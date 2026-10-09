/**
 * \file ERF_SolveWithGMRES.cpp
 */
#include "ERF.H"
#include "ERF_Utils.H"
#include "ERF_TerrainPoisson.H"
#include "ERF_SolverUtils.H"

#include <AMReX_GMRES.H>
#include <AMReX_Reduce.H>

using namespace amrex;

/**
 * Solve the Poisson equation using FFT-preconditioned GMRES
 *
 * @param lev Level index for the solve
 * @param subdomain Box over which the solve is performed
 * @param rhs Right-hand side field for the Poisson solve
 * @param p Solution field to fill
 * @param fluxes Face-centered gradient fluxes to fill
 * @param ax_sub Terrain metric coefficient on x-faces
 * @param ay_sub Terrain metric coefficient on y-faces
 * @param az_sub Terrain metric coefficient on z-faces
 * @param znd_sub Node-centered physical height field
 *
 * The unnamed metric argument in the class declaration is the cell-centered
 * Jacobian determinant.
 */
int ERF::solve_with_gmres (int lev, Long projection_call, double l_time, double l_dt,
                            const Box& subdomain, MultiFab& rhs, MultiFab& phi,
                            Array<MultiFab,AMREX_SPACEDIM>& fluxes,
                            MultiFab& ax_sub, MultiFab& ay_sub, MultiFab& az_sub,
                            MultiFab& dJ_sub, MultiFab& znd_sub,
                            Real& solver_residual)
{
#ifdef ERF_USE_FFT
    BL_PROFILE("ERF::solve_with_gmres()");

    Real reltol = solverChoice.poisson_reltol;
    Real abstol = solverChoice.poisson_abstol;

    auto const dom_lo = lbound(Geom(lev).Domain());
    auto const dom_hi = ubound(Geom(lev).Domain());

    auto const sub_lo = lbound(subdomain);
    auto const sub_hi = ubound(subdomain);

    auto dx    = Geom(lev).CellSizeArray();

    Geometry my_geom;

    Array<int,AMREX_SPACEDIM> is_per; is_per[0] = 0; is_per[1] = 0; is_per[2] = 0;
    if (Geom(lev).isPeriodic(0) && sub_lo.x == dom_lo.x && sub_hi.x == dom_hi.x) { is_per[0] = 1;}
    if (Geom(lev).isPeriodic(1) && sub_lo.y == dom_lo.y && sub_hi.y == dom_hi.y) { is_per[1] = 1;}

    int coord_sys = 0;

    // If subdomain == domain then we pass Geom(lev) to the FFT solver
    if (subdomain == Geom(lev).Domain()) {
        my_geom.define(Geom(lev).Domain(), Geom(lev).ProbDomain(), coord_sys, is_per);
    } else {
        // else we create a new geometry based only on the subdomain
        // The information in my_geom used by the FFT routines is:
        //   1) my_geom.Domain()
        //   2) my_geom.CellSize()
        //   3) my_geom.isAllPeriodic() / my_geom.periodicity()
        RealBox rb( sub_lo.x   *dx[0],  sub_lo.y   *dx[1],  sub_lo.z   *dx[2],
                   (sub_hi.x+1)*dx[0], (sub_hi.y+1)*dx[1], (sub_hi.z+1)*dx[2]);
        my_geom.define(subdomain, rb, coord_sys, is_per);
    }

    amrex::GMRES<MultiFab, TerrainPoisson> gmsolver;

    TerrainPoisson tp(my_geom, Geom(lev), rhs.boxArray(), rhs.DistributionMap(), domain_bc_type,
                      stretched_dz_d[lev], ax_sub, ay_sub, az_sub, dJ_sub, &znd_sub,
                      solverChoice.use_real_bcs);

    gmsolver.define(tp);

    gmsolver.setVerbose(mg_verbose);

    gmsolver.setRestartLength(50);

    tp.usePrecond(true);

    gmsolver.solve(phi, rhs, reltol, abstol);

    const int solve_status = gmsolver.getStatus();
    solver_residual = gmsolver.getResidualNorm();
    if (mg_verbose > 0) {
        MultiFab independent_residual(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
        tp.apply(independent_residual, phi);
        independent_residual.minus(rhs, 0, 1, 0);
        MultiFab abs_independent_residual(
            independent_residual.boxArray(), independent_residual.DistributionMap(), 1, 0);
        for (MFIter mfi(independent_residual); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            const auto src = independent_residual.const_array(mfi);
            const auto dst = abs_independent_residual.array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                dst(i,j,k) = amrex::Math::abs(src(i,j,k));
            });
        }
        const IntVect residual_max_cell = abs_independent_residual.maxIndex(0);
        int residual_max_rank = ParallelDescriptor::NProcs();
        for (MFIter mfi(abs_independent_residual); mfi.isValid(); ++mfi) {
            if (mfi.validbox().contains(residual_max_cell)) {
                residual_max_rank = ParallelDescriptor::MyProc();
            }
        }
        ParallelDescriptor::ReduceIntMin(residual_max_rank);
        amrex::Print() << "ERF_TERRAIN_GMRES_EVENT"
                       << " step=" << istep[lev]
                       << " projection_call=" << projection_call
                       << " time=" << l_time
                       << " dt=" << l_dt
                       << " level=" << lev
                       << " rank=" << ParallelDescriptor::MyProc()
                       << " region=" << terrain_mlmg_box_string(subdomain)
                       << " status=" << solve_status
                       << " iterations=" << gmsolver.getNumIters()
                       << " initial_residual=" << gmsolver.getInitialResidualNorm()
                       << " reported_residual=" << gmsolver.getResidualNorm()
                       << " independent_residual_linf=" << independent_residual.norm0()
                       << " independent_residual_max_cell=" << residual_max_cell
                       << " independent_residual_max_rank=" << residual_max_rank
                       << " independent_residual_l2=" << independent_residual.norm2()
                       << " reltol=" << reltol
                       << " abstol=" << abstol << std::endl;
    }

    tp.getFluxes(phi, fluxes);

    auto report_normal_flux = [&] (char const* stage) {
        if (mg_verbose <= 1) { return; }
        const char* direction_name[AMREX_SPACEDIM] = {"x", "y", "z"};
        for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
            for (int side = 0; side < 2; ++side) {
                const int face_index = side == 0
                    ? subdomain.smallEnd(idir) : subdomain.bigEnd(idir) + 1;
                Box face_region = convert(
                    subdomain, IntVect::TheDimensionVector(idir));
                face_region.setRange(idir, face_index, 1);

                ReduceOps<ReduceOpMax,ReduceOpSum> reduce_op;
                ReduceData<Real,Long> reduce_data(reduce_op);
                using ReduceTuple = typename decltype(reduce_data)::Type;
                const MultiFab& face_flux = fluxes[idir];
                for (MFIter flux_mfi(face_flux); flux_mfi.isValid(); ++flux_mfi) {
                    const Box bx = face_region & flux_mfi.validbox();
                    if (bx.isEmpty()) { continue; }
                    const Array4<Real const> flux_arr = face_flux.const_array(flux_mfi);
                    reduce_op.eval(bx, reduce_data,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
                    {
                        return {amrex::Math::abs(flux_arr(i,j,k)), Long(1)};
                    });
                }
                Real max_abs = amrex::get<0>(reduce_data.value(reduce_op));
                Long samples = amrex::get<1>(reduce_data.value(reduce_op));
                ParallelDescriptor::ReduceRealMax(max_abs);
                ParallelDescriptor::ReduceLongSum(samples);
                amrex::Print() << "Terrain GMRES normal correction flux"
                               << " step=" << istep[lev]
                               << " time=" << l_time
                               << " level=" << lev
                               << " direction=" << direction_name[idir]
                               << " side=" << (side == 0 ? "lo" : "hi")
                               << " face_index=" << face_index
                               << " stage=" << stage
                               << " samples=" << samples
                               << " max_abs=" << max_abs << std::endl;
            }
        }
    };
    report_normal_flux("raw_getFluxes");

    for (MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        Box xbx = mfi.nodaltilebox(0);
        Box ybx = mfi.nodaltilebox(1);
        const Array4<Real      >& fx_ar = fluxes[0].array(mfi);
        const Array4<Real      >& fy_ar = fluxes[1].array(mfi);
        const Array4<Real const>& mf_ux = mapfac[lev][MapFacType::u_x]->const_array(mfi);
        const Array4<Real const>& mf_vy = mapfac[lev][MapFacType::v_y]->const_array(mfi);
        ParallelFor(xbx,ybx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fx_ar(i,j,k) *= mf_ux(i,j,0);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fy_ar(i,j,k) *= mf_vy(i,j,0);
        });
    } // mfi
    report_normal_flux("after_map_scaling");
#else
    amrex::ignore_unused(lev, projection_call, l_time, l_dt, subdomain, rhs, phi, fluxes,
                         ax_sub, ay_sub, az_sub, dJ_sub, znd_sub, solver_residual);
#endif

    // ****************************************************************************
    // Impose bc's on pprime
    // ****************************************************************************
    ImposeBCsOnPhi(lev, phi, subdomain);
#ifdef ERF_USE_FFT
    return solve_status;
#else
    return -1;
#endif
}
