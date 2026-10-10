/**
 * \file ERF_PoissonSolve.cpp
 */
#include "ERF.H"
#include "ERF_Utils.H"
#include "ERF_SolverUtils.H"

#include <AMReX_MLMG.H>
#if (AMREX_SPACEDIM == 3)
#include <AMReX_MLTerrainPoisson.H>
#endif

using namespace amrex;

/**
 * Solve the Poisson equation using the MLMG solver.
 *
 * @param lev Level index.
 * @param[in] rhs Right hand side of the Poisson equation.
 * @param[out] p Solution field.
 * @param[out] fluxes Gradient fluxes.
 * @param geom Geometry used for inverse cell spacing.
 * @param ref_ratio Refinement ratios.
 * @param l_domain_bc_type Domain boundary condition names.
 * @param mg_verbose Verbosity level for the solver.
 * @param reltol Relative tolerance.
 * @param abstol Absolute tolerance.
 */
void
solve_with_mlmg    (int lev,
                    Vector<amrex::MultiFab>& rhs, Vector<MultiFab>& p,
                    Vector<amrex::Array<MultiFab,AMREX_SPACEDIM>>& fluxes,
                    const Geometry& geom,
                    const amrex::Vector<amrex::IntVect>& ref_ratio,
                    Array<std::string,2*AMREX_SPACEDIM> l_domain_bc_type,
                    int mg_verbose, Real reltol, Real abstol);
/**
 * Solve the Poisson equation with embedded boundary using the MLMG solver.
 *
 * @tparam T EB data type.
 * @param lev Level index.
 * @param[in] rhs Right hand side of the Poisson equation.
 * @param[out] p Solution field.
 * @param[out] fluxes Gradient fluxes.
 * @param ebfact EB factory.
 * @param ebfact_u EB data for x-faces.
 * @param ebfact_v EB data for y-faces.
 * @param ebfact_w EB data for z-faces.
 * @param geom Geometry used for inverse cell spacing.
 * @param ref_ratio Refinement ratios.
 * @param l_domain_bc_type Domain boundary condition names.
 * @param mg_verbose Verbosity level for the solver.
 * @param reltol Relative tolerance.
 * @param abstol Absolute tolerance.
 */
template <typename T>
void
solve_with_EB_mlmg (int lev,
                    Vector<amrex::MultiFab>& rhs, Vector<MultiFab>& p,
                    Vector<amrex::Array<MultiFab,AMREX_SPACEDIM>>& fluxes,
                    EBFArrayBoxFactory const& ebfact,
                    T const& ebfact_u,
                    T const& ebfact_v,
                    T const& ebfact_w,
                    const Geometry& geom,
                    const amrex::Vector<amrex::IntVect>& ref_ratio,
                    Array<std::string,2*AMREX_SPACEDIM> l_domain_bc_type,
                    int mg_verbose, Real reltol, Real abstol);

/**
 * Project the single-level velocity field to enforce the anelastic constraint
 * Note that the level may or may not be level zero
 *
 * @param lev Level index for the velocity projection
 * @param time Time at which coarse data are registered
 * @param l_dt Time step used for coarse data registration
 */
void ERF::project_initial_velocity (int lev, double time, double l_dt)
{
    BL_PROFILE("ERF::project_initial_velocity()");
    // Impose FillBoundary on density since we use it in the conversion of velocity to momentum
    vars_new[lev][Vars::cons].FillBoundary(geom[lev].periodicity());

    const MultiFab* c_vfrac = nullptr;
    if (solverChoice.terrain_type == TerrainType::EB) {
        c_vfrac = &((get_eb(lev).get_const_factory())->getVolFrac());
    }

    VelocityToMomentum(vars_new[lev][Vars::xvel], IntVect{0},
                       vars_new[lev][Vars::yvel], IntVect{0},
                       vars_new[lev][Vars::zvel], IntVect{0},
                       vars_new[lev][Vars::cons],
                       rU_new[lev], rV_new[lev], rW_new[lev],
                       Geom(lev).Domain(), domain_bcs_type, c_vfrac);

    Vector<MultiFab> tmp_mom;

    tmp_mom.push_back(MultiFab(vars_new[lev][Vars::cons],make_alias,0,1));
    tmp_mom.push_back(MultiFab(rU_new[lev],make_alias,0,1));
    tmp_mom.push_back(MultiFab(rV_new[lev],make_alias,0,1));
    tmp_mom.push_back(MultiFab(rW_new[lev],make_alias,0,1));

    // If at lev > 0 we must first fill the velocities at the c/f interface -- this must
    //    be done *after* the projection at lev-1
    if (lev > 0) {
        int levc = lev-1;

        const MultiFab* c_vfrac_crse = nullptr;
        if (solverChoice.terrain_type == TerrainType::EB) {
            c_vfrac_crse = &((get_eb(levc).get_const_factory())->getVolFrac());
        }

        MultiFab& S_new_crse = vars_new[levc][Vars::cons];
        MultiFab& U_new_crse = vars_new[levc][Vars::xvel];
        MultiFab& V_new_crse = vars_new[levc][Vars::yvel];
        MultiFab& W_new_crse = vars_new[levc][Vars::zvel];

        VelocityToMomentum(U_new_crse, IntVect{0}, V_new_crse, IntVect{0}, W_new_crse, IntVect{0}, S_new_crse,
                           rU_new[levc], rV_new[levc], rW_new[levc],
                           Geom(levc).Domain(), domain_bcs_type, c_vfrac_crse);

        rU_new[levc].FillBoundary(geom[levc].periodicity());
        FPr_u[levc].RegisterCoarseData({&rU_new[levc], &rU_new[levc]}, {time, time+l_dt});

        rV_new[levc].FillBoundary(geom[levc].periodicity());
        FPr_v[levc].RegisterCoarseData({&rV_new[levc], &rV_new[levc]}, {time, time+l_dt});

        rW_new[levc].FillBoundary(geom[levc].periodicity());
        FPr_w[levc].RegisterCoarseData({&rW_new[levc], &rW_new[levc]}, {time, time+l_dt});
    }

    // Use the same time that was registered in the FillPatcher above so that the
    // FillSet assertion (time >= crse_times[0] && time <= crse_times[1]) is satisfied
    // when called at non-zero simulation time (restart or mid-run regrid).
    project_momenta(lev, time, l_dt, tmp_mom);

    MomentumToVelocity(vars_new[lev][Vars::xvel],
                       vars_new[lev][Vars::yvel],
                       vars_new[lev][Vars::zvel],
                       vars_new[lev][Vars::cons],
                       rU_new[lev], rV_new[lev], rW_new[lev],
                       Geom(lev).Domain(), domain_bcs_type, c_vfrac);
 }

/**
 * Project the single-level momenta to enforce the anelastic constraint
 * Note that the level may or may not be level zero
 *
 * @param lev Level index for the momentum projection
 * @param l_time Time used for coarse-fine momentum fills
 * @param l_dt Time step used in the projection update
 * @param vars Conserved density and face-centered momenta to project
 */
void ERF::project_momenta (int lev, double l_time, double l_dt_d, Vector<MultiFab>& mom_mf)
{
    BL_PROFILE("ERF::project_momenta()");

    // Stable call sequence identifies the projection stage when time-integration
    // invokes this routine multiple times within one coarse timestep.
    static Long projection_call_counter = 0;
    const Long projection_call = projection_call_counter++;

    const bool terrain_poisson_timing = terrain_poisson_timing_enabled();
    double projection_start_time = 0.0;
    if (terrain_poisson_timing) {
        terrain_poisson_timing_sync();
        projection_start_time = ParallelDescriptor::second();
    }
    double coarse_fine_ghost_fill_time = 0.0;

    // Keep the existing GMRES+FFT solver as the default. The experimental
    // multigrid path is opt-in through erf.terrain_poisson_solver = mlmg.
    static const std::string terrain_poisson_solver = [] () {
        std::string choice = terrain_poisson_solver_default();
        ParmParse pp("erf");
        pp.query("terrain_poisson_solver", choice);
        const std::string error = terrain_poisson_solver_choice_error(choice);
        if (!error.empty()) {
            amrex::Abort(error);
        }
        return choice;
    } ();

    if (terrain_poisson_solver == "mlmg") {
        const std::string error = terrain_mlmg_configuration_error(
            AMREX_SPACEDIM == 3,
            solverChoice.mesh_type == MeshType::VariableDz,
            solverChoice.terrain_type == TerrainType::StaticFittedMesh,
            solverChoice.use_real_bcs);
        if (!error.empty()) {
            amrex::Abort(error);
        }
        // The operator and post-solve flux scaling have not been validated for
        // nonunity horizontal map factors. Static fitted-terrain map factors
        // are immutable after a level is created, so inspect each new level
        // once rather than performing several global extrema reductions on
        // every projection call.
        static int checked_map_factor_levels = 0;
        for (int check_lev = checked_map_factor_levels;
             check_lev < static_cast<int>(subdomains.size()); ++check_lev) {
            Real max_abs_deviation = 0.0;
            if (mapfac[check_lev].size() < MapFacType::num) {
                amrex::Abort("erf.terrain_poisson_solver=mlmg cannot inspect the horizontal map factors "
                             "at level " + std::to_string(check_lev));
            }
            for (int imap = 0; imap < MapFacType::num; ++imap) {
                if (!mapfac[check_lev][imap]) {
                    amrex::Abort("erf.terrain_poisson_solver=mlmg found a missing horizontal map factor "
                                 "at level " + std::to_string(check_lev));
                }
                const MultiFab& scale = *mapfac[check_lev][imap];
                const Real min_scale = scale.min(0);
                const Real max_scale = scale.max(0);
                const Real low_deviation = amrex::Math::abs(min_scale - Real(1.0));
                const Real high_deviation = amrex::Math::abs(max_scale - Real(1.0));
                if (low_deviation > max_abs_deviation) { max_abs_deviation = low_deviation; }
                if (high_deviation > max_abs_deviation) { max_abs_deviation = high_deviation; }
            }
            const std::string map_factor_error =
                terrain_mlmg_map_factor_error(max_abs_deviation);
            if (!map_factor_error.empty()) {
                amrex::Abort(map_factor_error + "; level=" + std::to_string(check_lev));
            }
        }
        checked_map_factor_levels = static_cast<int>(subdomains.size());

        // Validate every currently constructed solve region before this call
        // can solve any level. Checking only `lev` lets a coarse projection
        // complete before a later fine-level configuration is rejected.
        for (int check_lev = 0; check_lev < static_cast<int>(subdomains.size()); ++check_lev) {
            const std::string scope_error = terrain_mlmg_scope_error(
                maxLevel(), check_lev, subdomains[check_lev].size());
            if (!scope_error.empty()) {
                amrex::Abort(scope_error);
            }
            if (subdomains[check_lev].empty()) {
                amrex::Abort("erf.terrain_poisson_solver=mlmg found no ERF solve regions at level " +
                             std::to_string(check_lev));
            }
            const Box level_domain = geom[check_lev].Domain();
            for (int isub = 0; isub < subdomains[check_lev].size(); ++isub) {
                const Box region(subdomains[check_lev][isub].minimalBox());
                BoxList region_box_list;
                for (int igrid = 0; igrid < grids[check_lev].size(); ++igrid) {
                    if (subdomains[check_lev][isub].intersects(grids[check_lev][igrid])) {
                        region_box_list.push_back(grids[check_lev][igrid]);
                    }
                }
                const BoxArray region_boxes(region_box_list);
                const std::string region_error = terrain_mlmg_region_error(
                    check_lev, isub, level_domain, region, region_boxes);
                if (!region_error.empty()) {
                    amrex::Abort(region_error);
                }
                const std::string periodic_error = terrain_mlmg_periodic_region_error(
                    check_lev, isub, geom[check_lev], region);
                if (!periodic_error.empty()) {
                    amrex::Abort(periodic_error);
                }
            }
        }
    }

#ifndef ERF_USE_FFT
    if (terrain_poisson_solver == "gmres_fft" &&
        solverChoice.mesh_type == MeshType::VariableDz) {
        amrex::Abort("erf.terrain_poisson_solver=gmres_fft requires an FFT-enabled build; "
                     "set erf.terrain_poisson_solver=mlmg for the experimental non-FFT solver");
    }
#endif

    // Print the chosen terrain solver once, when multigrid verbosity is enabled.
    static bool printed_terrain_solver = false;
    if (solverChoice.mesh_type == MeshType::VariableDz && mg_verbose > 0 &&
        !printed_terrain_solver) {
        amrex::Print() << "Terrain Poisson solver: " << terrain_poisson_solver << std::endl;
        printed_terrain_solver = true;
    }

    Real l_dt = static_cast<Real>(l_dt_d);
    //
    // If at lev > 0 we must first fill the momenta at the c/f interface with interpolated coarse values
    //
    if (lev > 0) {
        PhysBCFunctNoOp null_bc;
        FPr_u[lev-1].FillSet(mom_mf[IntVars::xmom], l_time, null_bc, domain_bcs_type);
        FPr_v[lev-1].FillSet(mom_mf[IntVars::ymom], l_time, null_bc, domain_bcs_type);
        FPr_w[lev-1].FillSet(mom_mf[IntVars::zmom], l_time, null_bc, domain_bcs_type);
    }

    // Make sure the solver only sees the levels over which we are solving
    Vector<BoxArray>            ba_tmp;   ba_tmp.push_back(mom_mf[Vars::cons].boxArray());
    Vector<DistributionMapping> dm_tmp;   dm_tmp.push_back(mom_mf[Vars::cons].DistributionMap());
    Vector<Geometry>          geom_tmp; geom_tmp.push_back(geom[lev]);

    Box domain = geom[lev].Domain();

    MultiFab r_hse(base_state[lev], make_alias, BaseState::r0_comp, 1);

    Vector<MultiFab> rhs;
    Vector<MultiFab> phi;

    if (solverChoice.terrain_type == TerrainType::EB)
    {
        rhs.resize(1); rhs[0].define(ba_tmp[0], dm_tmp[0], 1, 0, MFInfo(), EBFactory(lev));
        phi.resize(1); phi[0].define(ba_tmp[0], dm_tmp[0], 1, 1, MFInfo(), EBFactory(lev));
    } else {
        rhs.resize(1); rhs[0].define(ba_tmp[0], dm_tmp[0], 1, 0);
        phi.resize(1); phi[0].define(ba_tmp[0], dm_tmp[0], 1, 1);
    }

    MultiFab rhs_lev(rhs[0], make_alias, 0, 1);
    MultiFab phi_lev(phi[0], make_alias, 0, 1);

    auto dx    = geom[lev].CellSizeArray();
    auto dxInv = geom[lev].InvCellSizeArray();

    // Inflow on an x-face -- note only the normal velocity is used in the projection
    if (domain_bc_type[0] == "Inflow" || domain_bc_type[3] == "Inflow") {
        (*physbcs_u[lev])(vars_new[lev][Vars::xvel],vars_new[lev][Vars::xvel],vars_new[lev][Vars::yvel],
                        IntVect{1,0,0},t_new[lev],BCVars::xvel_bc,false);
    }

    // Inflow on a  y-face -- note only the normal velocity is used in the projection
    if (domain_bc_type[1] == "Inflow" || domain_bc_type[4] == "Inflow") {
        (*physbcs_v[lev])(vars_new[lev][Vars::yvel],vars_new[lev][Vars::xvel],vars_new[lev][Vars::yvel],
                          IntVect{0,1,0},t_new[lev],BCVars::yvel_bc,false);
    }

    if (domain_bc_type[0] == "Inflow" || domain_bc_type[3] == "Inflow" ||
        domain_bc_type[1] == "Inflow" || domain_bc_type[4] == "Inflow") {

        const MultiFab* c_vfrac = nullptr;
        if (solverChoice.terrain_type == TerrainType::EB) {
            c_vfrac = &((get_eb(lev).get_const_factory())->getVolFrac());
        }

        VelocityToMomentum(vars_new[lev][Vars::xvel], IntVect{0},
                           vars_new[lev][Vars::yvel], IntVect{0},
                           vars_new[lev][Vars::zvel], IntVect{0},
                           vars_new[lev][Vars::cons],
                           mom_mf[IntVars::xmom],
                           mom_mf[IntVars::ymom],
                           mom_mf[IntVars::zmom],
                           Geom(lev).Domain(),
                           domain_bcs_type, c_vfrac);
    }

    // If !fixed_density, we must convert (rho u) which came in
    // to (rho0 u) which is what we will project
    if (!solverChoice.fixed_density[lev]) {
        ConvertForProjection(mom_mf[Vars::cons], r_hse,
                             mom_mf[IntVars::xmom],
                             mom_mf[IntVars::ymom],
                             mom_mf[IntVars::zmom],
                             Geom(lev).Domain(),
                             domain_bcs_type);
    }

    // ERF's face-state contract assigns duplicate shared faces to the lower
    // BoxArray index. Reconcile the rho0-weighted projection inputs before
    // forming divergence so a split layout or restart cannot feed different
    // copies of one physical face to neighboring cells. Keep the established
    // default GMRES path unchanged.
    if (terrain_poisson_solver == "mlmg") {
        mom_mf[IntVars::xmom].OverrideSync(geom[lev].periodicity());
        mom_mf[IntVars::ymom].OverrideSync(geom[lev].periodicity());
        mom_mf[IntVars::zmom].OverrideSync(geom[lev].periodicity());
    }

    // OmegaFromW and WFromOmega need rho0-weighted x/y momentum one cell above
    // and below each valid w-face. Same-level FillBoundary cannot populate the
    // coarse/fine part of that stencil. Build rho0-weighted coarse sources from
    // the states bracketing l_time, then use ERF's ordinary face-centered
    // FillPatchTwoLevels interpolation, and copy only vertical ghost faces back
    // so valid fine momentum and projection corrections remain untouched.
    MultiFab coarse_xmom_old, coarse_ymom_old, coarse_zmom_old;
    MultiFab coarse_xmom_new, coarse_ymom_new, coarse_zmom_new;
    auto fill_vertical_projection_ghost = [&] (MultiFab& fine_mom,
                                               MultiFab& coarse_mom_old,
                                               MultiFab& coarse_mom_new,
                                               int bc_comp,
                                               char const* stencil_name)
    {
        MultiFab interp_fine(fine_mom.boxArray(), fine_mom.DistributionMap(),
                             1, IntVect(0,0,1));
        interp_fine.setVal(bogus_large_value);
        Vector<MultiFab*> coarse_data{&coarse_mom_old, &coarse_mom_new};
        Vector<MultiFab*> fine_data{&fine_mom, &fine_mom};
        Vector<Real> coarse_times{static_cast<Real>(t_old[lev-1]),
                                  static_cast<Real>(t_new[lev-1])};
        Vector<Real> data_times{static_cast<Real>(l_time), static_cast<Real>(l_time)};
        if (mg_verbose > 1) {
            amrex::Print() << "ERF_TERRAIN_PROJECTION_GHOST_TIME"
                           << " level=" << lev
                           << " fine_time=" << l_time
                           << " coarse_old_time=" << t_old[lev-1]
                           << " coarse_new_time=" << t_new[lev-1]
                           << " stencil=" << stencil_name << std::endl;
        }
        FillPatchTwoLevels(interp_fine, IntVect(0,0,1), IntVect(0,0,0),
                           static_cast<Real>(l_time),
                           coarse_data, coarse_times, fine_data, data_times,
                           0, 0, 1, geom[lev-1], geom[lev], refRatio(lev-1),
                           &face_cons_linear_interp, domain_bcs_type, bc_comp);

        const Box level_face_domain = convert(geom[lev].Domain(), fine_mom.boxArray().ixType());

        if (mg_verbose > 1) {
            // Strict diagnostic mode records the first failed interpolation and
            // synchronizes only when the user explicitly requests stencil checks.
            Gpu::DeviceScalar<int> d_bad(0), d_bad_i(-1), d_bad_j(-1), d_bad_k(-1), d_bad_box(-1);
            Gpu::DeviceScalar<Real> d_bad_value(zero);
            int* bad_ptr = d_bad.dataPtr();
            int* bad_i_ptr = d_bad_i.dataPtr();
            int* bad_j_ptr = d_bad_j.dataPtr();
            int* bad_k_ptr = d_bad_k.dataPtr();
            int* bad_box_ptr = d_bad_box.dataPtr();
            Real* bad_value_ptr = d_bad_value.dataPtr();

            for (MFIter mfi(fine_mom, false); mfi.isValid(); ++mfi) {
                const Box& valid = mfi.validbox();
                Box ghost_lo(valid);
                ghost_lo.setRange(2, valid.smallEnd(2)-1, valid.smallEnd(2)-1);
                ghost_lo &= interp_fine[mfi].box();
                ghost_lo &= level_face_domain;
                Box ghost_hi(valid);
                ghost_hi.setRange(2, valid.bigEnd(2)+1, valid.bigEnd(2)+1);
                ghost_hi &= interp_fine[mfi].box();
                ghost_hi &= level_face_domain;

                const Array4<Real> dst = fine_mom.array(mfi);
                const Array4<Real const> src = interp_fine.const_array(mfi);
                const int box_id = mfi.index();
                auto copy_and_check = [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    const Real value = src(i,j,k);
                    dst(i,j,k) = value;
                    if (terrain_projection_ghost_value_invalid(value, bogus_large_value) &&
                        Gpu::Atomic::CAS(bad_ptr,0,1) == 0) {
                        bad_i_ptr[0] = i; bad_j_ptr[0] = j; bad_k_ptr[0] = k;
                        bad_box_ptr[0] = box_id;
                        bad_value_ptr[0] = value;
                    }
                };
                if (!ghost_lo.isEmpty()) { ParallelFor(ghost_lo, copy_and_check); }
                if (!ghost_hi.isEmpty()) { ParallelFor(ghost_hi, copy_and_check); }
            }
            Gpu::synchronize();
            if (d_bad.dataValue() != 0) {
                const IntVect bad_index(d_bad_i.dataValue(), d_bad_j.dataValue(), d_bad_k.dataValue());
                const Box& this_box = fine_mom.boxArray()[d_bad_box.dataValue()];
                bool in_any_fine_box = false;
                for (int ib = 0; ib < fine_mom.boxArray().size(); ++ib) {
                    in_any_fine_box = in_any_fine_box || fine_mom.boxArray()[ib].contains(bad_index);
                }
                const char* where = in_any_fine_box ? "fine-fine ghost" : "coarse-fine ghost";
                AllPrint() << "Invalid interpolated rho0 momentum ghost before " << stencil_name
                        << ": ERF level=" << lev
                        << ", box=" << d_bad_box.dataValue() << ' ' << this_box
                        << ", component=" << ((bc_comp == BCVars::xvel_bc) ? "xmom" : "ymom")
                        << ", face=(" << bad_index << "), value=" << d_bad_value.dataValue()
                        << ", location=" << where << std::endl;
                Abort(std::string("coarse/fine momentum interpolation left an invalid ") +
                      stencil_name + " stencil value");
            }
        } else {
            // The normal path only copies interpolated values; it avoids
            // diagnostic atomics and a device-to-host synchronization.
            for (MFIter mfi(fine_mom, false); mfi.isValid(); ++mfi) {
                const Box& valid = mfi.validbox();
                Box ghost_lo(valid);
                ghost_lo.setRange(2, valid.smallEnd(2)-1, valid.smallEnd(2)-1);
                ghost_lo &= interp_fine[mfi].box();
                ghost_lo &= level_face_domain;
                Box ghost_hi(valid);
                ghost_hi.setRange(2, valid.bigEnd(2)+1, valid.bigEnd(2)+1);
                ghost_hi &= interp_fine[mfi].box();
                ghost_hi &= level_face_domain;

                const Array4<Real> dst = fine_mom.array(mfi);
                const Array4<Real const> src = interp_fine.const_array(mfi);
                auto copy = [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    dst(i,j,k) = src(i,j,k);
                };
                if (!ghost_lo.isEmpty()) { ParallelFor(ghost_lo, copy); }
                if (!ghost_hi.isEmpty()) { ParallelFor(ghost_hi, copy); }
            }
        }
    };

    auto validate_projection_stencil = [&] (MultiFab& xmom, MultiFab& ymom,
                                            char const* stencil_name)
    {
        if (mg_verbose <= 1) { return; }
        Gpu::DeviceScalar<int> d_bad(0), d_bad_box(-1), d_bad_component(-1);
        Gpu::DeviceScalar<int> d_src_i(-1), d_src_j(-1), d_src_k(-1);
        Gpu::DeviceScalar<int> d_w_i(-1), d_w_j(-1), d_w_k(-1);
        Gpu::DeviceScalar<Real> d_bad_value(zero);
        int* bad_ptr = d_bad.dataPtr();
        int* box_ptr = d_bad_box.dataPtr();
        int* component_ptr = d_bad_component.dataPtr();
        int* src_i_ptr = d_src_i.dataPtr();
        int* src_j_ptr = d_src_j.dataPtr();
        int* src_k_ptr = d_src_k.dataPtr();
        int* w_i_ptr = d_w_i.dataPtr();
        int* w_j_ptr = d_w_j.dataPtr();
        int* w_k_ptr = d_w_k.dataPtr();
        Real* value_ptr = d_bad_value.dataPtr();

        for (MFIter mfi(rhs_lev,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Array4<Real const> u = xmom.const_array(mfi);
            const Array4<Real const> v = ymom.const_array(mfi);
            const int box_id = mfi.index();
            const auto check_sample = [=] AMREX_GPU_DEVICE
                (Real value, int si, int sj, int sk, int component,
                 int wi, int wj, int wk) noexcept
            {
                if (terrain_projection_ghost_value_invalid(value, bogus_large_value) &&
                    Gpu::Atomic::CAS(bad_ptr,0,1) == 0) {
                    box_ptr[0] = box_id;
                    component_ptr[0] = component;
                    src_i_ptr[0] = si; src_j_ptr[0] = sj; src_k_ptr[0] = sk;
                    w_i_ptr[0] = wi; w_j_ptr[0] = wj; w_k_ptr[0] = wk;
                    value_ptr[0] = value;
                }
            };
            Box w_faces = mfi.nodaltilebox(2);
            ParallelFor(w_faces, [=] AMREX_GPU_DEVICE
                (int i, int j, int k) noexcept
            {
                const int klo = (k == 0) ? k : k-1;
                const int khi = (k == 0) ? k+1 : k;
                check_sample(u(i+1,j,klo), i+1,j,klo, 0, i,j,k);
                check_sample(u(i+1,j,khi), i+1,j,khi, 0, i,j,k);
                check_sample(u(i  ,j,klo), i  ,j,klo, 0, i,j,k);
                check_sample(u(i  ,j,khi), i  ,j,khi, 0, i,j,k);
                check_sample(v(i,j+1,klo), i,j+1,klo, 1, i,j,k);
                check_sample(v(i,j+1,khi), i,j+1,khi, 1, i,j,k);
                check_sample(v(i,j  ,klo), i,j  ,klo, 1, i,j,k);
                check_sample(v(i,j  ,khi), i,j  ,khi, 1, i,j,k);
            });
        }
        Gpu::synchronize();
        if (d_bad.dataValue() != 0) {
            const int box_id = d_bad_box.dataValue();
            const int component = d_bad_component.dataValue();
            const IntVect source(d_src_i.dataValue(), d_src_j.dataValue(), d_src_k.dataValue());
            const IntVect target(d_w_i.dataValue(), d_w_j.dataValue(), d_w_k.dataValue());
            const MultiFab& source_mom = component == 0 ? xmom : ymom;
            const Box& source_box = source_mom.boxArray()[box_id];
            const Box source_domain = convert(
                geom[lev].Domain(), source_mom.boxArray().ixType());
            const char* location = "coarse-fine ghost";
            if (!source_domain.contains(source)) {
                location = "physical-boundary ghost";
            } else if (source_box.contains(source)) {
                location = "valid face";
            } else {
                for (int ib = 0; ib < source_mom.boxArray().size(); ++ib) {
                    if (source_mom.boxArray()[ib].contains(source)) {
                        location = "fine-fine ghost";
                        break;
                    }
                }
            }
            AllPrint() << "Invalid " << stencil_name << " momentum stencil input"
                    << ": ERF level=" << lev << " projection_call=" << projection_call
                    << " rank=" << ParallelDescriptor::MyProc()
                    << " box=" << box_id << ' ' << source_box
                    << " target_wface=" << target
                    << " component=" << (component == 0 ? "xmom" : "ymom")
                    << " source_face=" << source
                    << " value=" << d_bad_value.dataValue()
                    << " location=" << location << std::endl;
            Abort(std::string(stencil_name) +
                  " read an invalid horizontal momentum stencil value");
        }
    };

    double omega_ghost_fill_start = 0.0;
    if (solverChoice.mesh_type == MeshType::VariableDz) {
        if (terrain_poisson_timing) {
            terrain_poisson_timing_sync();
            omega_ghost_fill_start = ParallelDescriptor::second();
        }
        mom_mf[IntVars::xmom].FillBoundary(IntVect(0,0,1), geom[lev].periodicity());
        mom_mf[IntVars::ymom].FillBoundary(IntVect(0,0,1), geom[lev].periodicity());
    }

    if (solverChoice.mesh_type == MeshType::VariableDz && lev > 0) {
        const int crse_lev = lev - 1;
        const auto& crse_u_old = rU_old[crse_lev];
        const auto& crse_v_old = rV_old[crse_lev];
        const auto& crse_w_old = rW_old[crse_lev];
        const auto& crse_u_new = rU_new[crse_lev];
        const auto& crse_v_new = rV_new[crse_lev];
        const auto& crse_w_new = rW_new[crse_lev];

        coarse_xmom_old.define(crse_u_old.boxArray(), crse_u_old.DistributionMap(), 1,
                               crse_u_old.nGrowVect());
        coarse_ymom_old.define(crse_v_old.boxArray(), crse_v_old.DistributionMap(), 1,
                               crse_v_old.nGrowVect());
        coarse_zmom_old.define(crse_w_old.boxArray(), crse_w_old.DistributionMap(), 1,
                               crse_w_old.nGrowVect());
        coarse_xmom_new.define(crse_u_new.boxArray(), crse_u_new.DistributionMap(), 1,
                               crse_u_new.nGrowVect());
        coarse_ymom_new.define(crse_v_new.boxArray(), crse_v_new.DistributionMap(), 1,
                               crse_v_new.nGrowVect());
        coarse_zmom_new.define(crse_w_new.boxArray(), crse_w_new.DistributionMap(), 1,
                               crse_w_new.nGrowVect());

        MultiFab::Copy(coarse_xmom_old, crse_u_old, 0, 0, 1, 0);
        MultiFab::Copy(coarse_ymom_old, crse_v_old, 0, 0, 1, 0);
        MultiFab::Copy(coarse_zmom_old, crse_w_old, 0, 0, 1, 0);
        MultiFab::Copy(coarse_xmom_new, crse_u_new, 0, 0, 1, 0);
        MultiFab::Copy(coarse_ymom_new, crse_v_new, 0, 0, 1, 0);
        MultiFab::Copy(coarse_zmom_new, crse_w_new, 0, 0, 1, 0);

        MultiFab crse_rho0(base_state[crse_lev], make_alias, BaseState::r0_comp, 1);
        ConvertForProjection(vars_old[crse_lev][Vars::cons], crse_rho0,
                             coarse_xmom_old, coarse_ymom_old, coarse_zmom_old,
                             geom[crse_lev].Domain(), domain_bcs_type);
        ConvertForProjection(vars_new[crse_lev][Vars::cons], crse_rho0,
                             coarse_xmom_new, coarse_ymom_new, coarse_zmom_new,
                             geom[crse_lev].Domain(), domain_bcs_type);

        coarse_xmom_old.FillBoundary(geom[crse_lev].periodicity());
        coarse_ymom_old.FillBoundary(geom[crse_lev].periodicity());
        coarse_xmom_new.FillBoundary(geom[crse_lev].periodicity());
        coarse_ymom_new.FillBoundary(geom[crse_lev].periodicity());
        fill_vertical_projection_ghost(mom_mf[IntVars::xmom],
                                       coarse_xmom_old, coarse_xmom_new,
                                       BCVars::xvel_bc, "OmegaFromW");
        fill_vertical_projection_ghost(mom_mf[IntVars::ymom],
                                       coarse_ymom_old, coarse_ymom_new,
                                       BCVars::yvel_bc, "OmegaFromW");
    }
    if (solverChoice.mesh_type == MeshType::VariableDz && terrain_poisson_timing) {
        terrain_poisson_timing_sync();
        coarse_fine_ghost_fill_time +=
            ParallelDescriptor::second() - omega_ghost_fill_start;
    }

    //
    // ****************************************************************************
    // Now convert the rho0w MultiFab to hold Omega rather than rhow
    // ****************************************************************************
    //
    if (solverChoice.mesh_type == MeshType::VariableDz)
    {
        validate_projection_stencil(mom_mf[IntVars::xmom],
                                    mom_mf[IntVars::ymom], "OmegaFromW");
        // OmegaFromW below averages (rho0 u) and (rho0 v) over the faces below and above
        // each w-face, so at the lowest and highest w-face of a box it reads one face in the
        // z-ghost layer. Same-level FillBoundary handles box splits; for partial fine
        // levels the explicit two-level fill above supplies coarse/fine values in rho0
        // units. At the physical bottom k=0 is set to zero and nothing below it is read;
        // the physical top remains governed by VelocityToMomentum's existing extrapolation.
        AMREX_ALWAYS_ASSERT(mom_mf[IntVars::xmom].nGrowVect()[2] >= 1 &&
                            mom_mf[IntVars::ymom].nGrowVect()[2] >= 1);
        for (MFIter mfi(rhs_lev,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Array4<Real const>& rho0u_arr = mom_mf[IntVars::xmom].const_array(mfi);
            const Array4<Real const>& rho0v_arr = mom_mf[IntVars::ymom].const_array(mfi);
            const Array4<Real      >& rho0w_arr = mom_mf[IntVars::zmom].array(mfi);
            const Array4<Real const>& z_nd = z_phys_nd[lev]->const_array(mfi);
            const Array4<Real const>& mf_u = mapfac[lev][MapFacType::u_x]->const_array(mfi);
            const Array4<Real const>& mf_v = mapfac[lev][MapFacType::v_y]->const_array(mfi);

            Box tbz = mfi.nodaltilebox(2);
            ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (k == 0) {
                    rho0w_arr(i,j,k) = zero;
                } else {
                    Real rho0w = rho0w_arr(i,j,k);
                    rho0w_arr(i,j,k) = OmegaFromW(i,j,k,rho0w,
                                                  rho0u_arr,rho0v_arr,
                                                  mf_u,mf_v,z_nd,dxInv);
                }
            });
        } // mfi
    }

    // ****************************************************************************
    // Allocate fluxes
    // ****************************************************************************
    Vector<Array<MultiFab,AMREX_SPACEDIM> > fluxes;
    fluxes.resize(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (solverChoice.terrain_type == TerrainType::EB) {
            fluxes[0][idim].define(convert(ba_tmp[0], IntVect::TheDimensionVector(idim)), dm_tmp[0], 1, 0, MFInfo(), EBFactory(lev));
        } else {
            fluxes[0][idim].define(convert(ba_tmp[0], IntVect::TheDimensionVector(idim)), dm_tmp[0], 1, 0);
        }
        //
        // A subdomain whose RHS is already below poisson_abstol is skipped below, so its
        // fluxes are never written by a solve. They are still added to the momenta after
        // the loop, so they must start at zero (FArrayBox data is uninitialized unless
        // AMREX_DEBUG is set).
        //
        fluxes[0][idim].setVal(0.0);
    }

    // ****************************************************************************
    // Initialize phi to 0
    // (It is essential that we do this in order to fill the corners; these are never
    //  used but the Saxpy requires the values to be initialized.)
    // ****************************************************************************
    phi_lev.setVal(0.0);

    // ****************************************************************************
    // Break into subdomains
    // ****************************************************************************

    std::map<int,int> index_map;

    BoxArray ba(grids[lev]);

    Vector<MultiFab> rhs_sub; rhs_sub.resize(1);
    Vector<MultiFab> phi_sub; phi_sub.resize(1);
    Vector<Array<MultiFab,AMREX_SPACEDIM>> fluxes_sub; fluxes_sub.resize(1);
    Vector<int> solve_status(subdomains[lev].size(), -1);
    Vector<Real> solve_residual(subdomains[lev].size(), Real(-1));
    Vector<Real> compatibility_mean(subdomains[lev].size(), zero);
    Real pre_projection_linf = zero;

    MultiFab ax_sub, ay_sub, az_sub, dJ_sub, znd_sub;
    MultiFab mfmx_sub, mfmy_sub, mfvx_sub, mfuy_sub;

    Array<MultiFab,AMREX_SPACEDIM> rho0_u_sub;
    Array<MultiFab const*, AMREX_SPACEDIM> rho0_u_const;

    // If we are going to solve with MLMG then we do not need to break this into subdomains
    bool will_solve_with_mlmg = false;
    if (solverChoice.mesh_type == MeshType::ConstantDz) {
        will_solve_with_mlmg = true;
#ifdef ERF_USE_FFT
        if (use_fft) {
            bool all_boxes_ok = true;
            for (int isub = 0; isub < subdomains[lev].size(); ++isub) {
                Box my_region(subdomains[lev][isub].minimalBox());
                bool boxes_make_rectangle = (my_region.numPts() == subdomains[lev][isub].numPts());
                if (!boxes_make_rectangle) {
                    all_boxes_ok = false;
                }
            } // isub
            if (all_boxes_ok) {
                will_solve_with_mlmg = false;
            }
        } // use_fft
#else
        if (use_fft) {
            amrex::Warning("You set use_fft=true but didn't build with USE_FFT = TRUE; defaulting to MLMG");
        }
#endif
    } // No terrain or grid stretching

    for (int isub = 0; isub < subdomains[lev].size(); ++isub)
    {
        Box my_region(subdomains[lev][isub].minimalBox());
        BoxList bl_sub;
        Vector<int> dm_sub;

        for (int j = 0; j < ba.size(); j++)
        {
            if (subdomains[lev][isub].intersects(ba[j]))
            {
                //
                // Note that bl_sub.size() is effectively a counter which is
                // incremented above
                //
                // if (ParallelDescriptor::MyProc() == j) {
                // }
                index_map[bl_sub.size()] = j;

                bl_sub.push_back(grids[lev][j]);
                dm_sub.push_back(dmap[lev][j]);
            } // intersects
        } // loop over ba (j)

        BoxArray ba_sub(bl_sub);

        BoxList bl2d_sub = ba_sub.boxList();
        for (auto& b : bl2d_sub) {
            b.setRange(2,0);
        }
        BoxArray ba2d_sub(std::move(bl2d_sub));

        // Define MultiFabs that hold only the data in this particular subdomain
        if (solverChoice.terrain_type == TerrainType::EB) {
            if (ba_sub != ba) {
                amrex::Print() << "EB Solves with multiple regions is not yet supported" << std::endl;
            }
            rhs_sub[0].define(ba_sub, DistributionMapping(dm_sub), 1, rhs_lev.nGrowVect(), MFInfo{}.SetAlloc(false), EBFactory(lev));
            phi_sub[0].define(ba_sub, DistributionMapping(dm_sub), 1, phi_lev.nGrowVect(), MFInfo{}.SetAlloc(false), EBFactory(lev));

            mfmx_sub.define(ba2d_sub, DistributionMapping(dm_sub), 1, mapfac[lev][MapFacType::m_x]->nGrowVect(), MFInfo{}.SetAlloc(false), EBFactory(lev));
            mfmy_sub.define(ba2d_sub, DistributionMapping(dm_sub), 1, mapfac[lev][MapFacType::m_y]->nGrowVect(), MFInfo{}.SetAlloc(false), EBFactory(lev));
            mfvx_sub.define(convert(ba2d_sub, IntVect(0,1,0)), DistributionMapping(dm_sub), 1, mapfac[lev][MapFacType::v_x]->nGrowVect(), MFInfo{}.SetAlloc(false), EBFactory(lev));
            mfuy_sub.define(convert(ba2d_sub, IntVect(1,0,0)), DistributionMapping(dm_sub), 1, mapfac[lev][MapFacType::u_y]->nGrowVect(), MFInfo{}.SetAlloc(false), EBFactory(lev));
              dJ_sub.define(ba_sub, DistributionMapping(dm_sub), 1, detJ_cc[lev]->nGrowVect(), MFInfo{}.SetAlloc(false), EBFactory(lev));

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                fluxes_sub[0][idim].define(convert(ba_sub, IntVect::TheDimensionVector(idim)), DistributionMapping(dm_sub), 1,
                                                   IntVect::TheZeroVector(), MFInfo{}.SetAlloc(false), EBFactory(lev));
            }
            rho0_u_sub[0].define(convert(ba_sub, IntVect::TheDimensionVector(0)), DistributionMapping(dm_sub), 1,
                                         mom_mf[IntVars::xmom].nGrowVect(), MFInfo{}.SetAlloc(false), EBFactory(lev));
            rho0_u_sub[1].define(convert(ba_sub, IntVect::TheDimensionVector(1)), DistributionMapping(dm_sub), 1,
                                         mom_mf[IntVars::ymom].nGrowVect(), MFInfo{}.SetAlloc(false), EBFactory(lev));
            rho0_u_sub[2].define(convert(ba_sub, IntVect::TheDimensionVector(2)), DistributionMapping(dm_sub), 1,
                                         mom_mf[IntVars::zmom].nGrowVect(), MFInfo{}.SetAlloc(false), EBFactory(lev));
        } else {
            rhs_sub[0].define(ba_sub, DistributionMapping(dm_sub), 1, rhs_lev.nGrowVect(), MFInfo{}.SetAlloc(false));
            phi_sub[0].define(ba_sub, DistributionMapping(dm_sub), 1, phi_lev.nGrowVect(), MFInfo{}.SetAlloc(false));

            mfmx_sub.define(ba2d_sub, DistributionMapping(dm_sub), 1, mapfac[lev][MapFacType::m_x]->nGrowVect(), MFInfo{}.SetAlloc(false));
            mfmy_sub.define(ba2d_sub, DistributionMapping(dm_sub), 1, mapfac[lev][MapFacType::m_y]->nGrowVect(), MFInfo{}.SetAlloc(false));
            mfvx_sub.define(convert(ba2d_sub, IntVect(0,1,0)), DistributionMapping(dm_sub), 1, mapfac[lev][MapFacType::v_x]->nGrowVect(), MFInfo{}.SetAlloc(false));
            mfuy_sub.define(convert(ba2d_sub, IntVect(1,0,0)), DistributionMapping(dm_sub), 1, mapfac[lev][MapFacType::u_y]->nGrowVect(), MFInfo{}.SetAlloc(false));
              dJ_sub.define(ba_sub, DistributionMapping(dm_sub), 1, detJ_cc[lev]->nGrowVect(), MFInfo{}.SetAlloc(false));

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                fluxes_sub[0][idim].define(convert(ba_sub, IntVect::TheDimensionVector(idim)), DistributionMapping(dm_sub), 1,
                                                   IntVect::TheZeroVector(), MFInfo{}.SetAlloc(false));
            }
            rho0_u_sub[0].define(convert(ba_sub, IntVect::TheDimensionVector(0)), DistributionMapping(dm_sub), 1,
                                         mom_mf[IntVars::xmom].nGrowVect(), MFInfo{}.SetAlloc(false));
            rho0_u_sub[1].define(convert(ba_sub, IntVect::TheDimensionVector(1)), DistributionMapping(dm_sub), 1,
                                         mom_mf[IntVars::ymom].nGrowVect(), MFInfo{}.SetAlloc(false));
            rho0_u_sub[2].define(convert(ba_sub, IntVect::TheDimensionVector(2)), DistributionMapping(dm_sub), 1,
                                         mom_mf[IntVars::zmom].nGrowVect(), MFInfo{}.SetAlloc(false));
        }

        // Link the new MultiFabs to the FABs in the original MultiFabs (no copy required)
        for (MFIter mfi(rhs_sub[0]); mfi.isValid(); ++mfi)
        {
            int orig_index = index_map[mfi.index()];
            rhs_sub[0].setFab(mfi, FArrayBox(rhs_lev[orig_index], amrex::make_alias, 0, 1));
            phi_sub[0].setFab(mfi, FArrayBox(phi_lev[orig_index], amrex::make_alias, 0, 1));

            mfmx_sub.setFab(mfi, FArrayBox((*mapfac[lev][MapFacType::m_x])[orig_index], amrex::make_alias, 0, 1));
            mfmy_sub.setFab(mfi, FArrayBox((*mapfac[lev][MapFacType::m_y])[orig_index], amrex::make_alias, 0, 1));
            mfvx_sub.setFab(mfi, FArrayBox((*mapfac[lev][MapFacType::v_x])[orig_index], amrex::make_alias, 0, 1));
            mfuy_sub.setFab(mfi, FArrayBox((*mapfac[lev][MapFacType::u_y])[orig_index], amrex::make_alias, 0, 1));

            fluxes_sub[0][0].setFab(mfi,FArrayBox(fluxes[0][0][orig_index], amrex::make_alias, 0, 1));
            fluxes_sub[0][1].setFab(mfi,FArrayBox(fluxes[0][1][orig_index], amrex::make_alias, 0, 1));
            fluxes_sub[0][2].setFab(mfi,FArrayBox(fluxes[0][2][orig_index], amrex::make_alias, 0, 1));

            rho0_u_sub[0].setFab(mfi,FArrayBox(mom_mf[IntVars::xmom][orig_index], amrex::make_alias, 0, 1));
            rho0_u_sub[1].setFab(mfi,FArrayBox(mom_mf[IntVars::ymom][orig_index], amrex::make_alias, 0, 1));
            rho0_u_sub[2].setFab(mfi,FArrayBox(mom_mf[IntVars::zmom][orig_index], amrex::make_alias, 0, 1));
        }

        rho0_u_const[0] = &rho0_u_sub[0];
        rho0_u_const[1] = &rho0_u_sub[1];
        rho0_u_const[2] = &rho0_u_sub[2];

        if (solverChoice.mesh_type != MeshType::ConstantDz) {
            ax_sub.define(convert(ba_sub,IntVect(1,0,0)), DistributionMapping(dm_sub), 1,
                          ax[lev]->nGrowVect(), MFInfo{}.SetAlloc(false));
            ay_sub.define(convert(ba_sub,IntVect(0,1,0)), DistributionMapping(dm_sub), 1,
                          ay[lev]->nGrowVect(), MFInfo{}.SetAlloc(false));
            az_sub.define(convert(ba_sub,IntVect(0,0,1)), DistributionMapping(dm_sub), 1,
                          az[lev]->nGrowVect(), MFInfo{}.SetAlloc(false));
            znd_sub.define(convert(ba_sub,IntVect(1,1,1)), DistributionMapping(dm_sub), 1,
                           z_phys_nd[lev]->nGrowVect(), MFInfo{}.SetAlloc(false));

            for (MFIter mfi(rhs_sub[0]); mfi.isValid(); ++mfi) {
                int orig_index = index_map[mfi.index()];
                ax_sub.setFab(mfi, FArrayBox((*ax[lev])[orig_index], amrex::make_alias, 0, 1));
                ay_sub.setFab(mfi, FArrayBox((*ay[lev])[orig_index], amrex::make_alias, 0, 1));
                az_sub.setFab(mfi, FArrayBox((*az[lev])[orig_index], amrex::make_alias, 0, 1));
                znd_sub.setFab(mfi, FArrayBox((*z_phys_nd[lev])[orig_index], amrex::make_alias, 0, 1));
                dJ_sub.setFab(mfi, FArrayBox((*detJ_cc[lev])[orig_index], amrex::make_alias, 0, 1));
            }
        }

        if (solverChoice.terrain_type == TerrainType::EB) {
            for (MFIter mfi(rhs_sub[0]); mfi.isValid(); ++mfi) {
                int orig_index = index_map[mfi.index()];
                dJ_sub.setFab(mfi, FArrayBox((*detJ_cc[lev])[orig_index], amrex::make_alias, 0, 1));
            }
        }

        // ****************************************************************************
        // Compute divergence which will form RHS
        // Note that we replace "rho0w" with the contravariant momentum, Omega
        // ****************************************************************************

        compute_divergence(lev, rhs_sub[0], rho0_u_const, mfmx_sub, mfmy_sub, mfvx_sub, mfuy_sub,
                           ax_sub, ay_sub, dJ_sub, geom_tmp[0]);

        Real rhsnorm;

        // Max norm over the entire MultiFab
        rhsnorm = rhs_sub[0].norm0();
        pre_projection_linf = std::max(pre_projection_linf, rhsnorm);

        if (mg_verbose > 0) {
            bool local = false;
            Real sum = volWgtSumMF(lev,rhs_sub[0],0,dJ_sub,mfmx_sub,mfmy_sub,false,local);
            Print() << "Max/L2 norm of divergence before solve in subdomain " << isub << " at level " << lev << " : " << rhsnorm << " " <<
                        rhs_sub[0].norm2() << " and volume-weighted sum " << sum << std::endl;
        }

        if (lev == 0 && solverChoice.use_real_bcs)
        {
            // We always use VariableDz if use_real_bcs is true
            AMREX_ALWAYS_ASSERT(solverChoice.mesh_type == MeshType::VariableDz);

            // Note that we always impose the projections one level at a time so this will always be a vector of length 1
            Array<MultiFab*, AMREX_SPACEDIM> rho0_u_vec =
               {&mom_mf[IntVars::xmom], &mom_mf[IntVars::ymom], &mom_mf[IntVars::zmom]};
            Array<MultiFab*, AMREX_SPACEDIM> area_vec = {ax[lev].get(), ay[lev].get(), az[lev].get()};
            //
            // Modify ax,ay,ax to include the map factors as used in the divergence calculation
            // We do this here so that it is seen in the call to enforceInOutSolvability
            //
            for (MFIter mfi(rhs_lev); mfi.isValid(); ++mfi)
            {
                Box xbx = mfi.nodaltilebox(0);
                Box ybx = mfi.nodaltilebox(1);
                Box zbx = mfi.nodaltilebox(2);
                const Array4<Real      >& ax_ar = ax[lev]->array(mfi);
                const Array4<Real      >& ay_ar = ay[lev]->array(mfi);
                const Array4<Real      >& az_ar = az[lev]->array(mfi);
                const Array4<Real const>& mf_uy = mapfac[lev][MapFacType::u_y]->const_array(mfi);
                const Array4<Real const>& mf_vx = mapfac[lev][MapFacType::v_x]->const_array(mfi);
                const Array4<Real const>& mf_mx = mapfac[lev][MapFacType::m_x]->const_array(mfi);
                const Array4<Real const>& mf_my = mapfac[lev][MapFacType::m_y]->const_array(mfi);
                ParallelFor(xbx,ybx,zbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    ax_ar(i,j,k) /= mf_uy(i,j,0);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    ay_ar(i,j,k) /= mf_vx(i,j,0);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    az_ar(i,j,k) /= (mf_mx(i,j,0)*mf_my(i,j,0));
                });
            } // mfi

            if (mg_verbose > 0) {
                Print() << "Calling enforceInOutSolvability" << std::endl;
            }
            enforceInOutSolvability(lev, rho0_u_vec, area_vec, geom[lev]);

            //
            // Return ax,ay,ax to their original definition
            //
            for (MFIter mfi(rhs_lev); mfi.isValid(); ++mfi)
            {
                Box xbx = mfi.nodaltilebox(0);
                Box ybx = mfi.nodaltilebox(1);
                Box zbx = mfi.nodaltilebox(2);
                const Array4<Real      >& ax_ar = ax[lev]->array(mfi);
                const Array4<Real      >& ay_ar = ay[lev]->array(mfi);
                const Array4<Real      >& az_ar = az[lev]->array(mfi);
                const Array4<Real const>& mf_uy = mapfac[lev][MapFacType::u_y]->const_array(mfi);
                const Array4<Real const>& mf_vx = mapfac[lev][MapFacType::v_x]->const_array(mfi);
                const Array4<Real const>& mf_mx = mapfac[lev][MapFacType::m_x]->const_array(mfi);
                const Array4<Real const>& mf_my = mapfac[lev][MapFacType::m_y]->const_array(mfi);
                ParallelFor(xbx,ybx,zbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    ax_ar(i,j,k) *= mf_uy(i,j,0);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    ay_ar(i,j,k) *= mf_vx(i,j,0);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    az_ar(i,j,k) *= (mf_mx(i,j,0)*mf_my(i,j,0));
                });
            } // mfi

            compute_divergence(lev, rhs_lev, rho0_u_const, *mapfac[lev][MapFacType::m_x],
                               *mapfac[lev][MapFacType::m_y], *mapfac[lev][MapFacType::v_x],
                               *mapfac[lev][MapFacType::u_y], *ax[lev], *ay[lev],
                               *detJ_cc[lev], geom_tmp[0]);

            // Re-define max norm over the entire MultiFab
            rhsnorm = rhs_lev.norm0();
            pre_projection_linf = std::max(pre_projection_linf, rhsnorm);

            if (mg_verbose > 0)
            {
                bool local = false;
                Real sum = volWgtSumMF(lev,rhs_sub[0],0,dJ_sub,mfmx_sub,mfmy_sub,false,local);
                Print() << "Max/L2 norm of divergence before solve at level " << lev << " : " << rhsnorm << " " <<
                            rhs_lev.norm2() << " and volume-weighted sum " << sum << std::endl;
            }
        } // lev 0 && use_real_bcs

        // *******************************************************************************************
        // Enforce solvability only for operators whose effective BCs are singular.
        TerrainMLMGRegionBCs terrain_region_bcs;
        bool is_singular = true;
        if (terrain_poisson_solver == "mlmg") {
            terrain_region_bcs = terrain_mlmg_region_bcs(geom[lev], my_region, domain_bc_type);
            is_singular = terrain_region_bcs.is_singular;
        } else if (lev == 0) {
            if ( (domain_bc_type[0] == "Outflow" || domain_bc_type[0] == "Open") && !solverChoice.use_real_bcs ) is_singular = false;
            if ( (domain_bc_type[1] == "Outflow" || domain_bc_type[1] == "Open") && !solverChoice.use_real_bcs ) is_singular = false;
            if ( (domain_bc_type[3] == "Outflow" || domain_bc_type[3] == "Open") && !solverChoice.use_real_bcs ) is_singular = false;
            if ( (domain_bc_type[4] == "Outflow" || domain_bc_type[4] == "Open") && !solverChoice.use_real_bcs ) is_singular = false;
            if ( (domain_bc_type[5] == "Outflow" || domain_bc_type[5] == "Open")                               ) is_singular = false;
        } else if ( (domain_bc_type[5] == "Outflow" || domain_bc_type[5] == "Open") &&
                    (my_region.bigEnd(2) == domain.bigEnd(2)) ) {
            is_singular = false;
        }

        if (is_singular)
        {
            bool local = false;
            Real sum = volWgtSumMF(lev,rhs_sub[0],0,dJ_sub,mfmx_sub,mfmy_sub,false,local);

            Real vol;
            if (terrain_poisson_solver == "mlmg") {
                // Match the compatibility integral's terrain Jacobian and
                // map-factor weighting. AMReX applies detJ to the operator RHS
                // internally, so the Jacobian is included here exactly once.
                MultiFab unit_weight(rhs_sub[0].boxArray(), rhs_sub[0].DistributionMap(), 1, 0);
                unit_weight.setVal(1.0);
                vol = volWgtSumMF(lev, unit_weight, 0, dJ_sub, mfmx_sub, mfmy_sub, false, false);
            } else if (solverChoice.mesh_type == MeshType::ConstantDz) {
                vol = rhs_sub[0].boxArray().numPts() * dx[0] * dx[1] * dx[2];
            } else {
                vol = dJ_sub.sum() * dx[0] * dx[1] * dx[2];
            }

            sum /= vol;
            compatibility_mean[isub] = sum;

            if (terrain_poisson_solver == "mlmg" &&
                terrain_mlmg_compatibility_mean_exceeds(sum)) {
                std::ostringstream message;
                message << "ERF_TERRAIN_MLMG_INCOMPATIBLE_COMPATIBILITY: "
                        << "singular homogeneous-normal solve region has incompatible net flux; "
                        << "lev=" << lev << ", subdomain=" << isub
                        << ", region=" << terrain_mlmg_box_string(my_region)
                        << ", weighted_compatibility_mean=" << sum
                        << ", tolerance=" << terrain_mlmg_compatibility_tolerance()
                        << ". Adjust the physical fluxes or rectangular solve region so the "
                        << "singular operator receives a compatible right-hand side.";
                amrex::Abort(message.str());
            }

            for (MFIter mfi(rhs_sub[0]); mfi.isValid(); ++mfi)
            {
                rhs_sub[0][mfi.index()].template minus<RunOn::Device>(sum);
            }
            if (mg_verbose > 0) {
                amrex::Print() << " Subtracting " << sum << " from rhs in subdomain " << isub << std::endl;

                sum = volWgtSumMF(lev,rhs_sub[0],0,dJ_sub,mfmx_sub,mfmy_sub,false,local);
                Print() << "Sum after subtraction " << sum << " in subdomain " << isub << std::endl;
            }

        } // if is_singular

        rhsnorm = rhs_sub[0].norm0();

        // ****************************************************************************
        // No need to build the solver if RHS == 0
        // ****************************************************************************
        if (rhsnorm <= solverChoice.poisson_abstol) {
            solve_status[isub] = 0;
            solve_residual[isub] = rhsnorm;
            continue; // this subdomain only
        }

        std::unique_ptr<MultiFab> rhs_solve_diag;
        if (mg_verbose > 1) {
            rhs_solve_diag = std::make_unique<MultiFab>(
                rhs_sub[0].boxArray(), rhs_sub[0].DistributionMap(), 1, 0);
            MultiFab::Copy(*rhs_solve_diag, rhs_sub[0], 0, 0, 1, 0);
        }

        double start_step = ParallelDescriptor::second();

        if (mg_verbose > 0) {
            amrex::Print() << " Solving in subdomain " << isub << " of " << subdomains[lev].size() << " bins at level " << lev << std::endl;
        }

        if (solverChoice.mesh_type == MeshType::VariableDz) {
            //
            // Modify ax,ay,ax to include the map factors as used in the divergence calculation
            // We do this here to set the coefficients used in the stencil -- the extra factor
            // of the mapfac comes from the gradient
            //
            for (MFIter mfi(rhs_sub[0]); mfi.isValid(); ++mfi)
            {
                Box xbx = mfi.nodaltilebox(0);
                Box ybx = mfi.nodaltilebox(1);
                Box zbx = mfi.nodaltilebox(2);
                const Array4<Real      >& ax_ar = ax_sub.array(mfi);
                const Array4<Real      >& ay_ar = ay_sub.array(mfi);
                const Array4<Real      >& az_ar = az_sub.array(mfi);
                const Array4<Real const>& mf_ux = mapfac[lev][MapFacType::u_x]->const_array(mfi);
                const Array4<Real const>& mf_uy = mapfac[lev][MapFacType::u_y]->const_array(mfi);
                const Array4<Real const>& mf_vx = mapfac[lev][MapFacType::v_x]->const_array(mfi);
                const Array4<Real const>& mf_vy = mapfac[lev][MapFacType::v_y]->const_array(mfi);
                const Array4<Real const>& mf_mx = mapfac[lev][MapFacType::m_x]->const_array(mfi);
                const Array4<Real const>& mf_my = mapfac[lev][MapFacType::m_y]->const_array(mfi);
                ParallelFor(xbx,ybx,zbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    ax_ar(i,j,k) *= (mf_ux(i,j,0) / mf_uy(i,j,0));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    ay_ar(i,j,k) *= (mf_vy(i,j,0) / mf_vx(i,j,0));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    az_ar(i,j,k) /= (mf_mx(i,j,0)*mf_my(i,j,0));
                });
            } // mfi
        }

        if (solverChoice.terrain_type != TerrainType::EB) {

        // ****************************************************************************
        // No terrain or grid stretching
        // ****************************************************************************
        if (solverChoice.mesh_type == MeshType::ConstantDz) {
            if (will_solve_with_mlmg) {
                solve_with_mlmg(lev, rhs_sub, phi_sub, fluxes_sub, geom[lev], ref_ratio, domain_bc_type,
                                mg_verbose, solverChoice.poisson_reltol, solverChoice.poisson_abstol);
            } else {
#ifdef ERF_USE_FFT
                solve_with_fft(lev, isub, my_region, rhs_sub[0], phi_sub[0], fluxes_sub[0]);
#endif
            }
        } // No terrain or grid stretching
        // ****************************************************************************
        // Grid stretching (flat terrain)
        // ****************************************************************************
        else if (solverChoice.mesh_type == MeshType::StretchedDz) {
#ifndef ERF_USE_FFT
            amrex::Abort("Rebuild with USE_FFT = TRUE so you can use the FFT solver");
#else
            bool boxes_make_rectangle = (my_region.numPts() == subdomains[lev][isub].numPts());
            if (!boxes_make_rectangle) {
                amrex::Abort("FFT won't work unless the union of boxes is rectangular");
            } else {
                if (!use_fft) {
                    amrex::Warning("Using FFT even though you didn't set use_fft to true; it's the best choice");
                }
                solve_with_fft(lev, isub, my_region, rhs_sub[0], phi_sub[0], fluxes_sub[0]);
            }
#endif
        } // grid stretching

        // ****************************************************************************
        // General terrain
        // ****************************************************************************
        else if (solverChoice.mesh_type == MeshType::VariableDz) {
            if (terrain_poisson_solver == "mlmg") {
#if (AMREX_SPACEDIM == 3)
                if (terrain_poisson_timing) { terrain_poisson_timing_sync(); }
                const double setup_start = terrain_poisson_timing ? ParallelDescriptor::second() : 0.0;
                Vector<Geometry> ml_geom{terrain_mlmg_region_geometry(Geom(lev), my_region)};
                Vector<BoxArray> ml_grids{rhs_sub[0].boxArray()};
                Vector<DistributionMapping> ml_dmap{rhs_sub[0].DistributionMap()};

                LPInfo lpinfo;
                MLTerrainPoisson terrain_op(ml_geom, ml_grids, ml_dmap, lpinfo);
                terrain_op.setDomainBC(terrain_region_bcs.lo, terrain_region_bcs.hi);
                terrain_op.setLevelBC(0, nullptr);
                terrain_op.setZPhys(0, znd_sub);
                Array<MultiFab const*,AMREX_SPACEDIM> terrain_areas{
                    &ax_sub, &ay_sub, &az_sub
                };
                terrain_op.setAreas(0, terrain_areas);
                terrain_op.setDetJ(0, dJ_sub);

                MLMG mlmg(terrain_op);
                mlmg.setMaxIter(200);
                mlmg.setVerbose(mg_verbose);
                mlmg.setBottomVerbose(0);
                if (terrain_poisson_timing) { terrain_poisson_timing_sync(); }
                const double setup_elapsed = terrain_poisson_timing
                    ? ParallelDescriptor::second() - setup_start : 0.0;
                if (mg_verbose > 0) {
                    amrex::Print() << "Terrain MLMG region: ERF level " << lev
                                   << ", subdomain " << isub << ", bounds "
                                   << terrain_mlmg_box_string(my_region)
                                   << ", operator AMR level 0" << std::endl;
                }
                if (mg_verbose > 1) {
                    auto bc_name = [] (LinOpBCType bc) {
                        if (bc == LinOpBCType::Periodic) return "periodic";
                        if (bc == LinOpBCType::Neumann) return "Neumann";
                        if (bc == LinOpBCType::Dirichlet) return "Dirichlet";
                        return "other";
                    };
                    amrex::Print() << "Terrain MLMG effective BCs: lo=("
                                   << bc_name(terrain_region_bcs.lo[0]) << ','
                                   << bc_name(terrain_region_bcs.lo[1]) << ','
                                   << bc_name(terrain_region_bcs.lo[2]) << "), hi=("
                                   << bc_name(terrain_region_bcs.hi[0]) << ','
                                   << bc_name(terrain_region_bcs.hi[1]) << ','
                                   << bc_name(terrain_region_bcs.hi[2]) << "), singular="
                                   << terrain_region_bcs.is_singular << std::endl;
                }

                if (terrain_poisson_timing) { terrain_poisson_timing_sync(); }
                const double solve_start = terrain_poisson_timing ? ParallelDescriptor::second() : 0.0;
                mlmg.solve(GetVecOfPtrs(phi_sub), GetVecOfConstPtrs(rhs_sub),
                           solverChoice.poisson_reltol, solverChoice.poisson_abstol);
                if (terrain_poisson_timing) { terrain_poisson_timing_sync(); }
                const double solve_elapsed = terrain_poisson_timing
                    ? ParallelDescriptor::second() - solve_start : 0.0;
                solve_residual[isub] = mlmg.getFinalResidual();
                const Real mlmg_tolerance = amrex::max(
                    solverChoice.poisson_reltol * mlmg.getInitResidual(),
                    solverChoice.poisson_abstol);
                solve_status[isub] =
                    amrex::Math::isfinite(solve_residual[isub]) &&
                    solve_residual[isub] <= mlmg_tolerance ? 0 : 1;

                // compFlux is inherited from MLCellLinOpT and dispatches to
                // MLTerrainPoisson::FFlux for this operator.
                Array<MultiFab*,AMREX_SPACEDIM> terrain_fluxes{
                    &fluxes_sub[0][0], &fluxes_sub[0][1], &fluxes_sub[0][2]
                };
                if (terrain_poisson_timing) { terrain_poisson_timing_sync(); }
                const double flux_start = terrain_poisson_timing ? ParallelDescriptor::second() : 0.0;
                terrain_op.compFlux(0, terrain_fluxes, phi_sub[0],
                                    MLTerrainPoisson::Location::FaceCenter);
                if (terrain_poisson_timing) { terrain_poisson_timing_sync(); }
                const double flux_elapsed = terrain_poisson_timing
                    ? ParallelDescriptor::second() - flux_start : 0.0;

                // Match the map-factor treatment used by the GMRES terrain path.
                for (MFIter mfi(phi_sub[0]); mfi.isValid(); ++mfi)
                {
                    Box xbx = mfi.nodaltilebox(0);
                    Box ybx = mfi.nodaltilebox(1);
                    const Array4<Real      >& fx_ar = fluxes_sub[0][0].array(mfi);
                    const Array4<Real      >& fy_ar = fluxes_sub[0][1].array(mfi);
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
                }

                // Report the actual normal correction flux on effective Neumann
                // faces. This diagnostic is intentionally verbosity-gated so normal
                // projections do not pay for boundary reductions.
                if (mg_verbose > 1) {
                    const char* axis_name[AMREX_SPACEDIM] = {"x", "y", "z"};
                    for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
                        for (int side = 0; side < 2; ++side) {
                            const LinOpBCType bc = side == 0
                                ? terrain_region_bcs.lo[idir] : terrain_region_bcs.hi[idir];
                            if (bc != LinOpBCType::Neumann) { continue; }

                            ReduceOps<ReduceOpMax,ReduceOpSum> reduce_op;
                            ReduceData<Real,Long> reduce_data(reduce_op);
                            using ReduceTuple = typename decltype(reduce_data)::Type;
                            MultiFab const& face_flux = fluxes_sub[0][idir];
                            Box face_region = convert(
                                my_region, IntVect::TheDimensionVector(idir));
                            const int face_index = side == 0
                                ? my_region.smallEnd(idir) : my_region.bigEnd(idir) + 1;
                            face_region.setRange(idir, face_index, 1);
                            for (MFIter flux_mfi(face_flux); flux_mfi.isValid(); ++flux_mfi) {
                                Box bx = face_region & flux_mfi.validbox();
                                if (mg_verbose > 2) {
                                    amrex::Print() << "Terrain MLMG flux diagnostic boxes: plane="
                                                   << face_region << " valid=" << flux_mfi.validbox()
                                                   << " intersection=" << bx << std::endl;
                                }
                                if (bx.isEmpty()) { continue; }
                                const Array4<Real const> flux_arr = face_flux.const_array(flux_mfi);
                                reduce_op.eval(bx, reduce_data,
                                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
                                {
                                    return {amrex::Math::abs(flux_arr(i,j,k)), Long(1)};
                                });
                            }
                            Real max_normal_flux = amrex::get<0>(reduce_data.value(reduce_op));
                            Long sample_count = amrex::get<1>(reduce_data.value(reduce_op));
                            ParallelDescriptor::ReduceRealMax(max_normal_flux);
                            ParallelDescriptor::ReduceLongSum(sample_count);
                            if (sample_count == 0) { max_normal_flux = 0.0; }
                            amrex::Print() << "Terrain MLMG Neumann normal correction flux max: level="
                                           << lev << " subdomain=" << isub
                                           << " dir=" << axis_name[idir]
                                           << " side=" << (side == 0 ? "lo" : "hi")
                                           << " face_index=" << face_index
                                           << " operator_domain=" << terrain_op.Geom(0).Domain()
                                           << " samples=" << sample_count
                                           << " max_abs=" << max_normal_flux << std::endl;
                        }
                    }
                }

                ImposeBCsOnPhi(lev, phi_sub[0], my_region);
                if (terrain_poisson_timing) {
                    Real max_setup = static_cast<Real>(setup_elapsed);
                    Real max_solve = static_cast<Real>(solve_elapsed);
                    Real max_flux = static_cast<Real>(flux_elapsed);
                    ParallelDescriptor::ReduceRealMax(max_setup);
                    ParallelDescriptor::ReduceRealMax(max_solve);
                    ParallelDescriptor::ReduceRealMax(max_flux);
                    amrex::Print() << "ERF_TERRAIN_TIMING solver=mlmg level=" << lev
                                   << " subdomain=" << isub
                                   << " setup_s=" << max_setup
                                   << " solve_s=" << max_solve
                                   << " compFlux_s=" << max_flux << std::endl;
                }
                if (mg_verbose > 0) {
                    amrex::Print() << "Terrain MLMG level " << lev << " subdomain " << isub
                                   << ": iterations=" << mlmg.getNumIters()
                                   << ", final residual=" << mlmg.getFinalResidual()
                                   << ", setup=" << setup_elapsed << " s"
                                   << ", solve=" << solve_elapsed << " s"
                                   << ", compFlux=" << flux_elapsed << " s" << std::endl;
                }
#else
                amrex::Abort("MLTerrainPoisson is available only in 3D builds");
#endif
            } else {
#ifdef ERF_USE_FFT
                bool boxes_make_rectangle = (my_region.numPts() == subdomains[lev][isub].numPts());
                if (!boxes_make_rectangle) {
                    amrex::Abort("FFT preconditioner for GMRES won't work unless the union of boxes is rectangular");
                } else {
                    solve_status[isub] = solve_with_gmres(
                        lev, projection_call, l_time, l_dt_d, my_region, rhs_sub[0], phi_sub[0],
                        fluxes_sub[0], ax_sub, ay_sub, az_sub, dJ_sub, znd_sub,
                        solve_residual[isub]);
                }
#else
                amrex::Abort("Rebuild with USE_FFT = TRUE so you can use the FFT preconditioner for GMRES");
#endif
            }

            // Independently check that the correction fluxes satisfy the same ERF
            // terrain divergence used to form the solve RHS. This separates an
            // operator/flux mismatch from iterative convergence.
            if (mg_verbose > 1 && rhs_solve_diag) {
                MultiFab correction_divergence(
                    rhs_sub[0].boxArray(), rhs_sub[0].DistributionMap(), 1, 0);
                Array<MultiFab const*,AMREX_SPACEDIM> correction_fluxes{
                    &fluxes_sub[0][0], &fluxes_sub[0][1], &fluxes_sub[0][2]
                };
                compute_divergence(lev, correction_divergence, correction_fluxes,
                                   mfmx_sub, mfmy_sub, mfvx_sub, mfuy_sub,
                                   ax_sub, ay_sub, dJ_sub, geom_tmp[0]);
                MultiFab flux_operator_residual(
                    correction_divergence.boxArray(),
                    correction_divergence.DistributionMap(), 1, 0);
                MultiFab::Copy(flux_operator_residual, correction_divergence, 0, 0, 1, 0);
                MultiFab::Add(flux_operator_residual, *rhs_solve_diag, 0, 0, 1, 0);
                amrex::Print() << "Terrain projection flux/operator check"
                               << " level=" << lev
                               << " subdomain=" << isub
                               << " correction_div_Linf=" << correction_divergence.norm0()
                               << " solve_rhs_Linf=" << rhs_solve_diag->norm0()
                               << " flux_operator_residual_Linf=" << flux_operator_residual.norm0()
                               << " flux_operator_residual_L2=" << flux_operator_residual.norm2()
                               << std::endl;
            }

            //
            // Restore ax,ay,ax to their original definitions
            //
            for (MFIter mfi(rhs_lev); mfi.isValid(); ++mfi)
            {
                Box xbx = mfi.nodaltilebox(0);
                Box ybx = mfi.nodaltilebox(1);
                Box zbx = mfi.nodaltilebox(2);
                const Array4<Real      >& ax_ar = ax_sub.array(mfi);
                const Array4<Real      >& ay_ar = ay_sub.array(mfi);
                const Array4<Real      >& az_ar = az_sub.array(mfi);
                const Array4<Real const>& mf_ux = mapfac[lev][MapFacType::u_x]->const_array(mfi);
                const Array4<Real const>& mf_uy = mapfac[lev][MapFacType::u_y]->const_array(mfi);
                const Array4<Real const>& mf_vx = mapfac[lev][MapFacType::v_x]->const_array(mfi);
                const Array4<Real const>& mf_vy = mapfac[lev][MapFacType::v_y]->const_array(mfi);
                const Array4<Real const>& mf_mx = mapfac[lev][MapFacType::m_x]->const_array(mfi);
                const Array4<Real const>& mf_my = mapfac[lev][MapFacType::m_y]->const_array(mfi);
                ParallelFor(xbx,ybx,zbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    ax_ar(i,j,k) *= (mf_uy(i,j,0) / mf_ux(i,j,0));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    ay_ar(i,j,k) *= (mf_vx(i,j,0) / mf_vy(i,j,0));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    az_ar(i,j,k) *= (mf_mx(i,j,0)*mf_my(i,j,0));
                });
            } // mfi

        } // MeshType::VariableDz

        // ****************************************************************************
        // Print time in solve
        // ****************************************************************************
        double end_step = ParallelDescriptor::second();
        if (mg_verbose > 0) {
            amrex::Print() << "Time in solve " << end_step - start_step << std::endl;
        }

        } // not EB
    } // loop over subdomains (i)

    // ****************************************************************************
    // When using multigrid we can solve for all of the level at once, even if there
    //      are disjoint regions
    // ****************************************************************************
    if (solverChoice.terrain_type == TerrainType::EB) {
        double start_step_eb = ParallelDescriptor::second();
        solve_with_EB_mlmg(lev, rhs_sub, phi_sub, fluxes_sub,
                           *(get_eb(lev).get_const_factory()),
                           *(get_eb(lev).get_u_const_factory()),
                           *(get_eb(lev).get_v_const_factory()),
                           *(get_eb(lev).get_w_const_factory()),
                           geom[lev], ref_ratio, domain_bc_type,
                           mg_verbose, solverChoice.poisson_reltol, solverChoice.poisson_abstol);
        double end_step_eb = ParallelDescriptor::second();
        if (mg_verbose > 0) {
            amrex::Print() << "Time in solve " << end_step_eb - start_step_eb << std::endl;
        }
    }

    // ****************************************************************************
    // Subtract dt grad(phi) from the momenta (rho0u, rho0v, Omega)
    // ****************************************************************************
    // In diagnostic mode compare duplicate face values against ERF's
    // lower-BoxArray-index ownership rule before and after correction. Audit
    // momentum and the solver-produced flux separately to locate any mismatch;
    // an OverrideSync copy is used so this instrumentation never changes state.
    auto audit_shared_faces = [&] (MultiFab const& face_values,
                                    char const* component,
                                    char const* stage) {
        MultiFab owner_values(face_values.boxArray(),
                              face_values.DistributionMap(), 1, 0);
        MultiFab::Copy(owner_values, face_values, 0, 0, 1, 0);
        owner_values.OverrideSync(geom[lev].periodicity());
        MultiFab::Subtract(owner_values, face_values, 0, 0, 1, 0);
        const Real difference = owner_values.norm0();
        const Real magnitude = face_values.norm0();
        const Real tolerance = Real(64) * std::numeric_limits<Real>::epsilon() *
                               amrex::max(Real(1.0), magnitude);
        amrex::Print() << "ERF_TERRAIN_SHARED_FACE_AUDIT level=" << lev
                       << " component=" << component
                       << " stage=" << stage
                       << " difference_linf=" << difference
                       << " tolerance=" << tolerance << std::endl;
        return difference;
    };
    if (terrain_poisson_solver == "mlmg" && mg_verbose > 1) {
        audit_shared_faces(mom_mf[IntVars::xmom], "xmom", "pre_correction_momentum");
        audit_shared_faces(mom_mf[IntVars::ymom], "ymom", "pre_correction_momentum");
        audit_shared_faces(mom_mf[IntVars::zmom], "zmom", "pre_correction_momentum");
        audit_shared_faces(fluxes[0][0], "xmom", "correction_flux");
        audit_shared_faces(fluxes[0][1], "ymom", "correction_flux");
        audit_shared_faces(fluxes[0][2], "zmom", "correction_flux");
    }

    if (terrain_poisson_timing) { terrain_poisson_timing_sync(); }
    const double correction_start_time = terrain_poisson_timing
        ? ParallelDescriptor::second() : 0.0;
    MultiFab::Add(mom_mf[IntVars::xmom],fluxes[0][0],0,0,1,0);
    MultiFab::Add(mom_mf[IntVars::ymom],fluxes[0][1],0,0,1,0);
    MultiFab::Add(mom_mf[IntVars::zmom],fluxes[0][2],0,0,1,0);

    if (terrain_poisson_solver == "mlmg" && mg_verbose > 1) {
        audit_shared_faces(mom_mf[IntVars::xmom], "xmom", "post_correction_momentum");
        audit_shared_faces(mom_mf[IntVars::ymom], "ymom", "post_correction_momentum");
        audit_shared_faces(mom_mf[IntVars::zmom], "zmom", "post_correction_momentum");
    }

    // ****************************************************************************
    // Define gradp from fluxes -- note that fluxes is dt * change in Gp
    //   (weighted by map factor!)
    // ****************************************************************************
    MultiFab::Saxpy(gradp[lev][GpVars::gpx],-one/l_dt,fluxes[0][0],0,0,1,0);
    MultiFab::Saxpy(gradp[lev][GpVars::gpy],-one/l_dt,fluxes[0][1],0,0,1,0);
    MultiFab::Saxpy(gradp[lev][GpVars::gpz],-one/l_dt,fluxes[0][2],0,0,1,0);

    gradp[lev][GpVars::gpx].FillBoundary(geom_tmp[0].periodicity());
    gradp[lev][GpVars::gpy].FillBoundary(geom_tmp[0].periodicity());
    gradp[lev][GpVars::gpz].FillBoundary(geom_tmp[0].periodicity());
    if (terrain_poisson_timing) { terrain_poisson_timing_sync(); }
    const double correction_elapsed = terrain_poisson_timing
        ? ParallelDescriptor::second() - correction_start_time : 0.0;

    //
    // This call is only to verify the divergence after the solve
    // It is important we do this before computing the rho0w_arr from Omega back to rho0w
    //
    // ****************************************************************************
    // THIS IS SIMPLY VERIFYING THE DIVERGENCE AFTER THE SOLVE
    // ****************************************************************************
    //
    if (mg_verbose > 0)
    {
        rho0_u_const[0] = &mom_mf[IntVars::xmom];
        rho0_u_const[1] = &mom_mf[IntVars::ymom];
        rho0_u_const[2] = &mom_mf[IntVars::zmom];

        compute_divergence(lev, rhs_lev, rho0_u_const, *mapfac[lev][MapFacType::m_x],
                           *mapfac[lev][MapFacType::m_y], *mapfac[lev][MapFacType::v_x],
                           *mapfac[lev][MapFacType::u_y], *ax[lev], *ay[lev],
                           *detJ_cc[lev], geom_tmp[0]);

        bool local = false;
        Real sum = volWgtSumMF(lev,rhs_lev,0,*detJ_cc[lev],*mapfac[lev][MapFacType::m_x],*mapfac[lev][MapFacType::m_y],false,local);

        if (mg_verbose > 0) {
            const Real post_projection_linf = rhs_lev.norm0();
            MultiFab abs_projection_divergence(
                rhs_lev.boxArray(), rhs_lev.DistributionMap(), 1, 0);
            for (MFIter mfi(rhs_lev); mfi.isValid(); ++mfi) {
                const Box bx = mfi.validbox();
                const auto src = rhs_lev.const_array(mfi);
                const auto dst = abs_projection_divergence.array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    dst(i,j,k) = amrex::Math::abs(src(i,j,k));
                });
            }
            const IntVect post_max_cell = abs_projection_divergence.maxIndex(0);
            int post_max_rank = ParallelDescriptor::NProcs();
            for (MFIter mfi(abs_projection_divergence); mfi.isValid(); ++mfi) {
                if (mfi.validbox().contains(post_max_cell)) {
                    post_max_rank = ParallelDescriptor::MyProc();
                }
            }
            ParallelDescriptor::ReduceIntMin(post_max_rank);
            int all_solvers_status = 0;
            Real max_solver_residual = zero;
            Real max_compatibility_mean = zero;
            for (int isub = 0; isub < subdomains[lev].size(); ++isub) {
                if (solve_status[isub] != 0) { all_solvers_status = 1; }
                max_solver_residual = amrex::max(
                    max_solver_residual, amrex::Math::abs(solve_residual[isub]));
                max_compatibility_mean = amrex::max(
                    max_compatibility_mean, amrex::Math::abs(compatibility_mean[isub]));
            }
            const Real relative_reduction = pre_projection_linf > zero
                ? (pre_projection_linf - post_projection_linf) / pre_projection_linf
                : zero;
            amrex::Print() << "ERF_TERRAIN_PROJECTION_EVENT"
                           << " step=" << istep[lev]
                           << " projection_call=" << projection_call
                           << " time=" << l_time
                           << " dt=" << l_dt_d
                           << " level=" << lev
                           << " solver=" << terrain_poisson_solver
                           << " regions=" << subdomains[lev].size()
                           << " pre_Linf=" << pre_projection_linf
                           << " post_Linf=" << post_projection_linf
                           << " post_max_cell=" << post_max_cell
                           << " post_max_rank=" << post_max_rank
                           << " relative_reduction=" << relative_reduction
                           << " solve_status=" << all_solvers_status
                           << " solver_reltol=" << solverChoice.poisson_reltol
                           << " solver_abstol=" << solverChoice.poisson_abstol
                           << " residual=" << max_solver_residual
                           << " compatibility_mean_max_abs=" << max_compatibility_mean
                           << std::endl;
            Print() << "Max/L2 norm of divergence after  solve at level " << lev << " : " << post_projection_linf << " " <<
                        rhs_lev.norm2() << " and volume-weighted sum " << sum << std::endl;
        }

#if 0
         // FOR DEBUGGING ONLY
         for ( MFIter mfi(rhs_lev,TilingIfNotGPU()); mfi.isValid(); ++mfi)
         {
            const Array4<Real const>& rhs_arr = rhs_lev.const_array(mfi);
            Box bx = mfi.validbox();
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (std::abs(rhs_arr(i,j,k)) > Real(1.e-10)) {
                    amrex::AllPrint() << "RHS after solve at " <<
                                          IntVect(i,j,k) << " " << rhs_arr(i,j,k) << std::endl;
                }
            });
         } // mfi
#endif

    } // mg_verbose

    //
    // ****************************************************************************
    // Now convert the rho0w MultiFab back to holding (rho0w) rather than Omega
    // ****************************************************************************
    //
    if (solverChoice.mesh_type == MeshType::VariableDz)
    {
        // WFromOmega reads the same z-ghost faces of (rho0 u) and (rho0 v) as OmegaFromW above,
        // but the fluxes were added to the valid faces only. Fill the ghost faces again, or the
        // lowest and highest w-face of a box inside the domain (a BoxArray split in z) would be
        // converted back with the horizontal momenta from before the projection.
        double w_omega_ghost_fill_start = 0.0;
        if (terrain_poisson_timing) {
            terrain_poisson_timing_sync();
            w_omega_ghost_fill_start = ParallelDescriptor::second();
        }
        mom_mf[IntVars::xmom].FillBoundary(IntVect(0,0,1), geom[lev].periodicity());
        mom_mf[IntVars::ymom].FillBoundary(IntVect(0,0,1), geom[lev].periodicity());
        if (lev > 0) {
            fill_vertical_projection_ghost(mom_mf[IntVars::xmom],
                                           coarse_xmom_old, coarse_xmom_new,
                                           BCVars::xvel_bc, "WFromOmega");
            fill_vertical_projection_ghost(mom_mf[IntVars::ymom],
                                           coarse_ymom_old, coarse_ymom_new,
                                           BCVars::yvel_bc, "WFromOmega");
        }
        if (terrain_poisson_timing) {
            terrain_poisson_timing_sync();
            coarse_fine_ghost_fill_time +=
                ParallelDescriptor::second() - w_omega_ghost_fill_start;
        }
        validate_projection_stencil(mom_mf[IntVars::xmom],
                                    mom_mf[IntVars::ymom], "WFromOmega");

        for (MFIter mfi(mom_mf[Vars::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
             Box tbz = mfi.nodaltilebox(2);
             const Array4<Real      >& rho0u_arr = mom_mf[IntVars::xmom].array(mfi);
             const Array4<Real      >& rho0v_arr = mom_mf[IntVars::ymom].array(mfi);
             const Array4<Real      >& rho0w_arr = mom_mf[IntVars::zmom].array(mfi);
             const Array4<Real const>&      z_nd = z_phys_nd[lev]->const_array(mfi);
             const Array4<Real const>&      mf_u =  mapfac[lev][MapFacType::u_x]->const_array(mfi);
             const Array4<Real const>&      mf_v =  mapfac[lev][MapFacType::v_y]->const_array(mfi);
             ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                 Real omega = rho0w_arr(i,j,k);
                 rho0w_arr(i,j,k) = WFromOmega(i,j,k,omega,
                                               rho0u_arr,rho0v_arr,
                                               mf_u,mf_v,z_nd,dxInv);
             });
        } // mfi
    }

    // If !fixed_density, we must convert (rho0 u) back
    // to (rho0 u) which is what we will pass back out
    if (!solverChoice.fixed_density[lev]) {
        ConvertForProjection(r_hse, mom_mf[Vars::cons],
                             mom_mf[IntVars::xmom],
                             mom_mf[IntVars::ymom],
                             mom_mf[IntVars::zmom],
                             Geom(lev).Domain(),
                             domain_bcs_type);
    }

    // ****************************************************************************
    // Update pressure variable with phi -- note that phi is dt * change in pressure
    // ****************************************************************************
    MultiFab::Saxpy(pp_inc[lev], one/l_dt, phi_lev,0,0,1,1);

    if (terrain_poisson_timing) {
        terrain_poisson_timing_sync();
        Real max_ghost_fill = static_cast<Real>(coarse_fine_ghost_fill_time);
        Real max_correction = static_cast<Real>(correction_elapsed);
        Real max_projection = static_cast<Real>(
            ParallelDescriptor::second() - projection_start_time);
        ParallelDescriptor::ReduceRealMax(max_ghost_fill);
        ParallelDescriptor::ReduceRealMax(max_correction);
        ParallelDescriptor::ReduceRealMax(max_projection);
        amrex::Print() << "ERF_TERRAIN_TIMING solver=" << terrain_poisson_solver
                       << " level=" << lev
                       << " projection_call=" << projection_call
                       << " projection_momentum_ghost_fill_s=" << max_ghost_fill
                       << " correction_application_s=" << max_correction
                       << " total_projection_s=" << max_projection
                       << std::endl;
    }
}
