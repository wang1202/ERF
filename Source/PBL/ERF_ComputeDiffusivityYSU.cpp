#include "ERF_SurfaceLayer.H"
#include "ERF_DirectionSelector.H"
#include "ERF_Diffusion.H"
#include "ERF_Constants.H"
#include "ERF_TurbStruct.H"
#include "ERF_PBLModels.H"

using namespace amrex;

void
ComputeDiffusivityYSU (const MultiFab& xvel,
                       const MultiFab& yvel,
                       const MultiFab& cons_in,
                       MultiFab& eddyViscosity,
                       const Geometry& geom,
                       const TurbChoice& turbChoice,
                       std::unique_ptr<SurfaceLayer>& SurfLayer,
                       bool use_terrain_fitted_coords,
                       bool /*use_moisture*/,
                       int level,
                       const BCRec* bc_ptr,
                       bool /*vert_only*/,
                       const std::unique_ptr<MultiFab>& z_phys_nd,
                       const MoistureComponentIndices& moisture_indices)
{
    {
        /*
          YSU PBL initially introduced by S.-Y. Hong, Y. Noh, and J. Dudhia, MWR, 2006 [HND06]

          Further Modifications from S.-Y. Hong, Q. J. R. Meteorol. Soc., 2010 [H10]

          Implementation follows WRF as of early 2024 with some simplifications
        */

        // Extract zref BEFORE the MFIter loop (get_zref calls MultiFab::min internally)
        const Real most_zref = SurfLayer->get_zref(level);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(eddyViscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            // Pull out the box we're working on, make sure it covers full domain in z-direction
            const Box &bx = mfi.growntilebox(1);
            const Box &dbx = geom.Domain();
            Box sbx(bx.smallEnd(), bx.bigEnd());
            sbx.grow(2,-1);
            AMREX_ALWAYS_ASSERT(sbx.smallEnd(2) == dbx.smallEnd(2) && sbx.bigEnd(2) == dbx.bigEnd(2));

            // Get some data in arrays
            const auto& cell_data = cons_in.const_array(mfi);
            const auto& uvel = xvel.const_array(mfi);
            const auto& vvel = yvel.const_array(mfi);

            const auto& z0_arr        = SurfLayer->get_z0(level)->const_array(mfi);
            const auto& ws10av_arr    = SurfLayer->get_mac_avg(level,5)->const_array(mfi);
            const auto& t10av_arr     = SurfLayer->get_mac_avg(level,2)->const_array(mfi);
            const auto& t_surf_arr    = SurfLayer->get_t_surf(level)->const_array(mfi);
            const auto& over_land_arr = (SurfLayer->get_lmask(level)) ? SurfLayer->get_lmask(level)->const_array(mfi) :
                                                                      Array4<int> {};
            const auto& u_star_arr = SurfLayer->get_u_star(level)->const_array(mfi);
            const auto& l_obuk_arr = SurfLayer->get_olen(level)->const_array(mfi);
            const auto& t_star_arr = SurfLayer->get_t_star(level)->const_array(mfi);
            const Array4<Real const> z_nd_arr = z_phys_nd->array(mfi);

            // Note: MOST zref (~10m) is used by the surface layer model to compute ws10av/t10av.
            // With terrain, zref adjusts to cell center height above local terrain.
            // The YSU scheme uses MOST-provided surface averages, so exact zref doesn't matter here.

            // create flattened boxes to store PBL height
            const GeometryData gdata = geom.data();
            const Box xybx = PerpendicularBox<ZDir>(bx, IntVect{0,0,0});
            FArrayBox pbl_height(xybx,1);
            IArrayBox pbl_index(xybx,1);
            const auto& pblh_arr = pbl_height.array();
            const auto& pbli_arr = pbl_index.array();

            // -- Diagnose PBL height - starting out assuming non-moist --
            // loop is only over i,j in order to find height at each x,y
            const Real f0 = turbChoice.pbl_ysu_coriolis_freq;
            const bool force_over_water = turbChoice.pbl_ysu_force_over_water;
            const Real land_Ribcr = turbChoice.pbl_ysu_land_Ribcr;
            const Real unst_Ribcr = turbChoice.pbl_ysu_unst_Ribcr;

            // Additional arrays for convective quantities (per-column)
            FArrayBox pbl_sfcflg(xybx,1);  // 1 = unstable, 0 = stable
            FArrayBox pbl_wstar3(xybx,1);
            FArrayBox pbl_hgamt(xybx,1);   // countergradient for theta
            FArrayBox pbl_hgamq(xybx,1);   // countergradient for moisture
            pbl_sfcflg.setVal<RunOn::Device>(0.0);
            pbl_wstar3.setVal<RunOn::Device>(0.0);
            pbl_hgamt.setVal<RunOn::Device>(0.0);
            pbl_hgamq.setVal<RunOn::Device>(0.0);
            const auto& sfcflg_arr = pbl_sfcflg.array();
            const auto& wstar3_arr = pbl_wstar3.array();
            const auto& hgamt_arr  = pbl_hgamt.array();
            const auto& hgamq_arr  = pbl_hgamq.array();

            ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                // Reconstruct a surface bulk Richardson number from the surface layer model
                // In WRF, this value is supplied to YSU by the MM5 surface layer model
                const Real t_surf = t_surf_arr(i,j,0);
                const Real t_layer = t10av_arr(i,j,0);
                const Real ws_layer = ws10av_arr(i,j,0);
                const Real Rib_layer = CONST_GRAV * most_zref / (ws_layer*ws_layer) * (t_layer - t_surf)/(t_layer);

                // Determine stability flag (WRF: br > 0 => stable)
                bool is_unstable = (Rib_layer < 0.0);
                sfcflg_arr(i,j,0) = is_unstable ? 1.0 : 0.0;

                // Surface virtual heat flux (WRF: sflux = hfx/rho/cp + ep1*th1*qfx/rho)
                // We approximate using t_star and u_star: sflux ~ -u* * t_star (kinematic)
                // For buoyancy flux, we use t_star as proxy
                const Real rho_sfc = cell_data(i,j,0,Rho_comp);
                const Real theta_sfc = cell_data(i,j,0,RhoTheta_comp) / rho_sfc;

                // Thermal for PBL height computation
                Real thermal = theta_sfc; // base theta for Rib calculation

                // For unstable: add thermal excess to virtual potential temperature
                // Following WRF YSU: thermal = thvx(1) + vpert
                if (is_unstable) {
                    // Compute preliminary wscale and thermal excess (WRF approach)
                    // hol1 = zol1 * hpbl / zl1 * sfcfrac (first guess)
                    constexpr Real sfcfrac = 0.1;
                    constexpr Real bfac = 6.8;
                    constexpr Real gamcrt = 3.0;
                    constexpr Real gamcrq = 2.0e-3;
                    constexpr Real phifac = 8.0;
                    constexpr Real aphi16 = 16.0;
                    constexpr Real aphi5 = 5.0;
                    constexpr Real h1 = 1.0/3.0;

                    // First guess PBL height using Rib = brcr_ub (unstable Ri_cr = 0)
                    // Just do a quick scan
                    Real hpbl_guess = 500.0; // first guess
                    int kpbl_guess = 0;
                    bool found = false;
                    for (int kk = 1; kk < bx.bigEnd(2); kk++) {
                        if (!bx.contains(i,j,kk)) break;
                        const Real zval_k = use_terrain_fitted_coords ?
                            Compute_Zrel_AtCellCenter(i,j,kk,z_nd_arr) : gdata.ProbLo(2) + (kk + 0.5)*gdata.CellSize(2);
                        const Real ws2_k = 0.25*( (uvel(i,j,kk)+uvel(i+1,j,kk))*(uvel(i,j,kk)+uvel(i+1,j,kk))
                                                + (vvel(i,j,kk)+vvel(i,j+1,kk))*(vvel(i,j,kk)+vvel(i,j+1,kk)) );
                        const Real theta_k = cell_data(i,j,kk,RhoTheta_comp) / cell_data(i,j,kk,Rho_comp);
                        const Real Rib_k = (theta_k - thermal) / thermal * CONST_GRAV * zval_k / amrex::max(ws2_k, 1.0);
                        if (Rib_k >= 0.0) {
                            hpbl_guess = zval_k;
                            kpbl_guess = kk;
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        hpbl_guess = use_terrain_fitted_coords ?
                            Compute_Zrel_AtCellCenter(i,j,bx.bigEnd(2)-1,z_nd_arr) :
                            gdata.ProbLo(2) + (bx.bigEnd(2)-1 + 0.5)*gdata.CellSize(2);
                    }

                    // Compute buoyancy flux and wstar
                    // sflux ~ -u* * t_star (MOST relation, kinematic theta flux)
                    const Real sflux = amrex::max(-u_star_arr(i,j,0) * t_star_arr(i,j,0), 0.0);
                    const Real govrth = CONST_GRAV / theta_sfc;
                    const Real ws3 = govrth * sflux * hpbl_guess;
                    wstar3_arr(i,j,0) = ws3;
                    const Real wstar_val = std::pow(amrex::max(ws3, 0.0), h1);

                    // zol1 = z1/L (stability parameter at first level)
                    const Real zl1 = use_terrain_fitted_coords ?
                        Compute_Zrel_AtCellCenter(i,j,0,z_nd_arr) : gdata.ProbLo(2) + 0.5*gdata.CellSize(2);
                    Real zol1 = zl1 / l_obuk_arr(i,j,0);
                    zol1 = amrex::min(zol1, -1.0e-6); // ensure negative for unstable

                    const Real hol1 = zol1 * hpbl_guess / zl1 * sfcfrac;
                    const Real phim = std::pow(1.0 - aphi16*hol1, -0.25);
                    const Real phih = std::pow(1.0 - aphi16*hol1, -0.5);

                    const Real ust3 = u_star_arr(i,j,0) * u_star_arr(i,j,0) * u_star_arr(i,j,0);
                    const Real wscale = std::pow(ust3 + phifac*KAPPA*ws3*0.5, h1);
                    const Real wscale_lim = amrex::max(amrex::min(wscale, u_star_arr(i,j,0)*aphi16),
                                                        u_star_arr(i,j,0)/aphi5);

                    // Countergradient terms (WRF)
                    const Real gamfac = bfac / rho_sfc / wscale_lim;
                    // hfx/rho/cp ~ -u* * t_star  (kinematic heat flux)
                    const Real hfx_kinematic = amrex::max(-u_star_arr(i,j,0) * t_star_arr(i,j,0), 0.0);
                    hgamt_arr(i,j,0) = amrex::min(gamfac * hfx_kinematic, gamcrt);
                    hgamt_arr(i,j,0) = amrex::max(hgamt_arr(i,j,0), 0.0);

                    // Moisture countergradient
                    // qfx/rho ~ -u* * q_star
                    // hgamq = min(gamfac * ep1 * theta * qfx/rho, gamcrq)
                    // For now, set to 0 (conservative)
                    hgamq_arr(i,j,0) = 0.0;

                    // Add thermal excess (WRF: vpert)
                    constexpr Real afac = 6.8; // same as bfac in WRF
                    const Real vpert = hgamt_arr(i,j,0) / bfac * afac;
                    thermal = thermal + amrex::max(vpert, 0.0) * amrex::min(zl1 / (sfcfrac * hpbl_guess), 1.0);
                }

                // PBL Height: Using (possibly enhanced) thermal
                Real Rib_cr;
                bool over_land = (over_land_arr) ? over_land_arr(i,j,0) : 1;
                if (is_unstable) {
                    // Unstable: use 0 as critical Ri (WRF: brcr_ub = 0)
                    Rib_cr = 0.0;
                } else if (over_land && !force_over_water) {
                    Rib_cr = land_Ribcr;
                } else { // over water
                    const Real z0 = z0_arr(i,j,0);
                    const Real Rossby = ws_layer/(f0*z0);
                    Rib_cr = min(0.16*std::pow(1.0e-7*Rossby,-0.18),0.3);
                }

                bool above_critical = false;
                int kpbl = 0;
                Real Rib_up = Rib_layer, Rib_dn;
                while (!above_critical and bx.contains(i,j,kpbl+1)) {
                    kpbl += 1;
                    const Real zval = use_terrain_fitted_coords ?
                                      Compute_Zrel_AtCellCenter(i,j,kpbl,z_nd_arr) : gdata.ProbLo(2) + (kpbl + 0.5)*gdata.CellSize(2);
                    const Real ws2_level = amrex::max(
                        0.25*( (uvel(i,j,kpbl)+uvel(i+1,j  ,kpbl))*(uvel(i,j,kpbl)+uvel(i+1,j  ,kpbl))
                             + (vvel(i,j,kpbl)+vvel(i  ,j+1,kpbl))*(vvel(i,j,kpbl)+vvel(i  ,j+1,kpbl)) ),
                        1.0);
                    const Real theta = cell_data(i,j,kpbl,RhoTheta_comp) / cell_data(i,j,kpbl,Rho_comp);
                    Rib_dn = Rib_up;
                    Rib_up = (theta - thermal) / thermal * CONST_GRAV * zval / ws2_level;
                    above_critical = Rib_up >= Rib_cr;
                }

                Real interp_fact;
                if (Rib_dn >= Rib_cr) {
                    interp_fact = 0.0;
                } else if (Rib_up <= Rib_cr)
                    interp_fact = 1.0;
                else {
                    interp_fact = (Rib_cr - Rib_dn) / (Rib_up - Rib_dn);
                }

                const Real zval_up = use_terrain_fitted_coords ?
                                     Compute_Zrel_AtCellCenter(i,j,kpbl,z_nd_arr) : gdata.ProbLo(2) + (kpbl + 0.5)*gdata.CellSize(2);
                const Real zval_dn = use_terrain_fitted_coords ?
                                     Compute_Zrel_AtCellCenter(i,j,amrex::max(kpbl-1,0),z_nd_arr) : gdata.ProbLo(2) + (amrex::max(kpbl-1,0) + 0.5)*gdata.CellSize(2);
                pblh_arr(i,j,0) = zval_dn + interp_fact*(zval_up-zval_dn);

                const Real zval_0 = use_terrain_fitted_coords ?
                                     Compute_Zrel_AtCellCenter(i,j,0,z_nd_arr) : gdata.ProbLo(2) + (0.5)*gdata.CellSize(2);
                const Real zval_1 = use_terrain_fitted_coords ?
                                     Compute_Zrel_AtCellCenter(i,j,1,z_nd_arr) : gdata.ProbLo(2) + (1.5)*gdata.CellSize(2);
                if (pblh_arr(i,j,0) < 0.5*(zval_0+zval_1) ) {
                    kpbl = 0;
                }
                // Minimum PBL height for unstable (WRF uses hpbl >= zval(1))
                if (is_unstable) {
                    pblh_arr(i,j,0) = amrex::max(pblh_arr(i,j,0), zval_0);
                }
                pbli_arr(i,j,0) = kpbl;

                // Recompute wstar3 with final PBL height for unstable
                if (is_unstable) {
                    const Real sflux = amrex::max(-u_star_arr(i,j,0) * t_star_arr(i,j,0), 0.0);
                    const Real govrth = CONST_GRAV / theta_sfc;
                    wstar3_arr(i,j,0) = govrth * sflux * pblh_arr(i,j,0);
                }
            });

            // -- Compute diffusion coefficients --

            const Array4<Real      > &K_turb = eddyViscosity.array(mfi);

            // Dirichlet flags to switch derivative stencil
            bool c_ext_dir_on_zlo = ( (bc_ptr[BCVars::cons_bc].lo(2) == ERFBCType::ext_dir) );
            bool c_ext_dir_on_zhi = ( (bc_ptr[BCVars::cons_bc].lo(5) == ERFBCType::ext_dir) );
            bool u_ext_dir_on_zlo = ( (bc_ptr[BCVars::xvel_bc].lo(2) == ERFBCType::ext_dir) );
            bool u_ext_dir_on_zhi = ( (bc_ptr[BCVars::xvel_bc].lo(5) == ERFBCType::ext_dir) );
            bool v_ext_dir_on_zlo = ( (bc_ptr[BCVars::yvel_bc].lo(2) == ERFBCType::ext_dir) );
            bool v_ext_dir_on_zhi = ( (bc_ptr[BCVars::yvel_bc].lo(5) == ERFBCType::ext_dir) );

            const auto& dxInv = geom.InvCellSizeArray();
            const Real dz_inv = geom.InvCellSize(2);
            const int izmin = geom.Domain().smallEnd(2);
            const int izmax = geom.Domain().bigEnd(2);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const Real zval = use_terrain_fitted_coords ?
                                  Compute_Zrel_AtCellCenter(i,j,k,z_nd_arr) : gdata.ProbLo(2) + (k + 0.5)*gdata.CellSize(2);
                const Real rho = cell_data(i,j,k,Rho_comp);
                const Real met_h_zeta = use_terrain_fitted_coords ? Compute_h_zeta_AtCellCenter(i,j,k,dxInv,z_nd_arr) : 1.0;
                const Real dz_terrain = met_h_zeta/dz_inv;
                const bool is_unstable = (sfcflg_arr(i,j,0) > 0.5);

                if (k < pbli_arr(i,j,0)) {
                    // -- Compute diffusion coefficients within PBL
                    constexpr Real zfacmin = 1e-8; // value from WRF
                    constexpr Real phifac = 8.0; // value from H10 and WRF
                    constexpr Real pfac = 2.0; // profile exponent
                    constexpr Real aphi16 = 16.0;
                    constexpr Real aphi5 = 5.0;
                    constexpr Real sfcfrac = 0.1;
                    constexpr Real h1 = 1.0/3.0;

                    const Real zfac = std::min(std::max(1 - zval / pblh_arr(i,j,0), zfacmin ), 1.0);
                    const Real ust3 = u_star_arr(i,j,0) * u_star_arr(i,j,0) * u_star_arr(i,j,0);

                    if (is_unstable) {
                        // Unstable case (WRF YSU): use full wscale with wstar
                        const Real ws3 = wstar3_arr(i,j,0);
                        Real wscalek = std::pow(ust3 + phifac * KAPPA * ws3 * (1.0 - zfac), h1);

                        // Compute Prandtl number profile for unstable (WRF approach)
                        const Real zl1 = use_terrain_fitted_coords ?
                            Compute_Zrel_AtCellCenter(i,j,0,z_nd_arr) : gdata.ProbLo(2) + 0.5*gdata.CellSize(2);
                        Real zol1 = zl1 / l_obuk_arr(i,j,0);
                        zol1 = amrex::min(zol1, -1.0e-6);
                        const Real hol1 = zol1 * pblh_arr(i,j,0) / zl1 * sfcfrac;
                        const Real phim = std::pow(1.0 - aphi16*hol1, -0.25);
                        const Real phih = std::pow(1.0 - aphi16*hol1, -0.5);

                        // prnum0 = phih/phim + bfac*kappa*sfcfrac (WRF)
                        constexpr Real bfac_pr = 6.8;
                        const Real prnum0 = phih/phim + bfac_pr*KAPPA*sfcfrac;
                        // Prandtl number with height correction
                        constexpr Real conpr = 6.8 * KAPPA * sfcfrac; // bfac*vk*sfcfrac
                        const Real prfac2 = 15.9*ws3/ust3/(1.0 + 4.0*KAPPA*ws3/ust3);
                        const Real prnumfac = -3.0 * std::pow(amrex::max(zval - sfcfrac*pblh_arr(i,j,0), 0.0), 2.0)
                                            / (pblh_arr(i,j,0)*pblh_arr(i,j,0));
                        const Real prnum = 1.0 + (prnum0 - 1.0) * std::exp(prnumfac)
                                         + conpr * std::exp(prnumfac)
                                         + prfac2 * std::exp(prnumfac);
                        const Real prnum_lim = amrex::max(prnum, 0.25);

                        K_turb(i,j,k,EddyDiff::Mom_v) = rho * wscalek * KAPPA * zval * std::pow(zfac, pfac);
                        K_turb(i,j,k,EddyDiff::Theta_v) = K_turb(i,j,k,EddyDiff::Mom_v) / prnum_lim;
                    } else {
                        // Stable case (original code)
                        // Not including YSU top down PBL term (not in H10, added to WRF later)
                        const Real phi_term = 1 + 5 * zval / l_obuk_arr(i,j,0);
                        const Real wscalek = std::max(u_star_arr(i,j,0) / phi_term, 0.001);
                        K_turb(i,j,k,EddyDiff::Mom_v) = rho * wscalek * KAPPA * zval * std::pow(zfac, pfac);
                        K_turb(i,j,k,EddyDiff::Theta_v) = K_turb(i,j,k,EddyDiff::Mom_v);
                    }
                } else {
                    // -- Compute coefficients in free stream above PBL
                    constexpr Real lam0 = 30.0;
                    constexpr Real min_richardson = -100.0;
                    constexpr Real prandtl_max = 4.0;
                    Real dthetadz, dudz, dvdz;
                    ComputeVerticalDerivativesPBL(i, j, k,
                                                  uvel, vvel, cell_data, izmin, izmax, 1.0/dz_terrain,
                                                  c_ext_dir_on_zlo, c_ext_dir_on_zhi,
                                                  u_ext_dir_on_zlo, u_ext_dir_on_zhi,
                                                  v_ext_dir_on_zlo, v_ext_dir_on_zhi,
                                                  dthetadz, dudz, dvdz, moisture_indices);
                    const Real shear_squared = dudz*dudz + dvdz*dvdz + 1.0e-9; // 1.0e-9 from WRF to avoid divide by zero
                    const Real theta = cell_data(i,j,k,RhoTheta_comp) / cell_data(i,j,k,Rho_comp);
                    Real richardson = CONST_GRAV / theta * dthetadz / shear_squared;
                    const Real lambdadz = std::min(std::max(0.1*dz_terrain , lam0), 300.0); // in WRF, H10 paper just says use lam0
                    const Real lengthscale = lambdadz * KAPPA * zval / (lambdadz + KAPPA * zval);
                    const Real turbfact = lengthscale * lengthscale * std::sqrt(shear_squared);

                    if (richardson < 0) {
                        richardson = max(richardson, min_richardson);
                        Real sqrt_richardson = std::sqrt(-richardson);
                        K_turb(i,j,k,EddyDiff::Mom_v) = rho * turbfact * (1.0 - 8.0 * richardson / (1.0 + 1.746 * sqrt_richardson));
                        K_turb(i,j,k,EddyDiff::Theta_v) = rho * turbfact * (1.0 - 8.0 * richardson / (1.0 + 1.286 * sqrt_richardson));
                    } else {
                        const Real oneplus5ri = 1.0 + 5.0 * richardson;
                        K_turb(i,j,k,EddyDiff::Theta_v) = rho * turbfact / (oneplus5ri * oneplus5ri);
                        const Real prandtl = std::min(1.0+2.1*richardson, prandtl_max); // limit from WRF
                        K_turb(i,j,k,EddyDiff::Mom_v) = K_turb(i,j,k,EddyDiff::Theta_v) * prandtl;
                    }
                }

                // limit both diffusion coefficients - from WRF, not documented in papers
                constexpr Real ckz = 0.001;
                constexpr Real Kmax = 1000.0;
                const Real rhoKmin = ckz * dz_terrain * rho;
                const Real rhoKmax = rho * Kmax;
                K_turb(i,j,k,EddyDiff::Mom_v) = std::max(std::min(K_turb(i,j,k,EddyDiff::Mom_v) ,rhoKmax), rhoKmin);
                K_turb(i,j,k,EddyDiff::Theta_v) = std::max(std::min(K_turb(i,j,k,EddyDiff::Theta_v) ,rhoKmax), rhoKmin);
                K_turb(i,j,k,EddyDiff::Turb_lengthscale) = pblh_arr(i,j,0);
            });

            // HACK set bottom ghost cell to 1st cell
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k==-1) {
                    K_turb(i,j,k,EddyDiff::Mom_v) = K_turb(i,j,0,EddyDiff::Mom_v);
                    K_turb(i,j,k,EddyDiff::Theta_v) = K_turb(i,j,0,EddyDiff::Theta_v);
                }
            });
        }
    }
}
