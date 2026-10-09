/**
 * \file ERF_SolverUtils.H
 */
#ifndef ERF_SOLVER_UTILS_H_
#define ERF_SOLVER_UTILS_H_

#include "AMReX_LO_BCTYPES.H"

#include <string>

using namespace amrex;

/**
 * Default solver for the terrain Poisson equation. Keep the established
 * GMRES+FFT implementation as the default until the MLMG path is validated
 * for broader ERF configurations.
 */
inline constexpr char const* terrain_poisson_solver_default () noexcept
{
    return "gmres_fft";
}

/**
 * Return an actionable error for an unsupported terrain Poisson solver name.
 */
inline std::string
terrain_poisson_solver_choice_error (std::string const& choice)
{
    if (choice != "gmres_fft" && choice != "mlmg") {
        return "erf.terrain_poisson_solver must be gmres_fft or mlmg";
    }
    return {};
}

/**
 * Validate runtime conditions required by the experimental terrain MLMG
 * projection. Boolean arguments keep this check independent of ERF's solver
 * configuration types and make its diagnostics directly unit-testable.
 */
inline std::string
terrain_mlmg_configuration_error (int lev, bool is_3d, bool variable_dz,
                                  bool static_fitted_terrain, bool use_real_bcs,
                                  bool full_domain)
{
    if (!is_3d) {
        return "erf.terrain_poisson_solver=mlmg requires a 3D ERF build (AMREX_SPACEDIM=3)";
    }
    if (!variable_dz) {
        return "erf.terrain_poisson_solver=mlmg requires erf.mesh_type=VariableDz; "
               "flat ConstantDz/StretchedDz meshes are unsupported";
    }
    if (!static_fitted_terrain) {
        return "erf.terrain_poisson_solver=mlmg requires erf.terrain_type=StaticFittedMesh; "
               "moving terrain, EB, and immersed terrain are unsupported";
    }
    if (use_real_bcs) {
        return "erf.terrain_poisson_solver=mlmg does not support erf.use_real_bcs=true; "
               "use homogeneous projection boundary conditions";
    }
    if (lev != 0) {
        return "erf.terrain_poisson_solver=mlmg supports level 0 only; set amr.max_level=0 "
               "or select erf.terrain_poisson_solver=gmres_fft";
    }
    if (!full_domain) {
        return "erf.terrain_poisson_solver=mlmg requires exactly one ERF subdomain covering "
               "the complete level-0 domain; use a full-domain rectangular layout";
    }
    return {};
}

/**
 * Define the domain boundary conditions for the (optional) Poisson solve
 * if we want to enforce that the initial conditions satisfy the constraint
 *
 * @param lev_geom Geometry used to identify periodic low-side boundaries
 * @param l_domain_bc_type Domain boundary-condition names
 * @return Low-side linear-operator boundary condition in each direction
 */
inline Array<LinOpBCType,AMREX_SPACEDIM>
get_lo_projection_bc (Geometry const& lev_geom,
                      Array<std::string,2*AMREX_SPACEDIM> l_domain_bc_type)
{
    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if ( lev_geom.isPeriodic(dir) )
        {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc_type = l_domain_bc_type[Orientation(dir,Orientation::low)];
            if (bc_type == "Outflow" || bc_type == "Open") {
                r[dir] = LinOpBCType::Dirichlet;
            } else {
                r[dir] = LinOpBCType::Neumann;
            }
        }
    }
    return r;
}

/**
 * Define low-side projection BCs for the terrain MLMG operator. ERF's legacy
 * terrain projection imposes zero correction flux at zlo even if the physical
 * velocity boundary is outflow, so preserve that terrain-specific convention.
 */
inline Array<LinOpBCType,AMREX_SPACEDIM>
get_lo_terrain_projection_bc (Geometry const& lev_geom,
                              Array<std::string,2*AMREX_SPACEDIM> l_domain_bc_type)
{
    auto r = get_lo_projection_bc(lev_geom, l_domain_bc_type);
#if (AMREX_SPACEDIM == 3)
    if (!lev_geom.isPeriodic(2)) {
        r[2] = LinOpBCType::Neumann;
    }
#endif
    return r;
}

/**
 * Define high-side domain boundary conditions for the projection solve.
 *
 * @param lev_geom Geometry used to identify periodic high-side boundaries
 * @param l_domain_bc_type Domain boundary-condition names
 * @return High-side linear-operator boundary condition in each direction
 */
inline Array<LinOpBCType,AMREX_SPACEDIM>
get_hi_projection_bc (Geometry const& lev_geom,
                      Array<std::string,2*AMREX_SPACEDIM> l_domain_bc_type)
{
    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if ( lev_geom.isPeriodic(dir) )
        {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc_type = l_domain_bc_type[Orientation(dir,Orientation::high)];
            if (bc_type == "Outflow" || bc_type == "Open") {
                r[dir] = LinOpBCType::Dirichlet;
            } else {
                r[dir] = LinOpBCType::Neumann;
            }
        }
    }
    return r;
}

#ifdef ERF_USE_FFT

/**
 * Define FFT boundary conditions on the specified solve box.
 *
 * @param lev_geom Geometry for the full level domain
 * @param l_domain_bc_type Domain boundary-condition names
 * @param bounding_box Box covered by the FFT solve
 * @param use_real_bcs Whether real lateral boundary conditions are active
 * @return Low/high FFT boundary-condition pair in each direction
 */
inline Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM>
get_fft_bc (Geometry const& lev_geom,
            Array<std::string,2*AMREX_SPACEDIM> l_domain_bc_type,
            Box const& bounding_box, bool use_real_bcs) noexcept
{
    Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM> r;

    FFT::Boundary first, second;

    for (int dir = 0; dir <= 2; dir++)
    {
        auto bc_type_lo = l_domain_bc_type[Orientation(dir,Orientation::low)];
        if ( lev_geom.isPeriodic(dir) &&
            (lev_geom.Domain().smallEnd(dir) == bounding_box.smallEnd(dir)) )
        {
            first = FFT::Boundary::periodic;
             // amrex::Print() << "SETTING FIRST OF " << dir << " TO PERIODIC " << std::endl;
        }

        else if ( (lev_geom.Domain().smallEnd(dir) == bounding_box.smallEnd(dir)) &&
                  (bc_type_lo == "Outflow" || bc_type_lo == "Open") && (dir == 2 || !use_real_bcs) )
        {
            first = FFT::Boundary::odd;
             // amrex::Print() << "SETTING FIRST OF " << dir << " TO ODD " << std::endl;
        }

        else
        {
            first = FFT::Boundary::even;
            // amrex::Print() << "SETTING FIRST OF " << dir << " TO EVEN " << std::endl;
        }

        auto bc_type_hi = l_domain_bc_type[Orientation(dir,Orientation::high)];
        if ( lev_geom.isPeriodic(dir) &&
            (lev_geom.Domain().bigEnd(dir) == bounding_box.bigEnd(dir)) )
        {
            second = FFT::Boundary::periodic;
             // amrex::Print() << "SETTING SECOND OF " << dir << " TO PERIODIC " << std::endl;
        }

        else if ( (lev_geom.Domain().bigEnd(dir) == bounding_box.bigEnd(dir)) &&
                  (bc_type_hi == "Outflow" || bc_type_hi == "Open") && (dir == 2 || !use_real_bcs) )
        {
            second = FFT::Boundary::odd;
             // amrex::Print() << "SETTING SECOND OF " << dir << " TO ODD " << std::endl;
        }

        else
        {
            second = FFT::Boundary::even;
            // amrex::Print() << "SETTING SECOND OF " << dir << " TO EVEN " << std::endl;
        }

        r[dir] = std::make_pair(first,second);

    } // dir

    return r;
}
#endif
#endif
