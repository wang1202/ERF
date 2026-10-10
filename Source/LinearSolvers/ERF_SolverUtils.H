/**
 * \file ERF_SolverUtils.H
 */
#ifndef ERF_SOLVER_UTILS_H_
#define ERF_SOLVER_UTILS_H_

#include "AMReX_BoxArray.H"
#include "AMReX_Geometry.H"
#include "AMReX_Gpu.H"
#include "AMReX_LO_BCTYPES.H"
#include "AMReX_ParmParse.H"
#ifdef ERF_USE_FFT
#include <AMReX_FFT.H>
#endif

#include <limits>
#include <sstream>
#include <string>

using namespace amrex;

/**
 * Opt-in synchronized timing for terrain Poisson phase comparisons.  The
 * normal production path does not introduce GPU synchronization or output.
 */
inline bool terrain_poisson_timing_enabled ()
{
    static const bool enabled = [] () {
        bool value = false;
        ParmParse pp("erf");
        pp.query("terrain_poisson_timing", value);
        return value;
    } ();
    return enabled;
}

inline void terrain_poisson_timing_sync ()
{
#ifdef AMREX_USE_GPU
    Gpu::streamSynchronize();
#endif
}

/**
 * Default solver for the terrain Poisson equation. Keep the established
 * GMRES+FFT implementation as the default until the MLMG path is validated
 * for broader ERF configurations.
 */
inline constexpr char const* terrain_poisson_solver_default () noexcept
{
    return "gmres_fft";
}

inline constexpr Real terrain_mlmg_compatibility_tolerance () noexcept
{
    return sizeof(Real) == sizeof(float) ? Real(1.e-4) : Real(1.e-6);
}

inline bool terrain_mlmg_compatibility_mean_exceeds (Real mean) noexcept
{
    return !amrex::Math::isfinite(mean) ||
           amrex::Math::abs(mean) > terrain_mlmg_compatibility_tolerance();
}

/**
 * The terrain MLMG operator has not been verified with horizontal map
 * factors. Keep that configuration from running silently until its
 * coefficient/flux identity has independent numerical coverage.
 */
inline std::string
terrain_mlmg_map_factor_error (Real max_abs_deviation)
{
    const Real roundoff_tolerance = Real(64) * std::numeric_limits<Real>::epsilon();
    if (amrex::Math::isfinite(max_abs_deviation) &&
        max_abs_deviation <= roundoff_tolerance) {
        return {};
    }
    std::ostringstream os;
    os << "erf.terrain_poisson_solver=mlmg currently supports only unit horizontal map factors; "
       << "max_abs_deviation_from_one=" << max_abs_deviation
       << ", roundoff_tolerance=" << roundoff_tolerance;
    return os.str();
}

/**
 * Detect an invalid value in a terrain-projection momentum ghost.
 *
 * Fill buffers are initialized with ERF's bogus_large_value, whose magnitude
 * differs substantially between SINGLE and DOUBLE builds. Comparing against
 * a fixed 1.e100 threshold accepts the SINGLE sentinel. Reject values at or
 * above half the actual fill sentinel while allowing physically valid values.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool terrain_projection_ghost_value_invalid (Real value, Real fill_sentinel) noexcept
{
    return !amrex::Math::isfinite(value) ||
           amrex::Math::abs(value) >= Real(0.5) * amrex::Math::abs(fill_sentinel);
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
 * Restrict the experimental per-region MLMG operator to the integration
 * scope that has been exercised: ERF levels 0 and 1 and one connected region
 * on each active level.
 */
inline std::string
terrain_mlmg_scope_error (int max_level, int lev, int num_regions)
{
    if (max_level > 1) {
        return "erf.terrain_poisson_solver=mlmg currently supports at most two ERF levels "
               "(levels 0 and 1); configured max_level=" + std::to_string(max_level);
    }
    if (num_regions > 1) {
        return "erf.terrain_poisson_solver=mlmg currently requires one connected rectangular "
               "solve region per active level; level " + std::to_string(lev) +
               " has " + std::to_string(num_regions) + " regions";
    }
    return {};
}

/**
 * Validate runtime conditions required by the experimental terrain MLMG
 * projection. Each ERF level/subdomain is represented by a separate one-level
 * operator.
 */
inline std::string
terrain_mlmg_configuration_error (bool is_3d, bool variable_dz,
                                  bool static_fitted_terrain, bool use_real_bcs)
{
    if (!is_3d) {
        return "erf.terrain_poisson_solver=mlmg requires a 3D ERF build (AMREX_SPACEDIM=3)";
    }
    if (!variable_dz) {
        return "erf.terrain_poisson_solver=mlmg requires a VariableDz mesh derived from "
               "erf.terrain_type=StaticFittedMesh; flat ConstantDz/StretchedDz meshes are unsupported";
    }
    if (!static_fitted_terrain) {
        return "erf.terrain_poisson_solver=mlmg requires erf.terrain_type=StaticFittedMesh; "
               "moving terrain, EB, and immersed terrain are unsupported";
    }
    if (use_real_bcs) {
        return "erf.terrain_poisson_solver=mlmg does not support erf.use_real_bcs=true; "
               "use homogeneous projection boundary conditions";
    }
    return {};
}

inline std::string
terrain_mlmg_box_string (Box const& box)
{
    std::ostringstream os;
    os << '[';
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (d != 0) { os << ','; }
        os << box.smallEnd(d) << ':' << box.bigEnd(d);
    }
    os << ']';
    return os.str();
}

/**
 * Require valid solve boxes to be an exact, disjoint rectangular cover of a
 * regional operator domain. This rejects holes and L-shaped connected unions.
 */
inline std::string
terrain_mlmg_region_error (int lev, int isub, Box const& level_domain,
                           Box const& region, BoxArray const& region_boxes)
{
    bool boxes_inside_region = true;
    for (int i = 0; i < region_boxes.size(); ++i) {
        boxes_inside_region = boxes_inside_region && region.contains(region_boxes[i]);
    }

    const Long actual_cells = region_boxes.numPts();
    const Long expected_cells = region.numPts();
    const bool exact_rectangle = region.ok() && region.cellCentered() &&
        level_domain.contains(region) && region_boxes.size() > 0 &&
        region_boxes.isDisjoint() && boxes_inside_region &&
        region_boxes.minimalBox() == region && actual_cells == expected_cells;

    if (exact_rectangle) { return {}; }

    std::ostringstream os;
    os << "erf.terrain_poisson_solver=mlmg cannot solve non-rectangular or invalid "
       << "ERF region: lev=" << lev << ", isub=" << isub
       << ", region=" << terrain_mlmg_box_string(region)
       << ", actual_cells=" << actual_cells
       << ", expected_cells=" << expected_cells
       << ". Connected irregular/L-shaped regions are unsupported; use one rectangular "
       << "refinement region";
    return os.str();
}

/**
 * A partial region touching a periodic physical seam cannot be represented by
 * a one-domain operator with ordinary periodic or Neumann conditions.
 */
inline std::string
terrain_mlmg_periodic_region_error (int lev, int isub,
                                    Geometry const& level_geom, Box const& region)
{
    const Box& domain = level_geom.Domain();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (!level_geom.isPeriodic(d)) { continue; }
        const bool touches_lo = region.smallEnd(d) == domain.smallEnd(d);
        const bool touches_hi = region.bigEnd(d) == domain.bigEnd(d);
        const bool spans_periodic_axis = touches_lo && touches_hi;
        if ((touches_lo || touches_hi) && !spans_periodic_axis) {
            std::ostringstream os;
            os << "erf.terrain_poisson_solver=mlmg does not support partial periodic-seam "
               << "refinement: lev=" << lev << ", isub=" << isub
               << ", dir=" << d << ", region=" << terrain_mlmg_box_string(region)
               << ". Extend the rectangular region across the complete periodic direction";
            return os.str();
        }
    }
    return {};
}

/** Construct a one-level Geometry preserving the level's index/physical mapping. */
inline Geometry
terrain_mlmg_region_geometry (Geometry const& level_geom, Box const& region)
{
    const Box& global_domain = level_geom.Domain();
    Array<Real,AMREX_SPACEDIM> prob_lo{};
    Array<Real,AMREX_SPACEDIM> prob_hi{};
    Array<int,AMREX_SPACEDIM> periodic{};
    const auto dx = level_geom.CellSizeArray();

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        prob_lo[d] = level_geom.ProbLo(d) +
            (region.smallEnd(d) - global_domain.smallEnd(d)) * dx[d];
        prob_hi[d] = level_geom.ProbLo(d) +
            (region.bigEnd(d) + 1 - global_domain.smallEnd(d)) * dx[d];
        periodic[d] = level_geom.isPeriodic(d) &&
            region.smallEnd(d) == global_domain.smallEnd(d) &&
            region.bigEnd(d) == global_domain.bigEnd(d);
    }

    const RealBox physical_region(prob_lo, prob_hi);
    return Geometry(region, physical_region, level_geom.Coord(), periodic);
}

struct TerrainMLMGRegionBCs
{
    Array<LinOpBCType,AMREX_SPACEDIM> lo{};
    Array<LinOpBCType,AMREX_SPACEDIM> hi{};
    bool is_singular = true;
};

/**
 * Build effective regional BCs: physical faces retain ERF's projection map,
 * artificial subdomain faces are homogeneous Neumann, and terrain zlo is
 * always homogeneous Neumann. Periodicity is retained only for a full-axis
 * regional domain.
 */
inline TerrainMLMGRegionBCs
terrain_mlmg_region_bcs (Geometry const& level_geom, Box const& region,
                         Array<std::string,2*AMREX_SPACEDIM> const& domain_bc_names)
{
    TerrainMLMGRegionBCs result;
    const Box& domain = level_geom.Domain();
    auto physical_bc = [&] (int dir, auto side) {
        if (dir == 2 && side == Orientation::low) {
            return LinOpBCType::Neumann;
        }
        const auto& name = domain_bc_names[Orientation(dir, side)];
        return (name == "Outflow" || name == "Open")
            ? LinOpBCType::Dirichlet : LinOpBCType::Neumann;
    };

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const bool at_lo = region.smallEnd(d) == domain.smallEnd(d);
        const bool at_hi = region.bigEnd(d) == domain.bigEnd(d);
        if (level_geom.isPeriodic(d) && at_lo && at_hi) {
            result.lo[d] = LinOpBCType::Periodic;
            result.hi[d] = LinOpBCType::Periodic;
        } else {
            result.lo[d] = at_lo ? physical_bc(d, Orientation::low)
                                 : LinOpBCType::Neumann;
            result.hi[d] = at_hi ? physical_bc(d, Orientation::high)
                                 : LinOpBCType::Neumann;
        }
        if (result.lo[d] == LinOpBCType::Dirichlet ||
            result.hi[d] == LinOpBCType::Dirichlet) {
            result.is_singular = false;
        }
    }
    return result;
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
