#ifndef ERF_CLOUD_CHAMBER_WALL_FLUX_DIAG_H_
#define ERF_CLOUD_CHAMBER_WALL_FLUX_DIAG_H_

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>

#include <array>
#include <string>

/**
 * Interval-averaged Cloud Chamber wall fluxes (erf.wallflux_diagnostics).
 *
 * Every row of wallflux_lev0.txt holds, for each of the six domain walls,
 * the mean over one reporting window of
 *
 *   tau  the magnitude of the applied tangential wall stress [Pa],
 *   H    the applied sensible heat flux Cp_d * Exner * F(rho theta) [W/m2],
 *   E    the applied conservative water-vapor flux F(rho qv) [kg/m2/s],
 *
 * labeled with the window's END time.  H and E keep ERF's coordinate-flux
 * sign (positive toward +x, +y or +z on both the low and the high wall);
 * tau is a magnitude and is nonnegative.
 *
 * The values are the face fluxes the RHS applies.  Each RK stage contributes
 * stage_weight * sum_patches(area_patch * value_patch); an accepted step
 * multiplies its stage sum by its actual duration once; a window divides the
 * MPI-global integral by face_area * sum(accepted dt).  The stage weights
 * follow the anelastic integrators (ERF_MRI.H): RK2 (SSP-RK2) averages the
 * two stage sources (1/2, 1/2); MidPoint applies only the second stage (0, 1).
 *
 * This diagnostic is independent of CloudChamberBudget and output-only: it
 * reads the retained stress and flux arrays and never writes to them.
 */
namespace erf_wallflux_diag {

enum Channel : int { Stress = 0, Heat = 1, Vapor = 2, NumChannels = 3 };

inline constexpr int NumFaces = 2 * AMREX_SPACEDIM;
inline constexpr int NumStages = 2;
inline constexpr int NumColumns = 1 + NumChannels * NumFaces;

using FaceValues = std::array<amrex::Real, NumFaces>;
using ChannelFaceValues = std::array<FaceValues, NumChannels>;

enum class StageWeighting { Trapezoidal, MidPoint };

/** Weight of anelastic RK stage nrk in its step's applied source. */
inline amrex::Real stage_weight (StageWeighting weighting, int nrk) noexcept
{
    if (nrk < 0 || nrk >= NumStages) { return amrex::Real(0.0); }
    if (weighting == StageWeighting::MidPoint) {
        return (nrk == 1) ? amrex::Real(1.0) : amrex::Real(0.0);
    }
    return amrex::Real(0.5);
}

/**
 * Host-side window accumulator holding the interval math.
 *
 * Stage captures go into pending step integrals and are committed only when
 * the step is accepted.  Every required channel must be captured exactly once
 * in every stage of a step; capture() and commit() return a nonempty message
 * when that contract is broken and leave the accumulator unchanged.
 */
class WallFluxWindow
{
public:
    WallFluxWindow (StageWeighting weighting, bool has_vapor) noexcept
        : m_weighting(weighting), m_has_vapor(has_vapor) {}

    StageWeighting weighting () const noexcept { return m_weighting; }
    bool requires_channel (int channel) const noexcept
    {
        return channel == Stress || channel == Heat ||
               (channel == Vapor && m_has_vapor);
    }
    amrex::Real stage_weight (int nrk) const noexcept
    {
        return erf_wallflux_diag::stage_weight(m_weighting, nrk);
    }

    std::string capture (int channel, int nrk, const FaceValues& area_integral);
    std::string commit (amrex::Real dt);
    void reset () noexcept;

    const ChannelFaceValues& window_integral () const noexcept { return m_window; }
    amrex::Real elapsed () const noexcept { return m_elapsed; }
    int committed_steps () const noexcept { return m_steps; }
    bool has_pending () const noexcept;

private:
    StageWeighting m_weighting;
    bool m_has_vapor;
    ChannelFaceValues m_pending{};
    std::array<std::array<int, NumStages>, NumChannels> m_stage_count{};
    ChannelFaceValues m_window{};
    amrex::Real m_elapsed = amrex::Real(0.0);
    int m_steps = 0;
};

/** Window means: integral / (face_area * elapsed) for every channel and face. */
ChannelFaceValues interval_means (const ChannelFaceValues& integral,
                                  const FaceValues& face_area,
                                  amrex::Real elapsed);

/** The comment header line of wallflux_lev0.txt (no trailing newline). */
std::string header_line ();

/** One data row: time then FACE_tau FACE_H FACE_E for xlo xhi ylo yhi zlo zhi. */
std::string format_row (double time, const ChannelFaceValues& means);

/** Physical area of each domain wall, in output face order. */
FaceValues face_areas (const amrex::Geometry& geom);

using FaceReduceOps = amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum,
                                       amrex::ReduceOpSum, amrex::ReduceOpSum,
                                       amrex::ReduceOpSum, amrex::ReduceOpSum>;
using FaceReduceData = amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real,
                                         amrex::Real, amrex::Real, amrex::Real>;

/**
 * Rank-local six-face sum filled from inside an MFIter tile loop.  Construct
 * it outside any OpenMP parallel region; AMReX then keeps one partial sum per
 * thread (CPU) or accumulates on the device (GPU).
 */
class FaceReducer
{
public:
    FaceReducer () : m_data(m_ops) {}
    FaceReducer (const FaceReducer&) = delete;
    FaceReducer& operator= (const FaceReducer&) = delete;
    FaceReducer (FaceReducer&&) = delete;
    FaceReducer& operator= (FaceReducer&&) = delete;
    ~FaceReducer () = default;

    FaceReduceOps& ops () noexcept { return m_ops; }
    FaceReduceData& data () noexcept { return m_data; }
    FaceValues value ();

private:
    FaceReduceOps m_ops;
    FaceReduceData m_data;
};

/**
 * Add area * |tau_t| for the wall patches owned by the cell tile @p tile.
 *
 * The stress arrays are ERF's retained cross stresses after the resolved
 * closure and any Cloud Chamber traction override (DiffusionSrcForMom
 * convention: low walls store the traction on the fluid, high walls its
 * negative).  tau12 is nodal in x and y, tau13 in x and z, tau23 in y and z.
 * Each wall patch is owned by its adjacent valid cell.  Both tangential
 * components are reconstructed at the patch center from the two wall-edge
 * nodes that bracket it, and their magnitude is taken before any reduction.
 * A node on a domain edge, shared with a perpendicular wall, is excluded: the
 * normal velocities on both walls are fixed there, so no RHS applies it, and
 * the traction override writes it last-wall-wins.  A patch next to such an
 * edge takes its component from its single applied node.
 */
void add_wall_stress (FaceReducer& reducer,
                      const amrex::Box& tile,
                      const amrex::Box& domain,
                      const amrex::Array4<const amrex::Real>& tau12,
                      const amrex::Array4<const amrex::Real>& tau13,
                      const amrex::Array4<const amrex::Real>& tau23,
                      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx);

/**
 * Rank-local area integrals of a retained scalar face flux on the six walls.
 *
 * @p cell_layout is any cell-centered MultiFab on the flux arrays' BoxArray
 * and DistributionMapping; it decides which valid cell owns each wall patch.
 * With @p base_state the coordinate flux of rho*theta is converted locally to
 * Cp_d * (p_hse/p_0)^rdOcp * F at the wall-adjacent cell's hydrostatic
 * pressure, the pressure the Cloud Chamber wall thermodynamics uses; without
 * it the flux is integrated as stored.
 */
FaceValues scalar_face_integrals (const amrex::MultiFab& cell_layout,
                                  const amrex::MultiFab& xflux,
                                  const amrex::MultiFab& yflux,
                                  const amrex::MultiFab& zflux,
                                  int flux_comp,
                                  const amrex::MultiFab* base_state,
                                  amrex::Real rdOcp,
                                  const amrex::Geometry& geom);

} // namespace erf_wallflux_diag

/**
 * ERF-owned level-0 interval-mean wall-flux writer (wallflux_lev0.txt).
 * Allocated only when erf.wallflux_diagnostics is true.
 */
class CloudChamberWallFluxDiag
{
public:
    CloudChamberWallFluxDiag (int interval,
                              erf_wallflux_diag::StageWeighting weighting,
                              bool has_vapor,
                              const amrex::Geometry& geom,
                              std::string filename = "wallflux_lev0.txt");

    int interval () const noexcept { return m_interval; }
    const std::string& filename () const noexcept { return m_filename; }
    amrex::Real stage_weight (int nrk) const noexcept { return m_window.stage_weight(nrk); }
    bool requires_channel (int channel) const noexcept { return m_window.requires_channel(channel); }

    /** Add one stage's rank-local area integrals; aborts on a duplicate capture. */
    void capture_stage (int channel, int nrk,
                        const erf_wallflux_diag::FaceValues& area_integral);

    /** Commit the accepted step of duration dt; aborts on a missing capture. */
    void commit_step (amrex::Real dt);

    bool due (int step) const noexcept
    {
        return step > 0 && (step % m_interval) == 0;
    }

    /**
     * When @p step closes a window: reduce, write one row at @p time (the
     * window end) and reset.  Collective; call on every rank.
     */
    void report (int step, double time);

    /** Abort if the output file already exists (fresh runs never append). */
    void require_fresh_output () const;

private:
    bool write_row (double time, const erf_wallflux_diag::ChannelFaceValues& means);

    int m_interval;
    std::string m_filename;
    erf_wallflux_diag::FaceValues m_face_area{};
    erf_wallflux_diag::WallFluxWindow m_window;
    bool m_header_written = false;
};

#endif
