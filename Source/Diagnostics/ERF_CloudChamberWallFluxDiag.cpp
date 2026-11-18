#include "Diagnostics/ERF_CloudChamberWallFluxDiag.H"

#include "ERF_Constants.H"
#include "ERF_IndexDefines.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <utility>

using namespace amrex;

namespace erf_wallflux_diag {

namespace {

const char* const face_names[NumFaces] = {"xlo", "xhi", "ylo", "yhi", "zlo", "zhi"};
const char* const channel_names[NumChannels] = {"tau", "H", "E"};

using FaceTuple = FaceReduceData::Type;

template <int F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
FaceTuple face_tuple (Real value) noexcept
{
    return {F == 0 ? value : Real(0.0), F == 1 ? value : Real(0.0),
            F == 2 ? value : Real(0.0), F == 3 ? value : Real(0.0),
            F == 4 ? value : Real(0.0), F == 5 ? value : Real(0.0)};
}

// One tangential stress component at a wall-patch center, from the two
// wall-edge nodes that bracket the patch.  A node on a domain edge is not
// applied by any RHS (see add_wall_stress) and is left out.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real patch_center_component (const Array4<const Real>& tau,
                             int i0, int j0, int k0,
                             int i1, int j1, int k1,
                             bool lo_applied, bool hi_applied) noexcept
{
    if (lo_applied && hi_applied) {
        return Real(0.5) * (tau(i0,j0,k0) + tau(i1,j1,k1));
    }
    if (lo_applied) { return tau(i0,j0,k0); }
    if (hi_applied) { return tau(i1,j1,k1); }
    return Real(0.0);
}

// The wall-adjacent cell layer of a box, or an empty box if the box does not
// reach the wall.
template <int DIR, bool HIGH>
Box wall_cells (const Box& bx, const Box& domain)
{
    const int cell = HIGH ? domain.bigEnd(DIR) : domain.smallEnd(DIR);
    Box cells = bx;
    if (bx.smallEnd(DIR) > cell || bx.bigEnd(DIR) < cell) {
        cells.setSmall(DIR, cell + 1);
        cells.setBig(DIR, cell);
        return cells;
    }
    cells.setSmall(DIR, cell);
    cells.setBig(DIR, cell);
    return cells;
}

// tau_a and tau_b hold the cross stresses whose tangent directions are TA and
// TB: tau12/tau13 on x walls, tau12/tau23 on y walls, tau13/tau23 on z walls.
template <int DIR, bool HIGH>
void add_stress_face (FaceReducer& reducer,
                      const Box& tile,
                      const Box& domain,
                      const Array4<const Real>& tau_a,
                      const Array4<const Real>& tau_b,
                      Real area)
{
    constexpr int TA = (DIR == 0) ? 1 : 0;
    constexpr int TB = (DIR == 2) ? 1 : 2;
    constexpr int F = 2 * DIR + (HIGH ? 1 : 0);

    const Box cells = wall_cells<DIR,HIGH>(tile, domain);
    if (cells.isEmpty()) { return; }

    const int node = HIGH ? domain.bigEnd(DIR) + 1 : domain.smallEnd(DIR);
    const int ta_lo = domain.smallEnd(TA);
    const int ta_hi = domain.bigEnd(TA);
    const int tb_lo = domain.smallEnd(TB);
    const int tb_hi = domain.bigEnd(TB);
    // DiffusionSrcForMom stores the wall-on-fluid traction on a low wall and
    // its negative on a high wall (ERF_CloudChamberWallStress.H).
    const Real to_traction = HIGH ? Real(-1.0) : Real(1.0);

    reducer.ops().eval(cells, reducer.data(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> FaceTuple
    {
        int n0[3] = {i, j, k};
        n0[DIR] = node;

        int a1[3] = {n0[0], n0[1], n0[2]};
        a1[TA] += 1;
        const Real comp_a = to_traction * patch_center_component(
            tau_a, n0[0], n0[1], n0[2], a1[0], a1[1], a1[2],
            n0[TA] != ta_lo, n0[TA] != ta_hi);

        int b1[3] = {n0[0], n0[1], n0[2]};
        b1[TB] += 1;
        const Real comp_b = to_traction * patch_center_component(
            tau_b, n0[0], n0[1], n0[2], b1[0], b1[1], b1[2],
            n0[TB] != tb_lo, n0[TB] != tb_hi);

        return face_tuple<F>(area * std::sqrt(comp_a * comp_a + comp_b * comp_b));
    });
}

template <int DIR, bool HIGH>
void add_scalar_face (FaceReducer& reducer,
                      const Box& valid,
                      const Box& domain,
                      const Array4<const Real>& flux,
                      int flux_comp,
                      const Array4<const Real>& base_state,
                      bool to_heat_flux,
                      Real rdOcp,
                      Real area)
{
    constexpr int F = 2 * DIR + (HIGH ? 1 : 0);

    const Box cells = wall_cells<DIR,HIGH>(valid, domain);
    if (cells.isEmpty()) { return; }

    reducer.ops().eval(cells, reducer.data(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> FaceTuple
    {
        // The wall face of cell (i,j,k), as in erf_cloud_chamber_wall_flux::apply_face
        const int fi = (HIGH && DIR == 0) ? i + 1 : i;
        const int fj = (HIGH && DIR == 1) ? j + 1 : j;
        const int fk = (HIGH && DIR == 2) ? k + 1 : k;
        Real value = flux(fi,fj,fk,flux_comp);
        if (to_heat_flux) {
            // Exner function at the hydrostatic pressure of the wall-adjacent
            // cell, the pressure the wall thermodynamics evaluates theta_wall at
            const Real exner = std::pow(base_state(i,j,k,BaseState::p0_comp) / p_0, rdOcp);
            value *= Cp_d * exner;
        }
        return face_tuple<F>(area * value);
    });
}

} // namespace

std::string
WallFluxWindow::capture (int channel, int nrk, const FaceValues& area_integral)
{
    if (channel < 0 || channel >= NumChannels) {
        return "invalid channel " + std::to_string(channel);
    }
    if (nrk < 0 || nrk >= NumStages) {
        return std::string("channel ") + channel_names[channel] +
               " captured at unsupported RK stage " + std::to_string(nrk);
    }
    if (!requires_channel(channel)) {
        return std::string("channel ") + channel_names[channel] +
               " captured although this configuration does not transport it";
    }
    if (m_stage_count[channel][nrk] != 0) {
        return std::string("duplicate capture of channel ") + channel_names[channel] +
               " at RK stage " + std::to_string(nrk);
    }
    m_stage_count[channel][nrk] = 1;
    const Real weight = stage_weight(nrk);
    for (int f = 0; f < NumFaces; ++f) {
        m_pending[channel][f] += weight * area_integral[f];
    }
    return {};
}

std::string
WallFluxWindow::commit (Real dt)
{
    if (!(dt > Real(0.0)) || !std::isfinite(static_cast<double>(dt))) {
        return "nonpositive or nonfinite step duration " + std::to_string(dt);
    }
    for (int c = 0; c < NumChannels; ++c) {
        if (!requires_channel(c)) { continue; }
        for (int s = 0; s < NumStages; ++s) {
            if (m_stage_count[c][s] != 1) {
                return std::string("missing capture of channel ") + channel_names[c] +
                       " at RK stage " + std::to_string(s);
            }
        }
    }
    for (int c = 0; c < NumChannels; ++c) {
        for (int f = 0; f < NumFaces; ++f) {
            m_window[c][f] += dt * m_pending[c][f];
        }
    }
    m_elapsed += dt;
    ++m_steps;
    m_pending = {};
    m_stage_count = {};
    return {};
}

void
WallFluxWindow::reset () noexcept
{
    m_window = {};
    m_elapsed = Real(0.0);
    m_steps = 0;
}

bool
WallFluxWindow::has_pending () const noexcept
{
    for (const auto& channel : m_stage_count) {
        for (const int count : channel) {
            if (count != 0) { return true; }
        }
    }
    return false;
}

ChannelFaceValues
interval_means (const ChannelFaceValues& integral,
                const FaceValues& face_area,
                Real elapsed)
{
    ChannelFaceValues means{};
    for (int c = 0; c < NumChannels; ++c) {
        for (int f = 0; f < NumFaces; ++f) {
            means[c][f] = integral[c][f] / (face_area[f] * elapsed);
        }
    }
    return means;
}

std::string
header_line ()
{
    // The column layout of wang1202/ERF_01272026 ERF_WriteWallFluxes.cpp,
    // with the comment marker at the start of the line.
    std::ostringstream out;
    out << std::left << std::setw(15) << "# time[s]" << std::right;
    for (int f = 0; f < NumFaces; ++f) {
        out << std::setw(18) << (std::string(face_names[f]) + "_tau[Pa]")
            << std::setw(18) << (std::string(face_names[f]) + "_H[W/m2]")
            << std::setw(18) << (std::string(face_names[f]) + "_E[kg/m2/s]");
    }
    return out.str();
}

std::string
format_row (double time, const ChannelFaceValues& means)
{
    std::ostringstream out;
    out << std::scientific << std::setprecision(8) << std::setw(15) << time;
    for (int f = 0; f < NumFaces; ++f) {
        for (int c = 0; c < NumChannels; ++c) {
            out << std::setw(18) << static_cast<double>(means[c][f]);
        }
    }
    return out.str();
}

FaceValues
face_areas (const Geometry& geom)
{
    const Real lx = geom.ProbLength(0);
    const Real ly = geom.ProbLength(1);
    const Real lz = geom.ProbLength(2);
    return {ly * lz, ly * lz, lx * lz, lx * lz, lx * ly, lx * ly};
}

FaceValues
FaceReducer::value ()
{
    const auto sums = m_data.value(m_ops);
    return {amrex::get<0>(sums), amrex::get<1>(sums), amrex::get<2>(sums),
            amrex::get<3>(sums), amrex::get<4>(sums), amrex::get<5>(sums)};
}

void
add_wall_stress (FaceReducer& reducer,
                 const Box& tile,
                 const Box& domain,
                 const Array4<const Real>& tau12,
                 const Array4<const Real>& tau13,
                 const Array4<const Real>& tau23,
                 const GpuArray<Real, AMREX_SPACEDIM>& dx)
{
    const Real area_x = dx[1] * dx[2];
    const Real area_y = dx[0] * dx[2];
    const Real area_z = dx[0] * dx[1];
    add_stress_face<0,false>(reducer, tile, domain, tau12, tau13, area_x);
    add_stress_face<0,true >(reducer, tile, domain, tau12, tau13, area_x);
    add_stress_face<1,false>(reducer, tile, domain, tau12, tau23, area_y);
    add_stress_face<1,true >(reducer, tile, domain, tau12, tau23, area_y);
    add_stress_face<2,false>(reducer, tile, domain, tau13, tau23, area_z);
    add_stress_face<2,true >(reducer, tile, domain, tau13, tau23, area_z);
}

FaceValues
scalar_face_integrals (const MultiFab& cell_layout,
                       const MultiFab& xflux,
                       const MultiFab& yflux,
                       const MultiFab& zflux,
                       int flux_comp,
                       const MultiFab* base_state,
                       Real rdOcp,
                       const Geometry& geom)
{
    AMREX_ALWAYS_ASSERT(cell_layout.ixType().cellCentered());
    AMREX_ALWAYS_ASSERT(flux_comp >= 0 && flux_comp < xflux.nComp() &&
                        flux_comp < yflux.nComp() && flux_comp < zflux.nComp());
    AMREX_ALWAYS_ASSERT(base_state == nullptr || base_state->nComp() > BaseState::p0_comp);

    const Box& domain = geom.Domain();
    const Real* dx = geom.CellSize();
    const Real area_x = dx[1] * dx[2];
    const Real area_y = dx[0] * dx[2];
    const Real area_z = dx[0] * dx[1];
    const bool to_heat = (base_state != nullptr);

    FaceReducer reducer;
    for (MFIter mfi(cell_layout); mfi.isValid(); ++mfi) {
        const Box& valid = mfi.validbox();
        const Array4<const Real> fx = xflux.const_array(mfi);
        const Array4<const Real> fy = yflux.const_array(mfi);
        const Array4<const Real> fz = zflux.const_array(mfi);
        const Array4<const Real> p0 = to_heat ? base_state->const_array(mfi) : Array4<const Real>{};
        add_scalar_face<0,false>(reducer, valid, domain, fx, flux_comp, p0, to_heat, rdOcp, area_x);
        add_scalar_face<0,true >(reducer, valid, domain, fx, flux_comp, p0, to_heat, rdOcp, area_x);
        add_scalar_face<1,false>(reducer, valid, domain, fy, flux_comp, p0, to_heat, rdOcp, area_y);
        add_scalar_face<1,true >(reducer, valid, domain, fy, flux_comp, p0, to_heat, rdOcp, area_y);
        add_scalar_face<2,false>(reducer, valid, domain, fz, flux_comp, p0, to_heat, rdOcp, area_z);
        add_scalar_face<2,true >(reducer, valid, domain, fz, flux_comp, p0, to_heat, rdOcp, area_z);
    }
    return reducer.value();
}

} // namespace erf_wallflux_diag

CloudChamberWallFluxDiag::CloudChamberWallFluxDiag (int interval,
                                                    erf_wallflux_diag::StageWeighting weighting,
                                                    bool has_vapor,
                                                    const Geometry& geom,
                                                    std::string filename)
    : m_interval(interval),
      m_filename(std::move(filename)),
      m_face_area(erf_wallflux_diag::face_areas(geom)),
      m_window(weighting, has_vapor)
{
    AMREX_ALWAYS_ASSERT(m_interval > 0);
}

void
CloudChamberWallFluxDiag::capture_stage (int channel, int nrk,
                                         const erf_wallflux_diag::FaceValues& area_integral)
{
    const std::string error = m_window.capture(channel, nrk, area_integral);
    if (!error.empty()) {
        Abort("erf.wallflux_diagnostics: " + error);
    }
}

void
CloudChamberWallFluxDiag::commit_step (Real dt)
{
    const std::string error = m_window.commit(dt);
    if (!error.empty()) {
        Abort("erf.wallflux_diagnostics: cannot commit the accepted step: " + error);
    }
}

void
CloudChamberWallFluxDiag::report (int step, double time)
{
    using namespace erf_wallflux_diag;
    if (!due(step)) { return; }
    AMREX_ALWAYS_ASSERT(!m_window.has_pending());
    AMREX_ALWAYS_ASSERT(m_window.committed_steps() > 0 && m_window.elapsed() > Real(0.0));

    // One packed reduction per report; the elapsed time is the same on every rank.
    std::array<Real, NumChannels * NumFaces> packed{};
    const ChannelFaceValues& local = m_window.window_integral();
    for (int c = 0; c < NumChannels; ++c) {
        for (int f = 0; f < NumFaces; ++f) {
            packed[c * NumFaces + f] = local[c][f];
        }
    }
    ParallelDescriptor::ReduceRealSum(packed.data(), static_cast<int>(packed.size()));
    ChannelFaceValues global{};
    for (int c = 0; c < NumChannels; ++c) {
        for (int f = 0; f < NumFaces; ++f) {
            global[c][f] = packed[c * NumFaces + f];
        }
    }
    const ChannelFaceValues means = interval_means(global, m_face_area, m_window.elapsed());

    int written = 1;
    if (ParallelDescriptor::IOProcessor()) {
        written = write_row(time, means) ? 1 : 0;
    }
    ParallelDescriptor::Bcast(&written, 1, ParallelDescriptor::IOProcessorNumber());
    if (written != 0) {
        m_window.reset();
    } else {
        // Keep the window: the next row then covers it, and its duration is
        // still the difference between its end time and the last row's.
        Warning("erf.wallflux_diagnostics: could not write " + m_filename +
                " at step " + std::to_string(step) + "; the window continues");
    }
}

void
CloudChamberWallFluxDiag::require_fresh_output () const
{
    int exists = 0;
    if (ParallelDescriptor::IOProcessor()) {
        exists = FileExists(m_filename) ? 1 : 0;
    }
    ParallelDescriptor::Bcast(&exists, 1, ParallelDescriptor::IOProcessorNumber());
    if (exists != 0) {
        Error("erf.wallflux_diagnostics: " + m_filename + " already exists; a fresh run "
              "refuses to overwrite or append to it. Move it away or run elsewhere.");
    }
}

bool
CloudChamberWallFluxDiag::write_row (double time,
                                     const erf_wallflux_diag::ChannelFaceValues& means)
{
    std::ofstream out;
    if (!m_header_written) {
        if (FileExists(m_filename)) { return false; }
        out.open(m_filename, std::ios::out | std::ios::trunc);
        if (!out.is_open()) { return false; }
        out << erf_wallflux_diag::header_line() << '\n';
        m_header_written = out.good();
    } else {
        out.open(m_filename, std::ios::out | std::ios::app);
        if (!out.is_open()) { return false; }
    }
    out << erf_wallflux_diag::format_row(time, means) << '\n';
    out.flush();
    return out.good();
}
