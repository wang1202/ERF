#include "ERF.H"
#include "ERF_EOS.H"
#include "ERF_MOSTStress.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_TerrainMetrics.H"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

using namespace amrex;

namespace {

bool
validNoahMPGroundStationSurfaceTemperature (Real t_sfc, Real t_flux, int land_mask)
{
    // The Noah-MP result scatter writes lsm_undefined for every field in an
    // unprocessed cell.  Use the temperature-flux validity as the coupling
    // gate, matching the scalar surface-flux path, and independently reject a
    // bad skin temperature.  Restrict this to land: Noah-MP does not provide
    // an ocean skin-temperature result, so ocean stations must use MOST's SST
    // based t_surf.
    return land_mask == 1 &&
           std::isfinite(t_flux) && t_flux > Real(-9990.0) &&
           t_flux < Real(0.5) * lsm_undefined &&
           std::isfinite(t_sfc) && t_sfc > Real(0.0) &&
           t_sfc < Real(0.5) * lsm_undefined;
}

} // anonymous namespace

std::string
ERF::formatGroundStationHeightLabel (Real height)
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << height;
    std::string label = oss.str();

    while (!label.empty() && label.back() == '0') { label.pop_back(); }
    if (!label.empty() && label.back() == '.') { label.pop_back(); }
    if (label.empty()) { label = "0"; }

    std::replace(label.begin(), label.end(), '.', 'p');
    return label + "m";
}

void
ERF::initializeGroundStationFile (GroundStationConfig& station)
{
    if (station.header_written) { return; }

    if (ParallelDescriptor::IOProcessor()) {
        bool write_header = true;
        {
            std::ifstream existing(station.file_name, std::ios::binary | std::ios::ate);
            write_header = (!existing.good() || existing.tellg() == std::streampos(0));
        }

        std::ofstream out(station.file_name, std::ios::out | std::ios::app);
        if (!out.good()) {
            Abort("Failed to open ground station output file");
        }

        if (write_header) {
            const Real xloc = geom[0].ProbLo(0) + (static_cast<Real>(station.i_loc) + Real(0.5)) * geom[0].CellSize(0);
            const Real yloc = geom[0].ProbLo(1) + (static_cast<Real>(station.j_loc) + Real(0.5)) * geom[0].CellSize(1);

            out << "# ERF ground station output\n";
            out << "# i_loc=" << station.i_loc
                << " j_loc=" << station.j_loc
                << " x_m=" << xloc
                << " y_m=" << yloc << "\n";
            out << "# heights_agl_m";
            for (const auto& height : station.heights) {
                out << " " << height;
            }
            out << "\n";

            out << "step time_s surface_pressure_pa accumulated_precip_mm";
            for (const auto& height : station.heights) {
                const std::string suffix = formatGroundStationHeightLabel(height);
                out << " u_" << suffix
                    << " v_" << suffix
                    << " theta_" << suffix
                    << " temp_" << suffix
                    << " qv_" << suffix
                    << " RH_" << suffix;
            }
            out << "\n";
        }
    }

    station.header_written = true;
}

int
ERF::findGroundStationLevel (const GroundStationConfig& station) const
{
    const Real xloc = geom[0].ProbLo(0) + (static_cast<Real>(station.i_loc) + Real(0.5)) * geom[0].CellSize(0);
    const Real yloc = geom[0].ProbLo(1) + (static_cast<Real>(station.j_loc) + Real(0.5)) * geom[0].CellSize(1);

    int lev_station = 0;
    for (int lev = finest_level; lev >= 0; --lev) {
        const int i_lev = static_cast<int>(std::floor((xloc - geom[lev].ProbLo(0)) * geom[lev].InvCellSize(0)));
        const int j_lev = static_cast<int>(std::floor((yloc - geom[lev].ProbLo(1)) * geom[lev].InvCellSize(1)));
        if (grids[lev].contains(IntVect(i_lev, j_lev, geom[lev].Domain().smallEnd(2)))) {
            lev_station = lev;
            break;
        }
    }

    return lev_station;
}

void
ERF::sampleGroundStationColumn (int lev,
                                const GroundStationConfig& station,
                                Vector<Real>& z_agl,
                                Vector<GroundStationSample>& samples,
                                Real& surface_pressure,
                                Real& accumulated_precip) const
{
    const Real xloc = geom[0].ProbLo(0) + (static_cast<Real>(station.i_loc) + Real(0.5)) * geom[0].CellSize(0);
    const Real yloc = geom[0].ProbLo(1) + (static_cast<Real>(station.j_loc) + Real(0.5)) * geom[0].CellSize(1);

    const int i_lev = static_cast<int>(std::floor((xloc - geom[lev].ProbLo(0)) * geom[lev].InvCellSize(0)));
    const int j_lev = static_cast<int>(std::floor((yloc - geom[lev].ProbLo(1)) * geom[lev].InvCellSize(1)));

    const Box& domain = geom[lev].Domain();
    const int klo = domain.smallEnd(2);
    const int khi = domain.bigEnd(2);
    const int nk = domain.length(2);

    Vector<Real> packed(2 + 7 * nk, Real(0.0));
    int found = 0;

    for (MFIter mfi(vars_new[lev][Vars::cons]); mfi.isValid(); ++mfi) {
        const Box& vbx = mfi.validbox();
        if (!vbx.contains(IntVect(i_lev, j_lev, klo))) {
            continue;
        }

        found = 1;

        const auto cons_arr = vars_new[lev][Vars::cons].const_array(mfi);
        const auto xvel_arr = vars_new[lev][Vars::xvel].const_array(mfi);
        const auto yvel_arr = vars_new[lev][Vars::yvel].const_array(mfi);
        const auto z_nd_arr = z_phys_nd[lev]->const_array(mfi);
        const auto z_cc_arr = z_phys_cc[lev]->const_array(mfi);

        const bool use_moisture = (solverChoice.moisture_type != MoistureType::None);
        const Real z_surf = Real(0.25) * (z_nd_arr(i_lev    , j_lev    , klo) +
                                          z_nd_arr(i_lev + 1, j_lev    , klo) +
                                          z_nd_arr(i_lev    , j_lev + 1, klo) +
                                          z_nd_arr(i_lev + 1, j_lev + 1, klo));

        // Diagnose pressure at the lowest cell center, then hydrostatically
        // extrapolate it downward to the actual terrain surface.
        const Real rho_low = cons_arr(i_lev, j_lev, klo, Rho_comp);
        const Real qv_low = use_moisture
                                ? cons_arr(i_lev, j_lev, klo, RhoQ1_comp) / rho_low
                                : Real(0.0);
        const Real p_low = getPgivenRTh(
            cons_arr(i_lev, j_lev, klo, RhoTheta_comp), qv_low);
        const Real dz_low_to_surface =
            amrex::max(z_cc_arr(i_lev, j_lev, klo) - z_surf, Real(0.0));
        const Real rho_moist_low = rho_low * (Real(1.0) + qv_low);

        packed[0] = p_low + rho_moist_low * CONST_GRAV * dz_low_to_surface;

        const int nprecip = std::min(3, static_cast<int>(qmoist[lev].size()));
        for (int n = 0; n < nprecip; ++n) {
            if (qmoist[lev][n] != nullptr) {
                const auto precip_arr = qmoist[lev][n]->const_array(mfi);
                packed[1] += precip_arr(i_lev, j_lev, klo);
            }
        }

        for (int k = klo; k <= khi; ++k) {
            const int base = 2 + 7 * (k - klo);

            const Real rho = cons_arr(i_lev, j_lev, k, Rho_comp);
            const Real qv = use_moisture ? cons_arr(i_lev, j_lev, k, RhoQ1_comp) / rho : Real(0.0);
            const Real rhotheta = cons_arr(i_lev, j_lev, k, RhoTheta_comp);
            const Real theta = rhotheta / rho;
            const Real temp = getTgivenRandRTh(rho, rhotheta, qv);
            const Real pressure = getPgivenRTh(rhotheta, qv);
            const Real u = Real(0.5) * (xvel_arr(i_lev, j_lev, k) + xvel_arr(i_lev + 1, j_lev, k));
            const Real v = Real(0.5) * (yvel_arr(i_lev, j_lev, k) + yvel_arr(i_lev, j_lev + 1, k));

            Real rh = Real(0.0);
            if (use_moisture) {
                Real qsat = Real(0.0);
                erf_qsatw(temp, pressure * Real(0.01), qsat);
                if (qsat > Real(0.0)) {
                    rh = Real(100.0) * qv / qsat;
                }
            }

            packed[base + 0] = z_cc_arr(i_lev, j_lev, k) - z_surf;
            packed[base + 1] = u;
            packed[base + 2] = v;
            packed[base + 3] = theta;
            packed[base + 4] = temp;
            packed[base + 5] = qv;
            packed[base + 6] = rh;
        }

        break;
    }

    ParallelDescriptor::ReduceIntSum(found);
    ParallelDescriptor::ReduceRealSum(packed.data(), static_cast<int>(packed.size()));

    if (found != 1) {
        Abort("Ground station location must map to exactly one valid box");
    }

    surface_pressure = packed[0];
    accumulated_precip = packed[1];

    z_agl.resize(nk);
    samples.resize(nk);
    for (int k = 0; k < nk; ++k) {
        const int base = 2 + 7 * k;
        z_agl[k] = packed[base + 0];
        samples[k].u = packed[base + 1];
        samples[k].v = packed[base + 2];
        samples[k].theta = packed[base + 3];
        samples[k].temp = packed[base + 4];
        samples[k].qv = packed[base + 5];
        samples[k].rh = packed[base + 6];
    }
}

ERF::GroundStationSample
ERF::makeMOSTGroundStationSample (const GroundStationSample& ref_sample,
                                  Real target_height,
                                  Real surface_pressure,
                                  Real ustar,
                                  Real tstar,
                                  Real qstar,
                                  Real olen,
                                  Real theta_surf,
                                  Real qv_surf,
                                  Real z0) const
{
    constexpr Real tiny = Real(1.0e-12);

    similarity_funs sfuns;

    const Real z_eval = amrex::max(target_height, amrex::max(z0 + tiny, Real(1.0e-3)));
    Real psi_m = Real(0.0);
    Real psi_h = Real(0.0);
    if (std::abs(olen) > tiny && std::abs(olen) < bogus_large_value * Real(0.5)) {
        const Real zeta = z_eval / olen;
        psi_m = sfuns.calc_psi_m2(zeta);
        psi_h = sfuns.calc_psi_h2(zeta);
    }

    const Real profile_factor = std::max(std::log(z_eval / amrex::max(z0, tiny)) - psi_m, tiny);
    const Real scalar_factor  = std::max(std::log(z_eval / amrex::max(z0, tiny)) - psi_h, tiny);

    const Real ref_speed = std::sqrt(ref_sample.u * ref_sample.u + ref_sample.v * ref_sample.v);
    const Real speed = ustar * profile_factor / KAPPA;
    const Real dir_u = (ref_speed > tiny) ? ref_sample.u / ref_speed : Real(1.0);
    const Real dir_v = (ref_speed > tiny) ? ref_sample.v / ref_speed : Real(0.0);

    GroundStationSample out;
    out.u = speed * dir_u;
    out.v = speed * dir_v;
    out.theta = theta_surf + tstar * scalar_factor / KAPPA;
    out.qv = amrex::max(Real(0.0), qv_surf + qstar * scalar_factor / KAPPA);
    out.temp = getTgivenPandTh(surface_pressure, out.theta, R_d / Cp_d);

    Real qsat = Real(0.0);
    erf_qsatw(out.temp, surface_pressure * Real(0.01), qsat);
    out.rh = (qsat > tiny) ? Real(100.0) * out.qv / qsat : Real(0.0);

    return out;
}

ERF::GroundStationSample
ERF::interpolateGroundStationSample (Real target_height,
                                     const Vector<Real>& z_agl,
                                     const Vector<GroundStationSample>& samples) const
{
    AMREX_ALWAYS_ASSERT(!z_agl.empty());
    AMREX_ALWAYS_ASSERT(z_agl.size() == samples.size());

    if (target_height <= z_agl.front()) {
        return samples.front();
    }
    if (target_height >= z_agl.back()) {
        return samples.back();
    }

    for (int k = 0; k < static_cast<int>(z_agl.size()) - 1; ++k) {
        if (target_height >= z_agl[k] && target_height <= z_agl[k + 1]) {
            const Real dz = z_agl[k + 1] - z_agl[k];
            const Real alpha = (dz > Real(0.0)) ? (target_height - z_agl[k]) / dz : Real(0.0);

            GroundStationSample out;
            out.u = (Real(1.0) - alpha) * samples[k].u + alpha * samples[k + 1].u;
            out.v = (Real(1.0) - alpha) * samples[k].v + alpha * samples[k + 1].v;
            out.theta = (Real(1.0) - alpha) * samples[k].theta + alpha * samples[k + 1].theta;
            out.temp = (Real(1.0) - alpha) * samples[k].temp + alpha * samples[k + 1].temp;
            out.qv = (Real(1.0) - alpha) * samples[k].qv + alpha * samples[k + 1].qv;
            out.rh = (Real(1.0) - alpha) * samples[k].rh + alpha * samples[k + 1].rh;
            return out;
        }
    }

    return samples.back();
}

void
ERF::writeGroundStationData (GroundStationConfig& station, Real time, int nstep)
{
    initializeGroundStationFile(station);

    const int lev = findGroundStationLevel(station);

    Vector<Real> z_agl;
    Vector<GroundStationSample> column_samples;
    Real surface_pressure = Real(0.0);
    Real accumulated_precip = Real(0.0);

    sampleGroundStationColumn(lev, station, z_agl, column_samples, surface_pressure, accumulated_precip);

    bool have_surface_layer = (m_SurfaceLayer != nullptr);
    Real ustar = Real(0.0), tstar = Real(0.0), qstar = Real(0.0);
    Real olen = bogus_large_value, theta_surf = Real(0.0), qv_surf = Real(0.0), z0 = Real(0.1);

    if (have_surface_layer) {
        const int i2d = station.i_loc;
        const int j2d = station.j_loc;
        const int mlev = 0;

        const MultiFab* noahmp_t_sfc = lsm.Has_Model()
                                         ? lsm.Get_Data_Ptr(mlev, "t_sfc")
                                         : nullptr;
        const MultiFab* noahmp_t_flux = nullptr;
        if (mlev < static_cast<int>(lsm_flux.size())) {
            for (int n = 0; n < static_cast<int>(lsm_flux_name.size()); ++n) {
                if (amrex::toLower(lsm_flux_name[n]) == "t_flux" &&
                    n < static_cast<int>(lsm_flux[mlev].size())) {
                    noahmp_t_flux = lsm_flux[mlev][n];
                    break;
                }
            }
        }
        const bool has_land_mask = mlev < static_cast<int>(lmask_lev.size()) &&
                                   !lmask_lev[mlev].empty();
        const iMultiFab* land_mask = has_land_mask ? lmask_lev[mlev][0].get() : nullptr;

        Vector<Real> surface_pack(7, Real(0.0));
        int found_surface = 0;

        for (MFIter mfi(*m_SurfaceLayer->get_u_star(mlev)); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            if (!bx.contains(IntVect(i2d, j2d, 0))) { continue; }

            found_surface = 1;

            const auto ustar_arr = m_SurfaceLayer->get_u_star(mlev)->const_array(mfi);
            const auto tstar_arr = m_SurfaceLayer->get_t_star(mlev)->const_array(mfi);
            const auto qstar_arr = m_SurfaceLayer->get_q_star(mlev)->const_array(mfi);
            const auto olen_arr  = m_SurfaceLayer->get_olen(mlev)->const_array(mfi);
            const auto tsurf_arr = m_SurfaceLayer->get_t_surf(mlev)->const_array(mfi);
            const auto qsurf_arr = m_SurfaceLayer->get_q_surf(mlev)->const_array(mfi);
            const auto z0_arr    = m_SurfaceLayer->get_z0(mlev)->const_array(mfi);
            const auto t_sfc_arr = noahmp_t_sfc ? noahmp_t_sfc->const_array(mfi)
                                                 : Array4<const Real> {};
            const auto t_flux_arr = noahmp_t_flux ? noahmp_t_flux->const_array(mfi)
                                                   : Array4<const Real> {};
            const auto lmask_arr = land_mask ? land_mask->const_array(mfi)
                                              : Array4<const int> {};

            surface_pack[0] = ustar_arr(i2d, j2d, 0);
            surface_pack[1] = tstar_arr(i2d, j2d, 0);
            surface_pack[2] = qstar_arr(i2d, j2d, 0);
            surface_pack[3] = olen_arr(i2d, j2d, 0);
            const int is_land = lmask_arr ? lmask_arr(i2d, j2d, 0) : 1;
            if (t_sfc_arr && t_flux_arr &&
                validNoahMPGroundStationSurfaceTemperature(
                    t_sfc_arr(i2d, j2d, 0), t_flux_arr(i2d, j2d, 0), is_land)) {
                // Noah-MP TSK/t_sfc is a physical skin temperature.  MOST
                // profiles use potential temperature, so convert it at the
                // station's diagnosed surface pressure before extrapolating.
                surface_pack[4] = getThgivenTandP(
                    t_sfc_arr(i2d, j2d, 0), surface_pressure, R_d / Cp_d);
            } else {
                surface_pack[4] = tsurf_arr(i2d, j2d, 0);
            }
            surface_pack[5] = qsurf_arr(i2d, j2d, 0);
            surface_pack[6] = z0_arr(i2d, j2d, 0);
            break;
        }

        ParallelDescriptor::ReduceIntSum(found_surface);
        ParallelDescriptor::ReduceRealSum(surface_pack.data(), static_cast<int>(surface_pack.size()));

        if (found_surface == 1) {
            ustar = surface_pack[0];
            tstar = surface_pack[1];
            qstar = surface_pack[2];
            olen = surface_pack[3];
            theta_surf = surface_pack[4];
            qv_surf = surface_pack[5];
            z0 = surface_pack[6];
        } else {
            have_surface_layer = false;
        }
    }

    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream out(station.file_name, std::ios::out | std::ios::app);
        if (!out.good()) {
            Abort("Failed to append ground station output");
        }

        out << nstep << " "
            << std::setprecision(timeprecision) << time << " "
            << std::setprecision(datprecision) << surface_pressure << " "
            << accumulated_precip;

        for (const auto& height : station.heights) {
            GroundStationSample sample;
            if (have_surface_layer && !z_agl.empty() && height <= z_agl.front()) {
                sample = makeMOSTGroundStationSample(column_samples.front(), height, surface_pressure,
                                                     ustar, tstar, qstar, olen, theta_surf, qv_surf, z0);
            } else {
                sample = interpolateGroundStationSample(height, z_agl, column_samples);
            }
            out << " "
                << sample.u << " "
                << sample.v << " "
                << sample.theta << " "
                << sample.temp << " "
                << sample.qv << " "
                << sample.rh;
        }
        out << "\n";
    }
}
