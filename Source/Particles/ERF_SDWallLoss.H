#ifndef ERF_SD_WALL_LOSS_H_
#define ERF_SD_WALL_LOSS_H_

#ifdef ERF_USE_PARTICLES

#include <AMReX_Algorithm.H>
#include <AMReX_Math.H>
#include <AMReX_Random.H>

#include <cmath>
#include <string>
#include <type_traits>

namespace amrex { class ParmParse; }

namespace sdwl {

using amrex::Real;
enum Model : int { SDWL_Constant, SDWL_Diffusion, SDWL_MeasuredRate };
enum Method : int { SDWL_Deterministic, SDWL_Stochastic };
inline constexpr int SDWL_MAX_BINS = 64;
inline constexpr int SDWL_MAX_PROFILE_POINTS = 64;
inline constexpr Real SDWL_kB = Real(1.380649e-23); // J/K

//! Parsed in container scope; lengths, velocities, times and rates are in SI.
struct Params {
    bool enable = false;
    Method method = SDWL_Deterministic;
    Model dep_vel_model = SDWL_Constant;
    bool deactivate_on_depletion = false;
    Real deactivate_threshold = Real(1.e-30);
    Real vd_side = Real(1.e-4), vd_top = Real(1.e-4), vd_bottom = Real(1.e-4);
    bool settling_correction = false;
    Real diff_temperature = Real(293), diff_viscosity = Real(1.8e-5), diff_delta = Real(1.e-3);
    bool diff_slip_correction = true;
    Real diff_mean_free_path = Real(6.6e-8), vd_floor = 0;
    int profile_size = 0;
    Real profile_diameter[SDWL_MAX_PROFILE_POINTS]{};
    Real profile_rate_start[SDWL_MAX_PROFILE_POINTS]{};
    Real profile_rate_end[SDWL_MAX_PROFILE_POINTS]{};
    bool profile_time_varying = false;
    Real profile_time_start = 0, profile_time_end = 0;
    int num_bins = 16;
    Real bin_rmin = Real(1.e-8), bin_rmax = Real(1.e-4);
    int diag_interval = 100;
};
static_assert(std::is_standard_layout_v<Params> && std::is_trivially_copyable_v<Params>);

void read_params (amrex::ParmParse& pp, Params& p);
std::string signature (const Params& p);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real sd_cunningham_slip_correction (Real radius, Real mean_free_path) noexcept
{
    const Real kn = mean_free_path / amrex::max(radius, Real(1.e-12));
    return kn > 0 ? 1 + kn * (Real(1.257) + Real(0.4) * std::exp(Real(-1.1) / kn)) : Real(1);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real sd_brownian_diffusivity (Real radius, const Params& p) noexcept
{
    const Real r = amrex::max(radius, Real(1.e-12));
    const Real slip = p.diff_slip_correction ? sd_cunningham_slip_correction(r, p.diff_mean_free_path) : Real(1);
    constexpr Real six_pi = 6 * amrex::Math::pi<Real>();
    return SDWL_kB * p.diff_temperature * slip / (six_pi * p.diff_viscosity * r);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real sd_vd_turbulent (Real radius, const Params& p, Real constant_face_velocity) noexcept
{
    const Real vd = p.dep_vel_model == SDWL_Diffusion
        ? sd_brownian_diffusivity(radius, p) / amrex::max(p.diff_delta, Real(1.e-12)) : constant_face_velocity;
    return amrex::max(vd, p.vd_floor);
}

//! log(b/a), for 0 < a <= b. Also resolves distinct nearby floating-point knots
//! whose separately rounded logarithms are equal, without overflowing b/a.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real log_span (Real a, Real b) noexcept
{
    return b - a <= a ? std::log1p((b - a) / a) : std::log(b) - std::log(a);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real interpolate (Real a, Real b, Real f) noexcept
{
    return f <= Real(0.5) ? a + f * (b-a) : b + (1-f) * (a-b);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real sd_wall_loss_profile_rate (Real diameter, const Params& p, bool end_profile) noexcept
{
    const int n = p.profile_size;
    if (n == 0) { return 0; }
    const Real* rates = end_profile ? p.profile_rate_end : p.profile_rate_start;
    if (diameter <= p.profile_diameter[0] || n == 1) { return rates[0]; }
    if (diameter >= p.profile_diameter[n-1]) { return rates[n-1]; }
    int lo = 0, hi = n-1;
    while (hi-lo > 1) {
        const int mid = (lo+hi)/2;
        if (diameter < p.profile_diameter[mid]) { hi = mid; } else { lo = mid; }
    }
    const Real f = log_span(p.profile_diameter[lo], diameter)
                 / log_span(p.profile_diameter[lo], p.profile_diameter[hi]);
    return interpolate(rates[lo], rates[hi], amrex::min(f, Real(1)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real sd_wall_loss_measured_rate (Real radius, Real time, const Params& p) noexcept
{
    const Real diameter = 2 * radius;
    if (!p.profile_time_varying || time <= p.profile_time_start) {
        return sd_wall_loss_profile_rate(diameter, p, false);
    }
    const Real end = sd_wall_loss_profile_rate(diameter, p, true);
    if (time >= p.profile_time_end) { return end; }
    const Real f = (time-p.profile_time_start) / (p.profile_time_end-p.profile_time_start);
    return interpolate(sd_wall_loss_profile_rate(diameter, p, false), end, f);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real sd_wall_loss_coefficient (Real vd, Real area, Real volume) noexcept
{
    return vd * (area/volume);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T sd_wall_loss_probability (T k, T dt) noexcept
{
    return -std::expm1(-k*dt);
}

//! Binomial approximations for real-valued multiplicity, as in the donor.
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T sd_wall_loss_stochastic_removed (T xi, T prob, const amrex::RandomEngine& engine) noexcept
{
    const T mean = xi * prob;
    if (mean <= 0) { return 0; }
    if (prob >= 1) { return xi; }
    const T draw = mean < 30 ? T(amrex::RandomPoisson(Real(mean), engine))
        : mean + std::sqrt(mean * (1-prob)) * T(amrex::RandomNormal(Real(0), Real(1), engine));
    return amrex::max(T(0), amrex::min(xi, draw));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int sd_wall_loss_bin_index (Real radius, const Params& p) noexcept
{
    if (radius <= p.bin_rmin) { return 0; }
    if (radius >= p.bin_rmax) { return p.num_bins-1; }
    const Real f = log_span(p.bin_rmin, radius) / log_span(p.bin_rmin, p.bin_rmax);
    return amrex::min(int(p.num_bins * f), p.num_bins-1);
}

} // namespace sdwl
#endif // ERF_USE_PARTICLES
#endif // ERF_SD_WALL_LOSS_H_
