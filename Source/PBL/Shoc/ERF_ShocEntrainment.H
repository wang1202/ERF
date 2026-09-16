#ifndef ERF_SHOC_ENTRAINMENT_H_
#define ERF_SHOC_ENTRAINMENT_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <cmath>

namespace shoc_entrainment
{

constexpr amrex::Real missing_value = amrex::Real(-999.0);

struct KinematicResult
{
    amrex::Real we_kinematic = missing_value;
    amrex::Real pblh_tendency = missing_value;
    amrex::Real pblh_hadv = missing_value;
};

struct FluxJumpResult
{
    amrex::Real we_flux_jump = missing_value;
    amrex::Real delta_theta_v = missing_value;
    amrex::Real wthv_at_pblh = missing_value;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool finite (amrex::Real value) noexcept
{
    return std::isfinite(value);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
KinematicResult kinematic (amrex::Real H,
                           amrex::Real H_previous,
                           amrex::Real dt,
                           amrex::Real u_at_pblh,
                           amrex::Real v_at_pblh,
                           amrex::Real w_at_pblh,
                           amrex::Real dHdx,
                           amrex::Real dHdy,
                           bool history_valid,
                           bool include_horizontal_advection) noexcept
{
    KinematicResult result;
    if (!history_valid || !(dt > amrex::Real(0)) ||
        !finite(H) || !finite(H_previous) || !finite(w_at_pblh) ||
        !finite(u_at_pblh) || !finite(v_at_pblh) ||
        !finite(dHdx) || !finite(dHdy)) {
        return result;
    }

    result.pblh_tendency = (H - H_previous) / dt;
    result.pblh_hadv = include_horizontal_advection
        ? u_at_pblh * dHdx + v_at_pblh * dHdy : amrex::Real(0.0);
    result.we_kinematic = result.pblh_tendency + result.pblh_hadv - w_at_pblh;
    return result;
}

inline FluxJumpResult flux_jump (const amrex::Vector<amrex::Real>& z,
                                 const amrex::Vector<amrex::Real>& dz,
                                 const amrex::Vector<amrex::Real>& theta_v,
                                 const amrex::Vector<amrex::Real>& wthv_sec,
                                 amrex::Real H,
                                 int jump_ncell,
                                 amrex::Real min_delta_theta_v)
{
    FluxJumpResult result;
    const int nlev = static_cast<int>(z.size());
    if (nlev < 2 || dz.size() != z.size() || theta_v.size() != z.size() ||
        wthv_sec.size() != z.size() || jump_ncell < 1 ||
        !(min_delta_theta_v > amrex::Real(0)) || !finite(H)) {
        return result;
    }

    int bracket = -1;
    for (int k = 0; k + 1 < nlev; ++k) {
        if (z[k] <= H && H <= z[k + 1]) {
            bracket = k;
            break;
        }
    }
    if (bracket < 0) {
        return result;
    }

    const amrex::Real dz_bracket = z[bracket + 1] - z[bracket];
    if (!(dz_bracket > amrex::Real(0)) || !finite(dz_bracket) ||
        !finite(wthv_sec[bracket]) || !finite(wthv_sec[bracket + 1])) {
        return result;
    }
    result.wthv_at_pblh = wthv_sec[bracket] +
        (wthv_sec[bracket + 1] - wthv_sec[bracket]) *
        (H - z[bracket]) / dz_bracket;

    amrex::Real below_sum = 0.0;
    amrex::Real below_weight = 0.0;
    amrex::Real above_sum = 0.0;
    amrex::Real above_weight = 0.0;
    int below_count = 0;
    int above_count = 0;
    for (int distance = 0; distance < nlev &&
         (below_count < jump_ncell || above_count < jump_ncell); ++distance) {
        const int kb = bracket - distance;
        if (below_count < jump_ncell && kb >= 0 && z[kb] < H) {
            if (!finite(z[kb]) || !finite(dz[kb]) || !finite(theta_v[kb]) ||
                !(dz[kb] > amrex::Real(0))) {
                return FluxJumpResult{};
            }
            below_sum += dz[kb] * theta_v[kb];
            below_weight += dz[kb];
            ++below_count;
        }
        const int ka = bracket + 1 + distance;
        if (above_count < jump_ncell && ka < nlev && z[ka] > H) {
            if (!finite(z[ka]) || !finite(dz[ka]) || !finite(theta_v[ka]) ||
                !(dz[ka] > amrex::Real(0))) {
                return FluxJumpResult{};
            }
            above_sum += dz[ka] * theta_v[ka];
            above_weight += dz[ka];
            ++above_count;
        }
    }
    if (below_count != jump_ncell || above_count != jump_ncell ||
        !(below_weight > amrex::Real(0)) || !(above_weight > amrex::Real(0))) {
        return FluxJumpResult{};
    }

    result.delta_theta_v = above_sum / above_weight - below_sum / below_weight;
    if (!finite(result.delta_theta_v) ||
        result.delta_theta_v <= min_delta_theta_v ||
        !finite(result.wthv_at_pblh)) {
        return FluxJumpResult{};
    }
    result.we_flux_jump = -result.wthv_at_pblh / result.delta_theta_v;
    return finite(result.we_flux_jump) ? result : FluxJumpResult{};
}

} // namespace shoc_entrainment

#endif
