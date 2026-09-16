#include "ERF_ShocDriver.H"
#include "ERF_ShocImplicit.H"

#include "ERF_IndexDefines.H"

#include <AMReX_BLProfiler.H>
#include <AMReX_Gpu.H>

#include <algorithm>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

using namespace amrex;

bool
shoc_boxarray_spans_full_height (const BoxArray& ba, const Box& domain)
{
    const int dom_klo = domain.smallEnd(2);
    const int dom_khi = domain.bigEnd(2);
    for (int ibox = 0; ibox < ba.size(); ++ibox) {
        const Box& bx = ba[ibox];
        if (bx.smallEnd(2) != dom_klo || bx.bigEnd(2) != dom_khi) {
            return false;
        }
    }
    return true;
}

namespace
{
    constexpr int k_shoc_vertical_diff_comp = EddyDiff::Mom_v;
    constexpr int k_shoc_vertical_diff_count = EddyDiff::Q_v - EddyDiff::Mom_v + 1;

    struct ShocBadColumnReport
    {
        amrex::Real score = 0.0_rt;
        std::string reason;
        int mfi_index = -1;
        int i = -1;
        int j = -1;
        int k = -1;
        int kk = -1;
        int ic = -1;
    };

    template <typename Fab>
    FArrayBox
    copy_fab_to_host (const Fab& src)
    {
#ifdef AMREX_USE_GPU
        FArrayBox dst(src.box(), src.nComp(), amrex::The_Pinned_Arena());
        const std::size_t nbytes = static_cast<std::size_t>(dst.size()) * sizeof(Real);
        if (src.arena()->isManaged() || src.arena()->isDevice()) {
            static bool printed_gpu_path_marker = false;
            if (!printed_gpu_path_marker) {
                amrex::Print() << "NATIVE_SHOC_DEBUG_COPY_HELPER_GPU_PATH\n";
                printed_gpu_path_marker = true;
            }
            amrex::Gpu::dtoh_memcpy_async(dst.dataPtr(), src.dataPtr(), nbytes);
            amrex::Gpu::streamSynchronize();
        } else {
            std::memcpy(dst.dataPtr(), src.dataPtr(), nbytes);
        }
        return dst;
#else
        FArrayBox dst(src.box(), src.nComp());
        dst.template copy<RunOn::Host>(src);
        return dst;
#endif
    }

    void print_shoc_debug_settings_once (const ShocRuntimeOptions& opts)
    {
        static bool printed = false;
        if (printed) {
            return;
        }
        printed = true;

        std::ostringstream os;
        os << std::boolalpha;
        os << "SHOC debug thresholds:\n"
           << "  theta_tend=" << opts.debug_bad_column_theta_tend_threshold << "\n"
           << "  q_tend=" << opts.debug_bad_column_q_tend_threshold << "\n"
           << "  brunt=" << opts.debug_bad_column_brunt_threshold << "\n"
           << "  min_dz=" << opts.debug_bad_column_min_dz << "\n"
           << "  scalar_moment=" << opts.debug_bad_column_scalar_moment_threshold << "\n"
           << "SHOC debug switches:\n"
           << "  debug_bad_column=" << opts.debug_bad_column << "\n"
           << "  debug_bad_column_abort=" << opts.debug_bad_column_abort << "\n"
           << "  debug_disable_pdf_cloud_increment=" << opts.debug_disable_pdf_cloud_increment << "\n"
           << "  debug_disable_theta_state_update=" << opts.debug_disable_theta_state_update << "\n"
           << "  debug_disable_moisture_state_update=" << opts.debug_disable_moisture_state_update << "\n"
           << "  debug_disable_tke_state_update=" << opts.debug_disable_tke_state_update << "\n";
        amrex::Print() << os.str();
    }

    void warn_if_shoc_debug_overrides_active_once (const ShocRuntimeOptions& opts)
    {
        static bool warned = false;
        if (warned) {
            return;
        }
        if (!opts.debug_disable_pdf_cloud_increment &&
            !opts.debug_disable_theta_state_update &&
            !opts.debug_disable_moisture_state_update &&
            !opts.debug_disable_tke_state_update) {
            return;
        }
        warned = true;

        std::ostringstream os;
        os << std::boolalpha;
        os << "WARNING: native SHOC debug override active:\n"
           << "  debug_disable_pdf_cloud_increment = " << opts.debug_disable_pdf_cloud_increment << "\n"
           << "  debug_disable_theta_state_update = " << opts.debug_disable_theta_state_update << "\n"
           << "  debug_disable_moisture_state_update = " << opts.debug_disable_moisture_state_update << "\n"
           << "  debug_disable_tke_state_update = " << opts.debug_disable_tke_state_update << "\n"
           << "These options are for debugging only and change native SHOC numerics.\n";
        amrex::Print() << os.str();
    }

    void require_full_height_shoc_boxes (const BoxArray& ba, const Box& domain)
    {
        // Tile-local SHOC packing and writeback still assume each MFIter tile
        // spans a complete vertical column.
        if (!shoc_boxarray_spans_full_height(ba, domain)) {
            amrex::Abort(
                "Native SHOC requires each BoxArray box on a SHOC-active level to span the full vertical domain. Z-split boxes are not supported. Increase the vertical max_grid_size/blocking size, or use full-column refined grids.");
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real weighted_linear_interp (Real x0, Real x1, Real y0, Real y1, Real x)
    {
        const Real denom = x1 - x0;
        if (amrex::Math::abs(denom) <= 1.0e-12_rt) {
            return 0.5_rt * (y0 + y1);
        }
        return y0 + (y1 - y0) * (x - x0) / denom;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void diagnose_profile_at_height (const Array4<const Real>& zt,
                                     const Array4<const Real>& dz,
                                     const Array4<const Real>& theta_v,
                                     const Array4<const Real>& wthv_sec,
                                     const Array4<const Real>& u,
                                     const Array4<const Real>& v,
                                     const Array4<const Real>& w,
                                     int ic,
                                     int nlev,
                                     Real H,
                                     int jump_ncell,
                                     Real min_delta_theta_v,
                                     Real& u_at_pblh,
                                     Real& v_at_pblh,
                                     Real& w_at_pblh,
                                     Real& wthv_at_pblh,
                                     Real& delta_theta_v,
                                     Real& we_flux_jump)
    {
        u_at_pblh = shoc_entrainment::missing_value;
        v_at_pblh = shoc_entrainment::missing_value;
        w_at_pblh = shoc_entrainment::missing_value;
        wthv_at_pblh = shoc_entrainment::missing_value;
        delta_theta_v = shoc_entrainment::missing_value;
        we_flux_jump = shoc_entrainment::missing_value;

        if (nlev < 2 || jump_ncell < 1 || !(min_delta_theta_v > 0.0_rt) ||
            !shoc_entrainment::finite(H)) {
            return;
        }

        int bracket = -1;
        for (int k = 0; k + 1 < nlev; ++k) {
            if (zt(ic,k,0) <= H && H <= zt(ic,k+1,0)) {
                bracket = k;
                break;
            }
        }
        if (bracket < 0) {
            return;
        }

        const Real z0 = zt(ic,bracket,0);
        const Real z1 = zt(ic,bracket+1,0);
        const Real denom = z1 - z0;
        if (!(denom > 0.0_rt) || !shoc_entrainment::finite(denom)) {
            return;
        }
        u_at_pblh = weighted_linear_interp(z0, z1, u(ic,bracket,0), u(ic,bracket+1,0), H);
        v_at_pblh = weighted_linear_interp(z0, z1, v(ic,bracket,0), v(ic,bracket+1,0), H);
        w_at_pblh = weighted_linear_interp(z0, z1, w(ic,bracket,0), w(ic,bracket+1,0), H);
        wthv_at_pblh = weighted_linear_interp(z0, z1,
                                              wthv_sec(ic,bracket,0),
                                              wthv_sec(ic,bracket+1,0), H);
        if (!shoc_entrainment::finite(u_at_pblh) ||
            !shoc_entrainment::finite(v_at_pblh) ||
            !shoc_entrainment::finite(w_at_pblh) ||
            !shoc_entrainment::finite(wthv_at_pblh)) {
            u_at_pblh = v_at_pblh = w_at_pblh =
                wthv_at_pblh = shoc_entrainment::missing_value;
            return;
        }

        Real below_sum = 0.0_rt, below_weight = 0.0_rt;
        Real above_sum = 0.0_rt, above_weight = 0.0_rt;
        int below_count = 0, above_count = 0;
        for (int distance = 0;
             distance < nlev && (below_count < jump_ncell || above_count < jump_ncell);
             ++distance) {
            const int kb = bracket - distance;
            if (below_count < jump_ncell && kb >= 0 && zt(ic,kb,0) < H) {
                const Real z = zt(ic,kb,0);
                const Real layer_dz = dz(ic,kb,0);
                const Real thv = theta_v(ic,kb,0);
                if (!shoc_entrainment::finite(z) || !shoc_entrainment::finite(layer_dz) ||
                    !shoc_entrainment::finite(thv) || !(layer_dz > 0.0_rt)) {
                    return;
                }
                below_sum += layer_dz * thv;
                below_weight += layer_dz;
                ++below_count;
            }
            const int ka = bracket + 1 + distance;
            if (above_count < jump_ncell && ka < nlev && zt(ic,ka,0) > H) {
                const Real z = zt(ic,ka,0);
                const Real layer_dz = dz(ic,ka,0);
                const Real thv = theta_v(ic,ka,0);
                if (!shoc_entrainment::finite(z) || !shoc_entrainment::finite(layer_dz) ||
                    !shoc_entrainment::finite(thv) || !(layer_dz > 0.0_rt)) {
                    return;
                }
                above_sum += layer_dz * thv;
                above_weight += layer_dz;
                ++above_count;
            }
        }
        if (below_count != jump_ncell || above_count != jump_ncell ||
            !(below_weight > 0.0_rt) || !(above_weight > 0.0_rt)) {
            u_at_pblh = v_at_pblh = w_at_pblh =
                wthv_at_pblh = shoc_entrainment::missing_value;
            return;
        }
        delta_theta_v = above_sum / above_weight - below_sum / below_weight;
        if (!shoc_entrainment::finite(delta_theta_v) ||
            delta_theta_v <= min_delta_theta_v) {
            delta_theta_v = wthv_at_pblh = we_flux_jump = shoc_entrainment::missing_value;
            return;
        }
        we_flux_jump = -wthv_at_pblh / delta_theta_v;
        if (!shoc_entrainment::finite(we_flux_jump)) {
            delta_theta_v = wthv_at_pblh = we_flux_jump = shoc_entrainment::missing_value;
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int face_neighbor_cell (int idx, int domlo, int domhi, bool periodic)
    {
        return periodic ? idx : shoc_clamp(idx, domlo, domhi);
    }

    void seed_carried_turbulence_impl (
        Box const& xy_box,
        ShocColumnLayout layout,
        Array4<const Real> const& rho_host,
        Array4<const Real> const& carried,
        Array4<const Real> const& host_diff,
        Array4<Real> const& tk,
        Array4<Real> const& tkh,
        bool prev_turb_valid)
    {
        ParallelFor(xy_box, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            const int ic = shoc_column_index(layout, i, j);
            for (int kk = 0; kk < layout.nlev; ++kk) {
                const int k = layout.kmin + kk;
                if (prev_turb_valid) {
                    tk(ic,kk,0) = carried(i,j,k,0);
                    tkh(ic,kk,0) = carried(i,j,k,1);
                } else {
                    const Real rho = amrex::max(rho_host(i,j,k,Rho_comp), 1.0e-12_rt);
                    tk(ic,kk,0) = amrex::max(0.0_rt, host_diff(i,j,k,EddyDiff::Mom_v) / rho);
                    tkh(ic,kk,0) = amrex::max(0.0_rt, host_diff(i,j,k,EddyDiff::Theta_v) / rho);
                }
            }
        });
    }

    void store_carried_turbulence_impl (
        Box const& xy_box,
        ShocColumnLayout layout,
        Array4<Real> const& carried,
        Array4<const Real> const& tk,
        Array4<const Real> const& tkh)
    {
        ParallelFor(xy_box, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            const int ic = shoc_column_index(layout, i, j);
            for (int kk = 0; kk < layout.nlev; ++kk) {
                const int k = layout.kmin + kk;
                carried(i,j,k,0) = tk(ic,kk,0);
                carried(i,j,k,1) = tkh(ic,kk,0);
            }
        });
    }

    void seed_carried_buoyancy_flux_impl (
        Box const& xy_box,
        ShocColumnLayout layout,
        Array4<const Real> const& carried,
        Array4<Real> const& wthv_sec)
    {
        ParallelFor(xy_box, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            const int ic = shoc_column_index(layout, i, j);
            for (int kk = 0; kk < layout.nlev; ++kk) {
                const int k = layout.kmin + kk;
                wthv_sec(ic,kk,0) = carried(i,j,k,0);
            }
        });
    }

    void store_carried_buoyancy_flux_impl (
        Box const& xy_box,
        ShocColumnLayout layout,
        Array4<Real> const& carried,
        Array4<const Real> const& wthv_sec)
    {
        ParallelFor(xy_box, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            const int ic = shoc_column_index(layout, i, j);
            for (int kk = 0; kk < layout.nlev; ++kk) {
                const int k = layout.kmin + kk;
                carried(i,j,k,0) = wthv_sec(ic,kk,0);
            }
        });
    }

    void sync_face_multifab_impl (MultiFab& mf, const Geometry& geom)
    {
        // OverrideSync reconciles duplicate values on overlapping face-centered
        // tiles before SHOC reads or writes them. It is a no-op for cell-centered
        // MultiFabs, so the helper can be reused for the face-centered state and
        // tendency fields without special casing the index type.
        mf.OverrideSync(geom.periodicity());
    }

    void apply_state_update_cons_impl (MultiFab& cons,
                                       const MultiFab& theta_tend_cc,
                                       const MultiFab& qv_tend_cc,
                                       const MultiFab& qc_tend_cc,
                                       const MultiFab& qi_tend_cc,
                                       const MultiFab& tke_tend_cc,
                                       int qv_comp,
                                       int qc_comp,
                                       int qi_comp,
                                       const ShocRuntimeOptions& opts,
                                       double dt_d)
    {
        Real dt = static_cast<Real>(dt_d);

        const bool has_qv = shoc_valid_comp(qv_comp, cons.nComp());
        const bool has_qc = shoc_valid_comp(qc_comp, cons.nComp());
        const bool has_qi = shoc_valid_comp(qi_comp, cons.nComp());

        for (MFIter mfi(cons, false); mfi.isValid(); ++mfi) {
            const auto rho = cons.const_array(mfi);
            auto cc = cons.array(mfi);
            const auto theta_tend = theta_tend_cc.const_array(mfi);
            const auto qv_tend = qv_tend_cc.const_array(mfi);
            const auto qc_tend = qc_tend_cc.const_array(mfi);
            const auto qi_tend = qi_tend_cc.const_array(mfi);
            const auto tke_tend = tke_tend_cc.const_array(mfi);

            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const Real rho_val = rho(i,j,k,Rho_comp);
                if (!opts.debug_disable_theta_state_update) {
                    cc(i,j,k,RhoTheta_comp) += rho_val * dt * theta_tend(i,j,k);
                }
                if (!opts.debug_disable_moisture_state_update) {
                    if (has_qv) {
                        cc(i,j,k,qv_comp) = amrex::max(
                            cc(i,j,k,qv_comp) + rho_val * dt * qv_tend(i,j,k), Real(0.));
                    }
                    if (has_qc) {
                        cc(i,j,k,qc_comp) = amrex::max(
                            cc(i,j,k,qc_comp) + rho_val * dt * qc_tend(i,j,k), Real(0.));
                    }
                    if (has_qi) {
                        cc(i,j,k,qi_comp) = amrex::max(
                            cc(i,j,k,qi_comp) + rho_val * dt * qi_tend(i,j,k), Real(0.));
                    }
                }
                if (!opts.debug_disable_tke_state_update) {
                    cc(i,j,k,RhoKE_comp) = amrex::max(
                        cc(i,j,k,RhoKE_comp) + rho_val * dt * tke_tend(i,j,k), Real(0.));
                }
            });
        }
    }

    void apply_state_update_face_velocity_impl (MultiFab& vel,
                                                const MultiFab& vel_tend,
                                                double dt_d)
    {
        Real dt = static_cast<Real>(dt_d);
        for (MFIter mfi(vel_tend, false); mfi.isValid(); ++mfi) {
            auto v = vel.array(mfi);
            const auto tend = vel_tend.const_array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                v(i,j,k) += dt * tend(i,j,k);
            });
        }
    }
}

ShocDriver::ShocDriver (int lev, const SolverChoice& solver_choice)
    : m_lev(lev),
      m_moisture_indices(solver_choice.moisture_indices)
{
    m_opts.transport_mode     = solver_choice.shoc_transport_mode;
    m_opts.momentum_transport = solver_choice.shoc_momentum_transport;

    read_shoc_runtime_options(m_opts);
    validate_shoc_runtime_options(m_opts);
    warn_if_shoc_debug_overrides_active_once(m_opts);
}

void
ShocDriver::ensure_storage (const MultiFab& cons,
                            const MultiFab& xvel,
                            const MultiFab& yvel,
                            const MultiFab& eddy_diffs)
{
    const bool storage_changed =
        !m_theta_tend_cc.isDefined() ||
        !m_u_tend_fc.isDefined() ||
        !m_v_tend_fc.isDefined() ||
        !m_eddy_coeffs_cc.isDefined() ||
        m_theta_tend_cc.boxArray() != cons.boxArray() ||
        m_theta_tend_cc.DistributionMap() != cons.DistributionMap() ||
        m_u_tend_fc.boxArray() != xvel.boxArray() ||
        m_u_tend_fc.DistributionMap() != xvel.DistributionMap() ||
        m_v_tend_fc.boxArray() != yvel.boxArray() ||
        m_v_tend_fc.DistributionMap() != yvel.DistributionMap() ||
        m_eddy_coeffs_cc.boxArray() != eddy_diffs.boxArray() ||
        m_eddy_coeffs_cc.DistributionMap() != eddy_diffs.DistributionMap() ||
        (m_opts.diagnose_entrainment && !m_H_cc.isDefined());

    if (storage_changed) {
        m_column_workspaces.clear();
        m_prev_turb_valid = false;
        m_theta_tend_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_qv_tend_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_qc_tend_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_qi_tend_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_tke_tend_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_u_tend_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 1);
        m_v_tend_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 1);
        m_u_tend_fc.define(xvel.boxArray(), xvel.DistributionMap(), 1, 0);
        m_v_tend_fc.define(yvel.boxArray(), yvel.DistributionMap(), 1, 0);
        m_eddy_coeffs_cc.define(eddy_diffs.boxArray(), eddy_diffs.DistributionMap(),
                                EddyDiff::NumDiffs, 0);
        m_prev_turb_cc.define(cons.boxArray(), cons.DistributionMap(), 2, 0);
        m_prev_wthv_sec_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_consumed_sens_flux_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_consumed_laten_flux_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_pblh_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_shoc_ustar_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_shoc_olen_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_shoc_cldfrac_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_shoc_ql_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_shoc_ql2_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_shoc_cond_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_w_sec_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_wqls_sec_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_wthv_sec_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_thl_sec_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_qw_sec_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_qwthl_sec_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_wthl_sec_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_wqw_sec_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_w3_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_brunt_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_isotropy_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_shear_prod_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_buoy_prod_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        m_diss_tke_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
        if (m_opts.diagnose_entrainment) {
            m_H_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 1);
            m_prev_H_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
            m_we_kinematic_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
            m_we_flux_jump_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
            m_pblh_tendency_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
            m_pblh_hadv_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
            m_w_at_pblh_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
            m_delta_theta_v_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
            m_wthv_at_pblh_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
            m_u_at_pblh_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
            m_v_at_pblh_cc.define(cons.boxArray(), cons.DistributionMap(), 1, 0);
            m_H_cc.setVal(shoc_entrainment::missing_value);
            m_prev_H_cc.setVal(shoc_entrainment::missing_value);
            m_we_kinematic_cc.setVal(shoc_entrainment::missing_value);
            m_we_flux_jump_cc.setVal(shoc_entrainment::missing_value);
            m_pblh_tendency_cc.setVal(shoc_entrainment::missing_value);
            m_pblh_hadv_cc.setVal(shoc_entrainment::missing_value);
            m_w_at_pblh_cc.setVal(shoc_entrainment::missing_value);
            m_delta_theta_v_cc.setVal(shoc_entrainment::missing_value);
            m_wthv_at_pblh_cc.setVal(shoc_entrainment::missing_value);
            m_u_at_pblh_cc.setVal(shoc_entrainment::missing_value);
            m_v_at_pblh_cc.setVal(shoc_entrainment::missing_value);
            m_prev_H_valid = false;
        }
        m_prev_turb_cc.setVal(0.0);
        m_prev_wthv_sec_cc.setVal(0.0);
    }
    m_theta_tend_cc.setVal(0.0);
    m_qv_tend_cc.setVal(0.0);
    m_qc_tend_cc.setVal(0.0);
    m_qi_tend_cc.setVal(0.0);
    m_tke_tend_cc.setVal(0.0);
    m_u_tend_cc.setVal(0.0);
    m_v_tend_cc.setVal(0.0);
    m_u_tend_fc.setVal(0.0);
    m_v_tend_fc.setVal(0.0);
    m_eddy_coeffs_cc.setVal(0.0);
    m_consumed_sens_flux_cc.setVal(0.0);
    m_consumed_laten_flux_cc.setVal(0.0);
    m_pblh_cc.setVal(0.0);
    m_shoc_ustar_cc.setVal(0.0);
    m_shoc_olen_cc.setVal(0.0);
    m_shoc_cldfrac_cc.setVal(0.0);
    m_shoc_ql_cc.setVal(0.0);
    m_shoc_ql2_cc.setVal(0.0);
    m_shoc_cond_cc.setVal(0.0);
    m_w_sec_cc.setVal(0.0);
    m_wqls_sec_cc.setVal(0.0);
    m_wthv_sec_cc.setVal(0.0);
    m_thl_sec_cc.setVal(0.0);
    m_qw_sec_cc.setVal(0.0);
    m_qwthl_sec_cc.setVal(0.0);
    m_wthl_sec_cc.setVal(0.0);
    m_wqw_sec_cc.setVal(0.0);
    m_w3_cc.setVal(0.0);
    m_brunt_cc.setVal(0.0);
    m_isotropy_cc.setVal(0.0);
    m_shear_prod_cc.setVal(0.0);
    m_buoy_prod_cc.setVal(0.0);
    m_diss_tke_cc.setVal(0.0);
    if (m_opts.diagnose_entrainment) {
        m_H_cc.setVal(shoc_entrainment::missing_value);
        m_we_kinematic_cc.setVal(shoc_entrainment::missing_value);
        m_we_flux_jump_cc.setVal(shoc_entrainment::missing_value);
        m_pblh_tendency_cc.setVal(shoc_entrainment::missing_value);
        m_pblh_hadv_cc.setVal(shoc_entrainment::missing_value);
        m_w_at_pblh_cc.setVal(shoc_entrainment::missing_value);
        m_delta_theta_v_cc.setVal(shoc_entrainment::missing_value);
        m_wthv_at_pblh_cc.setVal(shoc_entrainment::missing_value);
        m_u_at_pblh_cc.setVal(shoc_entrainment::missing_value);
        m_v_at_pblh_cc.setVal(shoc_entrainment::missing_value);
    }
}

void
ShocDriver::seed_carried_turbulence (ShocColumnData& col,
                                     const MFIter& mfi,
                                     const MultiFab& cons,
                                     const MultiFab& eddy_diffs) const
{
    const auto rho_host = cons.const_array(mfi);
    const auto carried = m_prev_turb_cc.const_array(mfi);
    const auto host_diff = eddy_diffs.const_array(mfi);
    auto tk = col.tk.array();
    auto tkh = col.tkh.array();
    const bool prev_turb_valid = m_prev_turb_valid;
    const auto layout = col.layout;
    const Box xy_box = amrex::makeSlab(mfi.validbox(), 2, layout.kmin);

    seed_carried_turbulence_impl(xy_box, layout, rho_host, carried, host_diff,
                                 tk, tkh, prev_turb_valid);
}

void
ShocDriver::store_carried_turbulence (const ShocColumnData& col,
                                      const MFIter& mfi)
{
    auto carried = m_prev_turb_cc.array(mfi);
    const auto tk = col.tk.const_array();
    const auto tkh = col.tkh.const_array();
    const auto layout = col.layout;
    const Box xy_box = amrex::makeSlab(mfi.validbox(), 2, layout.kmin);

    store_carried_turbulence_impl(xy_box, layout, carried, tk, tkh);
}

void
ShocDriver::seed_carried_buoyancy_flux (ShocColumnData& col,
                                        const MFIter& mfi) const
{
    if (!m_prev_turb_valid) {
        return;
    }

    const auto carried = m_prev_wthv_sec_cc.const_array(mfi);
    auto wthv_sec = col.wthv_sec.array();
    const auto layout = col.layout;
    const Box xy_box = amrex::makeSlab(mfi.validbox(), 2, layout.kmin);

    seed_carried_buoyancy_flux_impl(xy_box, layout, carried, wthv_sec);
}

void
ShocDriver::store_carried_buoyancy_flux (const ShocColumnData& col,
                                         const MFIter& mfi)
{
    auto carried = m_prev_wthv_sec_cc.array(mfi);
    const auto wthv_sec = col.wthv_sec.const_array();
    const auto layout = col.layout;
    const Box xy_box = amrex::makeSlab(mfi.validbox(), 2, layout.kmin);

    store_carried_buoyancy_flux_impl(xy_box, layout, carried, wthv_sec);
}

void
ShocDriver::advance (MultiFab& cons,
                     MultiFab& xvel,
                     MultiFab& yvel,
                     MultiFab& zvel,
                     MultiFab* tau13,
                     MultiFab* tau23,
                     MultiFab* hfx3,
                     MultiFab* qfx3,
                     MultiFab* eddy_diffs,
                     MultiFab& z_phys_nd,
                     const MultiFab& mapfac_mx,
                     const MultiFab& mapfac_my,
                     const Geometry& geom,
                     double dt)
{
    BL_PROFILE("SHOC::advance");

    require_full_height_shoc_boxes(cons.boxArray(), geom.Domain());

    m_cons_ptr = &cons;
    m_hfx3_ptr = hfx3;
    m_qfx3_ptr = qfx3;
    m_tau13_ptr = tau13;
    m_tau23_ptr = tau23;
    m_eddy_diffs_ptr = eddy_diffs;

    sync_face_multifab_impl(xvel, geom);
    sync_face_multifab_impl(yvel, geom);
    if (m_tau13_ptr) {
        sync_face_multifab_impl(*m_tau13_ptr, geom);
    }
    if (m_tau23_ptr) {
        sync_face_multifab_impl(*m_tau23_ptr, geom);
    }

    {
        BL_PROFILE("SHOC::advance::ensure_storage");
        ensure_storage(cons, xvel, yvel, *eddy_diffs);
    }

    std::string error_message;
    if (!shoc_driver_state_update_layout_supported(m_moisture_indices, cons.nComp(),
                                                    error_message)) {
        amrex::Abort(error_message.c_str());
    }

    // Keep one SHOC column workspace per local box so GPU writeback kernels
    // can finish reading a box's scratch before it is reused.
    const amrex::Long n_local = cons.local_size();
    AMREX_ALWAYS_ASSERT(n_local >= amrex::Long{0});
    const amrex::Long current_size = static_cast<amrex::Long>(m_column_workspaces.size());
    if (current_size < n_local) {
        const amrex::Long old_size = current_size;
        m_column_workspaces.resize(static_cast<std::size_t>(n_local));
        for (amrex::Long n = old_size; n < n_local; ++n) {
            m_column_workspaces[static_cast<std::size_t>(n)] =
                std::make_unique<ShocColumnWorkspace>();
        }
    }

    amrex::Long local_index = 0;
    for (MFIter mfi(cons, false); mfi.isValid(); ++mfi, ++local_index) {
        AMREX_ALWAYS_ASSERT(
            local_index < static_cast<amrex::Long>(m_column_workspaces.size()));
        ShocColumnWorkspace& workspace =
            *m_column_workspaces[static_cast<std::size_t>(local_index)];
        const Box& vbx = mfi.validbox();
        const ShocColumnLayout active_layout = make_shoc_layout(vbx, geom);
        {
            BL_PROFILE("SHOC::advance::define_column_data");
            workspace.ensure_capacity(active_layout, amrex::The_Async_Arena(),
                                      shoc::default_init_run_on());
        }
        ShocColumnData& col = workspace.col;
        {
            BL_PROFILE("SHOC::advance::preprocess");
            ShocPreprocess::fill_columns(col, mfi, cons, xvel, yvel, zvel,
                                         hfx3, qfx3, tau13, tau23, z_phys_nd, geom,
                                         m_moisture_indices);
        }
        {
            BL_PROFILE("SHOC::advance::seed_carried_buoyancy_flux");
            seed_carried_buoyancy_flux(col, mfi);
        }
        {
            BL_PROFILE("SHOC::advance::seed_carried_turbulence");
            seed_carried_turbulence(col, mfi, cons, *eddy_diffs);
        }
        const auto dx = geom.CellSizeArray();
        BL_PROFILE("SHOC::advance::cache_baseline_state");
        ShocImplicit::cache_baseline_state(col);
        ShocDiagnostics::diagnose_pre_implicit(col, m_opts, dx[0], dx[1], dt);
        BL_PROFILE("SHOC::advance::implicit");
        ShocImplicit::advance_implicit_state(col, m_opts, dt);
        ShocDiagnostics::diagnose_post_implicit(col, m_opts, dt);
        BL_PROFILE("SHOC::advance::finalize");
        ShocImplicit::finalize_from_pdf(col, m_opts, dt);
        BL_PROFILE("SHOC::advance::debug_bad_column");
        debug_check_bad_column(col, mfi, z_phys_nd, hfx3, qfx3, tau13, tau23, geom, dt);
        {
            BL_PROFILE("SHOC::advance::store_carried_buoyancy_flux");
            store_carried_buoyancy_flux(col, mfi);
        }
        {
            BL_PROFILE("SHOC::advance::store_carried_turbulence");
            store_carried_turbulence(col, mfi);
        }

        auto tk_arr = m_eddy_coeffs_cc.array(mfi);
        auto th_tend = m_theta_tend_cc.array(mfi);
        auto qv_tend = m_qv_tend_cc.array(mfi);
        auto qc_tend = m_qc_tend_cc.array(mfi);
        auto qi_tend = m_qi_tend_cc.array(mfi);
        auto tke_tend = m_tke_tend_cc.array(mfi);
        auto u_tend_cc = m_u_tend_cc.array(mfi);
        auto v_tend_cc = m_v_tend_cc.array(mfi);
        auto pblh_arr = m_pblh_cc.array(mfi);
        auto consumed_sens_flux_arr = m_consumed_sens_flux_cc.array(mfi);
        auto consumed_laten_flux_arr = m_consumed_laten_flux_cc.array(mfi);
        auto shoc_ustar_arr = m_shoc_ustar_cc.array(mfi);
        auto shoc_olen_arr = m_shoc_olen_cc.array(mfi);
        auto shoc_cldfrac_arr = m_shoc_cldfrac_cc.array(mfi);
        auto shoc_ql_arr = m_shoc_ql_cc.array(mfi);
        auto shoc_ql2_arr = m_shoc_ql2_cc.array(mfi);
        auto shoc_cond_arr = m_shoc_cond_cc.array(mfi);
        auto w_sec_arr = m_w_sec_cc.array(mfi);
        auto wqls_sec_arr = m_wqls_sec_cc.array(mfi);
        auto wthv_sec_arr = m_wthv_sec_cc.array(mfi);
        auto thl_sec_arr = m_thl_sec_cc.array(mfi);
        auto qw_sec_arr = m_qw_sec_cc.array(mfi);
        auto qwthl_sec_arr = m_qwthl_sec_cc.array(mfi);
        auto wthl_sec_arr = m_wthl_sec_cc.array(mfi);
        auto wqw_sec_arr = m_wqw_sec_cc.array(mfi);
        auto w3_arr = m_w3_cc.array(mfi);
        auto brunt_arr = m_brunt_cc.array(mfi);
        auto isotropy_arr = m_isotropy_cc.array(mfi);
        auto shear_prod_arr = m_shear_prod_cc.array(mfi);
        auto buoy_prod_arr = m_buoy_prod_cc.array(mfi);
        auto diss_tke_arr = m_diss_tke_cc.array(mfi);
        auto H_arr = m_opts.diagnose_entrainment ? m_H_cc.array(mfi) : Array4<Real>{};
        auto we_flux_jump_arr = m_opts.diagnose_entrainment ? m_we_flux_jump_cc.array(mfi) : Array4<Real>{};
        auto w_at_pblh_arr = m_opts.diagnose_entrainment ? m_w_at_pblh_cc.array(mfi) : Array4<Real>{};
        auto delta_theta_v_arr = m_opts.diagnose_entrainment ? m_delta_theta_v_cc.array(mfi) : Array4<Real>{};
        auto wthv_at_pblh_arr = m_opts.diagnose_entrainment ? m_wthv_at_pblh_cc.array(mfi) : Array4<Real>{};
        auto u_at_pblh_arr = m_opts.diagnose_entrainment ? m_u_at_pblh_cc.array(mfi) : Array4<Real>{};
        auto v_at_pblh_arr = m_opts.diagnose_entrainment ? m_v_at_pblh_cc.array(mfi) : Array4<Real>{};

        const auto shoc_mix = col.shoc_mix.const_array();
        const auto pblh = col.pblh.const_array();
        const auto surf_sens_flux = col.surf_sens_flux.const_array();
        const auto surf_lat_flux = col.surf_lat_flux.const_array();
        const auto ustar = col.ustar.const_array();
        const auto obklen = col.obklen.const_array();
        const auto shoc_cldfrac = col.shoc_cldfrac.const_array();
        const auto shoc_ql = col.shoc_ql.const_array();
        const auto shoc_ql2 = col.shoc_ql2.const_array();
        const auto shoc_cond = col.shoc_cond.const_array();
        const auto w_sec = col.w_sec.const_array();
        const auto wqls_sec = col.wqls_sec.const_array();
        const auto wthv_sec = col.wthv_sec.const_array();
        const auto thl_sec = col.thl_sec.const_array();
        const auto qw_sec = col.qw_sec.const_array();
        const auto qwthl_sec = col.qwthl_sec.const_array();
        const auto wthl_sec = col.wthl_sec.const_array();
        const auto wqw_sec = col.wqw_sec.const_array();
        const auto w3 = col.w3.const_array();
        const auto brunt = col.brunt.const_array();
        const auto isotropy = col.isotropy.const_array();
        const auto shear_prod = col.shear_prod.const_array();
        const auto buoy_prod = col.buoy_prod.const_array();
        const auto diss_tke = col.diss_tke.const_array();
        const auto col_u = col.u.const_array();
        const auto col_v = col.v.const_array();
        const auto col_w = col.w.const_array();
        const auto col_theta_v = col.theta_v.const_array();
        const auto col_dz = col.dz.const_array();
        const bool diagnose_entrainment = m_opts.diagnose_entrainment;
        const int entrainment_jump_ncell = m_opts.entrainment_jump_ncell;
        const Real entrainment_min_delta_theta_v = m_opts.entrainment_min_delta_theta_v;
        const auto tk = col.tk.const_array();
        const auto tkh = col.tkh.const_array();
        const auto rho = col.rho.const_array();
        const auto zt = col.zt.const_array();
        const auto zi = col.zi.const_array();
        const auto col_theta_tend = col.theta_tend.const_array();
        const auto col_qv_tend = col.qv_tend.const_array();
        const auto col_qc_tend = col.qc_tend.const_array();
        const auto col_qi_tend = col.qi_tend.const_array();
        const auto col_tke_tend = col.tke_tend.const_array();
        const auto col_u_tend = col.u_tend.const_array();
        const auto col_v_tend = col.v_tend.const_array();
        const auto layout = col.layout;
        const Box xy_box = amrex::makeSlab(vbx, 2, layout.kmin);
        {
            BL_PROFILE("SHOC::advance::writeback");
            ParallelFor(xy_box, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                const int ic = shoc_column_index(layout, i, j);
                const int k0 = layout.kmin;
                const Real rho_sfc = rho(ic,0,0);
                for (int kk = 0; kk < layout.nlev; ++kk) {
                    const int k = layout.kmin + kk;
                    const Real l = shoc_mix(ic,kk,0);
                    const Real km = tk(ic,kk,0);
                    const Real kh = tkh(ic,kk,0);

                    tk_arr(i,j,k,EddyDiff::Mom_v) = rho(ic,kk,0) * km;
                    tk_arr(i,j,k,EddyDiff::Theta_v) = rho(ic,kk,0) * kh;
                    tk_arr(i,j,k,EddyDiff::KE_v) = rho(ic,kk,0) * kh;
                    tk_arr(i,j,k,EddyDiff::Scalar_v) = rho(ic,kk,0) * kh;
                    tk_arr(i,j,k,EddyDiff::Q_v) = rho(ic,kk,0) * kh;
                    tk_arr(i,j,k,EddyDiff::Turb_lengthscale) = l;

                    th_tend(i,j,k) = col_theta_tend(ic,kk,0);
                    qv_tend(i,j,k) = col_qv_tend(ic,kk,0);
                    qc_tend(i,j,k) = col_qc_tend(ic,kk,0);
                    qi_tend(i,j,k) = col_qi_tend(ic,kk,0);
                    tke_tend(i,j,k) = col_tke_tend(ic,kk,0);
                    u_tend_cc(i,j,k) = col_u_tend(ic,kk,0);
                    v_tend_cc(i,j,k) = col_v_tend(ic,kk,0);
                }
                // Motivation: Native SHOC state_update consumes SurfaceLayer
                // fluxes before the host diffusion path clears the overlapping
                // SFS arrays. Preserve the consumed conservative fluxes here
                // so 2D diagnostics report the surface forcing SHOC actually
                // used, not the cleared host arrays.
                consumed_sens_flux_arr(i,j,k0) = rho_sfc * surf_sens_flux(ic,0,0);
                consumed_laten_flux_arr(i,j,k0) = rho_sfc * surf_lat_flux(ic,0,0);
                shoc_ustar_arr(i,j,k0) = ustar(ic,0,0);
                shoc_olen_arr(i,j,k0) = obklen(ic,0,0);
            });

            // Keep diagnostic writeback in a separate launch so the kernel
            // argument block fits CUDA's 4 KiB limit on pre-Volta devices.
            ParallelFor(xy_box, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                const int ic = shoc_column_index(layout, i, j);
                for (int kk = 0; kk < layout.nlev; ++kk) {
                    const int k = layout.kmin + kk;
                    // Native SHOC pblh is diagnosed as meters AGL and is
                    // copied directly into the plotfile diagnostic field.
                    pblh_arr(i,j,k) = pblh(ic,0,0);
                    shoc_cldfrac_arr(i,j,k) = shoc_cldfrac(ic,kk,0);
                    shoc_ql_arr(i,j,k) = shoc_ql(ic,kk,0);
                    shoc_ql2_arr(i,j,k) = shoc_ql2(ic,kk,0);
                    shoc_cond_arr(i,j,k) = shoc_cond(ic,kk,0);
                    w_sec_arr(i,j,k) = w_sec(ic,kk,0);
                    wqls_sec_arr(i,j,k) = wqls_sec(ic,kk,0);
                    wthv_sec_arr(i,j,k) = wthv_sec(ic,kk,0);
                    brunt_arr(i,j,k) = brunt(ic,kk,0);
                    isotropy_arr(i,j,k) = isotropy(ic,kk,0);
                    thl_sec_arr(i,j,k) = weighted_linear_interp(zi(ic,kk,0), zi(ic,kk+1,0),
                                                                thl_sec(ic,kk,0), thl_sec(ic,kk+1,0),
                                                                zt(ic,kk,0));
                    qw_sec_arr(i,j,k) = weighted_linear_interp(zi(ic,kk,0), zi(ic,kk+1,0),
                                                               qw_sec(ic,kk,0), qw_sec(ic,kk+1,0),
                                                               zt(ic,kk,0));
                    qwthl_sec_arr(i,j,k) = weighted_linear_interp(zi(ic,kk,0), zi(ic,kk+1,0),
                                                                  qwthl_sec(ic,kk,0), qwthl_sec(ic,kk+1,0),
                                                                  zt(ic,kk,0));
                    wthl_sec_arr(i,j,k) = weighted_linear_interp(zi(ic,kk,0), zi(ic,kk+1,0),
                                                                 wthl_sec(ic,kk,0), wthl_sec(ic,kk+1,0),
                                                                 zt(ic,kk,0));
                    wqw_sec_arr(i,j,k) = weighted_linear_interp(zi(ic,kk,0), zi(ic,kk+1,0),
                                                                wqw_sec(ic,kk,0), wqw_sec(ic,kk+1,0),
                                                                zt(ic,kk,0));
                    w3_arr(i,j,k) = weighted_linear_interp(zi(ic,kk,0), zi(ic,kk+1,0),
                                                           w3(ic,kk,0), w3(ic,kk+1,0),
                                                           zt(ic,kk,0));
                    shear_prod_arr(i,j,k) = shear_prod(ic,kk,0);
                    buoy_prod_arr(i,j,k) = buoy_prod(ic,kk,0);
                    diss_tke_arr(i,j,k) = diss_tke(ic,kk,0);
                }
                if (diagnose_entrainment) {
                    const Real H = zi(ic,0,0) + pblh(ic,0,0);
                    H_arr(i,j,k0) = H;
                    Real u_at, v_at, w_at, wthv_at, delta_thv, we_flux;
                    diagnose_profile_at_height(zt, col_dz,
                                               col_theta_v, wthv_sec, col_u, col_v, col_w,
                                               ic, layout.nlev, H,
                                               entrainment_jump_ncell,
                                               entrainment_min_delta_theta_v,
                                               u_at, v_at, w_at, wthv_at, delta_thv, we_flux);
                    u_at_pblh_arr(i,j,k0) = u_at;
                    v_at_pblh_arr(i,j,k0) = v_at;
                    w_at_pblh_arr(i,j,k0) = w_at;
                    wthv_at_pblh_arr(i,j,k0) = wthv_at;
                    delta_theta_v_arr(i,j,k0) = delta_thv;
                    we_flux_jump_arr(i,j,k0) = we_flux;
                    for (int kk = 1; kk < layout.nlev; ++kk) {
                        const int kt = layout.kmin + kk;
                        H_arr(i,j,kt) = H;
                        u_at_pblh_arr(i,j,kt) = u_at;
                        v_at_pblh_arr(i,j,kt) = v_at;
                        w_at_pblh_arr(i,j,kt) = w_at;
                        wthv_at_pblh_arr(i,j,kt) = wthv_at;
                        delta_theta_v_arr(i,j,kt) = delta_thv;
                        we_flux_jump_arr(i,j,kt) = we_flux;
                    }
                }
            });
        }

    }

    if (m_opts.diagnose_entrainment) {
        BL_PROFILE("SHOC::advance::entrainment_diagnostics");
        update_entrainment_diagnostics(mapfac_mx, mapfac_my, geom, dt);
    }

    {
        BL_PROFILE("SHOC::advance::tendency_interpolation");
        m_u_tend_cc.FillBoundary(geom.periodicity());
        m_v_tend_cc.FillBoundary(geom.periodicity());

        const auto dom = geom.Domain();
        const int ilo = dom.smallEnd(0);
        const int ihi = dom.bigEnd(0);
        const int jlo = dom.smallEnd(1);
        const int jhi = dom.bigEnd(1);
        const bool xper = geom.isPeriodic(0);
        const bool yper = geom.isPeriodic(1);

        for (MFIter mfi(m_u_tend_cc, false); mfi.isValid(); ++mfi) {
            const Box& cc_bx = mfi.validbox();
            const Box xface_bx = amrex::surroundingNodes(cc_bx, 0);
            const Box yface_bx = amrex::surroundingNodes(cc_bx, 1);
            const auto u_cc = m_u_tend_cc.const_array(mfi);
            const auto v_cc = m_v_tend_cc.const_array(mfi);
            auto u_fc = m_u_tend_fc.array(mfi);
            auto v_fc = m_v_tend_fc.array(mfi);

            ParallelFor(xface_bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int il = face_neighbor_cell(i - 1, ilo, ihi, xper);
                const int ir = face_neighbor_cell(i,     ilo, ihi, xper);
                u_fc(i,j,k) = 0.5_rt * (u_cc(il,j,k) + u_cc(ir,j,k));
            });

            ParallelFor(yface_bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int jb = face_neighbor_cell(j - 1, jlo, jhi, yper);
                const int jt = face_neighbor_cell(j,     jlo, jhi, yper);
                v_fc(i,j,k) = 0.5_rt * (v_cc(i,jb,k) + v_cc(i,jt,k));
            });
        }
    }

    sync_face_multifab_impl(m_u_tend_fc, geom);
    sync_face_multifab_impl(m_v_tend_fc, geom);

    BL_PROFILE("SHOC::advance::state_update");
    apply_state_update(cons, xvel, yvel, dt);

    if (uses_momentum_state_update()) {
        sync_face_multifab_impl(xvel, geom);
        sync_face_multifab_impl(yvel, geom);
    }

    m_prev_turb_valid = true;
    ++m_advance_calls;
    if (m_opts.debug_summary) {
        BL_PROFILE("SHOC::advance::debug_summary");
        print_debug_summary(dt);
    }
}

void
ShocDriver::update_entrainment_diagnostics (const MultiFab& mapfac_mx,
                                            const MultiFab& mapfac_my,
                                            const Geometry& geom,
                                            double dt)
{
    m_H_cc.FillBoundary(geom.periodicity());

    const Box dom = geom.Domain();
    const int ilo = dom.smallEnd(0);
    const int ihi = dom.bigEnd(0);
    const int jlo = dom.smallEnd(1);
    const int jhi = dom.bigEnd(1);
    const bool xper = geom.isPeriodic(0);
    const bool yper = geom.isPeriodic(1);
    const Real dxinv = geom.InvCellSize(0);
    const Real dyinv = geom.InvCellSize(1);
    const bool history_valid = m_prev_H_valid;
    const bool include_hadv = m_opts.entrainment_include_horizontal_advection;
    const Real dt_real = static_cast<Real>(dt);

    for (MFIter mfi(m_H_cc, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const auto H = m_H_cc.const_array(mfi);
        const auto Hprev = m_prev_H_cc.const_array(mfi);
        const auto mx = mapfac_mx.const_array(mfi);
        const auto my = mapfac_my.const_array(mfi);
        const auto u = m_u_at_pblh_cc.const_array(mfi);
        const auto v = m_v_at_pblh_cc.const_array(mfi);
        const auto w = m_w_at_pblh_cc.const_array(mfi);
        auto we = m_we_kinematic_cc.array(mfi);
        auto tendency = m_pblh_tendency_cc.array(mfi);
        auto hadv = m_pblh_hadv_cc.array(mfi);
        const int klo = dom.smallEnd(2);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (k != klo) {
                return;
            }
            const int il = xper ? i - 1 : amrex::max(ilo, i - 1);
            const int ir = xper ? i + 1 : amrex::min(ihi, i + 1);
            const int jb = yper ? j - 1 : amrex::max(jlo, j - 1);
            const int jt = yper ? j + 1 : amrex::min(jhi, j + 1);

            const Real dHdx = (i == ilo && !xper)
                ? mx(i,j,0) * (H(ir,j,k) - H(i,j,k)) * dxinv
                : (i == ihi && !xper)
                    ? mx(i,j,0) * (H(i,j,k) - H(il,j,k)) * dxinv
                    : mx(i,j,0) * (H(ir,j,k) - H(il,j,k)) * (0.5_rt * dxinv);
            const Real dHdy = (j == jlo && !yper)
                ? my(i,j,0) * (H(i,jt,k) - H(i,j,k)) * dyinv
                : (j == jhi && !yper)
                    ? my(i,j,0) * (H(i,j,k) - H(i,jb,k)) * dyinv
                    : my(i,j,0) * (H(i,jt,k) - H(i,jb,k)) * (0.5_rt * dyinv);

            const auto result = shoc_entrainment::kinematic(
                H(i,j,k), Hprev(i,j,k), dt_real, u(i,j,k), v(i,j,k), w(i,j,k),
                dHdx, dHdy, history_valid, include_hadv);
            we(i,j,k) = result.we_kinematic;
            tendency(i,j,k) = result.pblh_tendency;
            hadv(i,j,k) = result.pblh_hadv;
            for (int kk = klo + 1; kk <= dom.bigEnd(2); ++kk) {
                we(i,j,kk) = result.we_kinematic;
                tendency(i,j,kk) = result.pblh_tendency;
                hadv(i,j,kk) = result.pblh_hadv;
            }
        });
    }

    MultiFab::Copy(m_prev_H_cc, m_H_cc, 0, 0, 1, 0);
    m_prev_H_valid = true;
}

void
ShocDriver::set_eddy_diffs () const
{
    BL_PROFILE("SHOC::set_eddy_diffs");

    AMREX_ALWAYS_ASSERT(m_eddy_diffs_ptr != nullptr);
    m_eddy_diffs_ptr->setVal(0.0, k_shoc_vertical_diff_comp, k_shoc_vertical_diff_count,
                             m_eddy_diffs_ptr->nGrow());
    MultiFab::Copy(*m_eddy_diffs_ptr, m_eddy_coeffs_cc,
                   EddyDiff::Turb_lengthscale, EddyDiff::Turb_lengthscale, 1, 0);
}

void
ShocDriver::set_diff_stresses () const
{
    BL_PROFILE("SHOC::set_diff_stresses");

    if (!m_hfx3_ptr || !m_qfx3_ptr) {
        return;
    }

    for (MFIter mfi(*m_hfx3_ptr, false); mfi.isValid(); ++mfi) {
        const Box& vbx_cc = mfi.validbox();
        const Box& vbx_xz = convert(vbx_cc, IntVect(1,0,1));
        const Box& vbx_yz = convert(vbx_cc, IntVect(0,1,1));
        auto hfx = m_hfx3_ptr->array(mfi);
        auto qfx = m_qfx3_ptr->array(mfi);

        ParallelFor(vbx_cc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            hfx(i,j,k) = 0.0;
            qfx(i,j,k) = 0.0;
        });

        if (m_tau13_ptr && m_tau23_ptr) {
            auto tau13 = m_tau13_ptr->array(mfi);
            auto tau23 = m_tau23_ptr->array(mfi);
            ParallelFor(vbx_xz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                tau13(i,j,k) = 0.0;
            });
            ParallelFor(vbx_yz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                tau23(i,j,k) = 0.0;
            });
        }
    }
}

void
ShocDriver::add_fast_tend (Vector<MultiFab>& S_rhs) const
{
    AMREX_ALWAYS_ASSERT(m_cons_ptr != nullptr);

    amrex::ignore_unused(S_rhs);
}

void
ShocDriver::add_slow_tend (const MFIter& mfi,
                           const Box& tbx,
                           const Array4<Real>& cell_rhs) const
{
    AMREX_ALWAYS_ASSERT(m_cons_ptr != nullptr);
    amrex::ignore_unused(mfi, tbx, cell_rhs);
}

bool
ShocDriver::uses_momentum_state_update () const
{
    return shoc_uses_momentum_state_update(m_opts.momentum_transport);
}

bool
ShocDriver::disables_momentum_transport () const
{
    return shoc_disables_momentum_transport(m_opts.momentum_transport);
}

bool
ShocDriver::owns_scalar_surface_fluxes () const
{
    return true;
}

void
ShocDriver::apply_state_update (MultiFab& cons,
                                MultiFab& xvel,
                                MultiFab& yvel,
                                double dt) const
{
    AMREX_ALWAYS_ASSERT(dt > 0.0);

    apply_state_update_cons_impl(cons,
                                 m_theta_tend_cc,
                                 m_qv_tend_cc,
                                 m_qc_tend_cc,
                                 m_qi_tend_cc,
                                 m_tke_tend_cc,
                                 m_moisture_indices.qv,
                                 m_moisture_indices.qc,
                                 m_moisture_indices.qi,
                                 m_opts,
                                 dt);

    if (uses_momentum_state_update()) {
        apply_state_update_face_velocity_impl(xvel, m_u_tend_fc, dt);
        apply_state_update_face_velocity_impl(yvel, m_v_tend_fc, dt);
    }
}

void
ShocDriver::debug_check_bad_column (const ShocColumnData& col,
                                    const MFIter& mfi,
                                    const MultiFab& z_phys_nd,
                                    const MultiFab* hfx3,
                                    const MultiFab* qfx3,
                                    const MultiFab* tau13,
                                    const MultiFab* tau23,
                                    const Geometry& geom,
                                    double dt_d) const
{
    amrex::ignore_unused(geom, dt_d);

    if (!m_opts.debug_bad_column) {
        return;
    }

    AMREX_ALWAYS_ASSERT(m_cons_ptr != nullptr);
    print_shoc_debug_settings_once(m_opts);

    const auto cons_host = copy_fab_to_host((*m_cons_ptr)[mfi]);
    const auto z_host = copy_fab_to_host(z_phys_nd[mfi]);
    const auto theta_tend_host = copy_fab_to_host(col.theta_tend);
    const auto qv_tend_host = copy_fab_to_host(col.qv_tend);
    const auto qc_tend_host = copy_fab_to_host(col.qc_tend);
    const auto qi_tend_host = copy_fab_to_host(col.qi_tend);
    const auto tke_tend_host = copy_fab_to_host(col.tke_tend);
    const auto u_tend_host = copy_fab_to_host(col.u_tend);
    const auto v_tend_host = copy_fab_to_host(col.v_tend);
    const auto p_mid_host = copy_fab_to_host(col.p_mid);
    const auto p_int_host = copy_fab_to_host(col.p_int);
    const auto zt_host = copy_fab_to_host(col.zt);
    const auto zi_host = copy_fab_to_host(col.zi);
    const auto dz_host = copy_fab_to_host(col.dz);
    const auto rho_host = copy_fab_to_host(col.rho);
    const auto theta_host = copy_fab_to_host(col.theta);
    const auto thetal_host = copy_fab_to_host(col.thetal);
    const auto theta_base_host = copy_fab_to_host(col.theta_base);
    const auto thetal_base_host = copy_fab_to_host(col.thetal_base);
    const auto qv_base_host = copy_fab_to_host(col.qv_base);
    const auto qc_base_host = copy_fab_to_host(col.qc_base);
    const auto qi_base_host = copy_fab_to_host(col.qi_base);
    const auto tke_base_host = copy_fab_to_host(col.tke_base_state);
    const auto theta_v_host = copy_fab_to_host(col.theta_v);
    const auto qv_host = copy_fab_to_host(col.qv);
    const auto qc_host = copy_fab_to_host(col.qc);
    const auto qi_host = copy_fab_to_host(col.qi);
    const auto qw_host = copy_fab_to_host(col.qw);
    const auto tabs_host = copy_fab_to_host(col.tabs);
    const auto exner_host = copy_fab_to_host(col.exner);
    const auto host_dse_host = copy_fab_to_host(col.host_dse);
    const auto pblh_host = copy_fab_to_host(col.pblh);
    const auto obklen_host = copy_fab_to_host(col.obklen);
    const auto ustar_host = copy_fab_to_host(col.ustar);
    const auto shoc_mix_host = copy_fab_to_host(col.shoc_mix);
    const auto brunt_host = copy_fab_to_host(col.brunt);
    const auto isotropy_host = copy_fab_to_host(col.isotropy);
    const auto tk_host = copy_fab_to_host(col.tk);
    const auto tkh_host = copy_fab_to_host(col.tkh);
    const auto shear_prod_host = copy_fab_to_host(col.shear_prod);
    const auto buoy_prod_host = copy_fab_to_host(col.buoy_prod);
    const auto diss_tke_host = copy_fab_to_host(col.diss_tke);
    const auto tke_host = copy_fab_to_host(col.tke);
    const auto w_sec_host = copy_fab_to_host(col.w_sec);
    const auto wthv_sec_host = copy_fab_to_host(col.wthv_sec);
    const auto shoc_cldfrac_host = copy_fab_to_host(col.shoc_cldfrac);
    const auto shoc_ql_host = copy_fab_to_host(col.shoc_ql);
    const auto shoc_ql2_host = copy_fab_to_host(col.shoc_ql2);
    const auto shoc_cond_host = copy_fab_to_host(col.shoc_cond);
    const auto shoc_evap_host = copy_fab_to_host(col.shoc_evap);
    const auto wqls_sec_host = copy_fab_to_host(col.wqls_sec);
    const auto thl_sec_host = copy_fab_to_host(col.thl_sec);
    const auto qw_sec_host = copy_fab_to_host(col.qw_sec);
    const auto qwthl_sec_host = copy_fab_to_host(col.qwthl_sec);
    const auto wthl_sec_host = copy_fab_to_host(col.wthl_sec);
    const auto wqw_sec_host = copy_fab_to_host(col.wqw_sec);
    const auto uw_sec_host = copy_fab_to_host(col.uw_sec);
    const auto vw_sec_host = copy_fab_to_host(col.vw_sec);
    const auto wtke_sec_host = copy_fab_to_host(col.wtke_sec);
    const auto w3_host = copy_fab_to_host(col.w3);
    const auto surf_sens_flux_host = copy_fab_to_host(col.surf_sens_flux);
    const auto surf_lat_flux_host = copy_fab_to_host(col.surf_lat_flux);
    const auto surf_tau_u_host = copy_fab_to_host(col.surf_tau_u);
    const auto surf_tau_v_host = copy_fab_to_host(col.surf_tau_v);
    const auto hfx3_host = hfx3 ? copy_fab_to_host((*hfx3)[mfi]) : FArrayBox();
    const auto qfx3_host = qfx3 ? copy_fab_to_host((*qfx3)[mfi]) : FArrayBox();
    const auto tau13_host = tau13 ? copy_fab_to_host((*tau13)[mfi]) : FArrayBox();
    const auto tau23_host = tau23 ? copy_fab_to_host((*tau23)[mfi]) : FArrayBox();

    const auto z_arr = z_host.const_array();
    const auto cons_arr = cons_host.const_array();
    const auto theta_tend_arr = theta_tend_host.const_array();
    const auto qv_tend_arr = qv_tend_host.const_array();
    const auto qc_tend_arr = qc_tend_host.const_array();
    const auto qi_tend_arr = qi_tend_host.const_array();
    const auto tke_tend_arr = tke_tend_host.const_array();
    const auto u_tend_arr = u_tend_host.const_array();
    const auto v_tend_arr = v_tend_host.const_array();
    const auto dz_arr = dz_host.const_array();

    std::vector<ShocBadColumnReport> reports;
    reports.reserve(static_cast<std::size_t>(col.layout.ncell * col.layout.nlev));

    auto maybe_add = [&] (int i, int j, int k, int kk, int ic,
                          const char* name, Real value, Real threshold) {
        if (!std::isfinite(value)) {
            reports.push_back(ShocBadColumnReport{
                std::numeric_limits<Real>::infinity(),
                std::string(name) + " nonfinite",
                mfi.index(), i, j, k, kk, ic
            });
            return;
        }
        const Real abs_value = amrex::Math::abs(value);
        const Real ratio = abs_value / threshold;
        if (ratio > 1.0_rt) {
            std::ostringstream oss;
            oss << name << " |value|=" << abs_value << " threshold=" << threshold;
            reports.push_back(ShocBadColumnReport{ratio, oss.str(), mfi.index(), i, j, k, kk, ic});
        }
    };

    auto maybe_add_geom = [&] (int i, int j, int k, int kk, int ic, Real dz_val) {
        if (!std::isfinite(dz_val)) {
            reports.push_back(ShocBadColumnReport{
                std::numeric_limits<Real>::infinity(),
                "dz nonfinite",
                mfi.index(), i, j, k, kk, ic
            });
        } else if (dz_val <= 0.0_rt) {
            reports.push_back(ShocBadColumnReport{
                std::numeric_limits<Real>::infinity(),
                "dz <= 0",
                mfi.index(), i, j, k, kk, ic
            });
        } else if (dz_val < m_opts.debug_bad_column_min_dz) {
            std::ostringstream oss;
            oss << "dz below minimum |dz|=" << dz_val
                << " min_dz=" << m_opts.debug_bad_column_min_dz;
            reports.push_back(ShocBadColumnReport{
                m_opts.debug_bad_column_min_dz / dz_val,
                oss.str(),
                mfi.index(), i, j, k, kk, ic
            });
        }
    };

    for (int j = col.layout.jmin; j < col.layout.jmin + col.layout.ny; ++j) {
        for (int i = col.layout.imin; i < col.layout.imin + col.layout.nx; ++i) {
            const int ic = shoc_column_index(col.layout, i, j);
            for (int kk = 0; kk < col.layout.nlev; ++kk) {
                const int k = col.layout.kmin + kk;
                const Real dz_val = dz_arr(ic, kk, 0);
                maybe_add_geom(i, j, k, kk, ic, dz_val);

                maybe_add(i, j, k, kk, ic, "theta_tend", theta_tend_arr(ic, kk, 0), m_opts.debug_bad_column_theta_tend_threshold);
                maybe_add(i, j, k, kk, ic, "qv_tend", qv_tend_arr(ic, kk, 0), m_opts.debug_bad_column_q_tend_threshold);
                maybe_add(i, j, k, kk, ic, "qc_tend", qc_tend_arr(ic, kk, 0), m_opts.debug_bad_column_q_tend_threshold);
                maybe_add(i, j, k, kk, ic, "qi_tend", qi_tend_arr(ic, kk, 0), m_opts.debug_bad_column_q_tend_threshold);
                maybe_add(i, j, k, kk, ic, "brunt", brunt_host.const_array()(ic, kk, 0), m_opts.debug_bad_column_brunt_threshold);
                maybe_add(i, j, k, kk, ic, "thl_sec", thl_sec_host.const_array()(ic, kk, 0), m_opts.debug_bad_column_scalar_moment_threshold);
                maybe_add(i, j, k, kk, ic, "qw_sec", qw_sec_host.const_array()(ic, kk, 0), m_opts.debug_bad_column_scalar_moment_threshold);
                maybe_add(i, j, k, kk, ic, "qwthl_sec", qwthl_sec_host.const_array()(ic, kk, 0), m_opts.debug_bad_column_scalar_moment_threshold);
                maybe_add(i, j, k, kk, ic, "wthl_sec", wthl_sec_host.const_array()(ic, kk, 0), m_opts.debug_bad_column_scalar_moment_threshold);
                maybe_add(i, j, k, kk, ic, "wqw_sec", wqw_sec_host.const_array()(ic, kk, 0), m_opts.debug_bad_column_scalar_moment_threshold);

                // Treat NaN/Inf in the key state and diagnostic fields as bad.
                const Real key_values[] = {
                    rho_host.const_array()(ic, kk, 0),
                    theta_host.const_array()(ic, kk, 0),
                    thetal_host.const_array()(ic, kk, 0),
                    theta_v_host.const_array()(ic, kk, 0),
                    qv_host.const_array()(ic, kk, 0),
                    qc_host.const_array()(ic, kk, 0),
                    qi_host.const_array()(ic, kk, 0),
                    qw_host.const_array()(ic, kk, 0),
                    tabs_host.const_array()(ic, kk, 0),
                    exner_host.const_array()(ic, kk, 0),
                    p_mid_host.const_array()(ic, kk, 0),
                    host_dse_host.const_array()(ic, kk, 0),
                    pblh_host.const_array()(ic, 0, 0),
                    obklen_host.const_array()(ic, 0, 0),
                    ustar_host.const_array()(ic, 0, 0),
                    shoc_mix_host.const_array()(ic, kk, 0),
                    isotropy_host.const_array()(ic, kk, 0),
                    tk_host.const_array()(ic, kk, 0),
                    tkh_host.const_array()(ic, kk, 0),
                    shear_prod_host.const_array()(ic, kk, 0),
                    buoy_prod_host.const_array()(ic, kk, 0),
                    diss_tke_host.const_array()(ic, kk, 0),
                    w_sec_host.const_array()(ic, kk, 0),
                    wthv_sec_host.const_array()(ic, kk, 0),
                    shoc_cldfrac_host.const_array()(ic, kk, 0),
                    shoc_ql_host.const_array()(ic, kk, 0),
                    shoc_ql2_host.const_array()(ic, kk, 0),
                    shoc_cond_host.const_array()(ic, kk, 0),
                    shoc_evap_host.const_array()(ic, kk, 0),
                    wqls_sec_host.const_array()(ic, kk, 0),
                    thl_sec_host.const_array()(ic, kk, 0),
                    qw_sec_host.const_array()(ic, kk, 0),
                    qwthl_sec_host.const_array()(ic, kk, 0),
                    wthl_sec_host.const_array()(ic, kk, 0),
                    wqw_sec_host.const_array()(ic, kk, 0),
                    uw_sec_host.const_array()(ic, kk, 0),
                    vw_sec_host.const_array()(ic, kk, 0),
                    wtke_sec_host.const_array()(ic, kk, 0),
                    w3_host.const_array()(ic, kk, 0),
                    theta_tend_arr(ic, kk, 0),
                    qv_tend_arr(ic, kk, 0),
                    qc_tend_arr(ic, kk, 0),
                    qi_tend_arr(ic, kk, 0),
                    tke_tend_arr(ic, kk, 0),
                    u_tend_arr(ic, kk, 0),
                    v_tend_arr(ic, kk, 0)
                };
                for (Real value : key_values) {
                    if (!std::isfinite(value)) {
                        reports.push_back(ShocBadColumnReport{
                            std::numeric_limits<Real>::infinity(),
                            "key field nonfinite",
                            mfi.index(), i, j, k, kk, ic
                        });
                        break;
                    }
                }
            }
        }
    }

    if (reports.empty()) {
        return;
    }

    std::stable_sort(reports.begin(), reports.end(),
                     [] (const ShocBadColumnReport& a, const ShocBadColumnReport& b) {
                         return a.score > b.score;
                     });

    const int max_reports = std::min(m_opts.debug_bad_column_max_reports,
                                     static_cast<int>(reports.size()));

    const auto rho_arr = rho_host.const_array();
    const auto theta_arr = theta_host.const_array();
    const auto thetal_arr = thetal_host.const_array();
    const auto theta_base_arr = theta_base_host.const_array();
    const auto thetal_base_arr = thetal_base_host.const_array();
    const auto qv_base_arr = qv_base_host.const_array();
    const auto qc_base_arr = qc_base_host.const_array();
    const auto qi_base_arr = qi_base_host.const_array();
    const auto theta_v_arr = theta_v_host.const_array();
    const auto qv_arr = qv_host.const_array();
    const auto qc_arr = qc_host.const_array();
    const auto qi_arr = qi_host.const_array();
    const auto qw_arr = qw_host.const_array();
    const auto tabs_arr = tabs_host.const_array();
    const auto exner_arr = exner_host.const_array();
    const auto p_mid_arr = p_mid_host.const_array();
    const auto p_int_arr = p_int_host.const_array();
    const auto host_dse_arr = host_dse_host.const_array();
    const auto pblh_arr = pblh_host.const_array();
    const auto obklen_arr = obklen_host.const_array();
    const auto ustar_arr = ustar_host.const_array();
    const auto shoc_mix_arr = shoc_mix_host.const_array();
    const auto brunt_arr = brunt_host.const_array();
    const auto isotropy_arr = isotropy_host.const_array();
    const auto tk_arr = tk_host.const_array();
    const auto tkh_arr = tkh_host.const_array();
    const auto shear_prod_arr = shear_prod_host.const_array();
    const auto buoy_prod_arr = buoy_prod_host.const_array();
    const auto diss_tke_arr = diss_tke_host.const_array();
    const auto tke_state_arr = tke_host.const_array();
    const auto tke_base_arr = tke_base_host.const_array();
    const auto w_sec_arr = w_sec_host.const_array();
    const auto wthv_sec_arr = wthv_sec_host.const_array();
    const auto shoc_cldfrac_arr = shoc_cldfrac_host.const_array();
    const auto shoc_ql_arr = shoc_ql_host.const_array();
    const auto shoc_ql2_arr = shoc_ql2_host.const_array();
    const auto shoc_cond_arr = shoc_cond_host.const_array();
    const auto shoc_evap_arr = shoc_evap_host.const_array();
    const auto wqls_sec_arr = wqls_sec_host.const_array();
    const auto thl_sec_arr = thl_sec_host.const_array();
    const auto qw_sec_arr = qw_sec_host.const_array();
    const auto qwthl_sec_arr = qwthl_sec_host.const_array();
    const auto wthl_sec_arr = wthl_sec_host.const_array();
    const auto wqw_sec_arr = wqw_sec_host.const_array();
    const auto uw_sec_arr = uw_sec_host.const_array();
    const auto vw_sec_arr = vw_sec_host.const_array();
    const auto wtke_sec_arr = wtke_sec_host.const_array();
    const auto w3_arr = w3_host.const_array();
    const auto surf_sens_flux_arr = surf_sens_flux_host.const_array();
    const auto surf_lat_flux_arr = surf_lat_flux_host.const_array();
    const auto surf_tau_u_arr = surf_tau_u_host.const_array();
    const auto surf_tau_v_arr = surf_tau_v_host.const_array();
    const bool has_hfx3 = hfx3 && hfx3_host.box().ok();
    const bool has_qfx3 = qfx3 && qfx3_host.box().ok();
    const bool has_tau13 = tau13 && tau13_host.box().ok();
    const bool has_tau23 = tau23 && tau23_host.box().ok();

    for (int n = 0; n < max_reports; ++n) {
        const auto& rep = reports[static_cast<std::size_t>(n)];
        const int i = rep.i;
        const int j = rep.j;
        const int k = rep.k;
        const int kk = rep.kk;
        const int ic = rep.ic;

        const auto node_value = [&] (int iii, int jjj, int kkk) -> Real {
            const IntVect iv(iii, jjj, kkk);
            return z_host.box().contains(iv) ? z_arr(iii, jjj, kkk) : std::numeric_limits<Real>::quiet_NaN();
        };

        const Real z_nd_ijk = node_value(i, j, k);
        const Real z_nd_ip1jk = node_value(i + 1, j, k);
        const Real z_nd_ijp1k = node_value(i, j + 1, k);
        const Real z_nd_ip1jp1k = node_value(i + 1, j + 1, k);
        const Real z_nd_ijkp1 = node_value(i, j, k + 1);
        const Real z_nd_ip1jkp1 = node_value(i + 1, j, k + 1);
        const Real z_nd_ijp1kp1 = node_value(i, j + 1, k + 1);
        const Real z_nd_ip1jp1kp1 = node_value(i + 1, j + 1, k + 1);
        const Real four_node_zlo = 0.25_rt * (z_nd_ijk + z_nd_ip1jk + z_nd_ijp1k + z_nd_ip1jp1k);
        const Real four_node_zhi = 0.25_rt * (z_nd_ijkp1 + z_nd_ip1jkp1 + z_nd_ijp1kp1 + z_nd_ip1jp1kp1);
        const Real four_node_dz = four_node_zhi - four_node_zlo;
        const Real corner_dz = z_nd_ijkp1 - z_nd_ijk;
        const Real theta_base_val = theta_base_arr(ic, kk, 0);
        const Real theta_new_val = theta_arr(ic, kk, 0);
        const Real thetal_base_val = thetal_base_arr(ic, kk, 0);
        const Real thetal_new_val = thetal_arr(ic, kk, 0);
        const Real qv_base_val = qv_base_arr(ic, kk, 0);
        const Real qc_base_val = qc_base_arr(ic, kk, 0);
        const Real qi_base_val = qi_base_arr(ic, kk, 0);
        const Real qv_new_val = qv_arr(ic, kk, 0);
        const Real qc_new_val = qc_arr(ic, kk, 0);
        const Real qi_new_val = qi_arr(ic, kk, 0);
        const Real ql_base = qc_base_val + qi_base_val;
        const Real ql_new = qc_new_val + qi_new_val;
        const Real qw_base = qv_base_val + qc_base_val + qi_base_val;
        const Real qw_new = qv_new_val + qc_new_val + qi_new_val;
        const Real delta_theta = theta_new_val - theta_base_val;
        const Real delta_qv = qv_new_val - qv_base_val;
        const Real delta_qc = qc_new_val - qc_base_val;
        const Real delta_qi = qi_new_val - qi_base_val;
        const Real delta_ql = ql_new - ql_base;
        const Real delta_qw = qw_new - qw_base;

        const Real dt = static_cast<Real>(dt_d);

        const Real dt_theta_tend = theta_tend_arr(ic, kk, 0) * dt;
        const Real dt_qv_tend = qv_tend_arr(ic, kk, 0) * dt;
        const Real dt_qc_tend = qc_tend_arr(ic, kk, 0) * dt;
        const Real dt_qi_tend = qi_tend_arr(ic, kk, 0) * dt;
        const Real dt_tke_tend = tke_tend_arr(ic, kk, 0) * dt;
        const Real cond_dt = shoc_cond_arr(ic, kk, 0) * dt;
        const Real evap_dt = shoc_evap_arr(ic, kk, 0) * dt;
        const Real tke_base_val = tke_base_arr(ic, kk, 0);
        std::ostringstream msg;
        msg << "NATIVE_SHOC_BAD_COLUMN_BEGIN\n"
            << "  rank=" << ParallelDescriptor::MyProc()
            << " level=" << m_lev
            << " shoc_call=" << (m_advance_calls + 1)
            << " mfi_index=" << rep.mfi_index
            << " box_valid_lo=(" << mfi.validbox().smallEnd(0) << ","
            << mfi.validbox().smallEnd(1) << ","
            << mfi.validbox().smallEnd(2) << ")"
            << " box_valid_hi=(" << mfi.validbox().bigEnd(0) << ","
            << mfi.validbox().bigEnd(1) << ","
            << mfi.validbox().bigEnd(2) << ")"
            << " layout.nx=" << col.layout.nx
            << " layout.ny=" << col.layout.ny
            << " layout.ncell=" << col.layout.ncell
            << " layout.nlev=" << col.layout.nlev
            << " layout.imin=" << col.layout.imin
            << " layout.jmin=" << col.layout.jmin
            << " layout.kmin=" << col.layout.kmin
            << " layout.kmax=" << col.layout.kmax
            << "\n"
            << "  i=" << i << " j=" << j << " k=" << k
            << " kk=" << kk << " ic=" << ic
            << " score=" << rep.score
            << " reason=" << rep.reason
            << "\n"
            << "  geometry z_nd(i,j,k)=" << z_nd_ijk
            << " z_nd(i+1,j,k)=" << z_nd_ip1jk
            << " z_nd(i,j+1,k)=" << z_nd_ijp1k
            << " z_nd(i+1,j+1,k)=" << z_nd_ip1jp1k
            << " z_nd(i,j,k+1)=" << z_nd_ijkp1
            << " z_nd(i+1,j,k+1)=" << z_nd_ip1jkp1
            << " z_nd(i,j+1,k+1)=" << z_nd_ijp1kp1
            << " z_nd(i+1,j+1,k+1)=" << z_nd_ip1jp1kp1
            << "\n"
            << "  four_node_zlo=" << four_node_zlo
            << " four_node_zhi=" << four_node_zhi
            << " four_node_dz=" << four_node_dz
            << " corner_dz=" << corner_dz
            << " dz=" << dz_arr(ic, kk, 0)
            << "\n"
            << "  rho=" << rho_arr(ic, kk, 0)
            << " theta=" << theta_arr(ic, kk, 0)
            << " thetal=" << thetal_arr(ic, kk, 0)
            << " theta_v=" << theta_v_arr(ic, kk, 0)
            << " qv=" << qv_arr(ic, kk, 0)
            << " qc=" << qc_arr(ic, kk, 0)
            << " qi=" << qi_arr(ic, kk, 0)
            << " qw=" << qw_arr(ic, kk, 0)
            << " tabs=" << tabs_arr(ic, kk, 0)
            << " exner=" << exner_arr(ic, kk, 0)
            << " p_mid=" << p_mid_arr(ic, kk, 0)
            << " p_int_lower=" << p_int_arr(ic, kk, 0)
            << " p_int_upper=" << p_int_arr(ic, kk + 1, 0)
            << " host_dse=" << host_dse_arr(ic, kk, 0)
            << "\n"
            << "  pblh=" << pblh_arr(ic, 0, 0)
            << " obklen=" << obklen_arr(ic, 0, 0)
            << " ustar=" << ustar_arr(ic, 0, 0)
            << " shoc_mix=" << shoc_mix_arr(ic, kk, 0)
            << " Lturb=" << shoc_mix_arr(ic, kk, 0)
            << " brunt=" << brunt_arr(ic, kk, 0)
            << " isotropy=" << isotropy_arr(ic, kk, 0)
            << " tk=" << tk_arr(ic, kk, 0)
            << " tkh=" << tkh_arr(ic, kk, 0)
            << " shear_prod=" << shear_prod_arr(ic, kk, 0)
            << " buoy_prod=" << buoy_prod_arr(ic, kk, 0)
            << " diss_tke=" << diss_tke_arr(ic, kk, 0)
            << " tke=" << tke_state_arr(ic, kk, 0)
            << "\n"
            << "  w_sec=" << w_sec_arr(ic, kk, 0)
            << " wthv_sec=" << wthv_sec_arr(ic, kk, 0)
            << " shoc_cldfrac=" << shoc_cldfrac_arr(ic, kk, 0)
            << " shoc_ql=" << shoc_ql_arr(ic, kk, 0)
            << " shoc_ql2=" << shoc_ql2_arr(ic, kk, 0)
            << " shoc_cond=" << shoc_cond_arr(ic, kk, 0)
            << " shoc_evap=" << shoc_evap_arr(ic, kk, 0)
            << " wqls_sec=" << wqls_sec_arr(ic, kk, 0)
            << "\n"
            << "  thl_sec_lower=" << thl_sec_arr(ic, kk, 0)
            << " thl_sec_upper=" << thl_sec_arr(ic, kk + 1, 0)
            << " qw_sec_lower=" << qw_sec_arr(ic, kk, 0)
            << " qw_sec_upper=" << qw_sec_arr(ic, kk + 1, 0)
            << " qwthl_sec_lower=" << qwthl_sec_arr(ic, kk, 0)
            << " qwthl_sec_upper=" << qwthl_sec_arr(ic, kk + 1, 0)
            << " wthl_sec_lower=" << wthl_sec_arr(ic, kk, 0)
            << " wthl_sec_upper=" << wthl_sec_arr(ic, kk + 1, 0)
            << " wqw_sec_lower=" << wqw_sec_arr(ic, kk, 0)
            << " wqw_sec_upper=" << wqw_sec_arr(ic, kk + 1, 0)
            << " w3_lower=" << w3_arr(ic, kk, 0)
            << " w3_upper=" << w3_arr(ic, kk + 1, 0)
            << " uw_sec_lower=" << uw_sec_arr(ic, kk, 0)
            << " uw_sec_upper=" << uw_sec_arr(ic, kk + 1, 0)
            << " vw_sec_lower=" << vw_sec_arr(ic, kk, 0)
            << " vw_sec_upper=" << vw_sec_arr(ic, kk + 1, 0)
            << " wtke_sec_lower=" << wtke_sec_arr(ic, kk, 0)
            << " wtke_sec_upper=" << wtke_sec_arr(ic, kk + 1, 0)
            << "\n"
            << "  theta_tend=" << theta_tend_arr(ic, kk, 0)
            << " qv_tend=" << qv_tend_arr(ic, kk, 0)
            << " qc_tend=" << qc_tend_arr(ic, kk, 0)
            << " qi_tend=" << qi_tend_arr(ic, kk, 0)
            << " tke_tend=" << tke_tend_arr(ic, kk, 0)
            << " u_tend=" << u_tend_arr(ic, kk, 0)
            << " v_tend=" << v_tend_arr(ic, kk, 0)
            << "\n"
            << "  baseline theta=" << theta_base_val
            << " thetal=" << thetal_base_val
            << " qv=" << qv_base_val
            << " qc=" << qc_base_val
            << " qi=" << qi_base_val
            << " ql=" << ql_base
            << " qw=" << qw_base
            << " tke=" << tke_base_val
            << "\n"
            << "  updated  theta=" << theta_new_val
            << " thetal=" << thetal_new_val
            << " qv=" << qv_new_val
            << " qc=" << qc_new_val
            << " qi=" << qi_new_val
            << " ql=" << ql_new
            << " qw=" << qw_new
            << " tke=" << tke_state_arr(ic, kk, 0)
            << "\n"
            << "  deltas   dtheta=" << delta_theta
            << " dqv=" << delta_qv
            << " dqc=" << delta_qc
            << " dqi=" << delta_qi
            << " dql=" << delta_ql
            << " dqw=" << delta_qw
            << "\n"
            << "  tend_dt  theta=" << dt_theta_tend
            << " qv=" << dt_qv_tend
            << " qc=" << dt_qc_tend
            << " qi=" << dt_qi_tend
            << " tke=" << dt_tke_tend
            << "\n"
            << "  consistency dtheta_minus_tenddt=" << (delta_theta - dt_theta_tend)
            << " dqv_minus_tenddt=" << (delta_qv - dt_qv_tend)
            << " dqc_minus_tenddt=" << (delta_qc - dt_qc_tend)
            << " dqi_minus_tenddt=" << (delta_qi - dt_qi_tend)
            << "\n"
            << "  pdf_cloud shoc_ql=" << shoc_ql_arr(ic, kk, 0)
            << " shoc_cond_dt=" << cond_dt
            << " shoc_evap_dt=" << evap_dt
            << " delta_ql=" << delta_ql
            << " delta_ql_minus_cond_minus_evap=" << (delta_ql - (cond_dt - evap_dt))
            << "\n"
            << "  surf_sens_flux=" << surf_sens_flux_arr(ic, 0, 0)
            << " surf_lat_flux=" << surf_lat_flux_arr(ic, 0, 0)
            << " surf_tau_u=" << surf_tau_u_arr(ic, 0, 0)
            << " surf_tau_v=" << surf_tau_v_arr(ic, 0, 0)
            << " rho_sfc=" << cons_arr(i, j, k, Rho_comp)
            << "\n";

        if (has_hfx3) {
            msg << "  raw hfx3(i,j,klo)=" << hfx3_host.const_array()(rep.i, rep.j, rep.k) << "\n";
        }
        if (has_qfx3) {
            msg << "  raw qfx3(i,j,klo)=" << qfx3_host.const_array()(rep.i, rep.j, rep.k) << "\n";
        }
        if (has_tau13) {
            msg << "  raw tau13(i,j,klo)=" << tau13_host.const_array()(rep.i, rep.j, rep.k)
                << " raw tau13(i+1,j,klo)=" << tau13_host.const_array()(rep.i + 1, rep.j, rep.k)
                << "\n";
        }
        if (has_tau23) {
            msg << "  raw tau23(i,j,klo)=" << tau23_host.const_array()(rep.i, rep.j, rep.k)
                << " raw tau23(i,j+1,klo)=" << tau23_host.const_array()(rep.i, rep.j + 1, rep.k)
                << "\n";
        }
        msg << "NATIVE_SHOC_BAD_COLUMN_END\n";
        amrex::AllPrint() << msg.str() << std::flush;
    }

    if (m_opts.debug_bad_column_abort) {
        amrex::Abort("Native SHOC debug_bad_column abort: bad SHOC column detected before state update");
    }
}

void
ShocDriver::print_debug_summary (double dt) const
{
    BL_PROFILE("SHOC::print_debug_summary");

    auto print_minmax = [] (const char* name, const MultiFab& mf, int comp) {
        amrex::Print() << "    " << std::setw(12) << std::left << name
                       << " min=" << mf.min(comp)
                       << " max=" << mf.max(comp) << "\n";
    };

    amrex::Print() << "SHOC debug summary:"
                   << " level=" << m_lev
                   << " call=" << m_advance_calls
                   << " dt=" << dt
                   << " transport_mode=" << shoc_transport_mode_name(m_opts.transport_mode)
                   << " momentum_transport=" << shoc_momentum_transport_name(m_opts.momentum_transport)
                   << " state_update=on"
                   << "\n";
    print_shoc_debug_settings_once(m_opts);

    print_minmax("Kmv", m_eddy_coeffs_cc, EddyDiff::Mom_v);
    print_minmax("Khv", m_eddy_coeffs_cc, EddyDiff::Theta_v);
    print_minmax("KE_v", m_eddy_coeffs_cc, EddyDiff::KE_v);
    print_minmax("Q_v", m_eddy_coeffs_cc, EddyDiff::Q_v);
    print_minmax("l_turb", m_eddy_coeffs_cc, EddyDiff::Turb_lengthscale);
    print_minmax("pblh", m_pblh_cc, 0);

    print_minmax("th_tend", m_theta_tend_cc, 0);
    print_minmax("qv_tend", m_qv_tend_cc, 0);
    print_minmax("qc_tend", m_qc_tend_cc, 0);
    print_minmax("qi_tend", m_qi_tend_cc, 0);
    print_minmax("tke_tend", m_tke_tend_cc, 0);
    print_minmax("u_tend", m_u_tend_fc, 0);
    print_minmax("v_tend", m_v_tend_fc, 0);
    print_minmax("w_sec", m_w_sec_cc, 0);
    print_minmax("thl_sec", m_thl_sec_cc, 0);
    print_minmax("qw_sec", m_qw_sec_cc, 0);
    print_minmax("qwthl_sec", m_qwthl_sec_cc, 0);
    print_minmax("wthl_sec", m_wthl_sec_cc, 0);
    print_minmax("wqw_sec", m_wqw_sec_cc, 0);
    print_minmax("w3", m_w3_cc, 0);
    print_minmax("brunt", m_brunt_cc, 0);
    print_minmax("isotropy", m_isotropy_cc, 0);
    print_minmax("shear_prod", m_shear_prod_cc, 0);
    print_minmax("buoy_prod", m_buoy_prod_cc, 0);
    print_minmax("diss_tke", m_diss_tke_cc, 0);
    print_minmax("cldfrac", m_shoc_cldfrac_cc, 0);
    print_minmax("shoc_ql", m_shoc_ql_cc, 0);
    print_minmax("shoc_ql2", m_shoc_ql2_cc, 0);
    print_minmax("shoc_cond", m_shoc_cond_cc, 0);
    print_minmax("wqls_sec", m_wqls_sec_cc, 0);
    print_minmax("wthv_sec", m_wthv_sec_cc, 0);

    if (m_hfx3_ptr)  print_minmax("hfx3_in", *m_hfx3_ptr, 0);
    if (m_qfx3_ptr)  print_minmax("qfx3_in", *m_qfx3_ptr, 0);
    if (m_tau13_ptr) print_minmax("tau13_in", *m_tau13_ptr, 0);
    if (m_tau23_ptr) print_minmax("tau23_in", *m_tau23_ptr, 0);
}
