#include "ERF_SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

#include "ERF_SDWallLoss.H"
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParticleUtil.H>

#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>

namespace sdwl {

void read_params (amrex::ParmParse& pp, Params& p)
{
    p = Params{};
    pp.query("wall_loss.enable", p.enable);
    if (!p.enable) { return; }

    const std::string scope = pp.getPrefix().empty() ? "wall_loss" : pp.getPrefix() + ".wall_loss";
    auto remaining = amrex::ParmParse::getEntries(scope);
    remaining.erase(scope + ".enable");
    auto key = [&](const char* suffix) {
        remaining.erase(scope + "." + suffix);
        return std::string("wall_loss.") + suffix;
    };
    auto require = [&](bool ok, const char* message) {
        if (!ok) { amrex::Abort(scope + "." + message); }
    };
    auto scalar = [&](const char* suffix, Real& value, bool positive = false) {
        const bool present = pp.query(key(suffix), value);
        if (!amrex::Math::isfinite(value) || (positive ? value <= 0 : value < 0)) {
            amrex::Abort(scope + "." + suffix + (positive ? " must be finite and positive" : " must be finite and nonnegative"));
        }
        return present;
    };

    std::string method = "deterministic", model = "constant";
    pp.query(key("method"), method);
    if (method == "stochastic" || method == "binomial") { p.method = SDWL_Stochastic; }
    else { require(method == "deterministic" || method == "exponential", "method: expected deterministic/exponential or stochastic/binomial"); }
    pp.query(key("deposition_velocity_model"), model);
    if (model == "diffusion" || model == "stokes_einstein") { p.dep_vel_model = SDWL_Diffusion; }
    else if (model == "measured_rate" || model == "measured_profile") { p.dep_vel_model = SDWL_MeasuredRate; }
    else { require(model == "constant", "deposition_velocity_model: expected constant, diffusion/stokes_einstein or measured_rate/measured_profile"); }

    pp.query(key("deactivate_on_depletion"), p.deactivate_on_depletion);
    scalar("deactivate_threshold", p.deactivate_threshold);
    Real common = p.vd_side;
    if (scalar("deposition_velocity", common)) { p.vd_side = p.vd_top = p.vd_bottom = common; }
    scalar("vd_side", p.vd_side);
    scalar("vd_top", p.vd_top);
    scalar("vd_bottom", p.vd_bottom);
    scalar("vd_floor", p.vd_floor);
    pp.query(key("settling_correction"), p.settling_correction);
    scalar("diffusion_temperature", p.diff_temperature, true);
    scalar("diffusion_viscosity", p.diff_viscosity, true);
    scalar("diffusion_delta", p.diff_delta, true);
    pp.query(key("diffusion_slip_correction"), p.diff_slip_correction);
    scalar("diffusion_mean_free_path", p.diff_mean_free_path);

    pp.query(key("num_bins"), p.num_bins);
    require(p.num_bins >= 1 && p.num_bins <= SDWL_MAX_BINS, "num_bins must be in [1,64]");
    scalar("bin_rmin", p.bin_rmin, true);
    scalar("bin_rmax", p.bin_rmax, true);
    require(p.bin_rmax > p.bin_rmin, "bin_rmax must exceed bin_rmin");
    pp.query(key("diag_interval"), p.diag_interval);
    require(p.diag_interval >= 0, "diag_interval must be nonnegative (zero disables output)");

    if (p.dep_vel_model == SDWL_MeasuredRate) {
        amrex::Vector<Real> diam, first, last;
        // Check the lengths before allocating any input-sized arrays.
        const auto dk = key("profile_diameters"), fk = key("profile_rates"), lk = key("profile_rates_end");
        p.profile_size = pp.countval(dk);
        require(p.profile_size >= 1 && p.profile_size <= SDWL_MAX_PROFILE_POINTS, "profile_diameters must contain 1..64 diameters");
        require(pp.countval(fk) == p.profile_size, "profile_rates must match profile_diameters in length");
        p.profile_time_varying = pp.contains(lk);
        pp.getarr(dk, diam);
        pp.getarr(fk, first);
        if (p.profile_time_varying) {
            require(pp.countval(lk) == p.profile_size, "profile_rates_end must match profile_rates in length");
            pp.getarr(lk, last);
            const bool t0 = pp.query(key("profile_time_start"), p.profile_time_start);
            const bool t1 = pp.query(key("profile_time_end"), p.profile_time_end);
            require(t0 && t1 && amrex::Math::isfinite(p.profile_time_start) && amrex::Math::isfinite(p.profile_time_end)
                    && p.profile_time_start < p.profile_time_end && amrex::Math::isfinite(p.profile_time_end-p.profile_time_start),
                    "profile_rates_end requires finite, increasing profile_time_start/end in seconds");
        }
        std::string du = "m", ru = "s^-1";
        pp.query(key("profile_diameter_units"), du);
        pp.query(key("profile_rate_units"), ru);
        Real meters = 1, seconds = 1;
        if (du == "um" || du == "micron") { meters = Real(1.e-6); }
        else if (du == "nm") { meters = Real(1.e-9); }
        else { require(du == "m", "profile_diameter_units: expected m, um, micron or nm"); }
        if (ru == "h^-1" || ru == "h-1" || ru == "1/h") { seconds = 3600; }
        else { require(ru == "s^-1" || ru == "s-1" || ru == "1/s", "profile_rate_units: expected inverse seconds or inverse hours"); }

        auto rate = [&](Real input) {
            require(amrex::Math::isfinite(input) && input >= 0, "profile rates must be finite and nonnegative");
            const Real si = input / seconds;
            require(input == 0 || si > 0, "profile rate underflow after conversion to inverse seconds");
            return si;
        };
        for (int i = 0; i < p.profile_size; ++i) {
            const Real d = diam[i] * meters;
            require(amrex::Math::isfinite(d) && d > 0 && (i == 0 || d > p.profile_diameter[i-1]),
                    "profile_diameters must be finite, positive and strictly increasing in meters");
            p.profile_diameter[i] = d;
            p.profile_rate_start[i] = rate(first[i]);
            p.profile_rate_end[i] = p.profile_time_varying ? rate(last[i]) : p.profile_rate_start[i];
        }
        if (p.settling_correction) {
            amrex::Warning(scope + ".settling_correction is ignored for a measured whole-chamber loss rate");
        }
    }
    if (!remaining.empty()) { amrex::Abort("Unknown or inapplicable wall-loss parameter: " + *remaining.begin()); }
}

std::string signature (const Params& p)
{
    std::ostringstream s;
    s.imbue(std::locale::classic());
    s << std::setprecision(std::numeric_limits<Real>::max_digits10)
      << "sdwl1," << p.enable << ',' << int(p.method) << ',' << int(p.dep_vel_model)
      << ',' << p.deactivate_on_depletion << ',' << p.deactivate_threshold
      << ',' << p.vd_side << ',' << p.vd_top << ',' << p.vd_bottom << ',' << p.vd_floor
      << ',' << p.settling_correction << ',' << p.diff_temperature << ',' << p.diff_viscosity
      << ',' << p.diff_delta << ',' << p.diff_slip_correction << ',' << p.diff_mean_free_path
      << ',' << p.num_bins << ',' << p.bin_rmin << ',' << p.bin_rmax << ',' << p.diag_interval
      << ',' << p.profile_size << ',' << p.profile_time_varying << ',' << p.profile_time_start << ',' << p.profile_time_end;
    for (int i = 0; i < p.profile_size; ++i) {
        s << ',' << p.profile_diameter[i] << ',' << p.profile_rate_start[i] << ',' << p.profile_rate_end[i];
    }
    return s.str();
}

} // namespace sdwl

namespace {

using namespace amrex;
using Groups = GpuArray<Real,3>; // lateral walls, roof, floor

bool wall_bc (ERF_BC bc) noexcept
{
    return bc == ERF_BC::no_slip_wall || bc == ERF_BC::slip_wall || bc == ERF_BC::surface_layer;
}

struct Walls {
    Box domain;
    GpuArray<Real,3> plo, dxi;
    GpuArray<bool,3> low{}, high{};
    Groups area_fraction{}; // measured profile: reporting only

    explicit Walls (const Geometry& geom, const GpuArray<ERF_BC,6>& bc)
        : domain(geom.Domain()), plo(geom.ProbLoArray()), dxi(geom.InvCellSizeArray())
    {
        for (int d = 0; d < 3; ++d) {
            low[d] = !geom.isPeriodic(d) && wall_bc(bc[Orientation(d,Orientation::low)]);
            high[d] = !geom.isPeriodic(d) && wall_bc(bc[Orientation(d,Orientation::high)]);
        }
        if (!any()) { return; }
        // Dividing all wall areas by the same chamber volume leaves their
        // fractions unchanged: A_x/V = 1/Lx, etc., including one-cell domains.
        area_fraction[0] = (int(low[0])+int(high[0]))/geom.ProbLength(0)
                         + (int(low[1])+int(high[1]))/geom.ProbLength(1);
        area_fraction[1] = int(high[2])/geom.ProbLength(2);
        area_fraction[2] = int(low[2])/geom.ProbLength(2);
        const Real sum = area_fraction[0] + area_fraction[1] + area_fraction[2];
        for (int g = 0; g < 3; ++g) { area_fraction[g] /= sum; }
    }

    bool any () const noexcept { return low[0] || low[1] || low[2] || high[0] || high[1] || high[2]; }

    //! Cell-local A_b/V_cell, summed over every touched face (not else-if).
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    Groups area_per_volume (const IntVect& iv) const noexcept
    {
        Groups a{};
        for (int d = 0; d < 3; ++d) {
            if (low[d] && iv[d] == domain.smallEnd(d))  { a[d == 2 ? 2 : 0] += dxi[d]; }
            if (high[d] && iv[d] == domain.bigEnd(d)) { a[d == 2 ? 1 : 0] += dxi[d]; }
        }
        return a;
    }
};

// Capture only these tile pointers, geometry and a device pointer to Params;
// none of the profile arrays or unrelated species/ice data are kernel arguments.
struct TileLoss {
    SuperDropletPC::ParticleType* particles;
    int* active;
    ParticleReal* mult;
    const ParticleReal *radius, *mass, *vterm;
    const sdwl::Params* params;
    Real* totals;
    Walls walls;
    Real dt, midpoint;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void apply (int i, const RandomEngine* rng) const noexcept
    {
        auto& particle = particles[i];
        // Negative activity tags are pending recycled placements, not suspended SDs.
        if (particle.id() <= 0 || active[i] <= 0) { return; }
        using Calc = std::common_type_t<ParticleReal, Real>;
        const Calc old = mult[i];
        if (old < 0) { return; } // pending inflow multiplicity
        const auto& p = *params;
        const Real r = Real(radius[i]);
        Groups share;
        Real k;
        if (p.dep_vel_model == sdwl::SDWL_MeasuredRate) {
            k = sdwl::sd_wall_loss_measured_rate(r, midpoint, p);
            share = walls.area_fraction;
        } else {
            auto a = walls.area_per_volume(getParticleCell(particle, walls.plo, walls.dxi, walls.domain));
            if (a[0] == 0 && a[1] == 0 && a[2] == 0) { return; }
            // Stokes-Einstein is independent of face; evaluate it only once.
            const Real side = sdwl::sd_vd_turbulent(r, p, p.vd_side);
            const bool diff = p.dep_vel_model == sdwl::SDWL_Diffusion;
            const Real top = diff ? side : amrex::max(p.vd_top, p.vd_floor);
            const Real bottom = diff ? side : amrex::max(p.vd_bottom, p.vd_floor);
            const Real settle = p.settling_correction ? Real(std::abs(vterm[i])) : Real(0);
            share = {{a[0]*side, a[1]*amrex::max(top-settle, Real(0)), a[2]*(bottom+settle)}};
            k = share[0]+share[1]+share[2];
            if (k > 0) { for (int g = 0; g < 3; ++g) { share[g] /= k; } }
        }
        if (!(k > 0)) { return; }

        const Calc optical_depth = Calc(k)*Calc(dt);
        Calc next;
        if (rng != nullptr) {
            const Calc prob = sdwl::sd_wall_loss_probability(Calc(k), Calc(dt));
            next = old - sdwl::sd_wall_loss_stochastic_removed(old, prob, *rng);
        } else if (optical_depth < Calc(0.5)) {
            next = old - old * sdwl::sd_wall_loss_probability(Calc(k), Calc(dt));
        } else {
            // 1-P rounds to zero for large k*dt even with representable survivors.
            next = old * std::exp(-optical_depth);
        }
        ParticleReal stored = ParticleReal(next);
        if (p.deactivate_on_depletion && (stored == 0 || Calc(stored) < Calc(p.deactivate_threshold))) {
            stored = 0;
            active[i] = 0;
            particle.id().make_invalid(); // terminal loss, never recycled
        }
        // Tally the change actually representable in the particle SoA. In
        // particular a sub-ULP update must not deposit phantom number or mass.
        if (stored == mult[i]) { return; }
        mult[i] = stored;
        const Calc removed = old - Calc(stored);
        for (int g = 0; g < 3; ++g) {
            if (share[g] == 0) { continue; }
            const Calc number = removed * share[g];
            HostDevice::Atomic::Add(totals+g, Real(number));
            HostDevice::Atomic::Add(totals+3+g, Real(number * mass[i]));
        }
        HostDevice::Atomic::Add(totals+6+sdwl::sd_wall_loss_bin_index(r, p), Real(removed));
    }
};

} // namespace

void SuperDropletPC::WallLoss (int lev, amrex::Real dt, amrex::Real time, int iter, const BCTypeArr& bc)
{
    if (!m_wall_loss.enable) { return; }
    BL_PROFILE("SuperDropletPC::WallLoss()");
    using namespace amrex;
    using namespace SDPCDefn;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lev == 0 && m_gdb->maxLevel() == 0 && !m_variable_dz && Geom(0).IsCartesian(),
                                     "SD wall loss requires a single uniform Cartesian level");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Math::isfinite(dt) && dt >= 0 && Math::isfinite(time), "Invalid SD wall-loss time or dt");
    const Walls walls(Geom(lev), bc);
    const int n = 6 + m_wall_loss.num_bins;
    Vector<Real> global(n, Real(0));
    if (walls.any() && dt > 0) {
        // A host-to-device copy (not a fill kernel capturing the whole profile),
        // completed before ParIter can rotate through multiple GPU streams.
        const Gpu::DeviceVector<sdwl::Params> config{m_wall_loss};
        Gpu::DeviceVector<Real> totals(n, Real(0));
        const bool stochastic = m_wall_loss.method == sdwl::SDWL_Stochastic;
        const Real midpoint = time + Real(0.5)*dt;
        constexpr int roff = SuperDropletsRealIdx::ncomps, ioff = SuperDropletsIntIdx::ncomps;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
            auto& soa = pti.GetStructOfArrays();
            const TileLoss tile{pti.GetArrayOfStructs()().data(),
                soa.GetIntData(ioff + SuperDropletsIntIdxSoA_RT::active).data(),
                soa.GetRealData(roff + SuperDropletsRealIdxSoA_RT::multiplicity).data(),
                soa.GetRealData(roff + SuperDropletsRealIdxSoA_RT::radius).data(),
                soa.GetRealData(SuperDropletsRealIdx::mass).data(),
                soa.GetRealData(roff + SuperDropletsRealIdxSoA_RT::term_vel).data(),
                config.dataPtr(), totals.dataPtr(), walls, dt, midpoint};
            const int np = pti.GetArrayOfStructs().numParticles();
            if (stochastic) {
                ParallelForRNG(np, [=] AMREX_GPU_DEVICE (int i, const RandomEngine& rng) noexcept { tile.apply(i, &rng); });
            } else {
                ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept { tile.apply(i, nullptr); });
            }
        }
        Gpu::synchronize(); // all tile streams before copy and config destruction
        Gpu::copy(Gpu::deviceToHost, totals.begin(), totals.end(), global.begin());
    }
    ParallelDescriptor::ReduceRealSum(global.dataPtr(), n);
    m_wall_loss_history.record(m_name, m_wall_loss_signature, m_wall_loss.diag_interval, iter+1, time, dt, global);
}
#endif // ERF_USE_PARTICLES
