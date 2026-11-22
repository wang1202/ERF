#include <ERF_SDWallLoss.H>

#include <AMReX_OpenMP.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Random.H>
#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

namespace {

using amrex::Real;
constexpr Real eps = std::numeric_limits<Real>::epsilon();
// Relative allowance for a few products, divisions and logarithms in mesh precision.
constexpr Real rtol = 16 * eps;

TEST(SDWallLoss, StokesEinsteinWithAndWithoutSlip)
{
    sdwl::Params p;
    p.dep_vel_model = sdwl::SDWL_Diffusion;
    p.diff_temperature = 305;
    p.diff_viscosity = Real(2.e-5);
    p.diff_delta = Real(4.e-4);
    p.diff_mean_free_path = Real(8.e-8);
    const Real radius = 2 * p.diff_mean_free_path; // Kn = 1/2
    const double stokes = 1.380649e-23 * p.diff_temperature
        / (6.0 * std::acos(-1.0) * p.diff_viscosity * radius);
    const double cunningham = 1.0 + 0.5 * (1.257 + 0.4 * std::exp(-2.2));

    for (bool slip : {false, true}) {
        SCOPED_TRACE(slip);
        p.diff_slip_correction = slip;
        const double diffusivity = stokes * (slip ? cunningham : 1.0);
        const double velocity = diffusivity / p.diff_delta;
        EXPECT_NEAR(sdwl::sd_brownian_diffusivity(radius, p), diffusivity, rtol * diffusivity);
        EXPECT_NEAR(sdwl::sd_vd_turbulent(radius, p, Real(1)), velocity, rtol * velocity);
    }
}

TEST(SDWallLoss, VelocityFloorBoundsBothFaceModels)
{
    sdwl::Params p;
    const Real r = Real(1.e-6);
    p.vd_floor = Real(0.125);
    p.dep_vel_model = sdwl::SDWL_Constant;
    EXPECT_EQ(sdwl::sd_vd_turbulent(r, p, Real(0.0625)), p.vd_floor);
    EXPECT_EQ(sdwl::sd_vd_turbulent(r, p, Real(0.25)), Real(0.25));
    p.dep_vel_model = sdwl::SDWL_Diffusion;
    EXPECT_EQ(sdwl::sd_vd_turbulent(r, p, Real(0.25)), p.vd_floor);
}

TEST(SDWallLoss, FaceCoefficientAndSmallStepProbability)
{
    const Real k = sdwl::sd_wall_loss_coefficient(Real(0.5), Real(3), Real(12));
    ASSERT_EQ(k, Real(0.125)); // v_d A/V, in inverse seconds
    EXPECT_EQ(sdwl::sd_wall_loss_probability(k, Real(0)), Real(0));
    EXPECT_EQ(sdwl::sd_wall_loss_probability(Real(0), Real(8)), Real(0));
    EXPECT_NEAR(sdwl::sd_wall_loss_probability(k, Real(8) * std::log(Real(2))), Real(0.5), rtol);

    // At k*dt = eps/8, exp(-k*dt) rounds to 1, but the loss is k*dt to roundoff.
    const Real tiny = eps / 8;
    EXPECT_NEAR(sdwl::sd_wall_loss_probability(k, eps), tiny, rtol * tiny);
}

TEST(SDWallLoss, MeasuredRateInterpolatesLogDiameterAndTime)
{
    sdwl::Params p;
    const Real d = Real(1.e-6);
    p.profile_size = 3;
    p.profile_diameter[0] = d;  p.profile_diameter[1] = 4*d;  p.profile_diameter[2] = 16*d;
    p.profile_rate_start[0] = 2;  p.profile_rate_start[1] = 10;  p.profile_rate_start[2] = 4;
    p.profile_rate_end[0] = 14;   p.profile_rate_end[1] = 6;     p.profile_rate_end[2] = 18;
    p.profile_time_varying = true;
    p.profile_time_start = 10;
    p.profile_time_end = 30;

    // The public lookup takes a wet radius. Diameters 2d and 8d are the
    // geometric midpoints; the *rates* there are arithmetic means.
    auto rate = [&](Real radius, Real time) { return sdwl::sd_wall_loss_measured_rate(radius, time, p); };
    EXPECT_EQ(rate(d/4, 10), Real(2));     // diameter below the first knot
    EXPECT_EQ(rate(2*d, 10), Real(10));    // exact interior diameter knot
    EXPECT_EQ(rate(16*d, 30), Real(18));   // diameter above the last knot
    EXPECT_NEAR(rate(d, 10), 6, rtol * 6);
    EXPECT_NEAR(rate(4*d, 5), 7, rtol * 7);       // before the time window
    EXPECT_NEAR(rate(4*d, 35), 12, rtol * 12);    // after the time window
    EXPECT_NEAR(rate(4*d, 15), 8.25, rtol * 8.25); // 1/4 of the way from 7 to 12
    p.profile_time_varying = false;
    EXPECT_NEAR(rate(4*d, 35), 7, rtol * 7);
}

TEST(SDWallLoss, NearbyDiameterKnotsRemainDistinct)
{
    // At a micron, adjacent diameters can have identical rounded logarithms.
    // Three equally spaced representable diameters give a log midpoint to O(eps).
    const Real lo = std::ldexp(Real(1), -20);
    const Real mid = std::nextafter(lo, Real(1));
    const Real hi = std::nextafter(mid, Real(1));
    sdwl::Params p;
    p.profile_size = 2;
    p.profile_diameter[0] = lo;
    p.profile_diameter[1] = hi;
    p.profile_rate_start[0] = 2;
    p.profile_rate_start[1] = 6;
    EXPECT_NEAR(sdwl::sd_wall_loss_profile_rate(mid, p, false), 4, rtol * 4);
}

// ParmParse objects share a global table. Only claim an empty test prefix and
// erase its fully qualified entries even if a fatal assertion exits the test.
constexpr char profile_scope[] = "gtest_sd_wall_loss_profile";
struct ScopedProfileInputs {
    amrex::ParmParse pp{profile_scope};
    ~ScopedProfileInputs () {
        amrex::ParmParse root;
        for (const auto& key : amrex::ParmParse::getEntries(profile_scope)) { root.remove(key); }
    }
};

TEST(SDWallLoss, ParserConvertsMeasuredProfilesToSI)
{
    struct Units { const char* length; double d; const char* rate; double k; };
    // The same two-knot, two-time profile in default SI units and lab units.
    for (const Units u : {Units{nullptr, 1.e-6, nullptr, 0.01},
                          Units{"um", 1, "h^-1", 36}, Units{"nm", 1000, "1/h", 36}}) {
        SCOPED_TRACE(u.length ? u.length : "default SI");
        ASSERT_TRUE(amrex::ParmParse::getEntries(profile_scope).empty());
        ScopedProfileInputs inputs;
        auto& pp = inputs.pp;
        pp.add("wall_loss.enable", true);
        pp.add("wall_loss.deposition_velocity_model", "measured_rate");
        if (u.length) { pp.add("wall_loss.profile_diameter_units", u.length); }
        if (u.rate) { pp.add("wall_loss.profile_rate_units", u.rate); }
        pp.addarr("wall_loss.profile_diameters", std::vector<double>{u.d, 4*u.d});
        pp.addarr("wall_loss.profile_rates", std::vector<double>{u.k, 4*u.k});
        pp.addarr("wall_loss.profile_rates_end", std::vector<double>{2*u.k, 8*u.k});
        pp.add("wall_loss.profile_time_start", 10);
        pp.add("wall_loss.profile_time_end", 30);

        sdwl::Params p;
        sdwl::read_params(pp, p);
        ASSERT_TRUE(p.enable);
        ASSERT_EQ(p.dep_vel_model, sdwl::SDWL_MeasuredRate);
        ASSERT_EQ(p.profile_size, 2);
        ASSERT_TRUE(p.profile_time_varying);
        EXPECT_EQ(p.profile_time_start, Real(10));
        EXPECT_EQ(p.profile_time_end, Real(30));
        for (int i = 0; i < p.profile_size; ++i) {
            SCOPED_TRACE(i);
            const int scale = (i == 0) ? 1 : 4;
            EXPECT_NEAR(p.profile_diameter[i], scale * 1.e-6, rtol * scale * 1.e-6);
            EXPECT_NEAR(p.profile_rate_start[i], scale * 0.01, rtol * scale * 0.01);
            EXPECT_NEAR(p.profile_rate_end[i], scale * 0.02, rtol * scale * 0.02);
        }
    }
}

TEST(SDWallLoss, RadiusBinsAreLogSpacedAndClamped)
{
    sdwl::Params p;
    p.bin_rmin = Real(1.e-8);
    p.bin_rmax = Real(1.e-4);
    p.num_bins = 4; // one decade per bin
    auto bin = [&](Real r) { return sdwl::sd_wall_loss_bin_index(r, p); };
    EXPECT_EQ(bin(p.bin_rmin/2), 0);
    EXPECT_EQ(bin(p.bin_rmin), 0);
    Real radius = 3 * p.bin_rmin;
    for (int i = 0; i < 4; ++i, radius *= 10) {
        EXPECT_EQ(bin(radius), i) << "radius = " << radius;
    }
    EXPECT_EQ(bin(std::nextafter(p.bin_rmax, p.bin_rmin)), 3);
    EXPECT_EQ(bin(p.bin_rmax), 3);
    EXPECT_EQ(bin(2*p.bin_rmax), 3);
    p.num_bins = 1;
    EXPECT_EQ(bin(Real(1.e-6)), 0);
    EXPECT_EQ(bin(2*p.bin_rmax), 0);
}

TEST(SDWallLoss, StochasticRemovalHonorsEndpoints)
{
    const auto engine = amrex::getInvalidRandomEngine(); // host draws do not use a device engine
    const Real xi = Real(7.5);
    EXPECT_EQ(sdwl::sd_wall_loss_stochastic_removed(xi, Real(0), engine), Real(0));
    EXPECT_EQ(sdwl::sd_wall_loss_stochastic_removed(xi, Real(1), engine), xi);
    EXPECT_EQ(sdwl::sd_wall_loss_stochastic_removed(Real(0), Real(0.5), engine), Real(0));
}

TEST(SDWallLoss, StochasticRemovalRespectsBoundsAndMoments)
{
    // Draw only on the host. Leave the AMReX stream as found for shuffled/repeated suites.
    struct ScopedHostRng {
        std::stringstream state;
        ScopedHostRng () { amrex::SaveRandomState(state); }
        ~ScopedHostRng () { amrex::RestoreRandomState(state, amrex::OpenMP::get_max_threads(), 0); }
    } restore_rng;
    const auto engine = amrex::getInvalidRandomEngine();

    // Fractional multiplicities near depletion make upper clipping observable
    // in both the Poisson (xi*p < 30) and Gaussian branches.
    for (const Real xi : {Real(0.5), Real(31.5)}) {
        for (int i = 0; i < 128; ++i) {
            const Real removed = sdwl::sd_wall_loss_stochastic_removed(xi, Real(0.99), engine);
            ASSERT_GE(removed, 0) << "xi = " << xi;
            ASSERT_LE(removed, xi) << "xi = " << xi;
        }
    }

    // Far from the clipping limits: Poisson(16) has variance 16; the
    // large-mean approximation at xi=128, p=1/2 has binomial variance 32.
    struct Sample { Real xi, prob; double mean, variance; };
    for (const Sample s : {Sample{64, Real(0.25), 16, 16}, Sample{128, Real(0.5), 64, 32}}) {
        SCOPED_TRACE(s.mean);
        constexpr int n = 4096;
        double sum = 0, squared_deviations = 0;
        for (int i = 0; i < n; ++i) {
            const Real removed = sdwl::sd_wall_loss_stochastic_removed(s.xi, s.prob, engine);
            ASSERT_GE(removed, 0) << "draw " << i;
            ASSERT_LE(removed, s.xi) << "draw " << i;
            sum += removed;
            const double deviation = removed - s.mean;
            squared_deviations += deviation * deviation;
        }
        EXPECT_NEAR(sum/n, s.mean, 8 * std::sqrt(s.variance/n)); // eight standard errors
        // Relative standard error of this second moment is at most sqrt((2+1/16)/n).
        EXPECT_NEAR(squared_deviations/n, s.variance, 0.20 * s.variance);
    }
}

} // namespace
