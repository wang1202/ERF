#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <AMReX_Arena.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <gtest/gtest.h>

#include "Diagnostics/ERF_CloudChamberWallFluxDiag.H"
#include "../ERF_GTestTempDir.H"

namespace
{

using amrex::Real;
using namespace erf_wallflux_diag;

FaceValues face_scaled (Real value)
{
    // Distinct per-face values expose a face permutation.
    FaceValues result{};
    for (int f = 0; f < NumFaces; ++f) { result[f] = value * Real(f + 1); }
    return result;
}

// One accepted step: every required channel at both stages, then commit.
void take_step (WallFluxWindow& window, Real dt, Real stage0, Real stage1)
{
    for (int c = 0; c < NumChannels; ++c) {
        if (!window.requires_channel(c)) { continue; }
        ASSERT_TRUE(window.capture(c, 0, face_scaled(stage0)).empty());
        ASSERT_TRUE(window.capture(c, 1, face_scaled(stage1)).empty());
    }
    ASSERT_TRUE(window.commit(dt).empty());
}

FaceValues unit_area ()
{
    FaceValues area{};
    area.fill(Real(1.0));
    return area;
}

} // namespace

TEST(CloudChamberWallFluxDiag, StageWeightsFollowTheAnelasticIntegrators)
{
    EXPECT_EQ(stage_weight(StageWeighting::Trapezoidal, 0), Real(0.5));
    EXPECT_EQ(stage_weight(StageWeighting::Trapezoidal, 1), Real(0.5));
    EXPECT_EQ(stage_weight(StageWeighting::MidPoint, 0), Real(0.0));
    EXPECT_EQ(stage_weight(StageWeighting::MidPoint, 1), Real(1.0));
}

// dt = 1 and 3 with stage pairs (2, 6) and (10, 14): the trapezoidal mean is
// (1*(2+6)/2 + 3*(10+14)/2)/4 = 10, the MidPoint mean (1*6 + 3*14)/4 = 12.  A
// snapshot, an unweighted step mean or interval * nominal dt would all give
// other values.
TEST(CloudChamberWallFluxDiag, UnequalStepsWeightByActualDuration)
{
    for (const auto weighting : {StageWeighting::Trapezoidal, StageWeighting::MidPoint}) {
        WallFluxWindow window(weighting, true);
        take_step(window, Real(1.0), Real(2.0), Real(6.0));
        take_step(window, Real(3.0), Real(10.0), Real(14.0));
        EXPECT_EQ(window.committed_steps(), 2);
        EXPECT_EQ(window.elapsed(), Real(4.0));
        const auto means = interval_means(window.window_integral(), unit_area(), window.elapsed());
        const Real expected = (weighting == StageWeighting::Trapezoidal) ? Real(10.0) : Real(12.0);
        for (int c = 0; c < NumChannels; ++c) {
            for (int f = 0; f < NumFaces; ++f) {
                EXPECT_NEAR(means[c][f], expected * Real(f + 1), Real(1.0e-12) * expected * Real(f + 1))
                    << "channel " << c << " face " << f;
            }
        }
    }
}

TEST(CloudChamberWallFluxDiag, MeanDividesByFaceAreaAndElapsedTime)
{
    WallFluxWindow window(StageWeighting::MidPoint, false);
    take_step(window, Real(0.5), Real(0.0), Real(8.0));
    FaceValues area{};
    for (int f = 0; f < NumFaces; ++f) { area[f] = Real(2.0); }
    const auto means = interval_means(window.window_integral(), area, window.elapsed());
    EXPECT_NEAR(means[Stress][0], Real(4.0), Real(1.0e-14));
    // Not transported here: never captured, so exactly zero.
    EXPECT_EQ(means[Vapor][3], Real(0.0));
}

TEST(CloudChamberWallFluxDiag, MissingOrDuplicateStageCapturesAreRejected)
{
    WallFluxWindow window(StageWeighting::Trapezoidal, true);
    EXPECT_TRUE(window.capture(Stress, 0, face_scaled(1.0)).empty());
    EXPECT_FALSE(window.capture(Stress, 0, face_scaled(1.0)).empty());   // duplicate
    EXPECT_FALSE(window.capture(Stress, 2, face_scaled(1.0)).empty());   // no such stage
    EXPECT_FALSE(window.commit(Real(1.0)).empty());                      // incomplete step
    EXPECT_EQ(window.committed_steps(), 0);
    EXPECT_EQ(window.elapsed(), Real(0.0));

    WallFluxWindow dry(StageWeighting::Trapezoidal, false);
    EXPECT_FALSE(dry.capture(Vapor, 1, face_scaled(1.0)).empty());       // not transported

    WallFluxWindow midpoint(StageWeighting::MidPoint, false);
    // MidPoint stage 0 carries zero weight but must still be seen once.
    EXPECT_TRUE(midpoint.capture(Stress, 1, face_scaled(1.0)).empty());
    EXPECT_TRUE(midpoint.capture(Heat, 1, face_scaled(1.0)).empty());
    EXPECT_FALSE(midpoint.commit(Real(1.0)).empty());
}

TEST(CloudChamberWallFluxDiag, RowLayoutMatchesTheLegacyWriter)
{
    const std::string header = header_line();
    EXPECT_EQ(header.rfind("# time[s]", 0), 0u);
    EXPECT_EQ(header.size(), 15u + 18u * NumChannels * NumFaces);
    // Right-aligned 18-character columns after a 15-character time column.
    const std::string xlo_columns = std::string(7, ' ') + "xlo_tau[Pa]" +
                                    std::string(7, ' ') + "xlo_H[W/m2]" +
                                    std::string(4, ' ') + "xlo_E[kg/m2/s]";
    EXPECT_EQ(header.substr(0, 15 + 3 * 18), "# time[s]      " + xlo_columns);
    EXPECT_NE(header.find("zhi_E[kg/m2/s]"), std::string::npos);

    ChannelFaceValues means{};
    for (int c = 0; c < NumChannels; ++c) {
        for (int f = 0; f < NumFaces; ++f) { means[c][f] = Real(10 * f + c) - Real(20.0); }
    }
    const std::string row = format_row(0.25, means);
    EXPECT_EQ(row.size(), 15u + 18u * NumChannels * NumFaces);
    std::istringstream tokens(row);
    std::vector<std::string> fields;
    for (std::string token; tokens >> token;) { fields.push_back(token); }
    ASSERT_EQ(static_cast<int>(fields.size()), NumColumns);
    EXPECT_EQ(fields[0], "2.50000000e-01");
    // FACE_tau FACE_H FACE_E, face by face: xlo_H is column 2, yhi_E column 12.
    EXPECT_EQ(fields[2], "-1.90000000e+01");
    EXPECT_EQ(fields[12], "1.20000000e+01");
}

// interval 2, five steps: rows at the ends of steps 2 and 4 only, each the
// mean of its own window (reset after every row); the trailing one-step
// window writes nothing.
TEST(CloudChamberWallFluxDiag, CompleteWindowsOnlyAndResetAfterEachRow)
{
    const std::filesystem::path dir = erf_gtest_temp_path("erf_wallflux_diag");
    if (amrex::ParallelDescriptor::IOProcessor()) { std::filesystem::create_directories(dir); }
    amrex::ParallelDescriptor::Barrier();
    const std::string path = (dir / "wallflux_lev0.txt").string();

    const amrex::Box domain(amrex::IntVect(0), amrex::IntVect(3));
    const amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});
    const amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
    const amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

    CloudChamberWallFluxDiag diag(2, StageWeighting::MidPoint, false, geom, path);
    diag.require_fresh_output();
    const int nranks = amrex::ParallelDescriptor::NProcs();
    double time = 0.0;
    for (int step = 1; step <= 5; ++step) {
        // Each rank holds an equal share of the step's value step*10.
        const Real value = Real(10 * step) / Real(nranks);
        for (int c : {Stress, Heat}) {
            diag.capture_stage(c, 0, FaceValues{});
            FaceValues stage1{};
            stage1.fill(value);
            diag.capture_stage(c, 1, stage1);
        }
        const Real dt = Real(step);
        diag.commit_step(dt);
        time += dt;
        diag.report(step, time);
    }

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream in(path);
        std::string line;
        std::vector<std::string> lines;
        while (std::getline(in, line)) { lines.push_back(line); }
        ASSERT_EQ(lines.size(), 3u);
        EXPECT_EQ(lines[0], header_line());
        // Window 1: steps 1-2, (1*10 + 2*20)/3; window 2: steps 3-4, (3*30 + 4*40)/7.
        const auto field = [] (const std::string& text, int n) {
            std::istringstream tokens(text);
            std::string token;
            for (int i = 0; i <= n; ++i) { tokens >> token; }
            return std::stod(token);
        };
        EXPECT_NEAR(field(lines[1], 0), 3.0, 1.0e-12);
        EXPECT_NEAR(field(lines[1], 1), 50.0 / 3.0, 1.0e-6);
        EXPECT_NEAR(field(lines[2], 0), 10.0, 1.0e-12);
        EXPECT_NEAR(field(lines[2], 1), 250.0 / 7.0, 1.0e-6);
        EXPECT_EQ(field(lines[2], 3), 0.0);  // xlo_E: not transported
        std::filesystem::remove_all(dir);
    }
    amrex::ParallelDescriptor::Barrier();
}

// A face whose patch tractions alternate in sign has a zero mean vector but a
// unit mean magnitude: the magnitude must be taken per patch, before any
// average.  Wall-edge nodes, shared with a perpendicular wall and applied by
// no RHS, must not enter the reconstruction.
TEST(CloudChamberWallFluxDiag, StressMagnitudeIsLocalAndSkipsEdgeNodes)
{
    const int n = 4;
    const amrex::Box domain(amrex::IntVect(0), amrex::IntVect(n - 1));
    const amrex::GpuArray<Real, AMREX_SPACEDIM> dx{AMREX_D_DECL(Real(0.25), Real(0.25), Real(0.25))};

    const amrex::Box box12 = amrex::convert(domain, amrex::IntVect(AMREX_D_DECL(1, 1, 0)));
    const amrex::Box box13 = amrex::convert(domain, amrex::IntVect(AMREX_D_DECL(1, 0, 1)));
    const amrex::Box box23 = amrex::convert(domain, amrex::IntVect(AMREX_D_DECL(0, 1, 1)));
    amrex::FArrayBox tau12(box12, 1, amrex::The_Pinned_Arena());
    amrex::FArrayBox tau13(box13, 1, amrex::The_Pinned_Arena());
    amrex::FArrayBox tau23(box23, 1, amrex::The_Pinned_Arena());
    tau12.setVal<amrex::RunOn::Host>(Real(0.0));
    tau23.setVal<amrex::RunOn::Host>(Real(0.0));
    tau13.setVal<amrex::RunOn::Host>(Real(0.0));

    // zlo wall (node k = 0): x-traction +1 on even y rows and -1 on odd ones;
    // the x-edge nodes i = 0 and i = n hold a large value that must be ignored.
    const auto t13 = tau13.array();
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i <= n; ++i) {
            const bool edge = (i == 0 || i == n);
            t13(i, j, 0) = edge ? Real(1.0e6) : ((j % 2 == 0) ? Real(1.0) : Real(-1.0));
        }
    }

    erf_wallflux_diag::FaceReducer reducer;
    erf_wallflux_diag::add_wall_stress(reducer, domain, domain,
                                       tau12.const_array(), tau13.const_array(),
                                       tau23.const_array(), dx);
    const FaceValues sums = reducer.value();

    // zlo: 16 patches of area 1/16 with |tau| = 1, so the area integral is 1,
    // although the signed traction sums to zero.
    EXPECT_NEAR(sums[4], Real(1.0), Real(1.0e-12));
    // The xlo and xhi walls see tau13 only on their k = 0 edge, which is excluded.
    EXPECT_EQ(sums[0], Real(0.0));
    EXPECT_EQ(sums[1], Real(0.0));
    EXPECT_EQ(sums[5], Real(0.0));
}
