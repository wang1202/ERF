#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>

#include <gtest/gtest.h>

#include <limits>

#include "ERF_SolverUtils.H"
#include "ERF_Constants.H"

namespace {

#if (AMREX_SPACEDIM == 3)
amrex::Geometry
make_geometry (amrex::Array<int, AMREX_SPACEDIM> const& periodic)
{
    const amrex::Box domain(amrex::IntVect(10,20,30), amrex::IntVect(17,27,35));
    const amrex::RealBox real_box({-5.0, 2.0, 10.0}, {-1.0, 10.0, 13.0});
    return amrex::Geometry(domain, real_box, amrex::CoordSys::cartesian, periodic);
}

amrex::Array<std::string, 2*AMREX_SPACEDIM>
wall_boundary_names ()
{
    amrex::Array<std::string, 2*AMREX_SPACEDIM> names;
    names.fill("NoSlipWall");
    return names;
}
#endif

} // namespace

#if (AMREX_SPACEDIM == 3)
TEST(TerrainMLMGRegion, LocalGeometryPreservesIndexOriginSpacingAndPhysicalOffset)
{
    const auto global = make_geometry({0, 1, 0});
    const amrex::Box region(amrex::IntVect(12,22,32), amrex::IntVect(15,25,34));
    const auto local = terrain_mlmg_region_geometry(global, region);

    EXPECT_EQ(local.Domain(), region);
    EXPECT_DOUBLE_EQ(local.CellSize(0), 0.5);
    EXPECT_DOUBLE_EQ(local.CellSize(1), 1.0);
    EXPECT_DOUBLE_EQ(local.CellSize(2), 0.5);
    EXPECT_DOUBLE_EQ(local.ProbLo(0), -4.0);
    EXPECT_DOUBLE_EQ(local.ProbLo(1), 4.0);
    EXPECT_DOUBLE_EQ(local.ProbLo(2), 11.0);
    EXPECT_DOUBLE_EQ(local.ProbHi(0), -2.0);
    EXPECT_DOUBLE_EQ(local.ProbHi(1), 8.0);
    EXPECT_DOUBLE_EQ(local.ProbHi(2), 12.5);
    EXPECT_FALSE(local.isPeriodic(0));
    EXPECT_FALSE(local.isPeriodic(1));
    EXPECT_FALSE(local.isPeriodic(2));
}

TEST(TerrainMLMGRegion, FullDomainRetainsPeriodicAxes)
{
    const auto global = make_geometry({1, 1, 0});
    const auto local = terrain_mlmg_region_geometry(global, global.Domain());
    EXPECT_TRUE(local.isPeriodic(0));
    EXPECT_TRUE(local.isPeriodic(1));
    EXPECT_FALSE(local.isPeriodic(2));
    EXPECT_TRUE(terrain_mlmg_periodic_region_error(
        0, 0, global, global.Domain()).empty());
}

TEST(TerrainMLMGRegion, SplitBoxesCoverDisconnectedRectangularRegions)
{
    const auto global = make_geometry({0, 0, 0});
    const amrex::Box first(amrex::IntVect(11,21,31), amrex::IntVect(14,24,34));
    const amrex::Box second(amrex::IntVect(15,21,31), amrex::IntVect(16,24,34));
    amrex::BoxList boxes;
    boxes.push_back(first);
    boxes.push_back(second);
    const amrex::BoxArray ba(std::move(boxes));
    const amrex::Box region(amrex::IntVect(11,21,31), amrex::IntVect(16,24,34));

    EXPECT_TRUE(terrain_mlmg_region_error(1, 0, global.Domain(), region, ba).empty());
    EXPECT_TRUE(terrain_mlmg_periodic_region_error(1, 0, global, region).empty());

    const amrex::Box other_region(amrex::IntVect(11,25,31), amrex::IntVect(16,26,34));
    const amrex::BoxArray other_ba(other_region);
    EXPECT_TRUE(terrain_mlmg_region_error(2, 1, global.Domain(), other_region, other_ba).empty());
}

TEST(TerrainMLMGRegion, RejectsConnectedLShapeWithLevelAndCellCounts)
{
    const auto global = make_geometry({0, 0, 0});
    amrex::BoxList boxes;
    boxes.push_back(amrex::Box(amrex::IntVect(11,21,31), amrex::IntVect(12,22,34)));
    boxes.push_back(amrex::Box(amrex::IntVect(13,21,31), amrex::IntVect(14,22,34)));
    boxes.push_back(amrex::Box(amrex::IntVect(11,23,31), amrex::IntVect(12,24,34)));
    const amrex::BoxArray ba(std::move(boxes));
    const amrex::Box region(amrex::IntVect(11,21,31), amrex::IntVect(14,24,34));

    const auto error = terrain_mlmg_region_error(2, 3, global.Domain(), region, ba);
    EXPECT_NE(error.find("lev=2, isub=3"), std::string::npos);
    EXPECT_NE(error.find("region=[11:14,21:24,31:34]"), std::string::npos);
    EXPECT_NE(error.find("actual_cells="), std::string::npos);
    EXPECT_NE(error.find("expected_cells="), std::string::npos);
    EXPECT_NE(error.find("L-shaped"), std::string::npos);
}

TEST(TerrainMLMGRegion, RejectsPartialPeriodicSeamButAcceptsInteriorPatch)
{
    const auto global = make_geometry({1, 1, 0});
    const amrex::Box seam_patch(amrex::IntVect(10,22,31), amrex::IntVect(13,25,34));
    const auto error = terrain_mlmg_periodic_region_error(1, 0, global, seam_patch);
    EXPECT_NE(error.find("partial periodic-seam refinement"), std::string::npos);
    EXPECT_NE(error.find("lev=1, isub=0, dir=0"), std::string::npos);

    const amrex::Box interior_patch(amrex::IntVect(12,22,31), amrex::IntVect(15,25,34));
    EXPECT_TRUE(terrain_mlmg_periodic_region_error(
        2, 1, global, interior_patch).empty());
    const auto local = terrain_mlmg_region_geometry(global, interior_patch);
    EXPECT_FALSE(local.isPeriodic(0));
    EXPECT_FALSE(local.isPeriodic(1));
}

TEST(TerrainMLMGRegion, RejectsUnvalidatedLevelAndDisconnectedRegionScopes)
{
    EXPECT_TRUE(terrain_mlmg_scope_error(1, 1, 1).empty());
    const auto three_level_error = terrain_mlmg_scope_error(2, 2, 1);
    EXPECT_NE(three_level_error.find("at most two ERF levels"), std::string::npos);
    EXPECT_NE(three_level_error.find("max_level=2"), std::string::npos);
    const auto disconnected_error = terrain_mlmg_scope_error(1, 1, 2);
    EXPECT_NE(disconnected_error.find("one connected rectangular solve region"), std::string::npos);
    EXPECT_NE(disconnected_error.find("level 1 has 2 regions"), std::string::npos);
}

TEST(TerrainMLMGCompatibility, UsesPrecisionAppropriateDivergenceThreshold)
{
    const amrex::Real tolerance = terrain_mlmg_compatibility_tolerance();
    EXPECT_EQ(tolerance, sizeof(amrex::Real) == sizeof(float)
                            ? amrex::Real(1.e-4) : amrex::Real(1.e-6));
    EXPECT_FALSE(terrain_mlmg_compatibility_mean_exceeds(tolerance));
    EXPECT_TRUE(terrain_mlmg_compatibility_mean_exceeds(tolerance * amrex::Real(1.01)));
    EXPECT_TRUE(terrain_mlmg_compatibility_mean_exceeds(
        std::numeric_limits<amrex::Real>::quiet_NaN()));
    EXPECT_TRUE(terrain_mlmg_compatibility_mean_exceeds(
        std::numeric_limits<amrex::Real>::infinity()));
}

TEST(TerrainProjectionBC, PhysicalAndArtificialFacesUseRegionalConditions)
{
    const auto global = make_geometry({0, 1, 0});
    auto names = wall_boundary_names();
    names[amrex::Orientation(0, amrex::Orientation::low)] = "Outflow";
    names[amrex::Orientation(0, amrex::Orientation::high)] = "Outflow";
    names[amrex::Orientation(2, amrex::Orientation::low)] = "Outflow";
    names[amrex::Orientation(2, amrex::Orientation::high)] = "Outflow";

    const auto interior = terrain_mlmg_region_bcs(
        global, amrex::Box(amrex::IntVect(12,22,31), amrex::IntVect(15,25,34)), names);
    EXPECT_EQ(interior.lo[0], amrex::LinOpBCType::Neumann);
    EXPECT_EQ(interior.hi[0], amrex::LinOpBCType::Neumann);
    EXPECT_EQ(interior.lo[1], amrex::LinOpBCType::Neumann);
    EXPECT_EQ(interior.hi[1], amrex::LinOpBCType::Neumann);
    EXPECT_EQ(interior.lo[2], amrex::LinOpBCType::Neumann);
    EXPECT_EQ(interior.hi[2], amrex::LinOpBCType::Neumann);
    EXPECT_TRUE(interior.is_singular);

    const auto mixed = terrain_mlmg_region_bcs(
        global, amrex::Box(amrex::IntVect(10,22,30), amrex::IntVect(15,25,35)), names);
    EXPECT_EQ(mixed.lo[0], amrex::LinOpBCType::Dirichlet);
    EXPECT_EQ(mixed.hi[0], amrex::LinOpBCType::Neumann);
    EXPECT_EQ(mixed.lo[2], amrex::LinOpBCType::Neumann);
    EXPECT_EQ(mixed.hi[2], amrex::LinOpBCType::Dirichlet);
    EXPECT_TRUE(mixed.is_singular == false);
}

TEST(TerrainProjectionBC, PhysicalLateralWallAndTerrainLowerZAreNeumann)
{
    const auto global = make_geometry({0, 0, 0});
    auto names = wall_boundary_names();
    names[amrex::Orientation(0, amrex::Orientation::high)] = "Outflow";
    names[amrex::Orientation(2, amrex::Orientation::low)] = "Outflow";
    const auto bcs = terrain_mlmg_region_bcs(global, global.Domain(), names);

    EXPECT_EQ(bcs.lo[0], amrex::LinOpBCType::Neumann);
    EXPECT_EQ(bcs.hi[0], amrex::LinOpBCType::Dirichlet);
    EXPECT_EQ(bcs.lo[2], amrex::LinOpBCType::Neumann);
    EXPECT_TRUE(bcs.is_singular == false);
}

TEST(TerrainProjectionBC, WholeDomainPeriodicBCsRemainPeriodic)
{
    const auto global = make_geometry({1, 1, 0});
    auto names = wall_boundary_names();
    const auto bcs = terrain_mlmg_region_bcs(global, global.Domain(), names);
    EXPECT_EQ(bcs.lo[0], amrex::LinOpBCType::Periodic);
    EXPECT_EQ(bcs.hi[0], amrex::LinOpBCType::Periodic);
    EXPECT_EQ(bcs.lo[1], amrex::LinOpBCType::Periodic);
    EXPECT_EQ(bcs.hi[1], amrex::LinOpBCType::Periodic);
    EXPECT_TRUE(bcs.is_singular);
}
#endif

TEST(TerrainMLMGConfiguration, DefaultAndSupportedSolverNames)
{
    EXPECT_STREQ(terrain_poisson_solver_default(), "gmres_fft");
    EXPECT_TRUE(terrain_poisson_solver_choice_error("gmres_fft").empty());
    EXPECT_TRUE(terrain_poisson_solver_choice_error("mlmg").empty());
    EXPECT_NE(terrain_poisson_solver_choice_error("unknown").find("gmres_fft or mlmg"),
              std::string::npos);
}

TEST(TerrainMLMGConfiguration, AcceptsSupportedConfiguration)
{
    EXPECT_TRUE(terrain_mlmg_configuration_error(true, true, true, false).empty());
}

TEST(TerrainMLMGConfiguration, RejectsUnsupportedConfigurationsWithInputNames)
{
    EXPECT_NE(terrain_mlmg_configuration_error(false, true, true, false).find("AMREX_SPACEDIM=3"), std::string::npos);
    EXPECT_NE(terrain_mlmg_configuration_error(true, false, true, false).find("VariableDz mesh"), std::string::npos);
    EXPECT_NE(terrain_mlmg_configuration_error(true, true, false, false).find("erf.terrain_type=StaticFittedMesh"), std::string::npos);
    EXPECT_NE(terrain_mlmg_configuration_error(true, true, true, true).find("erf.use_real_bcs=true"), std::string::npos);
}

TEST(TerrainProjectionGhost, RejectsPrecisionSpecificFillSentinelsAndNonfiniteValues)
{
    EXPECT_TRUE(terrain_projection_ghost_value_invalid(bogus_large_value,
                                                        bogus_large_value));
    EXPECT_TRUE(terrain_projection_ghost_value_invalid(
        bogus_large_value * amrex::Real(0.75), bogus_large_value));
    EXPECT_TRUE(terrain_projection_ghost_value_invalid(
        std::numeric_limits<amrex::Real>::max(), bogus_large_value));
    EXPECT_TRUE(terrain_projection_ghost_value_invalid(
        std::numeric_limits<amrex::Real>::quiet_NaN(), bogus_large_value));
    EXPECT_TRUE(terrain_projection_ghost_value_invalid(
        std::numeric_limits<amrex::Real>::infinity(), bogus_large_value));
    EXPECT_FALSE(terrain_projection_ghost_value_invalid(amrex::Real(123.5),
                                                         bogus_large_value));
}

TEST(TerrainProjectionGhost, CatchesSingleSentinelMagnitudeInEitherMeshPrecision)
{
    const amrex::Real single_precision_sentinel = amrex::Real(1.e18);
    EXPECT_TRUE(terrain_projection_ghost_value_invalid(
        single_precision_sentinel, single_precision_sentinel));
    EXPECT_TRUE(terrain_projection_ghost_value_invalid(
        single_precision_sentinel * amrex::Real(0.75), single_precision_sentinel));
}
