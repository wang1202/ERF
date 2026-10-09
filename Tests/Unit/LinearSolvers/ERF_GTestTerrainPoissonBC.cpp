#include <AMReX_Geometry.H>
#include <AMReX_FFT.H>

#include <gtest/gtest.h>

#include "ERF_SolverUtils.H"

namespace {

#if (AMREX_SPACEDIM == 3)
amrex::Geometry
make_geometry (amrex::Array<int, AMREX_SPACEDIM> const& periodic)
{
    const amrex::Box domain(amrex::IntVect(0), amrex::IntVect(7));
    const amrex::RealBox real_box({0.0, 0.0, 0.0}, {8.0, 8.0, 8.0});
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
TEST(TerrainProjectionBC, LowerZOutflowRetainsLegacyNeumannAndOtherFaces)
{
    const amrex::Geometry geom = make_geometry({0, 1, 0});
    auto names = wall_boundary_names();
    names[amrex::Orientation(2, amrex::Orientation::low)] = "Outflow";
    names[amrex::Orientation(2, amrex::Orientation::high)] = "Outflow";
    names[amrex::Orientation(0, amrex::Orientation::low)] = "Outflow";

    const auto generic_lo = get_lo_projection_bc(geom, names);
    const auto lo = get_lo_terrain_projection_bc(geom, names);
    const auto hi = get_hi_projection_bc(geom, names);

    EXPECT_EQ(generic_lo[2], amrex::LinOpBCType::Dirichlet);
    EXPECT_EQ(lo[2], amrex::LinOpBCType::Neumann);
    EXPECT_EQ(hi[2], amrex::LinOpBCType::Dirichlet);
    EXPECT_EQ(lo[0], amrex::LinOpBCType::Dirichlet);
    EXPECT_EQ(lo[1], amrex::LinOpBCType::Periodic);
}

TEST(TerrainProjectionBC, PeriodicZIsNotOverridden)
{
    const amrex::Geometry geom = make_geometry({1, 1, 1});
    auto names = wall_boundary_names();
    names[amrex::Orientation(2, amrex::Orientation::low)] = "Outflow";

    const auto lo = get_lo_terrain_projection_bc(geom, names);
    EXPECT_EQ(lo[2], amrex::LinOpBCType::Periodic);
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
    EXPECT_TRUE(terrain_mlmg_configuration_error(
        0, true, true, true, false, true).empty());
}

TEST(TerrainMLMGConfiguration, RejectsUnsupportedConfigurationsWithInputNames)
{
    EXPECT_NE(terrain_mlmg_configuration_error(
        0, false, true, true, false, true).find("AMREX_SPACEDIM=3"), std::string::npos);
    EXPECT_NE(terrain_mlmg_configuration_error(
        0, true, false, true, false, true).find("erf.mesh_type=VariableDz"), std::string::npos);
    EXPECT_NE(terrain_mlmg_configuration_error(
        0, true, true, false, false, true).find("erf.terrain_type=StaticFittedMesh"), std::string::npos);
    EXPECT_NE(terrain_mlmg_configuration_error(
        0, true, true, true, true, true).find("erf.use_real_bcs=true"), std::string::npos);
    EXPECT_NE(terrain_mlmg_configuration_error(
        1, true, true, true, false, true).find("amr.max_level=0"), std::string::npos);
    EXPECT_NE(terrain_mlmg_configuration_error(
        0, true, true, true, false, false).find("complete level-0 domain"), std::string::npos);
}
