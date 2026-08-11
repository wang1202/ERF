# ERF Instrument Simulator Port — CLI Agent Implementation Plan

## Goal

Port the instrument simulator implementation from `wang1202/ERF_instrument` into the current public repository `wang1202/ERF`.

Create a new branch:

```bash
Test2
```

The implementation must support:

1. Ceilometer — Vaisala CL-31, 905 nm
2. Doppler lidar — Halo StreamLine XR+, 1550 nm
3. Microwave radiometer (MWR) — Radiometrics MP-3000A, K+V band

It must include the moving-column behavior used by the supplied `move_start_time`, `move_stop_time`, `move_speed_x`, and `move_speed_y` inputs.

## Repositories and legacy commits

Target:

```text
https://github.com/wang1202/ERF
```

Legacy source:

```text
https://github.com/wang1202/ERF_instrument
```

Relevant commits:

```text
52152ffa29fc53a7294c36410536c7c59b90d125
Add instrument simulators (ceilometer, Doppler lidar, MWR)
```

```text
6b9e71beeec7e067c3edd48ce25a536de8ffb2e4
Finish instrument simulator integration
```

Do not blindly cherry-pick the entire second commit; it includes unrelated changes. Port only simulator-specific files/hunks.

The current target default branch is `development`.

## 1. Create branch

```bash
git clone https://github.com/wang1202/ERF.git
cd ERF
git checkout development
git pull --ff-only
git checkout -b Test2

git remote add instrument https://github.com/wang1202/ERF_instrument.git
git fetch instrument development
```

Do not merge the old branch wholesale.

## 2. Port the completed simulator header

```bash
git show instrument/development:Source/IO/ERF_InstrumentSimulators.H \
  > Source/IO/ERF_InstrumentSimulators.H
```

The completed version should contain:

- `InstrumentSimUtil::SiteLocation`
- `InstrumentSimUtil::MoveSchedule`
- `InstrumentSimUtil::site_location_at_time`
- GPU-safe column extraction
- simple geometric-optics backscatter
- size-resolved extinction support
- `CeilometerSimulator`
- `DopplerLidarSimulator`
- `MWRSimulator`
- moving physical `(x,y)` and sampled `(i,j)` in outputs

The moving schedule arrays must have equal lengths. Intervals must not overlap.

For each site, physical position starts at the original cell center and accumulates piecewise-constant displacement. Convert the instantaneous `(x,y)` back to sampled `(i,j)` and clamp the indices to the domain.

## 3. Integrate in `Source/ERF.H`

Add:

```cpp
#include <ERF_InstrumentSimulators.H>
```

near the other IO/sampling headers.

Add members near the line/plane samplers:

```cpp
std::unique_ptr<CeilometerSimulator>   m_ceilometer_sim;
std::unique_ptr<DopplerLidarSimulator> m_doppler_lidar_sim;
std::unique_ptr<MWRSimulator>          m_mwr_sim;
```

Do not disturb unrelated current-branch members.

## 4. Initialize in `Source/ERF.cpp`

Near existing line/plane sampler setup, after `geom[0]` is valid:

```cpp
{
    amrex::ParmParse pp_ceil("erf.ceilometer");
    bool do_ceil = false;
    pp_ceil.query("do_sim", do_ceil);
    if (do_ceil) {
        m_ceilometer_sim = std::make_unique<CeilometerSimulator>();
        m_ceilometer_sim->init(geom[0]);
    }
}

{
    amrex::ParmParse pp_dl("erf.doppler_lidar");
    bool do_dl = false;
    pp_dl.query("do_sim", do_dl);
    if (do_dl) {
        m_doppler_lidar_sim = std::make_unique<DopplerLidarSimulator>();
        m_doppler_lidar_sim->init(geom[0]);
    }
}

{
    amrex::ParmParse pp_mwr("erf.mwr");
    bool do_mwr = false;
    pp_mwr.query("do_sim", do_mwr);
    if (do_mwr) {
        m_mwr_sim = std::make_unique<MWRSimulator>();
        m_mwr_sim->init(geom[0]);
    }
}
```

## 5. Add write path in `ERF::post_timestep`

Port the simulator-specific `post_timestep` integration from legacy commit `52152ffa...`, then apply the moving-location corrections from `6b9e71b...`.

Use level-0 state:

```cpp
vars_new[0][Vars::cons]
vars_new[0][Vars::xvel]
vars_new[0][Vars::yvel]
vars_new[0][Vars::zvel]
geom[0]
```

Required gating:

```cpp
const bool ceil_write =
    m_ceilometer_sim &&
    m_ceilometer_sim->is_write_step(nstep);

const bool dl_write =
    m_doppler_lidar_sim &&
    m_doppler_lidar_sim->is_write_step(nstep);

const bool need_sigma = ceil_write || dl_write;
```

When SDM is active, only compute particle-resolved extinction on actual simulator write steps.

Use current moving site locations:

```cpp
m_ceilometer_sim->ILoc(s, t_new[0], geom[0])
m_ceilometer_sim->JLoc(s, t_new[0], geom[0])
```

and likewise for Doppler lidar.

Do not revert to fixed `ILoc(s)` / `JLoc(s)` for SDM sampling.

Call each simulator's `write(...)` using the current API in the copied header. Adapt only for current ERF API/type changes.

## 6. Preserve microphysics behavior

Bulk schemes should use standard conserved-state moisture slots:

```text
RhoQ1_comp — vapor
RhoQ2_comp — cloud liquid
RhoQ3_comp — cloud ice
```

Bulk defaults:

```text
r_eff_default     = 1.0e-5 m
r_eff_ice_default = 3.0e-5 m
```

Preserve SBM size-resolved extinction if the current branch still exposes the required bin interface.

Preserve SDM particle-resolved extinction and compute it only on output steps.

## 7. Port unit test

```bash
mkdir -p Tests/Unit/IO

git show instrument/development:Tests/Unit/IO/ERF_GTestInstrumentSimulators.cpp \
  > Tests/Unit/IO/ERF_GTestInstrumentSimulators.cpp
```

The legacy test checks WSM6-style conserved-state compatibility and validates:

- ceilometer output
- Doppler vertical stare
- Doppler VAD
- Doppler backscatter
- MWR profile
- MWR LWP/PWV column output

## 8. Register unit test

In `Tests/Unit/CMakeLists.txt`, add:

```cmake
${CMAKE_CURRENT_SOURCE_DIR}/IO/ERF_GTestInstrumentSimulators.cpp
```

to the existing unit-test target.

## 9. Add focused moving-site tests

Add deterministic tests for `MoveSchedule` / `ILoc` / `JLoc`:

1. unchanged before first start time;
2. correct displacement inside an interval;
3. frozen displacement after stop;
4. later intervals accumulate;
5. `(x,y)` converts to correct `(i,j)`;
6. out-of-domain motion clamps indices;
7. mismatched schedule lengths assert.

## 10. Acceptance input

```text
# ---- Ceilometer ----
erf.ceilometer.do_sim              = true
erf.ceilometer.i_loc               = 159
erf.ceilometer.j_loc               = 75
erf.ceilometer.move_start_time     = 60. 1200000000.
erf.ceilometer.move_stop_time      = 345600. 4800000000.
erf.ceilometer.move_speed_x        = -6.3 0.
erf.ceilometer.move_speed_y        = -3.1 0.
erf.ceilometer.output_interval     = 1
erf.ceilometer.output_file         = ceil_morrison_
erf.ceilometer.backscatter_model   = simple
erf.ceilometer.r_eff_default       = 1.0e-5
erf.ceilometer.r_eff_ice_default   = 3.0e-5
erf.ceilometer.beta_thresh         = 3.0e-4

# ---- Doppler Lidar ----
erf.doppler_lidar.do_sim               = true
erf.doppler_lidar.i_loc                = 159
erf.doppler_lidar.j_loc                = 75
erf.doppler_lidar.move_start_time      = 60. 1200000000.
erf.doppler_lidar.move_stop_time       = 345600. 4800000000.
erf.doppler_lidar.move_speed_x         = -6.3 0.
erf.doppler_lidar.move_speed_y         = -3.1 0.
erf.doppler_lidar.output_interval      = 1
erf.doppler_lidar.output_file          = dl_morrison_
erf.doppler_lidar.backscatter_model    = simple
erf.doppler_lidar.r_eff_default        = 1.0e-5
erf.doppler_lidar.do_vertical_stare    = true
erf.doppler_lidar.do_vad               = true
erf.doppler_lidar.vad_elevation_angle  = 75.0
erf.doppler_lidar.do_backscatter       = true
erf.doppler_lidar.cnr_model            = normalized

# ---- MWR ----
erf.mwr.do_sim              = true
erf.mwr.i_loc               = 159
erf.mwr.j_loc               = 75
erf.mwr.move_start_time     = 60. 1200000000.
erf.mwr.move_stop_time      = 345600. 4800000000.
erf.mwr.move_speed_x        = -6.3 0.
erf.mwr.move_speed_y        = -3.1 0.
erf.mwr.output_interval     = 1
erf.mwr.output_file         = mwr_morrison_
erf.mwr.do_brightness_temp  = false
```

Do not hard-code the earlier 96x96 example domain. The parser must continue to support one or multiple sites.

## 11. Preserve final output schemas

Ceilometer:

```text
time x y i_loc j_loc height att_backscatter cloud_base_height
```

Doppler stare:

```text
time x y i_loc j_loc height w CNR
```

Doppler VAD:

```text
time x y i_loc j_loc height u v wind_speed wind_dir
```

Doppler backscatter:

```text
time x y i_loc j_loc height beta_att
```

MWR profile:

```text
time x y i_loc j_loc height T qv qc
```

MWR column:

```text
time x y i_loc j_loc LWP PWV
```

MWR brightness-temperature output, when enabled, must also include instantaneous site location.

## 12. Compatibility review against current ERF

Review current versions of:

```text
Source/ERF.H
Source/ERF.cpp
Source/IO/ERF_InstrumentSimulators.H
Source/ERF_IndexDefines.H
Source/Microphysics/*
Source/Particles/*
Tests/Unit/CMakeLists.txt
```

Check for:

- renamed `MoistureType` values;
- particle-container API changes;
- base-state pressure storage changes;
- `Vars::*` indexing changes;
- `post_timestep` signature changes;
- AMReX API changes;
- warnings promoted to errors;
- `M_PI` portability issues.

Prefer the header's `ERF_INSTR_PI` constant rather than introducing new `M_PI` uses.

## 13. Do not port unrelated legacy changes

Do not port unrelated changes from `6b9e71beeec7e067c3edd48ce25a536de8ffb2e4`, especially:

```text
Source/Microphysics/Morrison/ERF_AdvanceMorrison.cpp
Source/Particles/ERFPCUtils.cpp
Source/PhysicsInterfaces/Radiation/*
```

unless a compile failure demonstrates a direct simulator dependency.

## 14. Documentation

Create:

```text
Docs/instrument_simulators/README.md
```

Document:

- purpose of each simulator;
- runtime parameters;
- moving-column semantics;
- output formats;
- backscatter assumptions;
- SBM/SDM behavior;
- MWR TB option;
- known limitations.

The legacy LaTeX design document is a physics reference only; do not make the build depend on its PDF.

## 15. Build and test

Use the current repository's documented build path.

At minimum run the unit tests containing `InstrumentSimulators`.

Typical commands:

```bash
ctest --output-on-failure
```

or:

```bash
<unit-test-executable> --gtest_filter='InstrumentSimulators.*'
```

Also run a normal ERF build representative of Morrison microphysics and, if supported by the build environment, a particles-enabled build for the SDM compile path.

## 16. Smoke-test acceptance criteria

Confirm:

- all three simulators initialize;
- output honors `output_interval`;
- ceilometer output contains finite backscatter;
- Doppler stare/VAD/beta outputs exist;
- MWR profile/column outputs exist;
- rows contain current `(x,y,i,j)`;
- movement follows the configured schedule;
- outputs continue after the first move interval ends;
- `do_sim = false` disables output;
- boundary-reaching motion causes no out-of-range access;
- MPI produces one consistent output stream per site.

## 17. Suggested commits

```text
Add instrument simulator forward operators
Integrate instrument simulators into ERF timestep output
Add instrument simulator unit tests
Document instrument simulator runtime configuration
```

## 18. Final branch checks

```bash
git status
git diff development...HEAD --stat
git diff development...HEAD
```

Then:

```bash
git push -u origin Test2
```

Do not open a PR unless requested.

## 19. Expected changed files

Normally limited to:

```text
Source/IO/ERF_InstrumentSimulators.H
Source/ERF.H
Source/ERF.cpp
Tests/Unit/IO/ERF_GTestInstrumentSimulators.cpp
Tests/Unit/CMakeLists.txt
Docs/instrument_simulators/README.md
```

A separate focused moving-site test file is fine.

Any changes outside this set should be justified.

## 20. Use exact legacy diffs

For the large integration, inspect exact historical changes:

```bash
git show 52152ffa29fc53a7294c36410536c7c59b90d125 -- Source/ERF.cpp Source/ERF.H
git show 6b9e71beeec7e067c3edd48ce25a536de8ffb2e4 -- Source/ERF.cpp Source/IO/ERF_InstrumentSimulators.H
```

For the new header and existing unit test, copying the completed files from `instrument/development` is preferred, followed by compatibility fixes against current `development`.

## Definition of done

Complete when:

- branch `Test2` exists in `wang1202/ERF`;
- all three simulators compile on current ERF;
- supplied runtime input is accepted;
- moving instrument positions work;
- unit tests pass;
- a short Morrison smoke test produces expected outputs;
- no unrelated legacy changes are included.
