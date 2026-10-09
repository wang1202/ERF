# Terrain MLMG coarse-fine AMR validation

> **Historical evidence only:** All run results below predate the current rebased source and are from dirty executables. They do not validate the current branch. The authoritative current-source test matrix, blockers, and merge verdict are in [terrain_mlmg_merge_readiness.md](terrain_mlmg_merge_readiness.md).

Date: 2026-10-09

## Status

**Experimental; general multi-level AMR support is not yet claimed.** The
two-level partial-vertical TerrainHill case completes with finite fields on
one and two MPI ranks, including x/y/z box splits. A full-height two-level
control completed one timestep with both MLMG and GMRES+FFT. Its final state
fields are close but not identical, and later fine-level GMRES projections
left substantially more divergence than MLMG. The partial-height GMRES case
still reaches the known full-domain-z-metric FFT assertion. Restart/regrid,
disconnected regions, three levels, boundary-flux comparison between solvers,
CPU and no-FFT builds remain open; the final one-level default regression passed.

Branch: `terrain_optimization`. Reviewed base: `c82278faa770f13c43bebef45371d44e1962831c`.
AMReX submodule: `5e1ef8467df5f067be4d673c6efe45e1ffb3c9b1`.
Configuration: 3D CUDA, H100, MPI and FFT enabled. The current production
executable was rebuilt after the diagnostic probes were removed.

## First invalid Omega input

The pre-fix scanner examined the x/y face momenta read by `OmegaFromW`
and reported the first invalid value:

- ERF level 1, source box 0: `((16,8,0) (48,23,25))`
- Target w-face: `(16,18,26)`
- Read component and source face: `xmom(17,18,26)`
- Value: `1e150`; outside the valid face box, classified as a coarse-fine ghost.

This is the top vertical ghost of the partially refined level-1 patch. It is
not a physical-boundary ghost. Before the fix, the resulting level-1 divergence
was `1.942566943e148` in L_inf. The temporary exact-stencil scanner and
log are preserved at
`MyBuild/ERF_TerrainAMRValidation/hill-18968275.out` (job 18968275).

The cause was that `ERFFillPatcher::FillSet()` / `InterpFace()`
only interpolates the valid fine face box, while `FillBoundary()`
exchanges fine-fine data and cannot populate coarse-fine ghosts. The projection
now builds coarse x/y momentum in the same `rho0`-weighted representation
as the valid fine faces, using ERF's `ConvertForProjection()`, fills
same-level coarse data, and calls face-centered `FillPatchTwoLevels()`.
It copies only vertical ghost strips inside the fine level's global face domain.
This targeted fill runs immediately before `OmegaFromW()` and again
after the projection correction, after `FillBoundary()` and before
`WFromOmega()`. Interpolated values are checked for nonfinite and
sentinel-scale magnitudes; failures report level, box, component, face index,
value, and coarse-fine/fine-fine classification. Physical-boundary extrapolation
is left to the existing ERF path.

## Numerical results

| Run | Layout and backend | Result |
|---|---|---|
| Pre-fix, job 18968275 | 2 levels, partial vertical refinement, MLMG | Exact first invalid x-momentum ghost above; level-1 divergence reached 1.94e148. |
| Pre-fix, job 18967978 | Same case, GMRES+FFT | Same 1.94e148 level-1 divergence before the legacy FFT assertion. This establishes the ghost defect was shared, not caused by MLMG. |
| Post-fix, job 18968692 | Same case, GMRES+FFT | Finite level-1 pre-solve divergence (0.34296); legacy backend then aborts at AMReX FFT assertion `dz.size() == domain.length(2)` because the valid fine region is vertically partial. |
| Post-fix, job 18968920 | Same case, MLMG, one MPI rank, one H100 | Completed one full timestep; exit 0. |
| Post-fix, job 18968965 | Same case, x/y/z multi-box layout, two MPI ranks, one H100 | Completed one full timestep; exit 0. This used the Cray Shasta MPI plugin and a 1 GiB per-rank AMReX arena initial allocation. |
| Post-fix, job 18969191 | Full-height two-level patch, one rank, MLMG and GMRES+FFT | Both completed one timestep; initial post-solve divergence was below 1e-8 on both levels. Final field differences are quantified below; during the timestep GMRES fine-level post-solve divergence reached 2.15e-3 while MLMG remained below 1e-8. |

For job 18968920, the level-1 solve region was
`[16:47,8:23,0:25]`, with a partial vertical extent. Initial
projection divergence fell from `0.3429587494` to
`1.078199744e-9` (L_inf); level 0 fell from
`0.145697201` to `4.433911233e-9`. Both levels logged
`operator AMR level 0`, confirming the per-region one-level operator
index convention; source calls `compFlux(0, ...)`. The initial level-1
MLMG solve took 6 iterations with final residual `1.06152197e-9`.
The time step completed at `t=0.5`. Post-correction ghost validation
also completed without an invalid-value diagnostic.

Job 18968965 used `amr.max_grid_size_x=8`,
`amr.max_grid_size_y=8`, and `amr.max_grid_size_z=8`. The
log reports two MPI processes, level 0 and 1 MLMG solves, and post-projection
level-1 L_inf values from `1.08e-9` initially through
`5.71e-8` during the step; the final logged level-1 solve was
`4.75e-9`. This exercises multiple boxes and vertical box splits on
two ranks. The first two-rank attempt (job 18968960) was invalid because the
default Slurm launcher started independent singleton jobs that contended for
the GPU; it is excluded from the result.

The partial-height legacy run is not a solver-parity result: GMRES+FFT
receives finite input after the shared ghost fill, then its FFT asserts because
the shortened vertical region has no full-domain z metric vector. The default
remains `gmres_fft`; solver selection and the legacy backend are unchanged.

For the common-supported full-height control in job 18969191, both backends
started from the same input and completed one timestep. Initial-projection
L_inf divergence changed from 0.145697201 to 4.44e-9 (MLMG) and 1.70e-10
(GMRES) on level 0, and from 0.3429587494 to 9.76e-9 (MLMG) and 1.78e-9
(GMRES) on level 1. Later GMRES fine-level post-solve divergence reached
2.15e-3; MLMG's largest logged fine-level value was 9.76e-9. This difference
needs investigation and prevents a parity claim.

The final `plt00001` fields were compared on matching level grids. Density was
bitwise identical. The other fields had these max absolute and relative RMS
differences (RMS normalized by the GMRES field RMS):

| Level | Field | Max absolute difference | Relative RMS |
|---|---|---:|---:|
| 0 | x_velocity | 3.05e-3 | 2.31e-5 |
| 0 | y_velocity | 4.62e-4 | 3.08e-4 |
| 0 | z_velocity | 2.12e-3 | 1.55e-3 |
| 0 | theta | 2.73e-2 | 4.94e-6 |
| 1 | x_velocity | 5.07e-3 | 5.66e-5 |
| 1 | y_velocity | 9.22e-4 | 4.15e-4 |
| 1 | z_velocity | 3.03e-3 | 1.61e-3 |
| 1 | theta | 7.28e-2 | 3.12e-5 |

The run did not measure GMRES boundary-normal correction fluxes, so the
boundary-flux comparison is still open.

## Automated coverage added

- `ERF_GTestTerrainPoissonBC.cpp` contains geometry, rectangular
  region, disconnected region, irregular L-shape rejection, periodic-seam
  rejection, effective BC, singularity, and configuration tests. The focused
  test target passed all 11 terrain tests in H100 job 18967961.
- `TerrainMLMG_TwoLevelPartialVertical` is registered in
  `Tests/CTestList.cmake`. Its checker requires both ERF levels to
  reach MLMG, requires a completed timestep, rejects invalid ghost diagnostics,
  and bounds each initial post-projection L_inf divergence below 1e-6.
- The CMake integration runner passed on H100 in job 18969119. Its checker
  verified one completed step, finite coarse-fine ghosts, and post-projection
  divergence below 1e-6 on both levels.
- CTest discovery from the existing MyBuild tree fails because the unit-test
  directory expects a generated MPI CTest file that is absent. The integration
  runner was run directly through CMake in the H100 job.
- The same run sampled each effective Neumann face. Level 0 terrain zlo/zhi
  had 512 samples each; the level-1 rectangular region had 416 samples on
  each x face, 832 samples on each y face, and 512 samples on each z face.
  Every sampled maximum absolute normal correction flux was zero.

## Required test matrix

| ID | Coverage | Status |
|---|---|---|
| T0 | One-level default GMRES+FFT regression | **PASS**, job 18969221, exit 0. One level completed a full timestep at t=0.5 on H100 using the default solver selection. |
| T1 | One-level MLMG regression | Historical one-level runs passed; not rerun on the final source tree. |
| T2 | Two levels, partial rectangular fine region, MLMG, one rank | **PASS**, job 18968920. |
| T3 | Same case, x/y multi-box, at least two MPI ranks | **PASS**, job 18968965. |
| T4 | Same case with z-split boxes and at least two MPI ranks | **PASS**, combined with T3 in job 18968965. |
| T5 | Two disconnected rectangular fine regions | Not run. |
| T6 | Three levels, partial regions on levels 1 and 2 | Not run. |
| T7 | Restart/regrid with MLMG | Not run. |
| T8 | Same partial-vertical case with GMRES+FFT | Backend reaches level 1 with finite input after the fix, then fails at the known full-height FFT dz assertion; no parity comparison is possible for this layout. |
| T9 | CPU MLMG build and run | Not run. |
| T10 | No-FFT MLMG build and run | Not run. |
| T11 | Boundary variants and measured normal correction-flux norms | **Partial PASS**: job 18969119 sampled terrain zlo and all six effective level-1 Neumann faces; every MLMG maximum absolute normal correction flux was 0. Job 18969191 compared state fields on a full-height layout, but did not measure GMRES boundary-normal fluxes. The zlo=outflow input variant remains unrun. |
| T12 | Executable-level connected-irregular-region rejection | Unit helper rejects an L-shaped region; executable-level negative run not performed. |
| T13 | Executable-level unsupported solver/configuration diagnostics | Helper-level configuration tests exist; executable-level negative run not performed. |

## Remaining correctness checks

1. Measure GMRES boundary-normal correction fluxes on the full-height
   comparison layout, investigate its fine-level post-solve divergence up to
   2.15e-3, and run the zlo=outflow boundary variant.
2. Run end-to-end negative executable tests for an irregular region and an
   incompatible solver configuration.
3. Exercise disconnected regions, three levels, and restart/regrid.
4. Run CPU and FFT-disabled builds.
5. The one-level final-source regression passed in job 18969221.

Do not describe general multi-level AMR support as validated until the remaining
layout and lifecycle tests succeed. The current evidence demonstrates one
two-level partial-vertical case on one and two MPI ranks with x/y/z box splits.
No current performance claim is made.

Job 18969119 also checked the normal correction flux directly after the solve.
The sampled faces were the physical terrain zlo/zhi faces on level 0 and the
six effective Neumann faces bounding the level-1 rectangular solve region.
Every face had a nonzero sample count and `max_abs=0`. This confirms the
instrumented homogeneous-Neumann normal flux for this fixture; the outflow
input variant and other layouts remain untested.

## Artifacts and commands

Key shared logs:
- `MyBuild/ERF_TerrainAMRValidation/hill-18968920.out`
- `MyBuild/ERF_TerrainAMRValidation/multibox-18968965.out`
- `MyBuild/ERF_TerrainAMRValidation/gmres-18968692.out`
- `MyBuild/ERF_TerrainAMRValidation/parity-18969191.out`
- `MyBuild/ERF_TerrainAMRValidation/fullheight_state_diff-18969191.txt`
- `MyBuild/ERF_TerrainAMRValidation/onelevel-18969221.out`
- `MyBuild/ERF_TerrainAMRValidation/hill-18968275.out`

Production rebuild commands were the generated Makefile targets for
`ERF_PoissonSolve.cpp.o`, `erf_srclib`, and
`Exec/erf_exec`; no CMake cache reconfiguration was run. Run commands
and exact runtime options are recorded in the associated Slurm scripts under
`MyBuild/ERF_TerrainAMRValidation`. Build hash printed by ERF:
`040db7a7d280-dirty`. The task-instruction markdown and `MyBuild`
artifacts are not intended for the commit.
