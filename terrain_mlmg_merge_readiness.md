# Terrain MLMG merge readiness

Updated: 2026-10-10  
Feature branch: `terrain-mlmg-pr-ready`  
Published feature tip after requested pull: `4237b43010efd66607d31b196055ac0df21ffc39`  
Current upstream `development`: `b364748d49142d5d8b0b2f8fcdb5f0b2a3dac407` (includes terrain-diffusion PR #4231)  
Feature AMReX pin: `5e1ef8467df5f067be4d673c6efe45e1ffb3c9b1`

## Verdict

**NOT READY.** The feature branch has been fast-forwarded to its current published tip and remains unchanged apart from this local test edit and documentation. A separate validation worktree integrates current upstream in local merge commit `0ffee843a60b356ae4917e184189888cb7f6da83`; it does not alter the feature branch. The current-source CPU DOUBLE/FFT executable builds. The latest restart/regrid runner completes its simulations and projection checks, but its continuous/restart state comparison fails at `plt00002` by large velocity differences. The GMRES residual discrepancy also remains unlocalized and lacks the required same-state dual-solver replay.

Do not open for upstream merge review yet. Keep `erf.terrain_poisson_solver=gmres_fft` as the default. No push, issue publication, or PR action was taken.

## Source and integration state

The initial local feature ref was `75b1a3027`, 10 commits ahead and 2 behind its published branch. `git fetch origin terrain-mlmg-pr-ready` found remote tip `4237b43010efd66607d31b196055ac0df21ffc39`; `git pull --ff-only origin terrain-mlmg-pr-ready` fast-forwarded the feature branch to that SHA. A local review-only commit adds the restart runner layout correction and evidence reports, without changing solver C++. Untracked user artifacts `MyBuild/` and `terrain_mlmg_gmres_diagnostic_agent.md` were preserved.

The validation worktree is `/tmp/ERF_TerrainMLMG_UpstreamValidation`, based on the feature tip plus upstream `development`, with merge commit `0ffee843a60b356ae4917e184189888cb7f6da83`. Its AMReX submodule is initialized at the same feature pin. Merge whitespace check passed. A clean CPU DOUBLE/FFT CMake configure and `erf_exec` build completed. The focused PR #4231 CTests `TerrainStress_Ridge_ImplicitMetric`, `TerrainStress_Hill3D`, and `TerrainStress_Ridge_ImplicitMetricAgree` all passed after building the missing `amrex_fextrema` and `amrex_fcompare` helper targets.

The feature CPU DOUBLE/FFT executable was rebuilt from the refreshed `4237b430` source in `MyBuild/PRReady/final-clean-cpu-double-fft` with CMake 3.29.2, Cray GNU 12.2, FFTW 3.3.10.6, 3D, MPI off, GPU off, DOUBLE, Release. `erf_exec` linked and starts; its banner reports ERF `26.09-578-g4237b43010ef-dirty` and AMReX `26.10-34-g5e1ef8467df5`. The aggregate build exited 2 during GoogleTest discovery because CMake could not execute `Tests/Unit/erf_unit_tests` (`permission denied`); no CTest result is claimed from that build.

## GMRES discrepancy

See [terrain_mlmg_gmres_root_cause.md](terrain_mlmg_gmres_root_cause.md) for the quantitative event record and limitations. The full-height case was rerun on refreshed feature source `4237b430` using the CPU DOUBLE/FFT executable, serial launch, max grid size 8, and the full-height level-1 region. At its first failing event GMRES reports `7.2833e-9` while the immediate independent true residual is L2 `9.0853e-2`, Linf `5.3547e-3` at `(30,14,1)`. Artificial low/high-y raw correction-flux maxima are `2.7016e-3` / `3.3125e-3`. Inputs match the newer report hashes. Exact command and log are in the root-cause report and `/tmp/terrain-mlmg-gmres-current/run.log`. This reproduces the symptom on current feature source; it does not isolate the cause or establish whether current upstream on identical data behaves the same.

Existing clean-upstream/feature same-pin field A/B is evidence of no measured feature regression in those configurations. A new run on the feature-plus-current-upstream integration worktree also reproduces the gap at the first bad event: level 1, call 4, `t=0.25`, estimated residual `7.339462237e-9`, independent Linf `1.388556916e-3` at `(30,15,1)`, L2 `1.818949566e-2`; see `/tmp/terrain-mlmg-gmres-upstream4231/run.log`. This run used the same fixture inputs and solver options as the refreshed feature run, but it is not an immutable in-call snapshot or dual-solver replay. There has been no captured preprojection RHS/state, operator repeatability/ghost test, or two-backend replay. Root cause and same-state upstream status remain **unresolved**. Do not fold a GMRES/FFT correction into this MLMG change without that evidence; the separate issue draft is local and unpublished.

## G0–G17 acceptance matrix

Statuses distinguish current-source evidence from historical evidence. Logs and binaries under `MyBuild/` are local untracked artifacts.

| Gate | Status | Evidence and remaining work |
|---|---|---|
| G0 — default GMRES upstream/feature parity, including fluxes | **BLOCKED** | Current feature-plus-upstream integration run reproduces the first-event residual gap (`7.34e-9` estimated vs `1.39e-3` Linf / `1.82e-2` L2 true residual). No identical snapshot replay or correction-flux parity. |
| G1 — one-level MLMG | **NOT RUN on refreshed source** | Historical CPU DOUBLE runs passed divergence checks; repeat on exact `4237b430` final build. |
| G2 — two-level partial-height MLMG | **NOT RUN on refreshed source** | Historical H100 one-rank run passed; no final-tip rerun. |
| G3 — multi-rank partial-height and rank parity | **BLOCKED** | New `RunTerrainMLMGRankParity.cmake` compares cell/face plots, mass samples, and shared-face audit records, but has not been run on the refreshed source. The available fresh CPU build is serial. |
| G4 — full-height two-level GMRES and MLMG | **FAIL / ROOT CAUSE OPEN** | Current refreshed feature source reproduces first-event GMRES reported residual `7.2833e-9` versus true Linf `5.3547e-3`, L2 `9.0853e-2`; separate MLMG trajectory is not a same-state comparison. No same-state replay or post-#4231 control yet. |
| G5 — disconnected fine regions | **HISTORICAL PASS (scope rejection)** | Earlier executable guard test; rerun on refreshed source remains required. |
| G6 — three active levels | **HISTORICAL PASS (scope rejection)** | Earlier executable guard test; rerun on refreshed source remains required. |
| G7 — changed-layout restart/regrid | **PASS, partial current evidence** | Current CPU runner detected level-1 remake and passed BoxArray and every projection check before reaching plot comparison. Full restart test fails later at G8. |
| G8 — continuous/restart cell and face state at every shared time | **FAIL** | Current fresh CPU run: all legs completed, projection checks passed, `plt00001` passed, but `plt00002` failed. Maximum absolute/relative errors include level-0 x velocity `5.0903e-2` / `7.6016e-3`, level-1 x velocity `1.4892e-1` / `2.0212e-2`, and level-1 z velocity `2.3305e-1` / `7.5765e-2`. Restoring the original 8/8/8 continuous grid setting reproduces the same mismatch; the grid-size discrepancy alone is not the cause. Runner is currently corrected to 16/16/8 to match the restart leg. Logs: `/tmp/terrain-mlmg-regrid-current/terrain_mlmg_regrid_restart/{checkpoint,restart,continuous}/simulation.log`. |
| G9 — physical/artificial boundary variants | **BLOCKED** | Lower-z Outflow, lateral Outflow/Open/Dirichlet, and mixed boundary variants not run. |
| G10 — CPU DOUBLE/FFT | **PARTIAL** | Current feature executable configured, built, linked, and starts; aggregate build failed during GTest discovery with permission denied, and the restart regression fails G8. |
| G11 — CPU SINGLE/no-FFT | **NOT RUN on refreshed source** | Historical tests only. |
| G12 — no-FFT MLMG and GMRES diagnostic | **NOT RUN on refreshed source** | Historical tests only. |
| G13 — CUDA/H100, MPI | **NOT RUN on refreshed source** | Historical H100 direct-run evidence is from prior source. |
| G14 — executable negative matrix, including valid `use_real_bcs=true` fixture | **BLOCKED** | Historical guards cover several invalid cases; no current-source matrix or valid real-boundary fixture result. |
| G15 — nonunity map factors and operator/flux consistency | **BLOCKED** | No manufactured varying-map-factor operator/flux test or verified support restriction. |
| G16 — rank parity, face ownership, conservation | **BLOCKED** | New checker/test exists, but no rank run. No measured global conservation budget including source, boundary, and AMR accounting. |
| G17 — upstream regressions after PR #4231 | **PARTIAL PASS** | Isolated merge build passed. Focused CTests `TerrainStress_Ridge_ImplicitMetric`, `TerrainStress_Hill3D`, and `TerrainStress_Ridge_ImplicitMetricAgree` passed. The integrated `TerrainMLMG_TwoLevelRestartRegridsLevel1` still fails G8 at `plt00002` with level-1 z-velocity abs/rel errors `0.233068` / `0.0757703`. These are applicable focused regressions, not the full upstream suite. |

## Restart failure details

The new runner first had a continuous leg set to `8,8,8` while the restart leg requested `16,16,8`, contrary to the comment that the continuous run keeps the post-regrid tiling. The test was corrected locally to `16,16,8`, and both settings were run with the refreshed CPU DOUBLE/FFT executable. Both runs reproduce the same `plt00002` field mismatch. Logs show both trajectories end coarse steps at times `0.5`, `1.0`, and `1.5`; restart projection calls and continuous projection calls nevertheless differ in their event/stage sequences after the restart point. This is a lifecycle/trajectory failure requiring diagnosis, not a passing comparison. The runner stops at `plt00002`, so later shared times and face plot comparisons were not completed.

## Dependency, support, and performance

The current upstream ERF AMReX pin is `1de18774af5d8a368fa2a398f048e205a2c9f84f`; the feature needs `5e1ef8467df5f067be4d673c6efe45e1ffb3c9b1`, which contains `AMReX_MLTerrainPoisson.H`. The dependency strategy remains an explicit review decision. Advertised scope remains 3D StaticFittedMesh/VariableDz, ERF levels 0 and 1, and one connected rectangular solve region per level; this is not composite-AMR Poisson and does not claim `erf.use_real_bcs=true` support. `gmres_fft` remains the default.

Historical matched one-level benchmarks reported CPU timestep sum 18.0% lower for MLMG and H100 timestep sum 20.6% higher, based on five repetitions. They are not repeated multilevel measurements and support no universal speed claim. They do not validate refreshed source performance.

## Next actions

1. Diagnose the restart mismatch from identical time-aligned output and stage state; retain G8 FAIL until every shared cell/face time passes.
2. Capture the exact first failing GMRES preprojection context and replay both backends on independent copies; exercise operator repeatability, boundary ghosts, and flux/operator consistency.
3. Run current-source G0–G17 gates that are feasible, including 2-rank parity, boundary variants, map-factor handling, and applicable post-#4231 CTests. The three focused #4231 CTests passed; the two-rank gate remains unavailable in the MPI-off build.
4. Decide the AMReX dependency approach and update the PR description to match verified scope.

**Recommendation: NOT READY.** No maintainer disposition exists for the legacy GMRES issue, G8 currently fails, and several numerical and platform gates lack current-source evidence.
# Follow-up validation (2026-10-10)

## G8 continuous/restart control

The previously reported 0.233 m/s level-1 vertical-velocity difference came from comparing different regrid lifecycles. In the controlled no-regrid test, continuous and restart runs with the same BoxArray matched bit-for-bit at all common output times. In the changed-layout test, a restart from 8/8/8 to 16/16/8 (same physical times, refinement coverage, timestep settings, and plot schedules) matched an uninterrupted control that changed its BoxArray at the same time: all fields in `plt`, `pltU`, `pltV`, and `pltW` at t=0.5, 1.0, and 1.5 were bitwise identical (zero absolute and relative difference). The earlier continuous reference started at the 16/16/8 layout and therefore omitted the same `RemakeLevel` operation.

Field-level diagnostics on the uninterrupted runtime-reblock control found the first valid-cell change in `zvel` during `FillPatchFineLevel` in `RemakeLevel` at t=0.5 (Linf 0.2471481141, L2 1.336395862); valid `cons`, `xvel`, `yvel`, and `base_state` were unchanged at that point. The subsequent `project_initial_velocity` reduced divergence Linf from 0.0261532 to 1.265e-8. This is a common remake operation in both equivalent paths, not a continuous/restart discrepancy. No production correction is indicated. The checked-in restart regression still needs to be changed to compare against this lifecycle-equivalent runtime-reblock control; that regression repair was not completed in this round. Raw evidence: `/tmp/g8-runtime-E/run.log`, `/tmp/g8-runtime-E/compare-with-D.txt`.

## Remaining readiness gates

The immutable GMRES preprojection snapshot, same-data independent GMRES/MLMG solves, and clean-upstream reproduction did not complete in this round; no attribution to MLMG or upstream is justified. The GMRES issue therefore remains unresolved and requires follow-up before merge. The feature MPI build did not finish sufficiently to report the requested final rank-parity/regression matrix. GPU tests were unavailable on this node (no NVIDIA device/toolchain). PR #4253 was isolated and tested with and without its FillPatch density-ghost patch; the failing output was bitwise unchanged, so that PR does not contribute to G8 and no changes from it were imported. Verdict: G8 restart parity is cleared for equivalent lifecycle; overall merge readiness remains blocked on GMRES attribution and the incomplete final regression matrix.
