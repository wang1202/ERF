# GMRES true-residual discrepancy: evidence and status

Date: 2026-10-10  
Candidate ref: `terrain-mlmg-pr-ready` at `4237b43010efd66607d31b196055ac0df21ffc39`  
AMReX pin: `5e1ef8467df5f067be4d673c6efe45e1ffb3c9b1`

## Finding

The Arnoldi residual estimate is not a reliable proxy for the independently evaluated residual in the failing full-height, two-level TerrainHill run. The exact mathematical cause is **not isolated**. There is not yet an immutable preprojection snapshot or a shadow solve of GMRES and MLMG using the same RHS and coefficients, so this evidence does not establish that the discrepancy is wholly pre-existing, nor that the y-boundary fluxes cause it.

The failing case was rerun on the refreshed feature source using the current CPU DOUBLE/FFT executable, serial MPI-off launch, `max_step=1`, `amr.max_grid_size=8`, and `erf.box1.in_box_hi=1200 600 600`. The executable banner/provenance reports ERF `26.09-578-g4237b43010ef-dirty` (tracked source tip `4237b43010efd66607d31b196055ac0df21ffc39`) and AMReX `26.10-34-g5e1ef8467df5` (`5e1ef8467df5f067be4d673c6efe45e1ffb3c9b1`). Input and sounding hashes are `554cfbd008bfddcd4ce3fa3a36f00eaabadcfebf61f0d1fbd4d6c8464ab19748` and `71a5ebf7206ecfb707c8f5a1fdf151734dcdf8275fa5abc8c19be0da8be4a0a6`. Raw log: `/tmp/terrain-mlmg-gmres-current/run.log`; run directory also contains the exact copied inputs and first plotfile metadata.

At the first failing event in the refreshed-source raw log, level 1, `projection_call=4`, `time=0.25`, region `[16:47,8:23,0:47]`, GMRES returns status 0 after 39 iterations. It reports `7.283288539e-9` from an initial norm `0.07233900501`. Immediately after the solve, before `getFluxes` and `ImposeBCsOnPhi`, the existing diagnostic computes `TerrainPoisson::apply(phi)-rhs`: L-infinity `5.354735552e-3` at cell `(30,14,1)` and L2 `9.085317844e-2`. Post-projection divergence is L-infinity `5.35473556e-3`, L2 `9.085317844e-2`. The following GMRES event also reports success (`3.441324621e-9`) while its independent L2 residual is `6.656139088e-3`.

The same event has raw artificial-boundary correction-flux maxima of `2.701562828e-3` on low-y and `3.312516965e-3` on high-y (1,536 samples per side); x maxima are `9.00e-11` and `1.08e-10`, and z maxima are zero to roundoff. Map factors are unity in this case. This spatial clue is compatible with a boundary/operator mismatch, but there is no repeatability or ghost-sentinel experiment proving that explanation. The source's `TerrainPoisson::apply` fills ghosts through `apply_bcs`; a single residual evaluation after GMRES does not determine whether the Krylov matrix-vector operation was invariant over the solve.

The older artifact `MyBuild/PRReady/fullheight-current-source/gmres_fft/run.log` identifies source `44b945bc5c1b` and differs numerically; it is retained as historical evidence only. The current-source rerun above reproduces the values in the newer failure report. The current raw run's first failing event is therefore independently available, but still does not contain the requested immutable snapshot/replay.

## What the A/B evidence establishes

The October 9 readiness report records clean-upstream and feature default-GMRES controls at the same AMReX pin with matching output fields and matching divergence failure. That is evidence of no measured feature regression for those distinct runs. It is **not** a same-RHS solver replay and does not localize the cause. A separate integration worktree at merge commit `0ffee843a60b356ae4917e184189888cb7f6da83` combines candidate `4237b43010efd66607d31b196055ac0df21ffc39` with current upstream development `b364748d49142d5d8b0b2f8fcdb5f0b2a3dac407`, using AMReX `5e1ef8467df5f067be4d673c6efe45e1ffb3c9b1`. The same case and input/options reproduce a residual gap at level 1, call 4, time 0.25: estimated `7.339462237e-9`, independent Linf `1.388556916e-3` at `(30,15,1)`, L2 `1.818949566e-2`. Raw log: `/tmp/terrain-mlmg-gmres-upstream4231/run.log`. This confirms persistence after integrating #4231; since this tree contains the candidate feature, it is not a clean upstream-only replay and cannot prove the issue is wholly pre-existing.

No current test has yet recorded repeated `A(phi)` on identical valid input with changed ghost sentinels, linearity, post-`getFluxes` residual, valid/ghost and interface strip norms, flux-form versus production-operator agreement, preconditioner RHS/ghost mutation, Arnoldi residual at each iteration, or both backends acting on copies of one captured state. The root cause remains unresolved. No production solver correction is justified by the available measurements. In particular, the refreshed work did not complete the requested same-call state capture and independent GMRES/MLMG replay, so the conditional “entirely pre-existing” verdict is not established and an upstream issue should not be published on this evidence.

## Reproduction and next diagnostic

The executed one-rank CPU DOUBLE/FFT command from the repository root was:

```sh
erf_exec TerrainHill.i erf.anelastic=1 erf.terrain_poisson_solver=gmres_fft erf.mg_v=2 max_step=1 amr.max_grid_size=8 erf.box1.in_box_hi=1200 600 600 erf.v=2
```

It ran in `/tmp/terrain-mlmg-gmres-current`, with the fixture input and sounding hashes above. `erf_exec` was serial (`MPI OFF`), so no launcher was used. `MyBuild/PRReady/fullheight.sbatch` is an older max-grid-size-64 A/B asset; it is not the executed command.

The required next experiment is an opt-in diagnostic build that snapshots valid and ghost state, RHS before/after compatibility handling, coefficients/metrics/geometry/BCs and layout at the first bad solve entry; creates two isolated copies; and runs repeatability/ghost tests plus both backends on the same immutable bytes. Compare operator residuals by interior/boundary/interface strips and compare `A(phi)` with the signed flux-divergence form. No restart from a timestep checkpoint can stand in for this in-call state capture.

## Disposition

Leave the default `gmres_fft` unchanged. Keep this issue separate from the MLMG implementation unless a same-state experiment demonstrates feature interference or upstream maintainers explicitly request a combined scope. The issue draft in `terrain_gmres_upstream_issue_draft.md` is for user review only and has not been published.
