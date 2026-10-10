# Draft upstream issue: terrain GMRES reports convergence with a large true residual

**Status: draft for user review only. Not published.**

## Summary

In a one-rank CPU DOUBLE/FFT TerrainHill run with a full-height two-level refinement region, the terrain GMRES solver can return success and an Arnoldi-estimated residual near `1e-9` while an immediate independent evaluation of `TerrainPoisson::apply(phi)-rhs` remains between `1e-3` and `1e-2` in L-infinity/L2. The following projection event can also return success while its independent residual grows. The default solver's post-projection divergence therefore fails the configured `1e-6` criterion.

This draft does not claim a proven mathematical root cause or an exact same-state upstream/feature reproduction. Existing historical same-pin A/B logs and the current feature-plus-upstream integration run show persistence across tested configurations, but do not establish that the defect is entirely pre-existing or independent of MLMG. This draft is therefore evidence for a possible future issue, not ready to file as a confirmed upstream defect.

## Reproduction evidence currently available

The discrepancy was reproduced on refreshed feature source `4237b43010efd66607d31b196055ac0df21ffc39` (runtime build string `26.09-578-g4237b43010ef-dirty`), AMReX `5e1ef8467df5f067be4d673c6efe45e1ffb3c9b1`, with input SHA-256 `554cfbd008bfddcd4ce3fa3a36f00eaabadcfebf61f0d1fbd4d6c8464ab19748` (`TerrainHill.i`) and sounding SHA-256 `71a5ebf7206ecfb707c8f5a1fdf151734dcdf8275fa5abc8c19be0da8be4a0a6` (`input_sounding`). The serial CPU DOUBLE/FFT command was `erf_exec TerrainHill.i erf.anelastic=1 erf.terrain_poisson_solver=gmres_fft erf.mg_v=2 max_step=1 amr.max_grid_size=8 erf.box1.in_box_hi=1200 600 600 erf.v=2`. At level 1, projection call 4, time `0.25`, region `[16:47,8:23,0:47]`, status is 0 after 39 iterations:

| Quantity | Value |
|---|---:|
| Initial norm | `7.233900501e-2` |
| GMRES reported norm | `7.283288539e-9` |
| Independent true-residual L-infinity | `5.354735552e-3` at `(30,14,1)` |
| Independent true-residual L2 | `9.085317844e-2` |
| Post-projection divergence L-infinity / L2 | `5.35473556e-3` / `9.085317844e-2` |

The raw correction flux on artificial low-y/high-y faces is `2.701562828e-3` / `3.312516965e-3` over 1,536 samples per side; x is about `1e-10`, z is zero to roundoff. These values make the artificial y boundaries a useful place to inspect, but do not prove causality.

A newer user report names source `fe07312367af067ee2f04ca36b3b736e79312d4e`, Slurm job `18983861`, GMRES estimate `7.283288539e-9`, and independent residual L-infinity/L2 `5.354735552e-3` / `9.085317844e-2`. The raw log for that report should accompany a public issue; it is not in the current checkout.

After integrating current upstream development `b364748d49142d5d8b0b2f8fcdb5f0b2a3dac407` with the feature in isolated merge worktree `0ffee843a60b356ae4917e184189888cb7f6da83`, the same case reproduced at level 1, call 4, time `0.25`: estimated residual `7.339462237e-9`, independent Linf `1.388556916e-3` at `(30,15,1)`, L2 `1.818949566e-2`. Raw log: `/tmp/terrain-mlmg-gmres-upstream4231/run.log`. Because the integration tree includes the feature and the state was not snapshotted/replayed across solvers, this is not clean-upstream proof.

## Exact controls to include before filing

Capture the complete preprojection solve state and compatible RHS at the first failing call, then run GMRES and MLMG independently on copies of those same bytes. Record source and AMReX commits, input/sounding hashes, rank/launch mode, complete options, all valid arrays and required ghosts, coefficients, metrics, geometry, regional offsets and effective BCs. Test repeated operator application with different ghost sentinels, linearity, preconditioner input mutation, boundary/interface residual strips, and flux-divergence versus operator output. Save per-backend true residuals, potentials and correction fluxes. A restart checkpoint is not a substitute for this in-call capture.

## Scope request

Once a clean upstream-only same-state reproduction establishes the defect and a root-cause experiment identifies the violated invariant, please advise whether the residual estimate/operator mismatch is expected for terrain GMRES with these artificial regional boundaries. If confirmed as a pre-existing GMRES/FFT defect, it should be fixed and regression-tested in a separate change; the opt-in MLMG PR should keep `gmres_fft` as its default and should not absorb an unrelated rewrite of that backend. Do not file this draft until those controls exist.

## Artifacts

Current-source raw evidence: `/tmp/terrain-mlmg-gmres-current/run.log`; it includes the copied inputs and generated initial plotfile. Historical artifact: `MyBuild/PRReady/fullheight-current-source/gmres_fft/run.log` (source `44b945bc5c1b`). No issue has been filed.
