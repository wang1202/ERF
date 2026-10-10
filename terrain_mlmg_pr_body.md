# Draft upstream PR body

**Proposed title:** Add an opt-in MLMG terrain projection solver

## Summary

Adds an opt-in terrain projection backend selected by `erf.terrain_poisson_solver=mlmg`, using AMReX `MLTerrainPoisson`. The existing `gmres_fft` backend remains the default. The feature fills rho0-weighted coarse/fine horizontal momentum ghosts around projection and checks solver compatibility, true residual, divergence, and effective-Neumann correction flux.

## Supported scope and dependency

The intended scope is 3D `StaticFittedMesh`/`VariableDz`, ERF levels 0 and 1, and one connected rectangular solve region on each active level. This is not a composite-AMR Poisson solve. Disconnected regions, more than two active levels, unsupported mesh geometry, and `erf.use_real_bcs=true` are not supported and are rejected.

The implementation requires AMReX `5e1ef8467df5f067be4d673c6efe45e1ffb3c9b1`, which provides `AMReX_MLTerrainPoisson.H`; current upstream `development` pins `1de18774af5d8a368fa2a398f048e205a2c9f84f`. The AMReX dependency update needs explicit review.

## Validation status

Historical CPU DOUBLE/SINGLE, H100, partial-height and performance results are in the local merge-readiness report and are not all rerun on current feature tip. The refreshed current-source CPU DOUBLE/FFT executable builds and links. Its restart/regrid test verifies the changed BoxArray and all projection checks, but the continuous/restart cell-field comparison fails at `plt00002` with large velocity differences; face comparisons at later shared times therefore remain unverified.

The default GMRES full-height two-level case also has an unresolved residual discrepancy: the feature-only run reports `7.2833e-9` while independent application leaves Linf `5.3547e-3` / L2 `9.0853e-2`; the feature-plus-current-upstream integration run reports `7.3395e-9` while independent Linf is `1.3886e-3` / L2 `1.8189e-2`. These runs reproduce the symptom but do not provide the required identical preprojection snapshot, operator repeatability/ghost tests, or GMRES/MLMG shadow solve. No cause or maintainer disposition is claimed. The default is unchanged.

The focused current-upstream terrain-diffusion regressions `TerrainStress_Ridge_ImplicitMetric`, `TerrainStress_Hill3D`, and `TerrainStress_Ridge_ImplicitMetricAgree` pass. The integrated changed-layout restart/regrid regression still fails its continuous comparison at `plt00002` (level-1 z-velocity absolute error `0.2331`, relative error `0.07577`).

Matched historical one-level measurements showed MLMG about 18% faster in summed CPU timestep time and about 21% slower on H100. These are hardware-specific results, not a general performance claim.

## Review status

The branch is **not ready to open for upstream merge review**. G8 currently fails and several required current-source numerical, rank-parity, boundary, mapping, and full GMRES diagnostic gates remain open. This text is a local draft; no PR has been opened.
# G8 follow-up evidence (2026-10-10)

The initial ~0.233 m/s comparison was not lifecycle-equivalent. A no-regrid continuous/restart pair matched bit-for-bit. For the changed-BoxArray test, a continuous run that reblocked 8/8/8 to 16/16/8 at the same physical time as the restart remake matched all `plt`, `pltU`, `pltV`, and `pltW` fields at t=0.5, 1.0, and 1.5 bit-for-bit. The earlier continuous reference began with the final 16/16/8 layout and skipped `RemakeLevel`. No MLMG/terrain production correction is indicated by this test. The GMRES residual discrepancy remains unresolved: immutable snapshot, independent solver replay, and clean-upstream confirmation were not completed. Overall merge readiness remains pending that diagnosis and the final regression matrix.
