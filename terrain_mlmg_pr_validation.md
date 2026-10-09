# Terrain MLMG PR validation

Date: 2026-10-08

## Verdict

**Draft PR only.** The lower-z BC selection, fail-fast checks, documentation, focused normal-wall MLMG run, and matched short GMRES/MLMG runs are in place. Broader validation remains incomplete, and the focused GoogleTest target is still building.

## Changes made

- Source/LinearSolvers/ERF_SolverUtils.H keeps gmres_fft as the default, validates solver names and MLMG restrictions, and adds a terrain-specific low-side BC helper. For nonperiodic z it forces zlo to Neumann, preserves periodic z, and leaves generic projection BC handling unchanged.
- Source/LinearSolvers/ERF_PoissonSolve.cpp validates MLMG at the start of project_momenta, before coarse/fine momentum fills and the small-RHS skip; reports an actionable error for terrain GMRES without FFT; prints the chosen terrain solver once when mg verbosity is enabled; and uses the terrain-specific zlo helper only for MLMG.
- Tests/Unit/LinearSolvers/ERF_GTestTerrainPoissonBC.cpp checks lower-z outflow mapping, other faces, periodic z, default selection, and MLMG diagnostics. Tests/Unit/CMakeLists.txt registers it.
- Docs/sphinx_doc/Inputs.rst documents the default, opt-in keyword, FFT requirements, restrictions, lower-z convention, lateral Neumann difference, and selection logging.
- AMReX internals and the pinned submodule were not changed.

## Validation

### Build and tests

- Production 3D CUDA H100 build and relink passed with CMake; architecture 90, MPI and FFT enabled. Refreshed runtime ERF hash: 040db7a7d280-dirty; executable SHA-256: 71479bac58d65d5ffbf1c46832f27a74ebdc0821483fba64387800b1c014ec7e.
- Focused GoogleTest target built successfully. On H100, Tests/Unit/erf_unit_tests --gtest_filter=TerrainProjectionBC.*:TerrainMLMGConfiguration.* passed all 5 tests (2 BC and 3 configuration); Slurm job 18967744 exited 0.
- Restoring the pre-test CMake options with a configure-only attempt failed because this shell could not satisfy the MPI_CXX_WORKS probe. CMakeCache has both test options OFF; the incomplete reconfiguration did not rebuild or remove the completed artifacts.
- git diff --check passed.
- AMReX submodule: 5e1ef8467df5f067be4d673c6efe45e1ffb3c9b1.

### Short matched terrain runs

Both runs restarted from chk09433 at t=7650 s, used one MPI rank and one H100, a 360 x 360 x 200 domain, periodic x/y, zlo surface_layer, zhi SlipWall, one full-domain box, stop_time 7650.5 s, and the same executable. Two steps advanced to step 9435. The input omitted erf.poisson_reltol and erf.poisson_abstol; SolverChoice raises these tolerances to at least 1e-6.

| Run | Step/stage | Divergence before (L_inf, L2) | Divergence after (L_inf, L2) | MLMG iterations / final residual |
|---|---:|---:|---:|---|
| MLMG | 9434 / 1 | 2.077e-3, 1.425e-1 | 5.040e-10, 1.737e-7 | 3 / abs 4.179e-10, rel 3.110e-7 |
| MLMG | 9434 / 2 | 2.258e-3, 1.657e-1 | 5.617e-10, 1.803e-7 | 3 / abs 4.757e-10, rel 2.918e-7 |
| MLMG | 9435 / 1 | 2.751e-4, 2.034e-2 | 4.289e-9, 6.607e-7 | 2 / abs 3.207e-9, rel 1.900e-5 |
| MLMG | 9435 / 2 | 4.427e-5, 2.866e-3 | 5.474e-10, 1.636e-7 | 2 / abs 4.260e-10, rel 1.153e-5 |
| GMRES+FFT | 9434 / 1 | 2.077e-3, 1.425e-1 | 8.332e-10, 1.403e-7 | not reported |
| GMRES+FFT | 9434 / 2 | 2.258e-3, 1.657e-1 | 7.195e-10, 1.403e-7 | not reported |
| GMRES+FFT | 9435 / 1 | 2.751e-4, 2.034e-2 | 1.079e-9, 1.404e-7 | not reported |
| GMRES+FFT | 9435 / 2 | 4.426e-5, 2.866e-3 | 8.406e-10, 1.404e-7 | not reported |

All four distinct short-run configurations exited 0 at t=7650.5 s; the relinked post-test MLMG confirmation also exited 0. Slurm jobs/logs: 18966599 (normal-wall MLMG), 18966892 (explicit GMRES+FFT), 18966928 (keyword omitted), 18967101 (zlo outflow MLMG). The no-keyword case selected gmres_fft; explicit gmres_fft selected the same solver and matched default divergence values to rounding precision. Maximum post-projection L_inf divergence was 4.3e-9 for MLMG and 1.1e-9 for GMRES+FFT.

A yt comparison of final MLMG and GMRES+FFT plotfiles removed the mean pert_pres offset before comparing pressure. Relative RMS differences were 8.12e-11, 2.74e-11, and 7.15e-10 for x/y/z velocity; 2.32e-6 for pert_pres; and 1.92e-12 for theta. Density matched exactly. This is a two-step comparison, not a long-term LES equivalence test.

The normal-wall integration exercises the solver but does not measure lower-face correction-flux norms. The focused BC-selection test checks outflow-to-Neumann mapping. A separate restart with zlo.type=outflow also completed through step 9435 with exit code 0; its post-projection L_inf divergence ranged from 5.0e-10 to 4.3e-9, although the lower-face flux norm was not instrumented. Invalid settings are tested at helper level, not through executable negative-input tests. Multi-box/MPI, z-split, lateral Neumann/outflow variations, no-FFT build, explicit finite-field assertions, and the full existing default regression suite remain unverified.

Initial login-node and first Slurm test attempts loaded the CUDA stub and aborted before assertions; the corrected H100 job passed. The post-test relinked MLMG run (job 18967747, log MyBuild/ERF_TerrainMLMGPostFixRelink/slurm-18967747.out) reported ERF hash 040db7a7d280-dirty, reached step 9435/t=7650.5 s, and exited 0 with post-projection L_inf divergence at or below 4.29e-9.

## Performance

The short run is not a benchmark: MLMG reported 19.39 s total application time and 0.49/0.31 s timestep times; GMRES reported 21.12 s total and 0.92/0.58 s. Initialization and output dominate these two-step samples; they do not establish speedup.

Historical matched evidence in MyBuild/ERF_TerrainMLMGConfirm/results.md reports one 450-simulation-second pair on one H100: GMRES+FFT 8:03 executable wall and 375.298 s summed Poisson time versus MLMG 3:58 and 132.358 s. It measured 2.03x executable and 2.84x Poisson-solve speedup, with 1,114 solves and similar final divergence. This predates the review fixes and is one pair, with no setup/flux phase timing or repeat. It is context only, not a current performance claim.

## Known limitations

- MLMG remains level-0 only, homogeneous-BC only, VariableDz with StaticFittedMesh only, and requires one full-domain subdomain.
- AMReX MLMG sets the entire normal correction flux to zero on lateral Neumann faces; legacy terrain GMRES may retain terrain cross terms.
- No boundary-normal correction-flux norm was recorded.
- No fresh repeated full-timestep benchmark was run.
- No CPU build, no no-FFT build, and no multi-rank/multi-box validation were run.

## Dependency

- Branch terrain_optimization; checkout commit 040db7a7d28039e8f07b3c1aa0250b6ed0102973.
- AMReX is pinned at 5e1ef8467df5f067be4d673c6efe45e1ffb3c9b1 and was not edited.
- Local refs: origin/development...HEAD is 0 1 (feature branch one commit ahead and none behind); no fetch was performed, so remote freshness remains unchecked.
- Tested CMake configuration is 3D CUDA architecture 90 with MPI and FFT enabled. CPU and no-FFT configurations remain untested.
- No push, merge, PR, or remote branch update was performed.

## Remaining PR blockers

1. Measure and assert the lower-z zero-flux norm in the outflow integration.
2. Run executable-level negative tests and the existing default regression suite against upgraded AMReX.
3. Expand solver comparison to box/MPI decompositions and additional boundary variations; capture correction-flux norms.
4. Repeat end-to-end benchmarks with setup, solve, flux, full-timestep, and memory metrics before making a current speedup claim.
