# Current experiment status

Updated 2026-10-07 21:10 MDT.

## Remote branch

The work is on wang1202/ERF:terrain_opt, based on remote development at
`e50b1611bc987a5188160dfa5aa5061232080859`. The branch includes the requested
agent instructions, copied `MyBuildTerrainOpt/cmake.sh`, the AMReX
MLTerrainPoisson pin, and reproducible benchmark decks. The experimental
solver remains opt-in; GMRES+FFT is the default.

## Builds

The clean Release CUDA build succeeded. A follow-up rebuild also succeeded
after replacing the unsupported MLMG flux query with
`MLTerrainPoisson::compFlux`, which provides the terrain-aware face fluxes
needed by ERF's momentum correction.

- Control executable SHA256: `0c4b71b2441b1d0e4ed78216e8ce3e8e5aadab546a35d6380336ea82402d482b`
- Corrected candidate executable SHA256: `9a6eaf213c4f596c6301ed103e77fcc3ca1768e60e2e704f5dd1a7a7b11b2101`

## Benchmark results

Control job 18949177 completed successfully on one H100 (30m41s Slurm
elapsed), covering model time 7200 to 9000 s in 2195 steps. Mature-half
measurements from `benchmarks/summarize_logs.py`:

- 0.81736 s/step and 3632 simulated seconds per wall hour.
- Median 19 GMRES iterations; solver work averaged 0.67414 s/step, or 82.48%
of measured step time.
- Median relative GMRES residual 8.20e-9. Maximum post-projection divergence
was 1.79e-9 in L-infinity and 1.26e-7 in the unnormalized L2 norm.
- No warnings or errors were found in the completed control log.

Candidate job 18949428 failed after its first projection. MLMG converged in
two cycles (reported absolute residual 3.17e-9), then the generic
`MLMG::getFluxes` call reached AMReX's unimplemented `MLLinOp::getFluxes`
abort. That failed attempt's log, preflight, backtrace, and partial outputs
are preserved under `benchmarks/candidate_mlmg/attempts/18949428/`. It provides
no candidate speed or post-projection divergence result.

The candidate now calls `MLTerrainPoisson::compFlux`; its corrected binary is
ready for a clean rerun. Performance and numerical equivalence remain unknown
until that run and field/statistics checks complete. The solver's homogeneous
Neumann flux behavior still requires explicit boundary-flux validation.

## Next steps

1. Resubmit the corrected candidate against the completed control interval.
2. Compare mature step and solver time, residuals, cycles, and post-projection
divergence using `benchmarks/summarize_logs.py`.
3. Compare native profiles, surface diagnostics, plotfile fields, and boundary
flux behavior before accepting any speedup.
4. Extend the run for statistical validation only if the candidate is both
faster and numerically acceptable.
