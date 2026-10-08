# Current experiment status

Updated 2026-10-07 21:14 MDT.

## Remote branch

The work is on wang1202/ERF:terrain_opt, based on remote development at
`e50b1611bc987a5188160dfa5aa5061232080859`. The branch includes the requested
agent instructions, copied `MyBuildTerrainOpt/cmake.sh`, the AMReX
MLTerrainPoisson pin, and reproducible benchmark decks. The experimental
solver remains opt-in; GMRES+FFT is the default.

## Builds

The clean Release CUDA build succeeded. The corrected candidate also compiled
and linked incrementally after replacing the unsupported MLMG flux query
with `MLTerrainPoisson::compFlux`, which calculates the terrain-aware face
fluxes needed for ERF's momentum correction.

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

Initial candidate job 18949428 failed after its first projection. MLMG
completed two cycles with an absolute residual of 3.17e-9, then the generic
`MLMG::getFluxes` call reached AMReX's unimplemented
`MLLinOp::getFluxes` abort. No post-projection divergence or candidate timing
is available from that attempt. Its local log, preflight, backtrace, and
partial outputs are preserved in `benchmarks/candidate_mlmg/attempts/18949428/`.

The flux call is fixed in source commit `33e6b4c`, and corrected candidate job
18953921 has been resubmitted with the new executable. It is currently pending
on scheduler priority. Thus non-FFT performance and numerical equivalence are
still unknown. The solver's homogeneous Neumann flux behavior still requires
explicit boundary-flux validation.

## Next steps

1. Monitor corrected candidate job 18953921 through H100 preflight and
completion.
2. Compare mature step and solver time, residuals, cycles, and post-projection
divergence using `benchmarks/summarize_logs.py`.
3. Compare native profiles, surface diagnostics, plotfile fields, and boundary
flux behavior before accepting any speedup.
4. Extend the run for statistical validation only if the candidate is both
faster and numerically acceptable.
