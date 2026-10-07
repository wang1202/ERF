# Current experiment status

Updated 2026-10-07.

## Remote branch

The work is on wang1202/ERF:terrain_opt, based on remote development at
e50b1611bc987a5188160dfa5aa5061232080859. The branch includes the requested
agent instructions, copied MyBuildTerrainOpt/cmake.sh, the AMReX
MLTerrainPoisson pin, and reproducible benchmark decks. The experimental
solver is committed as 5960a55.

## Builds

The clean Release CUDA build completed successfully at 15:23 MDT. Build
provenance and the control executable hash are in build_provenance.md. ERF and
AMReX, including MLTerrainPoisson, compiled and linked.

The experimental MLTerrainPoisson path is opt-in through
erf.terrain_poisson_solver = mlmg. The default remains GMRES+FFT. The
incremental candidate compile and erf_exec link succeeded with the same build
configuration. Runtime guards require 3D, level 0, a full-domain single
subdomain, and homogeneous BCs; the current case satisfies these requirements.

## Benchmark status

The control executable is pinned in
benchmarks/baseline_gmres_fft/erf_exec with SHA256
0c4b71b2441b1d0e4ed78216e8ce3e8e5aadab546a35d6380336ea82402d482b.
The candidate executable is pinned in
benchmarks/candidate_mlmg/erf_exec with SHA256
a43448f05279913c37e1252b14632fc694d96e7a0f3bb27337cd30a32d46d5d2.

Both cases use the same read-only 7200 s checkpoint and stop at 9000 s with
relative and absolute Poisson tolerances of 1e-8. Control job 18949177 and
candidate job 18949428 are submitted to gpu-h100s. At the latest scheduler
check, both were pending: control due to unavailable/reserved H100 nodes and
candidate due to priority. Neither has started, so there are no new runtime,
speedup, or numerical-equivalence results yet.

## Findings and limits

The existing terrain path does not support disabling FFT by setting
erf.use_fft = 0; the new AMReX MLTerrainPoisson multigrid solver is the
non-FFT alternative being tested.

MLTerrainPoisson supports homogeneous boundary conditions and sets Neumann
face flux to zero. ERF's terrain boundary treatment can retain slope-related
cross-term flux at those faces, so boundary flux, post-projection divergence,
and LES-field differences remain required acceptance checks.

The older strict campaign measured 0.811332 s/step and median 19 GMRES
iterations on ERF commit b0123b8. That is historical context only; this
branch is being remeasured from the same checkpoint.

## Next steps

1. Monitor both jobs through H100 preflight and simulation completion.
2. Compare solver time and full-step time at identical simulated intervals.
3. Check residual/divergence, primary flow fields, turbulence statistics,
   and surface diagnostics before accepting any speedup.
4. If the multigrid path is numerically valid and faster, extend the run for
   statistical validation; if not, preserve GMRES+FFT and test another
   tolerance-preserving optimization.
