# Current experiment status

Updated 2026-10-07.

## Remote branch

The work is on wang1202/ERF:terrain_opt, based on remote development at
e50b1611bc987a5188160dfa5aa5061232080859. The branch includes the requested
agent instructions, the copied MyBuildTerrainOpt/cmake.sh, the AMReX
MLTerrainPoisson source pin, and reproducible benchmark decks.

## Build and baseline

The clean Release CUDA build completed successfully at 15:23 MDT using the
copied script. Build provenance is in build_provenance.md. AMReX, including
MLTerrainPoisson, and all ERF targets compiled and linked. The build emitted
non-fatal third-party CUDA/Fortran warnings.

The control executable is copied to
benchmarks/baseline_gmres_fft/erf_exec and SHA256-pinned as
0c4b71b2441b1d0e4ed78216e8ce3e8e5aadab546a35d6380336ea82402d482b.
Slurm job 18949177 is submitted to gpu-h100s and is currently pending priority.
It runs from the common read-only 7200 s checkpoint to 9000 s with both Poisson
tolerances at 1e-8.

## Experiment setup and findings

The candidate and control use the same physics, 360 x 360 x 200 grid,
checkpoint, stop time, and tolerances. Their only input difference is the
terrain_poisson_solver choice. Each run has its own executable copy and output
folder.

The existing terrain solver does not support disabling FFT by setting
erf.use_fft = 0. The new AMReX MLTerrainPoisson multigrid solver is the
non-FFT alternative being tested. No speedup or numerical-equivalence result
has been measured on this branch yet.

MLTerrainPoisson supports homogeneous boundary conditions and sets Neumann
face flux to zero. ERF terrain boundary treatment can retain slope-related
cross-term flux at those faces, so boundary flux, post-projection divergence,
and LES-field differences remain required acceptance checks.

The older strict campaign measured 0.811332 s/step and median 19 GMRES
iterations on ERF commit b0123b8. That is historical context only; the
requested development branch is being remeasured from the same checkpoint.

## Next steps

1. Monitor and complete the GMRES+FFT control run.
2. Add terrain_poisson_solver = gmres_fft|mlmg, with GMRES+FFT remaining the
   default, and guard unsupported layouts and boundary modes.
3. Rebuild and run the candidate from the same checkpoint.
4. Compare solver and full-step timing, residual/divergence, primary flow
   fields, turbulence statistics, and surface diagnostics before accepting
   any optimization.
