# Current experiment status

Updated 2026-10-07.

## Remote branch

The work is on wang1202/ERF:terrain_opt, based on remote development at
e50b1611bc987a5188160dfa5aa5061232080859. The branch is pushed through
ef0a779. It includes the requested agent instructions, the copied
MyBuildTerrainOpt/cmake.sh, the AMReX MLTerrainPoisson source pin, and
reproducible benchmark inputs and job script.

## Build

A clean Release CUDA build was started from MyBuildTerrainOpt using the
copied script. The last observed progress was 70%, compiling ERF sources after
AMReX including MLTerrainPoisson built successfully. No ERF compile failure has
appeared. The initial module-load warnings for Cray MPI/LibSci were followed by
successful CMake configuration using the explicit compiler and library paths
in the copied script. The build is still running; its log is ignored inside
the isolated build directory.

## Experiment setup

The two decks under benchmarks/ use the same 7200 s strict-tolerance
checkpoint, physics, 360 x 360 x 200 grid, and 9000 s stop time. Both retain
relative and absolute Poisson tolerances of 1e-8. The candidate deck adds only
erf.terrain_poisson_solver = mlmg; the baseline will use the existing
GMRES+FFT path. The checkpoint remains read-only. The one-H100 Slurm job has
not yet been submitted.

## Findings and limits

The current terrain path does not support disabling FFT by setting
erf.use_fft = 0; the new AMReX MLTerrainPoisson multigrid solver is the
non-FFT alternative being tested. No speedup or numerical-equivalence result
has been measured on this branch yet.

MLTerrainPoisson uses homogeneous boundary conditions and zero flux on
Neumann faces. ERF's terrain boundary treatment can retain slope-related
cross-term flux at those faces, so boundary flux, post-projection divergence,
and LES-field differences are acceptance checks, not assumptions.

The older strict campaign measured 0.811332 s/step and median 19 GMRES
iterations on ERF commit b0123b8. That is historical context only; the
requested remote branch differs and will be remeasured from the same checkpoint.

## Next steps

1. Finish the clean build and record compiler, executable hash, and build
   provenance.
2. Run the GMRES+FFT control from the common checkpoint.
3. Integrate MLTerrainPoisson behind the selector with GMRES+FFT still the
   default, rebuild, and run the candidate.
4. Compare solver time, full-step time, residual/divergence, primary flow
   fields, turbulence statistics, and surface diagnostics before accepting
   any optimization.
