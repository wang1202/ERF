# Clean baseline build provenance

Build completed successfully at 2026-10-07 15:23:31 MDT with the copied MyBuildTerrainOpt/cmake.sh.

- ERF source commit: 54baec97f031824dc84fc27beb78f01ef55246db
- ERF branch: terrain_opt
- AMReX submodule commit: 5e1ef8467df5f067be4d673c6efe45e1ffb3c9b1
- Build type: Release
- ERF dimension: 3
- ERF FFT: ON
- CMake: 3.29.2 requested by cmake.sh
- C++ compiler: /opt/cray/pe/mpich/8.1.28/ofi/gnu/10.3/bin/mpicxx (g++ (GCC) 12.2.1 20221121 (Red Hat 12.2.1-7))
- Fortran compiler: /opt/cray/pe/mpich/8.1.28/ofi/gnu/10.3/bin/mpifort (GNU Fortran (GCC) 12.2.1 20221121 (Red Hat 12.2.1-7))
- CUDA compiler: /nopt/cuda/12.9/bin/nvcc (nvcc: NVIDIA (R) Cuda compiler driver)
- Build log: MyBuildTerrainOpt/cmake_build_20261007_135016.log
- Executable: MyBuildTerrainOpt/Exec/erf_exec
- Executable SHA256: 0c4b71b2441b1d0e4ed78216e8ce3e8e5aadab546a35d6380336ea82402d482b
- Control copy: terrain_projection_optimization/benchmarks/baseline_gmres_fft/erf_exec
- Control-copy SHA256: 0c4b71b2441b1d0e4ed78216e8ce3e8e5aadab546a35d6380336ea82402d482b

The initial module load reported that cray-mpich/8.1.28 and cray-libsci/23.12.5 were unavailable by module name. CMake still configured with the explicit Cray wrapper and library paths in the copied script, and the complete build succeeded. Third-party CUDA/Fortran warnings were non-fatal.

The runtime GPU, driver, loaded libraries, and input/checkpoint hashes are captured by the Slurm preflight in the case directory.

## Experimental candidate build

The MLTerrainPoisson integration was added in ERF commit 5960a55 and built
incrementally against the same CMake configuration and AMReX commit. The
single changed ERF translation unit compiled, and the erf_exec target linked
successfully after restoring the CUDA 12.9 and Cray libfabric paths used by
the clean build environment.

- Candidate executable: terrain_projection_optimization/benchmarks/candidate_mlmg/erf_exec
- Candidate SHA256: a43448f05279913c37e1252b14632fc694d96e7a0f3bb27337cd30a32d46d5d2
- Candidate source commit: 5960a55
- Link check: nm confirms MLTerrainPoisson constructor, setZPhys, setAreas, and setDetJ symbols are present.

The login node has no GPU driver library, so its ldd output reports
libcuda.so.1 as unavailable. The H100 job performs the authoritative ldd and
nvidia-smi preflight before launching each case.


## Corrected candidate build and runtime attempt

The first candidate executable (SHA256
`a43448f05279913c37e1252b14632fc694d96e7a0f3bb27337cd30a32d46d5d2`) was
built from source commit `5960a55`. Job 18949428 ran on an NVIDIA H100 80GB
with driver 550.54.15, converged the first projection in two MLMG cycles, and
then aborted because this AMReX MLTerrainPoisson implementation does not
provide the generic `MLMG::getFluxes` method.

Source commit `33e6b4c` replaces that call with the operator's supported
`MLTerrainPoisson::compFlux` face-flux calculation. The incremental CUDA
rebuild and executable link succeeded. Corrected candidate executable SHA256:
`9a6eaf213c4f596c6301ed103e77fcc3ca1768e60e2e704f5dd1a7a7b11b2101`.
Retry job 18953921 was submitted with this executable and is pending; its
Slurm preflight will capture the runtime libraries, GPU, and executable hash.
