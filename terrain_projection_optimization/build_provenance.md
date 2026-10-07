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
