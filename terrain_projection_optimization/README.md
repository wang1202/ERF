# Terrain projection optimization

This work is isolated on the terrain_opt branch of wang1202/ERF, based on
the requested development branch. The scientific reference remains the
terrain LES with erf.poisson_reltol = 1.e-8 and erf.poisson_abstol = 1.e-8.

## Plan

1. Record the ERF, AMReX, compiler, build, executable, and input provenance.
2. Build and benchmark the existing GMRES+FFT path from a common 7200 s
   checkpoint. Add lightweight profiling around projection setup, solver,
   preconditioner, operator application, and flux recovery before optimizing.
3. Integrate AMReX MLTerrainPoisson behind an experimental runtime selector;
   keep the existing GMRES+FFT implementation as the default.
4. Benchmark both solvers with identical inputs, decomposition, checkpoint,
   and simulated interval. Compare solve residual, actual post-projection
   divergence, flow fields, turbulence statistics, and surface diagnostics.
5. Keep only isolated changes that improve measured whole-step and projection
   wall time without unacceptable numerical differences. Record rejected
   alternatives and remaining risks.

## Build

MyBuildTerrainOpt/cmake.sh is copied from the existing MyBuild/cmake.sh.
Run it from MyBuildTerrainOpt so all CMake caches, logs, and binaries stay in
this experiment build directory:

    cd MyBuildTerrainOpt
    bash cmake.sh

## Current state

- ERF base: e50b1611bc987a5188160dfa5aa5061232080859 (requested remote
  development).
- Original ERF solver remains the default.
- AMReX MLTerrainPoisson candidate: 5e1ef8467df5f067be4d673c6efe45e1ffb3c9b1.
- No optimization is accepted until the matched benchmark and physical checks
  are complete.
