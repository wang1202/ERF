# Baseline provenance

## Existing measured reference

The originating performance/statistical campaign measured the strict terrain
case over simulated hours 3–6 with:

- Poisson relative and absolute tolerances: 1e-8
- Median GMRES iterations: 19; mean 19.002; maximum 20
- Mature wall time: 0.811332 s/step
- Throughput: 3699.62 simulated seconds per wall hour
- Grid: 360 x 360 x 200
- Integrator: RK2; CFL: 0.6; terrain smoothing: 1
- Executable SHA256: 46d2d358cfb2105da3bb1dc45037848f198db168a1ca10f18c50cb914185bb93

Original provenance records ERF commit b0123b8c748c57896ce09a1f8057edd1b8ef9228.
The requested remote development branch starts at e50b1611bc987a5188160dfa5aa5061232080859.
The optimization benchmark will therefore remeasure the GMRES+FFT baseline and
candidate from one executable built at the requested remote branch, rather than
treating the older timing as a matched result.

## Common restart

The existing common checkpoint is
/kfs2/projects/erf/aaronwang/ERF/MyBuild/Exec/TerrainToleranceValidation/common_checkpoint.
Its job_info records artifact_time_seconds=7200. It is used read-only. All new
benchmark outputs will live under this branch's terrain_projection_optimization
tree. The existing checkpoint and original experiment directories will not be
modified.

## Current optimization build

- ERF commit: e50b1611bc987a5188160dfa5aa5061232080859.
- Original pinned AMReX commit: 1de18774af5d8a368fa2a398f048e205a2c9f84f.
- Candidate AMReX commit: 5e1ef8467df5f067be4d673c6efe45e1ffb3c9b1.
- Build script: MyBuildTerrainOpt/cmake.sh, copied from MyBuild/cmake.sh.
- Executable SHA, compiler/module versions, and runtime allocation will be
  appended after the build and preflight.
