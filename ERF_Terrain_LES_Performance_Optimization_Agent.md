# ERF Terrain LES Performance Optimization Task

## Objective

Optimize the terrain-enabled anelastic LES path in ERF, focusing on the terrain pressure projection, while preserving numerical correctness and the existing `1e-8` Poisson-tolerance scientific baseline.

Reference ERF commit:

```text
b0123b8c748c57896ce09a1f8057edd1b8ef9228
```

The current terrain case shows that the terrain projection is the dominant runtime cost. Relaxing the Poisson tolerance from `1e-8` to `3e-7` or `1e-6` improves throughput, but the current statistical validation shows persistent differences that still require scientific review. Therefore, **do not use tolerance relaxation as the main optimization strategy**.

The goal is to accelerate the `1e-8` case by improving implementation efficiency.

## Current measured baseline

| Poisson tolerance | Median GMRES iterations | s/step | Relative throughput |
|---:|---:|---:|---:|
| `1e-8` | 19 | 0.81133 | 1.000x |
| `3e-7` | 15 | 0.61309 | 1.319x |
| `1e-6` | 14 | 0.56494 | 1.434x |

The `1e-8` run is the numerical/scientific reference.

Approximate domain:

```text
360 x 360 x 200
```

About:

```text
25,920,000 cells
```

The case uses terrain-fitted coordinates and FFT-preconditioned GMRES for the anelastic projection.

# Main task

Investigate and implement code optimizations that accelerate terrain-enabled ERF LES without relaxing the `1e-8` Poisson tolerance.

Work in stages. Do not combine multiple optimization ideas into one unreviewable change.

For every optimization:

1. Profile first.
2. Implement one logically isolated change.
3. Verify correctness.
4. Benchmark against the exact baseline.
5. Keep or revert based on measured performance and correctness.
6. Record all results.

# Phase 0 — Reproduce and protect the baseline

Before modifying code:

1. Confirm the exact ERF commit.
2. Record the AMReX submodule commit.
3. Record compiler, GPU architecture, MPI layout, build options, and executable SHA256.
4. Confirm that the existing `1e-8` terrain case reproduces approximately:
   - median GMRES iterations: `19`
   - runtime: approximately `0.81 s/step`
5. Preserve the current executable and input files.
6. Create a new optimization branch.

Suggested branch:

```bash
git switch -c optimize-terrain-projection
```

Do not modify production result directories.

Create:

```text
terrain_projection_optimization/
├── baseline/
├── profiles/
├── benchmarks/
├── validation/
└── notes/
```

# Phase 1 — Fine-grained profiling of the current GMRES+FFT solver

Determine exactly where the terrain projection time is spent.

Relevant ERF files:

```text
Source/LinearSolvers/ERF_PoissonSolve.cpp
Source/LinearSolvers/ERF_SolveWithGMRES.cpp
Source/LinearSolvers/ERF_TerrainPoisson.cpp
Source/LinearSolvers/ERF_TerrainPoisson.H
Source/LinearSolvers/ERF_TerrainPoisson_3D_K.H
```

Relevant AMReX files:

```text
Src/LinearSolvers/AMReX_GMRES.H
Src/FFT/AMReX_FFT_Poisson.H
```

Add or use profiling regions to separate at least:

```text
ERF::project_momenta
  compute_divergence
  solvability correction
  coefficient/map-factor preparation
  solve_with_gmres
    TerrainPoisson construction/setup
    GMRES total
      TerrainPoisson::precond
        FFT forward transform
        PoissonHybrid::solve_z
        FFT inverse transform
      TerrainPoisson::apply
        apply_bcs
        terrpoisson_adotx
      GMRES Gram-Schmidt
      GMRES dot products / reductions
      GMRES vector operations
    getFluxes
  coefficient restoration
  momentum correction
  W/Omega conversion
```

Use lightweight timers or AMReX `BL_PROFILE`.

Run at least 100 timesteps after warm-up if practical.

Report:

```text
time per projection
time per timestep
percentage of total timestep
GMRES iteration count
time per GMRES iteration
time per preconditioner application
time per terrain operator application
time in Gram-Schmidt/reductions
```

Create:

```text
terrain_projection_optimization/profiles/baseline_profile.md
```

Do not optimize until this profile is recorded.

# Phase 2 — Test the new native AMReX MLTerrainPoisson

A new AMReX terrain Poisson operator was merged on October 6, 2026:

```text
AMReX commit:
5e1ef8467df5f067be4d673c6efe45e1ffb3c9b1

Title:
MLTerrainPoisson: cell-centered terrain-following Poisson operator (#6030)
```

The AMReX commit reports:

```text
200 x 200 x 176 WPS-like case

RTX 5070:
ERF GMRES+FFT     3.31 s / 53 iterations
MLMG              0.112 s / 3 V-cycles
GMRESMLMG         0.184 s / 3 iterations
```

This is highly relevant to the current `360 x 360 x 200` case.

## Important

Do this work in a separate experimental branch:

```bash
git switch -c test-amrex-mlterrainpoisson
```

Do not immediately replace the production solver.

### Tasks

1. Update AMReX to include at least:

```text
5e1ef8467df5f067be4d673c6efe45e1ffb3c9b1
```

2. Inspect the native AMReX `MLTerrainPoisson` API.
3. Determine the cleanest ERF integration path.
4. Reuse ERF's existing:
   - `z_phys_nd`
   - `ax`
   - `ay`
   - `az`
   - `detJ`
   - boundary-condition definitions
   - terrain projection RHS
   - pressure-gradient correction
5. Add a runtime-selectable experimental solver if practical:

```ini
erf.terrain_poisson_solver = gmres_fft
erf.terrain_poisson_solver = mlmg
```

Keep `gmres_fft` as the default.
6. Do not change the existing default answer.
7. Use:

```ini
erf.poisson_reltol = 1.e-8
```

for the primary comparison.

### Validation

Compare at minimum:

```text
solver residual
post-projection divergence
u
v
w
theta
pressure
resolved TKE
SGS TKE
uw
vw
w_theta
surface flux proxies
```

Reuse the existing statistical validation framework where practical.

Report whether differences are:

```text
bitwise
roundoff-level
tolerance-level
scientifically meaningful
```

### Benchmark

Measure:

```text
solver wall time
whole timestep wall time
iterations/V-cycles
GPU memory
setup cost
steady-state solve cost
```

Use the same GPU count and decomposition as the reference.

Do not accept a solver solely because it uses fewer iterations.

# Phase 3 — Optimize the existing GMRES+FFT path

Even if MLTerrainPoisson is promising, independently optimize the existing validated/reference path.

Each item below should be a separate benchmarkable change.

## Optimization A — Remove redundant preconditioner zeroing

Current code contains approximately:

```cpp
lhs.setVal(0.);
m_2D_fft_precond->solve(lhs, rhs_tmp, m_stretched_dz_d);
```

Inspect the pinned AMReX `PoissonHybrid::solve()` implementation.

The FFT solve appears to overwrite the destination through the backward transform.

If confirmed:

1. Remove `lhs.setVal(0.)`.
2. Add a comment explaining why the destination is fully overwritten.
3. Run regression tests.
4. Benchmark.

A `360 x 360 x 200` double-precision field is about 207 MB. At approximately 19 GMRES iterations and two projections per timestep, this may remove several GB of unnecessary writes per timestep.

Measure the actual effect.

## Optimization B — Cache TerrainPoisson / FFT preconditioner setup

Current `solve_with_gmres()` constructs approximately:

```cpp
amrex::GMRES<MultiFab, TerrainPoisson> gmsolver;
TerrainPoisson tp(...);
gmsolver.define(tp);
```

on every projection.

`TerrainPoisson` constructs an `FFT::PoissonHybrid` object.

For static terrain, fixed grids, fixed BCs, and no regridding, this structure is invariant.

Investigate caching by level.

Potential design:

```text
ERF
  terrain_poisson_cache[level]
      TerrainPoisson
      FFT::PoissonHybrid
```

Invalidate/rebuild only when required:

```text
grid changes
DistributionMapping changes
terrain changes
boundary-condition configuration changes
level creation/removal
```

Be careful with ownership and lifetime of:

```text
ax
ay
az
detJ
z_phys_nd
stretched_dz
```

Profile setup before and after.

## Optimization C — Avoid scaling and restoring static metric coefficients every projection

Current terrain projection modifies coefficients before the solve and restores them afterward, conceptually:

```cpp
ax *= mf_ux / mf_uy;
ay *= mf_vy / mf_vx;
az /= (mf_mx * mf_my);

solve();

ax *= mf_uy / mf_ux;
ay *= mf_vx / mf_vy;
az *= (mf_mx * mf_my);
```

These quantities are static for static terrain and static map factors.

Investigate persistent projection coefficients:

```text
ax_proj
ay_proj
az_proj
```

computed when terrain/metrics/map factors are built.

Requirements:

1. Do not change coefficients used by other physics.
2. Do not modify original metric arrays in place every solve.
3. Support non-unity map factors.
4. Add a fast path if all relevant map factors are one.
5. Verify bitwise or roundoff agreement.

Benchmark memory cost versus runtime benefit.

## Optimization D — Precompute static terrain metric ratios/stencil geometry

The current `terrpoisson_adotx` repeatedly reconstructs terrain geometry from `z_phys_nd`.

It recomputes terms such as:

```text
h_xi
h_eta
h_zeta
h_xi / h_zeta
h_eta / h_zeta
1 / h_zeta
```

during every matrix-vector product.

For static terrain, these do not change.

Investigate precomputing a compact set of terrain metric coefficients.

Do not immediately store a complete 15-point matrix.

Prefer a compact representation that reduces:

```text
z_phys loads
subtractions
divisions
repeated metric interpolation
```

without excessive memory traffic.

Candidate:

```text
precompute face/edge metric ratios
reuse them in terrpoisson_adotx and flux kernels
```

Benchmark:

```text
kernel time
GMRES iteration time
whole projection time
GPU memory use
register use/occupancy if available
```

Explicitly consider the GPU tradeoff:

```text
more arithmetic vs more memory traffic
```

Keep only measured improvements.

## Optimization E — Optimize FFT::PoissonHybrid vertical tridiagonal solve

This is high priority.

The pinned AMReX code itself contains a TODO approximately saying:

```cpp
// TODO: explore how to optimize this.
// Maybe use cusparse.
// Maybe make z-direction the unit stride direction.
```

The GPU path creates workspace and reconstructs tridiagonal coefficients for every preconditioner call.

For this LES:

```text
static geometry
static dz
static FFT operator
```

much of this structure is invariant.

### E1. Persistent workspace

Avoid repeatedly allocating:

```cpp
FArrayBox tridiag_workspace(box,4);
```

if profiling shows measurable cost.

### E2. Precompute tridiagonal coefficients/factors

The lower/upper coefficients are fixed by `dz`.

Parts of the diagonal depend on horizontal spectral wavenumber but are also invariant for a fixed grid.

Determine whether the Thomas/LU factors can be precomputed for each spectral `(i,j)` mode.

Potential design:

```text
initialization:
  for every spectral (i,j):
      construct tridiagonal matrix
      factor once
      cache factors

preconditioner:
  FFT(rhs)
  apply cached tridiagonal solve
  inverse FFT
```

### E3. GPU memory layout

Assess whether the current `k` loop is inefficient because x is unit stride while each thread walks vertically.

Consider:

```text
transposed spectral workspace
batched tridiagonal solve
cuSPARSE/gtsv2 or suitable GPU batched tridiagonal routine
custom PCR/Thomas hybrid
```

Do not add an external dependency unless clearly justified.

### E4. Measure

Report:

```text
FFT forward time
solve_z time
FFT backward time
preconditioner total
```

before and after.

Keep AMReX-side changes isolated so they could be proposed upstream.

## Optimization F — Reduce GMRES synchronization/reduction overhead

AMReX GMRES uses two-pass unmodified Gram-Schmidt.

At about 19 iterations, dot products are roughly:

```text
2 * (1 + 2 + ... + 19) = 380
```

per projection.

Investigate whether Gram-Schmidt/global reductions are material.

### F1. Batched/fused dot products

Instead of sequential:

```text
dot(v, q0)
dot(v, q1)
...
```

investigate multiple inner products in one reduction.

### F2. Benchmark BiCGStab

AMReX has a generic right-preconditioned `BiCGStab` compatible with the same operator interface.

Create an experimental runtime option if practical:

```ini
erf.terrain_krylov_solver = gmres
erf.terrain_krylov_solver = bicgstab
```

Default remains `gmres`.

Compare:

```text
iterations
operator applications
preconditioner applications
reductions
wall time
final residual
post-projection divergence
```

Do not assume fewer vectors means faster convergence.

### F3. Restart length

Current ERF uses:

```cpp
gmsolver.setRestartLength(50);
```

Current solves converge near 19 iterations, so restart is not triggered.

Do not spend much time tuning restart length unless profiling identifies a real cost.

## Optimization G — Warm-start / temporal recycling

Current projection initializes:

```cpp
phi_lev.setVal(0.0);
```

and GMRES starts from zero.

Investigate using:

```text
previous RK-stage phi
previous timestep phi
scaled previous pressure increment
```

as an initial guess.

AMReX supports:

```cpp
gmsolver.setInitialGuessNonzero(true);
```

But relative tolerance is relative to the initial residual when using a nonzero guess.

Therefore simply enabling a warm start may unintentionally tighten the absolute convergence target.

Design the experiment so the final numerical accuracy is comparable to the zero-guess `1e-8` reference.

Possible approach:

1. Record the reference RHS norm.
2. Convert the desired convergence threshold to an absolute tolerance.
3. Use the warm initial guess while preserving the same effective stopping criterion.

Measure whether iterations decrease substantially.

Do not accept this optimization if the stopping criterion becomes less strict.

## Optimization H — Cache per-projection workspace

Inspect allocations/definitions inside `project_momenta()`.

Candidates:

```text
rhs
phi
fluxes
subdomain MultiFabs
temporary alias structures
BoxArray
DistributionMapping
index maps
```

Determine which can be safely persistent per level.

Prioritize only if profiling shows setup/allocation cost is material.

# Phase 4 — Multi-GPU scaling diagnosis

If the terrain case scales poorly beyond one GPU, determine whether the main cause is:

```text
FFT transpose/all-to-all communication
GMRES reductions
FillBoundary communication
domain decomposition
spectral decomposition
GPU under-occupancy
```

Run at least:

```text
1 GPU
2 GPUs
4 GPUs
```

using the same total problem size.

Record:

```text
whole timestep
projection
preconditioner
FFT forward/backward
GMRES reductions
operator apply
MPI communication if measurable
```

Calculate speedup and parallel efficiency.

Do not change science settings between GPU-count tests.

# Phase 5 — Validation requirements

Every accepted optimization must satisfy numerical validation.

Strongest goal:

```text
bitwise identical
```

when the mathematical algorithm is unchanged.

If not bitwise identical, explain why and quantify differences.

At minimum compare:

```text
solver residual
post-projection divergence
u
v
w
theta
pressure
resolved TKE
SGS TKE
uw
vw
total_uw
total_vw
w_theta
total_w_theta
surface quantities
canopy-region statistics
terrain-conditioned statistics
```

For longer statistical comparisons, reuse the existing analysis scripts and methodology.

Do not change these while benchmarking implementation changes unless explicitly required:

```text
Poisson tolerance
physical parameterizations
terrain smoothing
canopy properties
time step
advection scheme
SGS closure
boundary conditions
```

# Phase 6 — Benchmark methodology

For each candidate:

1. Build Release.
2. Record executable SHA256.
3. Start from the same checkpoint if possible.
4. Run the same simulated interval.
5. Discard initialization/warm-up from performance averages.
6. Run enough steps to reduce noise.
7. Prefer multiple repetitions where affordable.
8. Record:
   - mean s/step
   - median s/step
   - standard deviation
   - total projection time
   - solver iterations
   - time per solver iteration
   - GPU count
   - MPI ranks
   - OpenMP threads
   - node type
9. Compare against the untouched `1e-8` baseline.

Create:

| Variant | s/step | Projection s/step | Median iterations | Speedup | Numerical status |
|---|---:|---:|---:|---:|---|
| baseline | 0.8113 | ... | 19 | 1.000x | reference |
| no precond zero | ... | ... | 19 | ... | bitwise |
| cached FFT setup | ... | ... | 19 | ... | bitwise |
| cached metrics | ... | ... | 19 | ... | roundoff |
| warm start | ... | ... | ... | ... | validated |
| BiCGStab | ... | ... | ... | ... | validated |
| AMReX MLTerrainPoisson | ... | ... | ... | ... | validated |

# Phase 7 — Decision criteria

Prioritize:

1. Numerical correctness
2. Whole-LES wall-clock improvement
3. Projection wall-clock improvement
4. Robustness
5. Maintainability
6. Memory overhead
7. Code complexity

Guideline:

```text
<1% whole-run improvement:
    probably not worth merging unless extremely simple

1–3%:
    keep only if low-risk and simple

3–10%:
    worthwhile

>10%:
    high-value

>25%:
    major improvement
```

# Important constraints

Do not:

- Relax `erf.poisson_reltol` below the `1e-8` reference for the primary optimization.
- Modify physical parameterizations to get a faster result.
- Change the time step to make the benchmark look faster.
- Edit production outputs.
- Combine unrelated optimizations before benchmarking separately.
- Trust iteration counts alone; use wall-clock timing.
- Assume a new solver is better without validating LES fields.
- Remove correctness-related synchronization or boundary fills without dedicated decomposition/MPI tests.

# Relevant existing ERF work

Inspect:

```text
ERF PR #4189
Multigrid solver for the anelastic projection on terrain-fitted meshes
```

Lessons:

- An earlier ERF-specific MLMG terrain solver reduced iteration counts but was slower in several terrain tests.
- FFT-preconditioned GMRES remained faster in those tests.
- Caching multigrid setup alone gave little improvement.
- Large costs included:
  - column smoother
  - repeated 15-point stencil evaluation
  - bottom-solver operator applications
  - ghost filling
- The existing FFT preconditioner is already effective for terrain-following meshes.

Do not blindly reimplement PR #4189.

Instead compare against the newly merged native AMReX `MLTerrainPoisson`.

# Recommended order of work

```text
1. Baseline fine-grained profile

2. Test native AMReX MLTerrainPoisson
   - separate branch
   - exact 1e-8 comparison

3. Existing GMRES+FFT:
   a. remove redundant lhs.setVal(0)
   b. cache TerrainPoisson/PoissonHybrid setup
   c. cache projection metric coefficients
   d. profile and optimize PoissonHybrid::solve_z
   e. precompute compact terrain metric coefficients

4. GMRES algorithm:
   a. profile Gram-Schmidt reductions
   b. test fused reductions if worthwhile
   c. benchmark BiCGStab
   d. test warm start with equivalent convergence criterion

5. Multi-GPU scaling

6. Combine only individually validated winners
```

# Deliverables

Produce:

```text
terrain_projection_optimization/
├── README.md
├── baseline.md
├── profiles/
│   ├── baseline_profile.md
│   └── optimized_profile.md
├── benchmarks/
│   ├── benchmark_summary.csv
│   ├── benchmark_summary.md
│   └── raw/
├── validation/
│   ├── numerical_validation.md
│   └── statistical_validation.md
└── notes/
    ├── amrex_mlterrainpoisson.md
    ├── gmres_fft_optimization.md
    └── rejected_ideas.md
```

Final report must contain:

## 1. Executive summary

State:

```text
best implementation
whole-run speedup
projection speedup
numerical impact
recommended production configuration
```

## 2. Bottleneck breakdown

Give the percentage of terrain timestep spent in:

```text
GMRES
FFT
vertical tridiagonal solve
terrain operator
Gram-Schmidt/reductions
boundary communication
other projection work
non-projection LES physics
```

## 3. Optimization results

For each attempted optimization:

```text
hypothesis
code change
performance result
numerical result
keep/reject decision
```

## 4. Recommended patch set

List exact commits that should be retained.

## 5. Remaining opportunities

Rank future work by expected payoff.

# Git workflow

Use small commits.

Examples:

```text
perf: add detailed terrain projection profiling

perf: remove redundant FFT preconditioner zeroing

perf: cache terrain FFT preconditioner setup

perf: cache terrain projection metric coefficients

perf: prototype warm-started terrain GMRES

perf: add BiCGStab terrain projection experiment

perf: test AMReX MLTerrainPoisson for terrain projection

perf: optimize PoissonHybrid vertical tridiagonal solve
```

Do not squash until all results are documented.

For every performance commit, record the benchmark result in the commit message or optimization report.

# Final success target

The target is not merely fewer GMRES iterations.

The target is:

```text
substantially lower terrain LES wall time
at the original 1e-8 projection accuracy
with scientifically equivalent results
and maintainable code
```

Given that the terrain projection dominates the current runtime, focus on that path first.
