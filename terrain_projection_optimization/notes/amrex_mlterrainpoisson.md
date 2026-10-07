# AMReX MLTerrainPoisson candidate

## Source and scope

AMReX commit 5e1ef8467df5f067be4d673c6efe45e1ffb3c9b1 is titled
“MLTerrainPoisson: cell-centered terrain-following Poisson operator (#6030)”.
It adds an MLMG cell-centered operator for div(A grad phi)/J on a 3-D
terrain-following mesh. It supports one AMR level covering the domain and
homogeneous periodic, Neumann, and Dirichlet boundary conditions. It offers
MLMG and GMRESMLMG on CPU and GPU. This is a multigrid candidate, not a runtime
switch to unpreconditioned GMRES.

The AMReX commit message reports a WPS 200 x 200 x 176 case on an RTX 5070:
ERF GMRES+FFT 3.31 s / 53 iterations, MLMG 0.112 s / 3 V-cycles, and
GMRESMLMG 0.184 s / 3 iterations. These are upstream results on another
hardware/domain and are not assumed to transfer to the current H100 case.

## ERF integration contract

The API takes nodal physical height via setZPhys, face terrain areas via
setAreas, cell volume factor via setDetJ, and the same homogeneous domain
boundary types. ERF's make_areas produces the face factors expected by the new
operator. The experimental path must pass the unmodified area factors and
detJ, not the temporarily map-factor-scaled ax/ay/az arrays used in the
existing GMRES path.

The ERF projection RHS is reused. MLMG internally solves the J-scaled system;
verify its residual convention and pressure-gradient flux sign before applying
the correction to momentum. The original GMRES+FFT implementation remains
the default and remains selectable as the reference.

## Boundary-condition risk

The AMReX commit notes that Neumann domain-face flux is set to zero. ERF's
existing ghost-cell reflection can retain terrain cross-term flux at a
sloping Neumann face. Therefore boundary columns can differ. Compare the
operator and actual post-projection divergence on the full domain and inspect
boundary-column contributions; report any resulting physical difference
explicitly rather than calling it roundoff.

## Current case compatibility

The target LES is a 3-D, single-level, full-domain terrain solve with periodic
horizontal directions and a vertical boundary treatment within the documented
BC set. Confirm grid coverage, vertical box splitting, map factors, and exact
boundary names in the runtime preflight before benchmarking.
