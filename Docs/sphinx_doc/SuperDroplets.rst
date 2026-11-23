.. _sec:SuperDroplets:

Super-Droplets
==============

The super-droplet method represents the condensed phase with Lagrangian
particles rather than with bulk mixing ratios.  Each particle carries a
multiplicity -- the number of real droplets it stands for -- together with the
mass of each condensate species and of each aerosol it contains.  The method is
selected with ``erf.moisture_model = SuperDroplets``; see
:ref:`Microphysics` for where it sits among the other moisture models, and
:ref:`sec:Plotfile3DReference` for the per-species plotfile names it provides.

All of the inputs below use the prefix ``super_droplets_moisture``, which is the
name of the particle container (``ERF_SuperDropletsMoist.H``).  Inputs that the
container inherits from the generic ERF particle container use the same prefix;
they are listed under :ref:`inputs-particles` for tracer particles and behave
the same way here.

.. _superdroplets-anelastic-compatibility:

Anelastic compatibility
-----------------------

``SuperDroplets`` uses a thermodynamic pressure source selected independently
on each AMR level. On a level with ``erf.compressibility_type`` not ``FullyCompressible``, the pressure supplied
to the super-droplet microphysics is the hydrostatic base-state pressure at
each cell, rather than a pressure diagnosed from the compressible equation of
state. ERF diagnoses the temperature from dry potential temperature using

.. math::

   \theta = \frac{\rho_d\theta}{\rho_d}, \qquad
   p_{\mathrm{SDM}} = p_{\mathrm{base}}, \qquad
   T = \theta\left(\frac{p_{\mathrm{base}}}{p_{\mathrm{ref}}}\right)^{R_d/c_p}.

Here :math:`\rho_d` is dry-air density, :math:`p_{\mathrm{base}}` is the
cellwise hydrostatic base-state pressure, :math:`p_{\mathrm{ref}}=10^5` Pa is
the fixed reference pressure used to define potential temperature, and
:math:`c_p` is the configured ``erf.c_p``. These pressures are in pascals.
No additional pressure or Exner-function input is required.

During an anelastic SDM microphysics source step, the pressure remains the
base-state pressure while the temperature is refreshed from the updated
potential temperature after each phase-change substep. Latent heating uses
the fixed inverse Exner factor for every process, including sequential
processes within a substep:

.. math::

   \Delta\theta =
   \left(\frac{p_{\mathrm{ref}}}{p_{\mathrm{base}}}\right)^{R_d/c_p}
   \frac{L}{c_p}\,\Delta q.

Here :math:`\Delta q` is the signed, dry-air-referenced mass transfer
and :math:`L` is the corresponding latent heat. This pressure and temperature
are supplied consistently to the scheme's saturation, phase-change,
terminal-velocity, and collision calculations. On compressible levels, SDM
diagnoses pressure and temperature from the local equation of state and uses
the local potential-temperature-to-temperature ratio for latent heating.

For example, with a particles-enabled build:

.. code-block:: text

   erf.compressibility_type = Anelastic
   erf.moisture_model = SuperDroplets

Choose compatible mesh, solver, particle, and AMR options, and validate the
selected configuration before using results for quantitative cloud-microphysics
analysis.

.. warning::

   Anelastic ``SuperDroplets`` execution has not been validated on CUDA or HIP
   GPUs. Check the selected particle processes, restart behavior, and AMR
   configuration before using GPU results for scientific analysis; do not
   assume reproducibility between CPU and GPU runs.

Species transport and diffusion
-------------------------------

In coupled atmospheric simulations, each additional species specified by
``super_droplets_moisture.species`` has an Eulerian vapor mixing ratio and a
condensate mixing ratio diagnosed from the super-droplets.
ERF advects the vapor with the resolved flow and applies the
configured turbulent diffusivity. Vertical turbulent diffusion includes an
implicit contribution when ``erf.implicit_moisture_diffusion`` is enabled and
the stage factor in ``erf.vert_implicit_fac`` is nonzero. These additional
species have no molecular scalar diffusion, even when molecular diffusion is
enabled for water or passive scalars.

The condensate field is a diagnostic of particle-carried mass, not a second
mass reservoir. It is not independently advected, diffused, or positivity-clipped
by the atmospheric solver. Particle motion and microphysical processes determine
its evolution. In particular, vapor diffusion alone does not represent turbulent
dispersion of condensed particles. These rules apply to the additional species;
the existing water microphysics coupling is unchanged.

For vertical implicit diffusion, the supported species boundary conditions are
zero-gradient extrapolation and prescribed Neumann gradients. If a transported
species has an unsupported vertical boundary condition, ERF issues a warning
and disables the shared moisture implicit-diffusion option. Horizontal diffusion,
explicitly treated stages, and terrain mixed terms can still constrain the
timestep. An implicit vertical solve is not an unconditional positivity or
stability guarantee for the full time integration.

AMR and particle sampling
~~~~~~~~~~~~~~~~~~~~~~~~~

Particles are assigned to the finest available AMR level containing their
positions. Newly injected particles are assigned before phase change or
condensate projection, including when ``include_advection = false``.
This assignment does not require ``split_merge_amr = true``.

Initial populations are seeded only in cells not covered by finer grids, so
refining a region does not add another copy of its prescribed particle mass.
Density scaling and condensate-based attribute initialization are completed
on each level before redistribution. If initial phase-change relaxation is
enabled, it then runs on each particle's owning level and retains the associated
vapor and latent-heat changes in the initial atmospheric state. These
initialization operations are not repeated when restarting from a checkpoint.

A restart restores the condensate fields from the checkpoint, including
fine-level averages in covered coarse cells. The ground accumulations
``rain_accum``, ``snow_accum`` and ``precip_accum`` are checkpointed on every
level. They and the per-species and per-aerosol ground accumulations
(``accum_<species>``, ``accum_<aerosol>``) are kept when the grids of a level
change, a new fine cell taking the depth, or an equal share of the amount, of
the coarse cell under it. With one-way coupling (``erf.coupling_type = OneWay``)
the coarse level does not receive what reached the ground on a fine level, so
that amount is lost where the fine level moves away or is removed. The
per-species and per-aerosol accumulations are not checkpointed: they start again
from zero after a restart. Turbulent diffusivities, including ``Kmv``, are not
checkpointed. They can be zero in the initial restart plotfile until the next
turbulence update, even when the conserved state has been restored correctly.

For additional species, cloud-in-cell projection uses only cells on the
particle's level that are not covered by finer grids. The retained weights are
renormalized, so projection at a physical wall or coarse/fine boundary preserves
the particle's represented mass. Near these boundaries the stencil is one-sided;
it does not transfer condensate mass between levels independently of particle
motion. Diagnostic condensate is refreshed after initialization, time-step
redistribution, and regridding without advancing phase change or vapor transport.
For water, ice and the other projected quantities, the share of the stencil
outside a non-periodic domain boundary is added to the adjacent boundary cell.

With ``split_merge_amr = true``, ERF also adjusts the super-droplet sampling as
particles change levels. When empty fine cells are populated from a neighboring
donor, the donor and its daughters share the original multiplicity. Each retained
daughter represents at least one real droplet; insufficient multiplicity leaves
some cells empty rather than creating physical mass. The allocation uses a fixed
cell ordering for a fixed particle population and grid decomposition. This is
not a guarantee of identical particle distributions across decompositions.

Condensate output and mass budgets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The output fields ``qv_<species>`` and ``qc_<species>`` are mixing ratios per unit
dry-air mass (:math:`\mathrm{kg\,kg^{-1}}`), whereas
``super_droplets_moisture_species_mass_density_<species>`` is particle-derived
condensate mass per physical volume (:math:`\mathrm{kg\,m^{-3}}`). Thus the
condensate density corresponding to ``qc_<species>`` is
:math:`\rho q_{c,s}`, not :math:`q_{c,s}` alone; :math:`\rho` is ERF's dry-air
density (the ``density`` output field).

For a three-dimensional mass budget, sum only uncovered cells of the AMR
hierarchy, using their physical volumes :math:`V_c`:

.. math::

   M_s = \sum_{\mathrm{composite\ cells}} V_c
         \left(\rho q_{v,s} + \rho q_{c,s}\right),
   \qquad
   M_{c,s} = \sum_{\mathrm{active\ particles}} \xi_p m_{p,s}.

Here :math:`\xi_p` is the number of real droplets represented by super-droplet
:math:`p`, and :math:`m_{p,s}` is the mass of species :math:`s` in one real
droplet. The particle sum is an independent measure of condensate mass:
do not add it to the mesh condensate term, which represents the same mass.
Phase change exchanges vapor and condensate. Account separately for injection,
precipitation, and other physical sources, sinks, or boundary fluxes when
assessing conservation; an injection run is not a closed system.

With ``dimensionality = one_d_z``, species condensate density from the reference
column (the lowest horizontal indices in each grid box) is replicated
horizontally before conversion to mixing ratio. Both normal microphysics
updates and diagnostic refreshes use this convention. It is not a horizontal
average: a volume integral over the replicated field must be interpreted using
the represented column area, rather than compared directly with an unreplicated
particle-column sum. The ``three_d``, ``two_d_xz``, and ``two_d_yz`` species
projections do not apply this extra horizontal replication.

Inflow limitations
~~~~~~~~~~~~~~~~~~

There is no per-additional-species prescribed inflow interface. Ordinary inflow
with a specified density uses a zero external value for these species;
``qv`` and ``qc`` inflow settings prescribe water only.
When inflow density is not prescribed, the existing fallback extrapolates
non-ingested conserved fields from the interior. Boundary-plane ingestion
does not supply additional species. Do not interpret water-vapor boundary
settings as non-water inflow concentrations.

Model configuration
-------------------

These inputs select which processes the model includes and what it carries.
They are read by the moisture model itself (``ERF_SuperDropletsMoistInit.cpp``)
rather than by the particle container, but use the same prefix.

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **species**                                     | additional condensate species carried besides water, and | List of species names                            | none                                |
|                                                 | ice when cold processes are on.  Naming a species that   |                                                  |                                     |
|                                                 | is already modeled aborts the run, and at most five      |                                                  |                                     |
|                                                 | extra species may be given                               |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **aerosols**                                    | aerosol materials carried by the super-droplets          | List of species names                            | none                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **include_phase_change**                        | include condensation and evaporation                     | Boolean                                          | true                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **include_advection**                           | advect the super-droplets                                | Boolean                                          | true                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **include_coalescence**                         | include coalescence                                      | Boolean                                          | true                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **include_cold_processes**                      | include ice; adds ice to the species list                | Boolean                                          | true                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **distribution_type**                           | initialization of the particles                          | uniform                                          | uniform                             |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **radius_raindrop**                             | radius [m] above which a drop is counted as rain rather  | Real > 0                                         | 4.0e-5                              |
|                                                 | than cloud water                                         |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **rime_mass_ratio**                             | rime mass fraction above which a particle is treated as  | Real                                             | 0.3                                 |
|                                                 | rimed                                                    |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **kinematic_mode**                              | run with a prescribed flow field rather than a coupled   | Boolean                                          | false                               |
|                                                 | one                                                      |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **dimensionality**                              | dimensionality of the super-droplet dynamics             | one_d_z, two_d_xz, two_d_yz, three_d             | three_d                             |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **recycle_particles**                           | return deactivated particles to the domain instead of    | Boolean                                          | false                               |
|                                                 | dropping them                                            |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **num_substeps_phase_change**                   | sub-steps taken per time step for the phase-change       | Integer > 0                                      | 1                                   |
|                                                 | update                                                   |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **initial_phase_change_relaxation**             | relax the particles to an equilibrium size at            | Boolean                                          | false                               |
|                                                 | initialization                                           |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **initial_phase_change_relaxation_time**        | duration [s] of that relaxation                          | Real > 0                                         | 10.0                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **diagnostics_interval**                        | number of steps between writes of the distribution       | Integer > 0                                      | 1                                   |
|                                                 | diagnostics                                              |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+

The species and aerosol names given here determine the per-material keys used by
the initialization and injection blocks below, and the ``qv_<species>`` and
``qc_<species>`` plotfile names.

Solver and time integration
---------------------------

The condensational growth of each particle is an ODE in the droplet mass, solved
per particle per step.

The implicit Newton solve uses the larger of ``newton_solver_stol`` and the
particle arithmetic's machine epsilon as its relative step stopping criterion.
This avoids rejecting a physically resolved update when its remaining Newton
correction is at roundoff. Requesting a tolerance below machine epsilon does
not increase the arithmetic precision of the solve.

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **mass_change_ti_method**                       | time integrator for the mass-change ODE.  An             | backward_euler, crank_nicolson, dirk2, rk3bs,    | backward_euler                      |
|                                                 | unrecognized value aborts the run                        | rk4                                              |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **mass_change_cfl**                             | CFL-like limit on the mass-change sub-step               | Real > 0                                         | 1000.0                              |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **newton_solver_rtol**                          | relative tolerance of the Newton solve for the implicit  | Real > 0                                         | 1.0e-6                              |
|                                                 | integrators                                              |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **newton_solver_atol**                          | absolute tolerance of that solve                         | Real >= 0                                        | 1.0e-99 (0.0 in single precision)   |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **newton_solver_stol**                          | step tolerance of that solve                             | Real > 0                                         | 1.0e-12                             |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **newton_solver_maxits**                        | maximum Newton iterations per particle                   | Integer > 0                                      | 10                                  |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **mass_change_unconverged_log**                 | write a log entry for every particle whose mass-change   | Boolean                                          | false                               |
|                                                 | solve fails to converge                                  |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **mass_change_unconverged_log_filename**        | name of that log file                                    | String                                           | unconverged_superdroplets.log       |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+

Ventilation
-----------

The ventilation factor enhances the mass-change rate of a falling drop.  It
follows Bayley et al. (2025), Eq. 3, and is off by default.

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **mass_change_ventilation**                     | include the ventilation factor in the phase change       | Boolean                                          | false                               |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **ventilation_alpha1**                          | fit coefficient alpha_1                                  | Real                                             | 6.954e7                             |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **ventilation_beta1**                           | fit exponent beta_1                                      | Real                                             | 1.963                               |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **ventilation_alpha2**                          | fit coefficient alpha_2                                  | Real                                             | 1.069e3                             |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **ventilation_beta2**                           | fit exponent beta_2                                      | Real                                             | 0.702                               |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **ventilation_fcap**                            | upper bound applied to the ventilation factor            | Real > 0                                         | 20.0                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+

Coalescence
-----------

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **coalescence_kernel**                          | collision kernel used for coalescence.  An unrecognized  | sedimentation, golovin, Longs, Halls, geometric  | sedimentation                       |
|                                                 | value aborts the run.  Note that ``Longs`` and           |                                                  |                                     |
|                                                 | ``Halls`` are matched case-sensitively;                  |                                                  |                                     |
|                                                 | ``geometric`` is ``sedimentation`` with collision        |                                                  |                                     |
|                                                 | efficiency 1                                             |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **kernel_relative_velocity**                    | which relative velocity enters the kernel                | terminal_velocity, absolute_velocity,            | terminal_velocity                   |
|                                                 |                                                          | radial_velocity                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **include_brownian_coalescence**                | add the Brownian contribution to the kernel              | Boolean                                          | false                               |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **ice_aggregation_efficiency**                  | collection efficiency prefactor for ice-ice aggregation  | Real                                             | 0.1                                 |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **coalescence_bin_size**                        | size, in cells, of the bins within which candidate       | 3 Integers > 0                                   | 1 1 1                               |
|                                                 | collision pairs are drawn                                |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **terminal_velocity_model**                     | terminal-velocity model.  One value sets the water       | RogersYau, AtlasUlbrich, CloudRainShima,         | CloudRainShima IceBohm              |
|                                                 | model; a second, if given, sets the ice model.  Giving   | IceBohm                                          |                                     |
|                                                 | more than two values aborts the run                      |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+

Particle population and motion
------------------------------

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **nucleate_particles**                          | create new super-droplets from the vapor field           | Boolean                                          | false                               |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **density_scaling**                             | scale the initial number density with the local air      | Boolean                                          | false                               |
|                                                 | density                                                  |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **place_randomly_in_cells**                     | place new particles at random positions within each cell | Boolean                                          | true                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **prescribed_advection**                        | move the particles with a prescribed vertical velocity   | Boolean                                          | false                               |
|                                                 | instead of the resolved flow                             |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **prescribed_w**                                | the prescribed updraft [m/s] used in that mode; a value  | Real                                             | 0.0                                 |
|                                                 | of zero selects a sinusoidal pulse instead of a constant |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **split_merge_amr**                             | split and merge particles as they cross AMR level        | Boolean                                          | false                               |
|                                                 | boundaries                                               |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **sigma0**                                      | width parameter of the assumed mass distribution         | Real > 0                                         | 0.62                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+

Deactivation and recycling
--------------------------

A particle whose multiplicity falls below a threshold is deactivated and may be
returned to the domain at a new position.

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **inactive_threshold**                          | multiplicity below which a particle is deactivated       | Real > 0                                         | 0.01                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **write_inactive_plt**                          | include deactivated particles in the particle plotfiles  | Boolean                                          | false                               |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **recycle_xmin**, **recycle_xmax**              | x-bounds of the region into which recycled particles are | Real                                             | domain bounds                       |
|                                                 | returned                                                 |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **recycle_ymin**, **recycle_ymax**              | y-bounds of that region                                  | Real                                             | domain bounds                       |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **recycle_zmin**, **recycle_zmax**              | z-bounds of that region                                  | Real                                             | domain bounds                       |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+

.. _sec:SDWallLoss:

Super-droplet wall loss
-----------------------

``super_droplets_moisture.wall_loss.enable = true`` removes represented droplets
from suspension (default ``false``). All keys below use the prefix
``super_droplets_moisture.wall_loss.``. Use a flat, uniform Cartesian grid with
``amr.max_level = 0`` and ``erf.mesh_type = ConstantDz``; stretched grids,
terrain, embedded boundaries and configured refinement are not supported.
Use ``erf.moisture_tight_coupling = false``: the sink acts once per physical
step, after optional advection and boundary crossings, before coalescence and
recycling, not once per phase-change substep. It reduces active particles'
multiplicity, not their radius or mass per real droplet. Condensate projection
reflects that loss; it is not a vapor source, latent heat flux, precipitation
accumulation or aerosol re-emission. For chamber setup and avoiding duplicate
losses, see :ref:`CloudChamberParticleLoss`.

Choosing the loss rate
~~~~~~~~~~~~~~~~~~~~~~

``deposition_velocity_model`` selects one of three models. Only nonperiodic
no-slip, slip or surface-layer walls count; open and periodic faces never
contribute. With no eligible wall there is no sink.

``constant`` (default)
   Only particles in wall-adjacent cells lose multiplicity. For each wall face
   touching the cell, :math:`k_b=v_{d,b} A_b/V_{\rm cell}` in inverse seconds;
   add **all** such rates, including at edges and corners. Interior cells have
   zero rate. ``vd_side`` (all four lateral walls), ``vd_top`` and ``vd_bottom``
   each default to ``1e-4`` m/s. Optional ``deposition_velocity`` sets all three
   first; individual keys override it. Use nonnegative speeds. For one face the
   local decay timescale is :math:`\Delta n/v_d`, where :math:`\Delta n` is the
   cell width normal to the wall, not the chamber width. For example, a 0.0625 m
   cell at the default speed gives 625 s. A *well-mixed chamber* instead has
   rate :math:`\sum_b v_{d,b} A_{b,\rm chamber}/V_{\rm chamber}`.

``diffusion`` (alias ``stokes_einstein``)
   Uses the same wall cells and face-area factors, but replaces the constant
   speeds by :math:`D_B(r)/\delta`, where the Brownian diffusivity is
   :math:`D_B=k_B T C_c/(6\pi\mu r)` in m²/s. The prescribed gas temperature
   ``diffusion_temperature`` is 293 K by default, dynamic viscosity
   ``diffusion_viscosity`` is ``1.8e-5`` Pa s, and wall-layer thickness
   ``diffusion_delta`` is ``1e-3`` m. This thickness is not the grid spacing.
   ``diffusion_slip_correction = true`` uses the Cunningham factor :math:`C_c`
   with ``diffusion_mean_free_path = 6.6e-8`` m; ``false`` sets :math:`C_c=1`.

For **both face models**, ``vd_floor`` (default ``0`` m/s) is a lower bound
on the base speed at *every* wall, not the bottom-wall speed. After applying
that bound, ``settling_correction = true`` adds the magnitude of the particle's
terminal speed at the bottom, subtracts it at the top (minimum zero), and
leaves the sides unchanged. The correction is off by default. It uses the speed
stored by particle advection, not a separate fall-speed calculation. To refresh
that speed without moving the particles, use ``include_advection = true`` with
``advect_with_flow = false`` and ``advect_with_gravity = false`` and no prescribed
motion. If advection is skipped, a newly initialized particle has zero stored speed.

``measured_rate`` (alias ``measured_profile``)
   Applies a chamber-wide first-order rate :math:`k(D_p,t)` to **every**
   suspended particle, including in the interior. There is no face-velocity
   or cell-size factor. :math:`D_p=2r` is the current effective (wet) diameter,
   not the dry aerosol injection diameter. Supply ``profile_diameters`` (1–64
   positive, strictly increasing knots) and equally many nonnegative
   ``profile_rates``. Rates interpolate linearly in *log diameter*, not log
   rate; sizes outside the table use the nearest endpoint. A single knot gives
   a size-independent rate. ``profile_diameter_units`` accepts ``m`` (default),
   ``um`` (alias ``micron``) or ``nm``. ``profile_rate_units`` accepts ``s^-1``
   (default; ``s-1`` or ``1/s``) or ``h^-1`` (``h-1`` or ``1/h``).

   The profile is constant in time unless ``profile_rates_end`` supplies a
   second set of rates on the same knots in the same units. Then specify
   ``profile_time_start`` and ``profile_time_end`` in simulation seconds, end
   greater than start. ERF blends the two rates linearly at the step midpoint,
   holding the endpoint profiles outside the time range. ``settling_correction``
   is ignored: the measured rate can already include settling.

Multiplicity, output and restart
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With loss probability :math:`P=1-\exp(-k\Delta t)`, ``method = deterministic``
(default; alias ``exponential``) leaves :math:`\xi(1-P)` of the original
multiplicity :math:`\xi`. ``stochastic`` (alias ``binomial``) approximates
binomial thinning: a Poisson draw for :math:`\xi P<30`, otherwise a Gaussian
with mean :math:`\xi P` and variance :math:`\xi P(1-P)`, bounded to
:math:`[0,\xi]`. ``deactivate_on_depletion`` defaults to ``false``. If enabled
and the remaining multiplicity falls below ``deactivate_threshold`` (default
``1e-30`` real droplets), the residual is also counted as deposited and the
super-droplet is removed, not recycled.

``diag_interval`` sets the cadence for ``super_droplets_moisture_wall_loss.txt``.
The default (100) samples steps 100, 200, and so on; ``0`` disables logging,
not deposition. Each row reports domain-wide losses, not fluxes per unit area:

#. **Columns 1–2:** start time and duration :math:`\Delta t` of the sampled step (s).
#. **Columns 3–5:** represented droplet number lost on **that step** to sides,
   top and bottom, respectively; not the number of computational particles.
#. **Columns 6–8:** cumulative **total material mass** (kg) deposited on those
   same wall groups through the end of the step. Every step contributes removed
   multiplicity times the mass of one droplet, including aerosols and all
   condensates, not just liquid water.
#. **Columns 9 onward:** number lost on all walls in each **radius** bin divided
   by the sampled step's :math:`\Delta t` (droplets/s), smallest radius first.
   ``num_bins`` defaults to 16 (1–64); ``bin_rmin`` and ``bin_rmax`` default to
   ``1e-8`` and ``1e-4`` m. Bins are log-spaced with outlying radii assigned to
   the end bins. These are not the measured-profile diameter knots.

Face models attribute loss by :math:`k_b/k`; the measured model uses each wall
group's fraction of the **absorbing chamber area** for reporting only. The bin
rates summed and multiplied by that row's :math:`\Delta t` recover its number
loss. At cadence greater than one, neither number nor bin rates cover skipped
steps. Difference the cumulative mass columns for an interval mass, or use
``diag_interval = 1`` for a complete number history.

Native checkpoints preserve cumulative masses even between log writes. Restart
with compatible wall-loss settings: in the same run directory, the log must
end at the checkpoint's last recorded row and is appended without duplicating
it. In a fresh output directory the new log continues from checkpoint totals,
not zero. A stale, ahead or incompatible log is rejected rather than overwritten.
A checkpoint without valid wall-loss state cannot continue a loss history;
older checkpoints remain usable with wall loss off. Continuity of the totals
does not imply identical stochastic draws across restarts.

If an interrupted write leaves a ``.lock`` directory, verify that no writer is
active before inspecting or removing it. If the log is damaged or ahead of the
checkpoint, restart in a fresh output directory instead of truncating the log.

Diagnostics
-----------

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **distribution_grid_size**                      | number of points in the one-dimensional radius grid on   | Integer > 0                                      | 100                                 |
|                                                 | which size distributions are accumulated                 |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **distribution_rmin**                           | smallest radius [m] of that grid.  Read only when built  | Real > 0                                         | 1.0e-6                              |
|                                                 | with ``ERF_USE_ML_UPHYS_DIAGNOSTICS``                    |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **distribution_rmax**                           | largest radius [m] of that grid.  Read only when built   | Real > 0                                         | 5.0e-3                              |
|                                                 | with ``ERF_USE_ML_UPHYS_DIAGNOSTICS``                    |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+

Local droplet size distribution probes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A named probe counts liquid droplets inside a fixed sampling volume. It writes
one whitespace-delimited file, ``<name>.dsd``, in the run directory: one row per
sampled time. It needs SDM but not the optional ML diagnostics build or particle
plotfiles. No probe runs unless named. For example, for a 2 x 2 x 1 meter chamber:

::

   super_droplets_moisture.dsd_probe_names    = S03 CellA
   super_droplets_moisture.dsd_probe_interval = 40

   super_droplets_moisture.S03.center       = 1.0 1.0 0.5
   super_droplets_moisture.S03.side_lengths = 0.125
   super_droplets_moisture.S03.rmin         = 0
   super_droplets_moisture.S03.rmax         = 50.e-6
   super_droplets_moisture.S03.nbins        = 50

   super_droplets_moisture.CellA.cell       = 4 4 4
   super_droplets_moisture.CellA.r_edges    = 0 2.e-6 5.e-6 10.e-6 20.e-6
   super_droplets_moisture.CellA.write_mass = true

Give **one** geometry per probe: ``cell`` is a level-0 cell index (including
its refined descendants); ``box_lo`` and ``box_hi`` are opposite corners in physical
model coordinates; or give ``center`` and ``side_lengths`` (one length for a cube,
three for a cuboid). Boxes include the lower but not the upper faces. Each owned
particle is tested at its position, without deposition or interpolation. On a
terrain-fitted mesh, box heights are physical, not computational or above ground.
The normalizing volume is the **full** box volume, not clipped to the terrain; put
the whole box in the air for an air-volume concentration. An indexed cell uses
its physical volume, including its local layer thickness.

Specify strictly increasing ``r_edges``, or ``rmin``, ``rmax``, ``nbins`` and
``spacing`` (``linear``, the default, or ``log`` with positive ``rmin``). All
lengths and radii are in meters. Only active particles with positive liquid H2O
mass contribute. Radius is the liquid-water-equivalent radius
:math:`r=(3m_w/(4\pi\rho_w))^{1/3}`, not the diameter or total wet radius
including aerosol and ice. Out-of-range radii are omitted.

The header lists the edges and column units. After time in seconds, column
:math:`N_b=\sum_{r_b\le r<r_{b+1}}\xi/V` is physical droplet number per cubic
meter **in** bin :math:`b` (multiplicity :math:`\xi`, sample volume :math:`V`).
Optional mass columns follow the number columns and contain :math:`\sum\xi m_w/V`
in kg per cubic meter. Empty bins are zero. Divide by the linear or logarithmic
bin width for a differential spectrum.

The interval is in coarse timesteps (default 1), independent of
``super_droplets_moisture.diagnostics_interval``. At a fixed timestep of 0.025 s,
the example samples at 1 Hz. Set ``diagnostics_interval`` to ``0`` for probes
without the usual whole-domain distribution files. A fresh run samples step 0;
a restart preserves rows before the checkpoint time and replaces later rows,
rejecting an incompatible header. A new probe name starts a new series.

.. _sec:SuperDropletInitialization:

Initialization blocks
---------------------

The initial particle population is described by one or more numbered blocks.
``super_droplets_moisture.num_initializations`` sets how many there are, and
block ``i`` is read under the prefix ``super_droplets_moisture.<i>`` for
``i = 0 ... num_initializations-1``. A block-specific value overrides the same
key under ``super_droplets_moisture``; otherwise that common value is used.

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **num_initializations**                         | how many initialization blocks are read                  | Integer >= 0                                     | 1                                   |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **<i>.initial_distribution_type**               | shape of the seeded region                               | uniform, bubble                                  | uniform                             |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **<i>.initial_particles_per_cell**              | super-droplets created per cell in that region           | Integer > 0                                      | set from the number density         |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **particle_box_lo**, **particle_box_hi**        | corners of the seeded box; read when                     | 3 Reals                                          | domain corners                      |
|                                                 | ``initial_distribution_type`` = ``uniform``              |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **particle_bubble_center**,                     | center and radii of the seeded bubble; read when         | 3 Reals                                          | none                                |
| **particle_bubble_radius**                      | ``initial_distribution_type`` = ``bubble``               |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **initial_number_density**                      | number density of real droplets represented              | Real > 0                                         | none                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **initial_super_droplet_density**               | number density of super-droplets used to represent them  | Real > 0                                         | none                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **maximum_multiplicity**                        | cap on the multiplicity assigned to a super-droplet      | Real > 0                                         | none                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **multiplicity_type**                           | how multiplicity is distributed across the population    | see ``ERF_SDInitialization.cpp``                 | none                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **ice_apparent_density**                        | apparent density [kg/m^3] used for ice particles         | Real > 0                                         | none                                |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+

Each material has its own size or mass distribution. For example,
``initial_species_distribution_type_<species> = mass_constant`` uses
``initial_species_mean_mass_<species>`` in kg per real droplet, not per
super-droplet. Other supported distributions are ``mass_exponential``,
``radius_log_normal``, and ``radius_lognormal_autorange``. Their parameters use
the same prefix with ``min_mass``, ``max_mass``, ``min_radius``, ``mean_radius``,
or ``max_radius`` in place of ``mean_mass``, as appropriate. Radii are in meters.
For log-normal radii, specify either ``initial_species_std_radius_<species>``
(the standard deviation of log radius) or
``initial_species_geomstd_radius_<species>`` (the geometric standard deviation),
not both. Aerosols use ``initial_aerosol_`` in place of ``initial_species_``.
Material names must match the configured species and aerosol names.

Injection blocks
----------------

Particles may also be injected during the run.
``super_droplets_moisture.num_injections`` sets how many injection blocks are
read, and block ``i`` uses the prefix ``super_droplets_moisture.injection.<i>``.
Block-specific values override common injection values under
``super_droplets_moisture.injection``. The rate is a volumetric number source
within the injection region, not a total particle count per second for the domain.
The region may lie inside fine coverage or span a coarse/fine boundary;
particles are injected once and then assigned to their finest available level.

+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| Parameter                                       | Definition                                               | Acceptable Values                                | Default                             |
+=================================================+==========================================================+==================================================+=====================================+
| **num_injections**                              | how many injection blocks are read                       | Integer >= 0                                     | 0                                   |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **rate**                                        | real-droplet number source [m^-3 s^-1]                   | Real > 0                                         | 0 (no injection)                    |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **sd_rate**                                     | super-droplet number source [m^-3 s^-1]                  | Real > 0                                         | derived from rate and multiplicity  |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **min_multiplicity**                            | target minimum multiplicity in per-box injection mode    | Real > 0                                         | 100                                 |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **t_start**, **t_stop**                         | simulation times [s] between which injection is active   | Real                                             | whole run                           |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **domain_velocity**                             | translation velocity [m/s] of the injection region       | 3 Reals                                          | 0.0 0.0 0.0                         |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+
| **fractional_tol**                              | tolerance used when converting a fractional injection    | Real > 0                                         | 1.0e-3                              |
|                                                 | count into whole particles                               |                                                  |                                     |
+-------------------------------------------------+----------------------------------------------------------+--------------------------------------------------+-------------------------------------+

Within an injection block, use ``distribution_type = uniform`` for a box
specified by ``particle_box_lo`` and ``particle_box_hi``, or ``bubble`` with
``particle_bubble_center`` and ``particle_bubble_radius``. Material distribution
keys follow the initialization pattern without ``initial_``; for example,
``species_distribution_type_<species>`` and ``species_mean_mass_<species>``.

Specifying ``particles_per_cell`` selects per-cell injection. Otherwise, ERF
accumulates a super-droplet count over the injection region using ``sd_rate``
or ``rate / min_multiplicity``. If both are explicitly specified, the smaller
super-droplet rate is used. Physical-particle counts and integer multiplicities
are accumulated or rounded according to the selected mode, so low rates need
not produce an injection event every time step. Include the actual represented
injected mass in conservation diagnostics rather than assuming a continuously
realized source at every step.

**Exact per-cell releases.** Set ``super_droplets_moisture.injection.period_s``
to a positive interval in seconds, explicitly set ``particles_per_cell = 1``
and ``multiplicity_type = constant``, and omit ``sd_rate`` and
``min_multiplicity``. This opt-in mode releases one new super-droplet in
*every* cell at step starts ``t_start + n * period_s < t_stop``. It does not
set the resident population. For example, 600 to 3600 seconds with a 10-second
period gives 300 releases; the state at 600 seconds precedes the first one.
Without ``period_s`` the existing injection modes are unchanged.

The source must span the entire stationary, uniform, flat, single-level domain.
Use fixed time steps from a zero-origin simulation clock, with interval and
both time bounds on step boundaries. Each event supplies a full period's dose
even if the final interval ends early. The physical dose per event is ``rate``
times domain volume times period, and must average at least one real particle
per cell. The integer dose is the difference of consecutive rounded cumulative
doses. Floor/ceiling weights share the remainder between cells independently of
MPI layout, rotating the recipients over releases. The schedule and weights
restart from the saved coarse step, with no stochastic timing or replay. Only
after verifying the count and weights of newly created records does ERF log
``Cell pulse: event=... time=... cells=... represented=...``.
