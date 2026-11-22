# Quiescent 2 x 2 x 1 m six-wall chamber. Use a uniform base state and a
# zero-amplitude Bubble so even fluid gravity can be off and the periodic
# control can use the same deck. The named Cloud Chamber requires fluid
# gravity and six nonperiodic faces; the driver also exercises it separately.
erf.prob_name = Bubble
erf.init_type = Uniform
prob.T_0 = 292
prob.U_0 = 0
prob.T_pert = 0
prob.do_moist_bubble = false
erf.use_gravity = false
erf.fix_random_seed = 1
erf.les_type = None
erf.moisture_model = SuperDroplets
erf.timestepping_type = RK3_NoSubstepping
erf.fixed_dt = 0.0001220703125   # 1/8192 s
max_step = 6

geometry.prob_lo = 0 0 0
geometry.prob_hi = 2 2 1
geometry.is_periodic = 0 0 0
amr.n_cell = 16 8 8
amr.max_grid_size = 8
amr.max_grid_size_z = 8
amr.max_level = 0
xlo.type = NoSlipWall
xhi.type = NoSlipWall
ylo.type = NoSlipWall
yhi.type = NoSlipWall
zlo.type = NoSlipWall
zhi.type = NoSlipWall

erf.wallflux_diagnostics = false
erf.sum_interval = -1
erf.check_file = chk
erf.check_int = -1
erf.plot_file_1 = plt
erf.plot_int_1 = 1
erf.plot_vars_1 = super_droplets_moisture_number_density super_droplets_moisture_mass_density
particles.disable_plt = true   # Only omit the particle dump, not its mesh moments or checkpoint.

# 1 SD at the exact center of every cell, representing 128 physical drops.
# Both H2O and NaCl contribute appreciably to the total mass of 2e-12 kg/drop.
super_droplets_moisture.num_initializations = 1
super_droplets_moisture.num_injections = 0
super_droplets_moisture.initial_distribution_type = uniform
super_droplets_moisture.initial_particles_per_cell = 1
super_droplets_moisture.initial_number_density = 32768
super_droplets_moisture.multiplicity_type = constant
super_droplets_moisture.place_randomly_in_cells = false
super_droplets_moisture.aerosols = NaCl
super_droplets_moisture.initial_aerosol_distribution_type_NaCl = mass_constant
super_droplets_moisture.initial_aerosol_mean_mass_NaCl = 0.5e-12
super_droplets_moisture.initial_species_distribution_type_H2O = mass_constant
super_droplets_moisture.initial_species_mean_mass_H2O = 1.5e-12
super_droplets_moisture.density_scaling = false
super_droplets_moisture.initial_phase_change_relaxation = false
super_droplets_moisture.include_phase_change = false
super_droplets_moisture.include_cold_processes = false
super_droplets_moisture.include_advection = false
super_droplets_moisture.advect_with_flow = false
super_droplets_moisture.advect_with_gravity = false
super_droplets_moisture.include_coalescence = false
super_droplets_moisture.recycle_particles = false
super_droplets_moisture.diagnostics_interval = 0

# Intentionally no wall_loss keys: this deck alone tests the default-off path.
# The driver supplies synthetic velocities/rates, not a chamber calibration.
