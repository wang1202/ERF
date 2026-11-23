# A quiescent, flat box. Kilometre-sized cells allow a cheap fixed 2.5-s
# acoustic step without changing the ten-second source period. Start at t=0;
# the first advanced states containing new particles are 7.5, 17.5 and 27.5 s.
erf.prob_name = Bubble
erf.init_type = Uniform
prob.T_pert = 0
prob.U_0 = 0
prob.do_moist_bubble = false
erf.use_gravity = false
erf.les_type = None
erf.moisture_model = SuperDroplets
erf.timestepping_type = RK3_NoSubstepping
erf.fixed_dt = 2.5
max_step = 15

geometry.prob_lo = 0 0 0
geometry.prob_hi = 20000 10000 10000
geometry.is_periodic = 1 1 0
zlo.type = SlipWall
zhi.type = SlipWall
amr.n_cell = 8 4 4
amr.max_level = 0
amr.blocking_factor = 4
amr.max_grid_size = 8

erf.sum_interval = -1
erf.check_file = chk
erf.check_int = 3
erf.file_name_digits = 5
erf.plot_file_1 = plt
erf.plot_int_1 = 1
erf.plot_file_on_restart = true
erf.plot_vars_1 = density
particles.disable_plt = false
erf.fix_random_seed = 1

super_droplets_moisture.num_initializations = 0
super_droplets_moisture.num_injections = 1
super_droplets_moisture.include_advection = false
super_droplets_moisture.advect_with_flow = false
super_droplets_moisture.advect_with_gravity = false
super_droplets_moisture.include_coalescence = false
super_droplets_moisture.include_phase_change = false
super_droplets_moisture.include_cold_processes = false
super_droplets_moisture.wall_loss.enable = false
super_droplets_moisture.recycle_particles = false
super_droplets_moisture.density_scaling = false
super_droplets_moisture.place_randomly_in_cells = false
super_droplets_moisture.diagnostics_interval = 0

super_droplets_moisture.injection.distribution_type = uniform
super_droplets_moisture.injection.domain_velocity = 0 0 0
super_droplets_moisture.injection.multiplicity_type = constant
super_droplets_moisture.injection.species_mean_mass_H2O = 1e-12
super_droplets_moisture.injection.particles_per_cell = 1
super_droplets_moisture.injection.period_s = 10
super_droplets_moisture.injection.t_start = 5
super_droplets_moisture.injection.t_stop = 35
# 299.375 physical droplets/domain/pulse over 128 cells: each newly born SD
# has weight 2 or 3. Rounding each event in isolation fails on the second one.
super_droplets_moisture.injection.rate = 1.496875e-11
