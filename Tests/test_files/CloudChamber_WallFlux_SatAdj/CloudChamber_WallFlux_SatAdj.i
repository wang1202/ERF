# SatAdj Cloud Chamber regression for the interval-mean wall-flux output
# (wallflux_lev0.txt) against the native budget at the same cadence.  Every
# wall differs (wet or dry, temperature) so a face swap cannot go unnoticed;
# xhi and yhi are dry and must exchange no vapor.  25 steps with interval 10:
# the trailing five-step window writes no row.
erf.prob_name = "Cloud Chamber"
erf.init_type = ConstantDensity
erf.compressibility_type = Anelastic
erf.use_gravity = true
erf.moisture_model = SatAdj
erf.vert_implicit = false
erf.project_initial_velocity = 0
erf.molec_diff_type = ConstantAlpha
erf.dynamic_viscosity = 0.0
erf.alpha_T = 0.001
erf.alpha_C = 0.001
erf.les_type = Smagorinsky
erf.Cs = 0.1
erf.fixed_dt = 0.01
erf.sum_interval = -1
erf.check_int = -1
erf.cloud_chamber_budget_interval = 10
erf.wallflux_diagnostics = true
erf.wallflux_output_interval = 10
erf.plot_file_1 = plt
erf.plot_int_1 = 10
erf.plot_vars_1 = density theta temp pressure pres_hse qv qc x_velocity y_velocity z_velocity

max_step = 25
geometry.prob_lo = 0.0 0.0 0.0
geometry.prob_hi = 2.0 2.0 1.0
geometry.is_periodic = 0 0 0
amr.n_cell = 16 16 16
amr.max_level = 0
amr.max_grid_size = 8

prob.p_inf = 100000.0
prob.U_0 = 0.5
prob.T_0 = 292.0
prob.thermodynamic_initialization = physical_temperature_rh
prob.initial_temperature_bottom = 300.0
prob.initial_temperature_top = 284.0
prob.initial_relative_humidity = 0.95
prob.temperature_perturbation_amplitude = 0.02
prob.perturbation_mode = deterministic_sine

xlo.type = NoSlipWall
xhi.type = NoSlipWall
ylo.type = NoSlipWall
yhi.type = NoSlipWall
zlo.type = NoSlipWall
zhi.type = NoSlipWall
xlo.temperature = 296.0
xlo.moisture = wet
xlo.wall_transfer_model = resolved_molecular
xhi.temperature = 292.0
xhi.moisture = dry
xhi.wall_transfer_model = resolved_molecular
ylo.temperature = 288.0
ylo.moisture = wet
ylo.wall_transfer_model = resolved_molecular
yhi.temperature = 290.0
yhi.moisture = dry
yhi.wall_transfer_model = resolved_molecular
zlo.temperature = 300.0
zlo.moisture = wet
zlo.wall_transfer_model = resolved_molecular
zhi.temperature = 284.0
zhi.moisture = wet
zhi.wall_transfer_model = resolved_molecular
