cmake_minimum_required(VERSION 3.18)
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
foreach(arg NRANKS MODE WORKING_DIRECTORY TEST_EXE CHECKER)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunSDWallLoss: missing ${arg}")
    endif()
endforeach()
foreach(width IN ITEMS ${NRANKS} 1)
    erf_mpi_launcher_command(launch_${width} LAUNCHER "${MPIEXEC}"
        NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}" NRANKS ${width} PREFLAGS "${MPIEXEC_PREFLAGS}"
        CONTEXT "RunSDWallLoss")
endforeach()

set(deck CloudChamber_SDM_WallLoss.i)
set(input ${deck})
set(sd super_droplets_moisture)
set(w ${sd}.wall_loss)
set(history ${sd}_wall_loss.txt)
set(on ${w}.enable=true ${w}.num_bins=4 ${w}.bin_rmin=1e-7 ${w}.bin_rmax=1e-3 ${w}.diag_interval=2)
# Synthetic rates chosen for resolvable losses in a six-step acoustic-scale run.
set(constant ${on} ${w}.deposition_velocity_model=constant ${w}.vd_side=128 ${w}.vd_top=192 ${w}.vd_bottom=64)
set(measured ${on} ${w}.deposition_velocity_model=measured_rate
    "${w}.profile_diameters=1 100" ${w}.profile_diameter_units=um "${w}.profile_rates=1024 1024")
set(periodic "geometry.is_periodic=1 1 0" xlo.type=Periodic xhi.type=Periodic ylo.type=Periodic yhi.type=Periodic)

function(fresh name)
    set(dir "${WORKING_DIRECTORY}/${name}")
    file(REMOVE_RECURSE "${dir}")
    file(MAKE_DIRECTORY "${dir}")
    file(COPY "${WORKING_DIRECTORY}/${deck}" "${WORKING_DIRECTORY}/physical_chamber.i" DESTINATION "${dir}")
endfunction()

# A negative leg must exit nonzero AND name its own error: a launcher failure,
# timeout, absent checkpoint, or failure of another solver guard isn't a pass.
function(run name log expected)
    set(dir "${WORKING_DIRECTORY}/${name}")
    execute_process(COMMAND ${launch_${NRANKS}} "${TEST_EXE}" "${input}" amrex.call_addr2line=0 ${ARGN}
        WORKING_DIRECTORY "${dir}" OUTPUT_FILE "${dir}/${log}" ERROR_FILE "${dir}/${log}"
        RESULT_VARIABLE rc TIMEOUT 90)
    if(expected STREQUAL "OK")
        if(NOT rc STREQUAL "0")
            message(FATAL_ERROR "ERF failed (${rc}): ${dir}/${log}")
        endif()
    else()
        file(READ "${dir}/${log}" output)
        # A serial abort can be reported as a signal name, not an exit code.
        if(rc STREQUAL "0" OR rc MATCHES "[Tt]imeout" OR NOT output MATCHES "${expected}")
            message(FATAL_ERROR "Expected '${expected}', got '${rc}': ${dir}/${log}")
        endif()
    endif()
    if(EXISTS "${dir}/wallflux_lev0.txt")
        message(FATAL_ERROR "ERF unexpectedly enabled Eulerian wall fluxes in ${dir}")
    endif()
endfunction()

# Checker reads the *actual* native ERF plotfiles and text history on one rank.
function(check name model first last cadence)
    execute_process(COMMAND ${launch_1} "${CHECKER}" "${WORKING_DIRECTORY}/${name}" ${model}
        ${first} ${last} ${cadence} ${ARGN} WORKING_DIRECTORY "${WORKING_DIRECTORY}/${name}"
        RESULT_VARIABLE rc TIMEOUT 60)
    if(NOT rc STREQUAL "0")
        message(FATAL_ERROR "SD wall-loss check failed (${rc}): ${name} ${model}, steps ${first}..${last}")
    endif()
endfunction()

function(reject_without_touching name reason)
    set(path "${WORKING_DIRECTORY}/${name}/${history}")
    file(SHA256 "${path}" before)
    run(${name} reject.log "${reason}" ${ARGN})
    file(SHA256 "${path}" after)
    if(NOT before STREQUAL after)
        message(FATAL_ERROR "ERF changed ${path} on rejection")
    endif()
endfunction()

if(MODE STREQUAL "Disabled")
    fresh(default_off)
    run(default_off sim.log OK)
    check(default_off off 0 6 0)
    fresh(explicit_off)
    run(explicit_off sim.log OK ${w}.enable=false)
    check(explicit_off off 0 6 0 "${WORKING_DIRECTORY}/default_off")
    # The diagnostic switch must not also switch off the physics.
    fresh(no_text)
    run(no_text sim.log OK ${constant} ${w}.diag_interval=0)
    check(no_text constant 0 6 0)
elseif(MODE STREQUAL "Constant" OR MODE STREQUAL "Measured")
    set(input physical_chamber.i)
    string(TOLOWER "${MODE}" model)
    set(layout "")
    if(NRANKS EQUAL 1)
        set(layout amr.max_grid_size=16) # one box here, split in x in Constant_MPI
    endif()
    fresh(active)
    run(active sim.log OK ${${model}} ${layout})
    check(active ${model} 0 6 2)
elseif(MODE STREQUAL "Periodic")
    foreach(model constant measured)
        fresh(${model})
        run(${model} sim.log OK ${${model}} ${periodic})
        check(${model} ${model}_periodic 0 6 2)
    endforeach()
    # A positive whole-chamber rate must do nothing if there is no wall.
    fresh(no_walls)
    run(no_walls sim.log OK ${measured} ${periodic} "geometry.is_periodic=1 1 1"
        zlo.type=Periodic zhi.type=Periodic max_step=2)
    check(no_walls off 0 2 2) # unchanged mesh AND an enabled, all-zero loss row
elseif(MODE STREQUAL "Settling")
    # ERF permits Open on lateral faces only. The opposite x face is a wall.
    # Advection stores Atlas-Ulbrich vterm even with both motion flags off in
    # the deck. Use rain-sized drops for a well-resolved signal in float plots.
    function(settling_leg name correction top)
        fresh(${name})
        run(${name} sim.log OK ${on} xlo.type=Open xhi.type=SlipWall max_step=3
            ${sd}.include_advection=true ${sd}.terminal_velocity_model=AtlasUlbrich
            ${sd}.initial_species_mean_mass_H2O=2e-7 ${w}.settling_correction=${correction}
            ${w}.vd_side=8 ${w}.vd_top=${top} ${w}.vd_bottom=4 ${w}.diag_interval=1)
        check(${name} ${name} 0 3 1)
    endfunction()
    settling_leg(face_control false 2)
    settling_leg(settling true 2)       # vt > vd_top: no roof loss
    settling_leg(reduced_roof true 8)   # vt < vd_top: reduced, not zero
elseif(MODE STREQUAL "Threshold")
    # 128 drops per SD initially; wall cells cross 103 on different steps.
    # Keep interior SDs alive, and allow recycling of *any* inactive fraction.
    function(threshold_leg name deactivate)
        fresh(${name})
        run(${name} sim.log OK ${constant} ${sd}.recycle_particles=true ${sd}.inactive_threshold=0
            ${w}.deactivate_on_depletion=${deactivate} ${w}.deactivate_threshold=103 ${w}.diag_interval=1
            "erf.plot_vars_1=${sd}_number_density ${sd}_mass_density ${sd}_sd_number_density")
        check(${name} ${name} 0 6 1)
    endfunction()
    threshold_leg(constant false)
    threshold_leg(threshold true)
elseif(MODE STREQUAL "Diffusion")
    fresh(active)
    run(active sim.log OK ${on} ${w}.deposition_velocity_model=diffusion ${w}.diffusion_temperature=293
        ${w}.diffusion_viscosity=1.8e-8 ${w}.diffusion_delta=1e-11 ${w}.diffusion_slip_correction=false)
    check(active diffusion 0 6 2)
elseif(MODE STREQUAL "Stochastic")
    fresh(active)
    run(active sim.log OK ${measured} ${w}.method=stochastic)
    check(active stochastic 0 6 2)
elseif(MODE STREQUAL "Restart")
    set(input physical_chamber.i)
    # Steps 2,4,6 write. Checkpoint 3 contains real particles AND nonzero loss
    # since row 2. A restart may neither replay row 2 nor forget step 3.
    fresh(straight)
    run(straight sim.log OK ${constant} erf.check_int=3)
    check(straight constant 0 6 2)
    fresh(same)
    run(same checkpoint.log OK ${constant} max_step=3 erf.check_int=3)
    check(same constant 0 3 2)
    foreach(name same straight)
        if(NOT EXISTS "${WORKING_DIRECTORY}/${name}/chk00003/super_droplets_moisture/Header")
            message(FATAL_ERROR "No native particle checkpoint in ${name}/chk00003")
        endif()
    endforeach()
    set(resume "erf.restart=${WORKING_DIRECTORY}/same/chk00003")
    file(READ "${WORKING_DIRECTORY}/same/${history}" prefix)

    # Valid, same schema, same last step and values, but a different run's lineage.
    fresh(foreign)
    file(COPY "${WORKING_DIRECTORY}/same/${history}" DESTINATION "${WORKING_DIRECTORY}/foreign")
    reject_without_touching(foreign "foreign history|history content differs" ${constant}
        "erf.restart=${WORKING_DIRECTORY}/straight/chk00003")

    fresh(new)
    run(new restart.log OK ${constant} "${resume}" erf.plot_file_on_restart=true)
    check(new constant 3 6 2 "${WORKING_DIRECTORY}/straight")

    run(same restart.log OK ${constant} "${resume}" erf.plot_file_on_restart=false)
    file(READ "${WORKING_DIRECTORY}/same/${history}" appended)
    string(FIND "${appended}" "${prefix}" prefix_offset)
    if(NOT prefix_offset EQUAL 0)
        message(FATAL_ERROR "Same-directory continuation changed its pre-checkpoint history")
    endif()
    check(same constant 0 6 2 "${WORKING_DIRECTORY}/straight")
    # Genuine ahead log of this same lineage: try to rewind to checkpoint 3.
    reject_without_touching(same "ahead.*history" ${constant} "${resume}")
elseif(MODE STREQUAL "Reject")
    fresh(existing)
    file(WRITE "${WORKING_DIRECTORY}/existing/${history}" "unrelated user file\n")
    reject_without_touching(existing "refuses existing|already exists" ${constant})
    fresh(bad_rates)
    run(bad_rates reject.log "profile_rates.*length" ${measured} ${w}.profile_rates=1024)
    fresh(bad_velocity)
    run(bad_velocity reject.log "vd_top.*nonnegative" ${constant} ${w}.vd_top=-1)
    # Use the generic box, not the named chamber (which independently rejects
    # AMR). No tagged fine patch is needed: the configured hierarchy matters.
    fresh(amr)
    run(amr reject.log "[Ww]all[ _-]loss[^\n]*(level|AMR)" ${constant}
        amr.max_level=1 amr.ref_ratio=2 erf.coupling_type=TwoWay)
    # The same configuration must advance with wall loss disabled. Do not
    # mistake a general ERF input rejection for an enabled-only restriction.
    function(enabled_guard name reason)
        fresh(${name})
        run(${name} reject.log "${reason}" ${w}.enable=true max_step=1 ${ARGN})
        fresh(${name}_off)
        run(${name}_off sim.log OK ${w}.enable=false max_step=1 ${ARGN})
        check(${name}_off off 0 1 0)
    endfunction()
    set(mesh_guard "SD wall loss requires a flat, uniform Cartesian mesh")
    enabled_guard(tight "SD wall loss requires moisture_tight_coupling = false" erf.moisture_tight_coupling=true)
    enabled_guard(terrain "${mesh_guard}" erf.terrain_type=StaticFittedMesh erf.flat_terrain=true)
    enabled_guard(map "${mesh_guard}" erf.test_mapfactor_amp=0.125)
    foreach(name bad_rates bad_velocity amr tight terrain map)
        if(EXISTS "${WORKING_DIRECTORY}/${name}/${history}" OR EXISTS "${WORKING_DIRECTORY}/${name}/plt00001")
            message(FATAL_ERROR "${name} advanced or created a wall-loss history before rejecting inputs")
        endif()
    endforeach()
else()
    message(FATAL_ERROR "Unknown SD wall-loss test mode: ${MODE}")
endif()
