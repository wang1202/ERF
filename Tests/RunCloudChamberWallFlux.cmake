include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

# Runs for the interval-mean Cloud Chamber wall-flux output (wallflux_lev0.txt,
# erf.wallflux_diagnostics); Tests/CloudChamberWallFluxCheck.cpp inspects them.
#
#   MODE=dry     20 fixed steps of 0.025 s, interval 10, budget off: two rows,
#                one box on one rank against many boxes on NRANKS ranks, and
#                the fields with the diagnostic disabled.
#   MODE=satadj  SatAdj with distinct wet and dry walls, wall-flux and budget
#                at the same cadence, 25 steps (the short trailing window
#                writes no row), under both anelastic integrators.
#   MODE=sdm     SuperDroplets that never injects (the same walls): the wall
#                vapor flux must close the vapor mass change under both
#                integrators, no cloud may form, and the shared fields must
#                match MoistNoCondensation run without the diagnostic.
#   MODE=reject  an existing output file, a nonpositive interval and a restart
#                must each stop the run before it starts.

if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED NRANKS OR NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED MODE)
    message(FATAL_ERROR "RunCloudChamberWallFlux.cmake missing required argument")
endif()
if(NOT MODE STREQUAL "reject" AND NOT DEFINED CHECKER)
    message(FATAL_ERROR "RunCloudChamberWallFlux.cmake: MODE=${MODE} needs CHECKER")
endif()

erf_mpi_launcher_command(mpi_launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunCloudChamberWallFlux.cmake")
erf_mpi_launcher_command(mpi_launch_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunCloudChamberWallFlux.cmake")

# run_erf(<name> <ONE|MANY> <expect_success> <args...>): run ERF in a fresh
# WORKING_DIRECTORY/<name>, logging to simulation.log there.
function(run_erf name ranks expect_success)
    set(dir "${WORKING_DIRECTORY}/${name}")
    if(ranks STREQUAL "ONE")
        set(launch ${mpi_launch_one})
    else()
        set(launch ${mpi_launch})
    endif()
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${INPUT} ${ARGN}
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${dir}/simulation.log"
        ERROR_FILE "${dir}/simulation.log"
        RESULT_VARIABLE result)
    if(expect_success AND NOT result EQUAL 0)
        message(FATAL_ERROR "wall-flux run ${name} failed: ${result}")
    elseif(NOT expect_success AND result EQUAL 0)
        message(FATAL_ERROR "wall-flux run ${name} was expected to stop but succeeded")
    endif()
endfunction()

function(fresh_dir name)
    file(REMOVE_RECURSE "${WORKING_DIRECTORY}/${name}")
    file(MAKE_DIRECTORY "${WORKING_DIRECTORY}/${name}")
endfunction()

function(run_checker log)
    execute_process(
        COMMAND ${mpi_launch_one} ${CHECKER} ${ARGN}
        WORKING_DIRECTORY "${WORKING_DIRECTORY}"
        OUTPUT_FILE "${WORKING_DIRECTORY}/${log}"
        ERROR_FILE "${WORKING_DIRECTORY}/${log}"
        RESULT_VARIABLE result)
    if(NOT result EQUAL 0)
        file(READ "${WORKING_DIRECTORY}/${log}" text)
        message(FATAL_ERROR "wall-flux check failed (${log}):\n${text}")
    endif()
endfunction()

function(expect_log name pattern)
    file(READ "${WORKING_DIRECTORY}/${name}/simulation.log" text)
    if(NOT text MATCHES "${pattern}")
        message(FATAL_ERROR "wall-flux run ${name} did not report '${pattern}':\n${text}")
    endif()
endfunction()

if(MODE STREQUAL "dry")
    foreach(name multi one_box disabled)
        fresh_dir(${name})
    endforeach()
    run_erf(multi MANY TRUE)
    run_erf(one_box ONE TRUE amr.max_grid_size=64)
    run_erf(disabled MANY TRUE erf.wallflux_diagnostics=false)
    if(EXISTS "${WORKING_DIRECTORY}/disabled/wallflux_lev0.txt")
        message(FATAL_ERROR "the disabled run wrote wallflux_lev0.txt")
    endif()
    run_checker(checker_dry.log dry
        "${WORKING_DIRECTORY}/multi" "${WORKING_DIRECTORY}/one_box"
        "${WORKING_DIRECTORY}/disabled" plt00020 2 0.25)
elseif(MODE STREQUAL "satadj")
    foreach(scheme MidPoint RK2)
        fresh_dir(${scheme})
        run_erf(${scheme} MANY TRUE erf.timestepping_type=${scheme})
        run_checker(checker_${scheme}.log satadj
            "${WORKING_DIRECTORY}/${scheme}" plt00000 2 0.1 "xhi,yhi")
    endforeach()
elseif(MODE STREQUAL "sdm")
    foreach(scheme MidPoint RK2)
        fresh_dir(${scheme})
        run_erf(${scheme} MANY TRUE erf.timestepping_type=${scheme})
        run_checker(checker_${scheme}.log vapor_closure
            "${WORKING_DIRECTORY}/${scheme}" 2 0.1 10 "xhi,yhi")
    endforeach()
    # Without droplets SuperDroplets only transports vapor, like
    # MoistNoCondensation (run here without the diagnostic).
    fresh_dir(noncondensing)
    run_erf(noncondensing MANY TRUE erf.timestepping_type=RK2
        erf.moisture_model=MoistNoCondensation erf.wallflux_diagnostics=false)
    run_checker(checker_same_fields.log same_fields
        "${WORKING_DIRECTORY}/RK2/plt00020" "${WORKING_DIRECTORY}/noncondensing/plt00020" 1.0e-9)
elseif(MODE STREQUAL "reject")
    fresh_dir(existing)
    set(unrelated "unrelated user data\n")
    file(WRITE "${WORKING_DIRECTORY}/existing/wallflux_lev0.txt" "${unrelated}")
    run_erf(existing MANY FALSE)
    expect_log(existing "wallflux_lev0.txt already exists")
    file(READ "${WORKING_DIRECTORY}/existing/wallflux_lev0.txt" kept)
    if(NOT kept STREQUAL unrelated)
        message(FATAL_ERROR "the run modified an existing wallflux_lev0.txt")
    endif()

    fresh_dir(interval)
    run_erf(interval MANY FALSE erf.wallflux_output_interval=0)
    expect_log(interval "wallflux_output_interval must be a positive number of steps")

    fresh_dir(restart)
    run_erf(restart MANY FALSE amr.restart=chk_wallflux_missing)
    expect_log(restart "restart is unsupported")
else()
    message(FATAL_ERROR "RunCloudChamberWallFlux.cmake: unknown MODE ${MODE}")
endif()
