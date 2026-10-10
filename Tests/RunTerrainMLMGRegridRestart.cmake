cmake_minimum_required(VERSION 3.20)
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

foreach(arg TEST_EXE INPUT WORKING_DIRECTORY NRANKS PYTHON_EXECUTABLE DIVERGENCE_TOLERANCE)
    if(NOT DEFINED ${arg} OR "${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunTerrainMLMGRegridRestart.cmake: missing required argument ${arg}")
    endif()
endforeach()

set(run_root "${WORKING_DIRECTORY}/terrain_mlmg_regrid_restart")
file(REMOVE_RECURSE "${run_root}")
file(MAKE_DIRECTORY "${run_root}/checkpoint" "${run_root}/restart")
file(GLOB auxiliary_files "${WORKING_DIRECTORY}/*")
foreach(auxiliary_file IN LISTS auxiliary_files)
    if(NOT IS_DIRECTORY "${auxiliary_file}")
        file(COPY "${auxiliary_file}" DESTINATION "${run_root}/checkpoint")
        file(COPY "${auxiliary_file}" DESTINATION "${run_root}/restart")
    endif()
endforeach()

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS "${NRANKS}"
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "Terrain MLMG regrid/restart regression")

set(common_options
    erf.anelastic=1
    erf.terrain_poisson_solver=mlmg
    erf.mg_v=2)
set(checkpoint_options
    max_step=1
    erf.check_int=1
    erf.plot_int_1=-1
    amr.max_grid_size_x=8
    amr.max_grid_size_y=8
    amr.max_grid_size_z=8)
set(restart_options
    erf.restart=chk00001
    max_step=3
    erf.check_int=-1
    erf.plot_int_1=3
    amr.max_grid_size_x=16
    amr.max_grid_size_y=16
    amr.max_grid_size_z=8
    erf.regrid_int=1)

execute_process(
    COMMAND ${launch} "${TEST_EXE}" "${INPUT}" ${common_options} ${checkpoint_options}
    WORKING_DIRECTORY "${run_root}/checkpoint"
    OUTPUT_FILE "${run_root}/checkpoint/simulation.log"
    ERROR_FILE "${run_root}/checkpoint/simulation.log"
    TIMEOUT 300
    RESULT_VARIABLE checkpoint_result)
if(NOT "${checkpoint_result}" STREQUAL "0")
    message(FATAL_ERROR "Terrain MLMG checkpoint leg failed (${checkpoint_result}); see ${run_root}/checkpoint/simulation.log")
endif()
if(NOT EXISTS "${run_root}/checkpoint/chk00001/Header")
    message(FATAL_ERROR "Terrain MLMG checkpoint leg did not write chk00001")
endif()

execute_process(
    COMMAND ${launch} "${TEST_EXE}" "${INPUT}" ${common_options} ${restart_options}
    WORKING_DIRECTORY "${run_root}/restart"
    OUTPUT_FILE "${run_root}/restart/simulation.log"
    ERROR_FILE "${run_root}/restart/simulation.log"
    TIMEOUT 300
    RESULT_VARIABLE restart_result)
if(NOT "${restart_result}" STREQUAL "0")
    message(FATAL_ERROR "Terrain MLMG restart/regrid leg failed (${restart_result}); see ${run_root}/restart/simulation.log")
endif()
if(NOT EXISTS "${run_root}/restart/plt00003/Level_1/Cell_H")
    message(FATAL_ERROR "Terrain MLMG restart/regrid leg did not write the expected two-level plt00003")
endif()

file(READ "${run_root}/restart/simulation.log" restart_log)
if(NOT restart_log MATCHES "REMAKING WITH NEW BA AT LEVEL 1")
    message(FATAL_ERROR "Terrain MLMG restart did not recreate level-1 grids; see ${run_root}/restart/simulation.log")
endif()
execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_LIST_DIR}/check_terrain_mlmg_regrid.py"
        --checkpoint-cell-h "${run_root}/checkpoint/chk00001/Level_1/Cell_H"
        --final-cell-h "${run_root}/restart/plt00003/Level_1/Cell_H"
        --restart-log "${run_root}/restart/simulation.log"
    RESULT_VARIABLE box_check_result
    OUTPUT_VARIABLE box_check_stdout
    ERROR_VARIABLE box_check_stderr)
if(NOT "${box_check_result}" STREQUAL "0")
    message(FATAL_ERROR "Terrain MLMG level-1 BoxArray check failed: ${box_check_stdout}${box_check_stderr}")
endif()

set(check_command "${PYTHON_EXECUTABLE}"
    "${CMAKE_CURRENT_LIST_DIR}/check_terrain_mlmg_projection.py"
    --log "${run_root}/restart/simulation.log"
    --expected-levels 0,1
    --min-events 1
    --final-step 3
    --divergence-tolerance "${DIVERGENCE_TOLERANCE}"
    --required-neumann-faces "0:z:lo,0:z:hi,1:x:lo,1:x:hi,1:y:lo,1:y:hi,1:z:lo,1:z:hi"
    --neumann-flux-tolerance 1.0e-12)
execute_process(
    COMMAND ${check_command}
    RESULT_VARIABLE projection_result
    OUTPUT_VARIABLE projection_stdout
    ERROR_VARIABLE projection_stderr)
if(NOT "${projection_result}" STREQUAL "0")
    message(FATAL_ERROR "Terrain MLMG post-regrid projection validation failed: ${projection_stdout}${projection_stderr}")
endif()
message(STATUS "${box_check_stdout}${projection_stdout}")
