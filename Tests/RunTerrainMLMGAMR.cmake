cmake_minimum_required(VERSION 3.20)

include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
foreach(arg TEST_EXE INPUT WORKING_DIRECTORY LOG NRANKS EXPECTED_LEVELS
            MIN_EVENTS_PER_LEVEL FINAL_STEP DIVERGENCE_TOLERANCE PYTHON_EXECUTABLE)
    if(NOT DEFINED ${arg} OR "${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunTerrainMLMGAMR.cmake: missing required argument ${arg}")
    endif()
endforeach()
if(NOT EXISTS "${TEST_EXE}")
    message(FATAL_ERROR "Terrain MLMG executable is missing: ${TEST_EXE}")
endif()
if(NOT EXISTS "${INPUT}")
    message(FATAL_ERROR "Terrain MLMG input is missing: ${INPUT}")
endif()

erf_mpi_launcher_command(run_command
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS "${NRANKS}"
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "Terrain MLMG AMR regression")
separate_arguments(runtime_options UNIX_COMMAND "${RUNTIME_OPTIONS}")
list(APPEND run_command "${TEST_EXE}" "${INPUT}" ${runtime_options})
execute_process(
    COMMAND ${run_command}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_VARIABLE simulation_stdout
    ERROR_VARIABLE simulation_stderr
    RESULT_VARIABLE simulation_result)
set(combined_output "${simulation_stdout}\n${simulation_stderr}")
file(WRITE "${LOG}"
    "RESULT=${simulation_result}\nNRANKS=${NRANKS}\nINPUT=${INPUT}\n\n${combined_output}")
if(NOT "${simulation_result}" STREQUAL "0")
    message(FATAL_ERROR
        "Terrain MLMG AMR run failed with result ${simulation_result}; see ${LOG}\n${combined_output}")
endif()
if(combined_output MATCHES "Invalid interpolated rho0 momentum ghost|First invalid Omega input|First invalid W input")
    message(FATAL_ERROR "Terrain MLMG AMR run reported an invalid projection ghost; see ${LOG}")
endif()
if(DEFINED REQUIRED_EFFECTIVE_BCS AND NOT "${REQUIRED_EFFECTIVE_BCS}" STREQUAL "")
    string(FIND "${combined_output}" "${REQUIRED_EFFECTIVE_BCS}" effective_bc_position)
    if(effective_bc_position EQUAL -1)
        message(FATAL_ERROR
            "Expected effective terrain MLMG BCs '${REQUIRED_EFFECTIVE_BCS}' were not reported; see ${LOG}")
    endif()
endif()

set(check_command "${PYTHON_EXECUTABLE}"
    "${CMAKE_CURRENT_LIST_DIR}/check_terrain_mlmg_projection.py"
    "--log" "${LOG}"
    "--expected-levels" "${EXPECTED_LEVELS}"
    "--min-events" "${MIN_EVENTS_PER_LEVEL}"
    "--final-step" "${FINAL_STEP}"
    "--divergence-tolerance" "${DIVERGENCE_TOLERANCE}")
if(DEFINED REQUIRED_NEUMANN_FACES AND NOT "${REQUIRED_NEUMANN_FACES}" STREQUAL "")
    list(APPEND check_command "--required-neumann-faces" "${REQUIRED_NEUMANN_FACES}")
endif()
if(DEFINED NEUMANN_FLUX_TOLERANCE AND NOT "${NEUMANN_FLUX_TOLERANCE}" STREQUAL "")
    list(APPEND check_command "--neumann-flux-tolerance" "${NEUMANN_FLUX_TOLERANCE}")
endif()
execute_process(
    COMMAND ${check_command}
    RESULT_VARIABLE checker_result
    OUTPUT_VARIABLE checker_stdout
    ERROR_VARIABLE checker_stderr)
if(NOT "${checker_result}" STREQUAL "0")
    message(FATAL_ERROR
        "Terrain projection checks failed; see ${LOG}\n${checker_stdout}\n${checker_stderr}\n${combined_output}")
endif()
message(STATUS "${checker_stdout}")
