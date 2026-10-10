cmake_minimum_required(VERSION 3.20)
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
foreach(arg TEST_EXE INPUT WORKING_DIRECTORY LOG REQUIRED_MARKER)
    if(NOT DEFINED ${arg} OR "${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunTerrainMLMGRejectIncompatible.cmake: missing required argument ${arg}")
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
    CONTEXT "Terrain MLMG runtime rejection regression")
separate_arguments(runtime_options UNIX_COMMAND "${RUNTIME_OPTIONS}")
list(APPEND run_command "${TEST_EXE}" "${INPUT}" ${runtime_options})
execute_process(
    COMMAND ${run_command}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_VARIABLE simulation_stdout
    ERROR_VARIABLE simulation_stderr
    RESULT_VARIABLE simulation_result)
set(combined_output "${simulation_stdout}\\n${simulation_stderr}")
file(WRITE "${LOG}" "RESULT=${simulation_result}; NRANKS=${NRANKS}; INPUT=${INPUT}; ${combined_output}")
if("${simulation_result}" STREQUAL "0")
    message(FATAL_ERROR "MLMG input unexpectedly succeeded; see ${LOG}: ${combined_output}")
endif()
if(NOT combined_output MATCHES "${REQUIRED_MARKER}")
    message(FATAL_ERROR "MLMG failed without the required diagnostic; see ${LOG}: ${combined_output}")
endif()
if(DEFINED REJECT_BEFORE_PROJECTION AND REJECT_BEFORE_PROJECTION AND
   combined_output MATCHES "ERF_TERRAIN_PROJECTION_EVENT|Terrain MLMG level|Solving in subdomain")
    message(FATAL_ERROR "Unsupported MLMG scope was rejected only after a projection solve; see ${LOG}: ${combined_output}")
endif()
message(STATUS "MLMG configuration rejected with the expected diagnostic.")
