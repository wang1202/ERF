cmake_minimum_required(VERSION 3.20)
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
foreach(arg MPIEXEC MPIEXEC_NUMPROC_FLAG TEST_EXE INPUT SOUNDING WORKING_DIRECTORY
            FCOMPARE PYTHON_EXECUTABLE LOG NRANKS RTOL ATOL)
    if(NOT DEFINED ${arg} OR "${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunTerrainMLMGRankParity.cmake: missing ${arg}")
    endif()
endforeach()
if(NOT EXISTS "${TEST_EXE}" OR NOT EXISTS "${INPUT}" OR NOT EXISTS "${SOUNDING}")
    message(FATAL_ERROR "RunTerrainMLMGRankParity.cmake: executable or input fixture missing")
endif()
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")
erf_resolve_executable(FCOMPARE "${FCOMPARE}" CONFIG "${CONFIG}"
    CONTEXT "RunTerrainMLMGRankParity.cmake: fcompare")
set(ONE_DIR "${WORKING_DIRECTORY}/one_rank")
set(TWO_DIR "${WORKING_DIRECTORY}/two_ranks")
file(REMOVE_RECURSE "${ONE_DIR}" "${TWO_DIR}")
file(MAKE_DIRECTORY "${ONE_DIR}" "${TWO_DIR}")
separate_arguments(runtime_options UNIX_COMMAND "${RUNTIME_OPTIONS}")
erf_mpi_launcher_command(one_command
    LAUNCHER "${MPIEXEC}" NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1 PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "Terrain MLMG one-rank state parity")
erf_mpi_launcher_command(two_command
    LAUNCHER "${MPIEXEC}" NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS} PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "Terrain MLMG two-rank state parity")
execute_process(
    COMMAND ${one_command} "${TEST_EXE}" "${INPUT}" ${runtime_options}
            "erf.input_sounding_file=${SOUNDING}"
    WORKING_DIRECTORY "${ONE_DIR}"
    OUTPUT_FILE "${ONE_DIR}/simulation.log" ERROR_FILE "${ONE_DIR}/simulation.log"
    RESULT_VARIABLE one_result)
if(NOT one_result EQUAL 0)
    message(FATAL_ERROR "Terrain MLMG one-rank simulation failed: ${one_result}; see ${ONE_DIR}/simulation.log")
endif()
execute_process(
    COMMAND ${two_command} "${TEST_EXE}" "${INPUT}" ${runtime_options}
            "erf.input_sounding_file=${SOUNDING}"
    WORKING_DIRECTORY "${TWO_DIR}"
    OUTPUT_FILE "${TWO_DIR}/simulation.log" ERROR_FILE "${TWO_DIR}/simulation.log"
    RESULT_VARIABLE two_result)
if(NOT two_result EQUAL 0)
    message(FATAL_ERROR "Terrain MLMG two-rank simulation failed: ${two_result}; see ${TWO_DIR}/simulation.log")
endif()
foreach(stem "plt" "pltU" "pltV" "pltW")
    set(one_plot "${ONE_DIR}/${stem}00001")
    set(two_plot "${TWO_DIR}/${stem}00001")
    if(NOT EXISTS "${one_plot}" OR NOT EXISTS "${two_plot}")
        message(FATAL_ERROR "Terrain MLMG parity is missing ${stem}00001; check plotfile controls")
    endif()
    execute_process(
        COMMAND ${one_command} "${FCOMPARE}" --abort_if_not_all_found --allow_diff_grids
                --rel_tol "${RTOL}" --abs_tol "${ATOL}" "${one_plot}" "${two_plot}"
        WORKING_DIRECTORY "${WORKING_DIRECTORY}"
        OUTPUT_FILE "${WORKING_DIRECTORY}/${stem}-parity.log"
        ERROR_FILE "${WORKING_DIRECTORY}/${stem}-parity.log"
        RESULT_VARIABLE compare_result)
    if(NOT compare_result EQUAL 0)
        message(FATAL_ERROR "One-rank and two-rank ${stem} plotfiles differ; see ${WORKING_DIRECTORY}/${stem}-parity.log")
    endif()
endforeach()
execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_LIST_DIR}/check_terrain_mlmg_rank_parity.py"
            --one-log "${ONE_DIR}/simulation.log" --two-log "${TWO_DIR}/simulation.log"
            --relative-tolerance "${RTOL}" --absolute-tolerance "${ATOL}"
    RESULT_VARIABLE diagnostic_result
    OUTPUT_VARIABLE diagnostic_stdout
    ERROR_VARIABLE diagnostic_stderr)
if(NOT diagnostic_result EQUAL 0)
    message(FATAL_ERROR "Terrain MLMG rank parity diagnostics failed: ${diagnostic_stdout}${diagnostic_stderr}")
endif()
message(STATUS "${diagnostic_stdout}")
