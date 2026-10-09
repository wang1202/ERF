cmake_minimum_required(VERSION 3.20)

include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

foreach(arg TEST_EXE INPUT WORKING_DIRECTORY LOG)
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
    NRANKS 1
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
    "RESULT=${simulation_result}\nINPUT=${INPUT}\n\n${combined_output}")

if(NOT "${simulation_result}" STREQUAL "0")
    message(FATAL_ERROR
        "Terrain MLMG AMR run failed with result ${simulation_result}; see ${LOG}")
endif()
if(combined_output MATCHES "Invalid interpolated rho0 momentum ghost|First invalid Omega input")
    message(FATAL_ERROR "Terrain MLMG AMR run reported an invalid projection ghost; see ${LOG}")
endif()
if(NOT combined_output MATCHES "Coarse STEP 1 ends")
    message(FATAL_ERROR "Terrain MLMG AMR run did not complete one full timestep; see ${LOG}")
endif()

string(REPLACE "\n" ";" output_lines "${combined_output}")
set(neumann_flux_lines)
foreach(output_line IN LISTS output_lines)
    if(output_line MATCHES "^Terrain MLMG Neumann normal correction flux max:")
        list(APPEND neumann_flux_lines "${output_line}")
    endif()
endforeach()
if(NOT neumann_flux_lines)
    message(FATAL_ERROR "No Neumann boundary-normal flux measurements were reported; see ${LOG}")
endif()
foreach(flux_line IN LISTS neumann_flux_lines)
    if(NOT flux_line MATCHES "max_abs=([0-9.eE+-]+)")
        message(FATAL_ERROR "Cannot parse a Neumann flux measurement: ${flux_line}")
    endif()
    set(max_normal_flux "${CMAKE_MATCH_1}")
    if(NOT flux_line MATCHES "samples=([0-9]+)")
        message(FATAL_ERROR "Cannot parse Neumann face sample count: ${flux_line}")
    endif()
    set(sample_count "${CMAKE_MATCH_1}")
    if(NOT sample_count GREATER 0)
        message(FATAL_ERROR "Neumann face had no valid flux samples: ${flux_line}")
    endif()
    if("${max_normal_flux}" GREATER 1.0e-12)
        message(FATAL_ERROR
            "Neumann boundary-normal correction flux ${max_normal_flux} exceeds 1e-12: ${flux_line}")
    endif()
endforeach()

foreach(level 0 1)
    if(NOT combined_output MATCHES "Terrain MLMG level ${level} subdomain 0:")
        message(FATAL_ERROR "ERF level ${level} did not reach terrain MLMG; see ${LOG}")
    endif()
    string(REGEX MATCH
        "Max/L2 norm of divergence after[ \t]+solve at level ${level} : ([0-9.eE+-]+)"
        divergence_line "${combined_output}")
    if(NOT divergence_line)
        message(FATAL_ERROR "No post-projection divergence was reported for level ${level}; see ${LOG}")
    endif()
    set(divergence_linf "${CMAKE_MATCH_1}")
    if("${divergence_linf}" GREATER 1.0e-6)
        message(FATAL_ERROR
            "Level ${level} post-projection divergence ${divergence_linf} exceeds 1e-6; see ${LOG}")
    endif()
    message(STATUS "Terrain MLMG level ${level}: post-projection L_inf=${divergence_linf}")
endforeach()
