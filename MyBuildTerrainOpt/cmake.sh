#!/bin/bash

# --- Log file setup ---
LOG_FILE="cmake_build_$(date '+%Y%m%d_%H%M%S').log"
exec > >(tee -a "${LOG_FILE}") 2>&1
echo "Logging output to: ${LOG_FILE}"
echo "Build started at: $(date)"
echo "==========================================="

# --- Load a CMake new enough to know about nvcc 12.x C++20 (cuda_std_20) ---
# The system default (cmake 3.20.2) errors out with:
#   "The compiler feature 'cuda_std_20' is not known to CUDA compiler 'NVIDIA' version 12.9.41"
module load cmake/3.29.2

# # --- Ensure submodules required by the build are initialized ---
# # googletest is needed because ERF_ENABLE_TESTS=ON below;
# # Tests/Unit/CMakeLists.txt does add_subdirectory(Submodules/googletest)
# # and configuration fails if that directory is empty.
# ERF_SRC_DIR="$(cd .. && pwd)"
# for sm in googletest AMReX RRTMGP Noah-MP ekat WW3; do
#     if [ ! -f "${ERF_SRC_DIR}/Submodules/${sm}/CMakeLists.txt" ] \
#        && [ ! -f "${ERF_SRC_DIR}/Submodules/${sm}/Makefile" ]; then
#         echo "Initializing submodule: Submodules/${sm}"
#         (cd "${ERF_SRC_DIR}" && git submodule update --init --recursive "Submodules/${sm}")
#     fi
# done

# --- Clean the build directory first (IMPORTANT!) ---
rm -rf CMakeCache.txt CMakeFiles/ install

# --- Define paths for clarity ---
# HDF5_INSTALL_DIR="/home/wang1202/oracle/aaronwang/hdf5_install"
HDF5_INSTALL_DIR="/nopt/nrel/apps/cpu_stack/libraries-craympich/hdf5/hdf5-1.14.6"
NETCDF_INSTALL_DIR="/nopt/nrel/apps/cpu_stack/libraries-craympich/netcdf/netcdf-4.9.3"
CRAY_MPI_DIR="/opt/cray/pe/mpich/8.1.28/ofi/gnu/10.3"
CRAY_MPI_BIN_DIR="${CRAY_MPI_DIR}/bin"
CRAY_MPI_INC_DIR="${CRAY_MPI_DIR}/include"
CRAY_MPI_LIB_DIR="${CRAY_MPI_DIR}/lib"

# --- Explicitly define NetCDF path for global flags ---
NETCDF_DIR="${NETCDF_INSTALL_DIR}"
NETCDF_INCLUDE_PATH="${NETCDF_INSTALL_DIR}/include"

# --- Submodule googletest include (must beat conda's old gtest) ---
# The loaded `anaconda3` module places /nopt/.../conda_23/include on the
# system include path. CMake adds the submodule googletest with -isystem
# AFTER conda's -isystem path, so #include <gtest/gtest.h> resolves to
# conda's old gtest 1.x while libgtest.a is built from the fresh submodule.
# Result: undefined reference to testing::internal::MakeAndRegisterTestInfo
# (signature mismatch — old gtest takes const char*, new gtest takes std::string).
# We prepend the submodule path with -I so it wins over all -isystem entries.
ERF_SRC_DIR="$(cd .. && pwd)"
GTEST_INCLUDE_DIR="${ERF_SRC_DIR}/Submodules/googletest/googletest/include"

# --- Make sure runtime finds correct HDF5 first ---
export LD_LIBRARY_PATH="${HDF5_INSTALL_DIR}/lib:${LD_LIBRARY_PATH}"
export CMAKE_LIBRARY_PATH="${HDF5_INSTALL_DIR}/lib:${CMAKE_LIBRARY_PATH}"

# --- GPU ARCHITECTURE - Set for NVIDIA H100 (with correct case) ---
KOKKOS_GPU_ARCH="HOPPER90" # <-- THIS IS THE FIX (was "Hopper90")
CMAKE_CUDA_ARCH="90"

# --- CUDA Toolkit ---
CUDA_ROOT="/nopt/cuda/12.9"
export PATH="${CUDA_ROOT}/bin:${PATH}"
export LD_LIBRARY_PATH="${CUDA_ROOT}/lib64:${LD_LIBRARY_PATH}"

# Standard MPI without GPU-aware GTL transport
export MPICH_GPU_SUPPORT_ENABLED=0

# --- Let test discovery start the CUDA unit-test binary on this CPU node ---
# Tests/Unit/CMakeLists.txt registers GoogleTest cases with
# `gtest_discover_tests(... DISCOVERY_MODE POST_BUILD)`, which EXECUTES
# erf_unit_tests during `make` to enumerate test names. This is a CUDA build,
# so the binary's NEEDED entry is libcuda.so.1; the dynamic loader must resolve
# it before main() runs. This CPU build node has no GPU driver (no libcuda.so.1),
# only the CUDA toolkit stub libcuda.so (SONAME libcuda.so.1), so the process
# can't even start -> "error while loading shared libraries: libcuda.so.1".
# ERF_GTestMain.cpp already skips amrex::Initialize under --gtest_list_tests, so
# discovery never touches the GPU; it just needs *a* libcuda.so.1 to load. We
# point a symlink at the stub and put it on LD_LIBRARY_PATH for the build only.
# (Actual GPU test runs happen later via `ctest` on a GPU node in a fresh shell,
#  so this stub path does not shadow the real driver there.)
CUDA_STUB_LINK_DIR="$(pwd)/.cuda_stub"
mkdir -p "${CUDA_STUB_LINK_DIR}"
ln -sf "${CUDA_ROOT}/lib64/stubs/libcuda.so" "${CUDA_STUB_LINK_DIR}/libcuda.so.1"
export LD_LIBRARY_PATH="${CUDA_STUB_LINK_DIR}:${LD_LIBRARY_PATH}"

# --- Store all CMake arguments in a bash array ---
STD_LIBS="-L/opt/cray/pe/mpich/8.1.28/ofi/gnu/10.3/lib -lmpi_gnu_103 -L/nopt/cuda/12.9/lib64/stubs -lcuda"

CMAKE_ARGS=(
    "-DCMAKE_INSTALL_PREFIX:PATH=./install"
    "-DCMAKE_CUDA_STANDARD_LIBRARIES=${STD_LIBS}"
    "-DCMAKE_CXX_STANDARD_LIBRARIES=${STD_LIBS}"
    
    # --- Compiler and MPI Settings ---
    "-DCMAKE_C_COMPILER=${CRAY_MPI_BIN_DIR}/mpicc"
    "-DCMAKE_CXX_COMPILER=${CRAY_MPI_BIN_DIR}/mpicxx"
    "-DCMAKE_Fortran_COMPILER=${CRAY_MPI_BIN_DIR}/mpifort"
    "-DCMAKE_CUDA_HOST_COMPILER=${CRAY_MPI_BIN_DIR}/mpicxx"
    "-DMPIEXEC_PREFLAGS:STRING=--oversubscribe"

    # --- ERF Build Options ---
    "-DCMAKE_BUILD_TYPE:STRING=Release"
    "-DERF_DIM:STRING=3"
    "-DERF_ENABLE_MPI:BOOL=ON"
    "-DERF_ENABLE_CUDA:BOOL=ON"
    "-DERF_ENABLE_HDF5:BOOL=ON"
    "-DERF_ENABLE_NETCDF:BOOL=ON"
    "-DERF_ENABLE_FFT:BOOL=ON"
    "-DERF_ENABLE_KOKKOS:BOOL=ON"
    "-DERF_ENABLE_RRTMGP:BOOL=ON"
    "-DERF_ENABLE_NOAHMP:BOOL=ON"
    "-DERF_ENABLE_MORR_FORT:BOOL=ON"   # Exp 19a: build Fortran Morrison alongside C++ to test divergence (erf.use_morr_cpp_answer=false)
    "-DERF_ENABLE_P3=OFF"
    "-DERF_ENABLE_SHOC=OFF"
    "-DERF_ENABLE_EKAT:BOOL=ON"
    "-DERF_ENABLE_TESTS:BOOL=OFF"
    "-DERF_ENABLE_FCOMPARE:BOOL=ON"
    "-DERF_ENABLE_DOCUMENTATION:BOOL=OFF"
    "-DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON"
    "-DERF_ENABLE_PARTICLES=ON"
    
    # --- Fix EKAT MPI detection for Cray MPI (which is MPICH-based) ---
    "-DHAVE_MPICH=TRUE"
    "-DERF_CHECK_MODULES=OFF"
    "-DMPICH_SKIP_MPICXX=TRUE"

    # --- Pre-populate FindMPI cache vars to skip wrapper interrogation. ---
    # CMake's FindMPI probes `mpicc -showme:compile`, which the Cray MPICH wrapper
    # does not implement and falls through into an infinite re-exec loop.
    # Setting MPI_<lang>_LIB_NAMES + MPI_<lang>_INCLUDE_DIRS satisfies the guard
    # in FindMPI.cmake (~line 1481) so the whole interrogate block is skipped.
    "-DMPI_SKIP_COMPILER_WRAPPER=TRUE"
    "-DMPI_C_LIB_NAMES=mpi_gnu_103"
    "-DMPI_C_INCLUDE_DIRS=${CRAY_MPI_INC_DIR}"
    "-DMPI_CXX_LIB_NAMES=mpi_gnu_103"
    "-DMPI_CXX_INCLUDE_DIRS=${CRAY_MPI_INC_DIR}"
    "-DMPI_Fortran_LIB_NAMES=mpifort_gnu_103;mpi_gnu_103"
    "-DMPI_Fortran_INCLUDE_DIRS=${CRAY_MPI_INC_DIR}"
    "-DMPI_mpi_gnu_103_LIBRARY=${CRAY_MPI_LIB_DIR}/libmpi_gnu_103.so"
    "-DMPI_mpifort_gnu_103_LIBRARY=${CRAY_MPI_LIB_DIR}/libmpifort_gnu_103.so"

    # --- GPU Architecture Settings ---
    "-DCMAKE_CUDA_ARCHITECTURES=${CMAKE_CUDA_ARCH}"
    "-DCUDAToolkit_ROOT=${CUDA_ROOT}"
    "-DKokkos_ARCH_${KOKKOS_GPU_ARCH}=ON"
    
    # --- Forcing include path for stubborn submodules ---
    # The conda env has a NON-parallel HDF5 in its include path which gets picked
    # up before the Cray-MPICH parallel HDF5.  Using -I (not -isystem) for the
    # parallel HDF5 ensures it wins over the conda -isystem path.
    "-DCMAKE_CXX_FLAGS=-I${GTEST_INCLUDE_DIR} -I${HDF5_INSTALL_DIR}/include -I${NETCDF_INCLUDE_PATH}"
    "-DCMAKE_CUDA_FLAGS=-I${GTEST_INCLUDE_DIR} -I${HDF5_INSTALL_DIR}/include -I${CRAY_MPI_INC_DIR}"

    # HDF5 explicit paths
    "-DHDF5_ROOT=${HDF5_INSTALL_DIR}"
    "-DHDF5_DIR=${HDF5_INSTALL_DIR}"
    "-DHDF5_INCLUDE_DIR=${HDF5_INSTALL_DIR}/include"
    "-DHDF5_LIBRARY_DIR=${HDF5_INSTALL_DIR}/lib"
    "-DHDF5_HL_LIBRARY=${HDF5_INSTALL_DIR}/lib/libhdf5_hl.so"
    "-DHDF5_C_LIBRARY=${HDF5_INSTALL_DIR}/lib/libhdf5.so"
    "-DCMAKE_INSTALL_RPATH=${HDF5_INSTALL_DIR}/lib"
    "-DCMAKE_BUILD_WITH_INSTALL_RPATH=ON"
    
    # NetCDF path
    "-DNETCDF_DIR=${NETCDF_DIR}"
    "-DCMAKE_PREFIX_PATH=${HDF5_INSTALL_DIR};${NETCDF_INSTALL_DIR}"

    # --- Source Directory ---
    ".."
)

# --- Execute cmake with the arguments, then build ---
echo "Running CMake with the following arguments:"
printf "  %s\n" "${CMAKE_ARGS[@]}"
echo "-------------------------------------------"

cmake "${CMAKE_ARGS[@]}" && make -j8
BUILD_STATUS=$?

# --- Summary ---
echo "==========================================="
if [ ${BUILD_STATUS} -eq 0 ]; then
    echo "Build SUCCEEDED at: $(date)"
else
    echo "Build FAILED at: $(date) (exit code: ${BUILD_STATUS})"
fi
echo "Full log saved to: ${LOG_FILE}"

exit ${BUILD_STATUS}