#!/bin/bash
#SBATCH --account=erf
#SBATCH --partition=gpu-h100s
#SBATCH --job-name=terrainopt
#SBATCH --output=slurm-%j.out
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --exclude=x3101c0s33b0n0

set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
export LD_LIBRARY_PATH="/opt/cray/libfabric/1.15.2.0/lib64:$LD_LIBRARY_PATH"
export MPICH_GPU_SUPPORT_ENABLED=0
export CUDA_LAUNCH_BLOCKING=0
export KOKKOS_ENABLE_CUDA=1
export MPICH_SMP_SINGLE_COPY_MODE=NONE
export FI_CXI_DISABLE_HOST_REGISTER=1
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

EXECUTABLE="$SLURM_SUBMIT_DIR/erf_exec"

{
    date -Is
    hostname
    printf 'SLURM_JOB_ID=%s\n' "$SLURM_JOB_ID"
    printf 'SLURM_JOB_NODELIST=%s\n' "$SLURM_JOB_NODELIST"
    sha256sum "$EXECUTABLE"
    ldd "$EXECUTABLE"
    nvidia-smi --query-gpu=name,driver_version --format=csv
} > runtime_environment.txt 2>&1
if rg -q 'not found' runtime_environment.txt; then
    cat runtime_environment.txt
    exit 127
fi

srun --cpu-bind=cores -n 1 "$EXECUTABLE" inputs_canopy
