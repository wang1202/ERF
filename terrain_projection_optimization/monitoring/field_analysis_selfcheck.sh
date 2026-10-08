#!/bin/bash
#SBATCH --account=erf
#SBATCH --partition=shared
#SBATCH --job-name=TerrainFieldQA
#SBATCH --output=/kfs2/projects/erf/aaronwang/ERF/TerrainOptRemote/terrain_projection_optimization/monitoring/field-selfcheck-%j.out
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail
export OMP_NUM_THREADS=8
export MPLBACKEND=Agg
TEST_OUT=$(mktemp -d /tmp/terrainopt-field-selfcheck.XXXXXX)
trap 'rm -rf "$TEST_OUT"' EXIT
PYTHON=/home/wang1202/.conda-envs/metpy/bin/python3
ROOT=/kfs2/projects/erf/aaronwang/ERF/TerrainOptRemote/terrain_projection_optimization
BASELINE="$ROOT/benchmarks/baseline_gmres_fft"
"$PYTHON" "$ROOT/monitoring/compare_field_runs.py" --control "$BASELINE" --candidate "$BASELINE" --output "$TEST_OUT" --self-compare
echo "Control self-comparison passed; generated outputs in temporary directory and removed."
