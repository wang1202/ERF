#!/bin/bash
#SBATCH --account=erf
#SBATCH --partition=shared
#SBATCH --job-name=TerrainFieldCompare
#SBATCH --output=/kfs2/projects/erf/aaronwang/ERF/TerrainOptRemote/terrain_projection_optimization/monitoring/field-analysis-%j.out
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail
export OMP_NUM_THREADS=8
export MPLBACKEND=Agg
PYTHON=/home/wang1202/.conda-envs/metpy/bin/python3
ROOT=/kfs2/projects/erf/aaronwang/ERF/TerrainOptRemote/terrain_projection_optimization
BASELINE="$ROOT/benchmarks/baseline_gmres_fft"
CANDIDATE="$ROOT/benchmarks/candidate_mlmg"
OUTPUT="$ROOT/monitoring/job_18953921_fields"
"$PYTHON" "$ROOT/monitoring/compare_field_runs.py" --control "$BASELINE" --candidate "$CANDIDATE" --output "$OUTPUT"
echo "Matched native-field analysis completed at $(date -Is)"
