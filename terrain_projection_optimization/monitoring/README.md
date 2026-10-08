# Terrain LES monitoring and analysis

The detached tmux session `terrainopt-monitor` polls candidate Slurm job 18953921 every 1,800 seconds and appends snapshots to `job_18953921.jsonl`. When it reaches a terminal state, the watcher parses matched performance, solver, and post-projection divergence metrics against control job 18949177.

Native field comparisons use `compare_field_runs.py` and the existing terrain statistics code. They cover exact common support for native profiles and surface diagnostics, physical-volume AGL/terrain region means, PDFs, and spectra. The run is a 30-minute single trajectory, so these outputs are descriptive and cannot establish statistical equivalence. Cartesian SGS momentum recovery retains the documented cell-centered approximation.

The control-only pipeline check uses `field_analysis_selfcheck.sh` (Slurm job 18954682). Full candidate field analysis is queued as job 18954684 by `field_analysis_job.sh` with dependency `afterok:18953921`, on the shared CPU partition (8 CPUs, 64 GB). The detached `terrainopt-field-monitor` session polls that analysis job every 30 minutes, appends its state, and commits/pushes the report and CSV artifacts to `terrain_opt` when it terminates.

The cluster blocks user crontabs, so detached tmux sessions provide the recurring checks while the login environment remains available. Inspect them with `tmux capture-pane -pt terrainopt-monitor` and `tmux capture-pane -pt terrainopt-field-monitor`.
