# Terrain LES monitoring

The detached tmux session `terrainopt-monitor` runs `monitor_slurm_job.py` for corrected candidate job 18953921. It polls both `squeue` and `sacct` every 1,800 seconds, preserving each observation in `job_18953921.jsonl`. At job completion or failure it waits for the Slurm output to settle, runs `benchmarks/summarize_logs.py` against control job 18949177, and saves `job_18953921_analysis.json` plus a short Markdown readout. These automated metrics cover performance, solver residuals/iterations, and post-projection divergence; field-statistics and boundary-flux checks remain a separate scientific review.

Inspect the live session with `tmux capture-pane -pt terrainopt-monitor`. The cluster denies user crontabs, so the detached tmux session provides the recurring poll while the login environment remains available.
