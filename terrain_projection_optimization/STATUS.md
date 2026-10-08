# Current experiment status

Updated 2026-10-07 21:14 MDT.

## Remote branch

The work is on wang1202/ERF:terrain_opt, based on remote development at
`e50b1611bc987a5188160dfa5aa5061232080859`. The branch includes the requested
agent instructions, copied `MyBuildTerrainOpt/cmake.sh`, the AMReX
MLTerrainPoisson pin, and reproducible benchmark decks. The experimental
solver remains opt-in; GMRES+FFT is the default.

## Builds

The clean Release CUDA build succeeded. The corrected candidate also compiled
and linked incrementally after replacing the unsupported MLMG flux query
with `MLTerrainPoisson::compFlux`, which calculates the terrain-aware face
fluxes needed for ERF's momentum correction.

- Control executable SHA256: `0c4b71b2441b1d0e4ed78216e8ce3e8e5aadab546a35d6380336ea82402d482b`
- Corrected candidate executable SHA256: `9a6eaf213c4f596c6301ed103e77fcc3ca1768e60e2e704f5dd1a7a7b11b2101`

## Benchmark results

Control job 18949177 completed successfully on one H100 (30m41s Slurm
elapsed), covering model time 7200 to 9000 s in 2195 steps. Mature-half
measurements from `benchmarks/summarize_logs.py`:

- 0.81736 s/step and 3632 simulated seconds per wall hour.
- Median 19 GMRES iterations; solver work averaged 0.67414 s/step, or 82.48%
of measured step time.
- Median relative GMRES residual 8.20e-9. Maximum post-projection divergence
was 1.79e-9 in L-infinity and 1.26e-7 in the unnormalized L2 norm.
- No warnings or errors were found in the completed control log.

Initial candidate job 18949428 failed after its first projection. MLMG
completed two cycles with an absolute residual of 3.17e-9, then the generic
`MLMG::getFluxes` call reached AMReX's unimplemented
`MLLinOp::getFluxes` abort. No post-projection divergence or candidate timing
is available from that attempt. Its local log, preflight, backtrace, and
partial outputs are preserved in `benchmarks/candidate_mlmg/attempts/18949428/`.

The flux call is fixed in source commit `33e6b4c`, and corrected candidate job
18953921 has been resubmitted with the new executable. It is currently pending
on scheduler priority. Thus non-FFT performance and numerical equivalence are
still unknown. The solver's homogeneous Neumann flux behavior still requires
explicit boundary-flux validation.

## Next steps

1. Monitor corrected candidate job 18953921 through H100 preflight and
completion.
2. Compare mature step and solver time, residuals, cycles, and post-projection
divergence using `benchmarks/summarize_logs.py`.
3. Compare native profiles, surface diagnostics, plotfile fields, and boundary
flux behavior before accepting any speedup.
4. Extend the run for statistical validation only if the candidate is both
faster and numerically acceptable.


## 30-minute monitoring

The detached tmux session `terrainopt-monitor` polls Slurm every 1,800 seconds for corrected candidate job 18953921 and appends scheduler snapshots to `terrain_projection_optimization/monitoring/job_18953921.jsonl`. The first check at 2026-10-08 05:48 UTC found the job pending for priority. When the job reaches a terminal state, the watcher waits for its log to settle, runs the existing control-versus-candidate performance parser, and writes a JSON comparison and Markdown readout in the monitoring folder. Field statistics and boundary-flux validation still require review after the run. The cluster blocks user crontabs, so this persistent tmux session is used for the recurring checks.


## Field-analysis pipeline and current queue

At 2026-10-08 06:10 UTC, Slurm still reported corrected candidate job 18953921 as PENDING for Priority, with a provisional start estimate of 2026-10-08 12:25:43 UTC; the estimate can move. The primary tmux monitor remains active and polls every 30 minutes.

The control output coverage is verified: plotfiles at 7200, 7800.133, 8400.469, and 9000 s; native profiles cover 7212.7978983.311 s; surface diagnostics cover 72009000 s. Parsed input decks differ only by the intended `erf.terrain_poisson_solver` selection (default control versus `mlmg` candidate).

A new comparison tool reuses the existing terrain analysis code for time-weighted native profiles and surface proxies, physical-volume AGL/terrain field statistics, PDFs, and spectra. Its control self-check is Slurm job 18954682, currently RUNNING on a shared CPU node; the first plotfield completed. Final matched field analysis is job 18954684 on the shared partition (8 CPUs, 64 GB), submitted with `afterok:18953921`, so it will run only if the candidate succeeds. A second detached tmux watcher records and publishes its terminal analysis result.

Source inspection confirms the MLTerrainPoisson flux kernel explicitly sets Neumann boundary-face flux to zero; runtime projection/divergence results remain pending. The field comparison is limited to the 30-minute single trajectory and is descriptive, not a statistical-equivalence test.

## Monitoring update — 2026-10-08 06:23 UTC

- Corrected MLMG candidate job 18953921 is still PENDING for scheduler priority. The current provisional start estimate is 2026-10-08 12:25:43 UTC and may move.
- Control field-pipeline self-check 18954682 is RUNNING (20m31s elapsed at this check). Its log confirms plotfields at 7200, 7800.133, and 8400.469 s have been processed; the 9000 s snapshot and final self-comparison are pending.
- Candidate field comparison 18954684 remains PENDING on `afterok:18953921`.
- Both detached 30-minute tmux watchers are alive. Their latest polls were 06:12:24 UTC (candidate) and 06:11:28 UTC (field analysis); next polls are due around 06:42 and 06:41 UTC.
- No corrected-candidate performance, residual, divergence, or field-difference result exists yet. The self-check is still validating the analysis pipeline, not measuring the MLMG candidate.

### Candidate integrity and QA resource check — 2026-10-08 06:28 UTC

The queued candidate executable still hashes to `9a6eaf213c4f596c6301ed103e77fcc3ca1768e60e2e704f5dd1a7a7b11b2101`; the control executable hashes to `0c4b71b2441b1d0e4ed78216e8ce3e8e5aadab546a35d6380336ea82402d482b`. A fresh input diff shows only `erf.terrain_poisson_solver` differs (`mlmg` in the candidate). Slurm job 18953921 still points to the candidate directory and requests one H100, 16 CPUs, and a 4-hour limit. At 06:28 UTC, control self-check 18954682 had run 25m06s, remained RUNNING, and had logged three of four plot snapshots; its latest accounting showed about 8.9 GiB resident memory, within its 64 GiB allocation. Candidate performance and numerical results remain unavailable until the queued run executes.

### Control field self-check progress — 2026-10-08 06:32 UTC

The control-only field self-check 18954682 has now processed all four plot snapshots through 9000 s. At 06:32:58 UTC it was still RUNNING after 29m50s, so the profile/statistics self-comparison has not yet produced a terminal result. Slurm accounting showed active CPU and disk reads, with peak resident memory about 8.9 GiB against 64 GiB requested. Candidate job 18953921 remained PENDING for Priority with a provisional 12:25:43 UTC start; field comparison 18954684 remained PENDING on its after-success dependency. No candidate performance or numerical result is available yet.

### Self-check aggregation stage — 2026-10-08 06:35 UTC

At 06:35:40 UTC, control self-check 18954682 remained RUNNING at 32m35s elapsed. Its log shows the control-side field aggregation reached all four plot times and the second aggregation (the identical control input supplied as the candidate for the self-test) has begun at 7200 s. The self-check has not yet returned its zero-difference summaries or exited. Candidate 18953921 remains PENDING for Priority; no candidate performance or numerical metrics are available.

### 30-minute scheduler and field-QA checkpoint — 2026-10-08 06:42 UTC

The candidate watcher recorded job 18953921 PENDING for Priority at 06:42:25 UTC; its provisional start estimate moved earlier from 12:25:43 to 12:10:00 UTC. The field-analysis watcher recorded job 18954684 PENDING on its `afterok:18953921` dependency at 06:41:28 UTC. Control self-check 18954682 remained RUNNING at 39m46s; its second, identical-input field aggregation had processed snapshots at 7200 and 7800.133 s. Slurm accounting showed active CPU and disk reads, about 9.1 GiB peak resident memory, and no error/exit status. Candidate performance, residual, divergence, and matched-field results remain pending.

### Control self-check progress — 2026-10-08 06:48 UTC

At 06:48:15 UTC, QA job 18954682 was RUNNING at 45m14s. Its identical-input aggregation pass had processed the 7200, 7800.133, and 8400.469 s plot snapshots; the final 9000 s snapshot and self-comparison summaries remain pending. Slurm accounting showed about 1h03m CPU time, 9.1 GB peak RSS, and increasing disk reads, within its 8-CPU/64-GB allocation. Corrected candidate 18953921 remained PENDING for Priority with a 12:10 UTC provisional start; dependent field analysis 18954684 remained pending. No candidate numerical or performance result is available.

### Field-analysis self-check failure and fix — 2026-10-08 07:01 UTC

Control pipeline QA job 18954682 ended FAILED after 51 minutes, having processed both passes over all four control plot snapshots. The traceback exposed a schema error in the new comparison wrapper: `field_statistics.aggregate()` returns time-aggregated spectra with `shell_power` and `high_wavenumber_fraction`, not raw-snapshot `total_power` and `high_wavenumber_power`. Commit `de8c24a` fixes this by summing shell power for total power, using the provided high-wavenumber fraction, and checking bin lengths. A focused identical-spectrum/schema-mismatch fixture and Python compilation pass. Full QA retry job 18954895 was submitted and is PENDING on the shared partition as of 07:01 UTC. Candidate 18953921 remains pending for priority; candidate field comparison 18954684 remains gated on its success.

### Self-check mode correction and retry — 2026-10-08 07:05 UTC

The first retry, job 18954895, failed in 4 seconds because the self-check intentionally supplies the same control directory twice while production comparison validation requires the control/MLMG solver-setting difference. Commit `3c52355` adds an explicit `--self-compare` mode that requires identical parsed inputs; the self-check launcher now uses it, while production mode still requires exactly the expected solver-setting difference. Synthetic tests pass for both modes and the corrected aggregate-spectrum summary, including rejection of schema mismatches; `py_compile`, `bash -n`, and `git diff --check` pass. Full QA retry 18954908 started on the shared node x1008c0s0b0n1 and was RUNNING at 25 seconds as of 07:05:20 UTC. Candidate 18953921 remains PENDING for Priority, with a provisional 12:40 UTC start; dependent field comparison 18954684 remains pending.

### 30-minute poll and corrected QA retry — 2026-10-08 07:12 UTC

The candidate watcher recorded job 18953921 PENDING for Priority at 07:12:26 UTC; its current provisional start estimate is 12:40 UTC. The field watcher recorded job 18954684 PENDING on its `afterok:18953921` dependency at 07:11:28 UTC. Corrected self-check 18954908 was RUNNING at 7m31s on x1008c0s0b0n1, with no immediate error; its log shows the first control plotfield at 7200 s entering field processing. Previous self-check failures and fixes are detailed above. No candidate performance or numerical output is available yet.

### Corrected field QA progress — 2026-10-08 07:19 UTC

Self-check retry 18954908 was RUNNING at 13m59s as of 07:19:11 UTC. The corrected log has advanced through plot snapshots at 7200 and 7800.133 s, with no reported error; Slurm accounting showed about 19m15s CPU time and 8.9 GB peak RSS within the 8-CPU/64-GB allocation. Candidate 18953921 remains PENDING for Priority with a 12:40 UTC provisional start. Field comparison 18954684 remains dependency-pending. The full self-check and all candidate performance/numerical results are still pending.

### Corrected field self-check progress — 2026-10-08 07:26 UTC

Self-check retry 18954908 was RUNNING at 21m21s on x1008c0s0b0n1. Its corrected identical-input pass had processed 7200, 7800.133, and 8400.469 s snapshots; the final 9000 s snapshot and zero-difference outputs are pending. Slurm accounting showed about 29m21s CPU, 9.1 GB peak RSS, and 9.4 GB disk reads within the allocation. Candidate 18953921 remains priority-pending with a provisional 12:40 UTC start; field comparison 18954684 remains dependency-pending. A terminal-state watcher is being added for this self-check so its final QA outcome is published automatically.

### Corrected field self-check: first aggregate pass complete — 2026-10-08 07:32 UTC

At 07:32:48 UTC, retry job 18954908 was RUNNING at 27m36s and had logged all four control snapshots (7200, 7800.133, 8400.469, and 9000 s) in the first identical-input aggregation pass; the second pass was not yet complete. Final regional/PDF/spectral comparisons and the success marker had not yet appeared. Slurm accounting showed active CPU and disk reads, about 9.1 GB peak RSS, and no reported error. Candidate 18953921 remained PENDING for Priority; field comparison 18954684 remained dependency-pending. The detached self-check watcher will publish its terminal result.


### Corrected field self-check: second pass started — 2026-10-08 07:38 UTC

At 07:38:35 UTC, job 18954908 remained RUNNING at 33m24s. The log confirms the first aggregate pass processed all four snapshots, then the second pass started at 7200 s; the remaining three second-pass snapshots and final self-comparison are pending. Slurm accounting showed active CPU and disk reads, about 9.1 GB peak RSS, with no reported error. This corrects the initial 07:32 wording, which incorrectly described both passes as complete. Candidate 18953921 remains priority-pending and field comparison 18954684 remains dependency-pending.

### 30-minute poll and second self-check pass — 2026-10-08 07:43 UTC

Candidate watcher poll at 07:42:26 UTC found job 18953921 still PENDING for Priority; the current provisional start remains 12:40 UTC. Field watcher poll at 07:41:28 UTC found job 18954684 still PENDING on its after-success dependency. Corrected self-check 18954908 was RUNNING at 38m06s and has begun its second identical-input aggregate pass, processing the 7200 s plotfield; the other second-pass snapshots and final zero-difference outputs remain pending. Slurm accounting showed about 53m49s CPU, 9.1 GB peak RSS, and 15.1 GB disk reads, with no reported error. The previously overstated 07:32 pass count is corrected in the prior status section.


### Control field self-check job 18954908

State: COMPLETED; exit code: 0:0; elapsed: 00:50:31; recorded 2026-10-08T07:57:21+00:00 UTC.
Self-comparison pass marker present: True.
Slurm output: `terrain_projection_optimization/monitoring/field-selfcheck-18954908.out`.

Final log excerpt:

```text
Fields processed: 7200.0 /kfs2/projects/erf/aaronwang/ERF/TerrainOptRemote/terrain_projection_optimization/benchmarks/baseline_gmres_fft/plt08876
Fields processed: 7800.133192763082 /kfs2/projects/erf/aaronwang/ERF/TerrainOptRemote/terrain_projection_optimization/benchmarks/baseline_gmres_fft/plt09615
Fields processed: 8400.468817454577 /kfs2/projects/erf/aaronwang/ERF/TerrainOptRemote/terrain_projection_optimization/benchmarks/baseline_gmres_fft/plt10344
Fields processed: 9000.0 /kfs2/projects/erf/aaronwang/ERF/TerrainOptRemote/terrain_projection_optimization/benchmarks/baseline_gmres_fft/plt11071
Fields processed: 7200.0 /kfs2/projects/erf/aaronwang/ERF/TerrainOptRemote/terrain_projection_optimization/benchmarks/baseline_gmres_fft/plt08876
Fields processed: 7800.133192763082 /kfs2/projects/erf/aaronwang/ERF/TerrainOptRemote/terrain_projection_optimization/benchmarks/baseline_gmres_fft/plt09615
Fields processed: 8400.468817454577 /kfs2/projects/erf/aaronwang/ERF/TerrainOptRemote/terrain_projection_optimization/benchmarks/baseline_gmres_fft/plt10344
Fields processed: 9000.0 /kfs2/projects/erf/aaronwang/ERF/TerrainOptRemote/terrain_projection_optimization/benchmarks/baseline_gmres_fft/plt11071
/tmp/terrainopt-field-selfcheck.aQmeij/field_analysis.md
Control self-comparison passed; generated outputs in temporary directory and removed.
```


### Self-check terminal record and publisher fix — 2026-10-08 07:58 UTC

Slurm job 18954908 completed at 07:56:46 UTC with exit code 0 after 50m31s. Its log confirms both identical-input passes processed all four plot snapshots and printed the pipeline success marker after writing a temporary `field_analysis.md`; the self-check trap removed temporary analysis files, so this verifies end-to-end pipeline completion but does not retain numeric zero-difference CSVs. The terminal watcher initially staged the event and status but did not commit because it passed an ignored `field-selfcheck-*.out` path to ordinary `git add`. The watcher now force-adds that explicit output file, and the terminal event, status, and Slurm log are being committed directly.

### Monitoring checkpoint — 2026-10-08 08:04 UTC

Corrected candidate job 18953921 remains PENDING for scheduler priority; Slurm's provisional start estimate is 12:25:43 UTC and may move. Dependent native-field comparison 18954684 remains PENDING on `afterok:18953921`. Both detached 30-minute tmux watchers are active; their latest recorded polls are 07:42:26 UTC and 07:41:28 UTC, with the next checks due around 08:12 and 08:11 UTC. Control pipeline self-check 18954908 completed successfully (exit 0, 50m31s); its temporary self-comparison outputs were removed after the end-to-end pass. Candidate performance, residual, divergence, and matched-field results remain unavailable.
