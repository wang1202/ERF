# Controlled terrain Poisson comparison

Both cases start from the same read-only strict-tolerance checkpoint at model time 7200 s and stop at 9000 s. The common input preserves the 360 x 360 x 200 terrain LES setup, the 1e-8 relative and absolute Poisson tolerances, and its existing physics and forcing. The only input difference is the terrain_poisson_solver = mlmg setting in the candidate deck; the baseline uses the default GMRES+FFT path.

The checkpoint provenance is copied from the existing validation campaign. The checkpoint itself remains at its original path and is not modified. Each case writes outputs to its own directory. Submit from the desired case directory with sbatch ../job.sh.

The job script requests one H100 and uses the same MPI/CUDA runtime settings as the established validation runs. Both executable copies are pinned and hashed; the candidate was built from the same CMake tree and AMReX revision. The cases write to separate directories and can run on separate H100 nodes.

Each case runs a private executable copy, with its SHA256 recorded by the Slurm preflight. When both logs exist, summarize matched interval throughput, solver time, iteration/cycle counts, residuals, and post-projection divergence with:

```bash
python3 summarize_logs.py baseline_gmres_fft/slurm-18949177.out candidate_mlmg/slurm-18953921.out --output log_metrics.json
```

The parser uses the later half of each run for mature timing, matching the established terrain-validation metrics. Native `surf`, `mean`, `flux`, and `subgrid` diagnostics plus plotfiles support field and LES-statistics comparisons. A measured speedup remains provisional until divergence, boundary-flux behavior, velocity and thermodynamic fields, turbulence statistics, and surface proxies are reviewed.

The first candidate attempt, job 18949428, exposed an unsupported generic flux
call after the first solve. The source now calls MLTerrainPoisson's supported
`compFlux` method; the corrected rerun is job 18953921. The failed attempt's
artifacts are retained locally under `candidate_mlmg/attempts/18949428/`.
