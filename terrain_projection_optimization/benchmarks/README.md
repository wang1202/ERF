# Controlled terrain Poisson comparison

Both cases start from the same read-only strict-tolerance checkpoint at model time 7200 s and stop at 9000 s. The common input preserves the 360 x 360 x 200 terrain LES setup, the 1e-8 relative and absolute Poisson tolerances, and its existing physics and forcing. The only input difference is the terrain_poisson_solver = mlmg setting in the candidate deck; the baseline uses the default GMRES+FFT path.

The checkpoint provenance is copied from the existing validation campaign. The checkpoint itself remains at its original path and is not modified. Each case writes outputs to its own directory. Submit from the desired case directory with sbatch ../job.sh; do not run both cases concurrently on the same H100 reservation.

The job script uses one H100 and the same MPI/CUDA runtime settings as the established validation runs. After the baseline completes, record its executable hash and step timings before building the experimental selector. Compare the candidate using the same executable revision and otherwise identical setup. Boundary flux differences and LES-field validation are required before interpreting any speedup.

Each case runs a private copy of its executable, with the SHA256 recorded by the Slurm script. This keeps the control binary fixed while the candidate is rebuilt.
