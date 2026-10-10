# Terrain MLMG final test matrix (2026-10-10)

| Area | Evidence from this round | Status |
|---|---|---|
| G8 continuous/restart, unchanged layout | One-rank CPU DOUBLE/FFT TerrainHill; continuous versus restart; common plot outputs bitwise equal | Pass |
| G8 changed BoxArray | Restart from 8/8/8 to 16/16/8 compared with uninterrupted runtime 8→16 reblock at matching physical time, coverage, dt and plot schedules; `plt`, `pltU`, `pltV`, `pltW` at t=0.5, 1.0, 1.5 all bitwise equal | Pass |
| G8 first state event | `RemakeLevel`/`FillPatchFineLevel` changes valid `zvel` first (Linf 0.2471481141, L2 1.336395862); initial velocity projection reduces divergence Linf 0.0261532→1.265e-8 | Diagnosed; no restart correction indicated |
| PR #4253 control | Isolated source patch on/off gave exact same `plt00002` fields; failing difference unchanged | Excluded as cause |
| MPI build | Fresh 3D DOUBLE CPU/FFT MPI build completed; `erf_exec`, `fcompare`, `fextrema` linked; two-rank hostname launch had passed | Build pass; registered CTest suite not run |
| GMRES residual | Existing feature/integrated runs show estimated residual ~7e-9 with true residual ~1e-3–1e-2; instrumented same-data snapshot build did not finish | Blocked; cause and clean-upstream status unestablished |
| Rank parity, boundary, map-factor, conservation | Registered cases available in MPI build, but CTest launch could not proceed because the shared host repeatedly failed sandbox namespace creation (`ENOSPC`) | Not run this round |
| GPU | No NVIDIA device/toolchain available on node | Not available |
| CPU single/no-FFT, upstream development regression, unit tests | No completed fresh results in this round | Not run |

Overall verdict: G8 restart parity is cleared when the continuous control follows the same remake/regrid lifecycle. Merge readiness remains blocked by unresolved GMRES attribution and the unrun final regression matrix. No numeric tolerances or defaults were changed; GMRES+FFT remains the default.
