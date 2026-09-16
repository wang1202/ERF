# ERF instrument simulators

ERF provides column-based forward operators for three ground-based remote-sensing instruments. They are enabled independently with `do_sim = true` beneath the corresponding `erf.*` input prefix.

| Input prefix | Instrument | Output |
| --- | --- | --- |
| `erf.ceilometer` | Vaisala CL-31, 905 nm | attenuated backscatter and cloud-base height |
| `erf.doppler_lidar` | Halo StreamLine XR+, 1550 nm | vertical-stare velocity/CNR, VAD wind retrieval, and attenuated backscatter |
| `erf.mwr` | Radiometrics MP-3000A, K+V band | thermodynamic/moisture profiles, LWP/PWV, and optional brightness temperatures |

## PBL diagnostic recorder

`erf.pbl_recorder.do_record = true` writes moving-column native-SHOC
diagnostics to `pbl_diag_00_summary.txt` and `pbl_diag_00_profile.txt` (one
pair per configured site). It requires native SHOC and
`erf.shoc.diagnose_entrainment = true`. Configure `level`, matching zero-based
`i_loc`/`j_loc`, `summary_output_interval`, `write_profiles`,
`profile_output_interval`, and `output_file`; the `move_*` schedule uses the
same physical trajectory convention as the instrument simulators above.

The recorder samples Eulerian diagnostics at the instantaneous moving cell;
it never differentiates boundary-layer height along the moving trajectory.
Files are self-describing whitespace-delimited scientific output. Summary
and profile rows use physical `z_phys` heights and `-999` for unavailable or
guarded values. Restart initialization appends without duplicating the last
time row, and each completed write is flushed.

## Configuration

All simulators require matching, zero-based `i_loc` and `j_loc` arrays; multiple entries configure multiple sites. Common optional parameters are `output_interval` (a positive number of coarse steps) and `output_file` (the output-file prefix). The ceilometer additionally accepts `backscatter_model`, `r_eff_default = 1.0e-5`, `r_eff_ice_default = 3.0e-5`, and `beta_thresh`. The Doppler lidar accepts `backscatter_model` and `r_eff_default = 1.0e-5`. `simple` is the implemented geometric-optics backscatter model; `mie` currently falls back to `simple` with a diagnostic message.

The Doppler lidar supports `do_vertical_stare`, `do_vad`, `vad_elevation_angle`, `do_backscatter`, and either a qualitative `cnr_model = normalized` or the parameterized `physical` CNR model. The physical model accepts `pulse_energy`, `optics_efficiency`, `beam_radius`, `bandwidth`, and `focus_range`. Set `erf.mwr.do_brightness_temp = true` to write the 35-channel MP-3000A brightness-temperature approximation.

## Moving columns

Each prefix can specify equal-length `move_start_time`, `move_stop_time`, `move_speed_x`, and `move_speed_y` arrays. A schedule moves every configured site from its original cell center at the listed piecewise-constant velocities. Motion accumulates over later intervals; an interval's final displacement remains in effect after it stops. Intervals must be ordered and non-overlapping. The instantaneous physical `(x, y)` is converted to the sampled `(i_loc, j_loc)` cell and clamped at the domain boundary, so motion beyond the domain remains safe while its reported physical position continues to reflect the requested trajectory.

Every data row starts with `time x y i_loc j_loc`. Every profile `height` field
is the level's cell-centered physical `z_phys` value, not the computational
coordinate height. The remaining fields are:

- Ceilometer: `height att_backscatter cloud_base_height`
- Doppler stare: `height w CNR`; VAD: `height u v wind_speed wind_dir`; backscatter: `height beta_att`
- MWR profile: `height T qv qc`; MWR column: `LWP PWV`; MWR TB: `frequency TB`

## Microphysics and limitations

Bulk microphysics reads `RhoQ1_comp`, `RhoQ2_comp`, and `RhoQ3_comp` as vapor, cloud liquid, and cloud ice. The simple lidar calculation applies fixed effective-radius defaults and geometric-optics extinction. With Super-Droplet microphysics, particle radii and multiplicities produce column extinction only on actual lidar write steps. The generic size-resolved interface also accepts spectral-bin extinction data, but this ERF branch does not expose an SBM bin container to populate it.

The MWR TB calculation is a lightweight plane-parallel approximation of water-vapor, oxygen, and cloud-liquid absorption. It is intended for simulator diagnostics rather than a replacement for a validated operational radiative-transfer retrieval. Terrain-following vertical geometry, finite beam volume, horizontal beam steering, and instrument noise/calibration effects are likewise outside the current column operators.
