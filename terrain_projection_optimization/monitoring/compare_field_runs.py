#!/usr/bin/env python3
"""Describe matched native LES profiles and field statistics for two ERF runs."""
import argparse
import datetime as dt
import hashlib
import json
import math
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

ANALYSIS = Path("/kfs2/projects/erf/aaronwang/ERF/MyBuild/Exec/TerrainToleranceValidation")
sys.path.insert(0, str(ANALYSIS))
import field_statistics as fs
import statistics_core as sc
import validation_campaign as vc

def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def scrub(value):
    if isinstance(value, dict):
        return {k: scrub(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [scrub(v) for v in value]
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    return value

def write_csv(path, rows):
    if not rows:
        Path(path).write_text("no_rows\n")
        return
    sc.write_csv(path, scrub(rows))

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--control", type=Path, required=True)
    ap.add_argument("--candidate", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    control, candidate, out = args.control.resolve(), args.candidate.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    control_params = vc.params(control / "inputs_canopy")
    candidate_params = vc.params(candidate / "inputs_canopy")
    parameter_differences = {
        key: [control_params.get(key), candidate_params.get(key)]
        for key in sorted(control_params.keys() | candidate_params.keys())
        if control_params.get(key) != candidate_params.get(key)
    }
    expected_differences = {"erf.terrain_poisson_solver": [None, "mlmg"]}
    if parameter_differences != expected_differences:
        raise ValueError("Unexpected control/candidate input differences: " + str(parameter_differences))

    # Native profile diagnostics are independently time-weighted on their actual common support.
    cp, npf = sc.profiles([control]), sc.profiles([candidate])
    pstart = max(float(cp["t"][0]), float(npf["t"][0]))
    pend = min(float(cp["t"][-1]), float(npf["t"][-1]))
    if not pstart < pend:
        raise ValueError("Native profile support does not overlap")
    profile_rows, profile_summary = [], []
    if not set(cp["series"]).issubset(npf["series"]):
        raise ValueError("Native profile schemas differ")
    for quantity in cp["series"]:
        a = sc.average(cp["t"], cp["series"][quantity], pstart, pend)
        b = sc.average(npf["t"], npf["series"][quantity], pstart, pend)
        delta = b - a
        ref_rms = sc.rms(a)
        profile_summary.append({
            "quantity": quantity, "window_start_s": pstart, "window_end_s": pend,
            "reference_profile_rms": ref_rms,
            "candidate_profile_rms": sc.rms(b),
            "difference_profile_rms": sc.rms(delta),
            "relative_difference_profile_rms": sc.relative(sc.rms(delta), ref_rms),
            "maximum_absolute_height_difference": float(abs(delta).max()),
            "sampling_note": "One 1800 s trajectory; descriptive only, no independent block-variability estimate.",
        })
        for z, av, bv, dv in zip(cp["z"], a, b, delta):
            profile_rows.append({
                "quantity": quantity, "z_m": float(z), "reference_mean": float(av),
                "candidate_mean": float(bv), "signed_difference": float(dv),
                "absolute_difference": float(abs(dv)),
                "relative_difference": sc.relative(float(dv), float(av)),
                "window_start_s": pstart, "window_end_s": pend,
            })
    write_csv(out / "native_profile_comparison.csv", profile_rows)
    write_csv(out / "native_profile_summary.csv", profile_summary)

    # Surface proxies use the common native time support and do not represent W/m2 or vector stress.
    sstart = max(float(cp["surface_t"][0]), float(npf["surface_t"][0]))
    send = min(float(cp["surface_t"][-1]), float(npf["surface_t"][-1]))
    if not sstart < send:
        raise ValueError("Surface diagnostic support does not overlap")
    surface_rows, surface_summary = sc.compare_surface(cp, npf, sstart, send)
    write_csv(out / "surface_comparison.csv", surface_rows)
    write_csv(out / "surface_summary.csv", surface_summary)

    # Native AMReX plots: physical-volume canopy/terrain statistics, PDFs, and spectra.
    ci = fs.plotfiles([str(control)])
    ni = fs.plotfiles([str(candidate)])
    if not ci or not ni:
        raise ValueError("A run has no readable plotfiles")
    fstart = max(float(ci[0]["time_s"]), float(ni[0]["time_s"]))
    fend = min(float(ci[-1]["time_s"]), float(ni[-1]["time_s"]))
    if not fstart < fend:
        raise ValueError("Plotfile time support does not overlap")
    cache = Path(tempfile.mkdtemp(prefix="terrainopt-field-cache-"))
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            cf = fs.aggregate(ci, vc.params(control / "inputs_canopy"), cache / "control", fstart, fend)
            nf = fs.aggregate(ni, vc.params(candidate / "inputs_canopy"), cache / "candidate", fstart, fend)
            region_rows, region_blocks, region_summary = fs.compare(cf, nf, fstart, fend)
    finally:
        shutil.rmtree(cache, ignore_errors=True)
    write_csv(out / "field_region_comparison.csv", region_rows)
    write_csv(out / "field_region_blocks.csv", region_blocks)
    write_csv(out / "field_region_summary.csv", region_summary)

    pdf_rows, pdf_summary = [], []
    cpdf = {(x["region"], x["quantity"]): x for x in cf["pdfs"]}
    npdf = {(x["region"], x["quantity"]): x for x in nf["pdfs"]}
    if cpdf.keys() != npdf.keys():
        raise ValueError("PDF schema differs")
    for key, a in cpdf.items():
        b = npdf[key]
        if a["edges"] != b["edges"]:
            raise ValueError("PDF bins differ")
        pa, pb = a["probabilities"], b["probabilities"]
        for lo, hi, x, y in zip(a["edges"][:-1], a["edges"][1:], pa, pb):
            pdf_rows.append({"region": key[0], "quantity": key[1], "bin_low": lo,
                             "bin_high": hi, "reference_probability": x,
                             "candidate_probability": y, "difference": y - x})
        tvd = 0.5 * (sum(abs(x-y) for x, y in zip(pa, pb)) +
                     abs(a["underflow_probability"]-b["underflow_probability"]) +
                     abs(a["overflow_probability"]-b["overflow_probability"]))
        row = {"region": key[0], "quantity": key[1], "total_variation_distance": tvd,
               "reference_underflow": a["underflow_probability"],
               "candidate_underflow": b["underflow_probability"],
               "reference_overflow": a["overflow_probability"],
               "candidate_overflow": b["overflow_probability"]}
        for i, label in enumerate(("q05", "q50", "q95")):
            row["reference_" + label] = a["quantiles"][i]
            row["candidate_" + label] = b["quantiles"][i]
        pdf_summary.append(row)
    write_csv(out / "field_pdf_bins.csv", pdf_rows)
    write_csv(out / "field_pdf_summary.csv", pdf_summary)

    spectral_rows, spectral_summary = [], []
    cspec = {(x["agl_height_m"], x["quantity"]): x for x in cf["spectra"]}
    nspec = {(x["agl_height_m"], x["quantity"]): x for x in nf["spectra"]}
    if cspec.keys() != nspec.keys():
        raise ValueError("Spectrum schema differs")
    for key, a in cspec.items():
        b = nspec[key]
        if a["wavenumber_cycles_per_m"] != b["wavenumber_cycles_per_m"]:
            raise ValueError("Spectral bins differ")
        for k, x, y in zip(a["wavenumber_cycles_per_m"], a["shell_power"], b["shell_power"]):
            spectral_rows.append({"agl_height_m": key[0], "quantity": key[1],
                                  "wavenumber_cycles_per_m": k,
                                  "reference_shell_power": x, "candidate_shell_power": y,
                                  "difference": y - x})
        spectral_summary.append({
            "agl_height_m": key[0], "quantity": key[1],
            "reference_total_power": a["total_power"], "candidate_total_power": b["total_power"],
            "reference_high_wavenumber_fraction": a["high_wavenumber_power"] / a["total_power"] if a["total_power"] else 0.0,
            "candidate_high_wavenumber_fraction": b["high_wavenumber_power"] / b["total_power"] if b["total_power"] else 0.0,
        })
    write_csv(out / "field_spectral_bins.csv", spectral_rows)
    write_csv(out / "field_spectral_summary.csv", spectral_summary)

    inventory = []
    for case, items in (("control", ci), ("candidate", ni)):
        for item in items:
            inventory.append({"case": case, "time_s": item["time_s"], "path": item["path"],
                              "header_sha256": item["header_sha256"]})
    write_csv(out / "plotfile_inventory.csv", inventory)
    geometry = cf["snapshots"][0]["geometry"]
    metadata = {
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "control_dir": str(control), "candidate_dir": str(candidate),
        "control_inputs_sha256": sha(control / "inputs_canopy"),
        "candidate_inputs_sha256": sha(candidate / "inputs_canopy"),
        "parameter_differences": parameter_differences,
        "native_profile_common_window_s": [pstart, pend],
        "surface_common_window_s": [sstart, send],
        "field_common_window_s": [fstart, fend],
        "control_plot_times_s": [x["time_s"] for x in ci],
        "candidate_plot_times_s": [x["time_s"] for x in ni],
        "field_snapshot_count": [len(cf["snapshots"]), len(nf["snapshots"])],
        "field_geometry": geometry,
        "SGS_stress_conversion": geometry.get("SGS_stress_conversion"),
        "sampling_limitation": "This is a 30-minute single-trajectory comparison with 600 s plot snapshots. It is descriptive, not a statistical-equivalence test; there are too few 30-minute blocks to estimate sampling variability.",
        "interpretation": "Turbulent trajectories may decorrelate. Compare time-averaged profiles, physical-volume regional statistics, PDFs, spectra, and surface proxies; do not treat instantaneous gridpoint differences as an accuracy threshold.",
    }
    (out / "field_analysis.json").write_text(json.dumps(scrub(metadata), indent=2, allow_nan=False) + "\n")
    largest = sorted(profile_summary, key=lambda x: x["difference_profile_rms"], reverse=True)[:5]
    lines = ["# Terrain LES native-field comparison", "",
             "Control job 18949177 versus corrected MLMG candidate job 18953921.", "",
             "This matched 30-minute single-trajectory window is descriptive, not a statistical-equivalence test. Turbulent decorrelation is expected; no universal acceptance margin is applied.", "",
             f"- Native profiles common support: {pstart:.3f} to {pend:.3f} s.",
             f"- Surface diagnostics common support: {sstart:.3f} to {send:.3f} s.",
             f"- Plotfield common support: {fstart:.3f} to {fend:.3f} s.",
             f"- Matched plot snapshots: {len(cf['snapshots'])} control and {len(nf['snapshots'])} candidate.",
             f"- Input parameter differences: {parameter_differences}.",
             f"- Geometry checks: top boundary error {geometry['top_boundary_error_m']:.3g} m; minimum cell thickness {geometry['minimum_dz_m']:.6g} m.", "",
             "## Largest profile RMS differences", ""]
    lines += [f"- {r['quantity']}: RMS difference {r['difference_profile_rms']:.6g}; reference profile RMS {r['reference_profile_rms']:.6g}; relative {r['relative_difference_profile_rms']}." for r in largest]
    lines += ["", "Regional means, PDFs, spectra, surface proxies, and their detailed support are in the CSV files beside this report.",
              "The Cartesian SGS momentum conversion uses Tau31/Tau32 divided by density and cell-centered h_zeta; the existing analysis documents this as a second-order approximation. These SGS diagnostics do not directly verify the Poisson face-flux boundary condition.", ""]
    (out / "field_analysis.md").write_text("\n".join(lines))
    print(str(out / "field_analysis.md"), flush=True)

if __name__ == "__main__":
    main()
