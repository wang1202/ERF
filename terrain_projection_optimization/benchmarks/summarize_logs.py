#!/usr/bin/env python3
"""Summarize matched ERF terrain-projection logs using only Python's stdlib."""
import argparse
import json
import re
import statistics
from pathlib import Path

N = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
P = {
    "end": re.compile(rf"Coarse STEP (\d+) ends\. TIME = ({N}) DT = ({N})"),
    "wall": re.compile(rf"Timestep time = ({N}) seconds\."),
    "gmres": re.compile(rf"GMRES: iter = (\d+), residual = ({N}), ({N})"),
    "gmres_time": re.compile(rf"GMRES: Solve Time = ({N})"),
    "solve": re.compile(rf"Time in solve ({N})"),
    "mg_init": re.compile(rf"MLMG: Initial residual \(resid0\) = ({N})"),
    "mg_final": re.compile(rf"MLMG: Final Iter\. (\d+) resid, resid/(?:resid0|bnorm) = ({N}), ({N})"),
    "mg_cycles": re.compile(r"MLTerrainPoisson iterations: (\d+)"),
    "div_before": re.compile(rf"Max/L2 norm of divergence before solve(?: in subdomain \d+)? at level \d+ : ({N}) ({N})"),
    "div_after": re.compile(rf"Max/L2 norm of divergence after\s+solve at level \d+ : ({N}) ({N})"),
}


def parse(path):
    lines = Path(path).read_text(errors="replace").splitlines()
    steps, solves, gmres, mlmg, divs = [], [], [], [], []
    current = end = g_pending = d_pending = None
    mg_pending = {}
    for line in lines:
        if m := re.search(r"Coarse STEP (\d+) starts", line):
            current = int(m.group(1))
        if m := P["end"].search(line):
            end = (int(m.group(1)), float(m.group(2)), float(m.group(3)))
        elif m := P["wall"].search(line):
            if end:
                steps.append({"step": end[0], "time_s": end[1], "dt_s": end[2], "wall_s": float(m.group(1))})
                end = None
        if m := P["gmres"].search(line):
            g_pending = {"step": current, "iterations": int(m.group(1)),
                         "residual": float(m.group(2)), "relative_residual": float(m.group(3))}
        elif m := P["gmres_time"].search(line):
            if g_pending:
                gmres.append({**g_pending, "solve_s": float(m.group(1))})
                g_pending = None
        if m := P["solve"].search(line):
            solves.append({"step": current, "solve_s": float(m.group(1))})
        if m := P["mg_init"].search(line):
            mg_pending[current] = {"initial_residual": float(m.group(1))}
        if m := P["mg_final"].search(line):
            mg_pending.setdefault(current, {}).update({"cycles": int(m.group(1)),
                "final_residual": float(m.group(2)), "relative_residual": float(m.group(3))})
        if m := P["mg_cycles"].search(line):
            mlmg.append({"step": current, "cycles": int(m.group(1)), **mg_pending.pop(current, {})})
        if m := P["div_before"].search(line):
            d_pending = {"step": current, "before_linf": float(m.group(1)), "before_l2": float(m.group(2))}
        elif m := P["div_after"].search(line):
            if d_pending and d_pending["step"] == current:
                divs.append({**d_pending, "after_linf": float(m.group(1)), "after_l2": float(m.group(2))})
                d_pending = None
    if not steps:
        raise ValueError(f"No complete coarse-step timing records found: {path}")
    mature = steps[len(steps)//2:]
    ids = {r["step"] for r in mature}
    wall, sim = sum(r["wall_s"] for r in mature), sum(r["dt_s"] for r in mature)
    ss = [r for r in solves if r["step"] in ids]
    gg = [r for r in gmres if r["step"] in ids]
    mm = [r for r in mlmg if r["step"] in ids]
    dd = [r for r in divs if r["step"] in ids]
    med = lambda a: statistics.median(a) if a else None
    return {
        "log": str(path), "steps": len(steps),
        "model_time_start_s": steps[0]["time_s"]-steps[0]["dt_s"],
        "model_time_end_s": steps[-1]["time_s"], "mature_steps": len(mature),
        "mature_step_wall_s_mean": statistics.mean(r["wall_s"] for r in mature),
        "mature_step_wall_s_median": med([r["wall_s"] for r in mature]),
        "mature_mean_dt_s": statistics.mean(r["dt_s"] for r in mature),
        "mature_sim_seconds_per_wall_hour": sim/wall*3600 if wall else None,
        "solver": {
            "solve_calls": len(ss), "solve_wall_s_per_step": sum(r["solve_s"] for r in ss)/len(mature),
            "solve_time_fraction_of_step_wall": sum(r["solve_s"] for r in ss)/wall if wall else None,
            "gmres_calls": len(gg), "gmres_iterations_median": med([r["iterations"] for r in gg]),
            "gmres_relative_residual_median": med([r["relative_residual"] for r in gg]),
            "mlmg_calls": len(mm), "mlmg_cycles_median": med([r["cycles"] for r in mm]),
            "mlmg_relative_residual_median": med([r["relative_residual"] for r in mm if "relative_residual" in r]),
            "projection_pairs": len(dd),
            "post_projection_max_divergence_linf": max((r["after_linf"] for r in dd), default=None),
            "post_projection_max_divergence_l2": max((r["after_l2"] for r in dd), default=None),
        },
        "warnings_and_errors": [s.strip() for s in lines if re.search(
            r"\b(:NaN|Inf|infinity|WARNING|Warning|Abort|Segmentation|CUDA error)\b", s, re.I)],
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("baseline_log", type=Path)
    ap.add_argument("candidate_log", type=Path)
    ap.add_argument("--output", type=Path, help="also write JSON to this path")
    args = ap.parse_args()
    a, b = parse(args.baseline_log), parse(args.candidate_log)
    a_step, b_step = a["mature_step_wall_s_mean"], b["mature_step_wall_s_mean"]
    a_solve, b_solve = a["solver"]["solve_wall_s_per_step"], b["solver"]["solve_wall_s_per_step"]
    data = {
        "comparison": {
            "matched_model_interval_s": [a["model_time_start_s"], a["model_time_end_s"]],
            "same_model_interval": abs(a["model_time_start_s"]-b["model_time_start_s"]) < 1e-8 and
                                   abs(a["model_time_end_s"]-b["model_time_end_s"]) < 1e-8,
            "step_speedup_candidate_vs_baseline": a_step/b_step if b_step else None,
            "solve_speedup_candidate_vs_baseline": a_solve/b_solve if b_solve else None,
            "interpretation": "Performance only; review divergence and LES statistics before accepting a speedup.",
        }, "baseline": a, "candidate": b,
    }
    out = json.dumps(data, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(out)
    print(out, end="")


if __name__ == "__main__":
    main()
