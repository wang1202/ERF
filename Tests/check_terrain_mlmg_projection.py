#!/usr/bin/env python3
"""Strict validation for every terrain MLMG projection event in an ERF log."""
import argparse
import math
import re
import sys
import unittest
from dataclasses import dataclass

NUMBER = r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
EVENT_RE = re.compile(
    r"^ERF_TERRAIN_PROJECTION_EVENT step=(\S+) projection_call=([0-9]+) "
    r"time=(\S+) dt=(\S+) level=([0-9]+) solver=(\S+) regions=([0-9]+) "
    r"pre_Linf=(\S+) post_Linf=(\S+) post_max_cell=(\([^\s]+\)) "
    r"post_max_rank=([0-9]+) relative_reduction=(\S+) "
    r"solve_status=([+-]?[0-9]+) solver_reltol=(\S+) solver_abstol=(\S+) "
    r"residual=(\S+) compatibility_mean_max_abs=(\S+)$"
)
PRE_RE = re.compile(
    rf"^Max/L2 norm of divergence before\s+solve in subdomain [0-9]+ "
    rf"at level ([0-9]+) : ({NUMBER})"
)
POST_RE = re.compile(
    rf"^Max/L2 norm of divergence after\s+solve at level ([0-9]+) : (\S+)"
)
FLUX_RE = re.compile(
    rf"^Terrain MLMG Neumann normal correction flux max: level=([0-9]+) "
    rf"subdomain=([0-9]+) dir=([xyz]) side=(lo|hi) face_index=(-?[0-9]+)"
    rf".* samples=([0-9]+) max_abs=(\S+)$"
)

@dataclass
class Event:
    step: str
    time: float
    dt: float
    level: int
    regions: int
    pre: float
    post: float
    relative: float
    status: int
    residual: float
    compatibility: float
    projection_call: int
    post_max_cell: str
    post_max_rank: int
    solver: str
    solver_reltol: float
    solver_abstol: float
    coarse_step: str | None = None

def finite_number(value, context):
    try:
        number = float(value)
    except ValueError as exc:
        raise ValueError(f"{context}: malformed number {value!r}") from exc
    if not math.isfinite(number):
        raise ValueError(f"{context}: nonfinite number {value!r}")
    return number

def parse_log(text, expected_levels, min_events, final_step, divergence_tolerance,
              required_faces=(), neumann_flux_tolerance=1.e-12):
    events = []
    human_posts = {level: [] for level in expected_levels}
    latest_pre = {}
    flux_faces = set()
    current_coarse_step = None
    for line_number, original in enumerate(text.splitlines(), 1):
        line = original.strip()
        coarse_start = re.match(r"Coarse STEP ([0-9]+) starts", line)
        if coarse_start:
            current_coarse_step = coarse_start.group(1)
        coarse_end = re.match(r"Coarse STEP ([0-9]+) ends", line)
        if "Max/L2 norm of divergence before" in line:
            match = PRE_RE.match(line)
            if match:
                level = int(match.group(1))
                latest_pre[level] = finite_number(match.group(2), f"line {line_number} pre L_inf")
        if "Max/L2 norm of divergence after" in line:
            match = POST_RE.match(line)
            if not match:
                raise ValueError(f"line {line_number}: malformed post-projection norm: {line}")
            level = int(match.group(1))
            value = finite_number(match.group(2), f"line {line_number} post L_inf")
            if level not in expected_levels:
                raise ValueError(f"line {line_number}: unexpected divergence level {level}")
            if value > divergence_tolerance:
                raise ValueError(
                    f"line {line_number}: post L_inf={value:g} exceeds {divergence_tolerance:g}")
            if level not in latest_pre:
                raise ValueError(f"line {line_number}: no pre-projection norm for level {level}")
            if value > latest_pre[level]:
                raise ValueError(
                    f"line {line_number}: projection increased L_inf at level {level}: "
                    f"pre={latest_pre[level]:g}, post={value:g}")
            human_posts[level].append(value)

        if line.startswith("ERF_TERRAIN_PROJECTION_EVENT"):
            match = EVENT_RE.fullmatch(line)
            if not match:
                raise ValueError(f"line {line_number}: malformed structured event: {line}")
            (step, projection_call, time, dt, level, solver, regions, pre, post,
             max_cell, max_rank, relative, status, reltol, abstol, residual,
             compatibility) = match.groups()
            level = int(level)
            if level not in expected_levels:
                raise ValueError(f"line {line_number}: unexpected projection level {level}")
            event = Event(
                step=step,
                time=finite_number(time, f"line {line_number} time"),
                dt=finite_number(dt, f"line {line_number} dt"),
                level=level,
                regions=int(regions),
                pre=finite_number(pre, f"line {line_number} pre L_inf"),
                post=finite_number(post, f"line {line_number} post L_inf"),
                relative=finite_number(relative, f"line {line_number} relative reduction"),
                status=int(status),
                residual=finite_number(residual, f"line {line_number} solver residual"),
                compatibility=finite_number(compatibility, f"line {line_number} compatibility mean"),
                coarse_step=current_coarse_step,
                projection_call=int(projection_call),
                post_max_cell=max_cell,
                post_max_rank=int(max_rank),
                solver=solver,
                solver_reltol=finite_number(reltol, f"line {line_number} solver reltol"),
                solver_abstol=finite_number(abstol, f"line {line_number} solver abstol"),
            )
            if event.regions < 1:
                raise ValueError(f"line {line_number}: no solve regions reported")
            if event.solver not in ("mlmg", "gmres_fft"):
                raise ValueError(f"line {line_number}: unknown solver {event.solver!r}")
            if event.solver_reltol <= 0 or event.solver_abstol < 0:
                raise ValueError(f"line {line_number}: invalid solver tolerances")
            if event.status != 0:
                raise ValueError(f"line {line_number}: solver status {event.status}")
            if event.post > divergence_tolerance:
                raise ValueError(
                    f"line {line_number}: structured post L_inf={event.post:g} exceeds "
                    f"{divergence_tolerance:g}")
            if event.post > event.pre:
                raise ValueError(
                    f"line {line_number}: structured projection increased divergence "
                    f"(pre={event.pre:g}, post={event.post:g})")
            expected_relative = (event.pre - event.post) / event.pre if event.pre else 0.0
            if not math.isclose(event.relative, expected_relative, rel_tol=2.e-5, abs_tol=2.e-7):
                raise ValueError(
                    f"line {line_number}: relative reduction {event.relative:g} does not match "
                    f"pre/post norms ({expected_relative:g})")
            events.append(event)

        if line.startswith("Terrain MLMG Neumann normal correction flux max:"):
            match = FLUX_RE.fullmatch(line)
            if not match:
                raise ValueError(f"line {line_number}: malformed Neumann flux diagnostic: {line}")
            level, _subdomain, direction, side, _face, samples, value = match.groups()
            value = finite_number(value, f"line {line_number} Neumann flux")
            if int(samples) <= 0:
                raise ValueError(f"line {line_number}: Neumann face sampled no values")
            if abs(value) > neumann_flux_tolerance:
                raise ValueError(
                    f"line {line_number}: Neumann flux {value:g} exceeds "
                    f"{neumann_flux_tolerance:g}")
            flux_faces.add(f"{level}:{direction}:{side}")
        if coarse_end:
            current_coarse_step = None

    if not re.search(rf"Coarse STEP {re.escape(str(final_step))} ends", text):
        raise ValueError(f"simulation did not complete coarse step {final_step}")
    if not events:
        raise ValueError("no structured terrain projection events found")
    call_ids = sorted(event.projection_call for event in events)
    if len(set(call_ids)) != len(call_ids):
        raise ValueError("duplicate projection_call identifier in structured events")
    if call_ids != list(range(call_ids[0], call_ids[-1] + 1)):
        raise ValueError("gap in structured projection_call identifiers")
    for level in expected_levels:
        level_events = [event for event in events if event.level == level]
        final_step_events = [event for event in level_events if event.coarse_step == str(final_step)]
        if len(final_step_events) < min_events:
            raise ValueError(f"expected at least {min_events} projection events at level {level} "
                             f"during coarse step {final_step}; found {len(final_step_events)}")
        if len(level_events) < min_events:
            raise ValueError(
                f"expected at least {min_events} projection events at level {level}; "
                f"found {len(level_events)}")
        if len(human_posts[level]) != len(level_events):
            raise ValueError(
                f"level {level}: {len(level_events)} structured events but "
                f"{len(human_posts[level])} human-readable norms")
        for index, (event, value) in enumerate(zip(level_events, human_posts[level]), 1):
            if not math.isclose(event.post, value, rel_tol=1.e-7, abs_tol=0.0):
                raise ValueError(
                    f"level {level} event {index}: structured post {event.post:g} differs "
                    f"from human norm {value:g}")
    missing_faces = sorted(set(required_faces) - flux_faces)
    if missing_faces:
        raise ValueError("required Neumann faces have no sampled flux: " + ", ".join(missing_faces))
    return events

def self_test():
    fixture = ["Coarse STEP 1 starts ..."]
    for level in (0, 1):
        for local_call in range(2):
            projection_call = 2 * level + local_call
            fixture.append(
                f"Max/L2 norm of divergence before solve in subdomain 0 at level {level} : "
                "1e-2 2e-2 and volume-weighted sum 0")
            fixture.append(
                f"ERF_TERRAIN_PROJECTION_EVENT step=1 projection_call={projection_call} time=0.5 dt=0.5 level={level} "
                "solver=mlmg regions=1 pre_Linf=1e-2 post_Linf=1e-8 post_max_cell=(0,0,0) "
                "post_max_rank=0 relative_reduction=0.999999 solve_status=0 "
                "solver_reltol=1e-6 solver_abstol=1e-10 residual=1e-9 "
                "compatibility_mean_max_abs=0")
            fixture.append(
                f"Max/L2 norm of divergence after  solve at level {level} : "
                "1e-8 2e-8 and volume-weighted sum 0")
    fixture.append("Coarse STEP 1 ends. TIME = 0.5 DT = 0.5")
    valid = "\n".join(fixture) + "\n"

    class CheckerTests(unittest.TestCase):
        def check_fails(self, text):
            with self.assertRaises(ValueError):
                parse_log(text, (0, 1), 2, 1, 1.e-6)

        def test_valid_multiple_events(self):
            self.assertEqual(len(parse_log(valid, (0, 1), 2, 1, 1.e-6)), 4)

        def test_late_bad_divergence_fails(self):
            late_bad = valid.rsplit(
                "post_Linf=1e-8 relative_reduction=0.999999", 1)
            late_bad = "post_Linf=2e-4 relative_reduction=0.98".join(late_bad)
            late_bad = late_bad.rsplit(
                " : 1e-8 2e-8 and volume-weighted sum 0", 1)
            late_bad = " : 2e-4 2e-2 and volume-weighted sum 0".join(late_bad)
            self.check_fails(late_bad)

        def test_nan_fails(self):
            self.check_fails(valid.replace("post_Linf=1e-8", "post_Linf=NaN", 1))

        def test_missing_level_fails(self):
            missing = "\n".join(
                line for line in valid.splitlines()
                if "level=1 " not in line and "at level 1 :" not in line
            ) + "\n"
            self.check_fails(missing)

        def test_missing_final_step_fails(self):
            self.check_fails(valid.replace("Coarse STEP 1 ends. TIME = 0.5 DT = 0.5\n", ""))

        def test_final_step_without_projection_fails(self):
            earlier = valid.replace("Coarse STEP 1 starts ...", "Coarse STEP 0 starts ...")
            earlier = earlier.replace("Coarse STEP 1 ends. TIME = 0.5 DT = 0.5",
                                      "Coarse STEP 0 ends. TIME = 0.5 DT = 0.5")
            self.check_fails(earlier + "Coarse STEP 1 starts ...\\n"
                             "Coarse STEP 1 ends. TIME = 1.0 DT = 0.5\\n")

        def test_small_float_threshold_is_numeric(self):
            flux = ("Terrain MLMG Neumann normal correction flux max: "
                    "level=1 subdomain=0 dir=x side=lo face_index=0 "
                    "operator_domain=[0:1,0:1,0:1] samples=3 max_abs=5e-9\n")
            self.assertEqual(
                len(parse_log(valid + flux, (0, 1), 2, 1, 1.e-6,
                              ("1:x:lo",), 1.e-8)),
                4)

        def test_small_flux_violation_fails_float_comparison(self):
            flux = ("Terrain MLMG Neumann normal correction flux max: "
                    "level=1 subdomain=0 dir=x side=lo face_index=0 "
                    "operator_domain=[0:1,0:1,0:1] samples=3 max_abs=2e-7\\n")
            with self.assertRaises(ValueError):
                parse_log(valid + flux, (0, 1), 2, 1, 1.e-6,
                          ("1:x:lo",), 1.e-8)

        def test_zero_flux_samples_fail(self):
            flux = ("Terrain MLMG Neumann normal correction flux max: "
                    "level=1 subdomain=0 dir=x side=lo face_index=0 "
                    "operator_domain=[0:1,0:1,0:1] samples=0 max_abs=0\\n")
            with self.assertRaises(ValueError):
                parse_log(valid + flux, (0, 1), 2, 1, 1.e-6,
                          ("1:x:lo",), 1.e-8)

    suite = unittest.defaultTestLoader.loadTestsFromTestCase(CheckerTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--log")
    parser.add_argument("--expected-levels", default="0,1")
    parser.add_argument("--min-events", type=int, default=2)
    parser.add_argument("--final-step", default="1")
    parser.add_argument("--divergence-tolerance", type=float, default=1.e-6)
    parser.add_argument("--required-neumann-faces", default="")
    parser.add_argument("--neumann-flux-tolerance", type=float, default=1.e-12)
    args = parser.parse_args()
    if args.self_test:
        return self_test()
    if not args.log:
        parser.error("--log is required unless --self-test is used")
    try:
        with open(args.log, encoding="utf-8") as stream:
            text = stream.read()
        levels = tuple(int(item) for item in args.expected_levels.split(",") if item)
        faces = tuple(item for item in args.required_neumann_faces.split(",") if item)
        events = parse_log(text, levels, args.min_events, args.final_step,
                           args.divergence_tolerance, faces,
                           args.neumann_flux_tolerance)
    except (OSError, ValueError) as exc:
        print(f"terrain projection validation failed: {exc}", file=sys.stderr)
        return 1
    for event in events:
        print(
            f"step={event.step} time={event.time:g} level={event.level} "
            f"pre_Linf={event.pre:g} post_Linf={event.post:g} "
            f"relative_reduction={event.relative:g} status={event.status} "
            f"residual={event.residual:g} compatibility_mean={event.compatibility:g}")
    print(f"Validated {len(events)} projection events across levels {levels}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
