#!/usr/bin/env python3
"""Check mass-integral parity and face-momentum ownership in rank-parity logs."""
import argparse
import math
import re
from pathlib import Path

MASS_RE = re.compile(r"TIME=\s*([+-]?[0-9.eE+-]+).*?MASS\s+SL/ML\s*=\s*([+-]?[0-9.eE+-]+)\s+([+-]?[0-9.eE+-]+)", re.S)
AUDIT_RE = re.compile(
    r"ERF_TERRAIN_SHARED_FACE_AUDIT level=([0-9]+) component=(xmom|ymom|zmom) "
    r"stage=(pre_correction_momentum|correction_flux|post_correction_momentum) "
    r"difference_linf=([+-]?[0-9.eE+-]+) tolerance=([+-]?[0-9.eE+-]+)"
)

def read_log(path):
    text = Path(path).read_text(encoding="utf-8")
    masses = [(float(t), float(ml)) for t, _sl, ml in MASS_RE.findall(text)]
    audits = [(int(lev), comp, stage, float(diff), float(tol))
              for lev, comp, stage, diff, tol in AUDIT_RE.findall(text)]
    return text, masses, audits

def close(a, b, rtol, atol):
    return abs(a-b) <= atol + rtol*max(abs(a), abs(b))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--one-log", required=True)
    parser.add_argument("--two-log", required=True)
    parser.add_argument("--relative-tolerance", type=float, required=True)
    parser.add_argument("--absolute-tolerance", type=float, required=True)
    args = parser.parse_args()
    one_text, one_mass, one_audits = read_log(args.one_log)
    two_text, two_mass, two_audits = read_log(args.two_log)
    if len(one_mass) < 2 or len(two_mass) < 2:
        raise SystemExit("need at least initial and final multilevel mass integrals in both logs")
    if len(one_mass) != len(two_mass):
        raise SystemExit(f"mass-integral sample counts differ: one={len(one_mass)} two={len(two_mass)}")
    for index, ((t1, m1), (t2, m2)) in enumerate(zip(one_mass, two_mass)):
        if not math.isclose(t1, t2, rel_tol=0.0, abs_tol=1.e-14):
            raise SystemExit(f"mass-integral times differ at sample {index}: {t1} vs {t2}")
        if not close(m1, m2, args.relative_tolerance, args.absolute_tolerance):
            raise SystemExit(f"multilevel mass differs at t={t1}: one={m1:.17g} two={m2:.17g}")
    for label, audits in (("one-rank", one_audits), ("two-rank", two_audits)):
        if not audits:
            raise SystemExit(f"{label} log has no shared-face ownership audit")
        for level, component, stage, difference, tolerance in audits:
            if not math.isfinite(difference) or difference > tolerance:
                raise SystemExit(f"{label} shared {component} level {level} at {stage} differs by {difference:g} > {tolerance:g}")
        required = {(level, component, stage)
                    for level in (0, 1)
                    for component in ("xmom", "ymom", "zmom")
                    for stage in ("pre_correction_momentum", "correction_flux", "post_correction_momentum")}
        present = {(level, component, stage) for level, component, stage, _, _ in audits}
        if not required.issubset(present):
            raise SystemExit(f"{label} ownership audit missing pairs: {sorted(required-present)}")
    print(f"rank parity passed: {len(one_mass)} multilevel mass samples agree; "
          f"ownership audits={len(one_audits)}/{len(two_audits)}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
