#!/usr/bin/env python3
"""Poll the dependency-gated native-field analysis and publish its terminal record."""
import datetime as dt
import json
import subprocess
import time
from pathlib import Path

REPO = Path("/kfs2/projects/erf/aaronwang/ERF/TerrainOptRemote")
JOB = "18954684"
INTERVAL = 1800
MON = REPO / "terrain_projection_optimization/monitoring"
EVENTS = MON / ("field-analysis-" + JOB + ".jsonl")
LOG = MON / ("field-analysis-" + JOB + ".out")
REPORT = MON / "job_18953921_fields/field_analysis.md"
STATUS = REPO / "terrain_projection_optimization/STATUS.md"
END = {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL",
       "OUT_OF_MEMORY", "PREEMPTED", "BOOT_FAIL", "DEADLINE", "REVOKED"}

def run(args, timeout=180):
    try:
        p = subprocess.run(args, text=True, capture_output=True, timeout=timeout)
        return p.returncode, p.stdout.strip(), p.stderr.strip()
    except Exception as e:
        return 127, "", type(e).__name__ + ": " + str(e)

def poll():
    rc, raw, err = run(["sacct", "-n", "-P", "-j", JOB,
                        "--format=JobID,State,Elapsed,ExitCode,NodeList"])
    state = None
    for line in raw.splitlines():
        f = line.split("|")
        if f and f[0].strip() == JOB and len(f) > 1:
            state = f[1].split()[0].split("+")[0]
            break
    return {"observed_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "job_id": JOB, "state": state, "sacct": raw, "sacct_rc": rc,
            "sacct_error": err}, state

def publish(event):
    for _ in range(60):
        try:
            st = LOG.stat()
            size_mtime = (st.st_size, st.st_mtime_ns)
        except FileNotFoundError:
            size_mtime = None
        if size_mtime is not None:
            break
        time.sleep(10)
    state = event.get("state")
    success = state == "COMPLETED" and REPORT.is_file()
    marker = "### Native-field analysis job " + JOB
    text = STATUS.read_text()
    if marker not in text:
        detail = (REPORT.read_text() if REPORT.is_file()
                  else "Field-analysis report was not created; inspect " + str(LOG))
        section = "\n\n" + marker + "\n\n" + "State: " + str(state) + ". Recorded UTC: " + event["observed_utc"] + ".\n\n" + detail + "\n"
        STATUS.write_text(text + section)
    paths = [
        str(EVENTS.relative_to(REPO)),
        str(LOG.relative_to(REPO)),
        str(STATUS.relative_to(REPO)),
    ]
    if success:
        paths.append(str(REPORT.parent.relative_to(REPO)))
    rc, out, err = run(["git", "-C", str(REPO), "add", "--", *paths])
    if rc:
        print("Staging final analysis failed: " + (err or out), flush=True)
        return
    _, staged_text, _ = run(["git", "-C", str(REPO), "diff", "--cached", "--name-only"])
    staged = set(staged_text.splitlines())
    if not staged.issubset(set(paths)):
        run(["git", "-C", str(REPO), "reset", "--", *paths])
        print("Skipped commit because unrelated paths were staged.", flush=True)
        return
    if staged:
        rc, out, err = run(["git", "-C", str(REPO), "commit", "-m",
                            "docs: record terrain LES field analysis"])
        if rc:
            print("Commit failed: " + (err or out), flush=True)
            return
    rc, out, err = run(["git", "-C", str(REPO), "push", "origin", "terrain_opt"])
    print("Field-analysis state=" + str(state) + " published; push exit=" + str(rc) +
          " " + (out or err), flush=True)

def main():
    MON.mkdir(parents=True, exist_ok=True)
    print("Monitoring field-analysis job " + JOB + " every 30 minutes.", flush=True)
    while True:
        event, state = poll()
        with EVENTS.open("a") as f:
            f.write(json.dumps(event, sort_keys=True) + "\n")
            f.flush()
        print(event["observed_utc"] + " state=" + str(state or "UNKNOWN"), flush=True)
        if state in END:
            publish(event)
            return
        time.sleep(INTERVAL)

if __name__ == "__main__":
    main()
