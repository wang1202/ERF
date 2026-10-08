#!/usr/bin/env python3
"""Poll the dependency-gated field analysis every 30 minutes and publish each poll."""
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


def delay_to_next_poll(path, interval):
    try:
        last = json.loads(path.read_text().splitlines()[-1])["observed_utc"]
        due = dt.datetime.fromisoformat(last) + dt.timedelta(seconds=interval)
        return max(0.0, (due - dt.datetime.now(dt.timezone.utc)).total_seconds())
    except (OSError, IndexError, KeyError, ValueError, json.JSONDecodeError):
        return 0.0


def poll():
    rc, raw, err = run(["sacct", "-n", "-P", "-j", JOB,
                        "--format=JobID,State,Elapsed,ExitCode,NodeList"])
    state = None
    for line in raw.splitlines():
        f = line.split("|")
        if f and f[0].strip() == JOB and len(f) > 1:
            state = f[1].split()[0].split("+")[0]
            break
    return ({"observed_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
             "job_id": JOB, "state": state, "sacct": raw, "sacct_rc": rc,
             "sacct_error": err}, state)


def stage_commit_push(paths, message, label):
    paths = [str(Path(p).relative_to(REPO)) for p in paths]
    rc, out, err = run(["git", "-C", str(REPO), "add", "--", *paths])
    if rc:
        print(label + " stage failed: " + (err or out), flush=True)
        return
    _, staged_text, _ = run(["git", "-C", str(REPO), "diff", "--cached", "--name-only"])
    staged = set(staged_text.splitlines())
    if not staged.issubset(set(paths)):
        run(["git", "-C", str(REPO), "reset", "--", *paths])
        print(label + " skipped commit because unrelated paths were staged.", flush=True)
        return
    if staged:
        rc, out, err = run(["git", "-C", str(REPO), "commit", "-m", message])
        if rc:
            print(label + " commit failed: " + (err or out), flush=True)
            return
        print(label + " committed: " + out, flush=True)
    rc, out, err = run(["timeout", "35s", "git", "-C", str(REPO), "push",
                        "origin", "terrain_opt"], timeout=45)
    print(label + " push exit=" + str(rc) + " " + (out or err), flush=True)


def publish_terminal(event):
    for _ in range(60):
        try:
            st = LOG.stat()
            if st.st_size >= 0:
                break
        except FileNotFoundError:
            pass
        time.sleep(10)
    state = event.get("state")
    success = state == "COMPLETED" and REPORT.is_file()
    marker = "### Native-field analysis job " + JOB
    text = STATUS.read_text()
    if marker not in text:
        detail = (REPORT.read_text() if REPORT.is_file()
                  else "Field-analysis report was not created; inspect " + str(LOG))
        section = ("\n\n" + marker + "\n\n" + "State: " + str(state) +
                   ". Recorded UTC: " + event["observed_utc"] + ".\n\n" + detail + "\n")
        STATUS.write_text(text + section)
    paths = [EVENTS, LOG, STATUS]
    if success:
        paths.append(REPORT.parent)
    stage_commit_push(paths, "docs: record terrain LES field analysis",
                      "Field-analysis terminal result state=" + str(state))


def main():
    MON.mkdir(parents=True, exist_ok=True)
    delay = delay_to_next_poll(EVENTS, INTERVAL)
    print("Monitoring field-analysis job " + JOB + " every 30 minutes.", flush=True)
    if delay:
        print("Preserving poll cadence; next check in " + str(round(delay)) + " seconds.", flush=True)
        time.sleep(delay)
    while True:
        event, state = poll()
        with EVENTS.open("a") as f:
            f.write(json.dumps(event, sort_keys=True) + "\n")
            f.flush()
        print(event["observed_utc"] + " state=" + str(state or "UNKNOWN"), flush=True)
        if state in END:
            publish_terminal(event)
            return
        stage_commit_push([EVENTS], "docs: record terrain LES field-analysis poll",
                          "Field-analysis poll " + event["observed_utc"])
        time.sleep(INTERVAL)


if __name__ == "__main__":
    main()
