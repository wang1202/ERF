#!/usr/bin/env python3
"""Poll one Slurm field self-check and publish its terminal outcome."""
import datetime as dt
import json
import subprocess
import sys
import time
from pathlib import Path

REPO = Path("/kfs2/projects/erf/aaronwang/ERF/TerrainOptRemote")
JOB = sys.argv[1] if len(sys.argv) > 1 else "18954908"
if not JOB.isdigit():
    raise SystemExit("job ID must be numeric")
INTERVAL = 1800
MON = REPO / "terrain_projection_optimization/monitoring"
EVENTS = MON / f"field-selfcheck-{JOB}.jsonl"
LOG = MON / f"field-selfcheck-{JOB}.out"
STATUS = REPO / "terrain_projection_optimization/STATUS.md"
END = {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL",
       "OUT_OF_MEMORY", "PREEMPTED", "BOOT_FAIL", "DEADLINE", "REVOKED"}


def run(args):
    try:
        p = subprocess.run(args, text=True, capture_output=True, timeout=180)
        return p.returncode, p.stdout.strip(), p.stderr.strip()
    except Exception as exc:
        return 127, "", type(exc).__name__ + ": " + str(exc)


def poll():
    rc, raw, err = run(["sacct", "-n", "-P", "-j", JOB,
                        "--format=JobID,State,Elapsed,ExitCode,NodeList"])
    state = None
    elapsed = exit_code = node = None
    for line in raw.splitlines():
        fields = line.split("|")
        if fields and fields[0].strip() == JOB and len(fields) >= 5:
            state = fields[1].split()[0].split("+")[0]
            elapsed, exit_code, node = fields[2], fields[3], fields[4]
            break
    if state is None:
        qrc, qraw, qerr = run(["squeue", "-h", "-j", JOB, "-o", "%i|%T|%M|%N|%R"])
        if qraw:
            fields = qraw.splitlines()[0].split("|")
            state = fields[1] if len(fields) > 1 else None
        else:
            qrc, qerr = qrc, qerr
    else:
        qrc, qraw, qerr = 0, "", ""
    return ({"observed_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
             "job_id": JOB, "state": state, "elapsed": elapsed,
             "exit_code": exit_code, "node": node,
             "sacct": raw, "sacct_rc": rc, "sacct_error": err,
             "squeue": qraw, "squeue_rc": qrc, "squeue_error": qerr}, state)


def publish(event):
    marker = f"### Control field self-check job {JOB}"
    current = STATUS.read_text()
    if marker not in current:
        try:
            output = LOG.read_text(errors="replace").splitlines()
        except FileNotFoundError:
            output = []
        passed = event.get("state") == "COMPLETED" and any(
            "Control self-comparison passed" in line for line in output)
        section = ["", "", marker, "",
                   f"State: {event.get('state')}; exit code: {event.get('exit_code')}; elapsed: {event.get('elapsed')}; recorded {event['observed_utc']} UTC.",
                   f"Self-comparison pass marker present: {passed}.",
                   f"Slurm output: `{LOG.relative_to(REPO)}`.", "", "Final log excerpt:", "", "```text"]
        section.extend(output[-30:] or ["No Slurm output was captured."])
        section.extend(["```", ""])
        STATUS.write_text(current + "\n".join(section))
    paths = [str(EVENTS.relative_to(REPO)), str(STATUS.relative_to(REPO))]
    rc, out, err = run(["git", "-C", str(REPO), "add", "--", *paths])
    if rc:
        print("Could not stage self-check result: " + (err or out), flush=True)
        return
    if LOG.is_file():
        log_path = str(LOG.relative_to(REPO))
        rc, out, err = run(["git", "-C", str(REPO), "add", "-f", "--", log_path])
        if rc:
            print("Could not stage self-check log: " + (err or out), flush=True)
            return
        paths.append(log_path)
    _, staged_text, _ = run(["git", "-C", str(REPO), "diff", "--cached", "--name-only"])
    staged = set(staged_text.splitlines())
    if not staged.issubset(set(paths)):
        run(["git", "-C", str(REPO), "reset", "--", *paths])
        print("Skipped self-check commit due to unrelated staged files.", flush=True)
        return
    rc, out, err = run(["git", "-C", str(REPO), "commit", "-m",
                        f"docs: record field self-check {JOB}"])
    if rc:
        print("Self-check commit failed: " + (err or out), flush=True)
        return
    rc, out, err = run(["git", "-C", str(REPO), "push", "origin", "terrain_opt"])
    print("Self-check result published; push exit=" + str(rc) + " " + (out or err), flush=True)


def main():
    MON.mkdir(parents=True, exist_ok=True)
    print(f"Monitoring field self-check {JOB} every {INTERVAL} seconds.", flush=True)
    while True:
        event, state = poll()
        with EVENTS.open("a") as stream:
            stream.write(json.dumps(event, sort_keys=True) + "\n")
            stream.flush()
        print(event["observed_utc"] + " state=" + str(state or "UNKNOWN"), flush=True)
        if state in END:
            publish(event)
            return
        time.sleep(INTERVAL)


if __name__ == "__main__":
    main()
