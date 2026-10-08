#!/usr/bin/env python3
"""Poll ERF job 18953921 every 30 minutes and save terminal performance results."""
import datetime as dt, hashlib, json, subprocess, sys, time
from pathlib import Path
R=Path("/kfs2/projects/erf/aaronwang/ERF/TerrainOptRemote")
B=R/"terrain_projection_optimization/benchmarks"
J="18953921"; I=1800; C=B/"candidate_mlmg"; L=C/("slurm-"+J+".out")
BASE=B/"baseline_gmres_fft/slurm-18949177.out"
O=R/"terrain_projection_optimization/monitoring"; O.mkdir(parents=True,exist_ok=True)
E=O/("job_"+J+".jsonl"); F=O/("job_"+J+"_analysis.json"); M=O/("job_"+J+"_analysis.md")
END={"COMPLETED","FAILED","CANCELLED","TIMEOUT","NODE_FAIL","OUT_OF_MEMORY","PREEMPTED","BOOT_FAIL","DEADLINE","REVOKED"}

def run(args):
 try:
  p=subprocess.run(args,text=True,capture_output=True,timeout=180)
  return p.returncode,p.stdout.strip(),p.stderr.strip()
 except Exception as e: return 127,"",str(e)

def poll():
 qr,q,qe=run(["squeue","-h","-j",J,"-o","%i|%T|%M|%N|%R"])
 ar,a,ae=run(["sacct","-n","-P","-j",J,"--format=JobID,State,Elapsed,ExitCode,NodeList"])
 state=None
 for line in a.splitlines():
  f=line.split("|")
  if f[0].strip()==J and len(f)>1: state=f[1].split()[0].split("+")[0]; break
 if state is None and q:
  f=q.splitlines()[0].split("|"); state=f[1] if len(f)>1 else None
 return {"observed_utc":dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),"job_id":J,"state":state,"squeue":q,"squeue_rc":qr,"squeue_error":qe,"sacct":a,"sacct_rc":ar,"sacct_error":ae},state

def finish(event):
 old=None; same=0
 for _ in range(80):
  try:
   s=L.stat(); now=(s.st_size,s.st_mtime_ns)
  except FileNotFoundError: now=None
  same=same+1 if now is not None and now==old else 0
  if same>=3: break
  old=now; time.sleep(15)
 r={"job_id":J,"final_state":event["state"],"final_scheduler_snapshot":event,"candidate_log":str(L),"baseline_log":str(BASE),"recorded_utc":dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}
 exe=C/"erf_exec"
 if exe.is_file():
  h=hashlib.sha256()
  with exe.open("rb") as f:
   for block in iter(lambda:f.read(8*1024*1024),b""): h.update(block)
  r["candidate_executable_sha256"]=h.hexdigest()
 env=C/"runtime_environment.txt"
 if env.is_file(): r["runtime_environment"]=env.read_text(errors="replace")
 if BASE.is_file() and L.is_file():
  rc,out,err=run([sys.executable,str(B/"summarize_logs.py"),str(BASE),str(L),"--output",str(F)])
  r["parser_exit_code"]=rc; r["parser_stderr"]=err
  if rc==0 and F.is_file(): r["performance_comparison"]=json.loads(F.read_text())
  else: r["performance_analysis_error"]=err or out or "Parser produced no output"
 else: r["performance_analysis_error"]="A benchmark log is missing"
 F.write_text(json.dumps(r,indent=2,allow_nan=False)+"\n")
 cmp=r.get("performance_comparison",{}).get("comparison",{}); cand=r.get("performance_comparison",{}).get("candidate",{}); sol=cand.get("solver",{})
 lines=["# Corrected terrain MLMG run "+J,"","- Final scheduler state: "+str(r["final_state"]),"","- Recorded UTC: "+r["recorded_utc"],"- Executable SHA256: "+str(r.get("candidate_executable_sha256","unavailable"))]
 if cand:
  lines += ["- Mature-window step time: "+str(cand.get("mature_step_wall_s_mean"))+" s/step.","- Step speedup vs control: "+str(cmp.get("step_speedup_candidate_vs_baseline"))+"x.","- Solver speedup vs control: "+str(cmp.get("solve_speedup_candidate_vs_baseline"))+"x.","- Median GMRES iterations: "+str(sol.get("gmres_iterations_median"))+".","- Median MLMG cycles: "+str(sol.get("mlmg_cycles_median"))+".","- Maximum post-projection divergence L-infinity: "+str(sol.get("post_projection_max_divergence_linf"))+".","- Maximum post-projection divergence unnormalized L2: "+str(sol.get("post_projection_max_divergence_l2"))+".","- Field statistics and boundary-flux checks remain required."]
 else: lines += ["","Performance parsing did not complete:",str(r.get("performance_analysis_error","unknown"))]
 M.write_text("\n".join(lines)+"\n")
 publish_candidate()

def publish_candidate():
 files=[
  "terrain_projection_optimization/monitoring/job_"+J+".jsonl",
  "terrain_projection_optimization/monitoring/job_"+J+"_analysis.json",
  "terrain_projection_optimization/monitoring/job_"+J+"_analysis.md",
 ]
 rc,out,err=run(["git","-C",str(R),"add","--",*files])
 if rc:
  print("Could not stage candidate results: "+(err or out),flush=True); return
 _,names,_=run(["git","-C",str(R),"diff","--cached","--name-only"])
 if not set(names.splitlines()).issubset(set(files)):
  run(["git","-C",str(R),"reset","--",*files])
  print("Skipped candidate result commit due to unrelated staged files.",flush=True); return
 rc,out,err=run(["git","-C",str(R),"commit","-m","docs: record corrected terrain MLMG run"])
 if rc:
  print("Candidate result commit failed: "+(err or out),flush=True); return
 rc,out,err=run(["git","-C",str(R),"push","origin","terrain_opt"])
 print("Candidate result push exit="+str(rc)+" "+(out or err),flush=True)

def main():
 print("Monitoring job "+J+" every 30 minutes; events: "+str(E),flush=True)
 while True:
  event,state=poll()
  with E.open("a") as f:
   f.write(json.dumps(event,sort_keys=True)+"\n")
   f.flush()
  print(event["observed_utc"]+" state="+str(state or "UNKNOWN"),flush=True)
  if state in END:
   finish(event)
   print("Final results written to "+str(F),flush=True)
   return
  time.sleep(I)

if __name__ == "__main__":
 main()
