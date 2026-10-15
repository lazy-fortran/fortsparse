#!/usr/bin/env python3
"""Optional stdlib-only fresh-process benchmark orchestration; not library code."""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import signal
import statistics
import subprocess
import sys
import time


def physical_cpus():
    cores = {}
    for cpu in sorted(os.sched_getaffinity(0)):
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        key = tuple((topology / name).read_text().strip() for name in
                    ("physical_package_id", "core_id"))
        cores.setdefault(key, cpu)
    return list(cores.values())


def sample(args, cpus, backend, family, level, threads, workload, trial):
    binary = args.bin_dir / ("bench_spd_pardiso" if backend == "mkl-pardiso"
                             else "bench_spd")
    # taskset allows Python's vfork/exec launch path. A Python preexec_fn would
    # fork the interpreter and contaminate the child's inherited RSS high-water.
    command = ["taskset", "--cpu-list", ",".join(map(str, cpus[:threads])),
               str(binary)]
    if backend != "mkl-pardiso":
        command.append(backend)
    command += [family, str(level), str(threads), workload,
                "heldout" if args.heldout else "training"]
    record = dict(backend=backend, family=family, level=level, threads=threads,
                  workload=workload, trial=trial, command=command)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS=str(threads), MKL_DYNAMIC="FALSE",
               MKL_INTERFACE_LAYER="LP64", MKL_THREADING_LAYER="GNU",
               OMP_DYNAMIC="FALSE", OMP_PROC_BIND="close", OMP_PLACES="cores")
    if not binary.exists():
        record.update(status="unavailable", reason="executable missing")
        return record
    begin = time.monotonic()
    process = subprocess.Popen(command, env=env, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               start_new_session=True)
    try:
        stdout, stderr = process.communicate(timeout=1.0)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        stdout, stderr = process.communicate()
        record.update(status="timeout", stdout=stdout, stderr=stderr)
        return record
    record["process_wall_s"] = time.monotonic() - begin
    if process.returncode:
        record.update(status="failed", returncode=process.returncode,
                      stdout=stdout, stderr=stderr)
        return record
    lines = stdout.splitlines()
    header = next((i for i, line in enumerate(lines)
                   if line.startswith("backend,family,")), None)
    if header is None:
        record.update(status="invalid_output", stdout=stdout, stderr=stderr)
        return record
    try:
        values = next(csv.DictReader(io.StringIO("\n".join(lines[header:]))))
        for key, value in values.items():
            if key in ("backend", "family", "reused"):
                record[key] = value
            else:
                record[key] = float(value)
        record["factor_solve_s"] = record["factor_s"] + record["solve_s"]
        record["workflow_cpu_s"] = record["parent_cpu_s"] + record["helper_cpu_s"]
        record["status"] = "ok"
    except (ValueError, StopIteration, TypeError) as error:
        record.update(status="invalid_output", reason=str(error), stdout=stdout)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--backends", default="superlu,umfpack,mumps,native,mkl-pardiso")
    parser.add_argument("--threads", default="1,2,4,8,16")
    parser.add_argument("--trials", type=int, default=5)
    parser.add_argument("--workloads", default="rhs1,rhs8,refactor,derivatives,warm")
    parser.add_argument("--families", default="band,grid2,grid3")
    parser.add_argument("--levels", default="1,2,3")
    parser.add_argument("--heldout", action="store_true")
    parser.add_argument("--calibrate-ladder", action="store_true")
    parser.add_argument("--scaling-corpus", type=Path)
    args = parser.parse_args()
    args.bin_dir = args.bin_dir.resolve()
    cpus = physical_cpus()
    threads = [int(t) for t in args.threads.split(",") if int(t) <= len(cpus)]
    if not threads or min(threads) < 1 or args.trials < 1:
        parser.error("positive available thread counts and trials required")
    backends = args.backends.split(",")
    cases = [(f, int(level)) for f in args.families.split(",")
             for level in args.levels.split(",")]
    if args.calibrate_ladder:
        backends = [b for b in backends if b != "native"]
        if not backends:
            parser.error("calibration requires external backends")
        cases = [(family, level) for family in ("grid2", "grid3")
                 for level in range(4, 8)]
    if args.scaling_corpus:
        cases += [tuple(case) for case in json.loads(
            args.scaling_corpus.read_text())["selected_scaling_cases"]]
    workloads = ["rhs1"] if args.calibrate_ladder else args.workloads.split(",")
    count = 1 if args.calibrate_ladder else args.trials
    records = []
    for family, level in cases:
        for backend in backends:
            for thread in threads:
                for workload in workloads:
                    for trial in range(count):
                        records.append(sample(args, cpus, backend, family,
                                              level, thread, workload, trial))
    selected = []
    if args.calibrate_ladder:
        for family in ("grid2", "grid3"):
            accepted = []
            for level in range(4, 8):
                group = [r for r in records if r["family"] == family
                         and r["level"] == level]
                if all(r["status"] == "ok" and r["process_wall_s"] < 0.5
                       for r in group):
                    accepted.append(level)
            if accepted:
                selected.append([family, max(accepted)])
    summaries = []
    keys = sorted({(r["backend"], r["family"], r["level"], r["threads"],
                    r["workload"]) for r in records})
    metrics = ("factor_s", "solve_s", "refactor_s", "total_s", "parent_cpu_s",
               "parent_rss_kib", "helper_rss_kib", "helper_cpu_s",
               "backward_error", "process_wall_s", "factor_solve_s",
               "workflow_cpu_s", "jvp_s", "vjp_s", "cold_factor_s",
               "factor_nnz", "factor_mem_mb", "analysis_s", "numeric_s")
    for key in keys:
        group = [r for r in records if (r["backend"], r["family"], r["level"],
                 r["threads"], r["workload"]) == key]
        good = [r for r in group if r["status"] == "ok"]
        summary = dict(zip(("backend", "family", "level", "threads", "workload"), key))
        summary.update(successes=len(good), samples=len(group))
        if good:
            for metric in metrics:
                values = [r[metric] for r in good]
                available = [v for v in values if v >= 0]
                summary[metric] = (dict(median=statistics.median(available),
                                        minimum=min(available), maximum=max(available))
                                   if available else None)
        summaries.append(summary)
    for summary in summaries:
        serial = next((s for s in summaries if s["backend"] == summary["backend"]
                       and s["family"] == summary["family"]
                       and s["level"] == summary["level"]
                       and s["workload"] == summary["workload"]
                       and s["threads"] == 1 and s["successes"] == s["samples"]), None)
        if serial and summary["successes"] == summary["samples"]:
            for metric in ("total_s", "factor_solve_s"):
                speedup = serial[metric]["median"] / summary[metric]["median"]
                summary[metric]["speedup"] = speedup
                summary[metric]["parallel_efficiency"] = speedup / summary["threads"]
    binaries = {}
    for name in ("bench_spd", "bench_spd_pardiso"):
        binary = args.bin_dir / name
        if binary.exists():
            binaries[name] = dict(sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                                 ldd=subprocess.run(["ldd", str(binary)], text=True,
                                                    capture_output=True).stdout)
    revision = subprocess.run(["git", "rev-parse", "HEAD"], text=True,
                              capture_output=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain"], text=True,
                           capture_output=True).stdout.strip()
    cpu_description = next((line.split(":", 1)[1].strip() for line in
                            Path("/proc/cpuinfo").read_text().splitlines()
                            if line.startswith("model name")), "unknown")
    result = dict(revision=revision, platform=platform.platform(),
                  physical_cpu_ids=cpus, binary_provenance=binaries,
                  selected_scaling_cases=selected, records=records, summaries=summaries,
                  calibration_only=args.calibrate_ladder, command=sys.argv,
                  working_tree_dirty=bool(dirty), cpu_model=cpu_description,
                  fixed_environment=dict(OPENBLAS_NUM_THREADS="1", MKL_DYNAMIC="FALSE",
                                         MKL_INTERFACE_LAYER="LP64", MKL_THREADING_LAYER="GNU",
                                         OMP_DYNAMIC="FALSE", OMP_PROC_BIND="close",
                                         OMP_PLACES="cores"))
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"{len(records)} samples; {sum(r['status'] == 'ok' for r in records)} passed; "
          f"results: {args.output}")
    if args.calibrate_ladder:
        print(f"Frozen scaling candidates accepted: {selected}")


if __name__ == "__main__":
    main()
