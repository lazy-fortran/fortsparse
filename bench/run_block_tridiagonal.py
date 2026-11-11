#!/usr/bin/env python3
"""Fresh-process trials of bench_block_tridiagonal.

Every sample is one process pinned to the first THREADS CPUs of the allowed
set. A process that exceeds --timeout is killed with its group and kept as a
TIMEOUT row; a backend that is not built is kept as one FAILED row. Writes raw.csv
(every trial), summary.csv (median, minimum, maximum) and manifest.txt.
"""
import argparse
import csv
import hashlib
import os
import platform
import resource
import signal
import statistics
import subprocess
import sys
import time

SHAPES = [(64, 104, 3), (256, 104, 1), (64, 312, 1)]  # blocks, width, repeats
THREADS = {
    "lapack": [1, 2, 4, 8, 16],
    "supernodal": [1, 2, 4, 8, 16],
    "blocked": [1, 8],
    "mumps": [1],
    "superlu": [1],
    "umfpack": [1],
}
MEASURED = ["factor_s", "factor_cpu_s", "refactor_s", "refactor_cpu_s",
            "analysis_s", "numeric_s", "solve_s", "peak_rss_kib",
            "process_wall_s", "process_cpu_s_full"]
KEY = ["mode", "threads", "blocks", "width"]


def run_one(binary, mode, threads, blocks, width, repeats, decades, timeout):
    cpus = sorted(os.sched_getaffinity(0))[:threads]
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="false",
               OMP_PLACES="cores", OMP_PROC_BIND="close",
               OPENBLAS_NUM_THREADS=str(threads))
    command = [binary, str(blocks), str(width), str(decades), mode, str(repeats)]
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.perf_counter()
    process = subprocess.Popen(
        command, env=env, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
        text=True, start_new_session=True,
        preexec_fn=lambda: os.sched_setaffinity(0, cpus))
    try:
        output, _ = process.communicate(timeout=timeout)
        status = "OK" if process.returncode == 0 else "FAILED"
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        output, _ = process.communicate()
        status = "TIMEOUT"
    wall = time.perf_counter() - start
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    row = {"mode": mode, "threads": str(threads), "blocks": str(blocks),
           "width": str(width)}
    lines = [line for line in output.splitlines() if "," in line]
    if status == "OK" and len(lines) >= 2:
        row.update(zip(lines[-2].split(","), (v.strip() for v in lines[-1].split(","))))
    row.update(status=status, process_wall_s=f"{wall:.4f}",
               # Completed-process user plus system time, loader included.
               process_cpu_s_full=f"{after.ru_utime - before.ru_utime + after.ru_stime - before.ru_stime:.4f}",
               cpus=" ".join(map(str, cpus)), command=" ".join(command))
    return row


def write_summary(rows, path):
    groups = {}
    for row in rows:
        groups.setdefault(tuple(row[name] for name in KEY), []).append(row)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        header = KEY + ["trials", "ok", "timeout", "failed"]
        for name in MEASURED:
            header += [f"{name}_median", f"{name}_min", f"{name}_max"]
        writer.writerow(header + ["negative", "expected", "reliability",
                                  "ordering", "active_threads", "backward_error_max"])
        for key, group in groups.items():
            good = [row for row in group if row["status"] == "OK"]
            line = list(key) + [len(group), len(good),
                                sum(row["status"] == "TIMEOUT" for row in group),
                                sum(row["status"] == "FAILED" for row in group)]
            for name in MEASURED:
                values = [float(row[name]) for row in good]
                if values:
                    line += [f"{statistics.median(values):.4e}",
                             f"{min(values):.4e}", f"{max(values):.4e}"]
                else:
                    line += ["", "", ""]
            for name in ["negative", "expected", "reliability", "ordering",
                         "active_threads"]:
                line.append("/".join(sorted({row[name] for row in good})))
            line.append(max((row["backward_error"] for row in good),
                            key=float, default=""))
            writer.writerow(line)


def text_of(command):
    try:
        return subprocess.run(command, capture_output=True, text=True).stdout.strip()
    except OSError:
        return "unavailable"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--modes", default=",".join(THREADS))
    parser.add_argument("--trials", type=int, default=5)
    parser.add_argument("--timeout", type=float, default=1.0)
    parser.add_argument("--decades", type=float, default=12.0)
    parser.add_argument("--select", default="",
                        help="comma list of MODE:THREADS:BLOCKS:WIDTH to run")
    parser.add_argument("--note", action="append", default=[],
                        help="NAME=VALUE line for the manifest")
    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)
    select = set(filter(None, args.select.split(",")))
    rows = []
    for blocks, width, repeats in SHAPES:
        for mode in args.modes.split(","):
            for threads in THREADS[mode]:
                if threads > len(os.sched_getaffinity(0)):
                    continue
                if select and f"{mode}:{threads}:{blocks}:{width}" not in select:
                    continue
                for trial in range(1, args.trials + 1):
                    row = run_one(args.binary, mode, threads, blocks, width,
                                  repeats, args.decades, args.timeout)
                    row["trial"] = trial
                    rows.append(row)
                    if row["status"] == "FAILED":
                        break  # not built, or a deterministic failure
    fields = []
    for row in rows:
        fields += [name for name in row if name not in fields]
    with open(os.path.join(args.output, "raw.csv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, restval="",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    write_summary(rows, os.path.join(args.output, "summary.csv"))
    with open(args.binary, "rb") as handle:
        digest = hashlib.sha256(handle.read()).hexdigest()
    with open(os.path.join(args.output, "manifest.txt"), "w") as handle:
        handle.write("\n".join(args.note + [
            f"command={' '.join(sys.argv)}",
            f"binary_sha256={digest}",
            f"trials={args.trials}",
            f"timeout_s={args.timeout}",
            f"allowed_cpus={' '.join(map(str, sorted(os.sched_getaffinity(0))))}",
            f"machine={platform.machine()} {platform.release()}",
            "environment=OMP_NUM_THREADS=THREADS OPENBLAS_NUM_THREADS=THREADS "
            "OMP_DYNAMIC=false OMP_PLACES=cores OMP_PROC_BIND=close",
            "[lscpu]", text_of(["lscpu"]),
            "[ldd]", text_of(["ldd", args.binary]), ""]))


if __name__ == "__main__":
    main()
