#!/usr/bin/env python3
"""Bounded fresh-process family measurements; raw attempts are never discarded."""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import signal
import statistics
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXTERNAL = ("superlu", "umfpack", "mumps", "intel-pardiso")
KEYS = ("backend", "family", "level", "threads", "workload", "variant", "operation")


def physical_cores():
    chosen = {}
    for cpu in sorted(os.sched_getaffinity(0)):
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        package = (topology / "physical_package_id").read_text().strip()
        core = (topology / "core_id").read_text().strip()
        chosen.setdefault((package, core), cpu)
    return list(chosen.values())


def parse_csv(stdout):
    lines = stdout.splitlines()
    for index, line in enumerate(lines):
        if not line.startswith("backend,family,level,threads,"):
            continue
        try:
            record = next(csv.DictReader(io.StringIO("\n".join(lines[index:index+2]))))
            for key, value in tuple(record.items()):
                if key in KEYS or key in ("status", "native_factor_family", "reused"):
                    continue
                record[key] = float(value)
            record["level"] = int(record["level"])
            record["threads"] = int(record["threads"])
            return record
        except (StopIteration, ValueError, TypeError):
            return None  # Partial/crashed output remains in the raw receipt.
    return None


def invoke(binary, args, threads, cores, timeout):
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS="1",
               OMP_PLACES="cores", OMP_PROC_BIND="close", MKL_NUM_THREADS=str(threads),
               MKL_DYNAMIC="FALSE", MKL_THREADING_LAYER="GNU")
    background_before = os.getloadavg()
    start = time.perf_counter()
    process = subprocess.Popen([str(binary), *args], stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, text=True, env=env,
                               start_new_session=True,
                               preexec_fn=lambda: os.sched_setaffinity(0, cores))
    timed_out = False
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        os.killpg(process.pid, signal.SIGKILL)
        stdout, stderr = process.communicate()
    receipt = dict(command=[str(binary), *args], affinity=cores, environment={
        key: env[key] for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                 "OMP_PLACES", "OMP_PROC_BIND", "MKL_NUM_THREADS",
                                 "MKL_DYNAMIC", "MKL_THREADING_LAYER")},
        process_wall_s=time.perf_counter()-start, returncode=process.returncode,
        background_load_before=background_before, background_load_after=os.getloadavg(),
        stdout=stdout, stderr=stderr, status="timeout" if timed_out else "failed")
    record = parse_csv(stdout)
    if record and not timed_out:
        record["process_wall_s"] = receipt["process_wall_s"]
        record["workflow_cpu_s"] = record["parent_cpu_s"]+record["helper_cpu_s"]
        receipt["record"] = record
        receipt["status"] = record["status"]
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mode", choices=("measure", "calibrate"), default="measure")
    parser.add_argument("--backends", default="native,"+",".join(EXTERNAL))
    parser.add_argument("--cases", default="")
    parser.add_argument("--threads", default="1,2,4,8,16")
    parser.add_argument("--trials", type=int, default=5)
    parser.add_argument("--timeout", type=float, default=0.9)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--levels-from", type=Path)
    parser.add_argument("--workloads", default="")
    parser.add_argument("--operations", default="N,T,H")
    parser.add_argument("--variants", default="training,heldout")
    parser.add_argument("--window-note", default="uncoordinated diagnostic window")
    args = parser.parse_args()
    if not 0 < args.timeout < 1 or args.trials < 1:
        parser.error("positive trials and timeout strictly below one second required")
    manifest = json.loads((HERE / "contract.json").read_text())
    if hashlib.sha256((HERE / "fixtures.f90").read_bytes()).hexdigest() != \
            manifest["fixture_source_sha256"]:
        parser.error("fixture source differs from the frozen contract")
    if hashlib.sha256((ROOT/"bench/fixtures/neo2-mini5/matrix.mtx").read_bytes()).hexdigest() != \
            manifest["mini5_matrix_sha256"]:
        parser.error("retained NEO-2 numeric fixture differs from frozen provenance")
    backends = args.backends.split(",")
    if any(b not in ("native", *EXTERNAL) for b in backends):
        parser.error("unknown backend")
    if args.mode == "calibrate":
        backends = [b for b in backends if b != "native"]
        if not backends or any(b not in EXTERNAL for b in backends):
            parser.error("calibration uses external comparators only")
    else:
        if args.levels_from:
            selected = json.loads(args.levels_from.read_text())
            if selected.get("mode") != "external_only" or selected.get(
                    "fixture_source_sha256") != manifest["fixture_source_sha256"]:
                parser.error("selection is not external-only calibration of these fixtures")
            manifest["performance_levels"] = selected["levels"]
        if manifest["performance_levels"] is None:
            parser.error("performance levels await external-only calibration; use --mode calibrate")
    cases = args.cases.split(",") if args.cases else list(manifest["cases"])
    threads = [int(v) for v in args.threads.split(",")]
    if args.mode == "measure" and any(
            manifest["performance_levels"].get(case) is None for case in cases):
        parser.error("requested case has no complete external-only size selection")
    if any(nt < 1 for nt in threads):
        parser.error("positive thread counts required")
    cores = physical_cores()
    if (args.output/"attempts.jsonl").exists():
        parser.error("output already contains attempts; use a new directory")
    args.output.mkdir(parents=True, exist_ok=True)
    provenance = dict(source_revision=subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        source_status=subprocess.check_output(
            ["git", "-C", str(ROOT), "status", "--porcelain"], text=True),
        contract=manifest, mode=args.mode, trials=args.trials, physical_cores=cores,
        timeout_s=args.timeout, window_note=args.window_note,
        background_load=os.getloadavg(), binaries={})
    for backend in backends:
        path = args.build.resolve() / "bin" / ("bench_native_families_pardiso" if
            backend == "intel-pardiso" else "bench_native_families")
        if path.exists():
            provenance["binaries"][backend] = dict(path=str(path), sha256=hashlib.sha256(
                path.read_bytes()).hexdigest(), ldd=subprocess.run(
                ["ldd", str(path)], capture_output=True, text=True).stdout)
    helper = args.build.resolve()/"bin/fortsparse_umfpack_helper"
    if helper.exists():
        provenance["helper"] = dict(path=str(helper), sha256=hashlib.sha256(
            helper.read_bytes()).hexdigest(), ldd=subprocess.run(
            ["ldd", str(helper)], capture_output=True, text=True).stdout)
    (args.output / "provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
    receipts = []
    raw_path = args.output / "attempts.jsonl"
    with raw_path.open("w") as output:
        for family in cases:
            level = manifest["initial_calibration_levels"][family] if args.mode == "calibrate" else \
                manifest["performance_levels"][family]
            for nt in threads:
                workloads = ["rhs8"] if args.mode == "calibrate" else (
                    args.workloads.split(",") if args.workloads else manifest["workloads"])
                variants = ["training"] if args.mode == "calibrate" else args.variants.split(",")
                operations = ["N"] if args.mode == "calibrate" else args.operations.split(",")
                for workload in workloads:
                    for variant in variants:
                        for operation in operations:
                            for trial in range(args.trials):
                                # Interleave peers at each matched point, rotating first position.
                                offset = trial % len(backends)
                                peers = backends[offset:]+backends[:offset]
                                for backend in peers:
                                    binary = args.build.resolve()/"bin"/("bench_native_families_pardiso"
                                        if backend == "intel-pardiso" else "bench_native_families")
                                    command = [family, str(level), str(nt), workload, operation,
                                               variant, str(ROOT/"bench/fixtures/neo2-mini5/matrix.mtx")]
                                    if backend != "intel-pardiso":
                                        command.insert(0, backend)
                                    if args.dry_run:
                                        print(str(binary), *command)
                                        continue
                                    if nt > len(cores) or not binary.exists():
                                        receipt = dict(status="unavailable", backend=backend,
                                            family=family, level=level, threads=nt, trial=trial,
                                            reason="physical cores or comparator binary unavailable")
                                    else:
                                        receipt = invoke(binary, command, nt, cores[:nt], args.timeout)
                                        receipt["trial"] = trial
                                    output.write(json.dumps(receipt)+"\n")
                                    output.flush()
                                    receipts.append(receipt)
    groups = {}
    for attempt in receipts:
        row = attempt.get("record")
        if row and attempt["status"] == "ok":
            groups.setdefault(tuple(row[k] for k in KEYS), []).append(row)
    summaries = []
    for key, rows in groups.items():
        summary = dict(zip(KEYS, key), successes=len(rows))
        for metric, value in rows[0].items():
            if metric in KEYS or not isinstance(value, (float, int)):
                continue
            values = [r[metric] for r in rows if r[metric] >= 0]
            if values:
                summary[metric] = dict(median=statistics.median(values),
                                      minimum=min(values), maximum=max(values))
        summaries.append(summary)
    if args.mode == "calibrate" and not args.dry_run:
        levels = {}
        for family in cases:
            eligible = []
            for level in [manifest["initial_calibration_levels"][family]]:
                sample = [r for r in receipts if r.get("record", {}).get("family") == family
                          and r.get("record", {}).get("level") == level]
                expected = len(backends)*len(threads)*args.trials
                if len(sample) == expected and all(r["status"] == "ok" and
                        r["process_wall_s"] < 0.5 for r in sample):
                    eligible.append(level)
            levels[family] = max(eligible) if eligible else None
        selection = dict(mode="external_only", fixture_source_sha256=manifest[
            "fixture_source_sha256"], levels=levels, backends=backends,
            threads=threads, trials=args.trials, attempts_path=str(raw_path.resolve()),
            attempts_sha256=hashlib.sha256(raw_path.read_bytes()).hexdigest(),
            admission="unknown; calibration establishes size only")
        (args.output / "calibration_selection.json").write_text(
            json.dumps(selection, indent=2)+"\n")
    (args.output / "summary.json").write_text(json.dumps(dict(
        attempts=len(receipts), successes=sum(r["status"] == "ok" for r in receipts),
        admission="unknown; no Pareto inference", groups=summaries), indent=2)+"\n")


if __name__ == "__main__":
    main()
