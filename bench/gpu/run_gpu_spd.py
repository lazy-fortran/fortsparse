#!/usr/bin/env python3
"""Optional matched GPU benchmark orchestration; never part of library runtime."""
import argparse
from collections import Counter
from copy import copy
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import random
import signal
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_spd import build_provenance, physical_cpus
from accounting_peak import AccountingPeak

EXECUTABLES = {"cudss": "bench_spd_cudss", "cholmod": "bench_spd_cholmod",
               "native": "bench_spd_native_gpu", "native-cuda": "bench_spd_native_gpu",
               "native-nd": "bench_spd_native_gpu", "native-cuda-nd": "bench_spd_native_gpu"}
STRINGS = {"backend", "family", "workload", "variant", "io_mode", "ordering",
           "ordering_hash", "solver_version", "execution_mode", "device_total_peak_source"}


def parse_metrics(stdout):
    lines = stdout.splitlines()
    first = next((i for i, line in enumerate(lines) if line.startswith("backend,")), None)
    if first is None:
        return {}
    row = next(csv.DictReader(io.StringIO("\n".join(lines[first:]))))
    result = {}
    for key, value in row.items():
        if key == "reused":
            logical = value.strip().upper().strip('.')
            if logical not in ("T", "F", "TRUE", "FALSE", "0", "1"):
                raise ValueError("invalid Fortran logical reused")
            result[key] = logical in ("T", "TRUE", "1")
        elif key in STRINGS:
            result[key] = value.strip()
        else:
            number = float(value)
            if not math.isfinite(number):
                raise ValueError("nonfinite benchmark metric: " + key)
            result[key] = number
    return result


def phase_objective(row):
    workload = row["workload"]
    terms = (("refactor_s", "solve_s") if workload == "refactor" else
             ("solve_s",) if workload == "repeated-rhs" else ("factor_s", "solve_s"))
    values = [row.get(key, -1) for key in terms]
    return sum(values) if all(v is not None and v >= 0 for v in values) else None


def unavailable(backend, io_mode, ordering):
    if io_mode == "resident" and backend != "cudss":
        return "device-resident solve is unsupported by this measured adapter"
    if ordering == "common" and backend.startswith("native"):
        return "native public API does not accept a supplied permutation"
    return None


def observe_accounting(args, record, pid):
    observer = getattr(args, "accounting_observer", None)
    if observer is not None:
        record["process_pid"] = pid
        record["vendor_device_total_peak_bytes"] = record.get("device_total_peak_bytes", -1)
        record.update(observer.query(pid))
    return record


def sample(args, cpus, backend, family, level, thread, workload, io_mode, ordering, trial):
    if backend in ("native-nd", "native-cuda-nd"):
        ordering = "nd"
    device = -1 if backend in ("native", "native-nd") else args.device
    repeats = args.repeats if workload == "repeated-rhs" else 1
    binary_dir = (getattr(args, "cpu_bin_dir", None) or args.bin_dir) if device < 0 else args.bin_dir
    command = ["taskset", "--cpu-list", ",".join(map(str, cpus[:thread])),
        str(binary_dir / EXECUTABLES[backend]), family, str(level), str(thread),
        workload, args.variant, str(device), io_mode, ordering,
        str(args.matrix) if args.matrix else "-", str(repeats),
        str(args.permutation) if args.permutation else "-"]
    record = dict(backend=backend, requested_backend=backend, family=family, level=level,
                  threads=thread, workload=workload, variant=args.variant, trial=trial,
                  requested_device=device, io_mode=io_mode, ordering=ordering,
                  repeats=repeats, command=command)
    record["pool_budget_bytes"] = int(os.environ.get(
        "CHOLMOD_GPU_POOL_BYTES" if backend == "cholmod" else
        "FSBENCH_CUDSS_POOL_LIMIT_BYTES", "268435456")) if backend in (
            "cholmod", "cudss") else 0
    reason = unavailable(backend, io_mode, ordering)
    if reason or not Path(command[3]).is_file():
        return dict(record, status="unavailable", reason=reason or "executable missing")
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS=str(thread), OMP_DYNAMIC="FALSE",
                       OMP_PROC_BIND="close", OMP_PLACES="cores", OPENBLAS_NUM_THREADS="1",
                       MKL_NUM_THREADS="1")
    start = time.monotonic()
    process = subprocess.Popen(command, env=environment, text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, start_new_session=True)
    try:
        stdout, stderr = process.communicate(timeout=1.0)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        stdout, stderr = process.communicate()
        record.update(status="timeout", process_wall_s=time.monotonic() - start,
                      stdout=stdout, stderr=stderr)
        try:
            record.update(parse_metrics(stdout))
        except (TypeError, ValueError, StopIteration):
            pass
        return observe_accounting(args, record, process.pid)
    record.update(process_wall_s=time.monotonic() - start, returncode=process.returncode)
    try:
        record.update(parse_metrics(stdout))
    except (TypeError, ValueError, StopIteration) as error:
        return observe_accounting(args, dict(record, status="invalid_output", reason=str(error),
                                  stdout=stdout, stderr=stderr), process.pid)
    if process.returncode or record.get("vendor_status", -1) != 0:
        record.update(status="failed", stdout=stdout, stderr=stderr)
    elif record.get("backward_error", -1) < 0 or record["backward_error"] > 1.e-12:
        record.update(status="accuracy_failed", stdout=stdout, stderr=stderr)
    else:
        record["status"] = "ok"
        record["primary_phase_s"] = phase_objective(record)
        record["workflow_cpu_s"] = record["parent_cpu_s"] + record.get("helper_cpu_s", 0.)
    return observe_accounting(args, record, process.pid)


def sample_block(args, cpus, backends, family, level, thread, workload,
                 io_mode, ordering, trial, randomizer):
    """Randomize candidate/reference processes together inside each trial block."""
    block = [(args, backend, "candidate") for backend in backends]
    reference_dir = getattr(args, "reference_bin_dir", None)
    if reference_dir:
        reference = copy(args)
        reference.bin_dir = reference_dir
        reference.cpu_bin_dir = reference_dir
        requested = getattr(args, "reference_backends", "native-cuda,native-cuda-nd").split(",")
        block += [(reference, backend, "reference") for backend in backends if backend in requested]
    randomizer.shuffle(block)
    return [dict(sample(configuration, cpus, backend, family, level, thread,
                        workload, io_mode, ordering, trial), series=series)
            for configuration, backend, series in block]


def build_sources(args):
    result = dict(build_provenance=build_provenance(args.bin_dir),
                  cpu_build_provenance=build_provenance(args.cpu_bin_dir)
                      if args.cpu_bin_dir else None,
                  reference_build_provenance=None)
    if args.reference_bin_dir:
        reference = build_provenance(args.reference_bin_dir)
        reference["source_hash_scope"] = (
            "Observed current source at the archived cache's CMAKE_HOME_DIRECTORY; "
            "not an attestation of the historical reference build. Executable SHA is authoritative.")
        result["reference_build_provenance"] = reference
    return result


def snapshot():
    result = subprocess.run(["nvidia-smi", "--query-gpu=index,uuid,pci.bus_id,name,driver_version,memory.total,memory.used",
                             "--format=csv,noheader,nounits"], text=True, capture_output=True)
    return dict(utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                nvidia_smi=result.stdout.strip(), status=result.returncode,
                loadavg=Path("/proc/loadavg").read_text().strip())


def file_identity(path):
    return dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--cpu-bin-dir", type=Path,
                        help="Use a CPU-only native build for CPU comparisons")
    parser.add_argument("--reference-bin-dir", type=Path,
                        help="Interleave frozen reference executables with candidates")
    parser.add_argument("--reference-backends", default="native-cuda,native-cuda-nd",
                        help="Reference copies of these active --backends (same backend names)")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--backends", default="native-cuda,cudss,cholmod")
    parser.add_argument("--cases", default="band:1,band:2,band:3,grid2:1,grid2:2,grid2:3,grid2:7,grid3:1,grid3:2,grid3:3,grid3:7")
    parser.add_argument("--workloads", default="rhs1,rhs8,warm,refactor,repeated-rhs")
    parser.add_argument("--io-modes", default="host")
    parser.add_argument("--orderings", default="tuned")
    parser.add_argument("--threads", default="1,2,4,8,16")
    parser.add_argument("--trials", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=8)
    parser.add_argument("--variant", choices=("training", "heldout"), default="training")
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--matrix", type=Path)
    parser.add_argument("--matrix-family", default="matrix-market")
    parser.add_argument("--permutation", type=Path)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--accounting-peak", action="store_true",
                        help="Read own terminated child PID's already-enabled NVML accounting peak")
    parser.add_argument("--calibration-only", action="store_true")
    parser.add_argument("--measurement-class", choices=("diagnostic-shared", "controlled"),
                        default="diagnostic-shared")
    parser.add_argument("--resource-note", default="No exclusive device lease asserted")
    args = parser.parse_args()
    args.bin_dir = args.bin_dir.resolve()
    if args.reference_bin_dir:
        args.reference_bin_dir = args.reference_bin_dir.resolve()
    backends = args.backends.split(",")
    if any(b not in EXECUTABLES for b in args.reference_backends.split(",")):
        parser.error("unknown reference backend")
    if any(b not in EXECUTABLES for b in backends):
        parser.error("unknown backend")
    if args.calibration_only:
        backends = [b for b in backends if not b.startswith("native")]
        if not backends:
            parser.error("external-only calibration requires an external backend")
    orderings = args.orderings.split(",")
    if "common" in orderings and not args.permutation:
        parser.error("common ordering requires the retained permutation file")
    if args.permutation and not args.permutation.is_file():
        parser.error("permutation file missing")
    cpus = physical_cpus()
    threads = [int(t) for t in args.threads.split(",")]
    if min(threads) < 1 or max(threads) > len(cpus) or args.trials < 1 or args.repeats < 1:
        parser.error("positive available thread/trial/repeat counts required")
    cases = [(f, int(n)) for f, n in (entry.split(":") for entry in args.cases.split(","))]
    if args.matrix:
        if not args.matrix.is_file() or args.variant == "heldout":
            parser.error("matrix file required; archived matrices have no held-out coefficients")
        args.matrix = args.matrix.resolve()
        cases = [(args.matrix_family, 1)]
    args.accounting_observer = AccountingPeak() if args.accounting_peak else None
    before = snapshot()
    records = []
    randomizer = random.Random(args.seed)
    try:
        for family, level in cases:
            for thread in threads:
                for workload in args.workloads.split(","):
                    for io_mode in args.io_modes.split(","):
                        for ordering in orderings:
                            for trial in range(args.trials):
                                records.extend(sample_block(args, cpus, backends, family, level,
                                    thread, workload, io_mode, ordering, trial, randomizer))
    finally:
        if args.accounting_observer is not None:
            args.accounting_observer.close()
    binaries = {}
    for binary in sorted({Path(r["command"][3]) for r in records}):
        if binary.is_file():
            binaries[str(binary)] = dict(file_identity(binary), ldd=subprocess.run(
                ["ldd", str(binary)], text=True, capture_output=True).stdout)
    for record in records:
        record["executable_sha256"] = binaries.get(record["command"][3], {}).get("sha256")
    result = dict(schema="fortsparse-gpu-benchmark-v1", trials=records,
        command=sys.argv, revision=subprocess.run(["git", "rev-parse", "HEAD"],
            text=True, capture_output=True).stdout.strip(),
        working_tree_dirty=bool(subprocess.run(["git", "status", "--porcelain"],
            text=True, capture_output=True).stdout.strip()), platform=platform.platform(),
        source_files=[file_identity(Path(__file__)), file_identity(Path(__file__).parents[1]/"spd_fixtures.f90")],
        binary_provenance=binaries, **build_sources(args),
        reference_bin_dir=str(args.reference_bin_dir) if args.reference_bin_dir else None,
        reference_backends=args.reference_backends.split(",") if args.reference_bin_dir else [],
        physical_cpu_ids=cpus, gpu_observation=dict(before=before, after=snapshot()),
        measurement_class=args.measurement_class, resource_note=args.resource_note,
        calibration_only=args.calibration_only, corpus=cases, seed=args.seed,
        accounting_observer=args.accounting_observer.metadata if args.accounting_observer else None,
        matrix_provenance=file_identity(args.matrix) if args.matrix else None,
        permutation_provenance=file_identity(args.permutation) if args.permutation else None,
        status_counts=dict(Counter(r["status"] for r in records)),
        timing_contract="primary_phase_s: cold/warm=factor+solve; refactor=refactor+solve; repeated-rhs=TOTAL retained-factor solves. total_s and workflow_cpu_s include setup/check/free and warmups.",
        memory_contract="Requested/used/reserved pools, wrapper bytes and phase-sampled process memory are separate. Optional device_total_peak_source=nvml_accounting is the driver's lifetime process maxMemoryUsage for one uniquely matched GPU; original vendor peak and parent observer cost remain separate. Unknown values never use sampled substitution.",
        environment={key: os.environ.get(key) for key in
            ("CUDA_VISIBLE_DEVICES", "CUDA_MODULE_LOADING", "OMP_WAIT_POLICY",
             "FSBENCH_CUDSS_POOL_LIMIT_BYTES", "CHOLMOD_GPU_POOL_BYTES")})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"{len(records)} trials: {result['status_counts']}; {args.output}")


if __name__ == "__main__":
    main()
