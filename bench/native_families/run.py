#!/usr/bin/env python3
"""Bounded fresh-process family measurements; raw attempts are never discarded."""
import argparse
import csv
import hashlib
import io
import json
import itertools
import math
import os
from pathlib import Path
import signal
import statistics
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXTERNAL = ("superlu", "umfpack", "mumps", "intel-pardiso")
NATIVE = ("native", "native-amd")
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
        record["workflow_cpu_s"] = (record["parent_cpu_s"]+record["helper_cpu_s"]
            if min(record["parent_cpu_s"], record["helper_cpu_s"]) >= 0 else -1.0)
        receipt["record"] = record
        receipt["status"] = record["status"]
        if process.returncode != 0 and record["status"] == "ok":
            receipt["status"] = "failed"
    return receipt



def attempt_passes(attempt, limit_s=0.5):
    """A completed successful external point below the strict wall-time gate."""
    try:
        wall = float(attempt["process_wall_s"])
    except (KeyError, TypeError, ValueError):
        return False
    return (attempt.get("status") == "ok" and attempt.get("returncode", 0) == 0
            and math.isfinite(wall) and 0 <= wall < limit_s)


def select_calibration_level(receipts, family, candidates, backends, threads,
                             trials, limit_s=0.5):
    """Largest complete passing prefix; never extrapolate past a failed level."""
    if not backends or any(backend not in EXTERNAL for backend in backends):
        raise ValueError("size selection requires external comparators only")
    if trials < 1 or not threads or any(nt < 1 for nt in threads):
        raise ValueError("positive trials and thread counts required")
    if len(set(backends)) != len(backends) or len(set(threads)) != len(threads):
        raise ValueError("duplicate calibration configurations")
    if sorted(set(candidates)) != list(candidates) or not candidates:
        raise ValueError("candidate levels must be nonempty and strictly increasing")
    expected = set(itertools.product(backends, threads, range(trials)))
    selected, outcomes, failed = None, [], False
    for level in candidates:
        if failed:
            outcomes.append(dict(level=level, status="unmeasured_after_failure"))
            continue
        sample = [row for row in receipts if row.get("family") == family
                  and row.get("level") == level and row.get("backend") in backends
                  and row.get("workload") == "rhs8" and row.get("operation") == "N"
                  and row.get("variant") == "training"]
        keys = [(row.get("backend"), row.get("threads"), row.get("trial"))
                for row in sample]
        complete = set(keys) == expected and len(keys) == len(expected)
        passed = complete and all(attempt_passes(row, limit_s) for row in sample)
        outcome = "pass" if passed else "incomplete_configurations"
        if not passed:
            statuses = {row.get("status") for row in sample}
            if statuses.intersection({"unavailable", "backend_unavailable"}):
                outcome = "comparator_unavailable"
            elif "timeout" in statuses:
                outcome = "timeout"
            elif any(row.get("status") != "ok" or row.get("returncode", 0) != 0
                     for row in sample):
                outcome = "solver_or_accuracy_failure"
            elif any(not attempt_passes(row, limit_s) for row in sample):
                outcome = "wall_gate_exceeded"
        outcomes.append(dict(level=level, status=outcome,
            attempts=len(sample), expected_attempts=len(expected), complete=complete,
            rejected_attempts=sum(not attempt_passes(row, limit_s) for row in sample)))
        if passed:
            selected = level
        else:
            failed = True
    return selected, outcomes



def validate_selected_levels(selection, manifest, families, receipts,
                             diagnostic=False):
    """Recheck claimed levels against hashed raw external-only observations."""
    if selection.get("mode") != "external_only" or selection.get(
            "fixture_source_sha256") != manifest["fixture_source_sha256"]:
        raise ValueError("selection does not describe these frozen external-only fixtures")
    peers, threads, trials = (selection.get("backends", []),
                              selection.get("threads", []), selection.get("trials", 0))
    if not peers or any(peer not in EXTERNAL for peer in peers):
        raise ValueError("native or unknown comparator in size selection")
    complete = (set(peers) == set(manifest["required_external_backends"])
                and set(threads) == set(manifest["thread_sweep"])
                and trials >= manifest["fresh_process_trials"])
    if not complete and not diagnostic:
        raise ValueError("selection omits declared external/thread/trial baselines; "
                         "only explicit --diagnostic use is allowed")
    for family in families:
        if family not in manifest["cases"]:
            raise ValueError("unknown selected family")
        expected, _ = select_calibration_level(receipts, family,
            manifest["cases"][family]["calibration_candidates"], peers, threads, trials)
        if expected is None or selection.get("levels", {}).get(family) != expected:
            raise ValueError("selected family level lacks complete passing raw evidence")
    return "complete_declared_baselines" if complete else "diagnostic_subset"


def matched_points(backends, threads, trials, workloads, variants, operations):
    for nt, workload, variant, operation, trial in itertools.product(
            threads, workloads, variants, operations, range(trials)):
        offset = trial % len(backends)
        for backend in backends[offset:]+backends[:offset]:
            yield backend, nt, workload, variant, operation, trial


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
    parser.add_argument("--diagnostic", action="store_true",
                        help="allow explicitly incomplete external-only selection evidence")
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
    if any(b not in (*NATIVE, *EXTERNAL) for b in backends):
        parser.error("unknown backend")
    selected = None
    selection_scope = "not_applicable"
    if args.mode == "calibrate":
        backends = [b for b in backends if b not in NATIVE]
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
    if any(nt < 1 for nt in threads) or len(set(threads)) != len(threads):
        parser.error("distinct positive thread counts required")
    if len(set(backends)) != len(backends):
        parser.error("distinct backends required")
    if any(case not in manifest["cases"] for case in cases):
        parser.error("unknown case")
    if selected is not None:
        raw = Path(selected["attempts_path"])
        if not raw.is_absolute():
            raw = args.levels_from.resolve().parent/raw
        raw_bytes = raw.read_bytes()
        if hashlib.sha256(raw_bytes).hexdigest() != selected.get("attempts_sha256"):
            parser.error("calibration raw receipt hash differs from the selection")
        raw_receipts = [json.loads(line) for line in raw_bytes.decode().splitlines() if line]
        try:
            selection_scope = validate_selected_levels(selected, manifest, cases,
                raw_receipts, diagnostic=args.diagnostic)
        except ValueError as error:
            parser.error(str(error))
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
        cpu_measurement_scope="parent_cumulative_rusage_self; helper_reaped_rusage_children",
        selection_scope=selection_scope, diagnostic=args.diagnostic,
        requested_backends=backends, requested_cases=cases, requested_threads=threads,
        requested_workloads=args.workloads or manifest["workloads"],
        requested_operations=args.operations, requested_variants=args.variants,
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
            levels = manifest["cases"][family]["calibration_candidates"] if args.mode == "calibrate" else \
                [manifest["performance_levels"][family]]
            for level in levels:
                workloads = ["rhs8"] if args.mode == "calibrate" else (
                    args.workloads.split(",") if args.workloads else manifest["workloads"])
                variants = ["training"] if args.mode == "calibrate" else args.variants.split(",")
                operations = ["N"] if args.mode == "calibrate" else args.operations.split(",")
                failed_candidate = False
                for backend, nt, workload, variant, operation, trial in matched_points(
                        backends, threads, args.trials, workloads, variants, operations):
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
                        receipt = dict(status="unavailable",
                            reason="physical cores or comparator binary unavailable")
                    else:
                        receipt = invoke(binary, command, nt, cores[:nt], args.timeout)
                    receipt.update(backend=backend, family=family, level=level, threads=nt,
                        workload=workload, variant=variant, operation=operation, trial=trial)
                    output.write(json.dumps(receipt)+"\n")
                    output.flush()
                    receipts.append(receipt)
                    if args.mode == "calibrate" and not attempt_passes(receipt):
                        failed_candidate = True
                        break  # The remaining configurations cannot make this level pass.
                if failed_candidate:
                    break  # Larger levels remain explicitly unmeasured for this family.
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
        levels, candidate_outcomes = {}, {}
        for family in cases:
            levels[family], candidate_outcomes[family] = select_calibration_level(
                receipts, family, manifest["cases"][family]["calibration_candidates"],
                backends, threads, args.trials)
        selection = dict(mode="external_only", fixture_source_sha256=manifest[
            "fixture_source_sha256"], levels=levels, backends=backends,
            threads=threads, trials=args.trials, candidate_outcomes=candidate_outcomes,
            declared_external_backends=manifest["required_external_backends"],
            calibration_scope=("complete_declared_baselines" if
                set(backends) == set(manifest["required_external_backends"]) and
                set(threads) == set(manifest["thread_sweep"]) and
                args.trials >= manifest["fresh_process_trials"] else "diagnostic_subset"),
            scan_scope="largest passing scanned prefix; larger levels after failure unmeasured",
            attempts_path=str(raw_path.resolve()),
            attempts_sha256=hashlib.sha256(raw_path.read_bytes()).hexdigest(),
            admission="unknown; calibration establishes size only")
        (args.output / "calibration_selection.json").write_text(
            json.dumps(selection, indent=2)+"\n")
    (args.output / "summary.json").write_text(json.dumps(dict(
        attempts=len(receipts), successes=sum(r["status"] == "ok" for r in receipts),
        admission="unknown; no Pareto inference", groups=summaries), indent=2)+"\n")


if __name__ == "__main__":
    main()
