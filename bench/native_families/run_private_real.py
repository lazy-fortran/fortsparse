#!/usr/bin/env python3
"""Explicit qualification and matched measurements of a frozen private fixture."""
import argparse
import collections
import json
import os
from pathlib import Path
import subprocess
import time

from real_fixture_audit import PrivateFixture, sha256
from run import (EXTERNAL, NATIVE, ROOT, binary_name, invoke, physical_cores, matched_points,
                 positive_i8, memory_environment, record_memory_configuration)

HERE = Path(__file__).resolve().parent


def binary_path(build, backend):
    return (build / "bin" / binary_name(backend)).resolve()


def prepare(fixture, output):
    output = output.resolve()
    if not output.is_relative_to(fixture.base):
        raise ValueError("numeric output must stay under the private fixture archive")
    parent = output.parent
    while not parent.exists():
        parent = parent.parent
    check = subprocess.run(["git", "-C", str(parent), "rev-parse",
                            "--is-inside-work-tree"], capture_output=True, text=True)
    if check.returncode == 0 and check.stdout.strip() == "true":
        raise ValueError("numeric output in a Git working tree is forbidden")
    output.mkdir(mode=0o700, parents=True, exist_ok=False)
    output.chmod(0o700)
    stage = output / "input"
    stage.mkdir(mode=0o700)
    paths = {"matrix.mtx": fixture.matrix,
             "physical-rhs.mtx": fixture.physical_rhs,
             "known-rhs.mtx": fixture.checked_file(fixture.reader["manufactured_rhs"]),
             "known-x.mtx": fixture.checked_file(fixture.reader["manufactured_known_x"]),
             "physical-reference.mtx": fixture.checked_file(
                 fixture.manifest["references"]["umfpack"]["physical"])}
    for name, path in paths.items():
        (stage / name).symlink_to(path)
    return stage


def verify_qualification(paths, fixture, build):
    """Reject stale, partial or merely asserted success before native execution."""
    if isinstance(paths, Path):
        paths = [paths]
    attempts, hashes = [], {}
    for path in paths:
        document = json.loads(path.read_text())
        if (document.get("mode") != "qualify" or document.get("status") != "pass"
                or document.get("manifest_sha256") != fixture.manifest_sha256
                or document.get("policy_sha256") != fixture.policy_sha256):
            raise ValueError("passing qualification of this exact fixture/policy required")
        attempts.extend(document.get("attempts", []))
        hashes[str(path.resolve())] = sha256(path)
    if len(attempts) != len(EXTERNAL) or {a["backend"] for a in attempts} != set(EXTERNAL):
        raise ValueError("all five declared external comparators must qualify")
    for entry in attempts:
        binary = binary_path(build, entry["backend"])
        if not binary.is_file() or sha256(binary) != entry.get("binary_sha256"):
            raise ValueError("external comparator binary changed after qualification")
        candidate = Path(entry["candidate_path"]).resolve()
        if not candidate.is_relative_to(fixture.base):
            raise ValueError("qualification candidate escapes private archive")
        if sha256(candidate) != entry.get("candidate_sha256"):
            raise ValueError("external qualification candidate changed")
        if entry.get("status") != "pass" or fixture.check(candidate)["status"] != "pass":
            raise ValueError("independent audit rejects external qualification")
    return hashes


def attempt(backend, build, stage, output, fixture, threads, cores, timeout, host_limit_bytes=None):
    binary = binary_path(build, backend)
    result = dict(backend=backend, binary=str(binary), status="unavailable")
    if not binary.is_file():
        result["reason"] = "comparator executable not built"
        return result
    result["binary_sha256"] = sha256(binary)
    deps = subprocess.run(["ldd", str(binary)], capture_output=True, text=True)
    result["dependency_resolution"] = dict(stdout=deps.stdout, stderr=deps.stderr,
                                            returncode=deps.returncode)
    candidate = output / (backend + "-candidate.mtx")
    args = ["real-mm", "1", str(threads), "rhs4", "N", "training", str(stage)]
    if backend not in ("intel-pardiso", "umfpack-direct"):
        args.insert(0, backend)
    receipt = invoke(binary, args, threads, cores, timeout,
                     extra_env={"FORTSPARSE_BENCH_SOLUTION": str(candidate),
                                **memory_environment(backend,host_limit_bytes)})
    record_memory_configuration(receipt,backend,host_limit_bytes)
    result["process_receipt"] = receipt
    if backend == "umfpack":
        result["aggregate_memory"] = "unknown; parent/helper RSS must not be summed"
    result["status"] = receipt["status"]
    if candidate.is_file():
        candidate.chmod(0o600)
        result.update(candidate_path=str(candidate), candidate_sha256=sha256(candidate))
        start = time.perf_counter()
        try:
            result["independent_audit"] = fixture.check(candidate)
        except (ValueError, OSError, StopIteration) as error:
            result["independent_audit"] = dict(status="audit_failure", reason=str(error))
        duration = time.perf_counter() - start
        result["independent_audit_including_reader_wall_s"] = duration
        result["solver_child_plus_audit_wall_s"] = receipt["process_wall_s"] + duration
        if receipt["status"] == "ok" and receipt["returncode"] == 0:
            result["status"] = result["independent_audit"]["status"]
    elif receipt["status"] == "ok":
        result["status"] = "missing_candidate_export"
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", required=True, type=Path)
    parser.add_argument("--policy-sha256", required=True)
    parser.add_argument("--build", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mode", choices=("prepare", "qualify", "native", "measure"),
                        default="prepare")
    parser.add_argument("--backend", choices=NATIVE)
    parser.add_argument("--native-host-limit-bytes", type=positive_i8,
                        help="required positive i8 cap for native-amd-capped")
    parser.add_argument("--qualification", type=Path, nargs="+")
    parser.add_argument("--backends", default=",".join(EXTERNAL),
                        help="external qualification subset; native admission requires all five")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--thread-sweep", default=None,
                        help="measure only; defaults to 1,2,4,8,16 physical cores")
    parser.add_argument("--trials", type=int, default=5,
                        help="fresh processes per measure configuration")
    parser.add_argument("--timeout", type=float, default=0.9)
    parser.add_argument("--window-note", default="")
    args = parser.parse_args()
    if args.threads < 1 or not 0 < args.timeout < 1:
        parser.error("positive threads and timeout strictly below one second required")
    if args.mode != "prepare" and not args.window_note:
        parser.error("execution requires the controller's coordinated window note")
    if args.mode in ("native", "measure") and (not args.backend or not args.qualification):
        parser.error("native execution requires explicit backend and passing external qualification")
    if args.mode in ("native", "measure") and args.backend=="native-amd-capped" and args.native_host_limit_bytes is None:
        parser.error("native-amd-capped requires --native-host-limit-bytes")
    if args.trials < 1:
        parser.error("positive fresh-process trial count required")
    if args.mode != "measure" and args.thread_sweep is not None:
        parser.error("thread-sweep requires explicit measure mode")
    peers = args.backends.split(",")
    if not peers or len(peers) != len(set(peers)) or any(p not in EXTERNAL for p in peers):
        parser.error("qualification peers must be distinct existing external backends")
    if args.mode == "measure" and set(peers) != set(EXTERNAL):
        parser.error("matched measurements require all five external peers")
    args.output = args.output.resolve()
    start = time.perf_counter()
    fixture = PrivateFixture(args.fixture, args.policy_sha256)
    build = args.build.resolve()
    available_cores = physical_cores()
    try:
        budgets = ([int(s) for s in (args.thread_sweep or "1,2,4,8,16").split(",")]
                   if args.mode == "measure" else [args.threads])
    except ValueError:
        parser.error("thread-sweep requires integer budgets")
    if not budgets or len(set(budgets)) != len(budgets) or min(budgets) < 1:
        parser.error("distinct positive physical-core budgets required")
    if max(budgets) > len(available_cores):
        parser.error("requested physical-core budget unavailable")
    qualification_hash = None
    if args.mode in ("native", "measure"):
        qualification_hash = verify_qualification(args.qualification, fixture, build)
    stage = prepare(fixture, args.output)
    document = dict(mode=args.mode, privacy="private-local-only; no numeric publication",
        source_revision=subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        source_status=subprocess.check_output(
            ["git", "-C", str(ROOT), "status", "--porcelain"], text=True),
        fixture_id=fixture.manifest["id"], manifest_sha256=fixture.manifest_sha256,
        policy_sha256=fixture.policy_sha256, window_note=args.window_note,
        qualification_sha256=qualification_hash,
        required_external_peers=list(EXTERNAL),
        comparison_scope=(list(EXTERNAL)+[args.backend] if args.mode == "measure"
                          else peers if args.mode == "qualify" else [args.backend]),
        requested_thread_budgets=budgets,
        requested_native_host_limit_bytes=args.native_host_limit_bytes,
        fresh_process_trials=args.trials if args.mode == "measure" else 1,
        comparison_policy="Numerical gates unchanged; existing five-peer family "
        "comparison scope includes direct and IPC UMFPACK separately.",
        script_sha256={p.name: sha256(p) for p in [Path(__file__),
                     HERE/"real_fixture_audit.py", HERE/"run.py"]},
        preparation_and_prerequisite_audit_wall_s=time.perf_counter()-start,
        attempts=[], timing_scope="child includes input/reference I/O, original-matrix "
        "checks, candidate export and teardown; exact rational audit is separate, "
        "explicitly timed. Child-only time is not certificate-inclusive workflow.",
        claims_unavailable=["forward-error certificate", "continuous accuracy",
                            "performance/frontier admission"])
    backends = peers if args.mode == "qualify" else (args.backend,)
    if args.mode == "measure":
        raw_path = args.output / "attempts.jsonl"
        with raw_path.open("x") as raw:
            for backend, nt, _, _, _, trial in matched_points(
                    list(EXTERNAL)+[args.backend], budgets, args.trials,
                    ["rhs4"], ["training"], ["N"]):
                folder = args.output / f"threads-{nt}-trial-{trial}-{backend}"
                folder.mkdir(mode=0o700)
                result = attempt(backend, build, stage, folder, fixture,
                                 nt, available_cores[:nt], args.timeout, args.native_host_limit_bytes)
                result.update(threads=nt, trial=trial)
                raw.write(json.dumps(result)+"\n")
                raw.flush()
                document["attempts"].append(result)
        raw_path.chmod(0o600)
        document["raw_attempts_sha256"] = sha256(raw_path)
    elif args.mode != "prepare":
        for backend in backends:
            document["attempts"].append(attempt(backend, build, stage, args.output,
                fixture, args.threads, available_cores[:args.threads], args.timeout, args.native_host_limit_bytes))
    document["status"] = ("prepared" if args.mode == "prepare" else "pass"
        if all(a["status"] == "pass" for a in document["attempts"]) else "failed")
    document["adapter_wall_s"] = time.perf_counter()-start
    target = args.output / "receipt.json"
    target.write_text(json.dumps(document, indent=2)+"\n")
    target.chmod(0o600)
    print(json.dumps(dict(status=document["status"], receipt=str(target),
                         attempts=len(document["attempts"]),
                         outcomes=dict(collections.Counter(
                             a["status"] for a in document["attempts"])),
                         backends={a["backend"]: a["status"]
                                   for a in document["attempts"]})))
    return 0 if document["status"] in ("prepared", "pass") else 2


if __name__ == "__main__":
    raise SystemExit(main())
