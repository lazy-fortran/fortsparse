"""Observed-range GPU comparisons with explicit offload and memory limits."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics

from offload_provenance import attach

CASE_KEYS = ("family", "level", "workload", "variant", "io_mode", "requested_device")
CONFIG_KEYS = ("requested_backend", "threads", "ordering", "repeats", "pool_budget_bytes",
               "executable_sha256")
OBJECTIVES = {
    "wall_host_rss": ("primary_phase_s", "parent_rss_kib"),
    "wall_device_total": ("primary_phase_s", "device_total_peak_bytes"),
    "wall_host_device": ("primary_phase_s", "parent_rss_kib", "device_total_peak_bytes"),
    "cpu_host_rss": ("workflow_cpu_s", "parent_rss_kib"),
}


def interval(records, field):
    values = [r.get(field) for r in records]
    if not values or any(v is None or not math.isfinite(v) or v < 0 for v in values):
        return None
    return dict(minimum=min(values), median=statistics.median(values), maximum=max(values))


def summarize(trials, required_trials=5, validated_trace=False):
    groups = {}
    for trial in trials:
        trial = dict(trial)
        trial.setdefault("requested_backend", trial["backend"])
        trial["executable_sha256"] = trial.get("executable_sha256") or ""
        trial.setdefault("pool_budget_bytes", trial.get("gpu_pool_budget_bytes",
                         trial.get("pool_limit_bytes", -1)))
        key = tuple(trial[k] for k in CASE_KEYS + CONFIG_KEYS)
        groups.setdefault(key, []).append(trial)
    result = []
    fields = set(sum((list(v) for v in OBJECTIVES.values()), [])) | {
        "factor_s", "refactor_s", "solve_s", "total_s", "backward_error",
        "active_threads", "requested_threads", "process_wall_s", "analysis_s", "numeric_s",
        "device_process_sampled_max_bytes", "wrapper_device_bytes", "pinned_host_bytes"}
    for key, records in sorted(groups.items()):
        good = [r for r in records if r["status"] == "ok"]
        row = dict(zip(CASE_KEYS + CONFIG_KEYS, key))
        row.update(samples=len(records), successes=len(good),
                   complete=len(good) == len(records) and len(good) >= required_trials,
                   offload_verified=bool(good) and all(
                       (r.get("gpu_offload_verified") == 1 and r.get("actual_device", -1) >= 0) or
                       (validated_trace and r.get("trace_verified") is True and
                        r.get("trace_actual_device") == r.get("requested_device") and
                        r.get("trace_actual_device", -1) >= 0)
                       for r in good),
                   trace_verified_samples=sum(validated_trace and r.get("trace_verified") is True
                                              for r in good),
                   statuses={s: sum(r["status"] == s for r in records)
                             for s in sorted({r["status"] for r in records})},
                   metrics={field: interval(good, field) for field in sorted(fields)})
        result.append(row)
    return result


def config(row):
    return {k: row[k] for k in CONFIG_KEYS}


def audit_candidate(native, externals, expected, objectives):
    """Unknown native upper bounds cannot establish memory-front membership."""
    if not native["complete"]:
        return dict(status="native_incomplete")
    if not native["offload_verified"]:
        return dict(status="native_offload_unverified")
    if not expected:
        return dict(status="no_external_comparators")
    missing = sorted(set(expected) - {r["requested_backend"] for r in externals})
    if missing:
        return dict(status="missing_comparators", backends=missing)
    incomplete = [config(r) for r in externals if not r["complete"]]
    if incomplete:
        return dict(status="external_incomplete", configurations=incomplete)
    unverified = [config(r) for r in externals if not r["offload_verified"]]
    if unverified:
        return dict(status="external_offload_unverified", configurations=unverified)
    unknown = [m for m in objectives if native["metrics"][m] is None]
    if unknown:
        return dict(status="unknown_native_objectives", metrics=unknown)
    possible, excluded, matches = [], [], []
    for external in externals:
        # An unknown external lower bound cannot exclude a possible dominator.
        separating = [m for m in objectives if external["metrics"][m] is not None
                      and native["metrics"][m]["maximum"] <
                      external["metrics"][m]["minimum"]]
        if separating:
            excluded.append(dict(config(external), separating_objectives=separating))
        else:
            possible.append(config(external))
        known = all(external["metrics"][m] is not None for m in objectives)
        if known and separating and all(native["metrics"][m]["maximum"] <=
                                        external["metrics"][m]["minimum"]
                                        for m in objectives):
            matches.append(config(external))
    return dict(status="possible_dominator" if possible else "observed_range_nondominated",
                possible_dominators=possible, excluded_dominators=excluded,
                same_comparator_match_or_beat=matches)


def analyze(document, required_trials=5, offload_provenance=None):
    if offload_provenance:
        document = attach(document, offload_provenance)
    trials = []
    for original in document["trials"]:
        row = dict(original)
        command = row.get("command", [])
        identity = document.get("binary_provenance", {}).get(
            command[3] if len(command) > 3 else "", {})
        known_hash = identity.get("sha256")
        if known_hash and row.get("executable_sha256", known_hash) != known_hash:
            raise ValueError("trial and document executable hashes disagree")
        row.setdefault("executable_sha256", known_hash)
        trials.append(row)
    rows = summarize(trials, required_trials, validated_trace=bool(offload_provenance))
    expected = sorted({r["requested_backend"] for r in rows
                       if not r["requested_backend"].startswith("native")})
    cases = sorted({tuple(r[k] for k in CASE_KEYS) for r in rows})
    results = []
    for case in cases:
        current = [r for r in rows if tuple(r[k] for k in CASE_KEYS) == case]
        external = [r for r in current if r["requested_backend"] in expected]
        candidates = []
        for row in current:
            if not row["requested_backend"].startswith("native"):
                continue
            candidates.append(dict(configuration=config(row), objectives={
                name: audit_candidate(row, external, expected, metrics)
                for name, metrics in OBJECTIVES.items()}))
        results.append(dict(zip(CASE_KEYS, case), candidates=candidates))
    return dict(schema="fortsparse-gpu-front-audit-v1", required_trials=required_trials,
                measurement_class=document.get("measurement_class", "unknown"),
                admission_eligible=document.get("measurement_class") == "controlled",
                limit="Observed ranges, not confidence intervals. Sampled GPU memory is not a total peak. "
                      "Uncontrolled runs cannot establish performance admission.",
                offload_provenance=document.get("offload_provenance") if offload_provenance else None,
                expected_external_backends=expected, summaries=rows, cases=results)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--required-trials", type=int, default=5)
    parser.add_argument("--offload-provenance", type=Path,
                        help="Revalidate exact-binary/configuration Nsight SQLite evidence offline")
    args = parser.parse_args()
    if args.required_trials < 1:
        parser.error("positive trial count required")
    result = analyze(json.loads(args.input.read_text()), args.required_trials, args.offload_provenance)
    result["source"] = dict(path=str(args.input.resolve()),
                            sha256=hashlib.sha256(args.input.read_bytes()).hexdigest())
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
