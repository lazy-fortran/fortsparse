#!/usr/bin/env python3
"""Combine raw runs into executable-consistent, auditable CPU admission evidence."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import sys

from benchmark_analysis import (CASE_KEYS, KEYS, add_objectives, analyze_admission,
                                summarize)

DEFAULT_CASES = ("band:1,band:2,band:3,grid2:1,grid2:2,grid2:3,grid2:7,"
                 "grid3:1,grid3:2,grid3:3,grid3:7,bcsstk09:1")
DEFAULT_BACKENDS = "superlu,umfpack,mumps,mkl-pardiso,native,native-nd"
DEFAULT_WORKLOADS = "rhs1,rhs8,refactor,derivatives,warm"
REPORT_METRICS = ("admission_s", "workflow_cpu_s", "parent_rss_kib",
                  "helper_rss_kib", "parent_cpu_s", "helper_cpu_s", "total_s",
                  "process_wall_s", "requested_threads", "active_threads",
                  "ordering_used", "actual_device", "factor_storage_bytes",
                  "workspace_bytes", "device_bytes", "backward_error")


def source_metadata(path, document, payload):
    metadata = {key: value for key, value in document.items()
                if key not in ("records", "summaries", "admission", "csv_exports")}
    metadata.update(path=str(path.resolve()), sha256=hashlib.sha256(payload).hexdigest(),
                    input_status_counts=dict(Counter(r["status"] for r in document["records"])))
    return metadata


def load_runs(paths):
    """Later inputs replace configurations; changed identity replaces whole series."""
    groups, identities, sources = {}, {}, []
    for path in paths:
        payload = path.read_bytes()
        document = json.loads(payload)
        source = source_metadata(path, document, payload)
        sources.append(source)
        current = {}
        for original in document["records"]:
            record = dict(original)
            record.setdefault("variant", "training")
            for key in ("level", "threads"):
                record[key] = int(record[key])
            if record["status"] == "ok":
                add_objectives(record)
            key = tuple(record[k] for k in KEYS)
            current.setdefault(key, []).append(record)
        for key, records in current.items():
            series = tuple(value for name, value in zip(KEYS, key) if name != "threads")
            executable = "bench_spd_pardiso" if key[0] == "mkl-pardiso" else "bench_spd"
            binary = source.get("binary_provenance", {}).get(executable, {}).get("sha256")
            # Missing binary identity cannot authorize merging separate inputs.
            identity = (binary or source["sha256"], json.dumps({
                field: source.get(field) for field in
                ("fixed_environment", "cpu_model", "platform", "physical_cpu_ids")},
                sort_keys=True), source.get("matrix_provenance", {}).get("matrix_sha256")
                if source.get("matrix_provenance") else None)
            if identities.get(series) != identity:
                for old in list(groups):
                    if tuple(value for name, value in zip(KEYS, old)
                             if name != "threads") == series:
                        del groups[old]
                identities[series] = identity
            groups[key] = dict(records=records, source_sha256=source["sha256"],
                               binary_sha256=binary, source_revision=source.get("revision"))
    records = [record for group in groups.values() for record in group["records"]]
    summaries = summarize(records)
    for row in summaries:
        group = groups[tuple(row[k] for k in KEYS)]
        row.update({key: group[key] for key in
                    ("source_sha256", "binary_sha256", "source_revision")})
        row["n"] = next((r.get("n") for r in group["records"] if r.get("n")), None)
    return summaries, sources


def expected_scope(cases, workloads, variants):
    # The archived structural matrix is a fixed case, with no invented variant.
    return [dict(family=family, level=level, workload=workload, variant=variant)
            for family, level in cases for workload in workloads for variant in variants
            if not (family == "bcsstk09" and variant != "training")]


def coverage(summaries, scope, backends, threads, trials):
    missing, trial_mismatches = [], []
    for case in scope:
        for backend in backends:
            rows = [r for r in summaries if r["backend"] == backend and
                    all(r[k] == case[k] for k in CASE_KEYS)]
            absent = sorted(set(threads) - {r["threads"] for r in rows})
            if absent:
                missing.append(dict(case, backend=backend, missing_threads=absent))
            for row in rows:
                if row["threads"] in threads and row["samples"] != trials:
                    trial_mismatches.append(dict(case, backend=backend,
                        threads=row["threads"], samples=row["samples"], expected=trials))
    present = {tuple(row[k] for k in CASE_KEYS) for row in summaries}
    missing_cases = [case for case in scope if tuple(case[k] for k in CASE_KEYS) not in present]
    return dict(complete=not missing and not trial_mismatches, expected_case_groups=len(scope),
                expected_backends=backends, expected_threads=threads, expected_trials=trials,
                missing_case_groups=missing_cases, missing_configurations=missing,
                trial_count_mismatches=trial_mismatches)


def describe(row):
    return {key: row.get(key) for key in (*KEYS, "n", "samples", "successes", "statuses",
            "source_sha256", "binary_sha256", "source_revision", *REPORT_METRICS)}


def build_report(summaries, sources, scope, backends, threads, trials):
    declared_cases = {tuple(case[k] for k in CASE_KEYS) for case in scope}
    scoped = [row for row in summaries if tuple(row[k] for k in CASE_KEYS) in declared_cases
              and row["backend"] in backends and row["threads"] in threads]
    admission = analyze_admission(scoped, backends)
    case_rows = []
    for selected in admission["best_native_by_case"]:
        key = tuple(selected[k] for k in CASE_KEYS)
        row = {k: selected[k] for k in CASE_KEYS}
        row.update(status=selected["status"], alternatives=selected["alternatives"])
        for field, role in (("selected_native", "native"), ("matched_comparator", "comparator")):
            choice = selected[field]
            match = next((r for r in summaries if tuple(r[k] for k in CASE_KEYS) == key
                          and choice and r["backend"] == choice["backend"]
                          and r["threads"] == choice["threads"]), None)
            row[role] = describe(match) if match else None
        case_rows.append(row)
    found = {tuple(r[k] for k in CASE_KEYS) for r in case_rows}
    for case in scope:
        if tuple(case[k] for k in CASE_KEYS) not in found:
            case_rows.append(dict(case, status="indeterminate_missing_case",
                                  alternatives=[], native=None, comparator=None))
    completeness = coverage(summaries, scope, backends, threads, trials)
    return dict(schema="fortsparse-spd-report-v1", mvp_claim="not_evaluated",
                interpretation=admission["interpretation"],
                native_selection_policy=admission["native_selection_policy"],
                scope_admission_status=(admission["best_native_overall"] if
                    completeness["complete"] else "indeterminate_incomplete_evidence"),
                coverage=completeness, sources=sources,
                out_of_scope_case_groups=[dict(zip(CASE_KEYS, case)) for case in sorted({
                    tuple(row[k] for k in CASE_KEYS) for row in summaries} - declared_cases)],
                out_of_scope_configurations=len(summaries) - len(scoped),
                selected_status_counts=dict(sum((Counter(r["statuses"]) for r in summaries),
                                                Counter())),
                case_status_counts=dict(Counter(r["status"] for r in case_rows)),
                admission=admission, cases=case_rows,
                configurations=[describe(row) for row in summaries],
                measurement_limits=[
                    "Elapsed admission phases and whole-workflow user+system CPU are distinct.",
                    "Whole-workflow CPU includes setup, warmups, checks, cleanup and helper CPU.",
                    "Parent/helper RSS are separate; aggregate UMFPACK peak memory is unknown.",
                    "Observed min/max ranges are not confidence intervals or proof of leadership.",
                    "Raw source archives retain superseded trials; their status counts remain in sources.",
                    "Shared-library/helper identities are limited to provenance supplied by each raw run."])


def flatten(row, prefix=""):
    result = {}
    for key, value in row.items():
        name = prefix + key
        if isinstance(value, dict) and key != "statuses":
            result.update(flatten(value, name + "_"))
        elif isinstance(value, (dict, list)):
            result[name] = json.dumps(value, separators=(",", ":"))
        else:
            result[name] = value
    return result


def write_report(report, output):
    exports = {}
    for field in ("cases", "configurations"):
        rows = [flatten(row) for row in report[field]]
        path = output.with_suffix(f".{field}.csv")
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=sorted({k for row in rows for k in row}),
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        exports[field] = dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    report["csv_exports"] = exports
    output.write_text(json.dumps(report, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path, help="raw run JSONs, oldest first")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cases", default=DEFAULT_CASES)
    parser.add_argument("--backends", default=DEFAULT_BACKENDS)
    parser.add_argument("--workloads", default=DEFAULT_WORKLOADS)
    parser.add_argument("--variants", default="training,heldout")
    parser.add_argument("--threads", default="1,2,4,8,16")
    parser.add_argument("--trials", type=int, default=5)
    args = parser.parse_args()
    cases = [(family, int(level)) for family, level in
             (entry.split(":") for entry in args.cases.split(","))]
    threads = [int(t) for t in args.threads.split(",")]
    if args.trials < 1 or not threads or min(threads) < 1:
        parser.error("positive thread counts and trial count required")
    scope = expected_scope(cases, args.workloads.split(","), args.variants.split(","))
    summaries, sources = load_runs(args.inputs)
    report = build_report(summaries, sources, scope, args.backends.split(","),
                          threads, args.trials)
    report["command"] = sys.argv
    report["analysis_source_sha256"] = {name: hashlib.sha256(
        Path(__file__).with_name(name).read_bytes()).hexdigest()
        for name in ("report_spd.py", "benchmark_analysis.py")}
    write_report(report, args.output)
    print(f"{len(summaries)} configurations; {len(report['cases'])} case groups; "
          f"scope admission: {report['scope_admission_status']}; output: {args.output}")


if __name__ == "__main__":
    main()
