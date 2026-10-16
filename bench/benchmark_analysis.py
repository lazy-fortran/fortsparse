"""Observed-range comparisons, not statistical confidence intervals or proof."""
import csv
import json
from pathlib import Path
import statistics

METRICS = (
    "factor_s", "solve_s", "refactor_s", "total_s", "parent_cpu_s",
    "parent_rss_kib", "helper_rss_kib", "helper_cpu_s", "backward_error",
    "process_wall_s", "factor_solve_s", "workflow_cpu_s", "jvp_s", "vjp_s",
    "cold_factor_s", "factor_nnz", "factor_mem_mb", "analysis_s", "numeric_s",
    "admission_s", "derivative_reuse_s", "requested_threads", "active_threads",
    "requested_device", "actual_device", "profile_iterations", "profile_s",
    "ordering_used", "symbolic_nnz", "factor_storage_bytes", "workspace_bytes",
    "device_bytes")
KEYS = ("backend", "family", "level", "threads", "workload", "variant")
CASE_KEYS = ("family", "level", "workload", "variant")
# These aliases select backend 4 directly in bench_spd; neither calls a fallback.
PRODUCTION_CPU_CONFIGS = ("native", "native-nd")


def add_objectives(record):
    """Keep costs explicit; never compare refactor or derivatives as primal-only."""
    record["factor_solve_s"] = record["factor_s"] + record["solve_s"]
    record["workflow_cpu_s"] = record["parent_cpu_s"] + record["helper_cpu_s"]
    record["derivative_reuse_s"] = -1.0
    workload = record["workload"]
    if workload in ("refactor", "profile-refactor"):
        record["admission_s"] = record["refactor_s"] + record["solve_s"]
    elif workload == "derivatives":
        record["derivative_reuse_s"] = record["jvp_s"] + record["vjp_s"]
        record["admission_s"] = (record["factor_solve_s"] +
                                  record["derivative_reuse_s"])
    else:
        record["admission_s"] = record["factor_solve_s"]


def summarize(records):
    groups = {}
    for record in records:
        record.setdefault("variant", "training")
        key = tuple(record[k] for k in KEYS)
        groups.setdefault(key, []).append(record)
    summaries = []
    for key, group in sorted(groups.items()):
        good = [r for r in group if r["status"] == "ok"]
        summary = dict(zip(KEYS, key), successes=len(good), samples=len(group))
        summary["statuses"] = {status: sum(r["status"] == status for r in group)
                               for status in sorted({r["status"] for r in group})}
        for metric in METRICS:
            values = [r[metric] for r in good if r.get(metric) is not None
                      and r[metric] >= 0]
            summary[metric] = (dict(median=statistics.median(values),
                                    minimum=min(values), maximum=max(values),
                                    available_samples=len(values))
                               if values else None)
        summaries.append(summary)
    for summary in summaries:
        serial = next((s for s in summaries if all(s[k] == summary[k]
                       for k in KEYS if k != "threads") and s["threads"] == 1
                       and s["successes"] == s["samples"]), None)
        if serial and summary["successes"] == summary["samples"]:
            for metric in ("total_s", "factor_solve_s", "admission_s"):
                if not serial[metric] or not summary[metric]:
                    continue
                denominator = summary[metric]["median"]
                if denominator <= 0:
                    continue
                speedup = serial[metric]["median"] / denominator
                summary[metric]["speedup"] = speedup
                summary[metric]["parallel_efficiency"] = speedup / summary["threads"]
    return summaries


def reference(row):
    return dict(backend=row["backend"], threads=row["threads"])


def interval_relation(native, comparator):
    if native["maximum"] <= comparator["minimum"]:
        return "no_worse_over_observed_ranges"
    if native["minimum"] > comparator["maximum"]:
        return "worse_over_observed_ranges"
    return "overlapping_observed_ranges"


def analyze_admission(summaries, expected_backends):
    cases = sorted({tuple(s[k] for k in ("family", "level", "workload", "variant"))
                    for s in summaries})
    external = [b for b in expected_backends if not b.startswith("native")]
    native_variants = [b for b in expected_backends if b.startswith("native")]
    results = []
    for case in cases:
        if case[2].startswith("profile-"):
            continue
        rows = [s for s in summaries if tuple(s[k] for k in
                ("family", "level", "workload", "variant")) == case]
        valid = [s for s in rows if s["successes"] == s["samples"]
                 and s["admission_s"] and s["parent_rss_kib"]]
        known = [s for s in valid if s["backend"] in external
                 and s["backend"] != "umfpack"]
        unknown = [s for s in valid if s["backend"] == "umfpack"]
        missing = [b for b in external if not any(s["backend"] == b for s in valid)]

        def dominates(a, b):
            at, bt = a["admission_s"]["median"], b["admission_s"]["median"]
            am, bm = a["parent_rss_kib"]["median"], b["parent_rss_kib"]["median"]
            return at <= bt and am <= bm and (at < bt or am < bm)

        frontier = [s for s in known if not any(dominates(other, s)
                    for other in known if other is not s)]
        protected = []
        uncertain = []
        for candidate in frontier:
            blockers = [u for u in unknown if not (
                candidate["admission_s"]["maximum"] < u["admission_s"]["minimum"]
                or candidate["parent_rss_kib"]["maximum"] <
                u["parent_rss_kib"]["minimum"])]
            if blockers:
                uncertain.append(dict(comparator=reference(candidate),
                                      unknown_memory_blockers=[reference(u) for u in blockers]))
            else:
                protected.append(candidate)
        for backend in native_variants:
            native = [s for s in valid if s["backend"] == backend]
            comparisons = []
            for n in native:
                for c in protected:
                    time = interval_relation(n["admission_s"], c["admission_s"])
                    memory = interval_relation(n["parent_rss_kib"], c["parent_rss_kib"])
                    comparisons.append(dict(native=reference(n), comparator=reference(c),
                                            time=time, memory=memory))
            passed = [c for c in comparisons if c["time"] ==
                      "no_worse_over_observed_ranges" and c["memory"] ==
                      "no_worse_over_observed_ranges"]
            possible = [c for c in comparisons if "worse_over_observed_ranges"
                        not in (c["time"], c["memory"])]
            if not native:
                status = "failed_or_unavailable"
            elif not external or missing or not protected:
                status = "indeterminate"
            elif passed:
                status = "pass_over_observed_ranges"
            elif possible:
                status = "indeterminate_overlap"
            else:
                status = "not_competitive_over_observed_ranges"
            cuda_verified = backend != "native-cuda" or any(
                n["actual_device"] is not None and
                n["actual_device"]["available_samples"] == n["samples"]
                for n in native)
            results.append(dict(zip(("family", "level", "workload", "variant"), case),
                backend=backend, status=status, production_path=backend != "native-reference",
                device_execution_verified=cuda_verified, missing_comparators=missing,
                known_frontier=[reference(c) for c in frontier],
                protected_frontier=[reference(c) for c in protected],
                uncertain_frontier=uncertain,
                unknown_aggregate_memory=[reference(c) for c in unknown],
                comparisons=comparisons))
    overall = {}
    for backend in native_variants:
        groups = [r for r in results if r["backend"] == backend]
        states = {r["status"] for r in groups}
        if "failed_or_unavailable" in states:
            status = "failed_or_unavailable"
        elif "not_competitive_over_observed_ranges" in states:
            status = "not_competitive_over_observed_ranges"
        elif not groups or any(s.startswith("indeterminate") for s in states):
            status = "indeterminate"
        else:
            status = "pass_over_observed_ranges"
        if any(not r["device_execution_verified"] for r in groups):
            status = "indeterminate_device_execution"
        overall[backend] = status
    best_by_case = best_native_by_case(summaries, results)
    return dict(overall=overall, cases=results,
                best_native_by_case=best_by_case,
                best_native_overall=aggregate_status(best_by_case),
                native_selection_policy=("Explicit measured production CPU ordering/thread "
                    "configuration per case and workload. This union is not an AUTO/default "
                    "performance claim; reference, CUDA and external fallbacks are excluded."),
                interpretation="Empirical ranges of the retained trials; not confidence intervals.",
                mvp_claim="not_evaluated")


def aggregate_status(groups):
    states = {r["status"] for r in groups}
    if "failed_or_unavailable" in states:
        return "failed_or_unavailable"
    if "not_competitive_over_observed_ranges" in states:
        return "not_competitive_over_observed_ranges"
    if not groups or any(s.startswith("indeterminate") for s in states):
        return "indeterminate"
    return "pass_over_observed_ranges"


def best_native_by_case(summaries, results):
    """Union strict SAME-comparator passes; retain every evaluated case group."""
    cases = sorted({tuple(s[k] for k in CASE_KEYS) for s in summaries
                    if not s["workload"].startswith("profile-")})
    selected = []
    for case in cases:
        rows = [s for s in summaries if tuple(s[k] for k in CASE_KEYS) == case]
        alternatives = [r for r in results if tuple(r[k] for k in CASE_KEYS) == case
                        and r["backend"] in PRODUCTION_CPU_CONFIGS]
        passes = []
        for result in alternatives:
            # A pair may numerically pass while its external corpus is incomplete.
            if result["status"] != "pass_over_observed_ranges":
                continue
            for pair in result["comparisons"]:
                if any(pair[axis] != "no_worse_over_observed_ranges"
                       for axis in ("time", "memory")):
                    continue
                native = next(s for s in rows if reference(s) == pair["native"])
                comparator = next(s for s in rows if reference(s) == pair["comparator"])
                passes.append((native, comparator, pair))
        entry = dict(zip(CASE_KEYS, case), alternatives=[
            dict(backend=r["backend"], status=r["status"]) for r in alternatives])
        if passes:
            native, comparator, pair = min(passes, key=lambda p: (
                p[0]["admission_s"]["median"], p[0]["parent_rss_kib"]["median"],
                p[0]["backend"], p[0]["threads"],
                p[1]["backend"], p[1]["threads"]))
            entry.update(status="pass_over_observed_ranges", selected_native=reference(native),
                         matched_comparator=reference(comparator), comparison=pair,
                         ordering_used=native.get("ordering_used"),
                         native_time_s=native["admission_s"],
                         native_parent_rss_kib=native["parent_rss_kib"],
                         comparator_time_s=comparator["admission_s"],
                         comparator_parent_rss_kib=comparator["parent_rss_kib"])
        else:
            states = {r["status"] for r in alternatives}
            if any(s.startswith("indeterminate") for s in states):
                status = "indeterminate"
            elif "not_competitive_over_observed_ranges" in states:
                status = "not_competitive_over_observed_ranges"
            else:
                status = "failed_or_unavailable"
            entry.update(status=status, selected_native=None, matched_comparator=None)
        selected.append(entry)
    return selected


def export_csv(result, output):
    """Two flat tables alongside JSON; every failure and raw trial is retained."""
    paths = {}
    for name, rows in (("trials", result["records"]), ("summaries", result["summaries"])):
        flattened = []
        for row in rows:
            flat = {}
            for key, value in row.items():
                if isinstance(value, dict) and key in METRICS:
                    flat.update({f"{key}_{field}": item for field, item in value.items()})
                elif isinstance(value, (dict, list)):
                    flat[key] = json.dumps(value, separators=(",", ":"))
                else:
                    flat[key] = value
            flattened.append(flat)
        path = Path(output).with_suffix(f".{name}.csv")
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=sorted({k for row in flattened for k in row}))
            writer.writeheader()
            writer.writerows(flattened)
        paths[name] = str(path)
    return paths
