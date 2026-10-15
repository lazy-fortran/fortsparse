#!/usr/bin/env python3
"""Render measured sparse-solver evidence; PNGs stay outside source history."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import statistics
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

# Okabe–Ito palette, with redundant markers and line styles.
STYLES = {
    "native": ("#0072B2", "*", "Native CPU"),
    "native-cuda": ("#D55E00", (8, 1, 0), "Native CUDA"),
    "native-reference": ("#777777", "x", "Native reference"),
    "superlu": ("#56B4E9", "s", "SuperLU"),
    "umfpack": ("#E69F00", "o", "UMFPACK parent"),
    "mumps": ("#009E73", "^", "MUMPS"),
    "mkl-pardiso": ("#CC79A7", "D", "Intel MKL PARDISO"),
}
KEYS = ("backend", "family", "level", "threads", "workload", "variant")


def phase_value(record):
    if record["workload"] == "refactor":
        return record["refactor_s"] + record["solve_s"]
    value = record["factor_s"] + record["solve_s"]
    if record["workload"] == "derivatives":
        value += record["jvp_s"] + record["vjp_s"]
    return value


def phase_label(workload):
    return {"rhs1": "Cold factor + solve", "rhs8": "Cold factor + 8 RHS",
            "warm": "Warm-runtime factor + solve", "refactor": "Refactor + solve",
            "derivatives": "Factor + primal + JVP + VJP"}.get(workload, workload)


def interval(values, scale=1.0):
    if not values:
        return None
    return tuple(scale*x for x in (statistics.median(values), min(values), max(values)))


def load_groups(paths):
    groups = {}
    sources = []
    series_binary = {}
    for path in paths:
        document = json.loads(path.read_text())
        source = dict(path=str(path.resolve()), name=path.stem,
                      sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                      revision=document.get("revision"),
                      dirty=document.get("working_tree_dirty"),
                      cpu_model=document.get("cpu_model"),
                      binary_hashes={name: info.get("sha256") for name, info in
                                     document.get("binary_provenance", {}).items()},
                      build_provenance=document.get("build_provenance"))
        sources.append(source)
        current = {}
        for record in document["records"]:
            record = dict(record)
            record.setdefault("variant", "training")
            key = tuple(record[k] for k in KEYS)
            current.setdefault(key, []).append(record)
        # Keep scaling curves within one executable. Replacing a binary
        # replaces its whole thread series; trials never pool across inputs.
        for key, records in current.items():
            series = tuple(value for index, value in enumerate(key) if KEYS[index] != "threads")
            binary_name = "bench_spd_pardiso" if key[0] == "mkl-pardiso" else "bench_spd"
            binary = source["binary_hashes"].get(binary_name) or source["sha256"]
            if series_binary.get(series) != binary:
                groups = {old: value for old, value in groups.items() if tuple(
                    item for index, item in enumerate(old) if KEYS[index] != "threads") != series}
                series_binary[series] = binary
            good = [r for r in records if r["status"] == "ok"]
            group = dict(zip(KEYS, key), records=records, good=good,
                         samples=len(records), successes=len(good), source=source)
            group["wall"] = interval([phase_value(r) for r in good], 1000)
            group["cpu"] = interval([r["parent_cpu_s"] + r["helper_cpu_s"]
                                     for r in good], 1000)
            group["rss"] = interval([r["parent_rss_kib"] for r in good], 1/1024)
            group["helper"] = interval([r["helper_rss_kib"] for r in good], 1/1024)
            group["device"] = interval([r["device_bytes"] for r in good
                                         if r.get("device_bytes", -1) >= 0], 1/2**20)
            group["active"] = sorted({int(r["active_threads"]) for r in good
                                       if r.get("active_threads", -1) >= 0})
            group["failures"] = {s: sum(r["status"] == s for r in records)
                                  for s in sorted({r["status"] for r in records})
                                  if s != "ok"}
            groups[key] = group
    return list(groups.values()), sources


def choose(groups, case, workload):
    return [g for g in groups if (g["family"], g["level"]) == case
            and g["workload"] == workload]


def case_label(groups, case):
    rows = [r for g in groups if (g["family"], g["level"]) == case for r in g["good"]]
    size = f", n={int(rows[0]['n']):,}" if rows else ""
    return {"grid2": "Anisotropic 2D", "grid3": "3D diffusion",
            "band": "Banded SPD"}.get(case[0], case[0]) + f" · level {int(case[1])}" + size


def errorbars(axis, x, y, backend, **kwargs):
    color, marker, _ = STYLES.get(backend, ("#000000", "+", backend))
    hollow = backend in ("umfpack", "native-cuda")
    axis.errorbar(x[0], y[0], xerr=[[x[0]-x[1]], [x[2]-x[0]]],
                  yerr=[[y[0]-y[1]], [y[2]-y[0]]], linestyle="none",
                  color=color, marker=marker, markersize=12 if backend.startswith("native") else 6,
                  markerfacecolor="white" if hollow else color, capsize=2, alpha=0.9, **kwargs)


def legend_handles(groups):
    result = []
    for backend, (color, marker, label) in STYLES.items():
        rows = [g for g in groups if g["backend"] == backend]
        if not rows:
            continue
        counts = "/".join(map(str, sorted({g["samples"] for g in rows})))
        result.append(Line2D([], [], color=color, marker=marker, linestyle="none",
                             markerfacecolor="white" if backend in ("umfpack", "native-cuda") else color,
                             markersize=10 if backend.startswith("native") else 6,
                             label=f"{label} ({counts} trials/config.)"))
    return result


def footer(figure, sources, note, rect=(0.04, 0.16, 0.99, 0.88)):
    source_text = "; ".join(f"{s['name']}@{str(s['revision'])[:7]}" +
                            ("+dirty" if s["dirty"] else "") for s in sources)
    text = textwrap.fill(note, 155) + "\n" + textwrap.fill("Inputs: " + source_text, 155)
    figure.text(0.04, 0.025, text, ha="left", va="bottom", fontsize=8, color="#444444")
    figure.tight_layout(rect=rect)


def pareto_plot(groups, sources, cases, workload, title, output):
    fig, axes = plt.subplots(len(cases), 2, figsize=(14.2, 4.3*len(cases)+1.6), squeeze=False)
    plotted = []
    for row, case in enumerate(cases):
        selected = choose(groups, case, workload)
        plotted += selected
        for column, metric in enumerate(("wall", "cpu")):
            axis = axes[row, column]
            points, helper, failures = [], [], []
            for group in selected:
                if group[metric] and group["rss"]:
                    errorbars(axis, group["rss"], group[metric], group["backend"])
                    if group["backend"] not in ("umfpack", "native-cuda") and not group["failures"]:
                        points.append((group["rss"][0], group[metric][0]))
                    if group["backend"] == "umfpack":
                        helper.append(group["helper"][0])
                if group["failures"]:
                    failures.append(f"{group['backend']} {int(group['threads'])}t: {group['failures']}")
            frontier = [p for p in points if not any(q[0] <= p[0] and q[1] <= p[1]
                        and q != p for q in points)]
            if frontier:
                frontier = sorted(set(frontier))
                axis.plot(*zip(*frontier), color="#333333", linestyle="--", linewidth=1)
            helper_note = (f"\nUMFPACK helper: {min(helper):.1f}–{max(helper):.1f} MiB, separate"
                           if helper else "")
            if failures:
                axis.text(0.02, 0.98, "\n".join(failures), transform=axis.transAxes,
                          va="top", fontsize=8, color="#AA3322")
            axis.set(xscale="log", yscale="log", title=case_label(groups, case)+helper_note,
                     xlabel="Peak parent RSS (MiB)", ylabel=(phase_label(workload)+" wall time (ms)"
                     if metric == "wall" else "Whole-workflow CPU time (ms)"))
            axis.grid(True, which="both", alpha=0.18)
    fig.suptitle(title+" · time / host-memory Pareto", fontsize=16, y=0.985)
    fig.legend(handles=legend_handles(plotted) + [Line2D([], [], color="#333333", linestyle="--",
               label="CPU single-process frontier")], loc="upper center", bbox_to_anchor=(0.5, 0.94), ncol=3, fontsize=9)
    footer(fig, sources, "Medians and observed min–max; each point is one requested thread configuration. "
           "Frontier uses CPU single-process measurements only. UMFPACK's aggregate memory is unknown. "
           "GPU host RSS excludes device allocation; CUDA is plotted separately from this CPU frontier. "
           "CPU time = parent + helper user/system time for the whole workflow, including warmup/setup/checks/free; "
           "it is not factor/solve phase CPU time.")
    fig.savefig(output, dpi=160)
    plt.close(fig)


def scaling_plot(groups, sources, cases, workload, title, output):
    fig, axes = plt.subplots(len(cases), 3, figsize=(16, 4.2*len(cases)+1.7), squeeze=False)
    plotted = []
    for row, case in enumerate(cases):
        selected = choose(groups, case, workload)
        plotted += selected
        for backend in STYLES:
            series = sorted([g for g in selected if g["backend"] == backend and g["wall"]],
                            key=lambda g: g["threads"])
            if not series:
                continue
            color, marker, _ = STYLES[backend]
            for column, metric in enumerate(("wall", "cpu", "rss")):
                axis = axes[row, column]
                values = [g[metric] for g in series]
                x = [g["threads"] for g in series]
                y = [v[0] for v in values]
                axis.errorbar(x, y, yerr=[[v[0]-v[1] for v in values], [v[2]-v[0] for v in values]],
                              color=color, marker=marker, linestyle="--" if backend == "umfpack" else "-",
                              markerfacecolor="white" if backend in ("umfpack", "native-cuda") else color,
                              markersize=10 if backend.startswith("native") else 5, capsize=2)
                if backend == "native" and column == 0:
                    serial = next((g for g in series if g["threads"] == 1), None)
                    if serial:
                        axis.plot(x, [serial[metric][0]/t for t in x], ":", color=color, alpha=0.65)
                        last = series[-1]
                        speedup = serial[metric][0]/last[metric][0]
                        axis.text(0.02, 0.03, f"CPU {int(last['threads'])}t: {speedup:.2f}×; "
                                  f"efficiency {100*speedup/last['threads']:.0f}%", transform=axis.transAxes,
                                  fontsize=8, color=color)
                    for g in series:
                        if g["active"] and min(g["active"]) < g["threads"]:
                            axis.annotate("active "+"/".join(map(str, g["active"])),
                                          (g["threads"], g[metric][0]), xytext=(0, 9),
                                          textcoords="offset points", ha="center", fontsize=7, color=color)
            if backend == "umfpack":
                axes[row, 2].plot([g["threads"] for g in series], [g["helper"][0] for g in series],
                                  ":", color=color, alpha=0.8)
                axes[row, 2].text(0.02, 0.95, "UMFPACK dotted = helper RSS, separate",
                                  transform=axes[row, 2].transAxes, fontsize=7)
        for column, axis in enumerate(axes[row]):
            axis.set(xscale="log", yscale="log", xticks=[1, 2, 4, 8, 16],
                     xticklabels=[1, 2, 4, 8, 16], xlabel="Requested physical CPU threads",
                     ylabel=(phase_label(workload)+" wall (ms)", "Whole-workflow CPU time (ms)",
                             "Peak parent RSS (MiB)")[column],
                     title=case_label(groups, case))
            axis.grid(True, which="both", alpha=0.18)
        failed = [f"{g['backend']} {int(g['threads'])}t {g['failures']}" for g in selected if g["failures"]]
        if failed:
            axes[row, 0].text(0.02, 0.98, "\n".join(failed), transform=axes[row, 0].transAxes,
                              va="top", fontsize=7, color="#AA3322")
    fig.suptitle(title+" · thread scaling", fontsize=16, y=0.99)
    fig.legend(handles=legend_handles(plotted), loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=4, fontsize=9)
    footer(fig, sources, "Observed medians/min–max. Dotted native line = ideal 1/t from the measured one-thread wall time. "
           "CPU time sums parent and helper user+system time over the whole invocation, including setup, checks, cleanup "
           "and warmup factors for warm runs; it is not phase CPU time. Missing thread counts are unmeasured. "
           "Parent/helper RSS are not summed; device allocations are excluded.", rect=(0.05, 0.17, 0.99, 0.89))
    fig.savefig(output, dpi=160)
    plt.close(fig)


def best_group(groups, case, workload, backend):
    selected = [g for g in choose(groups, case, workload) if g["backend"] == backend]
    complete = [g for g in selected if g["wall"] and not g["failures"]]
    if complete:
        return min(complete, key=lambda g: g["wall"][0])
    partial = [g for g in selected if g["wall"]]
    return min(partial, key=lambda g: g["wall"][0]) if partial else next(iter(selected), None)


def cpu_gpu_plot(groups, sources, cases, title, output):
    cases = [case for case in cases if any(g["backend"] == "native-cuda" and
             (g["family"], g["level"]) == case for g in groups)] or cases[:2]
    entries = [(case, workload) for case in cases for workload in ("rhs1", "refactor")]
    fig, axes = plt.subplots(1, 3, figsize=(16, max(5.5, 0.65*len(entries)+2.8)), sharey=True)
    labels = []
    for index, (case, workload) in enumerate(entries):
        labels.append(case_label(groups, case)+"\n"+phase_label(workload))
        for backend, offset in (("native", -0.16), ("native-cuda", 0.16)):
            group = best_group(groups, case, workload, backend)
            color, marker, _ = STYLES[backend]
            for column, metric in enumerate(("wall", "rss", "device")):
                axis = axes[column]
                y = index + offset
                value = group[metric] if group else None
                if value:
                    if metric == "device" and value[0] == 0:
                        axis.text(0.02, y, "0 MiB", transform=axis.get_yaxis_transform(),
                                  fontsize=8, color=color, va="center")
                        continue
                    axis.errorbar(value[0], y, xerr=[[value[0]-value[1]], [value[2]-value[0]]],
                                  color=color, marker=marker, markersize=11, capsize=3, linestyle="none",
                                  markerfacecolor="white" if backend == "native-cuda" else color)
                    if column == 0:
                        axis.annotate(f" {int(group['threads'])}t", (value[0], y),
                                      textcoords="offset points", xytext=(5, 4), fontsize=7, color=color)
                else:
                    if group and group["failures"]:
                        text = str(group["failures"])
                    else:
                        text = "not measured" if column != 2 else "device bytes unreported"
                    axis.text(0.98, y, text, transform=axis.get_yaxis_transform(),
                              ha="right", va="center", fontsize=7, color=color)
            if group and group["failures"] and group["wall"]:
                axes[0].text(0.98, index+offset, str(group["failures"]),
                             transform=axes[0].get_yaxis_transform(), ha="right", fontsize=7, color="#AA3322")
    for axis in axes[:2]:
        axis.set_xscale("log")
    axes[2].set_xscale("symlog", linthresh=0.01)
    device_values = [g["device"][2] for g in groups if g["device"]]
    axes[2].set_xlim(left=0, right=max(device_values)*1.2 if device_values and max(device_values) > 0 else 1)
    axes[0].set(yticks=np.arange(len(entries)), yticklabels=labels)
    axes[0].invert_yaxis()
    for axis, label in zip(axes, ("Elapsed solver phases (ms)", "Peak parent RSS (MiB)",
                                 "Attributed device allocation (MiB)")):
        axis.set_xlabel(label)
        axis.grid(True, axis="x", which="both", alpha=0.18)
    fig.suptitle(title+" · native CPU versus CUDA", fontsize=16, y=0.985)
    fig.legend(handles=[Line2D([], [], color=STYLES[b][0], marker=STYLES[b][1], linestyle="none",
               markerfacecolor="white" if b == "native-cuda" else STYLES[b][0],
               markersize=11, label=STYLES[b][2]) for b in ("native", "native-cuda")],
               loc="upper center", bbox_to_anchor=(0.6, 0.94), ncol=2, fontsize=10)
    footer(fig, sources, "Each path uses its fastest complete measured thread configuration for that case/workload; "
           "labels give requested host threads. Whiskers are observed min–max; trial counts and all failed/timeout configurations "
           "remain in the numeric input/manifest. Cold = fresh process/runtime, not a flushed filesystem cache. "
           "Host RSS and attributed device allocations are different quantities; neither is silently added to the other.",
           rect=(0.02, 0.19, 0.99, 0.88))
    fig.savefig(output, dpi=160)
    plt.close(fig)


def parse_cases(value):
    return [(family, int(level)) for family, level in
            (entry.split(":") for entry in value.split(","))]


def numeric_rows(groups):
    rows = []
    for group in groups:
        row = {key: group[key] for key in KEYS}
        row.update(samples=group["samples"], successes=group["successes"],
                   n=int(group["good"][0]["n"]) if group["good"] else "",
                   failures=json.dumps(group["failures"], separators=(",", ":")),
                   source_name=group["source"]["name"], source_sha256=group["source"]["sha256"],
                   source_revision=group["source"]["revision"],
                   active_threads="/".join(map(str, group["active"])))
        for metric in ("wall", "cpu", "rss", "helper", "device"):
            for index, statistic in enumerate(("median", "minimum", "maximum")):
                row[f"{metric}_{statistic}"] = group[metric][index] if group[metric] else ""
        rows.append(row)
    return rows


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=sorted({key for row in rows for key in row}))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path, help="JSON inputs; last configuration wins")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--workload", default="warm", choices=["rhs1", "rhs8", "warm", "refactor", "derivatives"])
    parser.add_argument("--cases", help="Pareto/GPU cases, e.g. grid2:3,grid3:3")
    parser.add_argument("--scaling-cases", help="default: largest measured 2D and 3D cases")
    parser.add_argument("--gpu-cases", help="default: cases with native CUDA measurements")
    parser.add_argument("--variant", default="training", choices=["training", "heldout"])
    parser.add_argument("--title", default="FortSparse measured checkpoint")
    parser.add_argument("--retain-csv", type=Path, help="Optional compact numeric evidence path in source tree")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    groups, sources = load_groups(args.inputs)
    groups = [g for g in groups if g["variant"] == args.variant]
    cases = parse_cases(args.cases) if args.cases else sorted({(g["family"], g["level"]) for g in groups
            if g["backend"] == "native" and g["workload"] == args.workload})[-2:]
    scaling = parse_cases(args.scaling_cases) if args.scaling_cases else [
        (family, max(g["level"] for g in groups if g["family"] == family))
        for family in ("grid2", "grid3") if any(g["family"] == family for g in groups)]
    if not cases:
        parser.error("no matching cases; provide --cases")
    if not scaling:
        scaling = cases[:2]
    gpu_cases = parse_cases(args.gpu_cases) if args.gpu_cases else sorted({
        (g["family"], g["level"]) for g in groups if g["backend"] == "native-cuda"}) or cases[:2]
    plt.rcParams.update({"font.size": 10, "axes.titlepad": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.facecolor": "white"})
    outputs = [args.outdir/name for name in ("pareto.png", "scaling.png", "cpu-vs-gpu.png")]
    pareto_plot(groups, sources, cases, args.workload, args.title, outputs[0])
    scaling_plot(groups, sources, scaling, args.workload, args.title, outputs[1])
    cpu_gpu_plot(groups, sources, gpu_cases, args.title, outputs[2])
    used_cases = set(cases+scaling+gpu_cases)
    rows = numeric_rows([g for g in groups if (g["family"], g["level"]) in used_cases
                         and g["workload"] in (args.workload, "rhs1", "refactor")])
    csv_path = args.outdir/"plot-inputs.csv"
    write_csv(csv_path, rows)
    if args.retain_csv:
        write_csv(args.retain_csv, rows)
    manifest = dict(sources=sources, workload=args.workload, cases=cases, scaling_cases=scaling, gpu_cases=gpu_cases,
                    variant=args.variant, configuration_selection=("Latest binary wins for a complete backend/case/workload "
                        "thread series; identical binaries may merge disjoint thread configurations; latest duplicate wins"),
                    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    matplotlib_version=matplotlib.__version__, numeric_input=str(csv_path),
                    numeric_input_sha256=hashlib.sha256(csv_path.read_bytes()).hexdigest(),
                    outputs={str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in outputs},
                    failed_configurations=[{k: g[k] for k in (*KEYS, "failures")} for g in groups if g["failures"]],
                    notes=["Empirical ranges, not confidence intervals.",
                           "UMFPACK aggregate peak memory unavailable; parent/helper shown separately.",
                           "CPU frontier excludes CUDA device memory and UMFPACK's unknown aggregate.",
                           "Whole-workflow CPU includes setup/checks/cleanup and warmup where applicable."])
    (args.outdir/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print("\n".join(map(str, outputs)))


if __name__ == "__main__":
    main()
