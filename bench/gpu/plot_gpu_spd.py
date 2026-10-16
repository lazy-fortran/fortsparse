"""Render matched host-workflow evidence; generated PNGs stay outside Git."""
import argparse
import hashlib
import json
from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from analyze_gpu_spd import summarize

# Okabe-Ito palette with redundant shapes and line styles.
STYLES = {
    "native": ("#777777", "*", "Native CPU", "-"),
    "native-nd": ("#777777", "P", "Native CPU ND", "--"),
    "native-cuda": ("#D55E00", "^", "Native CUDA", "-"),
    "native-cuda-nd": ("#D55E00", "D", "Native CUDA ND", "--"),
    "cudss": ("#0072B2", "o", "cuDSS", "-"),
    "cholmod": ("#009E73", "s", "CHOLMOD CUDA", "-"),
}


def metric(row, name, scale=1):
    value = row["metrics"].get(name)
    if value is None:
        return None
    return [value[k] * scale for k in ("minimum", "median", "maximum")]


def front(rows, xname, yname):
    known = [r for r in rows if r["complete"] and metric(r, xname) and metric(r, yname)]
    result = []
    for row in known:
        x, y = metric(row, xname)[1], metric(row, yname)[1]
        if not any(metric(other, xname)[1] <= x and metric(other, yname)[1] <= y and
                   (metric(other, xname)[1] < x or metric(other, yname)[1] < y)
                   for other in known):
            result.append(row)
    return sorted(result, key=lambda r: metric(r, xname)[1])


def finish(fig, path, title, note):
    fig.suptitle(title, fontsize=14)
    width = max(65, int(fig.get_figwidth() * 10))
    wrapped = "\n".join(textwrap.fill(line, width) for line in note.splitlines())
    bottom = max(.18, (12 * len(wrapped.splitlines()) + 18) / (72 * fig.get_figheight()))
    fig.text(.5, .015, wrapped, ha="center", va="bottom", fontsize=9)
    fig.tight_layout(rect=(0, bottom, 1, .93))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def completion_note(rows):
    failed = [f"{STYLES[r['requested_backend']][2]} ({r['threads']} cores): "
              f"{r['successes']}/{r['samples']} passed" for r in rows if not r['complete']]
    return "Incomplete configurations: " + "; ".join(failed) if failed else "All configurations completed."


def pareto(rows, path, title, note):
    gpu = [r for r in rows if "cuda" in r["requested_backend"] or
           r["requested_backend"] in ("cudss", "cholmod")]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, field, units, ylabel in zip(axes,
            ("parent_rss_kib", "device_process_sampled_max_bytes"),
            (1 / 1024, 1 / 2**20),
            ("Whole-process peak host RSS (MiB)", "Phase-sampled process VRAM (MiB)")):
        for row in gpu:
            wall, memory = metric(row, "primary_phase_s", 1000), metric(row, field, units)
            if not wall or not memory:
                continue
            color, marker, _, _ = STYLES[row["requested_backend"]]
            ax.errorbar(wall[1], memory[1],
                xerr=[[wall[1] - wall[0]], [wall[2] - wall[1]]],
                yerr=[[memory[1] - memory[0]], [memory[2] - memory[1]]],
                marker=marker, color=color, markersize=8, capsize=2, linestyle="none",
                markerfacecolor=color if field == "parent_rss_kib" and row["complete"] else "white",
                alpha=1 if row["complete"] else .45)
        frontier = front(gpu, "primary_phase_s", field)
        if len(frontier) > 1:
            ax.plot([metric(r, "primary_phase_s", 1000)[1] for r in frontier],
                    [metric(r, field, units)[1] for r in frontier],
                    color="#555555", linestyle=":", linewidth=1, zorder=0)
        ax.set(xlabel="Matched workflow phase wall time (ms)", ylabel=ylabel, xscale="log")
        ax.grid(alpha=.2)
    axes[0].set_title("Median frontier; observed min–max ranges")
    axes[1].set_title("Sampled memory only; total VRAM peak unknown")
    handles = [Line2D([], [], color=STYLES[b][0], marker=STYLES[b][1], linestyle="none",
                     label=STYLES[b][2]) for b in STYLES
               if any(r["requested_backend"] == b for r in gpu)]
    axes[0].legend(handles=handles, fontsize=9)
    failed = sum(not r["complete"] for r in gpu)
    finish(fig, path, title, note + f"\n{failed} incomplete configurations retained; dotted lines are median frontiers, not admission.")


def scaling(rows, path, title, note):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    groups = {}
    for row in rows:
        key = (row["requested_backend"], row["ordering"], row["pool_budget_bytes"])
        groups.setdefault(key, []).append(row)
    for (backend, order, pool), values in sorted(groups.items()):
        color, marker, label, line = STYLES[backend]
        values = sorted([r for r in values if r["complete"]], key=lambda r: r["threads"])
        if not values:
            continue
        if pool > 0:
            label += f" ({pool / 2**20:g} MiB pool)"
        for ax, field in zip(axes, ("primary_phase_s", "workflow_cpu_s")):
            usable = [r for r in values if metric(r, field)]
            if not usable:
                continue
            intervals = [metric(r, field, 1000) for r in usable]
            ax.errorbar([r["threads"] for r in usable], [v[1] for v in intervals],
                yerr=[[v[1] - v[0] for v in intervals], [v[2] - v[1] for v in intervals]],
                color=color, marker=marker, linestyle=line, label=label, capsize=2)
    for ax, ylabel in zip(axes, ("Matched phase wall time (ms)", "Whole-workflow process CPU time (ms)")):
        ax.set(xlabel="Allowed host physical cores (one GPU per GPU solver)", ylabel=ylabel,
               xscale="log", yscale="log")
        ax.set_xticks([1, 2, 4, 8, 16], [1, 2, 4, 8, 16])
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    finish(fig, path, title, note + "\nRequested host budgets shown; native GPU factorization uses one host thread. No multi-GPU scaling claim.")


def cpu_gpu(rows, path, title, note):
    rows = [r for r in rows if r["threads"] == 1]
    fig, ax = plt.subplots(figsize=(9, 5))
    labels = []
    for i, row in enumerate(sorted(rows, key=lambda r: r["requested_backend"])):
        color, marker, label, _ = STYLES[row["requested_backend"]]
        value = metric(row, "primary_phase_s", 1000)
        if row["pool_budget_bytes"] > 0:
            label += f"\n{row['pool_budget_bytes'] / 2**20:g} MiB pool"
        if not row['complete']:
            label += f"\n{row['successes']}/{row['samples']} passed"
        if value:
            ax.errorbar(value[1], i, xerr=[[value[1] - value[0]], [value[2] - value[1]]],
                        marker=marker, color=color, markersize=9, capsize=3,
                        markerfacecolor=color if row['complete'] else 'white')
            ax.text(value[2] * 1.04, i, f"{value[1]:.2f}", va="center", fontsize=9)
        else:
            ax.text(.02, i, "No successful timing", transform=ax.get_yaxis_transform(),
                    va="center", color=color, fontsize=9)
        labels.append(label)
    ax.set_yticks(range(len(labels)), labels)
    ax.set(xlabel="Host input to host solution: matched phase wall time (ms)", xscale="log")
    ax.grid(axis="x", alpha=.2)
    ax.margins(x=.2)
    finish(fig, path, title, note + "\nOne allowed host core; retained-analysis/factor boundaries follow the named workload. Ranges are observed min–max.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--case", required=True, help="FAMILY:LEVEL")
    parser.add_argument("--workload", default="warm")
    parser.add_argument("--variant", default="training")
    args = parser.parse_args()
    document = json.loads(args.input.read_text())
    family, level = args.case.split(":")
    rows = [r for r in summarize(document["trials"]) if r["family"] == family and
            r["level"] == int(level) and r["workload"] == args.workload and
            r["variant"] == args.variant and r["io_mode"] == "host"]
    if not rows:
        parser.error("no matching host-workflow evidence")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    title = f"{family}, level {level}: {args.workload}, {args.variant} coefficients"
    counts = "/".join(str(n) for n in sorted({r["samples"] for r in rows}))
    note = (f"{document.get('measurement_class', 'unknown')} • {counts} trials/configuration • "
            "FP64, backward error ≤ 1e-12\n" + document.get("resource_note", "Resource control unknown"))
    note += "\n" + completion_note(rows)
    for name, render in (("pareto", pareto), ("scaling", scaling), ("cpu-gpu", cpu_gpu)):
        render(rows, args.output_dir / f"{name}-{family}-{level}-{args.workload}.png", title, note)
    (args.output_dir / "manifest.json").write_text(json.dumps(dict(
        input=str(args.input.resolve()), sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
        case=args.case, workload=args.workload, variant=args.variant,
        semantics="Median frontiers are diagnostic; sampled VRAM is not total peak memory."), indent=2) + "\n")


if __name__ == "__main__":
    main()
