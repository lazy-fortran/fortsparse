#!/usr/bin/env python3
"""Observed CPU/RSS: a completed baseline and an explicitly failed native attempt."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=Path, default=Path(__file__).with_name("plot-data.json"))
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
data = json.loads(args.input.read_text())
baseline = data["baseline"]
native = data["native_failed_attempt"]
assert baseline["exit_status"] == 0 and native["exit_status"] != 0
assert data["successful_native_delta_available"] is False
# Okabe-Ito palette; redundant filled-circle/X encoding survives grayscale.
blue, vermillion = "#0072B2", "#D55E00"
plt.rcParams.update({"font.size": 11, "axes.labelsize": 12, "figure.dpi": 130})
fig, ax = plt.subplots(figsize=(9.2, 7.2))
fig.subplots_adjust(left=0.13, right=0.97, top=0.80, bottom=0.27)
for record, marker, color, label in [
    (baseline, "o", blue, "Legacy baseline: PASS"),
    (native, "X", vermillion, "Native: TIMEOUT (failed attempt)"),
]:
    ax.scatter(record["cpu_s"], record["maxrss_kib"] / 1048576,
               s=120 if marker == "o" else 160, marker=marker,
               color=color, edgecolors="black", linewidths=0.65,
               zorder=3, label=label)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(7, 1600)
ax.set_ylim(0.40, 12)
ax.xaxis.set_major_locator(FixedLocator([10, 100, 1000]))
ax.yaxis.set_major_locator(FixedLocator([0.5, 1, 2, 5, 10]))
formatter = FuncFormatter(lambda value, _: f"{value:g}")
ax.xaxis.set_major_formatter(formatter)
ax.yaxis.set_major_formatter(formatter)
ax.grid(which="major", color="0.85", linewidth=0.7)
ax.spines[["top", "right"]].set_visible(False)
ax.set_xlabel("Observed process CPU (user + system) [s, log scale]")
ax.set_ylabel("Observed process peak RSS [GiB, log scale]")
ax.annotate(f"PASS\n{baseline['cpu_s']:.2f} CPU s; "
            f"{baseline['maxrss_kib']/1048576:.3f} GiB",
            (baseline["cpu_s"], baseline["maxrss_kib"]/1048576),
            xytext=(16, 13), textcoords="offset points", color=blue,
            fontsize=11, ha="left", va="bottom")
ax.annotate(f"TIMEOUT at {data['watchdog_seconds']} s wall\n"
            f"{native['cpu_s']:.2f} CPU s; {native['maxrss_kib']/1048576:.3f} GiB\n"
            "Failed attempt only",
            (native["cpu_s"], native["maxrss_kib"]/1048576),
            xytext=(-16, -24), textcoords="offset points", color=vermillion,
            fontsize=11, ha="right", va="top")
ax.legend(loc="upper left", frameon=False, fontsize=10)
fig.suptitle("KIN6D original coupled Solovev P2: levels 8 / 16 / 32",
             fontsize=15, weight="bold", y=0.96)
fig.text(0.5, 0.91, "One attempt each · 4 requested OMP threads · BLAS 1 · physical CPUs 0–3",
         ha="center", fontsize=11)
fig.text(0.5, 0.865, "Native completion unavailable; successful CPU/RSS delta unavailable",
         ha="center", fontsize=11, color=vermillion)
fig.text(0.13, 0.175,
         "Native X is excluded from feasible/Pareto points. Only if the same execution continued,\n"
         "accumulated CPU and peak RSS could not decrease; this does not predict a repeat.",
         fontsize=10, va="top", linespacing=1.5)
fig.text(0.13, 0.105,
         "Separate correctness run: PASS 797.81 s wall under different load/affinity; not plotted.\n"
         "n = 1; shared host. RSS is per-process, not aggregate workflow memory.\n"
         "Baseline 73891c54 → native fc8b36a / FortSparse 886f7dc.",
         fontsize=9.5, va="top", linespacing=1.45, color="0.25")
args.output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.output, facecolor="white")
