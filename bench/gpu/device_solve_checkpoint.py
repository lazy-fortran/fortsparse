"""Compact evidence and plots for the bounded GPU device-solve checkpoint."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import textwrap

from analyze_gpu_spd import summarize


def export(rows, path):
    fields = sorted({name for row in rows for name in row["metrics"]})
    keys = ["family", "level", "workload", "variant", "requested_backend", "series",
            "executable_sha256", "threads", "ordering", "pool_budget_bytes",
            "samples", "successes", "complete", "offload_verified", "statuses"]
    columns = keys + [f"{name}_{stat}" for name in fields
                      for stat in ("minimum", "median", "maximum")]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, columns, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            record = {key: row[key] for key in keys}
            record["statuses"] = json.dumps(record["statuses"], sort_keys=True)
            for name in fields:
                for stat in ("minimum", "median", "maximum"):
                    record[f"{name}_{stat}"] = (row["metrics"][name] or {}).get(stat)
            writer.writerow(record)


def render(rows, cases, directory):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from plot_gpu_spd import STYLES, series_label, metric

    def identity(row):
        return (row["requested_backend"], row["series"], row["executable_sha256"])

    identities = {identity(row): row for row in rows}
    workloads = ("warm", "refactor")
    footer = ("Shared GPU with running model service; compiler/profiler jobs paused. "
              "Diagnostic only; no performance admission.\n"
              "Five fresh-process trials/configuration, observed min–max; FP64 backward error ≤1e-12. "
              "Timeouts retain the 1s whole-invocation limit; successful-only ranges are marked below.")

    def canvas(title):
        fig, axes = plt.subplots(2, 2, figsize=(14, 9))
        fig.suptitle(title, fontsize=15, y=.985)
        for ax, (family, level) in zip(axes.flat, cases):
            ax.set_title(f"{family}, level {level}")
            ax.grid(alpha=.2)
        return fig, axes

    def finish(fig, filename, selected, extra):
        keys = sorted({identity(row) for row in selected})
        handles = []
        for key in keys:
            backend, series, _ = key
            color, marker, _, _ = STYLES[backend]
            handles.append(Line2D([], [], color=color, marker=marker, linestyle="none",
                markerfacecolor="white" if series == "reference" else color,
                label=series_label(identities[key])))
        legend = fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5,.95),
                            ncol=3, fontsize=8)
        legend.set_zorder(20)
        caption = "\n".join(textwrap.fill(line,145) for line in (footer+"\n"+extra).splitlines())
        fig.text(.5,.012,caption,ha="center",va="bottom",fontsize=9)
        fig.tight_layout(rect=(.01,.16,.99,.87))
        fig.savefig(directory/filename,dpi=160)
        plt.close(fig)

    def completion(ax, selected):
        failed = [f"{STYLES[r['requested_backend']][2]} {r['series']} {r['workload']}: "
                  f"{r['successes']}/{r['samples']} passed" for r in selected if not r["complete"]]
        if failed:
            ax.text(.02,.02,"\n".join(failed),transform=ax.transAxes,fontsize=7,
                    bbox=dict(facecolor="white",edgecolor="none",alpha=.85),va="bottom")

    gpu = [r for r in rows if r["requested_device"] >= 0]
    fig, axes = canvas("GPU libraries: warm factor+solve and refactor+solve")
    gpu_keys = sorted({identity(row) for row in gpu})
    for ax, case in zip(axes.flat,cases):
        selected = [r for r in gpu if (r["family"],r["level"]) == tuple(case)]
        for row in selected:
            value = metric(row,"primary_phase_s",1000)
            if value is None:
                continue
            offset = (gpu_keys.index(identity(row))-(len(gpu_keys)-1)/2)*.09
            color,marker,_,_ = STYLES[row["requested_backend"]]
            ax.errorbar(workloads.index(row["workload"])+offset,value[1],
                yerr=[[value[1]-value[0]],[value[2]-value[1]]],color=color,
                marker=marker,linestyle="none",capsize=2,
                markerfacecolor="white" if row["series"] == "reference" or
                    not row["complete"] else color)
        ax.set_xticks([0,1],["Warm new factor+solve","Retained-analysis refactor+solve"])
        ax.set(yscale="log",ylabel="Host input to host answer phase time (ms)")
        completion(ax,selected)
    finish(fig,"gpu-libraries-final.png",gpu,
           "Reference = immutable prior DAG factorization with CPU solve. "
           "GPU-enabled comparators retain their reported offload evidence; no counter flag is invented.")

    candidate = [r for r in rows if r["series"] == "candidate"]
    for filename,field,scale,ylabel,title in (
        ("cpu-gpu-final.png","workflow_cpu_s",1000,"Whole-workflow process CPU time (ms)",
         "CPU and GPU: matched phase wall time versus whole-workflow CPU"),
        ("host-memory-final.png","parent_rss_kib",1/1024,"Whole-process peak host RSS (MiB)",
         "CPU and GPU: matched phase wall time versus host memory")):
        fig,axes = canvas(title)
        for ax,case in zip(axes.flat,cases):
            selected = [r for r in candidate if (r["family"],r["level"]) == tuple(case)]
            for row in selected:
                x,y = metric(row,"primary_phase_s",1000),metric(row,field,scale)
                if x is None or y is None:
                    continue
                color,marker,_,_ = STYLES[row["requested_backend"]]
                ax.errorbar(x[1],y[1],xerr=[[x[1]-x[0]],[x[2]-x[1]]],
                    yerr=[[y[1]-y[0]],[y[2]-y[1]]],color=color,marker=marker,
                    linestyle="none",capsize=2,markerfacecolor=color if
                        row["workload"] == "refactor" and row["complete"] else "white")
            ax.set(xscale="log",yscale="log",
                   xlabel="Matched factor/refactor+solve phase wall time (ms)",ylabel=ylabel)
            completion(ax,selected)
        finish(fig,filename,candidate,"Hollow = warm new factor; filled = retained-analysis refactor. "
            "One allowed physical host core for every configuration.\n"+
            ("Whole-workflow CPU includes warmup, setup, checks and cleanup; it is not phase CPU time."
             if field == "workflow_cpu_s" else
             "Total device peak is unknown. Phase-sampled VRAM and attributed allocations remain separate in the CSV."))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input",type=Path)
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--offload-provenance",type=Path,
                        help="Revalidate exact retained Nsight proof before exporting flags")
    args = parser.parse_args()
    document = json.loads(args.input.read_text())
    if args.offload_provenance:
        from offload_provenance import attach
        document = attach(document,args.offload_provenance)
    rows = summarize(document["trials"],validated_trace=bool(args.offload_provenance))
    if (len(document["corpus"]) != 4 or not rows or
            any(r["threads"] != 1 or r["io_mode"] != "host" or
                r["variant"] != "training" or r["workload"] not in ("warm","refactor")
                for r in rows)):
        parser.error("this fixed checkpoint recipe requires four training cases, host I/O, "
                     "one host core, and warm/refactor workloads")
    args.output_dir.mkdir(parents=True,exist_ok=True)
    export(rows,args.output_dir/"summary.csv")
    render(rows,document["corpus"],args.output_dir)
    manifest = dict(input=str(args.input.resolve()),
        input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        measurement_class=document["measurement_class"],status_counts=document["status_counts"],
        offload_provenance=document.get("offload_provenance"),
        semantics="Min/median/max of successful trials, with all status counts retained. "
                  "No pooling across executable SHA or reference/candidate. No scaling was sampled. "
                  "Whole-workflow CPU includes setup/check/free/warmup. Unknown device total peak stays unknown.")
    (args.output_dir/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")


if __name__ == "__main__":
    main()
