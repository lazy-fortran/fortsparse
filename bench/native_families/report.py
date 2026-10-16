#!/usr/bin/env python3
"""Summarize the initial 2026-10-16 receipts; render diagnostic PNGs externally.

Interpretation and CPU plot labels describe that checkpoint, not later runs.
"""
import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

KEYS = ('family', 'backend', 'level', 'threads', 'workload', 'variant', 'operation')
METRICS = ('process_wall_s', 'total_s', 'factor_s', 'solve_s', 'refactor_s',
           'parent_cpu_s', 'parent_main_cpu_s', 'helper_cpu_s', 'process_tree_cpu_s',
           'parent_rss_kib', 'helper_rss_kib', 'active_threads',
           'native_factor_serial', 'factor_nnz', 'host_peak_bytes',
           'backward_error', 'forward_error')
BACKENDS = ('native', 'native-amd', 'superlu', 'umfpack', 'mumps', 'intel-pardiso')
FAMILIES = ('neo2-mini5', 'drift-real', 'drift-complex', 'hpd-chain', 'ldlt-chain')
# Okabe-Ito palette; marker shapes redundantly identify backends.
COLORS = ('#000000', '#CC79A7', '#0072B2', '#E69F00', '#009E73', '#D55E00')
MARKERS = ('o', 'v', 's', '^', 'D', 'P')


def summarize(attempts):
    groups = collections.defaultdict(list)
    for attempt in attempts:
        groups[tuple(attempt[key] for key in KEYS)].append(attempt)
    rows = []
    for key, group in sorted(groups.items()):
        ok = [a for a in group if a['status'] == 'ok' and a.get('returncode', 0) == 0]
        row = dict(zip(KEYS, key), attempts=len(group), ok=len(ok),
                   timeouts=sum(a['status'] == 'timeout' for a in group),
                   failures=len(group)-len(ok)-sum(a['status'] == 'timeout' for a in group))
        for metric in METRICS:
            values = [a.get('record', {}).get(metric, -1) for a in ok]
            values = [v for v in values if isinstance(v, (int, float))
                      and math.isfinite(v) and v >= 0]
            # Wall measurements include failed/killed launches, kept separately.
            for suffix, value in zip(('count', 'min', 'median', 'max'),
                (len(values), min(values) if values else '',
                 statistics.median(values) if values else '', max(values) if values else '')):
                row[f'{metric}_{suffix}'] = value
        killed = [a['process_wall_s'] for a in group if a['status'] == 'timeout']
        row['timeout_observed_wall_min_s'] = min(killed) if killed else ''
        row['timeout_observed_wall_max_s'] = max(killed) if killed else ''
        rows.append(row)
    return rows


def render(rows, out, timeout):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.titlesize': 11, 'axes.labelsize': 10})
    available = [b for b in BACKENDS if any(r['backend'] == b for r in rows)]

    def figure(name, title, panels, variant='training', logarithmic=True):
        fig, axes = plt.subplots(2, 5, figsize=(18, 7.2), squeeze=False)
        for col, family in enumerate(FAMILIES):
            for rowindex, (workload, metric, ylabel) in enumerate(panels):
                ax = axes[rowindex, col]
                for backend in available:
                    index = BACKENDS.index(backend)
                    subset = [r for r in rows if r['family'] == family
                              and r['backend'] == backend and r['variant'] == variant
                              and r['workload'] == workload and r['operation'] == 'N']
                    subset.sort(key=lambda r: r['threads'])
                    points = [r for r in subset if r[metric+'_count']
                              and (not logarithmic or r[metric+'_median'] > 0)]
                    if points:
                        x = [r['threads'] for r in points]
                        y = [r[metric+'_median'] for r in points]
                        lo = [v-r[metric+'_min'] for v, r in zip(y, points)]
                        hi = [r[metric+'_max']-v for v, r in zip(y, points)]
                        ax.errorbar(x, y, yerr=[lo, hi], color=COLORS[index],
                                    marker=MARKERS[index], capsize=2, lw=1, ms=4,
                                    label=backend)
                    if metric == 'process_wall_s':
                        killed = [r for r in subset if r['timeouts']]
                        if killed:
                            ax.scatter([r['threads'] for r in killed],
                                       [timeout]*len(killed), color=COLORS[index],
                                       marker='x', s=36)
                if logarithmic:
                    ax.set_yscale('log')
                ax.set_xscale('log', base=2)
                ax.set_xticks([1, 2, 4, 8, 16], ['1', '2', '4', '8', '16'])
                ax.grid(alpha=0.25)
                ax.set_title(family if rowindex == 0 else workload)
                ax.set_xlabel('Requested thread budget')
                if col == 0:
                    ax.set_ylabel(ylabel)
        handles, labels = [], []
        for ax in axes.flat:
            h, l = ax.get_legend_handles_labels()
            for handle, label in zip(h, l):
                if label not in labels:
                    handles.append(handle)
                    labels.append(label)
        fig.legend(handles, labels, ncol=len(available), loc='lower center')
        fig.suptitle(title)
        fig.tight_layout(rect=(0, 0.06, 1, 0.94))
        fig.savefig(out/name, dpi=150)
        plt.close(fig)

    for variant in ('training', 'heldout'):
        figure(f'wall-{variant}.png',
               f'{variant}: subprocess wall median and full trial range; × = timeout ≥ {timeout:g}s',
               [('rhs8', 'process_wall_s', 'Cold factor + RHS8 workflow (s)'),
                ('refactor', 'process_wall_s', 'Prime + refactor workflow (s)')], variant)
    figure('phases-training.png', 'Training: timed solver calls, separate from subprocess wall',
           [('rhs8', 'factor_s', 'Cold factor call (s)'),
            ('refactor', 'refactor_s', 'Refactor call (s)')])
    figure('cpu-training.png',
           'Initial baseline CPU: parent main-entry interval and reaped-helper CPU (not whole workflow CPU)',
           [('rhs8', 'parent_cpu_s', 'Parent main-entry CPU interval (s)'),
            ('rhs8', 'helper_cpu_s', 'Reaped-helper lifetime CPU (s)')], logarithmic=False)
    figure('rss-training.png',
           'Training RHS8: separate memory diagnostics; UMFPACK aggregate RSS remains unknown',
           [('rhs8', 'parent_rss_kib', 'Parent VmHWM peak (KiB)'),
            ('rhs8', 'helper_rss_kib', 'Helper launch-lifecycle ru_maxrss (KiB)')], logarithmic=False)


def report(attempts, rows, provenance, source_hash, cpu_scope):
    counts = collections.Counter(a['status'] for a in attempts)
    timeouts = collections.Counter((a['family'], a['workload']) for a in attempts
                                  if a['status'] == 'timeout')
    loads = [a['background_load_before'][0] for a in attempts
             if 'background_load_before' in a]
    text = [f'# Native sparse-family initial baseline', '',
            f'Source revision: `{provenance.get("source_revision", "unknown")}`. '
            f'Raw attempts SHA256: `{source_hash}`.', '',
            f'{len(attempts)} attempts: '+', '.join(f'{n} {s}' for s, n in counts.items())+'.',
            f'CPU interpretation: **{cpu_scope}**. Prior raw fields are preserved.', '',
            'The rows below show training, normal operation, RHS8, one requested thread. '
            'Wall is complete subprocess elapsed time; factor is a timed library call. '
            'Ranges span all successful trials. Timeout rows have no factor/RSS sample.', '',
            '| Family | Backend | OK/all | Wall median [min,max] ms | Factor median ms | Parent peak MiB |',
            '|---|---|---:|---:|---:|---:|']
    for family in FAMILIES:
        for r in rows:
            if (r['family'] != family or r['threads'] != 1 or r['variant'] != 'training'
                    or r['workload'] != 'rhs8' or r['operation'] != 'N'):
                continue
            def fmt(metric, scale=1):
                return f'{r[metric+"_median"]*scale:.3f}' if r[metric+'_count'] else 'unavailable'
            wall = (f'{fmt("process_wall_s", 1000)} '
                    f'[{r["process_wall_s_min"]*1000:.3f},{r["process_wall_s_max"]*1000:.3f}]'
                    if r['ok'] else 'timeout')
            text.append(f'| {family} | {r["backend"]} | {r["ok"]}/{r["attempts"]} | '
                        f'{wall} | {fmt("factor_s", 1000)} | {fmt("parent_rss_kib", 1/1024)} |')
    text += ['', 'Timeout counts (all thread budgets, variants and trials):', '']
    text += [f'- {family}, {workload}: {count}.' for (family, workload), count in sorted(timeouts.items())]
    text += ['', 'Interpretation and limits:', '',
        '- Natural-order native LU is a clear bottleneck on the retained 2D drift cases. '
        'Complex drift times out for both workloads; real drift times out for refactor.',
        '- The HPD chain has a small native factor-call cost; total process cost is much '
        'larger. The chain and mini5 fixtures do not establish broad production leadership.',
        '- Native records report one active factor thread at every requested budget. '
        'External active teams are unknown (-1), not evidence that they used the requested team.',
        '- Parent VmHWM and helper RUSAGE_CHILDREN ru_maxrss remain separate. Helper RSS '
        'can include inherited/pre-exec highwater. Never sum them or use their maximum as '
        'a validated aggregate lower bound. UMFPACK aggregate memory is unknown; its parent '
        'VmHWM alone is a safe lower bound.',
        '- CPU in the initial receipts is the parent main-entry interval plus separately '
        'recorded helper lifetime CPU. These do not measure complete workflow CPU. '
        'New wait4 instrumentation is a separate measurement checkpoint.',
        '- Successful-only timing ranges do not replace timeout evidence. A timeout marks '
        'a censored observation at the configured deadline, not a measured factor time.',
        '- No global Pareto/admission claim follows: this is a bounded fixture corpus, '
        'several timings are censored, true helper aggregate memory is unknown, and the '
        'shared host has background load.']
    if loads:
        text += [f'- Observed pre-launch one-minute host load: {min(loads):.2f}–{max(loads):.2f}. '
                 'Core affinity reduces placement variation; it does not eliminate shared-host interference.']
    return '\n'.join(text)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='immutable measurement directory')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plots', action='store_true')
    parser.add_argument('--cpu-scope', required=True,
                        help='explicit interpretation; never silently relabel older receipts')
    args = parser.parse_args()
    source = (args.input/'attempts.jsonl').read_bytes()
    attempts = [json.loads(line) for line in source.splitlines() if line]
    provenance = json.loads((args.input/'provenance.json').read_text())
    rows = summarize(attempts)
    for row in rows:
        row["parent_cpu_scope"] = args.cpu_scope
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output/'summary.csv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (args.output/'report.md').write_text(report(attempts, rows, provenance,
        hashlib.sha256(source).hexdigest(), args.cpu_scope))
    (args.output/'recipe.json').write_text(json.dumps(dict(
        raw_path=str(args.input.resolve()), raw_sha256=hashlib.sha256(source).hexdigest(),
        cpu_scope=args.cpu_scope, plot_filter='N; per-variant wall; training phases/CPU/RSS',
        aggregation='median and full min/max range across successful trials; all timeouts retained',
        memory='separate parent VmHWM and helper launch-lifecycle peak; no sum/maximum aggregate',
        timeout_s=provenance.get('timeout_s', 0.9)), indent=2)+'\n')
    if args.plots:
        render(rows, args.output, provenance.get('timeout_s', 0.9))
    print(json.dumps(dict(attempts=len(attempts), groups=len(rows),
                          status=collections.Counter(a['status'] for a in attempts),
                          output=str(args.output))))


if __name__ == '__main__':
    main()
