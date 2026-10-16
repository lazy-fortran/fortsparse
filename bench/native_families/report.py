#!/usr/bin/env python3
"""Summarize immutable family checkpoints without conflating CPU/memory scopes."""
import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

KEYS = ('family', 'backend', 'level', 'threads', 'workload', 'variant', 'operation')
METRICS = ('n', 'nnz', 'nrhs', 'process_wall_s', 'total_s', 'factor_s', 'solve_s', 'refactor_s',
           'parent_cpu_s', 'parent_main_cpu_s', 'helper_cpu_s', 'process_tree_cpu_s',
           'parent_rss_kib', 'helper_rss_kib', 'analysis_s', 'numeric_s',
           'prime_factor_s', 'active_threads',
           'native_factor_serial', 'factor_nnz', 'host_peak_bytes',
           'backward_error', 'forward_error')
BACKENDS = ('native', 'native-amd', 'superlu', 'umfpack', 'mumps', 'intel-pardiso', 'umfpack-direct',
            'native-amd-low', 'native-amd-capped', 'umfpack-direct-noir')
FAMILIES = ('neo2-mini5', 'drift-real', 'drift-complex', 'hpd-chain', 'ldlt-chain')
# Existing seven styles retained; memory aliases add distinct colors and markers.
COLORS = ('#000000', '#CC79A7', '#0072B2', '#E69F00', '#009E73', '#D55E00', '#56B4E9', '#882255', '#999933', '#332288')
MARKERS = ('o', 'v', 's', '^', 'D', 'P', 'X', '*', 'h', '+')


def summarize(attempts):
    groups = collections.defaultdict(list)
    for attempt in attempts:
        groups[tuple(attempt[key] for key in KEYS)].append(attempt)
    rows = []
    for key, group in sorted(groups.items()):
        ok = [a for a in group if a['status'] == 'ok' and a.get('returncode', 0) == 0]
        configurations = {(a.get('record', {}).get('memory_policy', -1),
                           a.get('record', {}).get('host_memory_limit_bytes', -1)) for a in group
                          if a.get('record')}
        if len(configurations)>1:
            raise ValueError("mixed memory policy/cap in one benchmark group")
        policy,cap = next(iter(configurations),(-1,-1))
        row = dict(zip(KEYS, key), memory_policy=policy,host_memory_limit_bytes=cap,
                   attempts=len(group), ok=len(ok),
                   timeouts=sum(a['status'] == 'timeout' for a in group),
                   failures=len(group)-len(ok)-sum(a['status'] == 'timeout' for a in group),
                   reused_true=sum(a.get('record', {}).get('reused') == 'T' for a in ok))
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


def available_backends(rows):
    unknown = {r['backend'] for r in rows}-set(BACKENDS)
    if unknown:
        raise ValueError("plot backend has no declared style: "+", ".join(sorted(unknown)))
    return [b for b in BACKENDS if any(r['backend']==b for r in rows)]


def render(rows, out, timeout, cpu_model, front_note='no frontier admission'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.titlesize': 11, 'axes.labelsize': 10})
    available = available_backends(rows)
    families = [f for f in FAMILIES if any(r['family'] == f for r in rows)]
    native_teams = sorted({int(r['active_threads_'+bound]) for r in rows
                           if r['backend'].startswith('native') and
                           r.get('active_threads_count', 0) for bound in ('min', 'max')})
    team_label = '/'.join(map(str, native_teams)) or 'unknown'
    external_team_label = ('unknown' if not any(r.get('active_threads_count', 0)
        for r in rows if not r['backend'].startswith('native')) else 'see numeric CSV')

    def figure(name, title, panels, variant='training', logarithmic=True):
        fig, axes = plt.subplots(2, len(families),
                                 figsize=(max(9, 3.6*len(families)), 7.2), squeeze=False)
        for col, family in enumerate(families):
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
                sizes = sorted({int(r['n_median']) for r in rows
                                if r['family'] == family and r.get('n_count', 0)})
                size_label = ','.join(map(str, sizes)) or 'unknown'
                ax.set_title(f'{family} (n={size_label})' if rowindex == 0 else workload)
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
        fig.tight_layout(rect=(0, 0.06, 1, 0.94), w_pad=2)
        fig.savefig(out/name, dpi=150)
        plt.close(fig)

    censored = (f'× = timeout ≥ {timeout:g}s' if any(r['timeouts'] for r in rows)
                else 'all trials completed')
    for variant in ('training', 'heldout'):
        figure(f'wall-{variant}.png',
               f'{variant}: subprocess wall median and full trial range; {censored}',
               [('rhs8', 'process_wall_s', 'Cold factor + RHS8 workflow (s)'),
                ('refactor', 'process_wall_s', 'Prime + refactor workflow (s)')], variant)
    figure('phases-training.png', 'Training: timed solver calls, separate from subprocess wall',
           [('rhs8', 'factor_s', 'Cold factor call (s)'),
            ('refactor', 'refactor_s', 'Refactor call (s)')])
    figure('solve-training.png', 'Training: solve calls, separate from factor/refactor and subprocess wall',
           [('rhs8', 'solve_s', 'Solve call after cold factor (s)'),
            ('refactor', 'solve_s', 'Solve call after refactor (s)')])
    if cpu_model == 'initial':
        figure('cpu-training.png',
               'Initial CPU: parent main-entry interval and helper CPU; workflow CPU unavailable',
               [('rhs8', 'parent_cpu_s', 'Parent main-entry CPU interval (s)'),
                ('rhs8', 'helper_cpu_s', 'Reaped-helper lifetime CPU (s)')], logarithmic=False)
    else:
        figure('cpu-training.png',
               'Training: completed workflow CPU (wait4, including reaped helpers)',
               [('rhs8', 'process_tree_cpu_s', 'Cold factor + RHS8 workflow CPU (s)'),
                ('refactor', 'process_tree_cpu_s', 'Prime + refactor workflow CPU (s)')],
               logarithmic=False)
    figure('rss-training.png',
           'Training RHS8: separate memory diagnostics; UMFPACK aggregate RSS remains unknown',
           [('rhs8', 'parent_rss_kib', 'Parent VmHWM peak (KiB)'),
            ('rhs8', 'helper_rss_kib', 'Helper launch-lifecycle ru_maxrss (KiB)')], logarithmic=False)

    if cpu_model == 'completed':
        from matplotlib.lines import Line2D
        for variant in ('training', 'heldout'):
            fig, axes = plt.subplots(len(families), 2,
                                     figsize=(12, 4*len(families)), squeeze=False)
            for i, family in enumerate(families):
                for j, workload in enumerate(('rhs8', 'refactor')):
                    ax = axes[i, j]
                    for backend in available:
                        if backend == 'umfpack':
                            continue  # Its parent peak is not aggregate workflow RSS.
                        index = BACKENDS.index(backend)
                        points = [r for r in rows if r['family'] == family and
                                  r['backend'] == backend and r['variant'] == variant and
                                  r['workload'] == workload and r['operation'] == 'N' and
                                  r['parent_rss_kib_count'] and r['process_wall_s_count']]
                        points.sort(key=lambda r: r['threads'])
                        for r in points:
                            x = r['process_wall_s_median']*1000
                            y = r['parent_rss_kib_median']/1024
                            ax.errorbar(x, y,
                                xerr=[[x-r['process_wall_s_min']*1000],
                                      [r['process_wall_s_max']*1000-x]],
                                yerr=[[y-r['parent_rss_kib_min']/1024],
                                      [r['parent_rss_kib_max']/1024-y]],
                                color=COLORS[index], marker=MARKERS[index], ms=5,
                                capsize=2, alpha=0.8)
                    sizes = sorted({int(r['n_median']) for r in rows
                                    if r['family'] == family and r.get('n_count', 0)})
                    ax.set_title(f'{family}, n={",".join(map(str, sizes))}, {workload}')
                    ax.set_xlabel('Complete subprocess wall (ms)')
                    ax.set_ylabel('Single-process VmHWM peak (MiB)')
                    ax.grid(alpha=0.25)
            handles = [Line2D([], [], color=COLORS[BACKENDS.index(b)],
                             marker=MARKERS[BACKENDS.index(b)], ls='', label=b)
                       for b in available if b != 'umfpack']
            fig.legend(handles=handles, loc='lower center', ncol=len(handles))
            fig.suptitle(f'{variant}: time/RSS diagnostic; {front_note}\n'
                         f'Native factor team={team_label}; external teams {external_team_label}; '
                         'UMFPACK IPC aggregate RSS unknown: no point shown')
            fig.tight_layout(rect=(0, 0.05, 1, 0.93))
            fig.savefig(out/f'time-rss-{variant}.png', dpi=150)
            plt.close(fig)


def report(attempts, rows, provenance, source_hash, cpu_scope, cpu_model):
    counts = collections.Counter(a['status'] for a in attempts)
    timeouts = collections.Counter((a['family'], a['workload']) for a in attempts
                                  if a['status'] == 'timeout')
    loads = [a['background_load_before'][0] for a in attempts
             if 'background_load_before' in a]
    text = [f'# Native sparse-family {cpu_model} checkpoint', '',
            f'Source revision: `{provenance.get("source_revision", "unknown")}`. '
            f'Raw attempts SHA256: `{source_hash}`.', '',
            f'{len(attempts)} attempts: '+', '.join(f'{n} {s}' for s, n in counts.items())+'.',
            f'CPU interpretation: **{cpu_scope}**. Prior raw fields are preserved.', '',
            'The rows below show training, normal operation, RHS8, one requested thread. '
            'Wall is complete subprocess elapsed time; factor is a timed library call. '
            'Ranges span all successful trials. Timeout rows have no factor/RSS sample.', '',
            '| Family | Backend | OK/all | Wall median [min,max] ms | Factor median ms | Workflow CPU median ms | Parent peak MiB |',
            '|---|---|---:|---:|---:|---:|---:|']
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
                        f'{wall} | {fmt("factor_s", 1000)} | '
                        f'{fmt("process_tree_cpu_s", 1000)} | '
                        f'{fmt("parent_rss_kib", 1/1024)} |')
    text += ['', 'Timeout counts (all thread budgets, variants and trials):', '']
    text += [f'- {family}, {workload}: {count}.' for (family, workload), count in sorted(timeouts.items())]
    if not timeouts:
        text += ['- None.']
    text += ['', 'Interpretation and limits:', '']
    if cpu_model == 'initial':
        text += ['- Natural-order native LU is a clear bottleneck on these 2D drift cases. '
                 'Complex drift times out for both workloads; real drift refactor times out.',
                 '- HPD chain native factor-call cost is small relative to subprocess cost. '
                 'Chains and mini5 do not establish broad production leadership.',
                 '- Initial parent CPU is the main-entry interval; complete workflow CPU '
                 'is unavailable. New wait4 measurements form a separate checkpoint.']
    else:
        text += ['- AMD is an explicit native ordering option; default ordering is unchanged. '
                 'Native AMD still takes more total wall time than the best declared external '
                 'comparator. This checkpoint does not admit it to the performance frontier.',
                 '- UMFPACK here uses the helper/IPC path. A standalone direct UMFPACK '
                 'comparator is still required before judging best-library performance.',
                 '- Completed workflow CPU uses wait4 user+system usage, including reaped '
                 'helper CPU and termination. Helper CPU must not be added again. Parent '
                 'pre-CSV cumulative CPU and main interval remain separate in the CSV.',
                 '- Uncapped native AMD attributed peak bytes are unavailable (-1), never zero. '
                 'LOW/capped aliases retain their explicit policy/cap in CSV; owned payload '
                 'is separate from RSS. '
                 'Parent VmHWM is a process RSS measurement, not attributed factor memory.']
    text += [
        '- Native records report one active factor thread at every requested budget. '
        'External active teams are unknown (-1), not evidence that they used the requested team.',
        '- Parent VmHWM and helper RUSAGE_CHILDREN ru_maxrss remain separate. Helper RSS '
        'can include inherited/pre-exec highwater. Never sum them or use their maximum as '
        'a validated aggregate lower bound. UMFPACK aggregate memory is unknown; its parent '
        'VmHWM alone is a safe lower bound.',
        '- Successful-only timing ranges do not replace timeout evidence. A timeout marks '
        'a censored observation at the configured deadline, not a measured factor time.',
        '- No global Pareto/admission claim follows: this is a bounded fixture corpus, '
        'the initial checkpoint retains 150 censored timings, true helper aggregate '
        'memory is unknown, and the shared host has background load.']
    if loads:
        text += [f'- Observed pre-launch one-minute host load: {min(loads):.2f}–{max(loads):.2f}. '
                 'Core affinity reduces placement variation; it does not eliminate shared-host interference.']
    return '\n'.join(text)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='immutable measurement directory')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plots', action='store_true')
    parser.add_argument('--cpu-model', choices=('initial', 'completed'), default='initial')
    parser.add_argument('--cpu-scope', required=True,
                        help='explicit interpretation; never silently relabel older receipts')
    args = parser.parse_args()
    source = (args.input/'attempts.jsonl').read_bytes()
    attempts = [json.loads(line) for line in source.splitlines() if line]
    provenance = json.loads((args.input/'provenance.json').read_text())
    if args.cpu_model == 'completed' and any(a.get('status') == 'ok' and
            a.get('process_tree_cpu_s', -1) < 0 for a in attempts):
        parser.error('completed CPU requested but successful receipts lack wait4 CPU')
    rows = summarize(attempts)
    for row in rows:
        row["parent_cpu_scope"] = args.cpu_scope
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output/'summary.csv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (args.output/'report.md').write_text(report(attempts, rows, provenance,
        hashlib.sha256(source).hexdigest(), args.cpu_scope, args.cpu_model))
    (args.output/'recipe.json').write_text(json.dumps(dict(
        raw_path=str(args.input.resolve()), raw_sha256=hashlib.sha256(source).hexdigest(),
        cpu_scope=args.cpu_scope, cpu_model=args.cpu_model, plot_filter='N; per-variant wall; training phases/CPU/RSS',
        aggregation='median and full min/max range across successful trials; all timeouts retained',
        memory='separate parent VmHWM and helper launch-lifecycle peak; no sum/maximum aggregate',
        timeout_s=provenance.get('timeout_s', 0.9)), indent=2)+'\n')
    if args.plots:
        render(rows, args.output, provenance.get('timeout_s', 0.9), args.cpu_model)
    print(json.dumps(dict(attempts=len(attempts), groups=len(rows),
                          status=collections.Counter(a['status'] for a in attempts),
                          output=str(args.output))))


if __name__ == '__main__':
    main()
