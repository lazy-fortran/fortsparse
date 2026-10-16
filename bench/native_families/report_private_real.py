#!/usr/bin/env python3
"""Publish aggregate-only private-fixture metrics; never copy numeric exports."""
import argparse
import collections
import csv
import hashlib
import json
import itertools
from pathlib import Path
import statistics

from report import KEYS, METRICS, BACKENDS, COLORS, MARKERS, summarize
from report_front import PEERS, frontier

AUDIT_METRICS = ('audit_reader_wall_s', 'child_plus_audit_wall_s',
                 'audit_normwise_backward_max', 'audit_row_scaled_backward_max',
                 'empirical_forward_relative_inf_max')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalize(wrapper, family):
    """Keep scalar metrics only; wrapper/audit failure overrides child success."""
    child = wrapper.get('process_receipt', {})
    source = child.get('record', {})
    audit = wrapper.get('independent_audit', {})
    passed = (wrapper['status'] == 'pass' and child.get('status') == 'ok'
              and child.get('returncode') == 0 and audit.get('status') == 'pass')
    point = dict(backend=wrapper['backend'], family=family, level=source.get('level', 1),
                 threads=wrapper['threads'], workload=source.get('workload', 'rhs4'),
                 variant=source.get('variant', 'training'), operation=source.get('operation', 'N'))
    record = {k: source[k] for k in METRICS if k in source}
    record.update(point)
    normalized = dict(point, status='ok' if passed else
                      'timeout' if child.get('status') == 'timeout' else 'failed',
                      returncode=child.get('returncode', -1), record=record,
                      process_wall_s=child.get('process_wall_s', -1), trial=wrapper['trial'])
    safe = dict(point, trial=wrapper['trial'], wrapper_status=wrapper['status'],
                child_status=child.get('status', 'unavailable'),
                audit_status=audit.get('status', 'unavailable'),
                audit_rhs_count=len(audit.get('rhs', [])))
    safe.update({k: record.get(k, -1) for k in METRICS})
    safe['audit_reader_wall_s'] = wrapper.get('independent_audit_including_reader_wall_s', -1)
    safe['child_plus_audit_wall_s'] = wrapper.get('solver_child_plus_audit_wall_s', -1)
    rhs = audit.get('rhs', [])
    for metric, key in (('audit_normwise_backward_max', 'normwise_backward'),
                        ('audit_row_scaled_backward_max', 'row_scaled_backward')):
        safe[metric] = max((r[key] for r in rhs if key in r), default=-1)
    safe['empirical_forward_relative_inf_max'] = max((value for r in rhs
        for value in r.get('forward_relative_inf', {}).values()), default=-1)
    return normalized, safe


def audit_summary(safe):
    groups = collections.defaultdict(list)
    for attempt in safe:
        groups[tuple(attempt[k] for k in KEYS)].append(attempt)
    rows = []
    for key, group in sorted(groups.items()):
        row = dict(zip(KEYS, key), attempts=len(group))
        for field in ('wrapper_status', 'child_status', 'audit_status'):
            row[field+'_counts'] = json.dumps(collections.Counter(a[field] for a in group), sort_keys=True)
        row['audit_rhs_count_min'] = min(a['audit_rhs_count'] for a in group)
        row['audit_rhs_count_max'] = max(a['audit_rhs_count'] for a in group)
        for metric in AUDIT_METRICS:
            values = [a[metric] for a in group if a[metric] >= 0]
            for suffix, value in zip(('count', 'min', 'median', 'max'),
                (len(values), min(values) if values else '',
                 statistics.median(values) if values else '', max(values) if values else '')):
                row[metric+'_'+suffix] = value
        rows.append(row)
    return rows


def render(rows, audits, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    available = [b for b in BACKENDS if any(r['backend'] == b for r in rows)]

    def panel(ax, metric, data, label):
        for backend in available:
            index = BACKENDS.index(backend)
            points = sorted([r for r in data if r['backend'] == backend and
                             r.get(metric+'_count', 0)], key=lambda r: r['threads'])
            if not points:
                continue
            y = [r[metric+'_median']*1000 for r in points]
            ax.errorbar([r['threads'] for r in points], y,
                yerr=[[v-r[metric+'_min']*1000 for v, r in zip(y, points)],
                      [r[metric+'_max']*1000-v for v, r in zip(y, points)]],
                color=COLORS[index], marker=MARKERS[index], label=backend, capsize=3)
        ax.set_xscale('log', base=2)
        ax.set_xticks([1, 2, 4, 8, 16], ['1', '2', '4', '8', '16'])
        ax.set_xlabel('Requested physical thread budget')
        ax.set_ylabel(label)
        ax.grid(alpha=.25)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for ax, metric, label in zip(axes.flat,
        ('process_wall_s', 'process_tree_cpu_s', 'factor_s', 'solve_s'),
        ('Child subprocess wall (ms)', 'Completed child wait4 CPU (ms)',
         'Timed factor call (ms)', 'Timed four-RHS solve call (ms)')):
        panel(ax, metric, rows, label)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=3)
    shape = {metric: ','.join(str(int(r[metric+'_median'])) for r in rows
                              if r[metric+'_count']) for metric in ('n', 'nnz', 'nrhs')}
    shape = {k: ','.join(sorted(set(v.split(',')))) if v else 'unknown' for k, v in shape.items()}
    fig.suptitle(f'Private NEO-2 aggregates: n={shape["n"]}, '
                 f'nnz={shape["nnz"]}, RHS={shape["nrhs"]}\n'
                 'Median and full five-trial range; child includes I/O/check/export, exact audit separate')
    fig.tight_layout(rect=(0, .08, 1, .93))
    fig.savefig(output/'scaling-and-phases.png', dpi=150)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for ax, metric, label in zip(axes, AUDIT_METRICS[:2],
        ('Exact audit + candidate reader wall (ms)', 'Child + per-trial audit wall sum (ms)')):
        panel(ax, metric, audits, label)
    fig.legend(handles, labels, loc='lower center', ncol=3)
    fig.suptitle('Benchmark oracle costs separately timed; no forward certificate or whole-application timing')
    fig.tight_layout(rect=(0, .14, 1, .92))
    fig.savefig(output/'audit-costs.png', dpi=150)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, metric, label in zip(axes, ('process_wall_s', 'process_tree_cpu_s'),
                                 ('Child subprocess wall (ms)', 'Completed child CPU (ms)')):
        for backend in available:
            if backend == 'umfpack':
                continue  # Its parent RSS is not aggregate memory.
            index = BACKENDS.index(backend)
            for r in [r for r in rows if r['backend'] == backend and r[metric+'_count']
                      and r['parent_rss_kib_count']]:
                x, y = r[metric+'_median']*1000, r['parent_rss_kib_median']/1024
                ax.errorbar(x, y, xerr=[[x-r[metric+'_min']*1000], [r[metric+'_max']*1000-x]],
                    yerr=[[y-r['parent_rss_kib_min']/1024], [r['parent_rss_kib_max']/1024-y]],
                    marker=MARKERS[index], color=COLORS[index], capsize=3)
        ax.set_xlabel(label)
        ax.set_ylabel('Single-process VmHWM peak (MiB)')
        ax.grid(alpha=.25)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=COLORS[BACKENDS.index(b)], marker=MARKERS[BACKENDS.index(b)],
                     ls='', label=b) for b in available if b != 'umfpack']
    fig.legend(handles=handles, loc='lower center', ncol=3)
    fig.suptitle('Five requested budgets per backend; see observed-range assessments\n'
                 'UMFPACK IPC aggregate RSS unknown: no point; exact audit costs excluded')
    fig.tight_layout(rect=(0, .12, 1, .9))
    fig.savefig(output/'time-rss.png', dpi=150)
    plt.close(fig)


def describe(receipt, wrappers, rows, audits, fronts, raw):
    counts = collections.Counter(w['status'] for w in wrappers)
    text = ['# Private real NEO-2 aggregate benchmark', '',
        f'Fixture `{receipt["fixture_id"]}`, source `{receipt["source_revision"]}`. '
        f'{len(wrappers)} attempts: '+', '.join(f'{n} {s}' for s, n in counts.items())+'.', '',
        f'Raw private attempts SHA256 `{digest(raw)}`. Frozen accuracy-policy SHA256 '
        f'`{receipt["policy_sha256"]}`; manifest SHA256 `{receipt["manifest_sha256"]}`.', '',
        'One real NEO-2 matrix, n=4,856, nnz=28,898. Four RHS contain three physical '
        'columns (one duplicate) and one manufactured known-solution column. Five '
        'fresh trials for each backend and requested budget 1/2/4/8/16; all five '
        'external peers qualified with the same frozen executable set before native '
        'measurement. No matrix, RHS, candidate solution or per-RHS numeric audit is published.', '',
        'One requested thread: successful child-process medians [min,max]; timed '
        'factor/solve medians are separate library calls. RSS is process VmHWM.', '',
        '| Backend | Child wall ms | Completed child CPU ms | Factor ms | Four-RHS solve ms | Parent peak MiB |',
        '|---|---:|---:|---:|---:|---:|']
    def span(r, metric, scale=1):
        if not r.get(metric+'_count', 0):
            return 'unavailable'
        return (f'{r[metric+"_median"]*scale:.3f} '
                f'[{r[metric+"_min"]*scale:.3f},{r[metric+"_max"]*scale:.3f}]')
    def median(r, metric, scale=1):
        return f'{r[metric+"_median"]*scale:.3f}' if r.get(metric+'_count', 0) else 'unavailable'
    for r in rows:
        if r['threads'] == 1:
            text.append(f'| {r["backend"]} | {span(r,"process_wall_s",1000)} | '
                        f'{span(r,"process_tree_cpu_s",1000)} | {median(r,"factor_s",1000)} | '
                        f'{median(r,"solve_s",1000)} | {span(r,"parent_rss_kib",1/1024)} |')
    text += ['', 'Observed-range child-time/RSS assessments against every external peer/thread configuration:', '']
    for metric in ('process_wall_s', 'process_tree_cpu_s'):
        c = collections.Counter(r['status'] for r in fronts if r.get('metric') == metric)
        text.append(f'- {metric}: '+', '.join(f'{n} {status}' for status, n in c.items())+'.')
    if not all('metric' in r for r in fronts):
        text.append('- Failed or missing trials make frontier assessment indeterminate.')
    text += ['', 'Exact benchmark audit and timing limits:', '']
    for backend in sorted({a['backend'] for a in audits}):
        group = [r for r in audits if r['backend'] == backend]
        wall = [r['audit_reader_wall_s_median'] for r in group if r['audit_reader_wall_s_count']]
        norm = [r['audit_normwise_backward_max_max'] for r in group if r['audit_normwise_backward_max_count']]
        scaled = [r['audit_row_scaled_backward_max_max'] for r in group if r['audit_row_scaled_backward_max_count']]
        wall_text = f'{min(wall):.3f}–{max(wall):.3f}s' if wall else 'unavailable'
        norm_text = f'{max(norm):.3e}' if norm else 'unavailable'
        scaled_text = f'{max(scaled):.3e}' if scaled else 'unavailable'
        text.append(f'- {backend}: per-budget median audit+reader wall {wall_text}; '
                    f'maximum exact-audit normwise backward diagnostic {norm_text}, '
                    f'row-scaled {scaled_text}.')
    text += ['', '- Completed child CPU uses wait4 user+system usage including reaped helpers; '
        'helper CPU must not be added again. Child wall includes input/reference I/O, '
        'original-A checks, candidate export and teardown. Factor/solve calls occupy '
        'only part of it; the remainder is not an isolated I/O measurement.',
        '- Exact rational residual auditing runs separately after each child. It is a '
        'benchmark correctness oracle and empirical forward regression, not a rigorous '
        'forward-error or continuous certificate. Child timings exclude this audit. '
        'Audit CPU and peak memory were not measured. Child RSS is not audit-inclusive '
        'application memory. The child+audit wall sum excludes prerequisite audits, hashing, dependency '
        'inspection and other orchestration; it is not whole-application elapsed time.',
        '- Native active factor teams and external unknown-team reports are retained '
        'in the numeric CSV. Requested thread budgets are not actual measured team sizes.',
        '- Observed-range nondominance uses strict native max time < external min time '
        'OR native max RSS < external min RSS against all 25 external configurations. '
        'IPC aggregate RSS is unknown; only parent VmHWM is a safe lower bound for '
        'excluding a potential dominator. Helper lifecycle ru_maxrss is separate and '
        'may include inherited highwater. IPC has no aggregate point in the figures.',
        '- Nondominance is distinct from matching/beating the same comparator on both '
        'axes. No broad matrix-class admission, certificate-inclusive timing or whole '
        'application result follows from this one frozen private fixture.',
        '- Every wrapper, child and exact-audit status is preserved in aggregate-safe '
        '`attempt-statuses.csv`; failed trial counts are retained. Empty numeric summaries '
        'with count zero mean unavailable, never zero-valued performance.']
    loads = [w['process_receipt']['background_load_before'][0] for w in wrappers
             if 'background_load_before' in w.get('process_receipt', {})]
    if loads:
        text += [f'- Pre-launch one-minute host load {min(loads):.3f}–{max(loads):.3f}; '
                 'controller quiescence and fixed affinity do not make the host exclusive.']
    return '\n'.join(text)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plots', action='store_true')
    args = parser.parse_args()
    raw = args.input/'attempts.jsonl'
    wrappers = [json.loads(l) for l in raw.read_text().splitlines() if l]
    receipt = json.loads((args.input/'receipt.json').read_text())
    if receipt.get('mode') != 'measure' or digest(raw) != receipt.get('raw_attempts_sha256'):
        parser.error('completed immutable measure receipt and matching raw hash required')
    qualification = json.loads((Path(__file__).parent/'private_real_qualification.json').read_text())
    if receipt['fixture_id'] != qualification['fixture_id'] or receipt['policy_sha256'] != qualification['policy_sha256']:
        parser.error('report requires the qualified private fixture and frozen policy')
    if (receipt['fresh_process_trials'] != 5 or
            set(receipt['requested_thread_budgets']) != {1, 2, 4, 8, 16} or
            set(receipt['comparison_scope']) != set(PEERS)|{'native-amd'}):
        parser.error('report requires the frozen six-peer, five-trial, five-budget scope')
    expected = set(itertools.product(receipt['comparison_scope'],
        receipt['requested_thread_budgets'], range(receipt['fresh_process_trials'])))
    observed = [(w['backend'], w['threads'], w['trial']) for w in wrappers]
    if not wrappers or len(observed) != len(set(observed)):
        parser.error('nonempty distinct measurement attempts required')
    pairs = [normalize(w, receipt['fixture_id']) for w in wrappers]
    normalized, safe = map(list, zip(*pairs))
    rows, audits = summarize(normalized), audit_summary(safe)
    complete = set(observed) == expected and all(r['ok'] == r['attempts'] == 5 and r['parent_rss_kib_count'] == 5
                   and r['process_tree_cpu_s_count'] == 5 for r in rows)
    fronts = frontier(rows) if complete else [dict(status='indeterminate_failed_or_missing_trials')]
    args.output.mkdir(parents=True, exist_ok=True)
    for name, data in [('summary.csv', rows), ('audits.csv', audits), ('fronts.csv', fronts),
                       ('attempt-statuses.csv', safe)]:
        with (args.output/name).open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    allowed = ('mode', 'source_revision', 'fixture_id', 'manifest_sha256', 'policy_sha256',
               'window_note', 'qualification_sha256', 'required_external_peers',
               'comparison_scope', 'requested_thread_budgets', 'fresh_process_trials',
               'script_sha256', 'preparation_and_prerequisite_audit_wall_s', 'timing_scope',
               'claims_unavailable', 'status', 'adapter_wall_s', 'raw_attempts_sha256')
    public = {k: receipt[k] for k in allowed if k in receipt}
    public['report_script_sha256'] = {name: digest(Path(__file__).parent/name)
        for name in ('report_private_real.py', 'report.py', 'report_front.py')}
    public['missing_attempts'] = len(expected-set(observed))
    public['privacy'] = 'aggregate metrics and hashes only; numeric inputs/candidates remain private'
    public['binaries'] = [dict(backend=b, path=next(w['binary'] for w in wrappers if w['backend'] == b),
        sha256=next(w.get('binary_sha256', 'unavailable') for w in wrappers if w['backend'] == b))
        for b in sorted({w['backend'] for w in wrappers})]
    public['raw_private_path'] = str(raw.resolve())
    public['receipt_private_path'] = str((args.input/'receipt.json').resolve())
    public['receipt_sha256'] = digest(args.input/'receipt.json')
    (args.output/'recipe.json').write_text(json.dumps(public, indent=2)+'\n')
    text = describe(receipt, wrappers, rows, audits, fronts, raw)
    (args.output/'report.md').write_text(text)
    if args.plots:
        render(rows, audits, args.output)
    print(json.dumps(dict(attempts=len(wrappers), groups=len(rows),
        wrapper_status=collections.Counter(w['status'] for w in wrappers),
        assessments=collections.Counter(r['status'] for r in fronts))))


if __name__ == '__main__':
    main()
