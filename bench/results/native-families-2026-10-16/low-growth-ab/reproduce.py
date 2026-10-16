#!/usr/bin/env python3
"""Regenerate frozen ordinary LOW-growth A/B evidence without running solvers."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path
import re
import statistics
import sys

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO/'bench/native_families'))
from report import KEYS, summarize
from run import parse_csv
BASE = Path('/home/ert/code/sparse-dev/benchmarks/2026-10-16-native-families')
RAW_SHA = '7292b8d45799c621f12e06d9a066450625aad6e73e94407948bd0181c2495832'
MANIFEST_SHA = '7fc0172f97052a1f2be318afaa3b3421c566435d6534ef36eb146739a7158697'
DIAGNOSTIC_MANIFEST_SHA = 'abd99975d46efdab434eb484efadd25d00393a32975bc78ace6a9e1f7a86127c'
DIAGNOSTIC_FILES = ('README.md', 'before/replay.csv', 'before/replay.stderr',
    'candidate/final-replay.csv', 'candidate/final-replay.stderr',
    'candidate/final-dense-exact.txt', 'candidate/final-dense-under.txt',
    'profile/before-report.txt')
POLICIES = ('native-amd', 'native-amd-capped', 'native-amd-low')
VERSIONS = ('baseline', 'candidate')
EXTRA = ('host_current_bytes', 'memory_policy', 'host_memory_limit_bytes')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def aggregate(raw):
    normalized = []
    for a in raw:
        point = dict(a['requested_case'], threads=1, operation='N')
        normalized.append(dict(a, **point))
    rows = []
    for version in VERSIONS:
        sample = [a for a in normalized if a['version'] == version]
        for r in summarize(sample):
            group = [a for a in sample if all(a[k] == r[k] for k in KEYS)]
            row = dict(version=version, **r)
            for metric in EXTRA:
                values = [a['record'][metric] for a in group if a['status'] == 'ok'
                          and a.get('returncode') == 0 and a['record'].get(metric, -1) >= 0]
                for suffix, value in zip(('count', 'min', 'median', 'max'),
                    (len(values), min(values) if values else '', statistics.median(values)
                     if values else '', max(values) if values else '')):
                    row[metric+'_'+suffix] = value
            rows.append(row)
    return rows


def diagnostic(path):
    result = []
    for version, stem in [('baseline', 'before/replay'), ('candidate', 'candidate/final-replay')]:
        text = (path/(stem+'.csv')).read_text()
        growth, copied, requested = map(int, re.search(
            r'DIAGNOSTIC_GROWTH_COPY_REQUEST\s+(\d+)\s+(\d+)\s+(\d+)', text).groups())
        process = re.search(r'allocation_calls=(\d+) allocation_requested_bytes=(\d+)',
                            (path/(stem+'.stderr')).read_text())
        metrics = parse_csv(text)
        result.append(dict(version=version, growth_events=growth,
            cumulative_factor_prefix_copied_bytes=copied,
            cumulative_factor_capacity_requested_bytes=requested,
            cumulative_process_allocation_calls=int(process[1]),
            cumulative_process_allocation_requested_bytes=int(process[2]),
            retained_array_payload_bytes=int(metrics['host_current_bytes']),
            phase_maximum_array_payload_bytes=int(metrics['host_peak_bytes']),
            factor_nnz=int(metrics['factor_nnz'])))
    return result


def render(rows, out, watchdog):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors, markers = ('#0072B2', '#D55E00'), ('o', 's')
    labels = ('AUTO\nno cap', 'AUTO\n64 MiB cap', 'LOW\nno cap')
    for variant in ('training', 'heldout'):
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        for i, family in enumerate(('drift-real', 'drift-complex')):
            for j, workload in enumerate(('rhs8', 'refactor')):
                ax = axes[i, j]
                for version, color, marker, offset in zip(VERSIONS, colors, markers, (-.1, .1)):
                    for k, policy in enumerate(POLICIES):
                        r = next(r for r in rows if r['version'] == version and
                                 r['family'] == family and r['workload'] == workload and
                                 r['variant'] == variant and r['backend'] == policy)
                        if r['process_wall_s_count']:
                            value = r['process_wall_s_median']*1000
                            ax.errorbar(k+offset, value,
                                yerr=[[value-r['process_wall_s_min']*1000],
                                      [r['process_wall_s_max']*1000-value]],
                                marker=marker, color=color, capsize=4)
                        if r['timeouts']:
                            ax.scatter(k+offset, watchdog*1000, marker='x', color=color)
                            ax.annotate(f'{r["timeouts"]} timeouts', (k+offset, watchdog*1000),
                                        xytext=(-12, -18), textcoords='offset points', fontsize=9)
                ax.set_xticks(range(3), labels)
                ax.set_yscale('log')
                ax.set_ylabel('Complete subprocess wall (ms)')
                ax.set_title(f'{family}, n={16384 if i == 0 else 4096}, {workload}')
                ax.grid(axis='y', alpha=.25)
        from matplotlib.lines import Line2D
        handles = [Line2D([], [], color=c, marker=m, ls='', label=v)
                   for c, m, v in zip(colors, markers, VERSIONS)]
        fig.legend(handles=handles, loc='lower center', ncol=2)
        fig.suptitle(f'{variant}: ordinary interleaved growth-policy A/B, CPU0\n'
                     'Median/full five-trial range; × = censored at 0.9s, not a completed solve time')
        fig.tight_layout(rect=(0, .05, 1, .93))
        fig.savefig(out/f'runtime-{variant}.png', dpi=150)
        plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for i, family in enumerate(('drift-real', 'drift-complex')):
        for j, (metric, scale, ylabel) in enumerate([
            ('host_current_bytes', 1e-6, 'Retained owned array payload (decimal MB)'),
            ('parent_rss_kib', 1/1024, 'Whole-process VmHWM peak (binary MiB)')]):
            ax = axes[i, j]
            for version, color, marker, offset in zip(VERSIONS, colors, markers, (-.1, .1)):
                for k, policy in enumerate(POLICIES):
                    r = next(r for r in rows if r['version'] == version and r['family'] == family
                             and r['variant'] == 'training' and r['workload'] == 'rhs8'
                             and r['backend'] == policy)
                    if not r[metric+'_count']:
                        continue  # No memory sample from a killed baseline LOW solve.
                    value = r[metric+'_median']*scale
                    ax.errorbar(k+offset, value,
                        yerr=[[value-r[metric+'_min']*scale], [r[metric+'_max']*scale-value]],
                        marker=marker, color=color, capsize=4)
            ax.set_xticks(range(3), labels)
            ax.set_ylabel(ylabel)
            ax.set_title(f'{family}, training RHS8')
            ax.grid(axis='y', alpha=.25)
    fig.legend(handles=handles, loc='lower center', ncol=2)
    fig.suptitle('Ordinary memory diagnostics: retained payload differs from process peak RSS\n'
                 'Candidate real LOW owns fewer bytes than AUTO but has higher RSS; baseline real LOW has no sample')
    fig.tight_layout(rect=(0, .05, 1, .93))
    fig.savefig(out/'owned-and-rss.png', dpi=150)
    plt.close(fig)


def describe(rows, manifest, raw, counters):
    text = ['# LOW growth: ordinary native A/B and separate diagnostics', '',
        f'Baseline `{manifest["sources"]["baseline"]}`, candidate `{manifest["sources"]["candidate"]}`. '
        '240 interleaved ordinary attempts on CPU0: baseline 100 passed/20 timed out; '
        'candidate 120 passed. Each version/policy/case/workload/variant has five trials. '
        f'Raw SHA256 `{digest(raw)}`.', '',
        'Policies are AUTO without a cap, AUTO with a 64 MiB cap, and LOW without a cap. '
        'The optional positive LOW-cap campaign is not present here. Native factorization '
        'is serial. No external-frontier or thread-scaling claim follows.', '',
        'LOW complete process wall milliseconds, median [minimum,maximum]:', '',
        '| Family | Workload | Variant | Baseline LOW | Candidate LOW |',
        '|---|---|---|---:|---:|']
    for family in ('drift-real', 'drift-complex'):
        for workload in ('rhs8', 'refactor'):
            for variant in ('training', 'heldout'):
                cells = []
                for version in VERSIONS:
                    r = next(r for r in rows if r['family'] == family and r['workload'] == workload
                             and r['variant'] == variant and r['version'] == version
                             and r['backend'] == 'native-amd-low')
                    cells.append(f'{r["process_wall_s_median"]*1000:.3f} '
                        f'[{r["process_wall_s_min"]*1000:.3f},{r["process_wall_s_max"]*1000:.3f}]'
                        if r['process_wall_s_count'] else '5 timeouts ≥900 ms')
                text.append(f'| {family} | {workload} | {variant} | {cells[0]} | {cells[1]} |')
    text += ['', 'Memory and scope:', '',
        '- Candidate real training RHS8 LOW: retained 10,611,660 bytes, median 33,440 KiB '
        '(32.66 MiB) RSS. AUTO: retained 13,891,596 bytes, median 29,948 KiB (29.25 MiB) '
        'RSS. Lower retained owned payload does not establish lower RSS; this LOW RSS '
        'is higher. Allocator behavior is investigated separately, not attributed here.',
        '- `summary.csv` retains all 48 version/configuration groups, timeouts, min/median/max '
        'wall and completed wait4 CPU, factor/solve/refactor, payload, process RSS and '
        'numerical errors. Timeout CPU is incomplete and remains unavailable. No finite '
        'completed-time speedup is assigned to the real baseline LOW timeouts.',
        '- Aggregation reads only raw receipts, never the archived summary.json. '
        'Its corrected filter required no new solver invocations. The candidate ELF '
        'also links METIS while the baseline does not; library identities are retained. '
        'Small AUTO timing differences are descriptive, not an isolated algorithm gain.',
        '- The source and binary identities are the measured snapshots above. Later '
        'production integrations f4f329d/9354c60 do not change this campaign’s identity. '
        'AUTO+capped64MiB is not LOW+capped64MiB.', '',
        'Separate instrumented profile/counter evidence:', '',
        '| Quantity | Exact-growth diagnostic | Small-slack diagnostic |',
        '|---|---:|---:|']
    for key in counters[0]:
        if key != 'version':
            text.append(f'| {key} | {counters[0][key]:,} | {counters[1][key]:,} |')
    text += ['', '- The 54,871,019,988 copied bytes (~54.87 decimal GB) are cumulative traffic, '
        'not peak memory. Diagnostic phase-maximum payload rises from 14,899,228 to '
        '15,187,408 bytes; retained payload also rises slightly versus old exact LOW.',
        '- The historical perf sample uses source 3b9384b: 740 samples, zero lost, reserve '
        '64.99% own cycles and visit 5.23%. It is separate from the ordinary 7e1f123 A/B '
        'snapshot. Instrumented diagnostic runtimes and allocations are not ordinary '
        'performance measurements. Whole-process allocation counters also include runtime work.',
        '- Dense80 diagnostic evidence retains one exact retry after the hinted growth '
        'attempt: success with peak 144,812 bytes at that cap, failure one byte below, '
        'required 38,880 versus available 38,879. Actual peak is max across sequential '
        'attempts, never their sum. These diagnostics are not the pending ordinary '
        'positive LOW-cap campaign. Source-only final cleanup is distinguished in the '
        'diagnostic archive from the final replay’s instrumented algorithm.']
    return '\n'.join(text)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=BASE/'low-growth-ab')
    parser.add_argument('--diagnostic', type=Path, default=BASE/'low-growth-diagnostic')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plots', action='store_true')
    args = parser.parse_args()
    raw = args.input/'attempts.jsonl'
    if digest(raw) != RAW_SHA or digest(args.input/'manifest.json') != MANIFEST_SHA:
        parser.error('frozen ordinary evidence identity differs')
    if digest(args.diagnostic/'manifest.json') != DIAGNOSTIC_MANIFEST_SHA:
        parser.error('frozen instrumented diagnostic identity differs')
    diagnostic_manifest = json.loads((args.diagnostic/'manifest.json').read_text())
    expected = {r['path']: r['sha256'] for r in diagnostic_manifest['files']}
    if any(digest(args.diagnostic/name) != expected[name] for name in DIAGNOSTIC_FILES):
        parser.error('diagnostic source/counter evidence differs from its manifest')
    if args.output.resolve().is_relative_to(REPO):
        parser.error('generated plots must stay outside source history')
    data = [json.loads(l) for l in raw.read_text().splitlines() if l]
    manifest = json.loads((args.input/'manifest.json').read_text())
    rows, counters = aggregate(data), diagnostic(args.diagnostic)
    if len(data) != 240 or len(rows) != 48 or any(r['attempts'] != 5 for r in rows):
        parser.error('complete frozen five-trial A/B groups required')
    args.output.mkdir(parents=True, exist_ok=True)
    for filename, values in [('summary.csv', rows), ('diagnostic-counters.csv', counters)]:
        with (args.output/filename).open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(values[0]))
            writer.writeheader(); writer.writerows(values)
    (args.output/'report.md').write_text(describe(rows, manifest, raw, counters))
    (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    files = [raw, args.input/'manifest.json', args.input/'run_ab.py',
        args.diagnostic/'manifest.json']+[args.diagnostic/name for name in DIAGNOSTIC_FILES]
    (args.output/'recipe.json').write_text(json.dumps(dict(
        ordinary_sources=manifest['sources'], input=str(args.input), diagnostic=str(args.diagnostic),
        files=[dict(path=str(p),sha256=digest(p)) for p in files],
        recipe_sha256=digest(Path(__file__)),
        report_sources={name: digest(REPO/'bench/native_families'/name)
                        for name in ('report.py', 'run.py')}, coverage='native-only policy A/B, CPU0; no frontier admission',
        omitted_campaign='ordinary LOW with positive cap not measured here',
        units='payload bytes; cumulative copy/allocation bytes are traffic, not peak; RSS KiB'
        ), indent=2)+'\n')
    if args.plots:
        render(rows, args.output, manifest['watchdog_s'])
    print(json.dumps(dict(attempts=len(data), groups=len(rows),
        statuses=collections.Counter((r['version']+':'+r['status']) for r in data))))


if __name__ == '__main__':
    main()
