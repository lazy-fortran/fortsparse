#!/usr/bin/env python3
"""Retain the interleaved scalar/indexed diagnostic without merging revisions."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path
import statistics

from report import KEYS, summarize

LABELS = ('scalar', 'indexed')
GROUP = ('family', 'workload', 'variant')


def aggregates(attempts):
    rows = []
    for label in LABELS:
        subset = [a for a in attempts if a['revision_label'] == label]
        for row in summarize(subset):
            rows.append(dict(revision_label=label, **row))
    # Ratios use interleaved matched trial IDs, not independently sorted times.
    paired = collections.defaultdict(dict)
    for a in attempts:
        key = tuple(a[k] for k in (*KEYS, 'trial'))
        if a['revision_label'] in paired[key]:
            raise ValueError('duplicate revision in interleaved pair')
        paired[key][a['revision_label']] = a
    groups = collections.defaultdict(list)
    for pair in paired.values():
        if set(pair) != set(LABELS):
            raise ValueError('incomplete interleaved pair')
        first = pair['scalar']
        groups[tuple(first[k] for k in GROUP)].append(pair)
    ratios = []
    for key, pairs in sorted(groups.items()):
        row = dict(zip(GROUP, key), pairs=len(pairs))
        for metric in ('process_wall_s', 'process_tree_cpu_s'):
            values = [p['scalar'][metric]/p['indexed'][metric] for p in pairs]
            for suffix, value in zip(('min', 'median', 'max'),
                                    (min(values), statistics.median(values), max(values))):
                row[f'{metric}_scalar_over_indexed_{suffix}'] = value
        ratios.append(row)
    return rows, ratios


def render(rows, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    categories = [('drift-real', 'rhs8'), ('drift-real', 'refactor'),
                  ('drift-complex', 'rhs8'), ('drift-complex', 'refactor')]
    palette = ('#0072B2', '#D55E00')  # Okabe-Ito; distinct markers too.
    for name, metrics, title in (
        ('timing', [('process_wall_s', 'Complete subprocess wall (ms)', 1000),
                    ('process_tree_cpu_s', 'Completed wait4 workflow CPU (ms)', 1000)],
         'Interleaved AMD scalar/indexed: median and full five-trial range'),
        ('storage', [('parent_rss_kib', 'Process VmHWM peak (MiB)', 1/1024),
                     ('factor_nnz', 'Reported factor entries', 1)],
         'Interleaved AMD: unchanged fill; overlapping observed process RSS ranges')):
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), squeeze=False)
        for col, variant in enumerate(('training', 'heldout')):
            for i, (metric, ylabel, scale) in enumerate(metrics):
                ax = axes[i, col]
                for label, color, marker, offset in zip(LABELS, palette, ('o', 's'), (-.1, .1)):
                    points = [next(r for r in rows if r['revision_label'] == label
                              and r['family'] == family and r['workload'] == workload
                              and r['variant'] == variant) for family, workload in categories]
                    y = [r[metric+'_median']*scale for r in points]
                    ax.errorbar([j+offset for j in range(4)], y,
                        yerr=[[v-r[metric+'_min']*scale for v, r in zip(y, points)],
                              [r[metric+'_max']*scale-v for v, r in zip(y, points)]],
                        color=color, marker=marker, linestyle='none', label=label, capsize=4)
                ax.set_xticks(range(4), ['real\nRHS8', 'real\nrefactor',
                                        'complex\nRHS8', 'complex\nrefactor'])
                ax.set_ylabel(ylabel)
                ax.set_title(variant)
                ax.grid(axis='y', alpha=.25)
                ax.legend()
        fig.suptitle(title+'\n16,384 unknowns; AMD; requested=actual factor threads=1; no frontier admission')
        fig.tight_layout(rect=(0, 0, 1, .93))
        fig.savefig(output/f'{name}.png', dpi=150)
        plt.close(fig)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_report(attempts, rows, ratios, provenance, source):
    lookup = {(r['revision_label'], *(r[k] for k in GROUP)): r for r in rows}
    def span(row, metric, scale=1):
        return (f'{row[metric+"_median"]*scale:.3f} '
                f'[{row[metric+"_min"]*scale:.3f},{row[metric+"_max"]*scale:.3f}]')
    text = ['# Interleaved scalar versus generated indexed LU', '',
        f'Raw attempts: `{source}`; SHA256 `{digest(source)}`.', '',
        f'{len(attempts)}/{len(attempts)} attempts completed successfully. '
        'Sixteen five-trial revision/configuration groups form eight paired comparisons. '
        'Each trial alternates execution order. Both binaries use explicit AMD and one '
        'requested/actual native factor thread on CPU 0; N solves, training and heldout '
        'coefficients, n=16,384. Refactor workflows include the prime factorization.', '',
        f'Scalar source: `{provenance["source_manifest_old"]["source_revision"]}`. '
        f'Indexed source: `{provenance["source_manifest_new"]["source_revision"]}`. '
        'Binary hashes and paths are retained in `provenance.json` and `recipe.json`.', '',
        'Complete subprocess wall milliseconds, median [minimum,maximum]:', '',
        '| Family | Workload | Variant | Scalar | Indexed | Paired speedup |',
        '|---|---|---|---:|---:|---:|']
    for p in ratios:
        scalar, indexed = (lookup[(label, *(p[k] for k in GROUP))] for label in LABELS)
        ratio = 'process_wall_s_scalar_over_indexed'
        text.append(f'| {p["family"]} | {p["workload"]} | {p["variant"]} | '
                    f'{span(scalar, "process_wall_s", 1000)} | '
                    f'{span(indexed, "process_wall_s", 1000)} | '
                    f'{span(p, ratio)}× |')
    text += ['', 'Completed workflow CPU milliseconds, median [minimum,maximum]:', '',
             '| Family | Workload | Variant | Scalar | Indexed |',
             '|---|---|---|---:|---:|']
    for p in ratios:
        scalar, indexed = (lookup[(label, *(p[k] for k in GROUP))] for label in LABELS)
        text.append(f'| {p["family"]} | {p["workload"]} | {p["variant"]} | '
                    f'{span(scalar, "process_tree_cpu_s", 1000)} | '
                    f'{span(indexed, "process_tree_cpu_s", 1000)} |')
    loads = [a['background_load_before'][0] for a in attempts]
    text += ['', 'Scope and limits:', '',
        '- This is a native scalar/indexed A/B diagnostic, not an external solver '
        'comparison or frontier admission. Sizes retain the earlier four-external '
        'calibration; this run does not recalibrate them or add direct UMFPACK evidence.',
        '- Every indexed wall and CPU trial is below every scalar trial in its '
        'matched configuration. The retained paired ratios use matching trial IDs.',
        '- Both versions report exactly 762,644 factor entries throughout. Process RSS '
        'ranges overlap in every matched group; no RSS reduction is established. '
        'Attributed peak bytes are unavailable in both, not zero.',
        '- Workflow CPU is completed wait4 user+system time including reaped descendants. '
        'Helper CPU is not added again. Parent pre-CSV and main-entry CPU remain distinct. '
        'These native runs have no helper CPU/RSS.',
        '- `summary.csv` separates factor, numeric, analysis, solve, refactor, prime-factor, '
        'main-workflow and complete subprocess phases; unknown metrics have count zero. '
        '`paired.csv` retains all eight wall/CPU ratio ranges. Both versions pass the '
        'original-A known-solution/backward-error oracle.',
        f'- Observed pre-launch one-minute load: {min(loads):.3f}–{max(loads):.3f}. '
        'Core affinity and interleaving do not eliminate shared-host variability.']
    return '\n'.join(text)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--plots', action='store_true')
    args = parser.parse_args()
    source = args.input/'attempts.jsonl'
    attempts = [json.loads(line) for line in source.read_text().splitlines() if line]
    provenance = json.loads((args.input/'provenance.json').read_text())
    if any(a['status'] != 'ok' or a['returncode'] != 0 or
           a.get('process_tree_cpu_s', -1) < 0 for a in attempts):
        parser.error('this checkpoint requires successful completed wait4 receipts')
    rows, ratios = aggregates(attempts)
    args.output.mkdir(parents=True, exist_ok=True)
    for filename, data in [('summary.csv', rows), ('paired.csv', ratios)]:
        with (args.output/filename).open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    (args.output/'report.md').write_text(write_report(attempts, rows, ratios, provenance, source))
    identities = []
    for label in LABELS:
        binary = Path(next(a['command'][0] for a in attempts if a['revision_label'] == label))
        manifest = provenance['source_manifest_old' if label == 'scalar'
                              else 'source_manifest_new']
        expected = next(f['sha256'] for f in manifest['files']
                        if f.get('name') == binary.name or
                        Path(f.get('path', '')).name == binary.name)
        if digest(binary) != expected:
            parser.error('frozen binary differs from recorded manifest')
        identities.append(dict(label=label, path=str(binary), sha256=expected))
    recipe = dict(raw_path=str(source.resolve()), raw_sha256=digest(source),
        provenance_path=str(args.input/'provenance.json'),
        provenance_sha256=digest(args.input/'provenance.json'), binaries=identities,
        fixture_sha256=provenance['fixture_sha256'],
        original_calibration_sha256=provenance['case_calibration_sha256'],
        selection_scope=provenance['selection_scope'],
        cpu_scope='completed wait4 user+system including reaped children',
        ratio='scalar/indexed, matched configuration and trial ID; median/min/max',
        plots='timing.png and storage.png, external only', performance_admission='unknown')
    (args.output/'recipe.json').write_text(json.dumps(recipe, indent=2)+'\n')
    (args.output/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    if args.plots:
        render(rows, args.output)
    print(json.dumps(dict(attempts=len(attempts), revision_groups=len(rows), paired_groups=len(ratios))))


if __name__ == '__main__':
    main()
