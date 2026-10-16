#!/usr/bin/env python3
"""Summarize retained receipts without running benchmarks or exporting private data."""
import argparse
import csv
import hashlib
import json
import pathlib
import re
import statistics


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def triples(values):
    return dict(zip(('median', 'min', 'max'),
                    (statistics.median(values), min(values), max(values))))


def write_csv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def summarize(args):
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    public = [json.loads(line) for line in (args.public / 'attempts.jsonl').read_text().splitlines()]
    # Private receipts contain numeric candidates and per-RHS exact fractions.
    # Only explicit aggregate scalar fields below are copied into output rows.
    private = [json.loads(line) for line in (args.private / 'attempts.jsonl').read_text().splitlines()]
    assert len(public) == 80 and len(private) == 10
    assert all(r['status'] == 'ok' and r['returncode'] == 0 for r in public)
    assert all(r['status'] == 'pass' and r['process_receipt']['status'] == 'ok'
               and r['independent_audit']['status'] == 'pass' for r in private)
    groups = {}
    for receipt in public:
        record = receipt['record']
        key = ('public', record['family'], record['level'], record['variant'],
               record['workload'], receipt['version'])
        groups.setdefault(key, []).append((receipt, record))
    for receipt in private:
        key = ('private-aggregate', 'neo2-lorentz-tag6-4856-v1', '', 'frozen',
               'rhs4', 'candidate' if receipt['revision_label'] == 'pruned' else 'baseline')
        child = receipt['process_receipt']
        groups.setdefault(key, []).append((child, child['record']))
    rows = []
    record_metrics = ['total_s', 'factor_s', 'refactor_s', 'solve_s', 'parent_rss_kib',
                      'factor_nnz', 'host_current_bytes', 'host_peak_bytes']
    for key, pairs in groups.items():
        assert len(pairs) == 5
        row = dict(zip(('scope', 'family', 'level', 'variant', 'workload', 'version'), key))
        assert all(r['active_threads'] == 1 for _, r in pairs)
        row.update(attempts=5, passed=5, failures=0, requested_threads=1, active_threads=1)
        for field in ('n', 'nnz', 'nrhs', 'memory_policy', 'host_memory_limit_bytes'):
            values = {r[field] for _, r in pairs if field in r}
            assert len(values) <= 1
            row[field] = next(iter(values), '')
        metrics = {name: [r[name] for _, r in pairs if name in r and r[name] >= 0]
                   for name in record_metrics}
        metrics['process_wall_s'] = [r['process_wall_s'] for r, _ in pairs]
        metrics['completed_cpu_s'] = [r['wait4_reaped_cpu_s'] for r, _ in pairs]
        if key[0] == 'private-aggregate':
            selected = [r for r in private if (r['revision_label'] == 'pruned') == (key[-1] == 'candidate')]
            for name in ('independent_audit_including_reader_wall_s', 'solver_child_plus_audit_wall_s'):
                metrics[name] = [r[name] for r in selected]
        for metric, values in metrics.items():
            for stat in ('median', 'min', 'max'):
                row[f'{metric}_{stat}'] = ''
            for stat, value in triples(values).items() if values else []:
                row[f'{metric}_{stat}'] = value
        rows.append(row)
    for receipt in public:
        if receipt['version'] != 'baseline':
            continue
        r = receipt['record']
        candidate = next(c['record'] for c in public if c['version'] == 'candidate'
                         and c['trial'] == receipt['trial']
                         and all(c['record'][k] == r[k] for k in ('family', 'variant', 'workload')))
        assert candidate['factor_nnz'] == r['factor_nnz']
        assert candidate['host_current_bytes'] - r['host_current_bytes'] == 4 * r['n']
    write_csv(out / 'summary.csv', rows)
    profiles, profile_sources = [], []
    cases = ('drift-real-128-rhs8', 'drift-real-128-refactor',
             'drift-complex-64-rhs8', 'drift-complex-128-refactor')
    for version, directory in [('baseline', 'indexed-profiles'), ('candidate', 'witnessed-pruning-profiles')]:
        folder = args.public.parent / directory
        manifest = json.loads((folder / 'manifest.json').read_text())
        if version == 'baseline':
            frozen = json.loads((args.public.parent / 'indexed-openblas-baseline' / 'binary-manifest.json').read_text())
            profile_binary = frozen['binaries']['bench_native_families']['sha256']
        else:
            profile_binary = manifest['binary_sha256']
        for case in cases:
            report = folder / f'{case}-report.txt'
            stderr = folder / f'{case}.stderr'
            replay = folder / f'{case}-replay.sh'
            text = report.read_text()
            shares = {}
            for percent, symbol in re.findall(r'^\s*(\d+\.\d+)%\s+\[\.\]\s+(\S+)', text, re.M):
                shares[symbol] = shares.get(symbol, 0) + float(percent)
            row = dict(version=version, case=case, source=manifest['source'],
                       binary_sha256=profile_binary,
                       replay_count=5 if version == 'candidate' or 'complex-128' in case else 10,
                       samples=int(re.search(r'\((\d+) samples\)', stderr.read_text())[1]),
                       approximate_event_count=int(re.search(r'Event count \(approx\.\): (\d+)', text)[1]))
            for label, suffix in [('dfs', '_MOD_visit'), ('prune', '_MOD_lu_prune_column'),
                                  ('scatter', '_scatter_update')]:
                row[f'{label}_self_percent'] = sum(v for s, v in shares.items() if suffix in s)
            profiles.append(row)
            profile_sources.append(dict(version=version, case=case,
                                        files=[dict(path=str(p), sha256=sha(p)) for p in
                                               (report, stderr, replay, folder / f'{case}.data')]))
    write_csv(out / 'profiles.csv', profiles)
    pubmeta = json.loads((args.public / 'manifest.json').read_text())
    privmeta = json.loads((args.private / 'provenance.json').read_text())
    recipe = dict(public_manifest=pubmeta,
                  public_raw=dict(path=str(args.public / 'attempts.jsonl'), sha256=sha(args.public / 'attempts.jsonl')),
                  public_runner=dict(path=str(args.public / 'run_ab.py'), sha256=sha(args.public / 'run_ab.py')),
                  private_aggregate=dict(fixture_id='neo2-lorentz-tag6-4856-v1',
                      raw_path=str(args.private / 'attempts.jsonl'), raw_sha256=sha(args.private / 'attempts.jsonl'),
                      **{k: privmeta[k] for k in ('script_sha256', 'frozen_policy_sha256', 'fixture_manifest_sha256')},
                      qualification_receipt_sha256={k: list(v.values()) for k, v in privmeta['qualifications'].items()}),
                  profile_sources=profile_sources,
                  hardware=dict(cpu='AMD Ryzen 9 5950X', physical_cores=16, smt=2, numa_nodes=1),
                  ordinary_recipe=dict(backend='native-amd', families_levels=[['drift-real',128],['drift-complex',64]],
                      variants=['training','heldout'], workloads=['rhs8','refactor'], operation='N', trials=5,
                      affinity=[0], watchdog_s=0.9, environment=public[0]['environment'],
                      interleave='Zero-based even trials baseline,candidate; odd trials candidate,baseline; cases in listed Cartesian-product order.',
                      command='bench_native_families native-amd FAMILY LEVEL 1 WORKLOAD N VARIANT bench/fixtures/neo2-mini5/matrix.mtx',
                      compiler_flags='-O3 -march=znver3; unchanged mathematical/compiler arithmetic policy',
                      cmake_package_pins=dict(AMD_DIR='/usr/lib/cmake/AMD', SuiteSparse_config_DIR='/usr/lib/cmake/SuiteSparse_config')),
                  profile_recipe=dict(event='cycles:u', frequency_hz=1999, call_graph='dwarf,8192', affinity=[0],
                      kind='Nonexclusive instrumented diagnostic; replay counts differ; not ordinary timing.',
                      command='perf record -F 1999 -e cycles:u -g --call-graph dwarf,8192 -o DATA -- bash REPLAY.sh'),
                  acceptance=dict(attempts=90, passed=90, public_backward_error_max=max(r['record']['backward_error'] for r in public),
                                  public_forward_error_max=max(r['record']['forward_error'] for r in public)),
                  pruning_source='7db6a1aef86d2c65b21c95d39bbe5251fd0c1c13',
                  integrated_source='7e1f123', fortsym_revision='7e8ef7abdd90c49bb5f75e2c119e7d181656bcc6',
                  timing_scope='Public and private child wall covers process lifetime; completed CPU is wait4 user+system. Private exact audit timed separately.',
                  unavailable='Blank metric values are unavailable (including negative peak sentinel); not zero.',
                  reproduction='See report.md for commands, interleaving, runtime and privacy scope.')
    (out / 'recipe.json').write_text(json.dumps(recipe, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--public', type=pathlib.Path, required=True)
    parser.add_argument('--private', type=pathlib.Path, required=True)
    parser.add_argument('--output', type=pathlib.Path, required=True)
    summarize(parser.parse_args())
