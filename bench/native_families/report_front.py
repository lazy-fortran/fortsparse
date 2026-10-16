#!/usr/bin/env python3
"""Combine completed five-peer checkpoints and assess observed time/RSS ranges."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path

from report import summarize, render
from run import record_umfpack_configuration
from helper_memory import ENVIRONMENT as HELPER_METRICS_ENV, parse_helper_memory

PEERS = ('superlu', 'umfpack', 'umfpack-direct', 'mumps', 'intel-pardiso')
SIX_PEERS = PEERS + ('umfpack-direct-noir',)
MATCH = ('family', 'level', 'workload', 'variant', 'operation')


def assess(candidate, peers, metric, use_helper_lower_bound=False):
    """Exclude every possible dominator by strict time OR safe RSS separation."""
    unresolved, dominated, by_time, by_memory, by_helper = [], [], [], [], []
    for peer in peers:
        label = f'{peer["backend"]}:{peer["threads"]}'
        if candidate[metric+'_max'] < peer[metric+'_min']:
            by_time.append(label)
        elif candidate['parent_rss_kib_max'] < peer['parent_rss_kib_min']:
            # For IPC this uses only a valid parent peak lower bound.
            by_memory.append(label)
        elif (use_helper_lower_bound and peer['backend'] == 'umfpack' and
              peer.get('helper_postexec_hwm_kib_count', 0) == peer.get('ok', -1) and
              peer.get('helper_postexec_hwm_kib_count', 0) > 0 and
              candidate['parent_rss_kib_max'] < peer['helper_postexec_hwm_kib_min']):
            by_helper.append(label)
        else:
            unresolved.append(label)
        if (peer['backend'] != 'umfpack' and
                peer[metric+'_max'] < candidate[metric+'_min'] and
                peer['parent_rss_kib_max'] < candidate['parent_rss_kib_min']):
            dominated.append(label)
    status = ('dominated' if dominated else 'robust_nondominated' if not unresolved
              else 'indeterminate')
    result = dict(status=status, excluded_by_time=';'.join(by_time),
                  excluded_by_memory=';'.join(by_memory), possible_dominators=';'.join(unresolved),
                  definite_dominators=';'.join(dominated))
    if use_helper_lower_bound:
        result['excluded_by_helper_lower_bound'] = ';'.join(by_helper)
    return result


def frontier(rows, peer_backends=PEERS, use_helper_lower_bound=False):
    result = []
    for candidate in rows:
        if not candidate['backend'].startswith('native'):
            continue
        peers = [r for r in rows if r['backend'] in peer_backends and
                 all(r[k] == candidate[k] for k in MATCH)]
        if {(r['backend'], r['threads']) for r in peers} != {
                (peer, nt) for peer in peer_backends for nt in (1, 2, 4, 8, 16)}:
            raise ValueError(f'incomplete {len(peer_backends)}-comparator/thread evidence')
        for metric in ('process_wall_s', 'process_tree_cpu_s'):
            result.append(dict((k, candidate[k]) for k in (*MATCH, 'backend', 'threads')) |
                          dict(metric=metric, **assess(candidate, peers, metric, use_helper_lower_bound)))
    return result


def validate_coverage(coverage, selection, provenance, attempts):
    peer_backends = PEERS if coverage == 'historical-five' else SIX_PEERS
    if set(selection['backends']) != set(peer_backends):
        raise ValueError(f'{coverage} selection required')
    if coverage == 'contract-v4-six':
        if (selection.get('calibration_scope') != 'complete_declared_baselines' or
                any(p.get('contract', {}).get('schema_version') != 4 or
                    set(p['contract']['required_external_backends']) != set(SIX_PEERS)
                    for p in provenance)):
            raise ValueError('contract-v4 complete six-comparator provenance required')
        for attempt in attempts:
            if attempt['status'] == 'ok':
                record_umfpack_configuration(attempt, attempt['backend'])
                if attempt['status'] != 'ok':
                    raise ValueError('direct UMFPACK controls not confirmed')
    return peer_backends


def validate_helper_lower_bounds(attempts):
    for attempt in attempts:
        if attempt['backend'] != 'umfpack':
            continue
        recorded = attempt.get('record', {}).get('helper_postexec_hwm_kib', -1)
        if recorded <= 0:
            continue
        parsed = parse_helper_memory(attempt.get('stderr', ''))
        if (parsed['helper_postexec_hwm_kib'] != recorded or
                attempt.get('environment', {}).get(HELPER_METRICS_ENV) != '1'):
            raise ValueError('helper lower bound lacks matching opt-in stderr evidence')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def describe(attempts, rows, fronts, calibration, inputs, revision, calibration_attempts,
             peer_backends=PEERS, use_helper_lower_bound=False):
    calibration_status = collections.Counter(a['status'] for a in calibration_attempts)
    rejections = [a for a in calibration_attempts if a['status'] != 'ok'
                  or a['returncode'] != 0 or a['process_wall_s'] >= .5]
    families = sorted({r['family'] for r in rows})
    operations = ','.join(sorted({r['operation'] for r in rows}))
    workloads = ','.join(sorted({r['workload'] for r in rows}))
    native_teams = sorted({int(a['record']['active_threads']) for a in attempts
                           if a['backend'].startswith('native')})
    unknown_peaks = sorted({r['backend'] for r in rows
                            if r['backend'].startswith('native') and not r['host_peak_bytes_count']})
    text = ['# Five-peer completed CPU/RSS checkpoint' if peer_backends == PEERS else
            '# Six-comparator completed CPU/RSS checkpoint', '',
        f'{len(attempts)} completed attempts; {len(rows)} five-trial configuration groups. '
        f'Source `{revision}`. Native backend labels: '+
        ', '.join(sorted({r['backend'] for r in rows if r['backend'].startswith('native')}))+'. '
        'All process CPU uses completed wait4 user+system counters including reaped helpers.', '',
        ('Five external peers: SuperLU, UMFPACK IPC, standalone direct UMFPACK, MUMPS '
         'and Intel oneMKL PARDISO. ' if peer_backends == PEERS else
         'Six external configurations: SuperLU, UMFPACK IPC, direct UMFPACK defaults, '
         'direct UMFPACK IRSTEP=0, MUMPS and Intel oneMKL PARDISO. ') +
        'Requested physical thread budgets are 1/2/4/8/16. '
        f'Native active factor teams reported: {native_teams}. External active teams '
        'are retained separately in the numeric CSV; negative values mean unknown.', '',
        f'External-only calibration retained {len(calibration_attempts):,} attempts: '+
        ', '.join(f'{count} {status}' for status, count in calibration_status.items())+
        f'. Rejected observations under the strict <0.5s process-wall/status gate: '
        f'{len(rejections)}. Rejected observations and unmeasured candidate tails '
        'are retained; their success is never inferred. Earlier differently sized '
        'diagnostics remain separate.', '',
        'The following counts compare each of five native requested-thread configurations '
        f'with all {len(peer_backends)*5} external peer/thread configurations at the same family, actual size, '
        'workload, variant and operation:', '',
        '| Family | Workload | Variant | Wall: robust / indeterminate / dominated | CPU: robust / indeterminate / dominated |',
        '|---|---|---|---:|---:|']
    groups = sorted({(r['family'], r['workload'], r['variant']) for r in fronts})
    for family, workload, variant in groups:
        entries = []
        for metric in ('process_wall_s', 'process_tree_cpu_s'):
            c = collections.Counter(r['status'] for r in fronts if r['family'] == family
                                    and r['workload'] == workload and r['variant'] == variant
                                    and r['metric'] == metric)
            entries.append(' / '.join(str(c[k]) for k in ('robust_nondominated', 'indeterminate', 'dominated')))
        text.append(f'| {family} | {workload} | {variant} | {entries[0]} | {entries[1]} |')
    text += ['', 'Actual measured sizes: '+', '.join(
        f'{family} n={int(next(r["n_median"] for r in rows if r["family"] == family))}'
        for family in sorted({r['family'] for r in rows}))+'.',
        '', 'Assessment rule and limits:', '',
        '- Robust observed-range nondominance excludes every external configuration '
        'by native maximum time < external minimum time OR native maximum process RSS '
        '< external minimum process RSS. Strict inequality treats observed overlap as '
        'indeterminate. These trial ranges are descriptive, not confidence intervals.',
        ('- IPC aggregate RSS is unknown. Only its parent VmHWM is used as a safe lower '
         'bound for excluding a potential dominator. Helper ru_maxrss can retain inherited '
         'highwater; no sum or maximum of parent/helper is an aggregate peak. IPC cannot '
         'be declared a definite dominator with these memory diagnostics.'
         if not use_helper_lower_bound else
         '- IPC aggregate RSS remains unknown. Parent VmHWM and the separately validated '
         'helper post-exec VmHWM are lower bounds only. The opt-in helper rule requires '
         'all five valid trial receipts and native RSS maximum < helper VmHWM minimum; '
         'its exclusions have a separate CSV field. No sum or maximum supplies aggregate '
         'RSS; neither bound can prove IPC domination.'),
        '- A definite dominator has strictly lower maximum time and process RSS than '
        'the native minimum values. Such a result requires known single-process RSS.',
        '- Nondominance is a time/memory tradeoff, distinct from matching or beating '
        'the same best comparator on both axes. This report does not promote a library '
        f'default or claim broad production/frontier completion. It covers {len(families)} '
        f'fixture families, operations {operations}, workloads {workloads} only.',
        '- CPU and elapsed front assessments are separate. The refactor workflow includes '
        'prime factorization. The numeric CSV separates timed analysis, factor, numeric, '
        'solve and refactor calls from complete subprocess wall and CPU.',
        f'- Native labels with unavailable attributed peak bytes: {unknown_peaks}. '
        'Process VmHWM is valid RSS evidence and must not be called attributed factor '
        'memory. Unknown summaries have count zero and must not be plotted as zero.',
        '- Five-peer selection differs from the older four-peer diagnostic. No '
        'cross-checkpoint speedup at different sizes or host loads is inferred.']
    text += ['', 'Calibration rejections:', '']
    text += [f'- {a["family"]}, level {a["level"]}, {a["backend"]}, '
             f'{a["threads"]} threads, trial {a["trial"]}: status {a["status"]}, '
             f'process wall {a["process_wall_s"]:.6f}s.' for a in rejections]
    if not rejections:
        text += ['- None.']
    loads = [a['background_load_before'][0] for a in attempts]
    text += [f'- Pre-launch one-minute host load: {min(loads):.3f}–{max(loads):.3f}. '
             'Affinity and controller quiescence retain ordinary resident services; '
             'this is not an exclusive host.', '', 'Reproduction and provenance:', '']
    text += [f'- `{path}`: raw SHA256 `{digest(path/"attempts.jsonl")}`.' for path in inputs]
    text += [f'- Calibration selection SHA256 `{digest(calibration)}`; all rejected '
             'observations remain in its sibling `attempts.jsonl`.']
    return '\n'.join(text)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', nargs='+', type=Path)
    parser.add_argument('--calibration', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--plots', action='store_true')
    parser.add_argument('--helper-memory-lower-bound', action='store_true',
                        help='use validated complete helper VmHWM ranges only as union lower bounds')
    parser.add_argument('--coverage', choices=('historical-five', 'contract-v4-six'),
                        default='historical-five')
    args = parser.parse_args()
    attempts, provenance = [], []
    for path in args.inputs:
        attempts += [json.loads(l) for l in (path/'attempts.jsonl').read_text().splitlines() if l]
        provenance.append(json.loads((path/'provenance.json').read_text()))
    if any(a['status'] != 'ok' or a['returncode'] != 0 or
           a.get('process_tree_cpu_s', -1) < 0 for a in attempts):
        parser.error('completed successful wait4 receipts required')
    if len({p['source_revision'] for p in provenance}) != 1:
        parser.error('combined inputs must share a frozen source revision')
    selection = json.loads(args.calibration.read_text())
    calibration_attempts = [json.loads(l) for l in
        (args.calibration.parent/'attempts.jsonl').read_text().splitlines() if l]
    try:
        peer_backends = validate_coverage(args.coverage, selection, provenance,
                                          attempts + calibration_attempts)
    except ValueError as error:
        parser.error(str(error))
    if any(a['level'] != selection['levels'][a['family']] for a in attempts):
        parser.error('actual level differs from frozen external-only selection')
    if args.helper_memory_lower_bound:
        try:
            validate_helper_lower_bounds(attempts)
        except ValueError as error:
            parser.error(str(error))
    rows = summarize(attempts)
    if any(r['ok'] != 5 or r['attempts'] != 5 or r['n_min'] != r['n_max']
           or r['parent_rss_kib_count'] != 5 or r['process_tree_cpu_s_count'] != 5
           for r in rows):
        parser.error('five successful trials per group required')
    fronts = frontier(rows, peer_backends, args.helper_memory_lower_bound)
    args.output.mkdir(parents=True, exist_ok=True)
    for filename, data in [('summary.csv', rows), ('fronts.csv', fronts)]:
        with (args.output/filename).open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    (args.output/'report.md').write_text(describe(attempts, rows, fronts, args.calibration, args.inputs,
        provenance[0]['source_revision'], calibration_attempts, peer_backends, args.helper_memory_lower_bound))
    (args.output/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    (args.output/'calibration_selection.json').write_text(json.dumps(selection, indent=2)+'\n')
    receipts = [p/name for p in args.inputs for name in ('attempts.jsonl', 'provenance.json')]
    receipts += [args.calibration, args.calibration.parent/'attempts.jsonl',
                 args.calibration.parent/'provenance.json']
    rejections = [{k: a[k] for k in ('family', 'backend', 'level', 'threads',
                                    'trial', 'status', 'process_wall_s')}
                  for a in calibration_attempts if a['status'] != 'ok'
                  or a['returncode'] != 0 or a['process_wall_s'] >= .5]
    recipe = dict(inputs=[str(p) for p in args.inputs], calibration=str(args.calibration),
        calibration_attempts=len(calibration_attempts), calibration_rejections=rejections,
        files=[dict(path=str(p), sha256=digest(p)) for p in receipts],
        cpu_scope='completed wait4 workflow, no helper double counting',
        memory_scope='single-process VmHWM; IPC unknown aggregate, parent-only safe lower bound',
        coverage=args.coverage, required_external_backends=list(peer_backends),
        helper_memory_lower_bound=args.helper_memory_lower_bound,
        helper_memory_scope="single-process lower bound only; aggregate IPC peak remains unknown",
        assessment=f'all {len(peer_backends)*5} external configs; strict min/max time OR memory exclusion')
    (args.output/'recipe.json').write_text(json.dumps(recipe, indent=2)+'\n')
    if args.plots:
        render(rows, args.output, .9, 'completed',
               'see range assessments; no global admission')
    print(json.dumps(dict(attempts=len(attempts), groups=len(rows),
                         assessments=collections.Counter(r['status'] for r in fronts))))


if __name__ == '__main__':
    main()
