#!/usr/bin/env python3
"""Read frozen pruned-cohort data and regenerate reports; never execute a solver."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

REPO = Path(__file__).resolve().parents[4]
SCRIPTS = REPO/'bench/native_families'
BASE = Path('/home/ert/code/sparse-dev/benchmarks/2026-10-16-native-families')
RAW_SHA256 = 'bbb373762d0b1347c5ee1c97c3c8de703dff4f9e1add4b11427c774d9f4acda1'
PROVENANCE_SHA256 = 'aeb76043d2c13d70e9035ad023bc5f51e0b1b7f7d97fa88fd3250a13c96ee2fa'
CALIBRATION_SHA256 = 'cb39c751f3832a279e3e378ef33b47fc3cf9c714ee5128e3aedf380e790a42dc'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=BASE/'pruned-five-peer-drifts')
    parser.add_argument('--calibration', type=Path,
        default=BASE/'openblas-calibration/calibration_selection.json')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plots', action='store_true')
    args = parser.parse_args()
    if args.output.resolve().is_relative_to(REPO):
        parser.error('generated artifacts must stay outside the source checkout')
    if digest(args.input/'attempts.jsonl') != RAW_SHA256 or digest(
            args.input/'provenance.json') != PROVENANCE_SHA256:
        parser.error('raw cohort or source/binary provenance differs from this checkpoint')
    selection = json.loads(args.calibration.read_text())
    if digest(args.calibration) != CALIBRATION_SHA256 or digest(
            args.calibration.parent/'attempts.jsonl') != selection['attempts_sha256']:
        parser.error('historical external-only size gate differs from this checkpoint')
    command = [sys.executable, str(SCRIPTS/'report_front.py'), str(args.input),
               '--calibration', str(args.calibration), '--output', str(args.output)]
    if args.plots:
        command.append('--plots')
    subprocess.run(command, check=True)
    with (args.output/'summary.csv').open() as handle:
        rows = list(csv.DictReader(handle))
    scaling = []
    keys = ('family', 'backend', 'level', 'workload', 'variant', 'operation')
    for row in rows:
        baseline = next(r for r in rows if r['threads'] == '1' and
                        all(r[k] == row[k] for k in keys))
        speedup = float(baseline['process_wall_s_median'])/float(row['process_wall_s_median'])
        scaling.append({k: row[k] for k in (*keys, 'threads')} | dict(
            wall_speedup_vs_own_1thread_median=speedup,
            requested_budget_wall_efficiency=speedup/int(row['threads']),
            active_factor_threads_min=row['active_threads_min'],
            active_factor_threads_max=row['active_threads_max'],
            active_factor_threads_samples=row['active_threads_count']))
    with (args.output/'requested-scaling.csv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(scaling[0]))
        writer.writeheader()
        writer.writerows(scaling)
    (args.output/'reproduction.json').write_text(json.dumps(dict(
        coverage='historical-five', solver_measurement_source_revision=
        '7e1f123ea21fd8412480ec6a6ddf0696893acb2a',
        input=str(args.input), raw_sha256=RAW_SHA256, provenance_sha256=PROVENANCE_SHA256,
        report_sources={name: digest(SCRIPTS/name) for name in ('report.py','report_front.py')},
        recipe_source_sha256=digest(Path(__file__)),
        scope='only pruned drift receipts; historical calibration is a size gate, not joined chain measurements',
        scaling='backend own t1median/tPmedian; efficiency /requested budget, not active-team efficiency'
        ), indent=2)+'\n')


if __name__ == '__main__':
    main()
