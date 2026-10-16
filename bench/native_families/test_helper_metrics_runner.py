"""Actual IPC runner diagnostic propagation; correctness smoke only."""
import os
from pathlib import Path
import sys
from helper_memory import ENVIRONMENT, parse_helper_memory
from run import invoke
from report import summarize

binary, mini5 = (Path(p).resolve() for p in sys.argv[1:])
core = min(os.sched_getaffinity(0))
for workload in ('rhs8', 'warm'):
    args = ['umfpack', 'drift-real', '3', '1', workload, 'T', 'heldout', str(mini5)]
    receipt = invoke(binary, args, 1, [core], .8)
    assert receipt['status'] == 'ok', receipt
    assert receipt['environment'][ENVIRONMENT] == '1', receipt
    row = receipt['record']
    assert row['backward_error'] < 1e-12 and row['forward_error'] < 1e-10, row
    if workload == 'rhs8':
        assert receipt['helper_memory_status'] == 'available', receipt
        assert receipt['helper_postexec_hwm_kib'] > 0, receipt
        assert row['helper_postexec_hwm_kib'] == parse_helper_memory(
            receipt['stderr'])['helper_postexec_hwm_kib'], receipt
    else:
        assert receipt['helper_memory_status'] == 'invalid_receipt', receipt
        assert row['helper_postexec_hwm_kib'] == -1, receipt
    receipt.update({key: row[key] for key in ('backend','family','level','threads',
                   'workload','variant','operation')})
    summary = summarize([receipt])[0]
    assert summary['helper_postexec_hwm_kib_count'] == (1 if workload == 'rhs8' else 0)
print('IPC postexec memory propagation, warm ambiguity and unchanged numerical checks pass')
