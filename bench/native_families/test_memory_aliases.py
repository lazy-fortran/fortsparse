#!/usr/bin/env python3
"""Tiny actual-driver CSV/configuration oracles; no performance admission."""
import os
from pathlib import Path
import subprocess
import sys
from run import parse_csv, record_memory_configuration

binary=Path(sys.argv[1])
mini=Path(sys.argv[2])
numeric=unavailable=0
for backend in ('native-amd-low','native-amd-capped'):
    for family in ('drift-real','drift-complex'):
        env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
        env.pop('FORTSPARSE_BENCH_HOST_LIMIT_BYTES',None)
        cap=((2**53+1 if family=='drift-real' else 2**63-1)
             if backend=='native-amd-capped' else None)
        if cap:
            env['FORTSPARSE_BENCH_HOST_LIMIT_BYTES']=str(cap)
        result=subprocess.run([str(binary),backend,family,'3','1','refactor','H','heldout',str(mini)],
                              capture_output=True,text=True,env=env,timeout=.4)
        row=parse_csv(result.stdout)
        if not row:
            raise AssertionError((result.returncode,result.stdout,result.stderr))
        if row['status']=='backend_unavailable':
            unavailable+=1
            continue  # An explicitly unavailable AMD build is not a passing numeric sample.
        numeric+=1
        assert result.returncode==0 and row['status']=='ok',(result.stdout,result.stderr)
        receipt=record_memory_configuration(dict(status='ok',record=row),backend,cap)
        assert receipt['status']=='ok',receipt
        assert row['backward_error']<1e-12 and row['forward_error']<1e-10
        assert row['host_peak_bytes']>0
        print(backend,family,'memory_policy='+str(row['memory_policy']),
              'host_memory_limit_bytes='+str(row['host_memory_limit_bytes']))
for cap in (None,'0','-1',str(2**63),'123junk'):
    env=dict(os.environ)
    env.pop('FORTSPARSE_BENCH_HOST_LIMIT_BYTES',None)
    if cap is not None:
        env['FORTSPARSE_BENCH_HOST_LIMIT_BYTES']=cap
    result=subprocess.run([str(binary),'native-amd-capped','drift-real','3','1',
                           'rhs8','N','training',str(mini)],capture_output=True,text=True,
                          env=env,timeout=.4)
    assert result.returncode!=0 and parse_csv(result.stdout) is None,'invalid cap accepted'
print(f'actual memory aliases: {numeric} numeric successes, {unavailable} unavailable; exact i8/configuration and invalid-cap guards pass')
