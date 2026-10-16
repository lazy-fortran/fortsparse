#!/usr/bin/env python3
"""Regenerate portable numeric tables from this checkpoint's immutable receipts."""
import argparse,collections,csv,hashlib,json,statistics
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('archive',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
here=Path(__file__).resolve().parent
provenance=json.loads((here/'provenance.json').read_text())
for name in ['attempts.jsonl','diagnostics.jsonl','manifest.json']:
    assert hashlib.sha256((a.archive/name).read_bytes()).hexdigest()==provenance['raw_sha256'][name],name
ordinary=[json.loads(line) for line in (a.archive/'attempts.jsonl').read_text().splitlines()]
diag=[json.loads(line) for line in (a.archive/'diagnostics.jsonl').read_text().splitlines()]
manifest=json.loads((a.archive/'manifest.json').read_text())
assert len(ordinary)==20 and len(diag)==2
assert {v:i['source_revision'] for v,i in manifest['builds'].items()}==provenance['measured_sources']
a.output.mkdir(parents=True,exist_ok=True)
def write(name,rows):
    with (a.output/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
trials=[]
for r in ordinary:
    z={k:r[k] for k in ['workload','version','trial','status','returncode','binary_sha256']}
    z.update({k:r['process'][k] for k in ['wall_s','cpu_s','user_s','system_s','maxrss_kib','exit_status','cpu_affinity']})
    z.update({k:r['record'][k] for k in ['backend','n','nnz','nrhs','factor_s','solve_s','refactor_s','reused','factor_nnz','host_current_bytes','host_peak_bytes','backward_error','forward_error','memory_policy','host_memory_limit_bytes','parent_rss_kib']})
    trials.append(z)
write('trials.csv',trials)
summary=[]
for workload in ['rhs8','refactor']:
    for version in ['baseline','candidate']:
        sample=[r for r in ordinary if r['workload']==workload and r['version']==version]
        good=[r for r in sample if r['status']=='ok']
        z={'workload':workload,'version':version,'source_revision':manifest['builds'][version]['source_revision'],
           'binary_sha256':manifest['builds'][version]['binary_sha256'],'attempts':len(sample),'ok':len(good),
           'n':16384,'nnz':81408,'threads_requested':1,'active_factor_threads':1,'host_cap_bytes':0}
        for metric in ['wall_s','cpu_s','maxrss_kib','factor_s','solve_s','refactor_s','parent_rss_kib','factor_nnz','host_current_bytes','host_peak_bytes']:
            vals=[float(r['process'][metric] if metric in r['process'] else r['record'][metric]) for r in good]
            vals=[v for v in vals if v>=0]
            values=[statistics.median(vals),min(vals),max(vals)] if vals else ['','','']
            for name,value in zip(['median','min','max'],values):z[metric+'_'+name]=value
        summary.append(z)
write('summary.csv',summary)
counts=[]
for r in diag:
    z={'workload':r['workload'],'status':r['status'],'returncode':r['returncode'],
       'diagnostic_binary_sha256':r['binary_sha256'],'watchdog_wall_s':3,'timing_metrics_excluded':True}
    z.update({k:v for k,v in r['advice'].items() if k!='errno_counts'})
    z['errno_counts']=json.dumps(r['advice']['errno_counts'],sort_keys=True)
    counts.append(z)
write('diagnostics.csv',counts)
print('Regenerated20ordinaryrows,4aggregategroupsand2excludeddiagnosticrows')
