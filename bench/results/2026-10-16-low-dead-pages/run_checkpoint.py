"""Frozen20-child ordinary comparison plus2 separately excluded diagnostics."""
from pathlib import Path
import argparse,collections,datetime,hashlib,importlib.util,json,os,statistics,subprocess,sys,time
p=argparse.ArgumentParser()
p.add_argument('--archive',type=Path,required=True)
p.add_argument('--fixture-source',type=Path,required=True)
p.add_argument('--go',action='store_true',required=True)
a=p.parse_args();out=a.archive.resolve();repo=a.fixture_source.resolve()
m=json.loads((out/'manifest.json').read_text())
for relative,digest in m['source_input_sha256'].items():
    assert hashlib.sha256((repo/relative).read_bytes()).hexdigest()==digest
sys.path.insert(0,str(repo/'bench/native_families'))
from run import parse_csv
h=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for v,info in m['builds'].items():
    info['binary']=str(out/'frozen'/v/'bench_native_families')
    assert h(info['binary'])==info['binary_sha256']
m['diagnostic']['binary']=str(out/'frozen/candidate/bench_native_families_diagnostic')
assert h(m['diagnostic']['binary'])==m['diagnostic']['binary_sha256']
assert h(out/'measure')==m['observer']['sha256']
env=os.environ.copy()
for name in ['LD_PRELOAD','FORTSPARSE_BENCH_HOST_LIMIT_BYTES']:env.pop(name,None)
env.update(MEASURE_CPUS='30',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',
           OMP_PLACES='cores',OMP_PROC_BIND='close',MKL_NUM_THREADS='1',MKL_DYNAMIC='FALSE')
logs=out/'logs';logs.mkdir(exist_ok=True)
start={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'loadavg':os.getloadavg(),
       'runner_sha256':h(__file__),'manifest_sha256':h(out/'manifest.json'),'environment':{k:env.get(k) for k in ['MEASURE_CPUS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','OMP_PLACES','OMP_PROC_BIND','LD_PRELOAD']}}
(out/'run-start.json').write_text(json.dumps(start,indent=2)+'\n')
rows=[]
def execute(version,workload,trial,diagnostic=False):
    binary=m['diagnostic']['binary'] if diagnostic else m['builds'][version]['binary']
    args=['native-amd-low','drift-real','128','1',workload,'N','training',str(repo/'bench/fixtures/neo2-mini5/matrix.mtx')]
    label=f'{workload}-{version}-{trial}';log=logs/(label+'.txt')
    deadline=3 if diagnostic else .9
    command=['/usr/bin/timeout','--signal=TERM','--kill-after=.2s',str(deadline)+'s',str(out/'measure'),binary,str(log),*args]
    before=os.getloadavg();result=subprocess.run(command,capture_output=True,text=True,env=env,start_new_session=True)
    text=log.read_text() if log.exists() else ''
    receipt={'version':version,'workload':workload,'trial':trial,'diagnostic':diagnostic,
        'command':command,'binary_sha256':h(binary),'watchdog_wall_s':deadline,'returncode':result.returncode,
        'observer_stdout':result.stdout,'observer_stderr':result.stderr,'output':text,'load_before':before,'load_after':os.getloadavg(),
        'status':'timeout' if result.returncode==124 else 'failed','record':parse_csv(text),
        'metric_scope':'INSTRUMENTED: ALL elapsed/CPU/RSS excluded' if diagnostic else 'ordinary native-only complete child; no helper',
        'output_sha256':h(log) if log.exists() else None}
    try:
        receipt['process']=json.loads(result.stdout)
        assert receipt['process']['cpu_affinity']=='30'
        if result.returncode==0 and receipt['process']['exit_status']==0 and receipt['record']['status']=='ok':
            r=receipt['record']
            assert r['n']==16384 and r['nnz']==81408 and r['memory_policy']==1 and r['host_memory_limit_bytes']==0
            assert r['backend']=='native-amd-low' and r['active_threads']==1
            receipt['status']='ok'
    except (ValueError,TypeError,KeyError,AssertionError) as error:
        receipt['validation_error']=str(error)
    if diagnostic:
        try:
            line=next(l for l in text.splitlines() if l.startswith('ADVICE_DIAGNOSTIC '))
            receipt['advice']=json.loads(line.split(' ',1)[1])
            assert receipt['advice']['advice_calls']>0 and receipt['advice']['invalid_interior_spans']==0
        except (StopIteration,ValueError,KeyError,AssertionError) as error:
            receipt['status']='diagnostic_invalid';receipt['diagnostic_error']=str(error)
    return receipt
with (out/'attempts.jsonl').open('x') as f:
    for trial in range(5):
        for workload in ['rhs8','refactor']:
            for version in (['baseline','candidate'] if trial%2==0 else ['candidate','baseline']):
                r=execute(version,workload,trial);rows.append(r);f.write(json.dumps(r)+'\n');f.flush()
        print('completed ordinary pair trial',trial+1,flush=True)
with (out/'diagnostics.jsonl').open('x') as f:
    for workload in ['rhs8','refactor']:
        r=execute('diagnostic',workload,0,True);f.write(json.dumps(r)+'\n');f.flush()
summary=[]
for workload in ['rhs8','refactor']:
    item={'workload':workload,'trials_per_variant':5,'versions':{},'unchanged_pairs':[]}
    for version in ['baseline','candidate']:
        sample=[r for r in rows if r['version']==version and r['workload']==workload]
        point={'attempts':len(sample),'statuses':dict(collections.Counter(r['status'] for r in sample))}
        for key in ['wall_s','cpu_s','maxrss_kib','factor_s','solve_s','refactor_s','parent_rss_kib','factor_nnz','host_current_bytes','host_peak_bytes']:
            vals=[float(r['process'][key] if key in r.get('process',{}) else r['record'][key]) for r in sample if r['status']=='ok']
            vals=[v for v in vals if v>=0]
            point[key]=dict(count=len(vals),median=statistics.median(vals),min=min(vals),max=max(vals)) if vals else None
        item['versions'][version]=point
    for trial in range(5):
        pair=[r for r in rows if r['workload']==workload and r['trial']==trial]
        keys=['factor_nnz','host_current_bytes','host_peak_bytes','backward_error','forward_error','reused']
        item['unchanged_pairs'].append({'trial':trial,'both_success':all(r['status']=='ok' for r in pair),
            'fields_equal':{k:len({r['record'][k] for r in pair if r['record']})==1 for k in keys}})
    summary.append(item)
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
(out/'run-end.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'statuses':dict(collections.Counter(r['status'] for r in rows)),'loadavg':os.getloadavg()},indent=2)+'\n')
print('ordinary statuses',dict(collections.Counter(r['status'] for r in rows)),flush=True)
