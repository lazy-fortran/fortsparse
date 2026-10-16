from pathlib import Path
import subprocess,shlex,shutil,json,hashlib,re,datetime,argparse
p=argparse.ArgumentParser()
for name in ['archive','baseline-source','candidate-source','baseline-build','candidate-build','observer','qualification']:
    p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();out=a.archive.resolve();out.mkdir(parents=True,exist_ok=True)
repo=a.candidate_source.resolve();baseline_repo=a.baseline_source.resolve()
builds={'baseline':a.baseline_build.resolve(),'candidate':a.candidate_build.resolve()}
shutil.copyfile(Path(__file__).resolve().parent/'advice_counter.c',out/'advice_counter.c')
h=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
m={'scope':'Frozen bounded ordinary20-child comparison; no run performed by freeze.py','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'groups':[dict(family='drift-real',level=128,variant='training',operation='N',backend='native-amd-low',workload=w,host_memory_limit_bytes=0) for w in ['rhs8','refactor']],'trials':5,'ordinary_deadline_wall_s':.9,'diagnostic_deadline_wall_s':3,'affinity':[30],'threads':{'OMP_NUM_THREADS':1,'OPENBLAS_NUM_THREADS':1},'builds':{}}
for variant,build in builds.items():
 d=out/'frozen'/variant;d.mkdir(parents=True,exist_ok=True)
 src=build/'bin/bench_native_families'
 if (d/src.name).exists(): (d/src.name).chmod(0o755)
 shutil.copyfile(src,d/src.name);(d/src.name).chmod(0o555)
 cache=(build/'CMakeCache.txt').read_text();shutil.copyfile(build/'CMakeCache.txt',d/'CMakeCache.txt')
 flags={line.split(':',1)[0]:line.split('=',1)[1] for line in cache.splitlines() if re.match(r'(CMAKE_(BUILD_TYPE|C_FLAGS_RELEASE|Fortran_FLAGS_RELEASE|C_COMPILER|Fortran_COMPILER)|FORTSPARSE_ENABLE_.*):',line)}
 ldd=subprocess.check_output(['ldd',str(src)],text=True)
 closure={str(Path(path).resolve()):h(Path(path).resolve()) for path in re.findall(r'(?:=> )?(/\S+)',ldd) if Path(path).is_file()}
 link=subprocess.check_output(['ninja','-C',str(build),'-t','commands','bench_native_families'],text=True).splitlines()[-1]
 core=build/'lib/libfortsparse.a'
 m['builds'][variant]={'source_revision':subprocess.check_output(['git','-C',str(repo if variant=='candidate' else baseline_repo),'rev-parse','HEAD'],text=True).strip(),'binary':str(d/src.name),'binary_sha256':h(src),'core_archive':str(core),'core_sha256':h(core),'cache_sha256':h(build/'CMakeCache.txt'),'flags':flags,'library_closure':closure,'ldd':ldd,'link_command':link}
 if variant=='candidate':
  obj=out/'advice_counter.o';subprocess.run(['cc','-O2','-c',str(out/'advice_counter.c'),'-o',str(obj)],check=True)
  command=shlex.split(link)[2:-2]
  command[command.index('-o')+1]=str(d/'bench_native_families_diagnostic')
  command=[x for x in command if not x.startswith('-Wl,--dependency-file=')]
  command.insert(command.index('lib/libfortsparse.a'),str(obj))
  command += ['-Wl,--wrap=fsparse_native_lu_discard_payload','-Wl,--wrap=madvise']
  subprocess.run(command,cwd=build,check=True)
  (d/'bench_native_families_diagnostic').chmod(0o555)
  m['diagnostic']={'binary':str(d/'bench_native_families_diagnostic'),'binary_sha256':h(d/'bench_native_families_diagnostic'),'core_sha256':h(core),'link_command':command,'scope':'same candidate driver objects and core archive, extra syscall counter object; all instrumented timing/CPU/RSS excluded','counter_source_sha256':h(out/'advice_counter.c'),'counter_object_sha256':h(obj)}
assert m['builds']['baseline']['flags']==m['builds']['candidate']['flags']
assert m['builds']['baseline']['library_closure']==m['builds']['candidate']['library_closure']
for relative in ['bench/native_families/fixtures.f90','bench/native_families/bench_native_families.f90','bench/native_families/records.f90','bench/process_metrics.c','bench/fixtures/neo2-mini5/matrix.mtx']:
 assert h(repo/relative)==h(baseline_repo/relative)
m['source_input_sha256']={r:h(repo/r) for r in ['bench/native_families/fixtures.f90','bench/native_families/bench_native_families.f90','bench/native_families/records.f90','bench/process_metrics.c','bench/fixtures/neo2-mini5/matrix.mtx']}
qualified=a.observer.resolve()
q=a.qualification.resolve()
assert h(qualified)==json.loads(q.read_text())['observer_sha256']
shutil.copyfile(qualified,out/'measure');(out/'measure').chmod(0o555)
shutil.copyfile(qualified.with_suffix('.c'),out/'measure.c');shutil.copyfile(q,out/'observer-qualification.json')
m['observer']={'binary':str(out/'measure'),'sha256':h(out/'measure'),'qualification_receipt_sha256':h(q),'qualification':'reused exact qualified ELF; child and observer affinity validated by observer=30; original qualification includes noninflated-RSS and watchdog/reaping oracle'}
m['freeze_script_sha256']=h(__file__)
(out/'manifest.json').write_text(json.dumps(m,indent=2)+'\n')
print('Frozen both Release ELFs, matched flags/closure, exact-core diagnostic ELF; no numerical launches')
