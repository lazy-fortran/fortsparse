from pathlib import Path
import csv,json,statistics,hashlib
root=Path(__file__).resolve().parent
fields=['experiment','source','control_source','case','pair','width','check_ok','timeout','outer_timeout_returncode',
        'wall_s','user_s','system_s','cpu_s','maxrss_kib','factor_wall_s','solve_wall_s','factor_nnz',
        'owned_current_bytes','owned_peak_bytes','normwise_be','componentwise_be','actual_order_sha256']
rows=[]
for exp,folder,source in [('arbitrary-mask','lu-panel8-ab','2926aad'),('uniform-full-mask','lu-panel8-fullmask-ab','0598590'),('contiguous-active-run','lu-panel8-range-ab','bfc607b')]:
    sourcepath=root.parent/folder/'ordinary-attempts.json'
    for r in json.loads(sourcepath.read_text()):
        rec={k:r[k] for k in fields if k in r}
        rec.update(experiment=exp,source=source,control_source='2926aad',case=r['label'])
        rows.append(rec)
with (root/'trials.csv').open('w') as stream:
    writer=csv.DictWriter(stream,fields);writer.writeheader();writer.writerows(rows)
summary=[]
for exp in ['arbitrary-mask','uniform-full-mask','contiguous-active-run']:
 for case in ['kin6d-5525','drift-real-32']:
  for width in [1,8]:
   group=[r for r in rows if (r['experiment'],r['case'],r['width'])==(exp,case,width)]
   for metric in ['wall_s','cpu_s','maxrss_kib','factor_wall_s','factor_nnz','owned_current_bytes','owned_peak_bytes']:
    values=[r[metric] for r in group]
    summary.append(dict(experiment=exp,case=case,width=width,metric=metric,n=len(values),
                        median=statistics.median(values),minimum=min(values),maximum=max(values)))
with (root/'summary.csv').open('w') as stream:
    writer=csv.DictWriter(stream,['experiment','case','width','metric','n','median','minimum','maximum'])
    writer.writeheader();writer.writerows(summary)
print('Retained',len(rows),'requested ordinary trials; passes',sum(r['check_ok'] for r in rows))
