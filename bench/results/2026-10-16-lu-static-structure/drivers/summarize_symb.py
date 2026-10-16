import sys, csv
for path in sys.argv[1:]:
    rows = [l.strip().split(',') for l in open(path) if l.startswith('CSV') or l.startswith('HDR')]
    hdr = rows[0]; data = [dict(zip(hdr[1:], r[1:])) for r in rows[1:]]
    print(f"## {path}")
    print("| matching | order | zero diag | nnz(L) | nnz(L+U) sym | LU flops | snodes given/post (exact) | snodes post relaxed | max w | max front | pad | payload MB | flops w>=8 | flops w>=16 | GEMM share | GEMM k>=8 | hist w1,2-3,4-7,8-15,16-31,32-63,64-128 |")
    print("|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    keys = []
    for d in data:
        k = (d['variant'], d['ordering'])
        if k not in keys: keys.append(k)
    for k in keys:
        g = {(d['order'], d['relaxed']): d for d in data if (d['variant'], d['ordering']) == k}
        pr = g[('postordered', 'T')]
        h = ','.join(pr[c] for c in ['w1','w2-3','w4-7','w8-15','w16-31','w32-63','w64-128'])
        print(f"| {k[0]} | {k[1]} | {pr['zero_diag']} | {int(pr['nnzL']):,} | {int(pr['nnzLU_sym']):,} | {float(pr['flopsLU_sym']):.3g} | "
              f"{g[('given','F')]['nsuper']}/{g[('postordered','F')]['nsuper']} | {pr['nsuper']} | {pr['maxwidth']} | {pr['maxfront']} | "
              f"{float(pr['pad_overhead']):.3f} | {int(pr['payload_bytes'])/1e6:.1f} | {float(pr['frac_flops_w>=8']):.3f} | {float(pr['frac_flops_w>=16']):.3f} | "
              f"{float(pr['gemm_frac']):.3f} | {float(pr['gemm_frac_k>=8']):.3f} | {h} |")
