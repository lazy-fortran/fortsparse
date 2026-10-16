import json, sys
for p in sys.argv[1:]:
    d = json.load(open(p))
    head = f"{d['matching']:10s} {d['ordering']:3s} {d['scaling']:4s} {d['mode']:10s}"
    for r in d['runs']:
        ref = r['refinement']
        s0 = ref[0]; s1 = ref[1] if len(ref) > 1 else ref[0]; s3 = ref[min(3, len(ref)-1)]
        print(f"{head} tau={r['tau']:.0e} pert={r['perturbed']} swaps={r['in_supernode_swaps']} viol0.1={r['threshold0.1_violations']} "
              f"growth={r['growth_maxU_over_maxM']:.3g} maxL={r['max_abs_L']:.3g} minUdiag={r['min_abs_Udiag']:.3g} "
              f"BE0(n/c)={s0['normwise_be']:.2g}/{s0['componentwise_be']:.2g} BE1={s1['normwise_be']:.2g}/{s1['componentwise_be']:.2g} "
              f"BE3={s3['normwise_be']:.2g}/{s3['componentwise_be']:.2g} fwd3={s3['forward_rel']:.2g}")
