# Final single-node SPD measurements, 2026-10-16

The retained final CPU corpus has 17,250 successful fresh-process trials:
11 synthetic cases in training and held-out variants, plus the fixed
HB/bcsstk09 structural stiffness matrix; five workloads, six backend/ordering
configurations, five physical thread counts and five trials per configuration.
Native CPU measurements use `d91affa`; retained external and CUDA measurements
use `19af66a`. Full source/binary/build identities are in the plot manifest.
The eight original raw JSONs are archived outside Git under
`~/code/sparse-dev/benchmarks/2026-10-16-native-mvp/`.

For all 115 case/workload/variant groups, at least one explicitly measured
production native CPU ordering/thread configuration is robustly nondominated
in elapsed solver phases versus parent peak RSS over the observed ranges.
A separate analysis using whole-workflow CPU time versus RSS also finds
115/115 robust front members. The selected configuration can differ between
these two objectives. Whole-workflow CPU includes parent plus helper user/system
time, setup, warmup factors, checks and cleanup; it is not phase CPU time.

The stronger same-comparator match/beat gate remains distinct: 109 groups
pass, four are indeterminate (band level 1 refactor and grid2 level 7 cold
one-RHS, each in both variants), and two retain a time/memory tradeoff
(grid2 level 7 cold eight-RHS, both variants). These two strict-unmet groups
are nondominated; they are not reported as being off the Pareto front.
Case-specific ordering and thread tuning is required; this is not a universal
AUTO/default claim or evidence of a novel factorization algorithm.

UMFPACK aggregate peak memory remains unknown. Its parent/helper peaks are
separate, and the parent is only a lower bound for dominance exclusion. The
300 CUDA trials all pass on actual device 0, but these bounded cases show no
GPU advantage. Host RSS and attributed device bytes are separate quantities.
The figures also show limited strong scaling, actual teams smaller than some
requested counts, and increasing workflow CPU consumption at high thread counts.
Observed ranges are not confidence intervals or a claim about other machines.

The figures plot training fixtures only; held-out results remain separately
identified in the complete case CSV and admission report.

Retained evidence:

- `2026-10-16-final-cases.csv`: all 115 groups, both gate statuses, selected
  native/comparator time, CPU and memory ranges, ordering, teams and identities.
- `2026-10-16-final-plot-inputs.csv`: 372 configurations used in the figures,
  with all plotted ranges, trial counts, identities and device verification.
- `2026-10-16-final-front-audit.json`: compact, separately labeled elapsed/RSS
  and whole-workflow CPU/RSS audit counts.
- `2026-10-16-final-plot-manifest.json`: raw/compact data and image hashes,
  source/build provenance, publication URLs and inspection/replay evidence.

Reproduce all three scientific figures directly from the retained numeric CSV:

```sh
python3 bench/plot_spd.py bench/results/2026-10-16-final-plot-inputs.csv \
  --outdir /var/tmp/fortsparse-mvp/plots/final-replay --workload warm \
  --cases grid2:7,grid3:7 --scaling-cases grid2:7,grid3:7 \
  --gpu-cases grid2:3,grid3:7,bcsstk09:1 \
  --title 'FortSparse final single-node SPD evidence'
```

Recompute the report from archived raw files, replacing only native series:

```sh
SPDBENCH_ARCHIVE="$HOME/code/sparse-dev/benchmarks/2026-10-16-native-mvp"
python3 bench/report_spd.py \
  "$SPDBENCH_ARCHIVE"/final-{training,heldout,bcsstk09}.json \
  "$SPDBENCH_ARCHIVE"/final-native-{training,heldout,bcsstk09}.json \
  --output /var/tmp/fortsparse-mvp/final-report.json
```

The PNGs are published through slopbox and expire three days after upload;
their opaque URLs and hashes are retained in the manifest. Images and large
raw/configuration-level report files are deliberately outside source history.
