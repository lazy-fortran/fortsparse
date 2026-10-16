# Witnessed LU reach-graph pruning, 2026-10-16

Conservative reach-graph pruning reduces native LU factor/refactor time on the
frozen drift cases. All 80 public paired attempts and ten private NEO-2 attempts
pass. The public median complete-process wall and completed CPU improve in all
eight workload/variant groups; some complex observed ranges overlap. This
checkpoint establishes a measured optimization, without external-frontier or
novel-algorithm admission. The separate full five-peer campaign owns that gate.

The candidate `0cb6c5c` contains pruning source `7db6a1a` and the same benchmark
schema as baseline `3b9384b`. Pruning was integrated as `7e1f123`. Both frozen
ELFs resolve the same system AMD/config/OpenBLAS libraries and have no native
RUNPATH or `LD_LIBRARY_PATH` override. Five fresh-process trials are interleaved
per configuration on CPU 0 of a Ryzen 9 5950X, with one actual factor thread,
identical thread environment and a 0.9 s child watchdog. Quiet windows exclude
other task builds and measurements; they do not imply cold filesystem caches.

## Public measurements

Times below are milliseconds, median [minimum, maximum]. Wall covers the complete
subprocess. The phase is factorization for RHS8 and numeric refactorization for
refactor; the latter complete workflow also includes initial factorization and
one solve. All operations are N. Real n=16,384; complex n=4,096. Training and
held-out coefficients use the frozen externally selected sizes.

| Case | Complete wall, baseline → pruned | Factor/refactor phase, baseline → pruned |
|---|---|---|
| Real training RHS8 | 66.96 [66.41,69.03] → 54.74 [53.68,58.30] | 47.52 [46.95,48.25] → 35.72 [34.18,37.91] |
| Real training refactor | 95.48 [93.42,98.83] → 69.09 [68.56,73.81] | 38.47 [38.20,40.10] → 25.85 [25.81,27.89] |
| Real held-out RHS8 | 67.55 [66.28,71.29] → 54.07 [53.24,56.18] | 47.11 [46.86,48.51] → 34.48 [33.86,35.93] |
| Real held-out refactor | 94.71 [93.62,97.56] → 70.02 [68.91,70.92] | 38.57 [38.45,41.16] → 26.35 [26.13,26.93] |
| Complex training RHS8 | 22.81 [22.31,24.36] → 21.38 [20.74,22.02] | 11.98 [11.58,12.41] → 10.57 [10.39,11.32] |
| Complex training refactor | 25.36 [25.23,25.93] → 23.59 [23.15,23.68] | 8.69 [8.57,8.76] → 7.72 [7.63,7.78] |
| Complex held-out RHS8 | 22.33 [21.31,23.49] → 21.50 [20.31,22.58] | 11.99 [11.42,12.61] → 10.60 [10.44,11.39] |
| Complex held-out refactor | 25.98 [25.14,26.23] → 23.75 [23.39,26.73] | 8.72 [8.52,8.94] → 7.73 [7.68,9.92] |

Completed CPU is reaped `wait4` user+system for the child lifetime, with no helper.
RSS below is the driver's pre-CSV Linux VmHWM sample, not a sampled aggregate or
an exact solver-allocation peak. Each entry is a median; all ranges remain in
[`summary.csv`](summary.csv), including the complex held-out refactor outlier.

| Case | Completed CPU ms, baseline → pruned | VmHWM KiB, baseline → pruned |
|---|---|---|
| Real training RHS8 | 65.88 → 53.70 | 30,040 → 29,636 |
| Real training refactor | 93.84 → 67.84 | 24,068 → 23,524 |
| Real held-out RHS8 | 66.21 → 52.86 | 29,668 → 29,872 |
| Real held-out refactor | 93.15 → 68.79 | 23,500 → 23,284 |
| Complex training RHS8 | 21.82 → 20.48 | 12,480 → 12,488 |
| Complex training refactor | 24.71 → 22.82 | 11,656 → 11,660 |
| Complex held-out RHS8 | 21.61 → 20.01 | 12,500 → 12,468 |
| Complex held-out refactor | 25.00 → 22.88 | 11,620 → 11,656 |

There is no demonstrated general RSS reduction. Retained factor fill is identical
in every matched pair: 762,644 real and 134,400 complex entries. The additional
signed graph endpoint costs exactly 4*n owned bytes: 65,536 and 16,384 bytes,
respectively. Cap preflight includes old/new metadata overlap during reuse;
reference execution retains no endpoint array. Negative/unavailable attributed
peak sentinels become blank CSV fields, never zero. Caller storage, runtime
pages and allocator behavior remain distinct from the attributed live payload.

Maximum public backward error is 1.051e-15 and forward error 2.679e-15, against
unchanged gates of 1e-12 and 1e-10. No failed attempts or outliers were discarded.

## Private consumer aggregate

The frozen `neo2-lorentz-tag6-4856-v1` fixture has n=4,856 and nnz=28,898, with
three physical RHS and one manufactured RHS. All five paired native-AMD trials
pass both the process checks and independent exact audit. Numeric matrix/vector
data, candidates and per-RHS fractions remain in the private archive.

| Metric | Baseline median [min,max] | Pruned median [min,max] |
|---|---|---|
| Factor ms | 10.106 [10.039,10.969] | 9.058 [8.967,9.928] |
| Complete child wall ms | 88.244 [87.886,91.389] | 89.023 [87.699,89.926] |
| Completed child CPU ms | 86.594 [86.351,89.849] | 87.252 [86.000,88.271] |
| VmHWM KiB | 11,416 [11,392,11,444] | 11,704 [11,416,11,716] |
| Child plus exact-audit wall ms | 623.014 [619.350,638.653] | 626.178 [621.298,645.602] |

Factor fill remains 192,971 entries; live owned storage rises from 3,533,268 to
3,552,692 bytes, exactly 4*n. The factor phase improves with separated ranges;
whole-child wall/CPU and RSS ranges overlap, and their medians increase. No
precise private end-to-end or memory improvement is established. Child wall
includes matrix/RHS input, checking, export and teardown; the independent exact
audit is separate. Only the last row includes that audit. The frozen fixture,
policy and fresh five-peer qualification receipt hashes are in the recipe.

## Mechanism, proof obligation and executable oracles

The numerical L span is retained in full. DFS visits a prefix of each column.
After pivot j, an old edge k→r is excluded from that prefix only when actual
retained edges k→j and j→r witness a length-two path. The witness edge k→j is
kept. Row/value pairs move together; numeric updates and solves use the full span.

Let G be the full retained numerical graph and H the DFS graph. Initially they
coincide. Adding a newly pivoted vertex's outgoing edges to both preserves equal
transitive closure. Removing k→r from H preserves its closure because k→j→r
remains in H. Induction gives TC(H)=TC(G), preserving reach sets under dynamic
pivoting. A canceled numerical coefficient is absent and cannot witness removal.
This is a conservative application of established sparse-LU pruning principles;
it is not a new-algorithm claim. See Demmel et al., *A Supernodal Approach to
Sparse Partial Pivoting*, [SIAM J. Matrix Anal. Appl. 20 (1999), §4](https://portal.nersc.gov/project/sparse/xiaoye-web/simax-29176.pdf),
and retained SuperLU `dpruneL.c`/`dcolumn_dfs.c` provenance below.

The actual generated helper is checked on all 1,024 edge sets with a fixed five-node
topological order against independent Boolean Floyd–Warshall closure after every partition, with permuted
original-row tags and aligned real/complex values. Additional oracles exercise
exact cancellation that invalidates blind symmetric-pair pruning, dynamic pivots,
row/column permutations, N/T/H block solves, real/complex implicit JVP/VJP,
refinement, failed-pivot recovery and exact metadata cap boundaries. Candidate
focused CTest passed 9/9 with bounds checking; four later benchmark-schema,
public-pruning and memory-alias gates passed. No tests were rerun for this
docs-only checkpoint.

The retained source is `gen/lu_factor.f90.inc`; the scalar FortSym-generated
primal/JVP/VJP leaf files and FortSym pin are unchanged. Source generation was
repeated with identical output. The simple unpruned reference remains available.

## Profile and reproduction scope

[`profiles.csv`](profiles.csv) records exact sample counts from `perf` stderr,
approximate event counts and self sample shares. The before/after DFS shares are
37.17→15.84% (real RHS8), 48.71→20.82% (real refactor), 18.82→8.13%
(complex n=4,096 RHS8), and 30.51→10.48% (complex n=16,384 refactor diagnostic).
Pruning itself takes 0.93–1.55% of sampled cycles. Nonexclusive instrumented
replays differ in count and source benchmark metadata; their counts and shares
support diagnosis, not ordinary timing or a cycle-ratio speedup claim.

[`recipe.json`](recipe.json) retains sources, ELF/library hashes, raw-receipt
hashes, environment, commands, affinity, interleaving and profile provenance.
The original A/B runner and replay scripts remain in the named external archive.
Reproduce ordinary runs only in a controller-approved quiet window. Preserve
system AMD/config package pins; do not override the direct UMFPACK helper's
separate library closure with a global `LD_LIBRARY_PATH`.

Rebuild these summaries without executing any solver or copying private data:

```sh
python3 bench/results/native-families-2026-10-16/witnessed-pruning/summarize.py \
  --public /home/ert/code/sparse-dev/benchmarks/2026-10-16-native-families/witnessed-pruning-ab \
  --private /mnt/storage/data2/sparse-dev-private/neo2-lorentz-tag6-2026-10-16/pruning-ab-3b9384b-0cb6c5c \
  --output /tmp/fortsparse-pruning-summary
```

Retained literature is under `sparse-dev/papers/core_solvers/` and its library
cards: `superlu_serial_1999` and `sparse_direct_survey_davis_2016`. Inspected
SuperLU source revision: `a42ca6c57b2888e005cf3e1c11d44dec4ad49668`.
