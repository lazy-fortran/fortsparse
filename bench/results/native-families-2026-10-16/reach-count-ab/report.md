# LU reach-count fusion, 2026-10-16

Counting nonzero unpivoted entries during the existing pivot/U-count pass removes
one complete reach traversal per column. A valid pivot contributes exactly one
such entry, so subtracting one gives the L count. Pivot comparisons, tie order,
scalar primal/JVP/VJP arithmetic, pruning, factor fill and owned bytes are
unchanged. The simpler generated schedule is retained with a bounded performance
claim: real training RHS8 complete-process wall improves in all five paired
trials. Other cases have mixed signs and overlapping ranges; this checkpoint
does not establish a uniform speedup, private end-to-end gain or external-front
admission.

Baseline `cac6656` and composed candidate `8848e2d` include identical accepted
LOW growth and optional-cap changes. Optimization commit `7fefa61` changes only
the retained LU template, its real/complex emission and focused range/recovery
oracles. Release flags are `-O3 -march=znver3`; both frozen native ELFs have no
RUNPATH and resolve identical library hashes. There is no `LD_LIBRARY_PATH`
override. Sources, ELF/library identities and durable receipts are in
[`recipe.json`](recipe.json).

## Ordinary paired measurements

All 80 public attempts pass: five fresh processes per version for each real
n=16,384/complex n=4,096, training/held-out, RHS8/refactor configuration. Each
trial alternates baseline/candidate order; affinity is CPU 0, with one
actual factor thread, fixed thread environment and a 0.9 s child watchdog.
The controller reserved an exclusive measurement window. A fresh process does
not imply a cold filesystem cache. Sizes and accuracy gates remain frozen.

The table gives complete subprocess wall and completed `wait4` user+system CPU,
in milliseconds, median [minimum, maximum]. Neither is the driver's shorter
main-entry interval. No helper executes in these native measurements. Refactor
workflow totals also include initial factorization and one RHS solve.

| Case | Complete wall, baseline → candidate | Completed CPU, baseline → candidate |
|---|---|---|
| Real training RHS8 | 57.086 [56.407,63.630] → 56.172 [54.529,60.042] | 55.155 [54.503,61.312] → 54.730 [52.923,58.353] |
| Real training refactor | 72.964 [70.171,73.652] → 71.234 [69.369,74.672] | 71.006 [68.677,71.774] → 69.405 [67.771,72.543] |
| Real held-out RHS8 | 58.517 [55.161,59.380] → 56.376 [56.113,57.626] | 56.601 [53.693,57.401] → 54.936 [53.717,55.463] |
| Real held-out refactor | 73.172 [72.006,77.054] → 72.096 [69.730,75.913] | 70.949 [70.155,72.489] → 70.134 [67.932,72.101] |
| Complex training RHS8 | 23.132 [21.755,24.626] → 22.494 [22.012,24.186] | 21.957 [20.768,22.238] → 21.131 [20.911,21.924] |
| Complex training refactor | 25.745 [24.281,26.174] → 25.387 [23.827,26.069] | 24.555 [23.253,25.050] → 24.059 [22.972,24.809] |
| Complex held-out RHS8 | 22.969 [22.357,23.157] → 22.535 [21.936,24.222] | 21.307 [21.238,22.027] → 21.226 [20.891,21.673] |
| Complex held-out refactor | 25.863 [25.522,26.932] → 26.162 [25.266,26.566] | 24.249 [23.924,25.217] → 24.382 [23.773,25.108] |

Real training RHS8 paired baseline/candidate wall ratios are 1.006–1.060,
median 1.029, with all five pairs improving. Completed CPU improves in four of
those five pairs, including one ratio 0.996. Paired ratios and ratios of group
medians differ because the median need not select the same trial. Across all
40 public pairs, wall improves in 27 and CPU in 32; seven of eight group wall
and CPU medians decrease. These counts describe the observed trials without
establishing a population confidence bound.

Complex held-out refactor is an explicit exception: median wall increases
1.16% and CPU 0.55%. Its paired signs are mixed and ranges overlap. The measured
effect is unresolved; the higher medians remain in the evidence. No additional
trials, outlier removal or favorable-subset timing was used for acceptance.

[`summary.csv`](summary.csv) retains every group's factor, refactor, solve,
main interval, complete wall/CPU and pre-CSV VmHWM ranges.
[`paired.csv`](paired.csv) retains matched-trial ratios and signs for every
available metric. Factor fill and live owned bytes are exactly equal in all
40 pairs: fill is 762,644 real and 134,400 complex entries. RSS varies between
processes and establishes no general memory improvement; unknown attributed
peak fields remain blank. Maximum backward/forward errors are 1.051e-15 and
2.679e-15 against unchanged 1e-12 and 1e-10 gates.

## Private aggregate regression

The frozen `neo2-lorentz-tag6-4856-v1` fixture has n=4,856, nnz=28,898 and four
RHS (three physical, one manufactured). All ten interleaved native-AMD attempts
and independent exact audits pass, following fresh five-peer qualification.
Only scalar aggregate metrics and provenance hashes leave the private archive;
no numeric inputs, candidates or per-RHS exact fractions are included here.

| Metric | Baseline median [min,max] | Candidate median [min,max] |
|---|---|---|
| Factor ms | 9.298 [9.224,9.415] | 9.014 [8.932,9.551] |
| Complete child wall ms | 90.052 [89.044,90.378] | 90.253 [88.821,92.851] |
| Completed child CPU ms | 87.524 [87.047,88.326] | 87.960 [86.852,90.920] |
| Pre-CSV VmHWM KiB | 11,616 [11,372,11,744] | 11,464 [11,412,11,704] |

Private factor median falls about 3%, with four of five pairs improving, but
its observed ranges overlap. Complete-child wall and CPU medians rise about
0.22% and 0.50%, with overlapping ranges and only two of five pairs improving.
No private end-to-end gain or precise regression is resolved. Fill remains
192,971 and live owned storage remains 3,552,692 bytes in every attempt.
Child wall includes input, checks, export and teardown; independent exact-audit
wall is separate. CSV columns retain that audit and child-plus-audit wall
without presenting child-only time as certificate-inclusive time. Completed
child CPU does not include the independent postprocess audit.

## Profile and correctness basis

Before editing, source-line sampling of the current pruned factorization found
the separate reach-count loop at generated real lines 319–321 and complex lines
1029–1031. It accounts for 10.84% of `factor_raw` local sampled period for real
RHS8, 8.54% for real refactor and 14.61% for complex refactor. Disassembly shows
its separate reach/pivot-map/work traversal and back edge. These are diagnostic
local-symbol shares, not percentages of workflow time or a predicted speedup.

Those nonexclusive diagnostics use source `e4295cf`, five replays per case,
`perf record -F 1999 -e cycles:u -g --call-graph dwarf,8192`, and identical
arithmetic flags plus `-g` for line attribution. Exact recorded sample counts are
709, 840 and 413, respectively. Raw traces, annotations, the retained baseline debug ELF
and replay/source hashes are retained separately from ordinary timing receipts.

Generation at the retained FortSym pin was repeated with identical output;
scalar generated leaves are unchanged. Eight focused bounds-checked native LU,
pruning, leaf, general API, refinement and AMD/cap gates pass. After composing
accepted LOW changes, LU/pruning/refinement pass again with bounds checking.
The added finite-input multiplier-overflow, failed-factor readiness,
same-pattern recovery and real/complex nonfinite-input oracle also passes
against the unchanged baseline. Existing exact cancellation, dynamic pivots,
N/T/H blocks, implicit derivative and cap oracles remain active. The graph
schedule and its independent closure oracle are unchanged.

The recipe records original runner and summary-script hashes. To reproduce
summaries from retained receipts without running a solver, execute the archived
`summary_script.path` from the repository root. Ordinary reruns require a
controller-approved quiet window and the recorded frozen ELF/runtime closure.
