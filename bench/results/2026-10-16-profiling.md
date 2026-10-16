# Profile-driven development checkpoint

This is intermediate evidence. The large cases have not met the MVP Pareto
criterion. All timings below use matched original fixtures and independent
backward-error checks. Ordinary invocations remain below one second.

## CPU profiles and the first optimization

Linux `perf` user-space cycle sampling works on this host. Bounded profiling
replays separate complete factorization and reused-pattern refactorization:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 0 \
  perf record -F 1999 -e cycles:u --call-graph dwarf -o /var/tmp/factor.data \
  /path/to/bin/bench_spd native grid2 7 1 profile-factor training
perf report --stdio --no-children --call-graph none -i /var/tmp/factor.data
```

The complete 16,384-unknown 2D profile collected 700 samples without lost
samples. METIS routines dominated setup; supernode construction used 7.1% of
cycles. A separate 13,824-unknown 3D refactor profile collected 857 samples:
27.4% of cycles were in the scalar generated `spd_update` function, with another
2.7% in division. Object-code inspection confirmed PLT calls inside generated
loops. PIC semantic interposition prevented the compiler from inlining those
mathematical leaves.

The GNU build now applies `-fno-semantic-interposition` to the generated module
only. This preserves PIC and mathematical semantics, permits inlining/SIMD, and
requires no hand edits to emitted mathematics. Independent generated value and
derivative checks, native solves/refactors, and supernode tests pass.

Five-trial warm-runtime factor-plus-solve medians:

| Case | Threads | Before (ms) | After (ms) | Before/after process CPU (ms) | Before/after parent peak RSS (MiB) |
|------|--------:|------------:|-----------:|-----------------------------:|---------------------------------:|
| 2D, n=16384 | 1 | 66.0 | 57.0 | 137.9 / 119.7 | 27.7 / 27.7 |
| 3D, n=13824 | 4 | 107.7 | 101.2 | 495.1 / 443.4 | 57.8 / 57.3 |

Process CPU includes both the first factorization and the measured warm
factorization, plus fixture/check/free work. It is not the CPU time of just
the warm phase. The smaller n=1024 2D and n=729 3D one-thread cases reach
0.59 ms and 0.94 ms respectively after this change.

The follow-up 3D profile collected 865 samples: scalar leaf calls disappeared;
mapped scatter now costs 22.6%, the dense kernel 16.8%, dense packing 6.5%, and
update orchestration 4.5%. This identifies compact amalgamation/scatter work
and ordering setup as the next targets. Sixteen threads currently add overhead;
requested thread counts alone are not successful scaling.

Raw traces and profiles live under `/var/tmp/fortsparse-mvp/perf-*` and are not
committed. The compact trial evidence accompanying this record retains CPU,
wall, memory and accuracy measurements.

## CUDA integration

The complete native sparse GPU path now retains all panels and updates on the
device, uploads assembled panel values once and downloads factors once for
CPU solves. Numeric refactorization reuses the context. A retained-scratch
regression caught and fixed an overwrite-versus-accumulate error before push.

Four registered CUDA kernel/binding/sparse-engine tests pass in 0.23–0.37 s
per case, and the public native backend CUDA fixture passes in 0.39 s. Thirty
fresh-process GPU trials (2D n=256/1024, one RHS/refactor/derivatives) all pass,
with actual device 0 reported. Cold factor-plus-solve is approximately
182–193 ms and parent peak RSS 527–536 MiB, including CUDA/vendor initialization.
These cases are slower and use more host memory than CPU execution. They are
correctness/crossover evidence, never native performance wins by fallback.

## Baseline repair and frozen scaling cases

The larger MUMPS baseline initially aborted because its bundled METIS header
and the system METIS runtime assigned different indices to the numbering
option, despite both reporting 5.2.1. A captured valid graph reproduced the
abort; changing only the option layout made the direct METIS call succeed.
The benchmark now explicitly selects MUMPS' matching METIS header and static
library. Static GKlib dependencies follow the selected static METIS target.

The registered larger MUMPS ordering regression passes in 0.05–0.06 s. All
160 external-only calibration invocations now pass, admitting 128x128 2D and
24x24x24 3D scaling cases. `../corpus.json` records that fixed selection.
The earlier initial table is historical evidence with its stated provider
configuration; use the repaired configuration for further admission runs.

## Compact factors and ordering checkpoint

Tree-aware amalgamation caps added structural zeros and releases temporary
simplicial analysis before numeric storage is touched. Linear-time symmetric
CSC permutation replaces triplet construction and sorting. Optional BSD AMD
is now the default for larger non-banded native problems; explicit METIS ND
remains available. Existing small and banded matrices retain natural ordering.
No external factorization is substituted for native execution.

Five-trial warm factor-plus-solve medians with GNU 16.2.1 `-O3 -march=znver3`:

| Case | Threads | Time (ms) | Parent peak RSS (MiB) |
|------|--------:|----------:|---------------------:|
| 2D, n=16384 | 1 | 17.16 | 20.90 |
| 2D, n=16384 | 4 | 38.90 | 27.25 |
| 3D, n=13824 | 1 | 79.79 | 38.92 |
| 3D, n=13824 | 4 | 63.54 | 49.01 |

These are intermediate results, not admission. Hardware ISA tuning alone
provided no substantial improvement over portable `-O3`. A separate METIS
refinement-count experiment saved setup time but increased fill; it was not
retained. AMD improved complete 3D setup more than that experiment. A bounded
16-thread 2D `perf` profile collected 3566 samples without losses: over 90% of
cycles were in libgomp task coordination. This motivates destination-panel
task coarsening before further thread comparisons. Portable CPU correctness,
explicit AMD/ND permutation/refactor checks and absent-ordering status tests
remain required. Numeric factor payload, retained auxiliary bytes and actual
process peaks are different measurements.

## Task granularity and direct dense updates

Destination-panel tasks reduced the 16-thread 2D warm time from 143.5 to
49.4 ms, but libgomp still consumed over 85% of sampled cycles. Batching by
estimated work reduced that case to 15.7 ms with five actual workers for a
16-thread request. Contiguous batches introduced extra dependency chains in
3D; complete elimination subtrees are the next schedule to measure.

A diagnostic link wrapper counted 10,563,508 mapped-scatter values in the
13,824-unknown 3D case: 27% belonged to entirely contiguous updates. Those
intersections now lower directly to lower-triangle BLAS SYRK and rectangular
GEMM on the owning destination panel. Irregular intersections retain generated
mapped scatter. No additional factor copy or scratch is introduced, and unused
upper triangles remain untouched. Independent native value/refactor/block
oracles and all five focused CPU tests pass.

Five fresh-process trials before/after, 3D one thread: warm factor-plus-solve
68.68 / 64.82 ms, changed-value refactor-plus-solve 43.68 / 39.66 ms, with
approximately 39.6 / 32.1 MiB process peaks in those respective workflows.
Warm memory includes two complete factorization lifetimes; refactor retains
one analysis. All 80 candidate trials pass their independent residual checks.
The complete external comparison has 5500/5500 passing trials; the frozen
HB/bcsstk09 structural stiffness case adds 500/500 external passes.

Two other bounded experiments were rejected: BLAS triangular-solve lowering
provided no end-to-end gain, and segmented generated scatter introduced 2D
regressions despite some 3D improvement. Their raw diagnostics remain in
`/var/tmp/fortsparse-mvp/`; their implementation changes were removed.
