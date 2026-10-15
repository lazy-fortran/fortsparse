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
