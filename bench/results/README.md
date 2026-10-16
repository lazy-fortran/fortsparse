# Native SPD development measurements

The [final 2026-10-16 evidence](2026-10-16-final-evidence.md) supersedes the
intermediate checkpoints below for current native CPU comparisons. It retains
all 115 case-group statuses and distinguishes nondominance from the stronger
same-comparator match/beat gate.

The [GPU device-solve checkpoint](2026-10-16-gpu-device-solve.md) retains
320 interleaved trials against cuDSS, standalone CHOLMOD, native CPU and the
prior GPU-factor/CPU-solve binary. It is shared-GPU diagnostic evidence with
319 passes and one retained timeout, rather than a performance admission.

The [general-LU target-panel experiments](2026-10-16-lu-panels/README.md)
retain all 60 internal A/B observations and portable plot recipes. All
correctness gates pass, but the captured case remains slower than scalar
execution; no default or external-frontier admission follows.

## 2026-10-16: simplicial reference checkpoint

This is an intermediate reference, not the competitive production MVP.
`2026-10-16-simplicial-reference.csv` retains 675 individual successful trials:
9 fixed cases, 5 libraries, 5 fresh-process trials, one physical core, with
cold one-RHS, warm-runtime factorization, and derivative workloads. Raw extended
JSON and binary provenance are retained in `/var/tmp/fortsparse-mvp/` on the
measurement host. Reproduce with `bench/run_spd.py` and the commands below.

Hardware: Ryzen 9 5950X, GNU Fortran/C 16.2.1, Release `-O3`, physical-core
pinning; OpenBLAS 0.3.34 with one BLAS thread, SuperLU 7.0.1, SuiteSparse 7.14.1,
sequential/OpenMP MUMPS from `/mnt/storage/fs1-envs/mumps-local`, Intel oneMKL
2026.0 PARDISO in a separate executable. Each invocation completed below one
second. Fresh process means cold numerical runtime, not cold filesystem cache.

Largest base cases, medians of five cold trials:

| Case | Library | Factor + solve (ms) | Process CPU, including helper (ms) | Parent peak RSS (MiB) | Helper peak RSS (MiB) |
|------|---------|-------------------:|----------------------------------:|----------------------:|----------------------:|
| band, n=2048 | native reference | 0.390 | 0.737 | 6.63 | — |
| band, n=2048 | SuperLU | 1.326 | 1.687 | 9.82 | — |
| band, n=2048 | MUMPS | 2.173 | 2.501 | 12.28 | — |
| band, n=2048 | UMFPACK helper | 3.665 | 4.261 | 6.50 | 6.48 |
| band, n=2048 | Intel PARDISO | 4.088 | 6.412 | 20.13 | — |
| grid2, n=1024 | native reference | 1.349 | 1.464 | 6.50 | — |
| grid2, n=1024 | SuperLU | 1.308 | 1.427 | 9.28 | — |
| grid2, n=1024 | MUMPS | 1.453 | 1.574 | 11.59 | — |
| grid2, n=1024 | UMFPACK helper | 3.074 | 3.372 | 6.18 | 6.03 |
| grid2, n=1024 | Intel PARDISO | 3.262 | 5.264 | 19.99 | — |
| grid3, n=729 | native reference | 3.834 | 3.927 | 6.50 | — |
| grid3, n=729 | SuperLU | 1.645 | 1.772 | 9.40 | — |
| grid3, n=729 | MUMPS | 1.562 | 1.674 | 11.65 | — |
| grid3, n=729 | UMFPACK helper | 3.744 | 4.057 | 6.21 | 6.11 |
| grid3, n=729 | Intel PARDISO | 3.418 | 5.483 | 20.18 | — |

CPU covers fixture construction, checking and destruction as well as the
numerical phases. Wall phase times exclude those other operations. Parent and
helper RSS include runtime/fixture/vector memory and are reported separately;
their sum would double-count shared pages. They are not factor-allocation
measurements. No memory Pareto admission is inferred from incomparable totals.
Five raw trials, rather than this rounded table, determine variability.

Warm-runtime factor + solve medians on those same cases, in milliseconds:

| Case | Native reference | SuperLU | MUMPS | UMFPACK helper | Intel PARDISO |
|------|-----------------:|--------:|------:|---------------:|--------------:|
| band | 0.346 | 1.154 | 1.768 | 3.590 | 2.338 |
| grid2 | 1.295 | 1.394 | 1.097 | 2.983 | 1.691 |
| grid3 | 3.715 | 1.516 | 1.181 | 3.419 | 1.630 |

The reference wins the banded time comparison but loses the 3D mesh comparison.
This motivates blocked updates and fill-reducing ordering. Its numerical
factorization remains serial; requested threads do not constitute scaling.

Correctness: registered native and generated-leaf tests pass, along with the
benchmark fixture test. Separate registered tests of CUDA leaves, retained
panel/update execution and Fortran bindings pass on an RTX 5060 Ti in 0.28,
0.30 and 0.22 seconds respectively. These are GPU kernel/interface results,
not evidence of a complete sparse GPU factorization.

Example reproduction, after configuring the optional benchmark targets:

```sh
python3 bench/run_spd.py --bin-dir /path/to/build/bin \
  --output /var/tmp/reference-small.json --threads 1 --workloads rhs1
python3 bench/run_spd.py --bin-dir /path/to/build/bin \
  --output /var/tmp/reference-warm-derivatives.json --threads 1 \
  --workloads warm,derivatives
```

See `../README.md` for compiler/library setup and the thread-scaling commands.
