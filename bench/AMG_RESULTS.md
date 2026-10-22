# AMG convergence and remote scaling

[V] Synthetic Dirichlet grids only. Measurements use d8d0417f51e4714b9424887c39abc5b2ed75c1b5.
Scaling job 3355222 ran on sCluster node18 (two AMD EPYC 9455 sockets,
48 physical cores each, SMT disabled) for 16 min 56 s and completed
with exit 0. It used a 16-core allocation with fixed OpenMP places (`cores`),
close binding, dynamic teams disabled, and one OpenBLAS thread. GNU Fortran
12.2.0 used Release `-O3 -fopenmp` and the project's GNU-12
`-finline-matmul-limit=0` portability option. FortNum is pinned at
8c1d82bb635af6381a6873db036d56f2a0d39fd8 and FortGen runtime at
68a38f1a2fc71e6a35c79653392687875dad7e7b.

Receipt job 3355322 on node18 verified the dynamic libraries. OpenBLAS SHA256:
3c5ecc47b35745ac84ecdd91fb4e422d00996678f2c11e7c22d6ca69c0ddfca9.
The complete library identities are in `receipt-3355322.log` beside the raw data.
An OpenMP team is started before the setup timer, so setup timings exclude
initial team creation. Full process receipts retain that cost.

The benchmark binary SHA256 is
76d4972e5b633f3cf68542716ae1f7c8ef3f911b6b1080a9586ab722a019de78.
The source snapshot is retained outside Git at
`scluster:~/work/amg1-20261022/source-d8d0417.tar.gz`.
Later changes add error handling, test discovery/coverage, CI and documentation;
the measured default numerical path is unchanged.

## Accuracy and interpretation

All iterative comparisons use the same true relative residual target 1e-12
and maximum 500 iterations. The manufactured solution is
`sin(0.271*i) + cos(0.173*i)`. The unit tests separately compare solutions
with native LDLT to maximum absolute difference 1e-10 for isotropic and
anisotropic 2D/3D cases. Million-row 2D AMG forward error is about 6.4e-10,
whereas the residual is below 1e-12. The benchmark does not claim 1e-10
forward agreement at that size, where LDLT is censored.

ILUTP-PCG is the requested empirical control, with fill 12/12, drop 1e-3 and
pivot tolerance zero. ILUTP is generally nonsymmetric. LDLT uses the same
frozen source pin, before the separate default-row-merge change in PR #28.
The thread column is the requested budget, checked against an OpenMP startup
probe. This frozen default LDLT factorization is serial despite that budget.
No current-upstream LDLT speed or sparse-solver Pareto admission is claimed.

Five fresh processes are retained for each completed configuration. A first
15-second timeout censors the remaining repetitions of that configuration.
There are 500 timing trials across 120 configurations, including 25 censored
LDLT configurations. Shared-host variation remains in the recorded ranges.
Setup excludes CSR generation; solve includes Krylov workspace creation and
the independently recomputed residual. Peak RSS covers the complete process,
including generation, caller arrays and runtime, and is not solver payload.
The Linux in-process VmHWM observation avoids inherited launcher high-water
floors. Cumulative user+system CPU is sampled before final output and cleanup;
full-lifetime GNU time receipts remain in the raw logs.

## h-refinement

One thread, fixed threshold 0.02, two Chebyshev steps, V-cycle, coarse size 64.
P1-like uses the 5/7-point graph. Q1 uses exact tensor stiffness/mass products.

| Dimension | Side | Rows | P1-like iterations | Q1 iterations |
|---|---:|---:|---:|---:|
| 2 | 16 | 256 | 16 | 15 |
| 2 | 32 | 1024 | 16 | 16 |
| 2 | 64 | 4096 | 16 | 16 |
| 2 | 128 | 16384 | 17 | 16 |
| 3 | 4 | 64 | 1 | 1 |
| 3 | 8 | 512 | 16 | 14 |
| 3 | 16 | 4096 | 16 | 16 |
| 3 | 32 | 32768 | 17 | 16 |

The 4^3 grids use the direct coarse solve. The timing sweep extends P1-like
grids to one million rows: 23 iterations in 2D and 22 in 3D.

## Grid-aligned anisotropy

Coefficient tensor is diag(1, eps) or diag(1, eps, eps). All thresholds
are fixed across epsilon. Cells give iterations / setup milliseconds.
The 2D side is 64 and the 3D side is 16. These are single-process
observations. Low iteration counts can accompany expensive coarse LDLT
when strong-graph coarsening stalls; setup costs are retained.

| Dimension | eps | threshold 0 | threshold 0.02 | threshold 0.25 |
|---|---:|---:|---:|---:|
| 2 | 0.01 | 50 / 2.19 | 16 / 6.63 | 16 / 124 |
| 2 | 0.001 | 104 / 2.1 | 16 / 9.41 | 16 / 114 |
| 2 | 0.0001 | 133 / 2.18 | 16 / 14.6 | 16 / 28.4 |
| 2 | 1e-05 | 138 / 2.16 | 16 / 14.2 | 16 / 28.7 |
| 2 | 1e-06 | 138 / 2.18 | 16 / 13.4 | 16 / 31.9 |
| 3 | 0.01 | 31 / 4.5 | 15 / 31 | 15 / 1.34e+03 |
| 3 | 0.001 | 39 / 4.12 | 14 / 61.5 | 14 / 1.32e+03 |
| 3 | 0.0001 | 41 / 4.31 | 14 / 331 | 14 / 1.34e+03 |
| 3 | 1e-05 | 41 / 4.27 | 14 / 328 | 14 / 1.32e+03 |
| 3 | 1e-06 | 41 / 4.3 | 13 / 326 | 13 / 1.47e+03 |

## Timings

Median setup / solve milliseconds. Raw per-trial commands and full
process receipts are retained at `scluster:~/work/amg1-20261022/`.
[amg_measurements.csv](amg_measurements.csv) contains every median/range,
CPU sample, RSS, forward error, residual, iteration count and command.
[amg_convergence.csv](amg_convergence.csv) contains the 46 convergence rows.

| Dimension | Rows | Threads | AMG setup / solve ms | ILUTP setup / solve ms | LDLT setup / solve ms |
|---|---:|---:|---:|---:|---:|
| 2 | 1024 | 1 | 0.743 / 0.536 | 1.05 / 0.443 | 57.2 / 0.151 |
| 2 | 1024 | 2 | 0.607 / 0.451 | 0.979 / 0.413 | 52.7 / 0.152 |
| 2 | 1024 | 4 | 0.836 / 1.71 | 1.08 / 0.678 | 52.7 / 0.163 |
| 2 | 1024 | 8 | 0.903 / 1.62 | 1.1 / 0.766 | 53.1 / 0.178 |
| 2 | 1024 | 16 | 1.42 / 3.08 | 1 / 0.912 | 57.6 / 0.231 |
| 2 | 16384 | 1 | 8.13 / 7.49 | 14 / 16.7 | 1.36e+04 / 11.9 |
| 2 | 16384 | 2 | 5.87 / 4.32 | 13.9 / 15.5 | 1.36e+04 / 15.1 |
| 2 | 16384 | 4 | 6.33 / 5.77 | 14.8 / 16.9 | 1.35e+04 / 15.6 |
| 2 | 16384 | 8 | 4.83 / 4.39 | 13.8 / 15.7 | 1.36e+04 / 15.5 |
| 2 | 16384 | 16 | 5.44 / 6.24 | 15 / 18.6 | 1.36e+04 / 15.7 |
| 2 | 262144 | 1 | 106 / 257 | 228 / 888 | censored >15 s |
| 2 | 262144 | 2 | 71.8 / 141 | 235 / 837 | censored >15 s |
| 2 | 262144 | 4 | 58.2 / 59 | 233 / 812 | censored >15 s |
| 2 | 262144 | 8 | 50.6 / 49.9 | 228 / 805 | censored >15 s |
| 2 | 262144 | 16 | 50.4 / 31.3 | 226 / 854 | censored >15 s |
| 2 | 1000000 | 1 | 429 / 1.19e+03 | 861 / 1.98e+03 | censored >15 s |
| 2 | 1000000 | 2 | 278 / 670 | 865 / 1.85e+03 | censored >15 s |
| 2 | 1000000 | 4 | 212 / 355 | 858 / 1.74e+03 | censored >15 s |
| 2 | 1000000 | 8 | 179 / 335 | 855 / 1.73e+03 | censored >15 s |
| 2 | 1000000 | 16 | 176 / 183 | 864 / 1.73e+03 | censored >15 s |
| 3 | 512 | 1 | 0.763 / 0.415 | 1.5 / 0.247 | 93.5 / 0.12 |
| 3 | 512 | 2 | 0.622 / 0.37 | 1.52 / 0.247 | 89.3 / 0.133 |
| 3 | 512 | 4 | 0.883 / 1.6 | 1.44 / 0.494 | 90.6 / 0.143 |
| 3 | 512 | 8 | 0.757 / 1.43 | 1.48 / 0.551 | 90.2 / 0.148 |
| 3 | 512 | 16 | 1.13 / 2.72 | 1.51 / 0.738 | 89.8 / 0.23 |
| 3 | 32768 | 1 | 25.2 / 21.6 | 61.6 / 37.9 | censored >15 s |
| 3 | 32768 | 2 | 15.9 / 11.5 | 59.5 / 36 | censored >15 s |
| 3 | 32768 | 4 | 13.4 / 9.9 | 63.9 / 39.1 | censored >15 s |
| 3 | 32768 | 8 | 11.5 / 8.44 | 59.6 / 35.9 | censored >15 s |
| 3 | 32768 | 16 | 13.7 / 10.1 | 59.8 / 40.2 | censored >15 s |
| 3 | 262144 | 1 | 253 / 355 | 471 / 448 | censored >15 s |
| 3 | 262144 | 2 | 167 / 184 | 473 / 415 | censored >15 s |
| 3 | 262144 | 4 | 114 / 86.6 | 473 / 407 | censored >15 s |
| 3 | 262144 | 8 | 85.3 / 73.6 | 476 / 397 | censored >15 s |
| 3 | 262144 | 16 | 81.1 / 42.3 | 474 / 404 | censored >15 s |
| 3 | 1000000 | 1 | 1.03e+03 / 1.59e+03 | 1.78e+03 / 1.89e+03 | censored >15 s |
| 3 | 1000000 | 2 | 681 / 875 | 1.77e+03 / 1.73e+03 | censored >15 s |
| 3 | 1000000 | 4 | 464 / 472 | 1.78e+03 / 1.63e+03 | censored >15 s |
| 3 | 1000000 | 8 | 347 / 420 | 1.79e+03 / 1.62e+03 | censored >15 s |
| 3 | 1000000 | 16 | 289 / 234 | 1.78e+03 / 1.6e+03 | censored >15 s |

## Million-row AMG scaling

Speedup and efficiency use the one-thread median for the same dimension.

| Dimension | Threads | Iterations | Setup speedup | Solve speedup | Solve efficiency | Peak RSS MiB |
|---|---:|---:|---:|---:|---:|---:|
| 2 | 1 | 23 | 1.00 | 1.00 | 100.0% | 333.1 |
| 2 | 2 | 23 | 1.54 | 1.78 | 88.9% | 333.2 |
| 2 | 4 | 23 | 2.02 | 3.36 | 83.9% | 330.0 |
| 2 | 8 | 23 | 2.39 | 3.55 | 44.4% | 337.7 |
| 2 | 16 | 23 | 2.44 | 6.49 | 40.6% | 337.7 |
| 3 | 1 | 22 | 1.00 | 1.00 | 100.0% | 448.6 |
| 3 | 2 | 22 | 1.51 | 1.82 | 91.0% | 446.8 |
| 3 | 4 | 22 | 2.22 | 3.37 | 84.3% | 446.8 |
| 3 | 8 | 22 | 2.97 | 3.80 | 47.5% | 451.4 |
| 3 | 16 | 22 | 3.56 | 6.81 | 42.6% | 461.9 |

## Reproduction

Use the source pin above or its retained archive. Build and run on a
compute allocation, not a login node or the local workstation:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DFORTSPARSE_ENABLE_SUPERLU=OFF -DFORTSPARSE_ENABLE_UMFPACK_IPC=OFF \
  -DFORTSPARSE_ENABLE_AMD=OFF -DFORTSPARSE_ENABLE_METIS=OFF \
  -DFORTSPARSE_BUILD_BENCHMARKS=ON
cmake --build build --target bench_amg -j4
export OPENBLAS_NUM_THREADS=1 OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close OMP_PLACES=cores OMP_NUM_THREADS=16
export AMG_REVISION=d8d0417f51e4714b9424887c39abc5b2ed75c1b5
timeout 15 build/bin/bench_amg 1000 2 1 1 .02 amg 500
```

The exact complete sweep is retained in `scaling.sh` with the raw logs
and `evidence.tar.gz`. The analysis command is
`python ~/work/amg1-20261022/summarize.py scaling-3355222.log`.
Pilot 3355200/node8 and final validation 3355287/node4 are separate from
the timing series. Six CTests passed with `-fcheck=all -fbacktrace`,
and fpm passed both new test targets with `-fopenmp -fcheck=all`.
See [API contracts](../docs/amg.md) for the test names and limitations.
