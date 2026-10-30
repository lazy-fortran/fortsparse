# Default ordering for the symmetric, Hermitian and capped SPD factors

`FORTSPARSE_ORDER_DEFAULT` now selects AMD (the rule of the general LU and
uncapped SPD paths: at least 4096 columns and an entry more than 12 positions
from the diagonal) for the real symmetric indefinite LDL^T, complex HPD
Cholesky, complex Hermitian LDL^H and complex symmetric LDL^T factors, and for
real SPD under a host cap, LOW_MEMORY or refinement. The indefinite factors
restrict their pivots to the order with delayed pivots, so the order survives
strongly indefinite matrices: on a 4900-unknown shifted Laplacian the new
default needs 3.5x to 4.5x fewer factor entries than natural order at every
shift from 0 to 3.5, where AMD with the previous unrestricted Bunch--Kaufman
pivoting needed up to 1.85x more than natural order. Inertia matched the exact
count in every run, and backward errors stayed below 1e-14.

## Method

Unrestricted Bunch--Kaufman takes its 2x2 partner, or a 1x1 replacement, from
anywhere in the remaining matrix and swaps positions with it. That moves a
vertex early in the AMD order to a late position, and a late vertex (often a
separator) forward, which undoes the order. The existing factorization keeps
the Schur complement as explicit sparse rows and chooses pivots one at a time,
so the restriction used by MA57 and MUMPS fits it without a supernodal
analysis: for the next vertex `k` in the order,

1. accept a 1x1 pivot if `|a_kk| >= u max_i |a_ik|`;
2. otherwise try the 2x2 pivot with the first remaining neighbour `p` of `k` in
   the order (its elimination-tree parent), accepted if all multipliers are at
   most `1/u` (`u |P^-1| (max |a_ik|, max |a_ip|)^T <= 1`);
3. otherwise delay `k` to the position right after `p` (the multifrontal
   "pass to the parent front"), keeping the relative order of all others;
4. a vertex delayed again before any elimination takes the unrestricted
   Bunch--Kaufman pivot, moved into place by the same order-preserving shift.

Step 4 guarantees a pivot and termination. Static pivoting with perturbations
was not chosen because it changes the factored matrix and needs iterative
refinement for every solve, and because the inertia, which the code reports,
would no longer be that of `A`. Natural order keeps the old unrestricted
pivoting unchanged (identical factors below). Threshold `u = 0.1`, chosen from
this sweep on the same matrix (`threshold.csv`, one run each):

| u | factor nnz, s = 0.3 / 1 / 2.01 / 3.5 | max backward error |
|---|---|---|
| 0.01 | 83,692 / 83,704 / 83,692 / 85,606 | 3.4e-13 |
| 0.1 | 84,008 / 84,758 / 96,249 / 87,975 | 9.6e-15 |
| 0.64 (Bunch--Kaufman alpha) | 108,574 / 184,127 / 224,925 / 346,903 | 3.5e-15 |

## Results

Native backend, CMake Release build, one thread, five fresh processes per
configuration unless noted; time and memory are medians [range]. "base" is
04db8da (main), "new" this branch; both linked against SuiteSparse AMD 3.3.4,
so the AMD orders are identical. Raw rows: `trials.csv`.

Real LDL^T (`SYM_INDEFINITE`), 70 x 70 five-point Laplacian with diagonal
`4 - s` (n = 4900). Natural order gives the same factors on both builds (one
confirming run on new).

| s | negative eigenvalues | factor nnz natural / base AMD / new default | factor s natural / base AMD / new default | peak RSS MiB natural / base AMD / new | backward error natural / new |
|---|---|---|---|---|---|
| 0 | 0 | 343,069 / 83,692 / 83,692 | 0.24 / 0.069 / 0.066 | 14.2 / 12.0 / 12.0 | 6.1e-16 / 5.2e-16 |
| 0.3 | 110 | 343,120 / 107,337 / 84,008 | 0.23 / 0.080 / 0.049 | 14.4 / 12.6 / 12.1 | 3.7e-15 / 6.4e-15 |
| 1 | 408 | 343,318 / 219,423 / 84,758 | 0.22 / 0.31 / 0.066 | 14.6 / 14.9 / 12.1 | 4.1e-15 / 6.4e-15 |
| 2.01 | 901 | 338,829 / 460,170 / 96,249 | 0.28 / 0.83 / 0.067 | 14.9 / 19.4 / 12.1 | 1.4e-14 / 5.2e-15 |
| 3.5 | 1902 | 397,017 / 735,210 / 87,975 | 0.29 / 1.39 / 0.053 | 15.6 / 22.4 / 12.1 | 1.9e-15 / 9.6e-15 |

Every run reported the exact inertia (negative count as listed, no zero
pivots). Time ranges per configuration are in `trials.csv`.

Real saddle point `[[L, B^T], [B, 0]]`, L the m x m Laplacian, B coupling
vertex pairs (inertia m^2 positive, m^2/2 negative, reported exactly by all):

| m (n) | configuration | factor nnz | factor s | peak RSS MiB |
|---|---|---|---|---|
| 40 (2400) | natural (= new default, n < 4096) | 1,025,239 | 5.7 [5.3-6.2] base, 6.8 [6.5-8.9] new | 35.1 |
| 40 (2400) | base AMD | 155,310 | 0.24 | 11.6 |
| 70 (7350) | natural, one run | 9,350,494 | 198 | 223 |
| 70 (7350) | base AMD | 968,867 | 2.7 [2.5-3.0] | 25.2 |
| 70 (7350) | new default (AMD) | 77,452 | 0.053 [0.040-0.054] | 12.8 |

Complex factors, base default (natural; these classes had no ordering) versus
new default (AMD). HPD: unitary gauge `D L D^H`, `D = diag(exp(0.3 i v))`;
LDL^H: `D (L - s I) D^H`, same spectrum as the real case; complex symmetric:
`(1 + 0.2i)(L - s I)`.

| class | m (n) | s | factor nnz base / new | factor s base / new | peak RSS MiB base / new | backward error base / new |
|---|---|---|---|---|---|---|
| HPD | 150 (22,500) | 0 | 3,375,149 / 540,630 | 1.43 / 0.16 | 77.8 / 24.5 | 1.1e-15 / 1.6e-15 |
| LDL^H | 70 (4900) | 0.3 | 343,143 / 84,022 | 1.07 / 0.18 | 17.3 / 12.8 | 4.9e-15 / 4.7e-15 |
| LDL^H | 70 (4900) | 3.5 | 418,608 / 99,015 | 1.13 / 0.29 | 19.1 / 13.4 | 3.5e-15 / 1.4e-14 |
| complex LDL^T | 70 (4900) | 3.5 | 398,413 / 87,975 | 1.09 / 0.19 | 18.7 / 12.8 | 3.1e-15 / 2.5e-14 |

LDL^H inertia was exact on both builds.

Real SPD, 200 x 200 Laplacian (n = 40,000), DEFAULT ordering. `plan` is the
`sparse_analyze` peak; `ledger` is the factorization's reported host peak.

| mode | factor nnz base / new | factor s base / new | peak RSS MiB base / new | plan = ledger peak, MB base / new |
|---|---|---|---|---|
| LOW_MEMORY | 8,000,199 / 1,081,911 | 0.56 / 0.12 | 106.0 / 32.2 | 98.2 / 18.9 |
| LOW_MEMORY + refinement | 8,000,199 / 1,081,911 | 0.52 / 0.072 | 109.1 / 35.4 | 101.4 / 22.1 |
| 8 GB cap (AUTO, blocked) | 8,618,155 / 1,102,860 | 0.34 / 0.085 | 87.2 / 33.9 | ledger 110.5 / 27.7 (no plan) |

The plan and the ledger peak agree to the byte in both LOW modes on both
builds; the allocator-level check of that equality is in
`test-cmake/native/test_spd_amd_bounded.f90`.

## Definitions and conditions

- factor nnz: `sparse_factor_stats%factor_nnz` (entries of L and D).
- factor s: wall time around `sparse_factor` (ordering, analysis, numeric).
- peak RSS: the process's `VmHWM` via `bench/process_metrics.c`, one process
  per row, no helpers. Process CPU (`process_cpu_s`) is in `trials.csv`.
- backward error: `max|Ax - b| / (||A||_inf max|x| + max|b|)` with the
  original `A`.
- AMD Ryzen 9 5950X shared with other users, 1-minute load 16-38 during the
  runs (column `load1`); runs pinned to core 5, `OMP_NUM_THREADS =
  OPENBLAS_NUM_THREADS = 1`, active threads 1. Times are indicative, not
  admission-grade; factor counts and memory do not depend on load.
- gfortran 16.2.1, `-O3` (CMake Release) for library and driver, OpenBLAS
  0.3.34, SuiteSparse 7.14.1 (AMD 3.3.4), METIS 5.2.1, Linux 7.2.11. Each
  run capped at 600 s and 16 GiB without swap; none hit a cap.
- The KKT m = 40 "new default" row is natural order because n < 4096.
- One grid shape per class; no held-out or consumer matrices.

## Reproduce

Configure both revisions with `cmake -G Ninja -DCMAKE_BUILD_TYPE=Release
-DFORTSPARSE_BUILD_TESTING=ON -DFORTSPARSE_USE_SYSTEM_DEPS=ON`, build the
`fortsparse` target, then compile `drivers/bench_symmetric_ordering.f90` with
`bench/process_metrics.c` against each build's `include/` and
`lib/libfortsparse.a` (plus SuperLU, fortnum, OpenBLAS, METIS, AMD) and run
`drivers/run_bench.sh BASE_BINARY NEW_BINARY trials.csv`. The threshold rows
come from rebuilding with `delay_threshold` set to 0.01 and to `alpha` in
`src/native/fortsparse_native_ldlt.f90`.
