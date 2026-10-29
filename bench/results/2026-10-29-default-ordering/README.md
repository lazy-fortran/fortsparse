# Default native ordering: natural versus AMD

`FORTSPARSE_ORDER_DEFAULT` for native general LU changed from natural order to
AMD (approximate minimum degree) for matrices with at least 4096 columns and an
entry more than 12 positions from the diagonal, the rule real SPD already used.
On the two row-numbered 2D meshes below (one size each), factor entries fall
7.0x and 10.3x and process peak memory 4.4x and 6.5x, and both residuals stay
at roundoff level. Factor time fell 6.7x and 15.9x on a shared host, so the
times show the size of the effect but are not admission-grade measurements.

## Results

Native general LU (`FORTSPARSE_SYM_GENERAL`), one thread
(`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=1`), pinned to core 2 (`taskset -c 2`),
one fresh process per trial. Default (AMD): five trials, median [range].
Natural order: one trial each, because a single run takes 18 s and 98 s.
Raw rows: `trials.csv`.

| Matrix | n | nnz(A) | Ordering | factor nnz | factor wall s | process CPU s | peak RSS MiB | rel. residual |
|---|---|---|---|---|---|---|---|---|
| Laplacian, real | 100,489 | 501,177 | natural | 63,710,658 | 17.6 | 16.2 | 755 | 1.4e-15 |
| Laplacian, real | 100,489 | 501,177 | default (AMD) | 6,207,474 | 1.11 [1.08-1.20] | 1.10 [1.06-1.18] | 116 [116-116] | 1.6e-15 |
| Elasticity type, complex symmetric | 200,978 | 2,803,556 | natural | 255,440,502 | 98.1 | 93.2 | 4,968 | 3.2e-16 |
| Elasticity type, complex symmetric | 200,978 | 2,803,556 | default (AMD) | 36,625,446 | 14.7 [14.0-14.9] | 14.1 [13.6-14.1] | 1,136 [1,136-1,136] | 7.7e-16 |

- Matrices, both on 317x317 grid vertices numbered row by row and generated in
  `drivers/bench_default_ordering.f90`:
  - Laplacian: five-point stencil, diagonal 4, off-diagonals -1.
  - Elasticity type: two unknowns per vertex on a structured triangulation
    with linear (P1) elements. Each triangle adds the 6x6 block
    `(1+0.2i) kron(L, C)`, with a positive semidefinite 3x3 `L` and SPD 2x2 `C`;
    `(0.05+0.3i) I` is added. The matrix equals its transpose (not its
    conjugate transpose) and is nonsingular.
- factor nnz: entries of L and U reported by `sparse_factor_stats`.
- factor wall: in-process timer around `sparse_factor`, which includes
  ordering, symbolic work and numerical factorization. The AMD ordering step
  alone (`stats%analysis_s`) took 0.027 s and 0.10 s (medians).
- process CPU: user plus system time of the whole process, sampled before
  output, so it includes matrix generation.
- peak RSS: the process's own resident high-water mark (`VmHWM` from
  `/proc/self/status`, which resets at exec) via `bench/process_metrics.c`;
  one process per row, no helper processes.
- rel. residual: `max|Ax-b| / (c max|x| + max|b|)` with `c` = 8 for the
  Laplacian and 20 for the elasticity matrix, fixed scales in the driver.
- `stats%ordering_used` was 1 (AMD) for every default row and 0 (natural) for
  every natural row.
- Shared host: AMD Ryzen 9 5950X (16 cores), 1-minute load average 44-50 from
  other users. An earlier campaign at load 75 gave the same factor nnz and
  1.4-1.9x longer times. Fill and memory do not depend on load; times do, and
  natural order has one trial per matrix.
- Toolchain: FortSparse ff261c0 built by fpm 0.13.0 with `--profile release`
  (gfortran 16.2.1, `-O3 -funroll-loops -fopenmp`), OpenBLAS 0.3.34, Linux
  7.2.11. The fpm build contains the project-owned AMD, not SuiteSparse. Every
  run was capped at 600 s and 16 GiB without swap; none reached a cap.

## AMD quality and nested dissection

Raw outputs of the three diagnostic drivers are in `diagnostics.txt`.
`drivers/amd_quality.c` compares Cholesky fill of the symmetric pattern under
the project-owned AMD and SuiteSparse AMD 3.3.4. The ratio is 0.99-1.01 on
randomly renumbered 2D grids up to n=160,000, 1.03 on a 40^3 3D grid and
0.99-1.00 on random graphs. With the minimum workspace (`MEMORY_LOW`), the
core compacted its graph 191 times over the run (counted in an instrumented
copy) and returned the same orders.

`drivers/amd_vs_nd.c` on the elasticity pattern above (symbolic Cholesky of
the symmetric pattern, n=200,978):

| Ordering | nnz(L) | flops (sum of squared column counts) |
|---|---|---|
| natural | 127.3 M | 8.1e10 |
| project-owned AMD | 19.3 M | 8.1e9 |
| SuiteSparse AMD 3.3.4 | 18.5 M | 7.0e9 |
| METIS 5.2.1 nested dissection | 16.1 M | 4.6e9 |

On this one 2D mesh, nested dissection (ND) needs 17% less fill and 43% fewer
flops than the project-owned AMD (13% and 34% fewer than SuiteSparse AMD),
against 6.6x less fill and 10x fewer flops for AMD over natural order. ND
usually gains more on 3D meshes. It stays available through optional METIS;
a native ND was not attempted.

## Why the scalar real LDLT keeps natural order by default

The scalar (non-supernodal) real symmetric indefinite factorization
(`FORTSPARSE_SYM_INDEFINITE`) uses Bunch-Kaufman pivoting without restricting
the pivot candidates. `drivers/ldlt_probe.f90` factors the 70x70 five-point
Laplacian with diagonal `4-s`, which is indefinite for s above about 0.004:

| s | factor nnz, natural | factor nnz, AMD | factor s, natural | factor s, AMD |
|---|---|---|---|---|
| 0 | 343,069 | 80,696 | 0.25 | 0.05 |
| 0.05 | 343,076 | 80,696 | 0.29 | 0.06 |
| 0.3 | 343,120 | 100,671 | 0.27 | 0.09 |
| 1.0 | 343,318 | 241,581 | 0.30 | 0.48 |
| 2.0137 | 337,880 | 787,870 | 0.23 | 5.37 |
| 3.5 | 397,017 | 2,014,465 | 0.31 | 25.9 |

AMD wins when the matrix is close to definite and loses badly when it is
strongly indefinite (from s=1 in time, from s=2 in fill). That a default order
cannot be chosen safely without knowing the pivots is why this factorization
stays natural by default; the mechanism (pivot interchanges with distant rows
undoing the order) is a hypothesis, as pivot counts were not recorded. This is
one grid size.

## Reproduce

In a scratch directory next to a FortSparse checkout at ff261c0 named
`fortsparse`: copy `drivers/fpm.toml` and `drivers/run_bench.sh`, put
`drivers/*.f90` and `fortsparse/bench/process_metrics.c` in `app/`, run
`fpm build --profile release`, then `./run_bench.sh` (writes `results.csv`)
and `fpm run --profile release ldlt_probe`. From the FortSparse checkout the C
drivers build with
`gcc -O2 -Isrc/native -I/usr/include/suitesparse <driver>.c
src/native/fsparse_native_amd_core.c src/native/fsparse_native_amd_bounded.c
-lamd -lmetis` and run as `./amd_quality` and `./amd_vs_nd 317`.
