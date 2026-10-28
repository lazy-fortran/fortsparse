# Backends

Every backend is reached through the same calls (`sparse_factor`,
`sparse_solve`, `sparse_solve_block`, `sparse_refactor`, `sparse_free`). Select
one by setting `solver%backend_id` before factoring:

| Constant | Backend |
| --- | --- |
| `FORTSPARSE_BACKEND_AUTO` (default) | Chosen per factorization; see [automatic selection](#automatic-selection) |
| `FORTSPARSE_BACKEND_SUPERLU` | SuperLU (BSD), in process |
| `FORTSPARSE_BACKEND_SUPERLU_ILU` | SuperLU incomplete LU, real only |
| `FORTSPARSE_BACKEND_MUMPS` | MUMPS (CeCILL-C), in process, optional |
| `FORTSPARSE_BACKEND_UMFPACK_IPC` | UMFPACK (GPL) in a separate helper process |
| `FORTSPARSE_BACKEND_NATIVE` | Project-owned native solvers |
| `FORTSPARSE_BACKEND_NATIVE_REFERENCE` | Serial simplicial native reference, for comparisons |

A backend that was not built returns `FORTSPARSE_BACKEND_UNAVAILABLE`.

## Native backend

The native backend is written in this repository and needs only BLAS/LAPACK;
AMD and METIS orderings are optional.

For real SPD matrices set `solver%options%symmetry = FORTSPARSE_SYM_SPD` and
supply full symmetric CSC with sorted, unique row indices and finite values
(`csc_from_triplet` canonicalizes duplicates and ordering). The solver uses a
simplicial path for low fill and blocked OpenMP tasks for larger fronts. AMD is
the default ordering for larger non-banded systems when available, with METIS
nested dissection as an alternative; `FORTSPARSE_ORDER_AMD` and
`FORTSPARSE_ORDER_ND` require the corresponding library. A nonpositive or
nonfinite pivot returns `FORTSPARSE_NOT_POSITIVE_DEFINITE`. The analysis is
kept for same-pattern refactorization. Statistics report factor storage and
the requested and actual number of threads; they are not a forward-error bound.

Other matrix classes (real and complex general LU, complex HPD Cholesky, real
symmetric, complex Hermitian and complex symmetric indefinite LDL^T) are
described in [native CPU matrix classes](native-cpu.md). Memory limits and the
low-memory mode are described in [memory](memory.md); fp32 factor storage with
fp64 refinement in [precision](precision.md#native-direct-factors-with-adaptive-storage).

### CUDA

Build with `FORTSPARSE_ENABLE_NATIVE_CUDA=ON` and set
`solver%options%device = 0` to factor real SPD matrices on CUDA device 0. Panels
and Schur updates stay on the device during factorization, and triangular
solves reuse the device factors. Host-array interfaces transfer right-hand
sides and results. The default device is `-1` (CPU); a GPU request in a CPU-only
build returns `FORTSPARSE_BACKEND_UNAVAILABLE`. Device statistics count the
library's own allocations and support an optional device memory cap; they do
not include driver-internal pools.

## Implicit derivatives

Real solves provide exact derivatives without differentiating the backend. For
\(A x = b\), `sparse_solve_jvp` reuses the primal factorization to solve
\(A \dot{x} = \dot{b} - \dot{A} x\). `sparse_solve_vjp` solves
\(A^T \bar{b} = \bar{x}\) with a factorization of `transpose(A)` and returns the
CSC-value cotangents \(\bar{A}_{ij} = -\bar{b}_i x_j\):

```fortran
call sparse_solve_jvp(solver, A_dot, x, b_dot, x_dot, status)

call csc_transpose(A, transpose_A, status)
call sparse_factor(transpose_solver, transpose_A, status)
call sparse_solve_vjp(transpose_solver, A, x, x_bar, &
                      b_bar, A_values_bar, status)
```

This works with every real direct backend.

## SuperLU

SuperLU factors the full supplied CSC matrix for every symmetry hint. The hints
`FORTSPARSE_SYM_SPD`, `FORTSPARSE_SYM_INDEFINITE` and `FORTSPARSE_SYM_HERMITIAN`
select SuperLU's symmetric mode: minimum-degree ordering on the structure of
`A + transpose(A)` and diagonal threshold pivoting (`DiagPivotThresh = 0.001`).
The default `FORTSPARSE_SYM_GENERAL` uses COLAMD and ordinary partial pivoting.
The hints are not checked; the caller must declare the structure correctly.

`FORTSPARSE_BACKEND_SUPERLU_ILU` selects an equilibrated real SuperLU
incomplete LU (drop tolerance `1e-4`, fill factor `10`) used as an approximate
inverse. It is never chosen by AUTO and gives no direct-solve residual guarantee.
Complex matrices, refinement and JVP/VJP requests return
`FORTSPARSE_BACKEND_UNAVAILABLE`.

## Automatic selection

`FORTSPARSE_BACKEND_AUTO` decides at each factorization from the number of
stored entries: MUMPS for `nnz >= 40000`, SuperLU below
(`sparse_auto_backend(nnz, is_complex)` returns the choice). MUMPS is chosen only
when it is built, for complex matrices only with ZMUMPS, and for an MPI build only
while the application has initialized MPI; AUTO never initializes MPI itself.
Without MUMPS, AUTO uses SuperLU if enabled and otherwise the native backend. A
same-pattern `sparse_refactor` resolves to the same backend, so the MUMPS
analysis is reused.

The crossover comes from `bench_fortsparse_backends` (built with
`FORTSPARSE_ENABLE_MUMPS=ON`) on finite-element matrices: 2D P1/P2/P3 SPD
stiffness, `hcurl` (complex edge plus nodal unknowns, indefinite), `blk3` (P2
with three coupled fields, nonsymmetric LU), `ideal` (the same, symmetric
indefinite LDL^T) and 3D P1 tetrahedra with one (`tet1`, SPD) and three (`tet3`,
LU) fields. Times are factor plus one solve plus destroy, i.e. the cost of
`sparse_solve_once`; the last column is one solve with a retained factor.
Measured on an AMD Ryzen 9 5950X with MUMPS 5.9.1 (OpenMP, METIS) and serial
SuperLU 7.0.1, reference BLAS unless noted, with other load on the machine.

| kind | n | nnz | SuperLU | MUMPS 1 thr | MUMPS 8 thr | SuperLU, OpenBLAS | MUMPS 8 thr, OpenBLAS | solve SuperLU / MUMPS 8 thr |
| --- | --: | --: | --: | --: | --: | --: | --: | --: |
| p1 | 49 | 289 | 0.0636 ms | 0.0829 ms | 0.0955 ms | 0.0884 ms | 0.0983 ms | 0.001 / 0.0138 ms |
| p1 | 441 | 2921 | 0.459 ms | 0.503 ms | 0.522 ms | 0.502 ms | 0.558 ms | 0.0126 / 0.0758 ms |
| p1 | 3721 | 25561 | 6.18 ms | 4.78 ms | 4.66 ms | 6.08 ms | 4.32 ms | 0.153 / 0.546 ms |
| p1 | 40401 | 281201 | 170 ms | 98.8 ms | 75.6 ms | 129 ms | 45.9 ms | 5.12 / 3.76 ms |
| p1 | 361201 | 2523601 | 4.62 s | 2.08 s | 1.49 s | 2.99 s | 389 ms | 138 / 37.4 ms |
| p2 | 81 | 801 | 0.0862 ms | 0.123 ms | 0.138 ms | 0.112 ms | 0.172 ms | 0.00217 / 0.0186 ms |
| p2 | 625 | 6817 | 0.751 ms | 0.787 ms | 0.843 ms | 0.816 ms | 0.78 ms | 0.0197 / 0.107 ms |
| p2 | 6561 | 74241 | 12.5 ms | 10.2 ms | 12.8 ms | 12.7 ms | 9.68 ms | 0.334 / 0.768 ms |
| p2 | 58081 | 664321 | 283 ms | 250 ms | 165 ms | 205 ms | 66.4 ms | 7.83 / 5.88 ms |
| p3 | 100 | 1468 | 0.103 ms | 0.139 ms | 0.137 ms | 0.129 ms | 0.151 ms | 0.00301 / 0.019 ms |
| p3 | 625 | 10033 | 0.718 ms | 0.786 ms | 0.918 ms | 0.75 ms | 0.887 ms | 0.0239 / 0.104 ms |
| p3 | 5776 | 96376 | 11.2 ms | 9.83 ms | 11.2 ms | 10.1 ms | 11.8 ms | 0.432 / 0.541 ms |
| p3 | 58081 | 981601 | 270 ms | 257 ms | 182 ms | 190 ms | 74.4 ms | 9.98 / 8.32 ms |
| hcurl | 81 | 801 | 0.152 ms | 0.213 ms | 0.238 ms | 0.15 ms | 0.183 ms | 0.00459 / 0.0179 ms |
| hcurl | 625 | 6817 | 1.48 ms | 1.48 ms | 1.73 ms | 1.7 ms | 1.44 ms | 0.0386 / 0.153 ms |
| hcurl | 6561 | 74241 | 52 ms | 28.9 ms | 25.2 ms | 37.5 ms | 18.3 ms | 2.19 / 1.1 ms |
| hcurl | 58081 | 664321 | 2.29 s | 654 ms | 464 ms | 1.39 s | 160 ms | 26.8 / 12.1 ms |
| blk3 | 147 | 4167 | 0.187 ms | 0.338 ms | 0.4 ms | 0.209 ms | 0.369 ms | 0.00486 / 0.0226 ms |
| blk3 | 867 | 27657 | 1.96 ms | 2.45 ms | 2.72 ms | 2.14 ms | 2.42 ms | 0.0469 / 0.133 ms |
| blk3 | 7803 | 262359 | 57.7 ms | 47.1 ms | 43.3 ms | 46.5 ms | 26 ms | 2.05 / 0.879 ms |
| blk3 | 59643 | 2038689 | 2.28 s | 962 ms | 571 ms | 1.43 s | 212 ms | 64.8 / 11.4 ms |
| ideal | 147 | 4167 | 0.212 ms | 0.347 ms | 0.355 ms | 0.205 ms | 0.419 ms | 0.00626 / 0.0268 ms |
| ideal | 867 | 27657 | 1.98 ms | 2.38 ms | 2.49 ms | 2.02 ms | 2.43 ms | 0.0463 / 0.144 ms |
| ideal | 7803 | 262359 | 34.3 ms | 42.3 ms | 34.7 ms | 26.8 ms | 19.2 ms | 0.96 / 0.897 ms |
| ideal | 59643 | 2038689 | 735 ms | 689 ms | 487 ms | 454 ms | 133 ms | 22.8 / 10.1 ms |
| tet1 | 125 | 1333 | 0.176 ms | 0.187 ms | 0.192 ms | 0.151 ms | 0.231 ms | 0.00573 / 0.0246 ms |
| tet1 | 1000 | 12718 | 4.41 ms | 3.63 ms | 3.8 ms | 4.02 ms | 2.58 ms | 0.0674 / 0.234 ms |
| tet1 | 3375 | 45403 | 47 ms | 37.2 ms | 43.7 ms | 28.8 ms | 24.5 ms | 0.556 / 0.901 ms |
| tet1 | 9261 | 128581 | 409 ms | 193 ms | 165 ms | 188 ms | 71.1 ms | 2.91 / 1.91 ms |
| tet1 | 29791 | 424171 | 5.45 s | 1.56 s | 1.2 s | 2.59 s | 277 ms | 21.7 / 10.5 ms |
| tet3 | 192 | 5598 | 0.484 ms | 0.509 ms | 0.669 ms | 0.458 ms | 0.581 ms | 0.0107 / 0.0336 ms |
| tet3 | 1536 | 56142 | 19.8 ms | 19.6 ms | 19.8 ms | 13.4 ms | 7.98 ms | 0.256 / 0.388 ms |
| tet3 | 5184 | 203454 | 253 ms | 196 ms | 144 ms | 140 ms | 41.1 ms | 3.01 / 1.72 ms |
| tet3 | 12288 | 499374 | 1.7 s | 985 ms | 755 ms | 898 ms | 123 ms | 9.22 / 6.71 ms |

SuperLU is faster below about 3e4 entries, where its solves are also 3 to 15
times faster. From about 5e4 entries MUMPS is as fast or faster, increasingly so
with size, threads and an optimized BLAS (with OpenBLAS, 2.5 to 9 times faster at
1e5 entries and above). Only serial SuperLU is built.

MUMPS's OpenMP tree parallelism (ICNTL(48)) costs about 1 ms per factorization
on small problems, so fortsparse disables it below n = 5000. The default
ordering analyses with AMD and switches to METIS nested dissection when the AMD
estimate exceeds 3000 elimination flops per entry (three-dimensional or strongly
coupled patterns). `FORTSPARSE_ORDER_AMD` and `FORTSPARSE_ORDER_ND` force either.

## MUMPS

`FORTSPARSE_BACKEND_MUMPS` targets large SPD finite-element systems. It runs a
multifrontal factorization with AMD or METIS ordering, threaded by OpenMP over
the elimination tree and by BLAS inside the fronts. Iterative refinement against
the original matrix makes results independent of the thread count up to
rounding.

```fortran
solver%backend_id = FORTSPARSE_BACKEND_MUMPS
solver%refine = .true.                          ! iterative refinement
solver%options%symmetry = FORTSPARSE_SYM_SPD    ! or _INDEFINITE (LDL^T), _GENERAL (LU)
solver%options%nthreads = 48                    ! 0: OMP_NUM_THREADS
call sparse_factor(solver, A, status)           ! full symmetric CSC; lower half used
call sparse_solve(solver, b, x, status)
call sparse_solve_block(solver, B, X, status)   ! B(n, nrhs)
call sparse_factor_stats(solver, stats, status) ! fill, memory, flops, times
```

An indefinite matrix declared `FORTSPARSE_SYM_SPD` fails with
`FORTSPARSE_SINGULAR`. Sizes are limited to 32-bit MUMPS integers.

Complex `csc_z_t` systems use ZMUMPS with the same calls:

| `options%symmetry` | Meaning for complex matrices | Factorization |
| --- | --- | --- |
| `FORTSPARSE_SYM_GENERAL` | any matrix | LU |
| `FORTSPARSE_SYM_HERMITIAN` | `A = A^H`, indefinite allowed | LU of the expanded matrix |
| `FORTSPARSE_SYM_SPD` | Hermitian positive definite | LU of the expanded matrix |
| `FORTSPARSE_SYM_INDEFINITE` | complex symmetric, `A = A^T` | LDL^T of the lower triangle |

For the Hermitian options the lower triangle is read, the upper triangle is
rebuilt by conjugation and the diagonal is taken as real. MUMPS has no Hermitian
factorization, so these cost a general LU. On a real matrix
`FORTSPARSE_SYM_HERMITIAN` means real symmetric (LDL^T). The complex path needs
a MUMPS installation with the `z` component; otherwise complex factorizations
return `FORTSPARSE_BACKEND_UNAVAILABLE`.

MUMPS keeps process-global state, so all MUMPS calls run under one process-wide
lock: independent solver objects may be used from different threads but execute
one at a time. For parallel speed-up use `options%nthreads` on one solver. An
MPI-enabled MUMPS can share the application's MPI runtime; the application must
initialize MPI before factoring and free every solver before finalizing MPI.

Building MUMPS without MPI, with OpenMP and METIS, using the `scivision/mumps`
CMake wrapper (add `-DBUILD_DOUBLE=ON -DBUILD_COMPLEX16=ON -DBUILD_SINGLE=OFF
-DBUILD_COMPLEX=OFF` for the complex path):

```sh
cmake -S mumps-cmake -B build-mumps -DMUMPS_parallel=OFF -DMUMPS_scalapack=OFF \
      -DMUMPS_openmp=ON -DMUMPS_metis=ON -DCMAKE_INSTALL_PREFIX=$PREFIX
cmake --build build-mumps && cmake --install build-mumps
cmake -S fortsparse -B build -DFORTSPARSE_ENABLE_MUMPS=ON \
      -DFORTSPARSE_USE_SYSTEM_DEPS=ON -DCMAKE_PREFIX_PATH=$PREFIX
```

Link an OpenMP-threaded BLAS (the OpenBLAS `openmp` variant) so BLAS calls inside
the parallel elimination tree stay serial per thread; this speeds up MUMPS by 2 to
6 times on large systems compared with reference BLAS.
`bench_fortsparse_mumps <p3|lap> <size> <threads>` measures factor and solve
times; `bench_fortsparse_backends <kind> <size> <threads> <superlu|mumps|auto>`
compares backends.

## UMFPACK

UMFPACK runs in a separate GPL helper process that communicates with the library
over shared memory and named semaphores (see [design](design.md) and
[licensing](licensing.md)).

A negative `solver%options%pivot_tolerance` keeps UMFPACK's defaults; a value in
`[0,1]` sets both the ordinary and the symmetric diagonal pivot thresholds.

`solver%options%memory_budget_bytes` defaults to zero (no limit). With a
positive budget the helper runs the symbolic analysis and rejects a larger
UMFPACK peak estimate with `FORTSPARSE_RESOURCE_LIMIT` before numeric
factorization. This is a check on UMFPACK's estimate, not a hard allocation cap,
and it excludes the caller's matrices, the IPC mapping and runtime overhead.
`sparse_factor_stats` reports `symbolic_peak_bytes`, `estimated_peak_bytes`,
`numeric_peak_bytes`, `ipc_bytes` and `numeric_started` also after a rejection;
unknown values are `-1`. These options are ignored by the other backends.
