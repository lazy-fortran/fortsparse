# fortsparse

A backend-pluggable sparse direct solver for Fortran. It solves square real or
complex double systems `A x = b` through one public API, independent of the
library that does the factorization. Factor a matrix once, solve any number of
right-hand sides against that factorization, then free it. The backend is an
integer tag on the solver; switching it never changes the calling code.

A small [native C/C++ API](docs/native-c-api.md) exposes owned general-LU
contexts with real/complex N/T/H block solves and explicit failure status.

The [roadmap](ROADMAP.md) defines the native single-node solver development
track, generated CPU/GPU kernels and performance/memory admission gates.
The next requirements are native CPU coverage across audited ITPplasma matrix
classes and predictable, user-limited memory with RAM-only low-memory modes.
Native CPU references and natural-order SPD blocked execution support strict
host payload caps; exact pre-factor planning initially covers compact real SPD.
General LU also supports bounded explicit AMD through an audited caller-owned
workspace. See [memory limits and cluster sizing](docs/memory.md). Production
planning, SPD ordering/GPU caps and broader performance admission remain planned.
Development continues on `main`; GPU usability precedes CPU coverage, with full GPU Pareto
admission retained as a later gate. MPI and disk-backed factors are out of scope.

The [measured real-SPD MVP](bench/results/2026-10-16-mvp.md) is on the observed
CPU time/memory Pareto front for all 115 frozen case/workload groups, including
held-out coefficients and HB/bcsstk09. This uses explicit measured ordering and
thread configurations. CUDA correctness is verified; GPU performance is not
admitted on this corpus. The report retains ranges, CPU costs and limitations.

`FORTSPARSE_BACKEND_NATIVE` selects the developing project-owned real SPD
Cholesky solver. Set `solver%options%symmetry = FORTSPARSE_SYM_SPD` and
supply full symmetric CSC with strictly sorted, unique row indices and finite
values; `csc_from_triplet` canonicalizes duplicates and ordering. The native
backend uses a simplicial path for low fill and blocked OpenMP tasks for larger fronts.
Optional BSD AMD provides default ordering for larger non-banded systems,
with METIS nested dissection as an alternative. Explicit
`FORTSPARSE_ORDER_AMD` and `FORTSPARSE_ORDER_ND` require the corresponding
ordering library. The solver returns
`FORTSPARSE_NOT_POSITIVE_DEFINITE` for a
nonpositive or nonfinite pivot. It retains analysis for same-pattern
refactorization and supports block solves and the existing implicit JVP/VJP.
It currently reports factor/storage statistics, not a rigorous forward-error
certificate. AUTO prefers eligible external backends when they are built.
`FORTSPARSE_BACKEND_NATIVE_REFERENCE` explicitly selects the serial simplicial
reference for independent comparisons. Statistics distinguish requested and
actual active factorization threads.

Build with `FORTSPARSE_ENABLE_NATIVE_CUDA=ON` to enable explicit native GPU
factorization, and set `solver%options%device = 0` for CUDA device 0. Panels
and Schur updates remain resident during factorization. Generated triangular
solves reuse device factors, with RHS batches of eight and no mandatory factor
download. Host-array interfaces transfer RHS and results; the context and
factors are retained for solves and numeric refactorization. The default device is `-1` (CPU); a GPU request in a CPU build
returns backend unavailable. Device statistics cover attributed allocations,
excluding unreported driver/vendor internal pools. GPU correctness is tested
on real hardware; performance admission remains case dependent.

Native CPU also provides real/complex general sparse LU, complex HPD Cholesky,
and pivoted real symmetric, complex Hermitian and complex symmetric indefinite
factors as bounded-memory reference families.
See [matrix classes and transpose/adjoint solves](docs/native-cpu.md). Their
production performance admission and parallel factorization remain open.
The CPU references support retained original-matrix refinement; the published CPU frontier claim applies only to the frozen real-SPD corpus.

The standalone native real CSR [`sparse_ilutp_t` preconditioner](docs/ilutp.md)
provides explicit lower/upper fill caps, column pivoting, a host array-payload
limit and allocation-free retained normal apply. Krylov remains in FortNum.
SPEC application migration and preconditioner performance admission remain pending.

The default backend, `FORTSPARSE_BACKEND_AUTO`, chooses per matrix between
SuperLU (BSD) for small systems and MUMPS for large ones when the library is
built with a sequential MUMPS (see [Backend selection](#backend-selection)).
With SuperLU enabled and MUMPS disabled, AUTO uses SuperLU through a small MIT
C shim. When neither eligible in-process external backend is built, AUTO selects the
native backend. This capability fallback does not extend the performance
admission beyond the declared corpus. UMFPACK (GPL) is reached through a second backend that never links into
the library: it runs in a standalone GPL helper process driven over shared
memory and named semaphores. The MIT library stays free of any GPL link. See
[License](#license) for why the process boundary holds.

## Build

CMake is the primary build, with CTest. By default it fetches and builds the
solver libraries it needs (SuperLU, and the SuiteSparse/UMFPACK chain) from
source, so the only system requirement is a BLAS/LAPACK. OpenBLAS or MKL both
work; on a cluster, `module load` it first.

CMake selects BLAS/LAPACK before loading solver packages. With no parent or
`BLA_VENDOR` choice it tries OpenBLAS, sequential MKL, then CMake's normal
system search (including BLIS); `-DBLA_VENDOR=Generic` retains reference
BLAS. Fetched solvers share this choice. A prebuilt system solver may already
depend on another BLAS: inspect `ldd` on the final executable and rebuild that
solver against the selected BLAS, or use the default fetched dependencies, to
avoid mixing implementations. For pthread OpenBLAS with OpenMP MUMPS, set
`OPENBLAS_NUM_THREADS=1` and control sparse parallelism with `OMP_NUM_THREADS`;
this prevents nested BLAS teams. The fpm path links the system dense BLAS/LAPACK directly; their compiled
threading configuration remains a property of those libraries.

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Options, all set with `-D`:

| Option | Default | Effect |
| --- | --- | --- |
| `FORTSPARSE_ENABLE_SUPERLU` | `ON` | Build the in-process SuperLU backend. `OFF` drops SuperLU entirely; native execution remains available. |
| `FORTSPARSE_ENABLE_UMFPACK_IPC` | `ON` | Provision the separate GPL UMFPACK helper. `OFF` performs no UMFPACK/SuiteSparse search or fetch and builds/installs no helper. The IPC protocol binding remains available for an explicitly supplied helper. |
| `FORTSPARSE_USE_SYSTEM_DEPS` | `OFF` | `ON` finds SuperLU and UMFPACK with `find_package` instead of fetching them. BLAS/LAPACK always come from the system. |
| `FORTSPARSE_ENABLE_MUMPS` | `OFF` | Build the in-process MUMPS backend (threaded supernodal Cholesky, LDL^T and LU). Needs an OpenMP-enabled MUMPS, found by `find_package(MUMPS CONFIG)` or else `cmake/FindMUMPS.cmake` (distribution packages); see [MUMPS backend](#mumps-backend). |
| `FORTSPARSE_INSTALL_HELPER` | `ON` | Install the UMFPACK helper to `bin/`. See [Use in your project](#use-in-your-project). |
| `FORTSPARSE_TEST_OOM` | `OFF` | Also register `test_fortsparse_superlu_status_oom` (label `oom`), which forces a real SuperLU allocation failure under `ulimit -v 250000`. Manual: the working limit is toolchain dependent. Run with `ctest -L oom`. |

With `FORTSPARSE_ENABLE_UMFPACK_IPC=ON`, the UMFPACK helper builds whenever
the SuiteSparse chain is available (fetched, or found under
`FORTSPARSE_USE_SYSTEM_DEPS`). Optional benchmark builds can also produce
a separate GPL UMFPACK comparator; neither target links into the MIT library.

A native-only CMake build needs BLAS/LAPACK and no sparse solver, helper or CAS:

```sh
cmake -S . -B build-native -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DFORTSPARSE_ENABLE_SUPERLU=OFF -DFORTSPARSE_ENABLE_UMFPACK_IPC=OFF \
  -DFORTSPARSE_ENABLE_MUMPS=OFF -DFORTSPARSE_ENABLE_AMD=OFF \
  -DFORTSPARSE_ENABLE_METIS=OFF
cmake --build build-native --target test_fortsparse_native_only -j 2
ctest --test-dir build-native -R '^test_fortsparse_native_only$' --output-on-failure
```

AMD/METIS ordering and native CUDA can be enabled separately when available.
The `test_fortsparse_native_only_build` check (label `build`) separately checks
configuration, compilation, real/complex native AUTO solves, caps and linked
runtime dependencies. Its build time is not an ordinary numerical-test timing.

An `fpm.toml` lets `fo` and `fpm` build the native CPU library too. It links
system dense BLAS/LAPACK, without SuperLU, MUMPS, AMD/METIS or a UMFPACK helper.
AUTO therefore selects the native CPU backend. Select optional external sparse
backends, ordering libraries and CUDA through CMake. An explicitly requested
unconfigured backend returns `FORTSPARSE_BACKEND_UNAVAILABLE`.

The fpm test inventory remains enabled. Its three dedicated SuperLU programs
check the unavailable-backend status and report a skip when SuperLU is absent;
the enabled CMake configuration runs their complete numerical and status tests.
Standalone C test mains live under `cmake/tests` because fpm treats C sources
under `test` as shared support objects for each Fortran test executable.

```
fo
```

## Use in your project

Fetch fortsparse and link it:

```cmake
include(FetchContent)
FetchContent_Declare(fortsparse
  GIT_REPOSITORY https://github.com/lazy-fortran/fortsparse.git
  GIT_TAG main)
FetchContent_MakeAvailable(fortsparse)

target_link_libraries(my_app PRIVATE fortsparse)
```

The UMFPACK backend spawns its helper process and finds the binary next to the
running executable. Two ways to put it there:

- Install `my_app` to the standard `bin/` directory. With `FORTSPARSE_INSTALL_HELPER`
  on (the default), the helper installs to the same `bin/`, so `make install`
  co-locates them and the solver finds it. Nothing else to write.
- Install elsewhere, or run from the build tree: call `fortsparse_colocate_helper`
  once per executable that issues sparse solves.

```cmake
# Copy the helper next to my_app after each build, and install it to <dir>.
fortsparse_colocate_helper(my_app ${CMAKE_INSTALL_PREFIX})
```

`fortsparse_colocate_helper(<target> [<install_dir>])` adds a POST_BUILD copy
of the helper into the target's directory, and installs it to `<install_dir>`
when given. It is a no-op on a SuperLU-only build, so calling it
unconditionally is safe. The override `FORTSPARSE_UMFPACK_HELPER=<path>` names
the helper explicitly when neither layout applies.

## API

```fortran
use fortsparse
type(csc_t)               :: A
type(sparse_solver_t)     :: solver
type(fortsparse_status_t) :: status
real(dp) :: b(n), x(n)

! Build A from COO triplets; row and column indices are 1-based.
call csc_from_triplet(n, n, rows, cols, vals, A, status)

! Factor once, solve any number of right-hand sides, then free.
call sparse_factor(solver, A, status)
call sparse_solve(solver, b, x, status)
call sparse_free(solver)
```

Select a backend before factoring. The default, `FORTSPARSE_BACKEND_DEFAULT`, is
`FORTSPARSE_BACKEND_AUTO` (see [Backend selection](#backend-selection));
`FORTSPARSE_BACKEND_SUPERLU`, `FORTSPARSE_BACKEND_MUMPS` and the out-of-process
UMFPACK backend can be named explicitly:

```fortran
solver%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC
```

`csc_z_t` and the same calls handle complex systems. `sparse_solve_once(A, b, x,
status)` factors, solves, and frees in one call. Every routine reports through
`status`; check `status%code` against `FORTSPARSE_OK`.

Real solves provide exact implicit derivatives without differentiating the
backend. For \(A x=b\), `sparse_solve_jvp` reuses the primal factorization to
solve
\(A\dot{x}=\dot{b}-\dot{A}x\). `sparse_solve_vjp` solves
\(A^T\bar{b}=\bar{x}\) through a retained factorization of `transpose(A)` and
returns the active CSC-value cotangents
\(\bar{A}_{ij}=-\bar{b}_i x_j\). Thus the derivative path works with GCC and
every real direct backend:

```fortran
call sparse_solve_jvp(solver, A_dot, x, b_dot, x_dot, status)

call csc_transpose(A, transpose_A, status)
call sparse_factor(transpose_solver, transpose_A, status)
call sparse_solve_vjp(transpose_solver, A, x, x_bar, &
                      b_bar, A_values_bar, status)
```

## SuperLU symmetry hints

SuperLU factors the full supplied CSC matrix for every symmetry hint. The
`FORTSPARSE_SYM_SPD`, `FORTSPARSE_SYM_INDEFINITE` and
`FORTSPARSE_SYM_HERMITIAN` hints select symmetric mode, minimum-degree ordering
on the structure of `A + transpose(A)` and diagonal threshold pivoting
(`DiagPivotThresh = 0.001`). Off-diagonal pivots remain available for indefinite
systems. The default `FORTSPARSE_SYM_GENERAL` retains COLAMD and ordinary partial
pivoting. These hints do not check symmetry or positive definiteness; callers
must declare the matrix structure correctly.

## Backend selection

`FORTSPARSE_BACKEND_AUTO` resolves at each factorization from the number of
stored entries: MUMPS for `nnz >= 40000`, SuperLU below
(`sparse_auto_backend(nnz, is_complex)` returns the choice). MUMPS is chosen
only in a sequential build (an MPI MUMPS needs the application to own
`MPI_Init`) and, for complex matrices, only with ZMUMPS; a build without MUMPS
always uses SuperLU. A same-pattern `sparse_refactor` resolves to the same
backend, so the MUMPS analysis is reused.

The crossover comes from `bench_fortsparse_backends` (built with
`FORTSPARSE_ENABLE_MUMPS=ON`) on finite-element-shaped matrices: 2D P1/P2/P3
SPD stiffness, `hcurl` (complex edge + nodal unknowns, indefinite), `blk3` (P2
with three coupled fields, nonsymmetric LU), `ideal` (the same, symmetric
indefinite LDL^T) and 3D P1 tetrahedra with one (`tet1`, SPD) and three
(`tet3`, LU) fields. Times are factor + one solve + destroy, the cost of
`sparse_solve_once`; the last column is one solve with a retained factor.
AMD Ryzen 9 5950X, MUMPS 5.9.1 (OpenMP, METIS), SuperLU 7.0.1 (serial),
reference BLAS unless noted, other load on the machine:

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

SuperLU wins below about 3e4 entries, where its solves are also 3 to 15 times
faster; from about 5e4 entries MUMPS is as fast or faster, increasingly with
size, threads and an optimized BLAS (with OpenBLAS, 2.5 to 9 times faster at
1e5 entries and above).
SuperLU here is the serial library; no threaded SuperLU (SuperLU_MT or
SuperLU_DIST) is built.

MUMPS is tuned for these sizes. Its OpenMP tree parallelism (ICNTL(48), L0-OMP)
parks finished threads in `usleep(1000)` (MUMPS 5.9 `DMUMPS_PERFORM_COPIES`),
which cost about 1.1 ms per factorization with 4 threads (1.2 ms instead of
0.07 ms for n = 9); fortsparse enables it only from n = 5000. The default
ordering analyses with AMD and switches to METIS nested dissection when the
AMD estimate exceeds 3000 elimination flops per entry (three-dimensional or
strongly coupled patterns): METIS costs up to ten times the AMD analysis and
gains little fill on 2D meshes, but halves or quarters the factorization in
3D. `FORTSPARSE_ORDER_AMD` and `FORTSPARSE_ORDER_ND` still force either.

### UMFPACK pivot and memory controls

For `FORTSPARSE_BACKEND_UMFPACK_IPC`, negative
`solver%options%pivot_tolerance` keeps UMFPACK's native defaults. A value in
`[0,1]` sets both its ordinary and symmetric diagonal pivot thresholds; `1`
requests partial pivoting without favoring a smaller diagonal entry.

`solver%options%memory_budget_bytes` defaults to zero (no gate). With a
positive budget the helper runs symbolic analysis, then rejects a larger
UMFPACK peak estimate with `FORTSPARSE_RESOURCE_LIMIT` **before numeric
factorization**. This is an estimate gate, not a hard allocation cap. It
excludes caller matrices, CSC arrays, IPC mapping and operating-system/runtime
costs; callers must reserve these separately. A compatible helper is negotiated
before sending any factor request, so a stale helper cannot bypass the gate.

`sparse_factor_stats` also works after this rejection. Its
`symbolic_peak_bytes`, `estimated_peak_bytes` and `numeric_peak_bytes` are
UMFPACK's reported Unit counts converted to bytes; `ipc_bytes` is the separate
shared mapping allocation. `numeric_started` distinguishes rejected attempts;
`effective_pivot_tolerance` records the ordinary threshold passed to UMFPACK.
Unknown statistics remain `-1`. A rejected factor leaves the solver unfactored;
a subsequent call with a sufficient budget can reuse the solver normally.
These options are currently ignored by other backends.

## MUMPS backend

`FORTSPARSE_BACKEND_MUMPS` is the backend for large SPD finite-element systems
(the P3 Grad-Shafranov stiffness matrix). It runs a multifrontal factorization
with AMD or METIS nested-dissection ordering, threaded by OpenMP over the elimination
tree and by BLAS inside the fronts. One factorization serves any number of
right-hand sides, a block solve handles many at once, and iterative refinement
against the original matrix makes the result thread-invariant to rounding.

```fortran
solver%backend_id = FORTSPARSE_BACKEND_MUMPS
solver%refine = .true.                       ! iterative refinement
solver%options%symmetry = FORTSPARSE_SYM_SPD ! or _INDEFINITE (LDL^T), _GENERAL (LU)
solver%options%nthreads = 48                 ! 0: OMP_NUM_THREADS
call sparse_factor(solver, A, status)        ! full symmetric CSC; lower half used
call sparse_solve(solver, b, x, status)      ! reuse the factor
call sparse_solve_block(solver, B, X, status)! B(n, nrhs) in one call
call sparse_factor_stats(solver, stats, status) ! fill, memory, flops, times
```

An indefinite matrix declared `FORTSPARSE_SYM_SPD` fails with
`FORTSPARSE_SINGULAR`. Matrix size and entry count are limited to 32-bit MUMPS
integers.

### Complex matrices (ZMUMPS)

Complex `csc_z_t` systems take the same calls, with the same factor-once,
solve-many and block-solve shape: `sparse_factor(solver, Az, status)`,
`sparse_solve(solver, bz, xz, status)` (also in place) and
`sparse_solve_block(solver, Bz, Xz, status)`. Iterative refinement
(`solver%refine`) runs against the original complex matrix. The symmetry
option selects the factorization:

| `options%symmetry` | Complex meaning | Factorization |
| --- | --- | --- |
| `FORTSPARSE_SYM_GENERAL` | any matrix, as given | LU |
| `FORTSPARSE_SYM_HERMITIAN` | `A = A^H`, indefinite allowed | LU of the expanded matrix |
| `FORTSPARSE_SYM_SPD` | Hermitian positive definite (same as Hermitian) | LU of the expanded matrix |
| `FORTSPARSE_SYM_INDEFINITE` | complex symmetric, `A = A^T` (not `A^H`) | LDL^T of the lower triangle |

For the Hermitian options the lower triangle of the supplied full CSC matrix is
read, the upper triangle is rebuilt by conjugation and the diagonal is taken
real. MUMPS has no Hermitian factorization, so these options cost a general LU;
they remain valid for indefinite matrices. `FORTSPARSE_SYM_HERMITIAN` on a real
matrix means real symmetric (LDL^T). A real solve against a complex
factorization (or the reverse) fails with `FORTSPARSE_INTERNAL_ERROR`.

The complex path is built when the MUMPS installation provides the `z`
component (`MUMPS_BUILD_COMPLEX16=ON`, i.e. `-DBUILD_COMPLEX16=ON` in the
`scivision/mumps` wrapper, alongside `-DBUILD_DOUBLE=ON`). Against a
double-only installation the library still builds, and complex factorizations
return `FORTSPARSE_BACKEND_UNAVAILABLE`.

Thread safety: MUMPS keeps process-global state, so every MUMPS call (real and
complex, all solver objects) runs under one process-wide lock. Independent
solver objects may be used from different OpenMP threads without corrupting
anything, but they then execute one at a time. One solver object must not be
shared between threads without the caller ordering its factor and solve calls.
For parallel speed-up use `options%nthreads` (MUMPS' own OpenMP threads) on one
solver object rather than user threads.

An MPI-enabled MUMPS can share the application's MPI runtime. The application
must initialize MPI before factorization and free every solver before finalizing
MPI; inactive MPI returns `FORTSPARSE_BACKEND_UNAVAILABLE`. Use this build when
the application already links MPI, to avoid the sequential MPI shim's symbols.

Build MUMPS without MPI and with OpenMP and METIS, for example with the
`scivision/mumps` CMake wrapper (it fetches the MUMPS source and METIS), then
point fortsparse at the install prefix (add `-DBUILD_DOUBLE=ON
-DBUILD_COMPLEX16=ON -DBUILD_SINGLE=OFF -DBUILD_COMPLEX=OFF` for the complex
path):

```
cmake -S mumps-cmake -B build-mumps -DMUMPS_parallel=OFF -DMUMPS_scalapack=OFF \
      -DMUMPS_openmp=ON -DMUMPS_metis=ON -DCMAKE_INSTALL_PREFIX=$PREFIX
cmake --build build-mumps && cmake --install build-mumps
cmake -S fortsparse -B build -DFORTSPARSE_ENABLE_MUMPS=ON \
      -DFORTSPARSE_USE_SYSTEM_DEPS=ON -DCMAKE_PREFIX_PATH=$PREFIX
```

Link an OpenMP-threaded BLAS (OpenBLAS `openmp` variant) so BLAS calls inside
the parallel elimination tree stay serial per thread; it speeds MUMPS up 2 to 6
times on large systems over the reference BLAS. The ordering step is
single-threaded; it runs once per factorization pattern.
`bench_fortsparse_mumps <p3|lap> <size> <threads>` measures factor and solve
times; `bench_fortsparse_backends <kind> <size> <threads> <superlu|mumps|auto>`
compares the backends.

## Layout

- `src/`: library modules. Kinds, status, version, the CSC types and builders,
  the abstract backend and the pluggable solver, the in-process SuperLU backend
  (`src/superlu/`), and the out-of-process UMFPACK backend (`src/ipc/`). The MIT
  C shims sit beside their Fortran modules.
- `include/`: C headers shared by the library and the helper, the shared-memory
  protocol and the IPC ABI.
- `backends/mumps/`: the optional MUMPS binding (Fortran module and MIT C
  shims for real and complex, plus a link stub for double-only MUMPS). Outside `src/` so `fpm` and `fo` never compile it.
- `helper/`: `fortsparse_umfpack_helper.c`, the standalone GPL product helper.
  Like the optional GPL benchmark under `bench/native_families/umfpack_direct/`,
  it lives outside `src/` and is never linked into the MIT library.
- `test/`: CTest suite, grouped into `csc/`, `solver/`, and `ipc/`.
- `docs/`: `design.md`, design notes plus the planned autodiff and device work.
- `cmake/`: helper modules, `FindSuperLU.cmake`, `FindUMFPACK.cmake`.

## License

The fortsparse library is MIT, Copyright (c) 2025 lazy-fortran. It links only
permissive code, SuperLU (BSD) and the MIT C shims. It contains no GPL source,
includes no GPL header, and links no GPL library. Linking it, static or shared,
imposes no copyleft obligation.

The UMFPACK helper, `helper/fortsparse_umfpack_helper.c`, is a separate
executable. It links UMFPACK and is therefore GPL. The library talks to it
across a process boundary over shared memory. A separate process behind that
boundary is the FSF-recognized separation that keeps the GPL of UMFPACK
confined to the helper binary.

The optional MUMPS backend links MUMPS (CeCILL-C, a weak-copyleft licence in
the LGPL family) in-process. fortsparse contains no MUMPS source; only the MIT
shim includes its public header. MUMPS stays under CeCILL-C wherever a build
redistributes it; the MIT terms of fortsparse do not change. A build with
`FORTSPARSE_ENABLE_MUMPS=OFF` carries no such term.

Building the helper is the user's choice. A SuperLU-only build never produces
it. A user who wants UMFPACK builds the helper and accepts the GPL for that one
binary; the MIT library it talks to is unaffected. See `docs/design.md` for the
process-boundary argument.

The optional `bench_native_families_umfpack` comparator directly links UMFPACK
in its own GPL-2.0-or-later executable. Its Fortran/C driver sources and license
are kept under `bench/native_families/umfpack_direct/`; shared fixture and record
sources remain MIT. It does not link FortSparse, replace the IPC backend or
change the library's license. See that directory's README for ABI and memory
measurement limits.
