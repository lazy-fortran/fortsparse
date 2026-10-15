# fortsparse

A backend-pluggable sparse direct solver for Fortran. It solves square real or
complex double systems `A x = b` through one public API, independent of the
library that does the factorization. Factor a matrix once, solve any number of
right-hand sides against that factorization, then free it. The backend is an
integer tag on the solver; switching it never changes the calling code.

The default backend is MUMPS when the library is built with a sequential MUMPS
that provides real and complex arithmetic (see [MUMPS backend](#mumps-backend)),
otherwise SuperLU (BSD), linked in-process through a small MIT C shim. UMFPACK (GPL) is reached through a second backend that never links into
the library: it runs in a standalone GPL helper process driven over shared
memory and named semaphores. The MIT library stays free of any GPL link. See
[License](#license) for why the process boundary holds.

## Build

CMake is the primary build, with CTest. By default it fetches and builds the
solver libraries it needs (SuperLU, and the SuiteSparse/UMFPACK chain) from
source, so the only system requirement is a BLAS/LAPACK. OpenBLAS or MKL both
work; on a cluster, `module load` it first.

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Options, all set with `-D`:

| Option | Default | Effect |
| --- | --- | --- |
| `FORTSPARSE_ENABLE_SUPERLU` | `ON` | Build the in-process SuperLU backend. `OFF` drops SuperLU entirely; the library then needs only the UMFPACK helper. |
| `FORTSPARSE_USE_SYSTEM_DEPS` | `OFF` | `ON` finds SuperLU and UMFPACK with `find_package` instead of fetching them. BLAS/LAPACK always come from the system. |
| `FORTSPARSE_ENABLE_MUMPS` | `OFF` | Build the in-process MUMPS backend (threaded supernodal Cholesky, LDL^T and LU). Needs an OpenMP-enabled MUMPS, found by `find_package(MUMPS CONFIG)` or else `cmake/FindMUMPS.cmake` (distribution packages); see [MUMPS backend](#mumps-backend). |
| `FORTSPARSE_INSTALL_HELPER` | `ON` | Install the UMFPACK helper to `bin/`. See [Use in your project](#use-in-your-project). |
| `FORTSPARSE_TEST_OOM` | `OFF` | Also register `test_fortsparse_superlu_status_oom` (label `oom`), which forces a real SuperLU allocation failure under `ulimit -v 250000`. Manual: the working limit is toolchain dependent. Run with `ctest -L oom`. |

The UMFPACK helper builds whenever the SuiteSparse chain is available (fetched,
or found under `FORTSPARSE_USE_SYSTEM_DEPS`). It is the one GPL artifact and the
one place `umfpack.h` appears.

An `fpm.toml` lets `fo` and `fpm` build the library too. That build sees only
`src/`, so it links SuperLU alone and never compiles the helper. Selecting the
UMFPACK backend with no helper present returns `FORTSPARSE_BACKEND_UNAVAILABLE`,
the documented fallback, not a crash.

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
MUMPS in a sequential MUMPS build with ZMUMPS and SuperLU otherwise; an MPI
MUMPS build keeps SuperLU because the application must own `MPI_Init`:

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

## MUMPS backend

`FORTSPARSE_BACKEND_MUMPS` is the backend for large SPD finite-element systems
(the P3 Grad-Shafranov stiffness matrix). It runs a multifrontal factorization
with METIS nested-dissection ordering, threaded by OpenMP over the elimination
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
the parallel elimination tree stay serial per thread. The ordering step is a
single-threaded METIS call; it runs once per factorization pattern.
`bench_fortsparse_mumps <p3|lap> <size> <threads>` measures factor and solve
times.

## Layout

- `src/`: library modules. Kinds, status, version, the CSC types and builders,
  the abstract backend and the pluggable solver, the in-process SuperLU backend
  (`src/superlu/`), and the out-of-process UMFPACK backend (`src/ipc/`). The MIT
  C shims sit beside their Fortran modules.
- `include/`: C headers shared by the library and the helper, the shared-memory
  protocol and the IPC ABI.
- `backends/mumps/`: the optional MUMPS binding (Fortran module and MIT C
  shims for real and complex, plus a link stub for double-only MUMPS). Outside `src/` so `fpm` and `fo` never compile it.
- `helper/`: `fortsparse_umfpack_helper.c`, the standalone GPL binary, the only
  file that includes `umfpack.h`. It lives outside `src/` so `fpm` and `fo`
  never compile it.
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
