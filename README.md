# fortsparse

A sparse direct solver library for Fortran with interchangeable backends. It
solves square real or complex double-precision systems `A x = b` through one API:
factor a matrix once, solve any number of right-hand sides, then free it. The
backend is an integer tag on the solver object, so switching it does not change
the calling code.

Status: research code under active development; interfaces may change.

## Features

- One API over several factorization backends: SuperLU (in process), MUMPS
  (optional, in process), UMFPACK (in a separate GPL helper process) and a
  project-owned native backend.
- Native CPU factorizations: real SPD Cholesky, complex HPD Cholesky, real and
  complex general LU, and pivoted symmetric/Hermitian indefinite LDL^T. Built-in
  AMD fill-reducing ordering by default, optional METIS ordering, OpenMP, and
  optional CUDA for real SPD.
- Real and complex CSC matrices (`csc_t`, `csc_z_t`), built from 1-based COO
  triplets; block solves, transpose/adjoint solves and same-pattern refactorization.
- Exact implicit derivatives of real solves (`sparse_solve_jvp`, `sparse_solve_vjp`)
  that reuse the primal factorization.
- Host and device memory limits and a low-memory mode for the native solvers.
- Smoothed-aggregation AMG, CSR smoothers and an ILUTP preconditioner.
- A small C/C++ API for native general LU, and an optional FreeFem plugin.

## Building

CMake is the main build. By default it fetches and builds SuperLU and the
SuiteSparse/UMFPACK chain from source, and fetches
[fortnum](https://github.com/lazy-fortran/fortnum); the only system requirement is
BLAS/LAPACK (OpenBLAS or MKL).

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The `fpm.toml` builds the native CPU library only, linked against system
BLAS/LAPACK:

```sh
fpm build
fpm test
```

Build options (MUMPS, CUDA, AMD/METIS, system dependencies, the UMFPACK
helper) and how to link fortsparse from another CMake project are described in
[docs/building.md](docs/building.md).

## Usage

```fortran
use fortsparse
type(csc_t)               :: A
type(sparse_solver_t)     :: solver
type(fortsparse_status_t) :: status
real(dp) :: b(n), x(n)

! Build A from COO triplets; row and column indices are 1-based.
call csc_from_triplet(n, n, rows, cols, vals, A, status)

! Optional: choose a backend (default FORTSPARSE_BACKEND_AUTO).
solver%backend_id = FORTSPARSE_BACKEND_NATIVE
solver%options%symmetry = FORTSPARSE_SYM_SPD

! Factor once, solve any number of right-hand sides, then free.
call sparse_factor(solver, A, status)
call sparse_solve(solver, b, x, status)
call sparse_free(solver)
```

Every routine reports through `status`; compare `status%code` with
`FORTSPARSE_OK`. `sparse_solve_once(A, b, x, status)` factors, solves and frees
in one call. Complex systems use `csc_z_t` with the same calls.

`FORTSPARSE_BACKEND_AUTO` uses MUMPS for large matrices when it is built,
SuperLU for smaller ones, and the native backend when neither is available.

## Documentation

- [Building and linking](docs/building.md)
- [Backends, selection and benchmarks](docs/backends.md)
- [Design: backend abstraction, CSC format, UMFPACK IPC](docs/design.md)
- [Native CPU matrix classes](docs/native-cpu.md),
  [memory limits](docs/memory.md), [precision and adaptive storage](docs/precision.md)
- [AMG and smoothers](docs/amg.md), [ILUTP preconditioner](docs/ilutp.md)
- [Native C/C++ API](docs/native-c-api.md),
  [FreeFem plugin](plugins/freefem/README.md)
- [Licensing of the optional backends](docs/licensing.md)

## Related projects

- [fortnum](https://github.com/lazy-fortran/fortnum): numerical library used for
  Krylov methods and precision policies.

## License

MIT (see [LICENSE](LICENSE)). The library links only permissive code. The
UMFPACK helper is a separate GPL executable reached over shared memory, and the
optional MUMPS backend links CeCILL-C code; see
[docs/licensing.md](docs/licensing.md).
