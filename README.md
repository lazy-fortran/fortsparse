# fortsparse

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24223525.svg)](https://doi.org/10.5281/zenodo.24223525)

A sparse direct solver library for Fortran with interchangeable backends. It
solves square real or complex double-precision systems `A x = b` through one API:
factor a matrix once, solve any number of right-hand sides, then free it. The
backend is an integer tag on the solver object, so switching it does not change
the calling code.

Status: experimental research software; interfaces may change.
[Implementation status](docs/status.md) lists the independent checks and
unsupported paths.

## Features

- One API over several factorization backends: SuperLU (in process), MUMPS
  (optional, in process), UMFPACK (in a separate GPL helper process) and a
  project-owned native backend.
- Native CPU factorizations: real symmetric positive-definite (SPD) Cholesky,
  complex Hermitian positive-definite (HPD) Cholesky, real and
  complex general LU, and pivoted symmetric/Hermitian indefinite LDL^T/LDL^H.
  Approximate minimum degree (AMD) fill-reducing ordering, optional METIS
  ordering, OpenMP, and optional CUDA for real SPD. Default ordering depends
  on the matrix class and resource settings; see
  [native matrix classes](docs/native-cpu.md).
- Real and complex compressed sparse column (CSC) matrices (`csc_t`, `csc_z_t`),
  built from 1-based coordinate (COO) triplets; block solves, transpose/adjoint solves and same-pattern refactorization.
- Implicit tangent and adjoint derivatives of real and complex solves
  (`sparse_solve_jvp`, `sparse_solve_vjp`) that reuse retained factors.
- Host and device memory limits and a low-memory mode for the native solvers.
- Smoothed-aggregation algebraic multigrid (AMG), compressed sparse row (CSR)
  smoothers and an incomplete LU preconditioner with thresholding and pivoting
  (ILUTP).
- A small C/C++ API for native general LU, and an optional FreeFem plugin.

## Building

CMake is the main build. It needs Fortran 2018 and C compilers, CMake 3.22 or
newer, BLAS/LAPACK (OpenBLAS or MKL), and Ninja for the commands below. By
default it fetches and builds SuperLU and the SuiteSparse/UMFPACK chain from
source, and fetches [fortnum](https://github.com/lazy-fortran/fortnum).

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
```

The `fpm.toml` builds the native CPU library only, linked against system
BLAS/LAPACK. The compiler must pass the self-assignment lifetime probe;
CMake runs this check during configuration. For fpm:

```sh
fpm build
fpm test --target test_fortsparse_compiler_self_assignment
fpm test --target test_fortsparse_native_only
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

Check `status%code` against `FORTSPARSE_OK` after each call that accepts a
status argument. `sparse_solve_once(A, b, x, status)` factors, solves and frees
in one call. Complex systems use `csc_z_t` with the same calls.

`FORTSPARSE_BACKEND_AUTO` uses MUMPS for large matrices when it is built,
SuperLU for smaller ones, and the native backend when neither is available.

The [user guide](docs/user-guide.md) gives a complete first program and
explains status handling, factor ownership and derivative calls.

## Documentation

- [User guide](docs/user-guide.md), [module overview](docs/modules.md),
  [implementation status](docs/status.md)

- [Building and linking](docs/building.md)
- [Backends, selection and benchmarks](docs/backends.md)
- [Design: backend abstraction, CSC format, UMFPACK IPC](docs/design.md)
- [Native CPU matrix classes](docs/native-cpu.md),
  [memory limits](docs/memory.md), [precision and adaptive storage](docs/precision.md)
- [AMG and smoothers](docs/amg.md), [ILUTP preconditioner](docs/ilutp.md)
- [Native C/C++ API](docs/native-c-api.md),
  [FreeFem plugin](plugins/freefem/README.md)
- [Licensing of the optional backends](docs/licensing.md)
- [Contributing](CONTRIBUTING.md): checks and submission to the primary
  repository on [fortyard](https://code.tugraz.at/yard)

## Related projects

- [fortnum](https://github.com/lazy-fortran/fortnum): numerical library used for
  Krylov methods and precision policies.

## Citation

Cite fortsparse through its Zenodo record,
[doi:10.5281/zenodo.24223525](https://doi.org/10.5281/zenodo.24223525), which
always resolves to the latest release; release 0.1.0 is
[doi:10.5281/zenodo.24223526](https://doi.org/10.5281/zenodo.24223526).

Use [CITATION.cff](CITATION.cff) to cite the software version used in your work.

## License

MIT (see [LICENSE](LICENSE)). The UMFPACK helper is a separate GPL executable
reached over shared memory, and the
optional MUMPS backend links CeCILL-C code; see
[docs/licensing.md](docs/licensing.md).
