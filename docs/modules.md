# Module overview

Most applications import `fortsparse`, which re-exports the matrix, solver,
status and option interfaces. Import narrower modules when their ownership or
build boundary matters.

| Module | Responsibility |
| --- | --- |
| `fortsparse_kinds` | Real and integer kinds |
| `fortsparse_status` | Public status codes and status helpers |
| `fortsparse_csc`, `fortsparse_csr` | Column/row sparse storage, builders and matrix actions |
| `fortsparse_solver` | Retained factors, backend selection, solve/refactor and implicit derivatives |
| `fortsparse_backend` | Backend abstraction, options and factor statistics |
| `fortsparse_memory` | Solver-owned memory accounting and planning |
| `fortsparse_mixed_lu` | Mixed-precision factors and refinement |
| `fortsparse_amg`, `fortsparse_smoothers` | Algebraic multigrid and sparse smoothers |
| `fortsparse_ilutp` | Thresholded incomplete LU with pivoting |

Native factors live in `src/native/`. Optional SuperLU bindings live in
`src/superlu/`; the UMFPACK interprocess client is under `src/ipc/`, while its
separate GPL executable is under `helper/`. MUMPS is outside the fpm source tree
in `backends/mumps/`. The [licensing guide](licensing.md) explains these build
boundaries.

The C interface is declared by `include/fortsparse_native.h`; its layouts and
status contract are in [the native C/C++ guide](native-c-api.md).
