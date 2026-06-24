# fortsparse

A backend-pluggable sparse direct solver for Fortran. It solves square real or
complex double systems `A x = b` through one public API, independent of the
library that does the factorization. Factor a matrix once, solve any number of
right-hand sides against that factorization, then free it. The backend is an
integer tag on the solver; switching it never changes the calling code.

The default backend is SuperLU (BSD), linked in-process through a small MIT C
shim. UMFPACK (GPL) is reached through a second backend that never links into
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
| `FORTSPARSE_INSTALL_HELPER` | `ON` | Install the UMFPACK helper to `bin/`. See [Use in your project](#use-in-your-project). |

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

Select a backend before factoring; the default is SuperLU:

```fortran
solver%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC
```

`csc_z_t` and the same calls handle complex systems. `sparse_solve_once(A, b, x,
status)` factors, solves, and frees in one call. Every routine reports through
`status`; check `status%code` against `FORTSPARSE_OK`.

## Layout

- `src/`: library modules. Kinds, status, version, the CSC types and builders,
  the abstract backend and the pluggable solver, the in-process SuperLU backend
  (`src/superlu/`), and the out-of-process UMFPACK backend (`src/ipc/`). The MIT
  C shims sit beside their Fortran modules.
- `include/`: C headers shared by the library and the helper, the shared-memory
  protocol and the IPC ABI.
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

Building the helper is the user's choice. A SuperLU-only build never produces
it. A user who wants UMFPACK builds the helper and accepts the GPL for that one
binary; the MIT library it talks to is unaffected. See `docs/design.md` for the
process-boundary argument.
