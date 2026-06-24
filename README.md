# fortsparse

A unified, backend-pluggable sparse direct solver for Fortran. fortsparse
solves square real or complex double systems `A x = b` through one public API,
independent of which direct-solver library does the factorization. Client code
factors a matrix once, solves any number of right-hand sides reusing that
factorization, then frees it. The backend is an integer tag on the solver;
selecting a different one never changes the calling code.

The default backend is SuperLU (BSD), bound in-process through a small MIT C
shim. UMFPACK (GPL) is available through a separate backend, but it never
links into the library: it runs in a standalone GPL helper process that the
MIT library drives over shared memory and named semaphores. This keeps the
MIT library free of any GPL link. SuperLU is the default precisely because it
imposes no such obligation; the UMFPACK helper is opt-in.

The numerics are primal-first and derivative-ready: the value path is kept
explicit so that automatic differentiation (planned Enzyme integration) and
device-memory sharing can be added later without rewriting the solver.
fortsparse may later depend on the sibling library fortnum for shared kinds and
numerical utilities.

## Build

CMake is the primary build system with CTest integration. It builds everything,
including the GPL UMFPACK helper when SuiteSparse is present. On Debian/Ubuntu
install `libsuperlu-dev` for the default backend and `libsuitesparse-dev` for
the helper.

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
ctest --test-dir build --output-on-failure
```

An `fpm.toml` is provided so the `fo` tool and `fpm` work as well. This build
sees only `src/`, so it links SuperLU only and never touches UMFPACK; the
helper is not built. Run `fo` with no arguments for the full pipeline: static
analysis, build, test, lint, format check.

```
fo
```

```
fpm build
fpm test
```

The IPC round-trip against the UMFPACK helper runs only under the CMake build,
which is the only build that produces the helper binary. The fpm/fo build
reports `FORTSPARSE_BACKEND_UNAVAILABLE` if a client selects the UMFPACK
backend with no helper present; that is the documented behavior, not an error.

## Layout

- `src/`: library modules. Kinds, status, version, the CSC storage types and
  builders, the abstract backend, the backend-pluggable solver, the in-process
  SuperLU backend under `src/superlu/`, and the out-of-process UMFPACK backend
  under `src/ipc/`. The MIT C shims (SuperLU bridge, IPC primitives) sit beside
  their Fortran modules.
- `include/`: the C headers shared by the library and the helper, the shared
  memory protocol and the IPC primitive ABI.
- `helper/`: `fortsparse_umfpack_helper.c`, the standalone GPL binary. It is
  the only file that includes `umfpack.h`. It is outside `src/`, so fpm and fo
  never compile it.
- `test/`: CTest suite, with `csc/`, `solver/`, and `ipc/` subgroups.
- `docs/`: `design.md` for the design notes and the roadmap.
- `cmake/`: CMake helper modules, `FindSuperLU.cmake` and `FindUMFPACK.cmake`.

## License

The fortsparse library is MIT, Copyright (c) 2025 lazy-fortran. It links only
permissive code: SuperLU (BSD) and the MIT C shims. It contains no GPL source,
includes no GPL header, and links no GPL library. Linking the library, static
or shared, imposes no copyleft obligation.

The UMFPACK helper, `helper/fortsparse_umfpack_helper.c`, is a separate
standalone executable. It links UMFPACK and is therefore GPL. It is the only
GPL artifact in the tree, and it is the only place `umfpack.h` appears. The
library talks to it across a process boundary over shared memory; a separate
process behind that boundary is the FSF-recognized separation that keeps the
GPL of UMFPACK confined to the helper binary.

Building and using the helper is the user's choice. The default SuperLU build
never produces it. A user who wants UMFPACK builds the helper through CMake and
accepts the GPL for that binary; the MIT library it talks to is unaffected. See
`docs/design.md` for why the process boundary is clean.
