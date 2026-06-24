# Contributing to fortsparse

fortsparse is a backend-pluggable sparse direct solver. One public API factors
a matrix once, solves any number of right-hand sides, and frees the
factorization, independent of which library does the work. This document
describes how to contribute, what the milestone structure is, and what contract
a new backend or module must satisfy.

## The license rule

This is non-negotiable. The fortsparse library is MIT and links only permissive
code: SuperLU (BSD) and the MIT C shims.

**No GPL code or header may enter `src/` or the `fortsparse` target.** Do not
add `#include <umfpack.h>` or any GPL header to a file under `src/`. Do not link
`-lumfpack` or any GPL library into the `fortsparse` library target. Do not copy
GPL source into the tree.

UMFPACK, and any future GPL solver, is reached only from a separate standalone
helper executable under `helper/`, which links the GPL library and is itself
GPL. The library drives it across a process boundary over shared memory. A PR
that links a GPL library into `libfortsparse`, or that puts a GPL header under
`src/`, will not merge. `ldd` and `nm` on `libfortsparse` and the test binaries
must show no GPL solver symbol. See `docs/design.md` for the boundary.

## Milestones

| Milestone | Scope |
|-----------|-------|
| M0 | Infrastructure: CMake, fpm, CI, code style, the MIT/GPL split |
| M1 | CSC storage, status, the abstract backend and solver dispatch |
| M2 | SuperLU backend, in-process, the default |
| M3 | Out-of-process UMFPACK backend over shared-memory IPC and the GPL helper |
| M4 | Hardening: edge cases, large systems, portability across Linux/macOS/Windows |
| M5 | Own clean-room MIT LU to retire the UMFPACK dependency |
| M6 | Differentiability: Enzyme wiring, device-memory sharing |

Work in a feature branch. Name it `m<N>/<short-description>`, e.g.
`m2/superlu-refine`. Open a PR against `main` when tests pass.

## Build

CMake is the primary build system. It builds everything, including the GPL
helper when SuiteSparse is present, and runs the IPC round-trip:

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
ctest --test-dir build --output-on-failure
```

`fpm` also works via `fo`, and sees only `src/`, so it builds the MIT library
and SuperLU and never the helper. Run `fo` (no arguments) for the full
pipeline: static analysis, build, test, lint, format check.

## Code style

- Free source form. `implicit none` in every scoping unit.
- `use <mod>, only: ...` before `implicit none`.
- `real(dp)` via the `fortsparse_kinds` module.
- All dummy arguments have explicit `intent`. Declarations at start of scope.
- Derived-type names end in `_t`. `allocatable` over pointers; `move_alloc` for
  transfer.
- No module-level `save` or global mutable state. No alloc/dealloc in hot loops.
- `fprettify`: 88-column width, 4-space indent.
- Modules under 500 lines, procedures under 50 lines.
- Comments say why a choice was made, not what the next line does.
- `.and.` and `.or.` do not short-circuit in Fortran. Never rely on it: put a
  guard in its own `if`. `fo lint` flags short-circuit reliance; keep it at
  zero.

## Per-backend contract

Every backend must satisfy the following before a PR merges.

### 1. Implements the abstract type

A backend is a type that extends `sparse_backend_t` and implements all five
deferred bindings: `factor_real`, `factor_complex`, `solve_real`,
`solve_complex`, `free`. It holds the retained factorization in the type, so
repeated solves reuse it. Add the public `FORTSPARSE_BACKEND_<NAME>` id and the
factory `case` in `fortsparse_solver`.

### 2. Status mapping

Map the underlying library's native error codes onto the `fortsparse_status`
codes. A singular matrix sets `FORTSPARSE_SINGULAR`; an allocation or internal
failure sets `FORTSPARSE_INTERNAL_ERROR`. A backend that is not present in the
build reports `FORTSPARSE_BACKEND_UNAVAILABLE` cleanly, never a crash. Every
status carries a human-readable message.

### 3. In-process or out-of-process

A backend that links only permissive code (BSD, MIT) runs in-process, like
SuperLU. A backend that would link copyleft or proprietary code runs
out-of-process behind the IPC boundary, like UMFPACK, with the GPL part
confined to a separate helper under `helper/`.

### 4. Tests

Every backend requires:

- **Solver tests** under `test/solver/` covering a known small real and complex
  system, second-RHS reuse, the 1D Poisson reference problem, and the error
  paths (singular matrix, solve before factor, unknown backend id). Permissive
  in-process backends run these directly.
- For an out-of-process backend, an **IPC round-trip test** under `test/ipc/`
  registered through CMake only, with the helper path in the environment. It
  must not be an fpm auto-test that fails when the helper is absent. The
  round-trip validates the solution against the in-process backend to the same
  tolerance.

Tests are plain programs: a failure counter, check helpers, `write(error_unit,
...)` on failure, `stop 1` on failure, `"PASS"` and `stop 0` on success. They
must stay runnable without the helper, on the SuperLU backend. Each test runs
in under 120 seconds. Tests run via CTest and via fpm/fo. Do not merge if any
test fails.

## Pull requests

A PR must include:

- a reference to the issue it closes,
- for a new backend, the id, the status mapping, and in-process vs.
  out-of-process,
- real `ctest` output showing all tests pass, and `fo` output for the
  MIT-only build.

Do not open a PR that skips, weakens, or disables existing tests, that links a
GPL library into `libfortsparse`, or that puts a GPL header under `src/`.
