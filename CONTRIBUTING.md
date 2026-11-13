# Contributing to FortSparse

Start with [building and linking](docs/building.md), the
[module overview](docs/modules.md) and [implementation status](docs/status.md).
Keep each change focused and preserve existing matrix-class, ownership and
status-code contracts. A numerical change needs an analytical or independent
behavioral test in the existing `test/` tree.

CMake is the primary build and test authority:

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
```

The native CPU-only fpm build links BLAS/LAPACK. Use `fpm test --target NAME`
for affected tests. See [the build guide](docs/building.md) for optional
backends, system dependencies and compiler limitations. Run the configured
checks for every affected backend; a successful native-only build does not
verify an external backend.

The library is MIT licensed. Never add GPL source or headers to `src/`, or
link UMFPACK or another GPL solver into the `fortsparse` target. UMFPACK stays
in the separate helper under `helper/`, reached through interprocess
communication. Verify the link boundary when changing backend dependencies;
see [licensing](docs/licensing.md).

Each backend implements factor, solve and free operations for its advertised
real and complex matrix classes. Map native errors to the public status codes,
including singular matrices and unavailable backends. Tests cover known
solutions, original-matrix residuals, retained-factor reuse and error paths.
Solve derivatives are checked with central differences and adjoint identities.
Do not weaken these checks or silently fall back from an explicit backend.

Use explicit imports and argument intents, allocatable ownership and the
existing four-space, 88-column style. Put reusable numerical mathematics in
`gen/` and regenerate its CPU/CUDA outputs; do not hand-edit generated code.
Handwritten code owns graph structure, schedules, allocation and backend
orchestration. Keep `.and.` and `.or.` expressions independent of
short-circuit evaluation.

The primary repository is on [fortyard](https://code.tugraz.at/yard); GitHub is
its public mirror. Submit a focused branch to the primary repository, run its
configured checks, and land it there. Changes to the mirror are not a substitute
for landing in the primary repository.
