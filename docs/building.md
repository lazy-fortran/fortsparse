# Building and linking fortsparse

## CMake

CMake is the primary build, with CTest. By default it fetches and builds the
solver libraries it needs (SuperLU and the SuiteSparse/UMFPACK chain) and
[fortnum](https://github.com/lazy-fortran/fortnum) from source, so the only
system requirement is BLAS/LAPACK. OpenBLAS and MKL both work; on a cluster,
load the module first.

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

### Options

All options are set with `-D`.

| Option | Default | Effect |
| --- | --- | --- |
| `FORTSPARSE_ENABLE_SUPERLU` | `ON` | Build the in-process SuperLU backend. `OFF` drops SuperLU; the native backend remains available. |
| `FORTSPARSE_ENABLE_UMFPACK_IPC` | `ON` | Provision the separate GPL UMFPACK helper. `OFF` performs no UMFPACK/SuiteSparse search or fetch and installs no helper; the IPC binding remains for an explicitly supplied helper. |
| `FORTSPARSE_USE_SYSTEM_DEPS` | `OFF` | Find SuperLU and UMFPACK with `find_package` instead of fetching them. |
| `FORTSPARSE_ENABLE_MUMPS` | `OFF` | Build the in-process MUMPS backend; see [backends](backends.md#mumps). |
| `FORTSPARSE_ENABLE_AMD` | `ON` | Use SuiteSparse's BSD AMD for the native AMD ordering; otherwise the built-in implementation is used. |
| `FORTSPARSE_ENABLE_METIS` | `ON` | Use the optional METIS nested-dissection ordering for the native solvers. |
| `FORTSPARSE_ENABLE_NATIVE_CUDA` | `OFF` | Build the native CUDA kernels. |
| `FORTSPARSE_ENABLE_OPENMP` | `ON` | Thread the CSC matrix-vector product with OpenMP. |
| `FORTSPARSE_BUILD_FREEFEM_PLUGIN` | `OFF` | Build the optional FreeFem SDK module ([plugin README](../plugins/freefem/README.md)); requires native-only core linkage and AMD. |
| `FORTSPARSE_BUILD_BENCHMARKS` | `OFF` | Build the benchmark drivers under `bench/`. |
| `FORTSPARSE_INSTALL_HELPER` | `ON` | Install the UMFPACK helper to `bin/`. |
| `FORTSPARSE_TEST_OOM` | `OFF` | Register `test_fortsparse_superlu_status_oom` (label `oom`), which forces a SuperLU allocation failure under `ulimit -v 250000`. Run with `ctest -L oom`; the working limit depends on the toolchain. |

With `FORTSPARSE_ENABLE_UMFPACK_IPC=ON` the helper is built whenever the
SuiteSparse chain is available (fetched, or found with
`FORTSPARSE_USE_SYSTEM_DEPS`).

### BLAS/LAPACK selection

CMake selects BLAS/LAPACK before loading the solver packages. With no parent
project or `BLA_VENDOR` choice it tries OpenBLAS, sequential MKL, then CMake's
normal search (including BLIS); `-DBLA_VENDOR=Generic` selects reference BLAS.
Fetched solvers share this choice. A prebuilt system solver may depend on a
different BLAS: check `ldd` on the final executable and avoid mixing
implementations.

For pthread OpenBLAS together with OpenMP MUMPS, set `OPENBLAS_NUM_THREADS=1`
and control sparse parallelism with `OMP_NUM_THREADS` to avoid nested thread
teams.

### Native-only build

A build without any external sparse solver needs only BLAS/LAPACK:

```sh
cmake -S . -B build-native -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DFORTSPARSE_ENABLE_SUPERLU=OFF -DFORTSPARSE_ENABLE_UMFPACK_IPC=OFF \
  -DFORTSPARSE_ENABLE_MUMPS=OFF -DFORTSPARSE_ENABLE_AMD=OFF \
  -DFORTSPARSE_ENABLE_METIS=OFF
cmake --build build-native --target test_fortsparse_native_only -j 2
ctest --test-dir build-native -R '^test_fortsparse_native_only$' --output-on-failure
```

AMD/METIS ordering and CUDA can be enabled on top of this when available.

## fpm

`fpm.toml` builds the native CPU library with system BLAS/LAPACK and fetches
fortnum. It does not build SuperLU, MUMPS, SuiteSparse AMD, METIS or the UMFPACK
helper (the built-in AMD ordering is included), so
`FORTSPARSE_BACKEND_AUTO` selects the native backend. An explicitly requested
backend that is not configured returns `FORTSPARSE_BACKEND_UNAVAILABLE`.

```sh
fpm build
fpm test
```

The SuperLU test programs report a skip when SuperLU is absent. Standalone C
test programs live under `cmake/tests`, because fpm treats C sources under
`test/` as shared support objects.

## Using fortsparse from a CMake project

```cmake
include(FetchContent)
FetchContent_Declare(fortsparse
  GIT_REPOSITORY https://github.com/lazy-fortran/fortsparse.git
  GIT_TAG main)
FetchContent_MakeAvailable(fortsparse)

target_link_libraries(my_app PRIVATE fortsparse)
```

The UMFPACK backend starts its helper process and looks for the binary next to
the running executable. Either install `my_app` to the standard `bin/`
directory (with `FORTSPARSE_INSTALL_HELPER=ON` the helper is installed there
too), or call

```cmake
fortsparse_colocate_helper(my_app ${CMAKE_INSTALL_PREFIX})
```

which copies the helper next to `my_app` after each build and installs it to the
optional directory argument. It is a no-op when no helper is built. The
environment variable `FORTSPARSE_UMFPACK_HELPER=<path>` names the helper
explicitly when neither layout applies.

## Source layout

- `src/`: library modules (kinds, status, CSC types, backend abstraction and
  solver, SuperLU backend in `src/superlu/`, UMFPACK IPC client in `src/ipc/`,
  native solvers in `src/native/`), with MIT C shims beside their Fortran modules.
- `include/`: C headers for the IPC protocol and the native C API.
- `backends/mumps/`: optional MUMPS binding, outside `src/` so fpm does not
  compile it.
- `helper/`: the standalone GPL UMFPACK helper.
- `test/`, `bench/`, `cmake/`, `docs/`.
