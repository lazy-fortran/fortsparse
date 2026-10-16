# Standalone UMFPACK comparator

This optional benchmark-only executable directly links UMFPACK. Both driver
sources in this directory are **GPL-2.0-or-later**, with the license in COPYING.
It never links FortSparse. The shared fixtures, record writer and process metrics
remain permissively licensed. Neither the product IPC protocol nor the MIT
library's linkage/license changes. The benchmark uses general sparse LU even
when a fixture is Hermitian positive definite or symmetric indefinite.

The installed official `/usr/include/suitesparse/umfpack.h` inspected here is
UMFPACK 6.3.8 and identifies UMFPACK itself as GPL-2.0+. Its hash is retained in
`../contract.json`. The shim includes that header rather than reproducing vendor
prototypes. It uses 32-bit `umfpack_di_*`/`umfpack_zi_*`, matching native I4
fixtures; C static assertions check 32-bit int. Fortran converts CSC indices to
explicit `c_int32_t`, zero-based arrays. Real data uses `c_double`; complex data
uses `c_double_complex`, compatible with UMFPACK's documented packed layout.
Az, Bz and Xz are NULL for packed complex calls. Scalar-kind/layout guards make
the assumed interoperability explicit.

Header evidence: packed complex syntax in symbolic/numeric/solve documentation;
solve system table (`UMFPACK_A` normal, `UMFPACK_Aat` ordinary transpose,
`UMFPACK_At` conjugate transpose); memory Info documentation
(`UMFPACK_SIZE_OF_UNIT`, `UMFPACK_PEAK_MEMORY`); license/version header.
Real T and H both select `UMFPACK_At`. The public I4 prototypes are checked by
the C compiler against the official header. The product helper uses I8 APIs;
this comparator deliberately measures the I4 counterpart to native factors.

`umfpack-direct` retains `umfpack_*_defaults` without overrides (IRSTEP=2 in
the inspected version). `umfpack-direct-noir` changes only UMFPACK_IRSTEP to
zero, matching the IPC/native no-refinement solve workflow. Numeric refactor frees the
old Numeric object and retains Symbolic. RHS8 is eight retained-factor solves.
Warm workflow destroys the priming context before a fresh symbolic/numeric
factorization; total wall/CPU includes priming, conversion, checking and cleanup.
The usual known solution and original full-CSC residual verify N/T/H independently.

Vendor factor fill is LNZ+UNZ-n, with both vendor diagonals counted once.
Vendor peak memory is Info[PEAK_MEMORY]*Info[SIZE_OF_UNIT]/1e6 decimal MB,
cached after numeric factorization because solve may overwrite Info. These
numbers are distinct from process memory and exact owned payload remains
unavailable. The complete standalone process's Linux VmHWM is measured after
cleanup and resets at exec; helper fields are zero. The registered memory test
forces fork+exec from a large live launcher to reject inherited-launcher HWM.
Cumulative process CPU is sampled with the shared process-metrics implementation.

The system library used for initial validation resolves `/usr/lib/libumfpack.so.6`
and reference `/usr/lib/libblas.so.3.12.0`, with LAPACK. No FortSparse, MUMPS or
MKL objects/libraries are linked. Thus this is a usable isolated comparator,
not a best-tuned performance claim. The runner retains binary hashes and ldd;
future builds against a different BLAS are different explicitly recorded
comparison checkpoints. Requested threads are recorded, while actual active
factor teams remain unknown; the runner fixes OpenBLAS to one BLAS thread.

Build/test (without building the core library or changing installed libraries):

```sh
cmake -S . -B build-umf-direct -G Ninja \
  -DFORTSPARSE_USE_SYSTEM_DEPS=ON -DFORTSPARSE_ENABLE_SUPERLU=OFF \
  -DFORTSPARSE_ENABLE_UMFPACK_IPC=OFF -DFORTSPARSE_BUILD_BENCHMARKS=ON
cmake --build build-umf-direct --target bench_native_families_umfpack -j 2
ctest --test-dir build-umf-direct -j 2 \
  -R '^native_family_(umfpack_direct_|selector$)' --output-on-failure
```

The driver takes case, level, requested threads, workload, N/T/H, variant and
mini5 path, like the isolated PARDISO driver. The runner selects it under the
separate `umfpack-direct` label. No timing/profile/admission campaign accompanies
this implementation; smoke elapsed times are test execution evidence only.

The optional local provision recipe uses committed source from the retained
SuiteSparse clone, not its uncommitted research edits. It builds only UMFPACK
and its required SuiteSparse_config/AMD/CAMD/CCOLAMD/COLAMD/CHOLMOD dependencies,
with shared libraries, Release, LP64 OpenBLAS, no CUDA, no demos, no tests, no
static libraries and two build processes. GraphBLAS and other projects are not
selected. All source, build, install and provenance output stays in TASK_ROOT;
no downloads or system installation occur.

```sh
bash bench/native_families/umfpack_direct/provision_openblas.sh \
  /home/ert/code/sparse-dev/repos/SuiteSparse \
  /home/ert/code/sparse-dev/provisions/umfpack-openblas-34db4df
```

The pinned commit is `34db4dfbe21c53cea7b0e3c5ee2da37bfaef369a`; its UMFPACK
version is 6.3.8, matching the system header used above. OpenBLAS 0.3.34 uses
`/usr/lib/libopenblas.so`, not the separate ILP64 `libopenblas_64.so`. The recipe
records exact source archive/tree hashes, header/license/library hashes,
OpenBLAS runtime configuration, CMake options and dependency resolution in
`TASK_ROOT/provenance.json`. It rejects a runtime reporting USE64BITINT.
The LP64 compile/runtime probe is included in the build receipt; consumer
correctness checks have a separate receipt. No timing or performance conclusion follows from provisioning.
The original reference-BLAS executable and system libraries remain untouched.

For a fresh FortSparse build tree, the project-owned FindUMFPACK module accepts
these explicit cache entries (setting UMFPACK_DIR alone does not select it):

```sh
prefix=/home/ert/code/sparse-dev/provisions/umfpack-openblas-34db4df/prefix
cmake -S . -B build-umf-openblas -G Ninja \
  -DFORTSPARSE_USE_SYSTEM_DEPS=ON -DFORTSPARSE_ENABLE_SUPERLU=OFF \
  -DFORTSPARSE_ENABLE_UMFPACK_IPC=ON -DFORTSPARSE_BUILD_BENCHMARKS=ON \
  -DCMAKE_PREFIX_PATH="$prefix" \
  -DUMFPACK_INCLUDE_DIR="$prefix/include/suitesparse" \
  -DUMFPACK_LIBRARY="$prefix/lib/libumfpack.so" \
  -DSUITESPARSECONFIG_LIBRARY="$prefix/lib/libsuitesparseconfig.so"
```

UMFPACK notice from the retained source: UMFPACK, Copyright (c) 2005-2023,
Timothy A. Davis, All Rights Reserved. GPL-2.0-or-later; original version
available from <http://www.suitesparse.com>. Used by permission. The installed
header carries the updated 2005-2026 copyright. The exact source license and
all dependency license files remain in the retained source archive; OpenBLAS's
BSD license is `/usr/share/licenses/openblas/LICENSE` on the provisioned host.

Local provision validation on 2026-10-16 passed 26 focused registered tests:
22 standalone N/T/H/workflow/record checks, the selector check, and IPC
roundtrip, pivot-control and retained-operation checks. A compile/runtime
probe checked `sizeof(blasint)==4`, `sizeof(SUITESPARSE_BLAS_INT)==4`, runtime
UMFPACK/header version agreement and one exact 1x1 DGEMM. The standalone driver
and helper resolve only the local SuiteSparse closure plus OpenBLAS; no
reference BLAS/LAPACK, MUMPS, MKL or FortSparse library appears in their `ldd`.
The MIT archive has no unresolved UMFPACK/CHOLMOD symbols and the IPC consumer
has no direct GPL library dependency. Debug checks used `-fcheck=all`.

Exact hashes, link maps, CMake flags and check-log hashes are retained beside
the prefix in `provenance.json` and `consumer-validation.json`. The fresh
consumer build is `/var/tmp/fortsparse-umfpack-openblas-consumer`; these are local
smoke checks, not benchmark performance evidence. The recipe-generated C ABI
probe is GPL-2.0-or-later, like the standalone comparator.


The optional eighth CLI argument selects `default` or `no-refinement`; omission
keeps the original default mode. The latter emits the distinct
`umfpack-direct-noir` CSV label. Both initial and warm contexts receive the
selected control. Refactor retains it. The final `umfpack_controls,` JSON line
records the runtime version, IRSTEP index, independent original-default vector
and effective vector used by every vendor call. These bounded diagnostics
remain in raw receipts and are validated by the runners. The vendor routines
accept Control as const; the snapshot is taken after execution. Accuracy gates
and general-LU classification remain unchanged. Both configurations enter the
new six-comparator coverage contract; tuning/admission remains open pending
calibration and measurements. Historical default results are preserved.
