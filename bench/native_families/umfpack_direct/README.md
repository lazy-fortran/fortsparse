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

`umfpack_*_defaults` are retained without overrides. Numeric refactor frees the
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
