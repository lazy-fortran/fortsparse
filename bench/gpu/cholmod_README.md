# Standalone CHOLMOD CUDA comparator

This optional executable links the GPL-2.0+ CHOLMOD Supernodal/GPU modules.
It is built separately and never links into `fortsparse`. The fixture and
independent residual oracle are the same Fortran module used by other SPD
comparators. Build only against a CUDA-enabled private SuiteSparse prefix:

```sh
cmake -S bench/gpu/cholmod_standalone -B /var/tmp/fortsparse-gpu/cholmod-driver \
  -DFORTSPARSE_BENCH_CHOLMOD_ROOT=$HOME/code/sparse-dev/gpu-install/cholmod
cmake --build /var/tmp/fortsparse-gpu/cholmod-driver --target bench_spd_cholmod -j 4
```

The owning root may instead include `bench/gpu/cholmod.cmake` after creating
`fortsparse_spd_fixtures`. Enable only the SuiteSparse config, AMD, CAMD,
CCOLAMD, COLAMD and CHOLMOD projects in the dependency build; use
`SUITESPARSE_USE_CUDA=ON`, `CHOLMOD_USE_CUDA=ON`,
`SUITESPARSE_CUDA_ARCHITECTURES=120`, shared libraries and OpenBLAS.
Record exact source, compiler, flags, configuration and binary hashes.
`SUITESPARSE_USE_STRICT=ON` prevents an unnoticed CUDA-disabled build.

```
bench_spd_cholmod FAMILY LEVEL THREADS WORKLOAD VARIANT DEVICE IO_MODE ORDERING \
  [MATRIX_PATH] [REPEATS] [PERMUTATION_PATH]
```

Workloads are `rhs1`, `rhs8`, `warm`, `refactor`, `repeated-rhs`.
`VARIANT` is `training` or `heldout`; an unused matrix path is `-`.
`IO_MODE=host` is supported. `resident` returns unavailable (exit 77):
CHOLMOD factors remain on the host and its solves execute on the CPU.
FP64 supernodal factorization may offload dense operations to the selected GPU;
this is explicitly a hybrid comparator, not an entirely device-resident solve.
The adapter uses `cholmod_l_*`: CHOLMOD's GPU code requires 64-bit indices.

Ordering is `tuned` (CHOLMOD defaults), `amd`, `nd` (METIS), or `common`.
A common permutation file contains `n`, then `n` zero-based new-to-old indices.
The driver validates the bijection, disables postordering, and checks the actual
factor permutation. Set `CHOLMOD_PERMUTATION_OUT` to export its actual order.
`ordering_hash` is FNV-1a over little-endian int64 permutation entries; the
runner should also retain the permutation file's SHA256.

Set `CHOLMOD_GPU_POOL_BYTES` explicitly when calibrating memory; the adapter
default is 268435456 (256 MiB). Values below 64 MiB are rejected. CHOLMOD's
unbounded default otherwise requests nearly all available GPU memory. Retain
pool sweeps on external solvers before native tuning, including failures.
Report `cholmod_device_pool_peak_bytes` and `pinned_host_bytes` separately from
CHOLMOD's tracked host allocator peak, wrapper arrays and fixture bytes.
These attributed sizes do not establish a complete process/device peak: driver,
cuBLAS and allocator internals may be unreported. Unknown allocator/process
GPU peaks are `-1`; process RSS is measured by the common process-metrics helper.
The refactor fixture temporary's full matrix allocation is reported separately.

`factor_s` includes cold runtime setup, CSC conversion, analysis and synchronized
numeric factorization. `cold_factor_s` preserves that first cost. For `warm`,
`factor_s` measures the subsequent analysis/factorization with the CUDA runtime
and CHOLMOD pools already loaded; whole-process totals include both runs.
`solve_s` includes CHOLMOD solve and copying the host result. For repeated RHS,
it is the total over all declared repeats; fixture RHS construction is included
in total workflow time. Refactor uses the shared positive nonuniform diagonal
update and reconstructs the exact RHS. Every result is independently checked
against the original full matrix at backward error below 1e-12.

GPU operation counters are accumulated over all numerical factorizations.
The pinned CHOLMOD source updates these counters only when compiled with
`BLAS_TIMER`, which its ordinary build disables. Positive counts supply evidence;
zero counts alone do **not** establish CPU fallback. In that situation the driver
reports `actual_device=-1` (numerical execution unknown),
`execution_mode=offload-unknown`, `gpu_operation_calls=-1` and
`gpu_offload_verified=-1`. Preserve the raw individual counter values and attach
an external Nsight Systems evidence record keyed by exact binary hash, case,
ordering, pool and device. A kernel trace of the ordinary binary establishes
actual offload without instrumenting/rebuilding CHOLMOD. Do not use instrumented
or shared-service diagnostic timings as ordinary performance samples.

The internal update threshold requires at least 256 descendant rows and 32
columns; GPU POTRF/TRSM requires at least 512 columns. POTRF additionally
requires the supernode to have received GPU work. These are necessary shape
conditions, not a guarantee of offload. In a simple geometric nested-dissection
model, a 3D root separator has side-squared size (512 columns needs side at
least 23), while a 2D line separator reaches 512 at side 512. Actual graph
ordering, supernode splitting and descendant-tail sizes must still be checked.

Internal transfer and synchronization phase times are unknown (`-1`) because
the public API does not separate them from numeric factorization. Use profiling
for that breakdown. `gpu_operation_calls` counts CHOLMOD GPU BLAS operations
when available, not the profiler's number of CUDA kernels.
