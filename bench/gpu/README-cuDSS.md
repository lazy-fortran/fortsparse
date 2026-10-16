# Isolated cuDSS GPU comparator

`cudss.cmake` adds `bench_spd_cudss` only when
`FORTSPARSE_BENCH_CUDSS_ROOT` names the external CUDA 13 SDK prefix. It links
the shared Fortran fixtures, cuDSS and CUDA, with no FortSparse, OpenBLAS or GPL
backend dependency. NVIDIA SDK files/binaries are not distributed with this
repository. `cudss-sdk.json` pins the downloaded 0.8.0.10 archive and primary
sources; the required SDK notice is in `CUDSS-NOTICE.md` and wrapper comments.
Full archive, header, license, binary and documentation hashes are retained in
`~/code/sparse-dev/sdk/cudss/0.8.0.10-cuda13/provenance.json`, with a local
library card under `library_cards/alternatives/sdk-cudss-0.8.0-cuda13.md`.

For an isolated build, independent of the main library:

```sh
cmake -S test-cmake/cudss_benchmark -B /var/tmp/fortsparse-cudss-bench-build \
  -G Ninja -DCMAKE_BUILD_TYPE=Release -DCUDAToolkit_ROOT=/opt/cuda \
  -DFORTSPARSE_BENCH_CUDSS_ROOT="$HOME/code/sparse-dev/sdk/cudss/0.8.0.10-cuda13/libcudss-linux-x86_64-0.8.0.10_cuda13-archive"
cmake --build /var/tmp/fortsparse-cudss-bench-build -j 2
```

The driver CLI is `FAMILY LEVEL THREADS WORKLOAD VARIANT DEVICE IO_MODE ORDERING
[MATRIX_PATH] [REPEATS] [PERMUTATION_PATH]`. Use `-` for an unused path.
Workloads are `rhs1`, `rhs8`, `warm`, `refactor` and `repeated-rhs`.
`fixture_spd_update` applies a positive, nonuniform diagonal update, retaining
the graph; every repeated RHS uses the shared `fixture_known_rhs`. Each result
is checked against the original full matrix by the independent Fortran
residual and known solution (backward error <=1e-12; relative known-solution
error <=1e-8). Failure/status output is retained.

`host` factor timing includes conversion, runtime initialization, matrix
transfers, analysis and numeric factorization; host solve includes RHS upload
and solution download. `resident` reports device-input/device-output phase
boundaries, with explicit staging/retrieval measured separately. Warm mode
retains the initialized stream, private pool, library handle, configuration and
input-buffer capacity. It destroys/recreates only the factor-data object, redoes
CSR conversion and host uploads, and runs fresh analysis/factorization.
Whole-process time/CPU includes warmup and factor-data release. API-call counters
identify one runtime handle, two data objects and two analysis/factorization
calls. This replaces the earlier recreate-handle warm boundary; historical
measurements retain their original source identity and must not be pooled.
Refactor retains symbolic state.
Repeated solve time is the total over the declared reuse count, not an average.
The wrapper disables hybrid execution, hybrid memory and iterative refinement;
this is an explicitly GPU-only cuDSS configuration, not a claim that its tuning
is optimal. A separate hybrid-execute comparator can be added if required.

Ordering supports `tuned` (vendor default), `amd`, `nd`, and `common`. Common
files contain `n` followed by `n` zero-based new-to-old entries; the driver
checks the permutation, disables reordering and verifies the returned vendor
permutation matches exactly. `ordering_hash` is an opaque integer-sequence
checksum; the runner separately retains the permutation file SHA256. A tiny
registered MatrixMarket test verifies a non-self-inverse supplied permutation
and the independent solve residual. A separate diagonal test requires the exact
vendor `CUDSS_DATA_INFO=2` nonpositive minor; unrelated allocation/device failures
cannot pass it. Both fixtures and permutations are generated into the build
tree by the retained CMake recipe.

The custom allocator uses a private stream-ordered CUDA pool. Its default hard
limit is 256 MiB, configurable via `FSBENCH_CUDSS_POOL_LIMIT_BYTES`; freed memory
is retained until pool destruction and its reservation is measured. Requested,
pool-used and reserved highwater, wrapper buffers, vendor memory estimates and
NVML own-process phase-boundary samples are distinct. CUDA/vendor/context
allocations outside the handler are not capped by the private pool. Total
process GPU peak and vendor pinned-host allocations remain unknown. Sampling
overhead is in `memory_sampling_s` and whole-process costs, outside primary
phase timers. No background poller or performance-counter access is required.
`actual_device` proves context selection; `gpu_offload_verified=0` and unknown
kernel counts remain until an independent trace proves offload.

`run_gpu_spd.py` emits `fortsparse-gpu-benchmark-v1` JSON with flat `trials`,
raw commands/statuses/metrics, randomized interleaved backend blocks, exact
binary/matrix/permutation hashes and declared resources. It bounds each process
to one second, defaults to five independent trials, preserves partial metrics
and timeouts, and labels shared-device results diagnostic. `native`/`native-nd`
request CPU device -1; `native-cuda`/`native-cuda-nd` request the selected GPU.
Unsupported resident/common native paths are explicit unavailable records.
External-only calibration (`--calibration-only`) excludes native backends;
select and freeze larger cases before any native tuning.

Initial shared-device diagnostic correctness used physical GPU1 with
`CUDA_VISIBLE_DEVICES=1`, leaving its LLM service running. First SDK/code-loading
startup hit the unchanged one-second limit for host/refactor; the refactor
already emitted a correct residual before exit. Resident eight-RHS passed in
0.36s. Retrying only the failed registered tests after loading passed both in
0.36s each. These are correctness receipts, not admission timing. The follow-up
retained-runtime implementation passes all six registered GPU checks, including
common ordering, the exact negative-pivot status and the warm lifecycle oracle,
in 0.33–0.38s per invocation under a 128MiB private-pool limit. Five warm trials
each of externally selected grid2 level9 and grid3 level8 also pass the original
one-second bound with residuals below 6e-16. These shared-device diagnostic
results do not establish performance admission.

The retained primary API descriptions distinguish runtime handles and factor
data: [types](https://docs.nvidia.com/cuda/cudss/types.html) and
[data creation/destruction](https://docs.nvidia.com/cuda/cudss/functions.html).
The exact pages and hashes are archived with the SDK provenance. The warm
implementation uses their supported object lifetimes; it retains no symbolic
analysis or numeric factors from the warmup.
