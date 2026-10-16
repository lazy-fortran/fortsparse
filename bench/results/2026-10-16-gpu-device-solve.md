# GPU device-solve checkpoint

This is a usability and performance-visibility checkpoint under a shared GPU,
not a controlled performance admission or a GPU frontier claim. Of 320
fresh-process trials, 319 passed the unchanged accuracy/deadline checks and one
native CUDA ND large-2D warm invocation timed out at 1 s. The timeout remains
in the raw evidence and compact CSV. The maximum backward error among successes
was `1.26e-15`, below the `1e-12` gate. The nine affected actual-device CUDA
tests pass, including the public permutation/block/JVP/VJP contract; each ran
below its unchanged 1 s limit. Its complex-GPU negative test was corrected to
expect the explicit unsupported-backend status.

Native CUDA now retains numeric factors on device and uses generated triangular
solves for scalar and batched RHS. The integrated public block path tiles eight
RHS columns. The checkpoint measures one RHS: warm new analysis/factor+solve and
changed-value refactor+solve with retained analysis. Inputs/answers remain on
the host. The prior DAG reference retains GPU factorization and CPU solves.
Derivative and nine-RHS correctness are covered by the focused tests rather
than timed in this checkpoint.

Median phase wall times are in milliseconds. Each complete configuration has
five trials; `*` is the median of four successes plus one retained timeout.
Observed minima/maxima and every configuration's status count are in the CSV.

| Case | Workflow | CPU auto | CPU ND | GPU auto | GPU ND | cuDSS | CHOLMOD |
|---|---|---:|---:|---:|---:|---:|---:|
| 2D 256×256 | warm factor+solve | 84.94 | 236.16 | 129.63 | 286.43* | 179.98 | 78.39 |
| 2D 256×256 | refactor+solve | 32.32 | 33.49 | 55.57 | 58.83 | 7.80 | 61.13 |
| 3D 28×28×28 | warm factor+solve | 150.36 | 184.68 | 137.72 | 214.73 | 175.02 | 153.31 |
| 3D 28×28×28 | refactor+solve | 99.17 | 67.13 | 71.14 | 72.02 | 56.69 | 59.43 |

The four frozen cases are `grid2:3` (32×32), `grid3:7` (24³),
`grid2:9` (256×256), and `grid3:8` (28³), with training coefficients.
Sizes and comparator pools were frozen by earlier external-only calibration;
no failing case was shrunk. Every configuration has one allowed physical host
core and physical GPU 1 is exposed as logical device 0. The running model
service remains active. Controller/worker compiler and profiler jobs were
paused during the ordinary trials; unrelated host activity remains. The
observed load averages changed from `3.67/6.97/9.05` to `2.79/5.29/8.10`.
No scaling was sampled.

Interleaved, SHA-separated comparisons with the immutable prior DAG binary
show a mixed result. For large 2D, median solve time improves from 7.25 to
5.06 ms with default ordering and from 7.71 to 3.15 ms with ND. Refactor+solve
changes from 59.47 to 55.57 ms and from 65.02 to 58.83 ms respectively.
For large 3D, solve time increases from 4.89 to 8.41 ms with default ordering
and from 4.70 to 5.28 ms with ND. Refactor+solve is 70.86→71.14 ms and
73.91→72.02 ms respectively. Small-case GPU solves also regress. These are
descriptive medians and observed ranges, not confidence intervals or a claim
of an all-case speedup.

CPU and memory costs remain material. On large-3D refactor, default native GPU
uses 424.56 ms whole-workflow process CPU and 603.71 MiB peak host RSS; native
CPU ND uses 255.95 ms and 47.04 MiB. cuDSS uses 379.67 ms/553.18 MiB and
CHOLMOD 508.92 ms/719.55 MiB. Whole-workflow CPU includes initial factorization,
runtime warmup, setup, residual checks and cleanup; it is not phase CPU time.
Host RSS, attributed device allocations and phase-sampled device usage are
separate. Total device peak is unknown; sampled usage does not establish a
device-memory frontier. No NVML accounting modes or services were changed.

Five separate Nsight Systems invocations verify exact executable/configuration
GPU execution for cuDSS and CHOLMOD on both large refactor cases, plus native
default large-3D refactor. The strict validator checks retained SQLite hashes,
process paths/PIDs, GPU UUID/bus mapping and positive documented compute-kernel
counts. These proofs attach to 25 ordinary trials without changing vendor
flags, trial status or timing. Other cases/workflows remain untraced. The native
trace contains 35 forward and 35 backward generated solve kernels. Vendor
proofs establish offload during the complete invocation, not attribution to
every timed phase. Profiler runs took about 3.6–3.8 s including startup/export
and are excluded from ordinary timing. An earlier launcher-inside-profiler
capture was rejected by the strict process-path check; those raw receipts are
retained separately, and the successful captures put CPU affinity outside the
profiler. The validator was not weakened.

Immutable executable SHA-256 prefixes are native GPU `6f88affe`, native CPU
`ac9a6872`, cuDSS `39e0a811`, CHOLMOD `ae9d7df2`, and prior native DAG
`db885359`. Native numerical source is main `a3bd4a2`, with GPU execution label
`c6f3a4f`; later `ea99799` only corrects the negative test expectation. Runner
source/cache hashes describe the observed checkout at capture/completion;
the immutable binary SHA is authoritative. cuDSS 0.8.0.10/CUDA 13 and the
standalone CHOLMOD adapter remain isolated from the core vendor linkage.
Pools are 128 MiB for cuDSS and 64 MiB for CHOLMOD. Full identities, commands,
compiler/cache information and GPU observations are retained in the manifest
and durable raw archive.

Reproduce compact data and the three inspected PNGs with:

```sh
python3 bench/gpu/device_solve_checkpoint.py \
  /home/ert/code/sparse-dev/benchmarks/2026-10-16-gpu-mvp/device-solve-checkpoint.json \
  --offload-provenance /home/ert/code/sparse-dev/benchmarks/2026-10-16-gpu-mvp/device-solve-profiles/offload-provenance.json \
  --output-dir /var/tmp/fortsparse-gpu/device-solve-plots
```

The optional plotting recipe uses the existing benchmark summary contract and
Matplotlib. PNGs and full raw traces stay outside source history. The compact
CSV, exact offload receipt and manifest accompany this report. The ordinary
runner command/environment are retained verbatim in the manifest; each child
keeps its 1 s deadline and all five trials.
