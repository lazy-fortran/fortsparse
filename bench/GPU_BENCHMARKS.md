# Matched single-node GPU SPD benchmarks

This is the GPU admission contract. Existing CPU evidence remains a separate
corpus. Primary comparators are native CUDA, NVIDIA cuDSS and standalone
CHOLMOD CUDA. PaStiX Cholesky CUDA with MPI disabled is a later comparator when
its runtime is available. Record exact versions, licenses, build options and
binary identities. Vendor/GPL adapters do not become core dependencies.

## Inputs and resources

Use identical full canonical FP64 SPD matrices and known-solution RHS, with
one and eight RHS. Evaluate backward error against the original matrix using
an independent CPU implementation; the acceptance threshold is 1e-12. Check
known solutions, invalid inputs and nonpositive pivots. Refactorization changes
values nonuniformly while preserving SPD and the sparsity pattern.

Use one declared GPU with fixed physical host-core affinity and thread budget.
Record driver, device, available VRAM, background load and actual host teams.
Shared-device timings are diagnostic until a controlled measurement window is
available. Prove offload through library counters and a kernel trace; successful
CUDA configuration alone does not establish GPU execution.

For CPU/GPU comparisons, pass `--cpu-bin-dir` pointing to the CPU-only native
build. The default runs the CPU path of the CUDA-linked executable, whose host
RSS includes CUDA shared-library mappings. Retain that distinction for older
records; it does not establish standalone CPU memory use.

Retain the small latency corpus. Select larger, higher-fill and public SPD
cases using external GPU implementations alone, freeze that calibration, then
tune native execution. Ordinary invocations must finish below one second;
five independent trials are the default. Native failures never shrink the
corpus. Keep held-out coefficients and randomized, block-interleaved solver
order. Preserve every failure, timeout and unavailable comparison.

The GPU calibration ladder extends the existing levels without changing them:
2D sides 192/256/384/512 and 3D sides 28/32/40/48 are levels 8 through 11.
This candidate list precedes native GPU tuning. An external-library timeout or
absent numerical offload is retained, not converted into a GPU success. The
generic `matrix-market` fixture reads separately hashed public SPD inputs;
the existing `bcsstk09` fixture keeps its frozen dimension checks.

## Workflow boundaries

| Workflow | Retained state and measured boundary |
|----------|--------------------------------------|
| Cold | Fresh process, host matrix/RHS to host solution, including runtime initialization, conversion, analysis, transfers, factorization and solve. |
| Warm | Initialized runtime, fresh analysis/factors, host matrix/RHS to host solution. Record warmup and whole-process costs separately. |
| Refactor | Retained symbolic analysis, changed host values and RHS to host solution. |
| Repeated RHS | Retained numeric factors, host RHS to host solution. |
| Resident | Device matrix/RHS to device solution, with setup and ownership stated separately. |

Native currently downloads factors and solves on the CPU. It is eligible for
the host workflows only. CHOLMOD's GPU factorization likewise does not imply
a device-resident solve. Do not impose factor downloads on a competitor to
equalize implementation limitations. Mark unsupported resident paths unavailable.

Measure each solver's best supported ordering and declared tuning separately
from a common, retained AMD/ND permutation with further reordering disabled
where possible. The same permutation does not imply identical supernodes.
Initially disable mixed precision, compression and refinement that changes
the accuracy contract; later accuracy-matched cohorts must be separate.

## Timing, memory and admission

Synchronize before host wall-time boundaries. Device events supplement host
timing; launch time alone is insufficient. Report phase times, transfers,
initialization, synchronization, checks and total workflow. Report process CPU
user plus system time across all threads, and whole-process wall time.

Report host peak RSS and pinned storage separately from device memory. Device
accounting distinguishes live requested allocations, reserved pools, attributed
buffers and observed process/device high-water marks including vendor/context
costs. Custom allocator counters do not prove total device memory. Sampling
can miss short peaks: retain its interval and limitations. Unknown aggregate
memory blocks a claim requiring that quantity. Never sum RSS and VRAM into an
unqualified memory number. Declare and externally calibrate bounded library
pool settings, including CHOLMOD's paired device/pinned pools.

Retain raw trials and ranges. Audit time/VRAM, time/host-RSS and the joint
time/host-RSS/VRAM frontier; report host CPU cost alongside them. Robust
observed-range membership must exclude every measured external configuration
as a possible dominator. Preserve the stronger same-comparator match/beat gate
separately. Median-front membership is weaker evidence. No CPU-only comparator
or inaccessible library establishes GPU leadership.

Profile unresolved losses with Nsight Systems and supported counters. Record
instrumentation overhead and denied counters. Inspect compiled generated code,
not just symbolic operation counts. Repair reusable emission defects in the
owning generator and regenerate numerical leaves; keep sparse topology and
hardware schedules in FortSparse. Recheck CPU accuracy/performance after shared
generator changes. Publish plots through slopbox, retaining recipes and compact
numeric/provenance evidence in Git, never generated images or large traces.

GPU offload evidence may also come from an exact-binary, exact-configuration
Nsight capture, validated offline against retained SQLite process/device/kernel
events. Keep trace verification separate from vendor counter fields and retain
the capture assertions and their limits. See
[`gpu/OFFLOAD_PROVENANCE.md`](gpu/OFFLOAD_PROVENANCE.md) for the receipt schema
and analyzer command. Diagnostic trace timings never join ordinary trials.
