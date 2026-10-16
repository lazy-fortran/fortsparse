# Optional parent-observed device-memory high-water counter

Add `--accounting-peak` to the ordinary GPU runner to read the exited benchmark
child's NVML accounting record. This is optional monitoring code, outside the
scientific library/runtime. The harness does **not** enable, disable or clear
GPU accounting. The controller must establish a bounded accounting window and
restore its prior state outside this runner. Existing GPU services stay running.

The parent loads NVML and discovers already-enabled device handles once,
before cases begin. After `Popen.communicate()` returns, including the
kill/wait path for a timeout, it queries only that child PID. It never scans
another process for a substitute measurement, creates a CUDA context or
sums records across GPUs. Initialization/shutdown costs are in the report's
`accounting_observer` metadata. Per-query enclosing parent wall/process CPU
costs are `accounting_observer_wall_s` / `accounting_observer_cpu_s`; these are
excluded from `process_wall_s` and the child's workflow CPU metrics.

For one uniquely matched, terminated record with an available counter:

- `device_total_peak_bytes` receives `nvmlAccountingStats_t.maxMemoryUsage`;
- `device_total_peak_source` is `nvml_accounting`;
- `vendor_device_total_peak_bytes` preserves the adapter's original value;
- `accounting_pid`, GPU UUID, `accounting_is_running`, start timestamp and
  context time identify the observed driver record.

Missing/disabled accounting, unsupported APIs, driver errors, unavailable
`ULL_MAX` values, still-running records and multiple matching GPUs yield
`device_total_peak_bytes=-1`, source `unknown`, and an explicit reason. Phase
samples and wrapper allocation counts never replace this unknown value.
Benchmark failures/timeouts keep their original status. An unavailable binary
launches no process and triggers no per-case accounting query.

The quantity is NVIDIA's driver-reported process-lifetime memory maximum on
one GPU, including the ownership visible to the driver. It is not an exact
sum of wrapper allocations, a phase-specific maximum, an independent allocation
audit or total memory across GPUs. Driver accounting is distinct from sampled
process memory, requested/reserved allocator pools and host RSS. The driver
keeps a circular PID buffer and reports the latest terminated process on PID
collision; the harness queries its own child immediately after exit and
retains the reported start timestamp. It does not infer a timestamp epoch.

The ctypes binding follows the local CUDA 13.4 `nvml.h` declaration and the
[official NVML accounting structure documentation](https://docs.nvidia.com/deploy/nvml-api/latest/api/structnvmlAccountingStats__t.html).
The legacy unversioned `nvmlDeviceGetAccountingStats(device,pid,stats)` uses
56 bytes on the tested 64-bit C ABI: two unsigned ints, three unsigned long
longs, an unsigned int and five reserved unsigned ints. The newer `_v2`
structure has a different layout and is never cast to this one. Legacy init,
count and device-handle symbols are accepted when their `_v2` aliases are
absent; absent legacy stats support remains unknown.

The registered `gpu_spd_accounting_peak` test uses independent packed C-layout
bytes and mock NVML functions. It checks 64-bit values, legacy symbols,
disabled/unavailable/error/sentinel cases, other PIDs, multiple GPUs and
runner success/failure/timeout paths without loading NVML or changing devices.
Actual-driver validation belongs to the controller's explicitly granted probe
window; mock results alone do not establish hardware support.
