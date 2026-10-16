# Offline GPU offload provenance

Vendor counters may be disabled even when numerical kernels execute on the
GPU. An optional capture receipt can establish observed CUDA activity without
changing `gpu_offload_verified`, `actual_device`, execution-mode labels, trial
status, timing or memory. Diagnostic trace events are never timing samples.
This proves activity in the captured workflow, not that every timed phase is
on the GPU or that offload is fast. Hybrid factor/solve paths remain hybrid.

After capturing the **exact final executable and exact case**, export its
Nsight Systems report to SQLite. Retain the executable before rebuilding it,
the SQLite export, the original report and capture command/environment. The
receipt itself records executable SHA256, full benchmark argv, configuration,
input identities, environment, trace PID/hash, hardware mapping and documented
numerical kernel names/counts. SQLite does not attest the binary hash, argv,
pool setting or input checksums: those are capture assertions. The validator
matches them to ordinary-run provenance; it independently queries the process,
CUDA device mapping and kernel events. A later rebuild needs new evidence.

Use this JSON layout (`proofs` may contain multiple exact captures). Paths are
relative to the manifest, or absolute. Replace all illustrative values; no
wildcards or case/ordering/thread/pool aliases are expanded.

```json
{
  "schema": "fortsparse-gpu-offload-proof-v1",
  "proofs": [{
    "id": "cudss-grid3-8-warm-t1-p128",
    "executable": {"path": "retained/bench_spd_cudss", "sha256": "EXACT_BINARY_SHA256"},
    "capture": {
      "command": ["taskset", "--cpu-list", "16", "/original/bin/bench_spd_cudss",
        "grid3", "8", "1", "warm", "training", "0", "host", "tuned", "-", "1", "-"],
      "configuration": {
        "requested_backend": "cudss", "family": "grid3", "level": 8,
        "workload": "warm", "variant": "training", "io_mode": "host", "threads": 1,
        "ordering": "tuned", "repeats": 1, "pool_budget_bytes": 134217728,
        "requested_device": 0
      },
      "environment": {"CUDA_VISIBLE_DEVICES": "1", "FSBENCH_CUDSS_POOL_LIMIT_BYTES": "134217728"},
      "inputs": {"matrix_sha256": null, "permutation_sha256": null},
      "logical_device": 0,
      "physical_gpu": {"index": 1, "uuid": "GPU-EXACT_UUID", "pci_bus_id": "0000:08:00.0"}
    },
    "trace": {"path": "capture.sqlite", "sha256": "EXACT_SQLITE_SHA256", "pid": 12345},
    "kernels": [{"name": "EXACT_DEMANGLED_NUMERICAL_KERNEL_NAME", "count": 123,
      "purpose": "Numerical factorization update, identified from the retained kernel timeline"}]
  }]
}
```

The trace process path must equal `command[3]`, which can differ from the
archived executable path. Its hash must match the ordinary run's
`binary_provenance` entry. Configuration and full command must match exactly.
Record the relevant pool variable for cuDSS/CHOLMOD, including a JSON null if
unset, and `CUDA_VISIBLE_DEVICES`. Every recorded environment field must match
the ordinary-run environment. Matrix/permutation SHA values must match their
ordinary input provenance. Other unrelated environment variables need not be
recorded. Logical device IDs are resolved through `TARGET_INFO_CUDA_DEVICE`
to `TARGET_INFO_GPU` UUID, index and PCI bus; kernel events are restricted to
that logical device and the captured process, excluding other services.

The current supported export uses `PROCESSES`, `TARGET_INFO_CUDA_DEVICE`,
`TARGET_INFO_GPU`, `CUPTI_ACTIVITY_KIND_KERNEL` and `StringIds`. Each documented
kernel name must appear with exactly its recorded positive event count and
positive duration. Unsupported schemas, missing events, wrong mappings or
changed file hashes fail closed. Reusing a valid proof for another case,
binary, input, thread count, ordering, pool or visible-device setting does not
attach. Ambiguous matching receipts are rejected.

Revalidate evidence when producing an audit:

```sh
python3 bench/gpu/analyze_gpu_spd.py ordinary.json \
  --offload-provenance proofs.json --output audit.json
```

Or produce an annotated copy for plotting without changing the original:

```sh
python3 bench/gpu/offload_provenance.py ordinary.json \
  --offload-provenance proofs.json --output ordinary-with-trace-proof.json
```

Matched trials gain `trace_verified`, `trace_actual_device` and `offload_proof`
with manifest/trace identities and observed counts. Analyzer verification
accepts either the original backend metric or a freshly validated trace proof.
An annotated boolean alone is insufficient: pass the manifest again for a
later audit. Executable hashes separate summary series; instrumentation data
never creates additional ordinary trials. Unknown memory and all existing
accuracy/completeness/frontier gates remain unchanged.
