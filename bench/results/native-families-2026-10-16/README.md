# Native sparse-family performance checkpoints, 2026-10-16

These numeric receipts retain all declared configurations and five-trial ranges.
They are diagnostic performance evidence, with no frontier admission yet.

- `initial/`: source 42b8aca, 2,500 attempts, 2,350 completed and 150 timeouts.
  The old `parent_cpu_s` is a main-entry interval; complete CPU is unavailable.
- `amd/`: source 52dde03, 1,000 completed attempts, 200 matched groups, both
  16,384-unknown drift fixtures, training/heldout, normal operation, cold RHS8
  and same-pattern refactor, five requested thread budgets and five backends.
  Workflow CPU uses completed `wait4` user+system usage, including reaped helpers.

Each directory retains the report, numeric CSV, recipe, source/binary provenance
and paths/hashes for durable raw attempts and external-only calibration. Raw
receipts and rendered PNGs remain outside Git. The initial checkpoint's original
report and CPU interpretation are retained unchanged.

The AMD checkpoint improves completion relative to the earlier natural-order
snapshot, but it still trails the fastest declared external time on these cases.
The snapshots have different recorded host loads and CPU instrumentation; they
are not matched before/after trials. Native factor teams remain serial at every
requested budget. External active teams are unknown, rather than assumed.

The AMD parent RSS measurement is valid process VmHWM; attributed peak bytes
are unavailable and must never be treated as zero. UMFPACK currently uses the
helper/IPC comparator. Its parent and helper lifecycle RSS cannot establish
aggregate peak RSS. The time/RSS figures omit UMFPACK aggregate points explicitly.
A standalone direct UMFPACK comparison is required before performance admission.

Regenerate the AMD plots outside source history:

```sh
python3 bench/native_families/report.py \
  /home/ert/code/sparse-dev/benchmarks/2026-10-16-native-families/amd-measurement \
  --output /tmp/fortsparse-amd-report --plots --cpu-model completed \
  --cpu-scope 'completed wait4 workflow CPU including reaped descendants'
```

`summary.csv` keeps successful counts, minimum/median/maximum values and timeout
counts separately. Empty numeric summaries mean unavailable; metric count zero
must not be converted to a zero-valued measurement. Timed factor, analysis,
numeric, solve, refactor, prime-factor and main-workflow phases are distinct from
complete subprocess wall. Parent cumulative pre-CSV CPU, main-entry interval,
helper CPU and completed workflow CPU have separate columns.

`indexed-ab/` retains the subsequent interleaved scalar/generated-indexed native
LU diagnostic: 80 completed trials, five pairs per workload/variant, one actual
factor thread and the same AMD ordering/fill. It adds no external-frontier claim
and preserves the original four-external size-selection scope. Its standalone
recipe records both frozen source/binary identities and paired ratio ranges.
