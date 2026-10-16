# IPC helper post-exec memory diagnostic

Setting `FORTSPARSE_BENCH_HELPER_METRICS=1` requests one helper stderr line on
normal shutdown: `FORTSPARSE_BENCH_HELPER_VMHWM_KIB=<positive integer>`.
Linux reads `/proc/self/status` after releasing factors; VmHWM retains the
helper's post-exec resident high-water. Other platforms or a failed read emit
`=unavailable`. Default helper behavior emits no diagnostic. Neither the shared
memory protocol nor solver API changes.

`helper_memory.parse_helper_memory(stderr)` validates the receipt. Missing
historical receipts and explicit unavailable values retain `-1`. Malformed,
nonpositive, overflowing or duplicate receipts raise `ValueError`. Preserve the
raw stderr, requested environment and status with each new measurement.

This is one helper's process peak and a lower bound on the parent-plus-helper
union peak. Parent and helper peaks can occur at different times and can include
shared mappings: neither their sum nor their maximum establishes aggregate
peak RSS. Unknown aggregate memory remains unknown. A strictly larger valid
helper lower bound can exclude a candidate comparator as a possible memory
dominator, but cannot establish IPC domination or an aggregate-memory win.

`spd_helper_postexec_memory` checks the normal helper READY/SHUTDOWN path,
receipt validation and default silence. Its Linux fork/exec witness temporarily
holds 128 MiB in the parent: inherited wait4 highwater must retain those pages,
while the helper's own post-exec VmHWM must be below half that allocation.
This is a measurement-method regression check, not performance evidence.


The family runner requests the diagnostic and stores the parser fields beside
raw stderr, plus the separate `helper_postexec_hwm_kib` scalar in each record.
Missing/malformed/duplicate diagnostics keep this scalar unavailable; they do
not change the numerical outcome. Incomplete child lifetimes cannot supply it.
A warm workflow currently starts two helper lifetimes, so its duplicate lines
remain unavailable rather than being combined. Aggregate IPC RSS stays unknown.

`report_front.py --helper-memory-lower-bound` is an explicit opt-in. It checks
positive scalars against the exact stderr/opt-in environment, requires a valid
bound for every successful trial in a peer group, and permits only the strict
comparison candidate-RSS-max < helper-HWM-min. This excludes one possible
memory dominator and is recorded separately from the parent lower bound.
It never supplies an IPC memory upper bound or proves IPC domination. Without
the flag, existing historical frontier assessments retain their old behavior.
