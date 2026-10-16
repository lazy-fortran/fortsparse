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
