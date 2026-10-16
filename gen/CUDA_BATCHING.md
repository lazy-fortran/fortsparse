# Thin Schur update batches

The native CUDA adapter groups the updates from each source panel of width
at most eight. The arithmetic remains the FortSym-generated scalar update in
`spd_schur_scatter_entry`, retained in `spd_schedule.cuh.inc`. The adapter owns
only launch geometry and arena ownership. Wider updates retain cuBLAS.

For one source panel, canonical sorted rows visit destination panels in
increasing order, once each. Every destination owns its disjoint full arena
slice, after the source panel's slice. Within a destination, unique positive
row and column maps give each lower-triangle output one owner. Therefore
updates from this panel can execute concurrently without atomics. Updates
from distinct source panels can share destinations, so their batches remain
ordered on one stream. Setup checks sorted disjoint destination slices and
that every operand read ends before the first destination. Invalid setup
preserves the preceding plan; replacing scatter maps invalidates it.

Setup retains a device task array with cumulative 256-entry block counts.
Each launched block binary-searches its task's end count and calls the existing
generated per-output schedule. Zero-work tasks contribute no blocks, and empty
batches launch nothing. This preserves parallel tiles for large destinations
rather than serializing each task onto one block. A prior failed pivot makes
all subsequent batched writes no-ops, as in the single-task path.

For U thin updates in P source panels, retained storage is
`U*sizeof(DeviceTask)` device bytes and
`sizeof(Context) + batches.capacity()*sizeof(Batch)` host bytes, reported
separately. On the tested 64-bit ABI these task/batch sizes are 80/32 bytes.
Setup also temporarily owns U public task descriptors (72 bytes each),
P+1 host offsets, and U device-format host staging descriptors; these are freed
before numerical factorization. Ordinary process peak RSS includes this
setup peak. The existing maps, factor arena, pivot info and cuBLAS workspace
remain included in device attribution. No numerical call allocates metadata.

Each batch is bounded by INT_MAX blocks and each task by validated arena/map
ranges; the batch setter rejects overflowing prefixes. Setup costs O(U+P),
with O(log r) task lookup per 256-entry block for r updates in that panel.
Thin Schur launches fall from U to the number of nonempty thin source panels;
CUDA activation and C ABI dispatch fall by the same number. This is a launch
count prediction, not an elapsed-time claim. Transfer-inclusive comparisons
and profiler launch counts must be remeasured on the integrated revision.

The registered C ABI test independently computes host dot products for two
destinations with distinct strides/maps, a multi-block task, zero work,
repeated batches, and retained-plan reuse. It also covers alias/bounds errors,
failed pivots, map replacement, and stable host/device allocation counts.
The existing CUDA supernode fixture checks sparse factor/refactor results
against independent known solutions; CPU and generated derivative checks
remain separate.
