# Retained SPD triangular solve schedule

The optional CUDA supernode engine solves `L y = b` and `L^T x = y` in
factor numbering. Its existing dependency levels are reused for **all**
panels, including wide panels. `spd_schedule.cuh.inc` owns the indexed GPU
schedule; its only numerical operations are the retained FortSym-generated
update and divide leaves. Regenerate the committed CUDA header with the same
`gen/regenerate_spd.sh` recipe as the factor kernels. No differentiation through
triangular-solve iterations is introduced: the public implicit JVP/VJP still
solve their equation-level systems with these same factors.

For every panel, setup validates a width in `[1,128]`, disjoint factor storage,
and a consecutive partition of the `n` factor columns. Tail row indices are
strictly increasing, one-based, outside the panel, and owned by a strictly
later dependency level. A stable counting sort retains descriptors by level.
Setup temporary ownership/counting arrays are freed before factorization.

In an ascending level, each cooperative block solves its panel's diagonal
system in shared memory, then scatters the negative trailing product into
later rows. Different blocks own disjoint panel RHS rows. The strict-level
condition prevents a scatter from racing a panel RHS read/write in the same
level. Multiple predecessors can update a later row, so these FP64 scatters
use atomic addition. In descending levels, every block gathers already solved
later rows, applies the diagonal transpose solve, and writes only its owned
rows. Stream order provides the barriers between levels; no inter-block spin
barrier is needed. Atomic summation order is hardware dependent. Accuracy is
checked against the original matrix, without a bitwise-reproducibility claim.

`fsparse_spd_cuda_set_solve_plan` installs the plan at setup.
`fsparse_spd_cuda_solve(ctx, rhs_offset, ld, nrhs)` accepts zero through eight
RHS columns in the existing arena, rejects overlap with the retained factor
region, and requires a successful synchronized `factor_info` after the latest
factor upload/stage. A failed pivot invalidates solve readiness. Numerical
calls allocate no project-owned memory. The Fortran supernode `solve_block`
accepts any number of contiguous host RHS columns, processes tiles of eight,
and synchronizes each tile before its host buffer can expire. Scalar `solve`
uses the same path. The existing public backend owns original/factor numbering
and alias-safe staging. The CPU numerical solve path is preserved.

Factors remain on device after `factor_info`; factorization no longer downloads
them. The original host panel allocation remains the assembly input for the
next refactor. RHS input and answer output are the only solve transfers. This
is a host-workflow implementation, with no new public device-resident API.

For `P` panels, `T` packed tail row indices and `D` levels, retained additional
device metadata is `P*sizeof(SolvePanel) + T*sizeof(int)`. The implementation
reports these actual types through its allocation counter. Host launch
metadata is `capacity(solve_levels)*sizeof(SolveLevel)` plus the enlarged
context, also counted. On the current 64-bit ABI the descriptor sizes are
40 and 16 bytes. The arena reserves `max(existing_factor_scratch, 8*n)` doubles
after the factor region, rather than adding a separate RHS allocation. Each
block has a bounded `128*8`-double shared tile. CUDA runtime/vendor internal
allocation is separate from these attributed counts; no strict memory-budget
or complete process/device peak claim follows from the attributed counters.

The focused C ABI oracle constructs independent dense `A = L L^T`, supplies
known lower factors with poisoned upper entries/padding, and checks original-A
residuals and manufactured solutions for widths 1/8/32/128, branching levels
with overlapping forward destinations, repeated changed factors, RHS 0/1/8,
padded RHS strides, invalid plans/bounds, failed pivots, and stable retained
bytes. It does not use the factor kernel to manufacture its expected answer.
The production supernode fixtures factor original CSC matrices, check nine
RHS columns across the tile boundary, nonuniform SPD refactors and recovery.
`test_fortsparse_native_cuda_contract` runs public permutation, scalar/block,
implicit JVP/VJP pairing and failure oracles with device zero. Together these
are correctness evidence, not a rigorous enclosure or performance admission.
