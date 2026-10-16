# SPD numerical source

`gen_spd_leaves.f90` retains the scalar square-root, division and Schur-update
definitions. FortSym differentiates these expressions and emits CPU and native
CUDA values/JVP/VJP from the same mathematics. The Fortran elemental function
wrappers contain no independent numerical formula. The domains are positive
inputs for the square root and nonzero denominators for division; matrix and
pivot validation belong to the calling solver.

`spd_schedule.f90.inc` retains allocation-free CPU panel factorization,
right-side triangular solve and rectangular Schur-update schedules. The solver
owns task parallelism; the panel schedules use SIMD without nested teams.
Large contractions may use BLAS rather than these reference schedules.
The Schur schedule accepts an optional first-row index (default one), allowing
strided panel tails to share their original contiguous storage. CPU factors
below ten million estimated operations use the generated contraction and
indexed triangular solves to avoid BLAS initialization; CPU panels of at least
96 columns retain BLAS triangular solves.

From the repository root, with a built FortSym at the retained revision:

```sh
bash gen/regenerate_spd.sh /path/to/fortsym
git diff --exit-code -- src/generated/
```

For an out-of-tree FortSym build, pass the generator build directory as the
second argument and the FortSym build directory as the third. The retained
FortSym revision pins FortGen's reciprocal-aware emitter. Exact small negative
powers use typed reciprocals rather than general CUDA `pow`; this is scalar
source shaping, not a faithfully rounded complete quotient or an enclosure.

The generator uses FortSym's native simplifier. Ordinary library builds use
the committed output and require neither FortSym nor a CAS. Compile the
generated Fortran module into `fortsparse`; register
`test/solver/test_fortsparse_spd_leaves.f90` with the existing focused test
helper. The generated CUDA header also retains the per-output contraction and
mapped lower-scatter schedule from `spd_schedule.cuh.inc`. Launch geometry,
residency and stream ownership belong to the optional native CUDA adapter.

The focused test uses known values, independent central differences and
adjoint pairings, a constructed triangular-factor oracle, and dense `matmul`
residual/update oracles. Numerical probes establish behavior, not a proof or
an outward-rounded enclosure.

The optional adapter in `src/native/fsparse_spd_cuda.cu` owns one retained
device arena and nonblocking stream. The legacy C context creation also prepares
cuBLAS and its explicit 1 MiB workspace. Native sparse analysis instead creates
a deferred context, preparing BLAS during setup only when a source panel wider
than eight columns has trailing rows. Thin-only plans therefore require no
cuBLAS handle or workspace. Explicit BLAS operations on a deferred context
require the setup-only prepare call; numerical calls never allocate this state.
The reported allocation bytes count the workspace only after preparation;
CUDA context/cuBLAS internal memory needs separate process/device measurement.
The thin C interface uses arena offsets, retaining panels and operands.

The cooperative factor schedule is retained in `spd_schedule.cuh.inc` and uses
the generated numerical leaves. For source panels up to eight columns, the same
block factors the diagonal and solves every trailing row, including panels
longer than the block's 256 threads. Upper entries and leading-dimension padding
remain untouched. The original diagonal-only panel API preserves its behavior.
Source panels up to eight columns also use a fused generated negative
dot product and mapped scatter, avoiding product scratch and a second launch.
Native GPU setup computes dependency levels from every actual source-to-
destination update edge, rejecting non-forward edges. Stable panel IDs are
grouped by level. Each level factors all thin panels with one cooperative
block per panel, then accumulates their mapped updates in one tiled grid;
wide panels keep serial cuBLAS operations on the same stream. Every update
targets a strictly later level. The generated negative-dot increment feeds
both the original serial scatter and the level scatter. Concurrent thin
sources use FP64 atomic addition, so the reduction order is hardware dependent
and bitwise reproducibility is not promised. Independent original-matrix
residual checks remain the accuracy criterion; this is no mixed-precision
path. Each factor block uses shared local pivot status and atomically reports
failure after its uniform cooperative schedule. Other independent panels may
have executed when a failure is reported; all failed factor output is discarded.
Level descriptors and update tiles are retained at setup and counted in device
bytes, with host launch metadata and wide-panel ID lists counted separately.
Setup uses temporary counting-sort arrays; factor/refactor allocates none.
The original single-panel and non-atomic batch C interfaces remain available.
Wider Schur contractions and trailing triangular solves retain cuBLAS. The fused
schedule skips absent map entries and preserves unused destination cells;
its local reduction order can differ from the vendor GEMM reduction.
Upload/download calls enqueue copies; callers must synchronize before reusing
host buffers and check factor info before consuming a failed panel. Each
context belongs to one caller and one stream. The supernode engine's optional
device argument selects a resident sparse factorization pipeline: setup
retains compact integer scatter maps and an all-panel triangular-solve plan.
Each factorization uploads its assembled panel input once. Completed factors
stay on the device across solves and refactors; the host assembly allocation
is no longer a valid numeric factor after a GPU factorization. Host-workflow
solves transfer only RHS/answer tiles, with at most eight columns per tile.
The public backend owns permutation into/out of factor numbering. This remains
an experimental execution path; competitiveness requires transfer-inclusive
measurements on the declared workload. See [CUDA_SOLVES.md](CUDA_SOLVES.md) for
ownership, memory bounds and independent oracles.

For optional integration, define a CUDA option defaulting off, include
`cmake/FortSparseNativeCUDA.cmake` after creating the library target, and call
`fortsparse_add_cuda_tests()` inside the testing branch. The helper respects
`CMAKE_CUDA_ARCHITECTURES`. The two device fixtures test emitted derivatives
and resident panel/triangular/Schur operations against independent host
oracles; they require an actual CUDA device.

`fortsparse_native_cuda.F90` exposes ISO_C_BINDING interfaces and the
`native_cuda_available` build capability. CUDA interfaces are guarded by
`FORTSPARSE_HAVE_NATIVE_CUDA`; compiling the module without that definition
requires no CUDA compiler, headers or runtime. The optional CMake helper adds
the bindings and definition. The Fortran device fixture verifies pointer,
size and status interoperability through an independently known 2-by-2 factor.
The supernode fixture accepts `--cuda` to exercise sparse/block-boundary
factorization, retained-context refactorization and recovery after a failed
pivot on a real device. The CPU fixture exercises the same mathematical cases
without CUDA.

The CPU sparse schedule also recognizes exact contiguous update intersections.
It lowers their diagonal block to BLAS SYRK and their remaining rectangle to
GEMM directly in the destination factor. Irregular graph intersections retain
the generated mapped-scatter schedule. This is a schedule for the same Schur
update definitions; it does not introduce a second mathematical formulation.

Exact (zero-added-fill) CPU/GPU supernodes use a 64-column cap by default,
with a 128-column cap when mean symbolic column fill exceeds 128. The generated
CPU triangular schedule lowers panels with at least 96 columns and 32 tail
rows to BLAS TRSM; smaller panels retain the generated SIMD reference loops.
The independent larger triangular fixture checks a known solution and unequal
leading dimensions, including untouched padding. This bounded schedule choice
was measured end to end; wider panels trade some upper-triangle storage for
fewer update/packing calls. Matrix residual and implicit derivative contracts
are unchanged.

## Sparse LU update schedules

`gen_lu.f90` owns the Cartesian complex update definition and its scalar
primal/JVP/VJP. Its generic FortSym snippet emission inserts the same primal
roots directly into the allocation-free scatter/gather loops retained in
`lu_indexed_schedule.f90.inc`. Both output components are computed before the
complex value is stored. Scatter indices must be valid and unique; gather
accepts arbitrary valid indices and optionally conjugates coefficients for
adjoint solves. The real LU schedule reuses the generated SPD scatter/gather.
Empty sparse columns are skipped before forming array slices. Division,
scaling, pivot selection, factor storage and memory caps retain their contracts.

The LU leaf fixture checks indexed values against independent complex products,
untouched entries, empty batches, finite differences, and the retained symbolic
JVP/VJP using the ordinary real component pairing. LU factor fixtures also cover
one-element and diagonal matrices with empty L storage, N/T/H solves, column
permutations, re-pivoting and unchanged owned memory across solves.

Witnessed LU pruning keeps the complete numerical L while DFS traverses a
reduced adjacency prefix. An edge from an earlier pivot to row r is removed
only when an actual retained path through the newly completed pivot reaches r.
The witness edge to that pivot remains in the prefix. Pair-partitioning row
indices and coefficients keeps every numerical factor entry aligned; scatter
and triangular solves still use the full column. Exact cancellations cannot
justify a removal because membership uses only retained nonzero L entries.
Adding each new column and deleting these witnessed edges preserves the graph's
transitive closure. This is a conservative transitive reduction, not a new
factorization algorithm. General pruning motivation is described in Demmel et
al., *A Supernodal Approach to Sparse Partial Pivoting* (1999), section 4.2;
the retained paper is indexed in the sparse-dev solver library cards.

The internal low-level factor option defaults to the full graph. The production
native adapter enables the candidate schedule; NATIVE_REFERENCE selects the full
graph. Enabled pruning owns one integer endpoint per column and preflights that
allocation alongside retained factors; disabling it releases that metadata.
Every factor/refactor rebuilds graph state, including after a failed attempt.
The focused pruning fixture checks actual graph partitions against independent
Boolean graph closure, aligned real/complex coefficients, an exact-cancellation
counterexample, pivot permutations, implicit derivatives and allocation caps.
Performance retention requires the controller's matched end-to-end measurements.
