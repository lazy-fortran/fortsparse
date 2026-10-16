# Native memory planning and limits

`host_memory_limit_bytes` and `device_memory_limit_bytes` are 64-bit byte
limits; zero means unlimited. `FORTSPARSE_MEMORY_LOW` requests the RAM-only
low-memory schedule. Defaults preserve existing backend selection and execution.
UMFPACK's older `memory_budget_bytes` is an **estimate-only numeric gate** and
has different scope. It is not a hard allocation limit.

## First supported schedule

The strict host cap supports natural-order and explicit-AMD CPU real-SPD
simplicial and blocked Cholesky. Bounded AMD requires the audited optional
AMD 3.3.4 caller-workspace ABI; unsupported versions fail closed. AUTO attempts
bounded generated blocked execution, reduces task scratch
concurrency when needed and can use the compact reference if panels cannot fit.
Statistics report `compact_memory_used`, `budget_reduced_threads`, the actual
team and `native_factor_family`. LOW_MEMORY explicitly selects the compact
reference. These bounded schedules need separate performance admission.

The real/complex LU, complex HPD and real/complex indefinite references also
enforce host caps; see [native CPU families](native-cpu.md). General LU supports
explicit AMD under caps and LOW_MEMORY through the audited AMD 3.3.4
[caller-workspace path](native-amd-bounded.md), including ordering scratch,
refinement reservations and retained-factor overlap. General LU raw input is
borrowed and capped. Real SPD also uses that shared AMD workspace, accounting
permutations, reordered CSC, symbolic scratch, panels and refinement together.
Natural DEFAULT selection is unchanged. SPD nested dissection and external
backends remain unsupported under strict host caps. CUDA supports the separate
device-owned payload cap described below; CUDA host caps and LOW_MEMORY
remain unsupported.
Normal uncapped paths retain their existing execution choices.

Structured RAW CSC factor/refactor supports real SPD/indefinite and complex
HPD/Hermitian-indefinite/transpose-symmetric inputs with host caps and LOW_MEMORY.
The caller arrays stay borrowed. A temporary owned canonical CSC copy contributes
`4*(ncol+1+nnz) + value_bytes*nnz` bytes, with `value_bytes=8` for real and 16
for complex. Its allocation is admitted alongside all old retained factors,
checker arrays and user vectors before copying. The remaining cap bounds the CSC
backend's factor growth, ordering and refinement; the copy stays charged until
that call returns, then is released. Reported current payload excludes the released
copy; the peak includes its simultaneous lifetime. Factor options are restored
on success and failure, so subsequent user vectors can use that released space.
Malformed input, a refused cap or allocation failure leaves the attempted RAW
factorization unusable, with later factorization recovery supported. Existing
ordering/device restrictions still apply; uncapped SPD AUTO and opaque ordering
may report an unknown peak. This is array-payload accounting, not a process-RSS bound.

The standalone [real CSR ILUTP preconditioner](ilutp.md#payload-cap-and-recovery)
also bounds every owned array, including fixed factor capacities, setup scratch
and retained apply workspace. Its explicit fill caps give a closed setup-payload
formula. Rebuild frees the old arrays before allocating replacements; failure
leaves an unfactored context. Caller storage, fixed descriptors, allocator
overhead and process RSS are excluded. This does not establish SPEC application
memory or preconditioner performance admission.

The pre-factorization planning API describes compact real-SPD layouts. Natural
ordering supports LOW_MEMORY or NATIVE_REFERENCE. Explicit AMD planning requires
LOW_MEMORY, including with NATIVE_REFERENCE, and the audited bounded AMD ABI.
It uses the production permutation and independently counts its symbolic fill.
Uncapped reference AMD without LOW_MEMORY keeps its existing execution policy;
its opaque ordering workspace is not planned. Blocked, ND, GPU and
pivot-dependent planning remain unavailable.

```fortran
use fortsparse
type(sparse_solver_t) :: solver
type(sparse_memory_plan_t) :: plan
type(fortsparse_status_t) :: status
! A is caller-owned, full symmetric canonical CSC.
solver%backend_id = FORTSPARSE_BACKEND_NATIVE
solver%options%symmetry = FORTSPARSE_SYM_SPD
solver%options%memory_policy = FORTSPARSE_MEMORY_LOW
solver%options%host_memory_limit_bytes = 128_i8*1024_i8*1024_i8
call sparse_analyze(solver, A, plan, status, nrhs=8)
! Check status before using the plan or attempting factorization.
if (status_ok(status)) call sparse_factor(solver, A, status)
```

Analysis preserves any retained factorization. Planning scratch shares its cap
with retained budgeted factors; release them before planning when coexistence
cannot fit. Existing unbudgeted or other-family factors must be released before
strict planning. The current fit check conservatively includes retained payload
alongside the planned new layout. An ordered plan reports EXACT_LAYOUT for a
fresh factor and BOUND when this conservative old-plus-new fit is requested;
its byte fields describe the new factor, not the sum. It checks structure/symmetry and
counts fill; it does not prove positive definiteness. The numeric factorization
still checks every pivot. Planning itself is capped and may fail for insufficient
graph scratch before a full layout can be calculated. Cap changes apply on the
next factor/refactor attempt, and invalidate reuse when the schedule changes.

## CUDA device-owned payload cap

Native real-SPD CSC CUDA execution accepts `device_memory_limit_bytes` with
`memory_policy=FORTSPARSE_MEMORY_AUTO` and no host cap. Zero is unlimited.
This limits the sum of explicit device allocation requests owned by the handle:
the FP64 factor/RHS/update arena, pivot-info integer, optional 1 MiB cuBLAS
workspace, maps, Schur/level task metadata, level panels and solve panels/rows.
The arena retains eight RHS columns and the largest required wide-panel update
scratch. Existing generated numerical kernels and arithmetic are unchanged.

Every explicit device allocation checks signed 64-bit/`size_t` ranges and
reserves against the cap before calling CUDA. Transactional Schur, level and
solve-plan replacements include simultaneous old/new arrays. Map replacement
keeps its existing destructive behavior: synchronize, discard dependent plans
and old maps, then allocate new maps. Failure cleanup releases every successful
allocation; no CPU or external-factor fallback bypasses a device refusal.
`FORTSPARSE_RESOURCE_LIMIT` means budget refusal, while CUDA allocation failure
reports `FORTSPARSE_OUT_OF_MEMORY`. A failed public factor/refactor is unusable;
changing the cap and factoring again supports recovery. A changed cap rebuilds
the plan, while same-policy numeric refactorization reuses device allocations.

Statistics keep separate scopes: `device_bytes` is current explicit device
payload, `device_peak_bytes` is the successfully allocated live high-water,
and `device_budget_required_bytes`/`device_budget_available_bytes` describe
the refused allocation increment and remaining quota at that point.
`strict_device_memory_supported` identifies this device-payload accounting.
The existing `strict_memory_supported` continues to describe host accounting;
it is not set by a device cap. Refusal diagnostics and successful high-water
survive failed-build cleanup, whose current device payload is zero.

This is not a total-GPU-memory or process-RSS cap. CUDA context/stream and
cuBLAS internal allocations, driver pools, allocator granularity and all host
storage are excluded and need separate headroom. Retained device payload does
not measure those quantities. Host setup contains overlapping Fortran arrays,
C++ vectors and old/new plans that are not yet capped; its peak stays unknown.
CUDA host caps and LOW_MEMORY fail explicitly rather than silently ignoring
the requested policy. GPU complex/indefinite/general LU remain unsupported.
The pre-factor planning API does not yet predict CUDA layouts.

The optional C CUDA interface preserves existing entry points. Additive
`fsparse_spd_cuda_create_deferred_limited` accepts the limit and an optional
output `fsparse_spd_cuda_memory` report, including failed-create diagnostics;
`fsparse_spd_cuda_memory_stats` queries a live context. C status -2 is a cap
refusal and -3 is byte-count overflow. The existing general-LU C API and its
memory-report struct are unchanged; they do not expose GPU execution.

## What the bytes mean

Plan categories are disjoint: factor payload, other persistent array payload and
transient workspace. `peak_host_bytes` covers their peak simultaneous lifetime.
For the supported natural-order real reference, with n columns, m supplied CSC
entries and L structural Cholesky entries, the unrefined exact array layout is
`12*L + 4*m + 36*n + 8` bytes on the supported 32-bit-index/FP64 build.
With `solver%refine=.true.`, add `4*(n+1+m) + 8*(m+2*n)` for the retained
original CSC and two checker vectors. The returned plan already includes this
quota and checks its overlap with planning scratch; do not add it twice.
See [refinement](native-refinement.md) for the numerical contract.

For explicit AMD with LOW_MEMORY, the retained unrefined compact layout is
`12*L + 20*m + 48*n + 12` bytes. Add the same refinement quota above. The exact
peak also considers the audited AMD arena plus its permutation outputs, the
reordered CSC construction cursor, and the production symbolic-count scratch
used at n >= 4096. These phases have different lifetimes; `workspace_host_bytes`
is the peak increment above retained factor plus persistent payload, rather than
the sum of all transient allocations. The plan runs the actual production AMD
ordering under the cap, releases that temporary layout after counting, and
preserves the existing factors and their statistics even on failure.

The factor value array contributes `8*L`. Natural analysis needs three integer
n-vectors while planning; AMD analysis additionally retains a temporary ordered
CSC and permutation until counting finishes. Scalar and block solves
reuse retained scratch without per-call allocations; the real JVP uses its
caller-owned output as RHS. A backend-owned `sparse_vector` also consumes the
cap and returns a disassociated pointer if its allocation cannot fit.

Statistics report attributed current/peak payload and required/available bytes
on budget failure. Following an allocation failure, the retained attempted peak
is an admission bound; successfully allocated live payload is still accounted,
and freed analysis reports zero live bytes. RESOURCE_LIMIT is distinct from operating-system allocation
failure (OUT_OF_MEMORY). A failed factor attempt leaves the public solver
unfactored; callers may change limits and retry or destroy it. Different handles
have independent limits. The limit applies per handle, not to their combined sum.

The plan excludes caller-owned A/B/X, fixed object descriptors, allocator page
rounding, runtime/library code, stacks and process RSS. Exact layout means exact
array payload, not an exact RSS prediction. New pivoted paths must label unknown
fill as an estimate and cap every growth allocation, including old/new overlap.
Opaque vendor allocations cannot acquire a strict guarantee from an estimate.
A large symbolic-L array may use one private mmap; its requested payload counts
toward the cap. OS page rounding can add less than one system page for that
mapping and requires separate process headroom, just like heap rounding.
The symbolic count routine uses ordinary heap arrays with no private mapping.

## Cluster-node sizing

Add the solver peak to the application's other simultaneous memory: matrix
assembly and conversions, retained input/output and derivative arrays, other
solver handles, I/O and application state. For this CSC representation alone,
the input is `4*(n+1) + 12*m` bytes; each real n-by-r RHS/output is `8*n*r`.
Add measured runtime/allocator/thread overhead and explicit application-selected
headroom. Include concurrent processes and shared-memory mappings according to
the scheduler's accounting rules. Never assume summed RSS equals unique memory.

Validate predictions using representative complete jobs, cold and reused
lifetimes, and measured process/job peaks. Retain hardware, compiler, library,
thread count, cap and prediction quality. A solver payload cap supports capacity
planning but is not an OS process/job cap; use the site's existing scheduler or
cgroup limit when a total-job limit is required. No disk-backed fallback is
provided. A cap that cannot hold the factors produces a clear failure.
