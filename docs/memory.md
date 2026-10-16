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
Natural DEFAULT selection is unchanged. SPD nested dissection, CUDA, external
backends and structured raw CSC conversion remain unsupported under strict caps.
Normal uncapped paths retain their existing execution choices.

The standalone [real CSR ILUTP preconditioner](ilutp.md#payload-cap-and-recovery)
also bounds every owned array, including fixed factor capacities, setup scratch
and retained apply workspace. Its explicit fill caps give a closed setup-payload
formula. Rebuild frees the old arrays before allocating replacements; failure
leaves an unfactored context. Caller storage, fixed descriptors, allocator
overhead and process RSS are excluded. This does not establish SPEC application
memory or preconditioner performance admission.

The first exact pre-factorization planning API describes the compact real-SPD
layout; select LOW_MEMORY or NATIVE_REFERENCE. It does not claim to predict
blocked/ordered/pivot-dependent layouts. Production planning remains roadmap work.

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
alongside the planned new layout. It checks structure/symmetry and
counts fill; it does not prove positive definiteness. The numeric factorization
still checks every pivot. Planning itself is capped and may fail for insufficient
graph scratch before a full layout can be calculated. Cap changes apply on the
next factor/refactor attempt, and invalidate reuse when the schedule changes.

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

The factor value array contributes `8*L`. Analysis needs only three integer
n-vectors while planning, released before factorization. Scalar and block solves
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
