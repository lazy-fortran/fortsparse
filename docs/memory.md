# Native memory planning and limits

`host_memory_limit_bytes` and `device_memory_limit_bytes` are 64-bit byte
limits; zero means unlimited. `FORTSPARSE_MEMORY_LOW` requests the RAM-only
low-memory schedule. Defaults preserve existing backend selection and execution.
UMFPACK's older `memory_budget_bytes` is an **estimate-only numeric gate** and
has different scope. It is not a hard allocation limit.

## First supported schedule

The strict host cap currently supports natural-order CPU real-SPD simplicial
Cholesky. Selecting a host cap or LOW_MEMORY selects this native schedule.
It trades blocked execution and parallel scratch for predictable compact storage;
it has no production-performance admission yet. Explicit AMD/METIS ordering,
CUDA with strict caps, external backends with strict caps and strict-budget raw
CSC conversion return BACKEND_UNAVAILABLE rather than ignoring the request.
They remain roadmap work. Normal uncapped paths retain their existing behavior.

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

Analysis preserves any retained factorization. It checks structure/symmetry and
counts fill; it does not prove positive definiteness. The numeric factorization
still checks every pivot. Planning itself is capped and may fail for insufficient
graph scratch before a full layout can be calculated. Cap changes apply on the
next factor/refactor attempt, and invalidate reuse when the schedule changes.

## What the bytes mean

Plan categories are disjoint: factor payload, other persistent array payload and
transient workspace. `peak_host_bytes` covers their peak simultaneous lifetime.
For the supported natural-order real reference, with n columns, m supplied CSC
entries and L structural Cholesky entries, the exact array layout is
`12*L + 4*m + 36*n + 8` bytes on the supported 32-bit-index/FP64 build.
The factor value array contributes `8*L`. Analysis needs only three integer
n-vectors while planning, released before factorization. Scalar and block solves
reuse retained scratch without per-call allocations; the real JVP uses its
caller-owned output as RHS. A backend-owned `sparse_vector` also consumes the
cap and returns a disassociated pointer if its allocation cannot fit.

Statistics report attributed current/peak payload and required/available bytes
on budget failure. RESOURCE_LIMIT is distinct from operating-system allocation
failure (OUT_OF_MEMORY). A failed factor attempt leaves the public solver
unfactored; callers may change limits and retry or destroy it. Different handles
have independent limits. The limit applies per handle, not to their combined sum.

The plan excludes caller-owned A/B/X, fixed object descriptors, allocator page
rounding, runtime/library code, stacks and process RSS. Exact layout means exact
array payload, not an exact RSS prediction. New pivoted paths must label unknown
fill as an estimate and cap every growth allocation, including old/new overlap.
Opaque vendor allocations cannot acquire a strict guarantee from an estimate.

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
