# Compressed CSR and GMRES iterative refinement

`fortsparse_precision_csr:precision_csr_t` is an additive real square CSR
interface. Existing CSC/direct-solver contexts are unchanged. It reuses
[FortNum](https://github.com/lazy-fortran/fortnum)'s precision vectors, stage
policies, refinement controller and FGMRES implementation.

`%store(n,row,column,values,storage,ok)` accepts one-based CSR, real64 values
and a storage identity from `fortnum_precision`. It validates monotonic row
pointers and column bounds before replacing the previous matrix. Duplicate
columns are summed by the ordinary SpMV accumulation. `%matvec(x,y,compute,ok)`
accepts real64 boundaries and fp32/fp64 arithmetic; fp32 arithmetic normalizes
input and matrix values before accumulation. Nonfinite output fails explicitly.
`%bytes()` counts value and index payload, excluding object and allocator
metadata. Large row loops use OpenMP. Matrix conversion allocates during setup;
SpMV allocates no workspace.

`csr_gmres_ir(a,b,x,policy,tolerance,max_steps,inner_iterations,restart,`
`info,steps,residual)` requires original **fp64** CSR storage. The original
operator is retained for every outer residual and final acceptance. A separate
rounded matrix follows `STAGE_PRECONDITION`, and is rebuilt after escalation.
Each correction uses existing restarted flexible GMRES with a selectable
low-precision operator and fp64 basis/orthogonalization, then fp64 defect
correction. This is not a low-precision direct factorization or an entirely
fp16/fp32 Krylov basis. The inner solve targets a 5% correction residual;
the caller chooses the outer tolerance. See FortNum for return statuses.

CSR fp32 computation with original fp64 storage requires entries representable
in fp32; use narrowed, scaled storage for larger magnitudes. A single matrix
scale does not preserve arbitrarily small entries in a large dynamic range.
No native BF16 dot-product, native FP16 arithmetic, CUDA backend, hard allocation
cap or performance admission is claimed by this API. Existing capped direct
solvers retain their separate contracts. Tests use an independently evaluated
nonsymmetric tridiagonal system, both arithmetic kinds, all four storage
formats and scaled coefficients far above the binary16 range.

## Native direct factors with adaptive storage

The umbrella `fortsparse` module exports `mixed_lu_t` (real), `mixed_lu_z_t`
(complex), `mixed_solve_report_t`, `mixed_factor_batch`, `mixed_solve_batch`
and `block_jacobi_apply`. These additive CPU contexts use the native sparse
threshold-pivoted LU; external sparse solvers are unnecessary.

```fortran
use fortsparse
use fortnum_precision, only: FP32
! a is canonical square CSC; b and x are real(dp) or complex(dp).
type(mixed_lu_t) :: factor
type(mixed_solve_report_t) :: report
type(fortsparse_status_t) :: status
call factor%factor(a,status,storage=FP32)
call factor%solve(b,x,report,status)
call factor%free()
```

Always check `status%code`. Setup factorizes in fp64, estimates conditioning,
then optionally narrows retained L/U values to fp32. It does **not** perform
fp32 factorization. Indices, row scales, original CSC, residuals and triangular
arithmetic remain in their original kinds. One preallocated column buffer
widens retained values for the existing generated compensated solve leaves.
The default storage is FP64. Requested FP32 storage remains FP64 when either
estimated reciprocal condition is at most `100*epsilon(1.0_sp)` or factors
are outside the safe fp32 range. `report%storage` records the actual choice.
This conservative heuristic does not guarantee convergence; rounding can still
make otherwise usable factors unsuitable.

`factor(a,status,storage,column_order,device,memory_limit_bytes)` accepts an
optional valid column permutation, such as a fill-reducing order prepared by
the caller. Existing public `sparse_factor` retains its AMD/ND options. Device
`-1` selects CPU; other devices return `BACKEND_UNAVAILABLE`. FP16/BF16 direct
factor storage is unsupported. The optional positive cap reserves the original
CSC and all checker vectors alongside the low-level factor peak, including
narrowing overlap. `%bytes()` reports retained owned array payload, excluding
caller matrices, object descriptors, allocator and runtime overhead. The wrapper
has no pre-factor fill prediction or GPU execution; existing solver planning
and CPU LDLT/Cholesky interfaces keep their own contracts.

`solve(b,x,report,status,operation)` accepts normal (0), transpose (1) or adjoint
(2), reusing factors. It independently evaluates residuals against the original
fp64 CSC and performs at most ten fp64 defect corrections to a normwise
backward-error target of 32 fp64 epsilons. Stagnation, nonfinite errors or the
iteration limit return `SINGULAR` with an explicit refinement-failure message;
`report%converged` remains false. Refactor with FP64 storage after this failure.
This status describes unusable corrections, rather than proving exact matrix
singularity. No refactor or allocation occurs inside retained solves.

Reports include normwise and componentwise backward errors, correction count,
actual storage and a reciprocal condition **estimate**. A five-step Hager-style
inverse-norm estimator uses normal/adjoint solves at setup for each orientation.
The reported conditioning estimate concerns the original fp64 factors, before
narrowing; it is reused per solve. It can overestimate reciprocal conditioning.
Neither conditioning nor backward error is a rigorous forward-error certificate.

`solve_block(b,x,reports,statuses,operation)` returns one report/status per RHS.
`mixed_factor_batch(contexts,matrices,statuses,storage)` factors independent
real or complex CSC systems in parallel. Each context owns its scratch; there
is no aggregate batch cap. Set individual caps by calling `%factor` separately
when budgeting a batch. `mixed_solve_batch(contexts,offsets,b,x,reports,statuses)`
solves independent systems with one-based offsets into packed vectors, including
different dimensions per system. OpenMP parallelizes over contexts. A context
must never be used concurrently by multiple callers; block RHS execution reuses
its scratch sequentially.

For Krylov preconditioners, `%apply(rhs,status,operation)` applies the retained
approximate inverse without refinement. `block_jacobi_apply(contexts,offsets,`
`b,x,statuses)` applies independently factored diagonal blocks in parallel.
Callers assemble those blocks and own the outer operator and partition; block
precision adapts independently at factor setup. These are preconditioner seams,
not a new Krylov solver. Approximate apply has no exact-solve error guarantee.

Existing native symmetric LDLT remains the inertia interface. Compressed LU
never reports inertia. Use `sparse_inertia` with the documented reliability state;
a near-zero uncertainty must remain unresolved in a stability mode count.
