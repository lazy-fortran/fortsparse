# Native CPU matrix classes

Explicit `FORTSPARSE_BACKEND_NATIVE` now dispatches real SPD matrices to the
existing Cholesky engine and other classes to the following new CPU references:

| Value/class | Native factors | Current execution |
|---|---|---|
| Real general | Row-scaled threshold-pivoted sparse LU | Natural or explicit AMD order, retained CSC L/U, serial factorization |
| Complex general | Row-scaled threshold-pivoted sparse LU | Natural or explicit AMD order, retained CSC L/U, serial factorization |
| Complex Hermitian positive definite (`SYM_SPD`) | Sparse complex Cholesky | Natural-order simplicial reference |
| Real symmetric indefinite (`SYM_INDEFINITE` or real `SYM_HERMITIAN`) | Sparse Bunch--Kaufman LDLT | Unrestricted 1x1/2x2 pivots, serial reference |
| Complex Hermitian indefinite (`SYM_HERMITIAN`) | Sparse complex Bunch--Kaufman LDLH | Unrestricted 1x1/2x2 pivots, serial reference |
| Complex transpose-symmetric indefinite (`SYM_INDEFINITE`) | Sparse complex Bunch--Kaufman LDLT | Unrestricted 1x1/2x2 pivots, serial reference |

These are usable native reference families, not production performance admission.
The previous real-SPD CPU frontier results do not extend to these classes.
HPD Cholesky and complex LDLH validate full conjugate symmetry. Real and complex
LDLT validate full ordinary symmetry. Invalid declared structures return
`INVALID_MATRIX`; HPD Cholesky also checks positive pivots.

Existing CSC factor/refactor, in-place/vector/block solve and lifetime APIs are
preserved. LU raw CSC input is borrowed directly, including under a host cap;
no intermediate input value copy or split-complex conversion is needed.
Structured raw CSC factorization remains unsupported in this first adapter.
Default AUTO prefers eligible in-process external backends when built; a
native-only build uses native AUTO. This is capability selection, not additional
performance admission.

The new references enforce `host_memory_limit_bytes` across their owned arrays,
including factor growth replacement overlap, and support RAM-only LOW_MEMORY.
LU uses exact-needed growth instead of geometric spare capacity in LOW_MEMORY;
HPD already uses compact exact structural storage. Separate contexts own their
mutable work independently. Strict-budget and LOW_MEMORY execution require
natural ordering; these matrix classes currently execute on CPU. LU, HPD and
real/complex LDL support `solver%refine=.true.` with
retained original-matrix residuals, a 32-epsilon normwise backward-error target
and at most three corrections. Statistics report correction steps and final
backward error. See [native refinement](native-refinement.md) for failure,
allocation and refactor-cap semantics. Production blocking/parallel
factorization, conditioning diagnostics and ordering beyond the optional LU AMD
path remain planned.

`pivot_tolerance` selects LU's finite threshold in [0,1]; a negative value uses
the native default 1.0 (maximum-candidate row pivoting after scaling). A numerical
refactor always selects pivots anew. No derivative is taken through pivot choices.

## Optional general LU ordering

For uncapped real and complex general matrices, production native LU accepts
`solver%options%ordering=FORTSPARSE_ORDER_AMD` when SuiteSparse AMD is built.
AMD accepts the canonical unsymmetric CSC graph and orders its structural
`A+A^T`; numerical factorization borrows the original CSC values and applies the
column permutation through its retained sparse factors. DEFAULT remains natural
ordering. NATIVE_REFERENCE remains natural and rejects explicit AMD.

Value-only refactors reuse the retained column order when the canonical pattern
and requested ordering strategy match; numerical row pivots are selected anew.
A changed pattern recomputes AMD, and a changed ordering strategy invalidates
pattern-reuse reporting. Normal, transpose and adjoint solves, refinement and
implicit solve derivatives act on the original matrix.

AMD's private workspace has no solver-owned hard-cap contract. Explicit AMD
with a positive host cap or LOW_MEMORY returns BACKEND_UNAVAILABLE before
calling AMD; the native-only build without AMD returns the same availability
status. AMD statistics report the actual ordering and exact retained owned
payload, but complete host peak bytes and factor memory are unknown (`-1`), and
`strict_memory_supported` is false. Use measured process RSS for comparisons;
these changes make no performance-admission claim.

## Retained transpose and adjoint solves

```fortran
call sparse_solve_op(solver, b, x, FORTSPARSE_OP_TRANSPOSE, status)
call sparse_solve_op(solver, b, x, FORTSPARSE_OP_ADJOINT, status)
call sparse_solve_block_op(solver, Bmany, Xmany, FORTSPARSE_OP_ADJOINT, status)
```

NORMAL=0, TRANSPOSE=1 and ADJOINT=2. Real symmetric paths have the same action
for all three. Complex Hermitian factors, including HPD, have a distinct ordinary
transpose action; complex transpose-symmetric factors have a distinct adjoint
action.
Existing `sparse_solve` remains a normal solve. SuperLU, MUMPS and the current
UMFPACK helper also support these retained operations and operation-aware
refinement. An older factor-capable UMFPACK helper remains usable for normal
solves; capability negotiation rejects unsupported transpose/adjoint requests
with BACKEND_UNAVAILABLE while preserving the retained normal factors.

For equation-level reverse products, existing callers may still pass a separately
factored transpose/adjoint matrix. To reuse primal factors on a backend that supports the retained operation:

```fortran
call sparse_solve_vjp(solver, A, x, x_bar, b_bar, A_values_bar, status, &
                     use_primal_factors=.true.)
```

The default remains `.false.` for compatibility. Real and complex JVPs reuse the
caller output and retained solve workspace, avoiding tangent RHS allocations.
Independent public API tests exercise N/T/H, multiple RHS, original-matrix known
solutions, real-pairing adjoint identities, hard caps and the exact rational
oracle of the retained NEO-2 mini matrix. Those small tests establish algebraic
behavior, not a full NEO-2 kinetic run or performance parity.
