# Native CPU matrix classes

Explicit `FORTSPARSE_BACKEND_NATIVE` now dispatches real SPD matrices to the
existing Cholesky engine and other classes to the following new CPU references:

| Value/class | Native factors | Current execution |
|---|---|---|
| Real general | Row-scaled threshold-pivoted sparse LU | Natural order, retained CSC L/U, serial factorization |
| Complex general | Row-scaled threshold-pivoted sparse LU | Natural order, retained CSC L/U, serial factorization |
| Complex Hermitian positive definite (`SYM_SPD`) | Sparse complex Cholesky | Natural-order simplicial reference |
| Real symmetric indefinite (`SYM_INDEFINITE` or real `SYM_HERMITIAN`) | Sparse Bunch--Kaufman LDLT | Unrestricted 1x1/2x2 pivots, serial reference |
| Complex symmetric/Hermitian indefinite | General complex sparse LU | Full supplied CSC; specialized pivoted complex LDL remains planned |

These are usable native reference families, not production performance admission.
The previous real-SPD CPU frontier results do not extend to these classes.
Complex structural hints other than SPD select general LU, not a symmetry-specific
factorization. HPD Cholesky validates full conjugate symmetry and nonpositive
pivots. The real indefinite factorization validates full ordinary symmetry.

Existing CSC factor/refactor, in-place/vector/block solve and lifetime APIs are
preserved. LU raw CSC input is borrowed directly, including under a host cap;
no intermediate input value copy or split-complex conversion is needed.
Structured raw CSC factorization remains unsupported in this first adapter.
Default AUTO selection still uses the external backends.

The new references enforce `host_memory_limit_bytes` across their owned arrays,
including factor growth replacement overlap, and support RAM-only LOW_MEMORY.
LU uses exact-needed growth instead of geometric spare capacity in LOW_MEMORY;
HPD already uses compact exact structural storage. Separate contexts own their
mutable work independently. Their first adapters require natural ordering and
CPU execution; requesting refinement currently returns BACKEND_UNAVAILABLE
rather than silently ignoring it. Production blocking/parallel factorization,
refinement/conditioning diagnostics and broader ordering are remaining work.

`pivot_tolerance` selects LU's finite threshold in [0,1]; a negative value uses
the native default 1.0 (maximum-candidate row pivoting after scaling). A numerical
refactor always selects pivots anew. No derivative is taken through pivot choices.

## Retained transpose and adjoint solves

```fortran
call sparse_solve_op(solver, b, x, FORTSPARSE_OP_TRANSPOSE, status)
call sparse_solve_op(solver, b, x, FORTSPARSE_OP_ADJOINT, status)
call sparse_solve_block_op(solver, Bmany, Xmany, FORTSPARSE_OP_ADJOINT, status)
```

NORMAL=0, TRANSPOSE=1 and ADJOINT=2. Real symmetric paths have the same action
for all three; complex HPD still needs a distinct ordinary transpose action.
Existing `sparse_solve` remains a normal solve. External backends that do not
expose retained transpose operations report BACKEND_UNAVAILABLE for those calls.

For equation-level reverse products, existing callers may still pass a separately
factored transpose/adjoint matrix. To reuse the primal native factors instead:

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
