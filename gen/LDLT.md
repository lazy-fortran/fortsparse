# Native real symmetric indefinite reference

`fortsparse_native_ldlt` implements a serial sparse CPU reference for
`P^T A P = L D L^T`, with unit lower `L` and scalar/2-by-2 diagonal blocks.
The pivot decisions use unrestricted classical Bunch--Kaufman diagonal
pivoting with `(1+sqrt(17))/8`, following the lower-triangle algorithm described
in [LAPACK DSYTF2](https://www.netlib.org/lapack/double/dsytf2.f).
Symmetric permutations are maps on original vertex labels. Active Schur rows
and factor columns are sparse sorted index/value arrays; no global dense
matrix is materialized. The real transpose solve is the normal solve.

The reference accepts finite, full symmetric, strictly sorted unique CSC.
Missing diagonal entries mean zero. Symmetric stored nonzeros must have equal
values; no symmetrization or triangle completion is performed. Row maxima
include every remaining vertex; this is not restricted supernode pivoting,
unpivoted LDL or a renamed LU fallback. Structural fill can still become dense
when elimination demands it. The implementation has no parallel or performance
admission claim, does no equilibration/refinement, and reports no forward bound
or deterministic certificate. Nonfinite intermediate arithmetic returns an
explicit failure; numerical singular-pivot status does not prove mathematical
singularity of the exact input.

The API is `native_ldlt_t` with `factor`, `refactor`, `solve`, `solve_block`,
`free`, `permutation`, `retained_bytes`, `peak_bytes` and `factor_nnz` methods.
`factor(A,status,host_limit_bytes,low_memory)` accepts an optional int64 cap;
zero means unlimited. `refactor` checks the retained CSC pattern, drops the old
numeric factors and repivots new values while keeping the pattern. Invalid
input/pattern is rejected before dropping good factors; a numerical/refactor
failure leaves no usable factors. Recovery with the retained pattern is allowed.
`solve` and `solve_block` accept optional logical `transpose`; RHS columns are
processed using one retained workspace and allocate nothing.

The cap covers solver-owned dynamic array payload, including derived-type
array descriptors, copied pattern, permutations, `D`, sparse Schur/factor row
capacities, scratch and old/new row payload overlap during growth. Caller-owned
input/output arrays, the fixed-size solver object and allocator/runtime overhead
are outside that attributed-payload cap. Each allocation checks its complete
requested payload before allocation. LOW_MEMORY uses exact row growth; the
ordinary reference reserves bounded spare capacity. Numeric factorization may
grow sparse storage and is therefore not an allocation-free production kernel.
Refactor deliberately frees old numeric storage before building new factors;
there is no uncounted retained/new factor overlap. Successful factorization frees
Schur/scratch storage; reported retained payload includes the solve workspace.
`factor_nnz` counts stored strict-lower entries plus diagonal entries and one
entry for each 2-by-2 off-diagonal block. It counts structural stored zeros.

`gen_ldlt.f90` defines the normalized 2-by-2 determinant and solve, and the
rank-2 Schur subtraction. FortSym emits their primal/JVP/VJP leaves together.
The scale is the maximum absolute pivot-block entry; division/scalar updates
reuse the existing generated real leaves. Pivoting, indexing, sparse structure
and execution order are handwritten orchestration. Derivatives of the complete
solve belong at the equation level, not through discrete pivot choices.

Regenerate with the pinned FortSym version using:

```sh
bash gen/regenerate_ldlt.sh /path/to/fortsym
```

Ordinary builds compile the checked-in generated module without FortSym.
`test_fortsparse_native_ldlt` checks zero-diagonal 2-by-2 and nonadjacent pivots,
symmetric interchanges, same-pattern repivot/recovery, normal/transpose block
solves, deterministic mixed-sign matrices, known solutions and independent dense
residuals. It checks scaled 2-by-2 blocks at 1e-300 and 1e300, sparse fill growth,
cap success/failure at the actual recorded peak, low-memory equivalence and a
3000-row disconnected diagonal with a 3 MB cap. Generated pair derivatives are
checked against the implicit 2-by-2 equation and adjoint pairing.
