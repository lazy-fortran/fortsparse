# Complex sparse LDL references

`native_ldlt_z_t` supplies a CPU sparse Bunch--Kaufman reference for two
separately validated classes. Hermitian input satisfies `A=conjg(transpose(A))`
and is factored as `P^T A P=L D L^H`; complex transpose-symmetric input satisfies
`A=transpose(A)` and is factored as `P^T A P=L D L^T`. Symmetry is a declared
mathematical contract, not a positive-definiteness claim. Hermitian diagonal
entries must be real. Both require finite, full canonical sorted unique CSC.
Missing diagonals mean zero. The input is not symmetrized or triangle-completed.

Pivot selection uses unrestricted 1x1/2x2 complex Bunch--Kaufman diagonal
pivoting with `(1+sqrt(17))/8` and the complex component one-norm comparison
convention described by [ZHETF2](https://www.netlib.org/lapack/complex16/zhetf2.f)
and [ZSYTF2](https://www.netlib.org/lapack/complex16/zsytf2.f). The implementation
uses original vertex labels and permutation maps, sorted sparse Schur rows and
sparse L columns. It has no global dense matrix workspace. Both classes share
one retained indexed schedule, `ldlt_z_factor.f90.inc`; the generated native
module is not hand edited. The existing real LDL reference is unchanged.

For `D=[[a,b],[d,c]]` with `d=conjg(b)` for Hermitian input or `d=b` for
transpose-symmetric input, the multiplier row `[l1,l2]*D=[u,v]` has
`det=a*c-b*d`, `l1=(c*u-d*v)/det` and `l2=(a*v-b*u)/det`.
Block and right-hand-side scales use separate exact powers of two. FortSym
emits Cartesian determinant, numerator and rank-two update primal/JVP/VJP
leaves from `gen_ldlt_z.f90`; scaled complex division reuses the retained LU
leaves. The row-solve helper is also used for the diagonal-block column solve,
with conjugated b in the Hermitian case. Backward substitution applies the
appropriate transpose or adjoint to L. After Hermitian Schur updates the
imaginary diagonal roundoff is set to zero, as in Hermitian factorization
arithmetic. Original-input residuals are checked independently by the tests.

Complex leaf pairings are ordinary real Cartesian component pairings. The
complete solve is differentiated through its equation, not through pivot
selection or the normalization schedule. This component does not itself expose
a full solve derivative; the frontend composes its retained N/T/H solves.

The methods are `factor`, `refactor`, `solve`, `solve_inplace`, `solve_block`,
`free`, `permutation`, `retained_bytes`, `peak_bytes`, `required_bytes`,
`available_bytes` and `factor_nnz`, matching the real reference's queries.
`factor(A,status,hermitian,host_limit_bytes,low_memory)` defaults to Hermitian.
`refactor` retains the previous class when the optional flag is absent and
rejects a change of class before modifying usable factors. For a matching
canonical pattern it drops old numeric storage before building new factors and
repivots the new coefficients. Numerical failure makes the factors unavailable;
the retained pattern permits recovery. Invalid input is rejected before a
refactor destroys usable factors.

All solve methods accept optional integer `operation`: 0=N (default), 1=T,
2=H. Hermitian N/H and transpose-symmetric N/T share the direct solve. The
remaining operation conjugates the right-hand side before the same factor solve
and conjugates the result. Out-of-place and inplace calls share one private
workspace core and do not alias intent-in and intent-out arguments. Block RHS
columns use that same retained vector workspace and allocate nothing.

The int64 cap covers solver-owned dynamic payload: complex sparse-row capacity,
row descriptors, pattern, permutation, D, factor/Schur data, scratch and complete
old/new row growth overlap. Zero means unlimited. LOW_MEMORY uses exact growth
instead of spare capacity; it changes storage/scheduling, not arithmetic.
Caller-owned CSC/RHS/output arrays, the fixed solver object and allocator/runtime
overhead are outside the attributed cap. `required_bytes` includes a rejected
allocation's full transient requirement; `available_bytes` reports the cap.
The reference dynamically grows sparse numeric storage and is not an
allocation-free production factorization. Solve methods allocate nothing.

Explicit status distinguishes invalid class/input/operation, singular numerical
pivots, nonfinite arithmetic, OOM, index overflow and resource-limit exhaustion.
Finite complex values whose component one-norm would overflow are rejected as
an arithmetic-range failure before the threshold sum is formed. Nonfinite
intermediates likewise fail explicitly. This is not an unconditional accuracy
or backward-stability guarantee for every finite matrix. The reference has no
parallelism, performance/frontier admission, refinement or rigorous certificate
claim. General complex LU remains a distinct external/native alternative;
it is never silently used under the LDL label.

Regenerate checked-in leaves and indexed module with:

```sh
bash gen/regenerate_ldlt_z.sh /path/to/fortsym
```

Ordinary builds require neither FortSym nor a CAS. The focused test covers both
classes, imaginary zero-diagonal 2x2 blocks at scales 1e-300 and 1e300,
nonadjacent pivots, mixed complex matrices, independent dense N/T/H residuals and
known solutions, out/inplace/block RHS, repivot/refactor/recovery, class rejection,
low-memory equivalence, exact/failing growth caps and 3000-row diagonal families.
Leaf tests independently check complex block equations, polynomial values,
Cartesian tangents and real adjoint pairing. These are behavioral oracles, not
formal proofs or deterministic certificates.
