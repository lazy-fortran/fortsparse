# Sparse LU source and contract

`gen_lu.f90` retains Cartesian complex `a-b*c` and division by a nonzero
complex denominator, emitting primal/JVP/VJP with the ordinary real pairing
of real and imaginary components. The real path reuses the retained generated
SPD divide/update leaves; it does not assume positive pivots. Complex division
normalizes both operands by powers of two, evaluates the generated normalized
quotient, then restores the exponent difference. The normalization schedule is
in `lu_complex_schedule.f90.inc`; it avoids squaring the original denominator.

`lu_factor.f90.inc` is the shared indexed real/complex LU schedule, emitted as
two concrete modules in `src/native/fortsparse_native_lu.f90`. Ordinary builds
need neither FortSym nor generation. Regenerate from the repository root:

```sh
bash gen/regenerate_lu.sh /path/to/fortsym [GEN_BUILD [FORTSYM_BUILD]]
git diff --exit-code -- src/native/fortsparse_native_lu.f90 src/generated/fortsparse_lu_leaves.f90
```

The script checks the retained FortSym revision. Change the definitions or
indexed template and regenerate; never edit generated output.

Modules `fortsparse_native_lu` and `fortsparse_native_lu_z` export `lu_factor_t`
and `lu_factor_z_t`. Both expose:

- `factor(a,status,column_permutation,pivot_threshold,memory_limit_bytes,low_memory)`;
- `factor_raw(n,col_ptr,row_idx,values,status,...)`, with the same optional controls;
- generic `solve(rhs,status,operation)` for one vector or a block of RHS;
- `free()` and `memory_bytes()`;
- `n`, `factor_nnz`, `peak_memory_bytes`, `last_reused_pattern`,
  `failed_required_bytes`, and `failed_available_bytes` diagnostics.

Operation codes are `LU_NORMAL=0`, `LU_TRANSPOSE=1`, and `LU_ADJOINT=2`.
For real values the last two operations coincide. Operations and dimensions are
validated even for an empty RHS block. Objects own mutable retained workspaces;
one object must not be used concurrently. Distinct objects share no state.

Input is square one-based CSC with sorted unique row indices and finite
components. Raw input is borrowed during the call; no input values are copied
into another CSC container. Pattern indices are retained as owned metadata.
The optional column permutation is a one-based new-to-old bijection. Numeric
row pivots are always selected anew. `last_reused_pattern` means identical
input column pointers and row indices, not unchanged pivots or numeric factors.

Each row is divided by its largest component magnitude. Real pivot magnitude
is absolute value; complex pivot magnitude is the maximum absolute real or
imaginary component. Threshold defaults to 1 and must be finite in `[0,1]`.
Among remaining rows, the original diagonal row for the current permuted
column is preferred only when nonzero and at least `threshold*maximum`.
Otherwise the largest candidate is selected. Sparse iterative DFS reach over
existing L columns precedes each numerical left-looking elimination. L has
implicit unit diagonal; U stores its diagonal. No dense n-by-n storage exists.

Normal, transpose and adjoint solves apply row scaling and both permutations
consistently. Solve allocations are retained O(n) workspace. Refactorization
recomputes reach/pivots and reuses available capacity; changed pivots/fill can
require bounded storage growth. All old factors are invalidated at factor entry
and remain unusable after any failed factorization. Singular pivots, nonfinite
elimination and detected scaling/multiplier underflow return `SINGULAR` with a
diagnostic. Failed outputs are not solutions. Invalid inputs/options, allocation
failure and a memory cap return their separate named status codes.

The byte cap is an integer(i8), zero meaning unlimited; a negative cap is an
invalid option. `memory_bytes()` counts exact allocated-array payload, including
unused factor capacity, indices, permutations, retained input pattern, numeric
work and row scaling. It excludes caller input/RHS, derived-object descriptors
and allocator overhead. Peak payload includes both old and replacement arrays
during growth and partial successful allocations before an allocation failure.
The cap is checked before allocations, including this copy peak. On a budget
failure, `failed_required_bytes` records the additional allocation payload
requested and `failed_available_bytes` records `max(0, cap-current payload)`.
Both diagnostics reset at each factor attempt and on free. Default growth
is geometric; low-memory mode grows only to the required factor count. It can
trade more copying for less unused capacity. Process RSS remains a separate
measurement. `factor_nnz` counts L/U including both diagonals, even though L's
diagonal is implicit and consumes no value allocation.

The focused fixtures use independent original-matrix known solutions for
N/T/H, numerical row-pivot adversaries, same-pattern changed values, failure
recovery, direct raw input validation, sparse cycles with linear storage,
strongly unequal row scales and generated arithmetic derivative pairings.
These establish numerical behavior; they are not forward-error certificates,
formal proofs, threaded/GPU implementations or production Pareto admission.
