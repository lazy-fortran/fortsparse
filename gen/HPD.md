# Native complex HPD reference

`gen_hpd_leaves.f90` retains the continuous arithmetic for natural-order sparse
Cholesky, `A = L L^H`. It represents complex scalars by real and imaginary
components and derives primal/JVP/VJP together with FortSym. The component
pairing is the real Euclidean pairing, equivalently `Re(conjg(u)*v)` for complex
scalars. The leaves are positive square root, complex division by a real
positive diagonal, ordinary and conjugated multiply/subtract, and subtraction
of a squared complex norm. Their domains are checked by sparse orchestration.

Regenerate with `gen/regenerate_hpd.sh /path/to/fortsym`; the script checks the
pinned FortSym revision. Ordinary builds use the checked-in source and require
no symbolic package. `hpd_schedule.f90.inc` only packs/unpacks complex values
around generated component contractions.

`fortsparse_native_hpd` accepts full canonical Hermitian CSC, exact conjugate
pairs, real diagonals and finite values. Symbolic reach creates CSC fill for L;
numerical factorization retains that layout. Zero, negative and nonfinite
pivots return `FORTSPARSE_NOT_POSITIVE_DEFINITE`. Failed refactors invalidate
the factor; another valid same-pattern refactor can recover without allocation.
Normal (`operation=0`) and adjoint (`2`) solves coincide. A transpose solve
(`1`) conjugates the RHS and result around the normal solve. `solve_inplace`
provides an allocation-free single-vector API with readiness, shape, operation
and finite-input/output guards. Out-of-place vector solves use that same path.
No input matrix
or dense matrix/Jacobian is retained.

The owning payload is exactly `20*nnz(L) + 44*n + 4*nnz(A) + 8` bytes for the
supported 32-bit-index/complex128 build: complex factor values, sparse indices,
original pattern, reach metadata and one complex work vector. The implementation
uses `storage_size` when charging actual arrays. Every payload is reserved before
allocation through `memory_budget_t`. Current/peak payload excludes descriptors,
allocator/runtime overhead, process RSS and caller A/B/X buffers. Rejected
reservations do not inflate measured live/peak payload. A new factorization frees
old storage before allocating new storage; same-pattern refactor and solves
allocate no payload. A cap below an existing retained payload is rejected. Changed-pattern refactors
explicitly forward the selected cap and low-memory policy to the new analysis.
`memory_failure_bytes` exposes the last denied reservation and its available
capacity; tightening below live storage reports the retained total and proposed
cap instead.

Both AUTO/default and LOW_MEMORY use this same minimal reference workspace;
LOW_MEMORY makes no accuracy or precision change. Factorization is a serial
simplicial reference with SIMD update loops. Independent block RHS columns can
use OpenMP. There is no throughput, scaling, ordering or Pareto-admission claim.
Focused tests use independent known solutions, transpose/adjoint identities,
finite-difference and component-pairing derivative checks, invalid-input/pivot
failures, exact-cap boundaries and a 2000-row tridiagonal linear-storage case.
