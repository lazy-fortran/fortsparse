# Native real CSR ILUTP preconditioner

`sparse_ilutp_t` owns one real incomplete factorization with dual dropping
and column pivoting. It is an independently written native preconditioner,
not a new direct-solver backend. It invokes no external sparse factorization.
Krylov orchestration remains in FortNum. SPEC migration, application accuracy,
representative consumer extraction, threading/scaling and performance admission
remain unimplemented/unverified by this slice.

```fortran
use fortsparse, only: sparse_ilutp_t, sparse_ilutp_stats_t
! factor%setup(n, row_ptr, columns, values, lower_fill, upper_fill, &
!              drop_tolerance, pivot_tolerance, status, host_memory_limit_bytes)
! factor%apply(b, x, status)
! info = factor%stats()
! call factor%free()
```

## Frozen numerical/storage contract

- Input is square one-based CSR, default-integer indices (I4 build), binary64
  real values. Row pointers must span the supplied arrays exactly and be
  monotone. Column order within a row is unrestricted. Nonfinite values,
  invalid indices/shapes and negative fill or memory limits are rejected.
  Drop and pivot tolerances must be finite and in [0,1]. Input is borrowed
  read-only; setup does not permute or overwrite the caller's arrays.
- Duplicates are accumulated in their original input order. The original
  stored-row mean `sum(abs(input entries))/stored count` is evaluated BEFORE
  coalescing and includes stored zeros. A scaled reduction avoids overflow of
  the magnitude sum. Thus adding a stored zero or splitting an entry can
  change thresholding even when the represented matrix is unchanged. No
  claim of representation-invariant thresholding or bitwise SPARSKIT parity
  is made. Empty/all-zero stored input rows fail; cancellation between
  nonzero duplicates can instead lead to a reported regularized pivot.
- Lower elimination proceeds in increasing current column order. A
  multiplier with magnitude <= `drop_tolerance` is discarded before updating
  the row. Surviving lower multipliers are later truncated by `lower_fill`.
  Updated upper entries with magnitude <= `drop_tolerance*original_row_mean`
  are discarded, then at most `upper_fill` off-diagonals are retained.
  The diagonal is separate and does not consume either fill cap.
- Both fill caps count off-diagonal entries per row. Largest magnitude wins;
  equal magnitudes prefer the smaller current elimination column. Selection
  uses bounded heaps. A candidate is considered for column pivoting only
  after upper dropping/truncation. It must satisfy
  `abs(candidate)*pivot_tolerance > abs(diagonal)`; the largest qualifying
  candidate wins, with the same column tie rule. Pivot tolerance zero
  disables swaps. There are no row swaps or restricted pivot blocks.
- Column maps are context-owned: `perm(new)=old`. Prior upper entries retain
  original-column labels during setup; their final elimination labels are
  resolved after all swaps. With Q's columns equal to `e_perm(new)`, the
  approximation is `A*Q ~= L*U`; normal apply returns `Q*U^-1*L^-1*b`.
- An exact zero pivot is replaced by `(1e-4+drop_tolerance)*original_row_mean`
  and increments `regularized_pivots`. A nonfinite/zero replacement or
  nonfinite elimination fails. Successful setup means an applicable
  preconditioner, not a nonsingular original A, an exact factorization,
  a small original-system residual or a rigorous certificate.
- Apply accepts distinct, correctly sized real input/output vectors and
  supports normal application only. It uses the retained scratch vector;
  there is no allocation, backend dispatch or external numerical call.
  Nonfinite input/arithmetic returns explicit failure. An apply failure does
  not destroy a valid factor; the next apply overwrites its scratch state.
  Output on failure is unspecified. Concurrent use of one context is
  unsupported; separate contexts have no shared mutable state.

Drop/pivot decisions are discrete and are not differentiated. Scalar quotient
and subtract-product operations reuse the existing FortSym-generated
`spd_divide` and `spd_update` definitions, including their retained JVP/VJP
leaves in [gen_spd_leaves.f90](../gen/gen_spd_leaves.f90). Their algebra does
not require an SPD matrix. Indexed elimination, heaps, permutations and
storage/lifetime orchestration are independently handwritten. This does not
add derivatives of Krylov iterations or a claim about nonsmooth choices;
scientific solve sensitivities still differentiate the original equation.

## Payload cap and recovery

Setup destructively frees all old arrays BEFORE validation or new allocation.
Failure leaves the context unfactored and owning zero array bytes; setup may
subsequently recover. No old/new factor overlap is hidden or charged twice.
The caller's old factor is not retained if a replacement fails.

Let `kL=min(lower_fill,n-1)`, `kU=min(upper_fill,n-1)` and
`C(k)=k*n-k*(k+1)/2`. Flat arrays reserve capacities `CL=C(kL)`, `CU=C(kU)`.
Capacity products use checked-range 64-bit arithmetic; capacities whose
terminal pointer cannot fit the default integer fail before allocation.
For integer/real element byte widths I/R, the exact requested setup payload is

```text
I*(CL+CU+7*n+2) + R*(CL+CU+2*n).
```

It includes both factor value/index capacities, two row-pointer arrays,
two permutation arrays, the diagonal, retained apply scratch, and all three
integer setup scratch arrays. Setup releases the latter after success, so
retained payload is setup payload minus `3*I*n`. Unused reserved factor slots
remain allocated and counted. There are no per-row allocations, dynamic growth
or hidden automatic arrays proportional to n.

A positive `host_memory_limit_bytes` rejects the complete setup payload before
any allocation; zero means unlimited. Partial allocation failure also frees
every successfully allocated component. `stats()` reports current allocated
payload, observed peak allocated payload, requested setup bytes, required bytes
on failure, fill, pivot replacements, swaps and duplicate/drop diagnostics.
Fixed Fortran objects/descriptors, scalar locals, allocator overhead, caller
CSR/RHS/output storage and process/runtime RSS are excluded. This is an array
payload cap, consistent with [memory.md](memory.md), not a portable RSS cap.
The apply workspace is context-owned and INCLUDED in retained payload.

## Actual SPEC seam and limits

Source inspection used SPEC revision
`b23a2408fc3f1dfbc14896e16d73a57f7d1d9e88`:

- `src/spsmat.f90:19–36,550–570` describes the selected harmonic/constraint
  approximation and builds one-based CSR. `src/mp00ac.f90:428–438` forms
  `dMAS-mu*dMDS` and invokes ILUTP. `:273–285` chooses `maxfil` and storage.
- SPEC passes `permtol=0.1`, `mbloc=n`, `epsILU` (default 1e-12), `lfil=maxfil`.
  Its actual `src/ilut.f:683–725` U loop stores at most `maxfil-1`
  off-diagonals, despite header prose. Native explicit caps can express
  L=`maxfil`, U=`max(0,maxfil-1)` without changing generic cap meanings.
  This cap correspondence does not promise identical selected entries,
  rounding or Krylov histories.
- `src/ilut.f:751–753` comments out input-column mutation; native input is
  likewise read-only. `src/mp00ac.f90:900–914` applies the inverse column map
  after triangular solves; native apply owns that correction.
- SPEC constructs factors for its primal call and reuses them for derivative
  RHS (`src/mp00ac.f90:435–459`). Its GMRES is LEFT-preconditioned, restart 5,
  with a preconditioned-residual test (`:807–839`, `src/iters.f:1988–2017`).
  FortNum's existing FGMRES is RIGHT-preconditioned and its GMRES is
  unpreconditioned; changing this is a separate consumer/Krylov contract.

No SPEC executable or physical case was run or migrated. The public tiny
tests are deliberately structural, not claimed to be extracted SPEC data.
Application acceptance needs an actual retained consumer operator/RHS,
original-equation residual and requested observable checks. No private
numeric data belongs in source history.

## Algorithm provenance and verification

Yousef Saad, *Iterative Methods for Sparse Linear Systems*, second edition,
SIAM (2003), §§10.4.1 and 10.4.4 (printed pp.321–322 and 327–329), supplies
the public mathematical ILUT/column-pivoting description:
[author-hosted PDF](https://www-users.cse.umn.edu/~saad/IterMethBook_2ndEd.pdf).
The cited paper is Saad (1994), *ILUT: A dual threshold incomplete LU
factorization*, DOI [10.1002/nla.1680010405](https://doi.org/10.1002/nla.1680010405).
The book is retained locally outside Git; the paper was used only through
bibliographic metadata, not an unavailable full-text claim.

Viewed SPEC source includes Saad-attributed 1993/1996 SPARSKIT routines. SPEC's
root MIT license does not by itself resolve original third-party provenance.
No SPARSKIT implementation text, comments, quick-split routine or factor
storage code is copied into FortSparse. The new source uses its own bounded
heaps, split flat CSR factors and context lifetime. Existing FortSparse
generated arithmetic remains MIT. Provenance archives, hashes, reading depth
and library cards are retained in `~/code/sparse-dev/provenance/ilutp/` and
`~/code/sparse-dev/library_cards/ilutp/`.

Affected checks are `test_fortsparse_ilutp`, `test_ilutp_allocations` (supported
Linux link wrappers), and the reused `test_fortsparse_spd_leaves`. They cover
exact dyadic constraint-block inverses, column swaps, unordered CSR,
duplicates/stored zeros, drop equalities, independent fill caps/ties,
regularization, pointer/cap failures, recovery, every partial setup allocation
failure and allocation-free repeated normal apply. These are correctness and
memory-lifetime checks, not performance admission or a forward certificate.
