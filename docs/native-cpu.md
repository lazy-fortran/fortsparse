# Native CPU matrix classes

Explicit `FORTSPARSE_BACKEND_NATIVE` now dispatches real SPD matrices to the
existing Cholesky engine and other classes to the following new CPU references:

| Value/class | Native factors | Current execution |
|---|---|---|
| Real general | Row-scaled threshold-pivoted sparse LU | Natural, AMD or METIS ND order, retained CSC L/U, serial factorization |
| Complex general | Row-scaled threshold-pivoted sparse LU | Natural, AMD or METIS ND order, retained CSC L/U, serial factorization |
| Complex Hermitian positive definite (`SYM_SPD`) | Sparse complex Cholesky | Natural-order simplicial reference |
| Real symmetric indefinite (`SYM_INDEFINITE` or real `SYM_HERMITIAN`) | Sparse Bunch--Kaufman LDLT | Unrestricted 1x1/2x2 pivots, serial reference |
| Complex Hermitian indefinite (`SYM_HERMITIAN`) | Sparse complex Bunch--Kaufman LDLH | Unrestricted 1x1/2x2 pivots, serial reference |
| Complex transpose-symmetric indefinite (`SYM_INDEFINITE`) | Sparse complex Bunch--Kaufman LDLT | Unrestricted 1x1/2x2 pivots, serial reference |

Real symmetric-indefinite CSC factors accept explicit AMD ordering. The initial
symmetric order is distinct from numerical Bunch--Kaufman pivot swaps and is
retained for same-pattern refactors; each refactor chooses its pivots afresh.
Natural DEFAULT and NATIVE_REFERENCE behavior is preserved; explicit AMD is
unavailable for NATIVE_REFERENCE. The common audited AMD 3.3.4 caller-workspace
path supports host caps and LOW_MEMORY. Its permutation input, retained initial
order (`4*n` bytes), graph workspace, factor growth, retained factors and optional
original-matrix refinement storage are all charged. Opaque uncapped AMD reports
complete peak memory as unknown. Ordered pre-factor planning remains unavailable.
Complex structured AMD remains unavailable.

These are usable native reference families, not production performance admission.
The previous real-SPD CPU frontier results do not extend to these classes.
HPD Cholesky and complex LDLH validate full conjugate symmetry. Real and complex
LDLT validate full ordinary symmetry. Invalid declared structures return
`INVALID_MATRIX`; HPD Cholesky also checks positive pivots.

For an assembled operator known mathematically to be real self-adjoint,
`csc_normalize_symmetric_raw(n, col_ptr, row_idx, values, rtol, &
max_entry_correction, status)` provides an explicit caller-owned normalization.
It validates finite values, exact array extents, sorted unique rows and every
reflected off-diagonal entry before changing any value. Missing zero diagonals
are allowed; missing reflected entries are rejected even for a stored zero.
With finite `rtol` in `[0,1]`, it accepts only when the maximum entry-difference
bound is at most `rtol*max(abs(original_values))`. It then copies the original
lower triangle into the upper triangle, preserving lower and diagonal entries.
The helper allocates nothing and also accepts a symmetric CSR view (where its
lower CSC triangle corresponds to the original upper CSR triangle).

On success, `max_entry_correction` is an outward one-ulp bound on the largest
absolute change to a stored entry; zero means exact symmetry. The caller must
retain this perturbation separately from solver error. Rejection leaves the
input unchanged; the correction output is only a partial diagnostic for malformed
input and may be infinity for an unrepresentable difference. This helper is
opt-in: factorization still requires exact symmetry, and a tolerance does not
establish that a general assembled operator is mathematically self-adjoint.

Existing CSC factor/refactor, in-place/vector/block solve and lifetime APIs are
preserved. LU raw CSC input is borrowed directly, including under a host cap;
no intermediate input value copy or split-complex conversion is needed.
Structured raw CSC factor/refactor uses a budgeted owned conversion for real SPD
and indefinite, complex HPD, Hermitian indefinite and transpose-symmetric systems.
It dispatches to the same Cholesky/LDL family as the CSC interface, including
zero-diagonal 2x2 LDL pivots. Host caps and LOW_MEMORY include the temporary CSC
copy beside factors, retained checker/user arrays and factor growth; caller RAW
arrays are excluded. The copy is released after factor/refactor and its peak
lifetime remains in statistics. See [memory accounting](memory.md). Exact full
canonical symmetry requirements and ordering/device capability checks are shared
with the CSC interface.
Default AUTO prefers eligible in-process external backends when built; a
native-only build uses native AUTO. This is capability selection, not additional
performance admission.

Public solver contexts own mutable resources. `sparse_copy(dst, src, status)`
and the defined assignment `dst = src` create an independent deep copy: native
CPU factors, ordering, symbolic, supernodal panel and refinement state are
duplicated into storage owned by `dst` (the symbolic-row mapping and panel
values are reallocated, never aliased), so either context can solve
(normal/transpose/adjoint, block RHS, derivatives), refactor or be destroyed
alone. `sparse_vector` slots are not copied. The destination keeps its
`host_memory_limit_bytes` and is charged for the copy (see
[memory accounting](memory.md)). SuperLU, UMFPACK IPC and MUMPS factors and
native CUDA device factors are not copied: `dst` receives the backend selection
and options with null handles, is unfactored, and the status is
FORTSPARSE_NOT_FACTORED. Assignment cannot return a status, so a failed
assignment leaves `dst` a valid unfactored context; use `sparse_copy` when the
outcome matters. Elemental array assignment is supported. Self-assignment is
a no-op without aliasing or double release, but gfortran's self-assignment
temporary deep-copies backend arrays that it does not release; avoid `a = a`.
Assigning a derived type that merely contains a `sparse_solver_t` component is
not supported: gfortran implements that through intrinsic component copies,
which can alias backend storage. Copy such components explicitly.

The new references enforce `host_memory_limit_bytes` across their owned arrays,
including factor growth replacement overlap, and support RAM-only LOW_MEMORY.
LOW_MEMORY LU grows factor buffers by at most 1/16 of their old capacity
beyond the requested entries, amortizing prefix copies; buffers below 16 entries
still grow exactly. At each growth to required size R, capacity is at most
R + floor(R/16), compared with AUTO's doubling policy. Every allocation checks
old/new overlap against a positive cap and first falls back to an exact request.
If spare capacity created by the current attempt still prevents growth, LOW
releases its already invalid numerical factors and retries once with exact
growth, retaining the pattern/order workspace. It does not retry singular or
OS allocation failures. Peak accounting includes both sequential attempts;
successful retry reports pattern reuse false. Tightly capped LOW can therefore
trade more copying for minimum padding. This retry adds no entry-cap bypass:
the existing low-level LU entry policy may release storage when it already
exceeds a newly requested cap. Public adapter reservations for ordering,
refinement and user vectors remain part of the budget before LU runs.
On Linux, real LOW LU advises strictly interior pages of dead factor buffers
at least 1 MiB long with `MADV_DONTNEED` immediately before deallocation. Both
replacement prefix copies finish first. Caller arrays, live factors, correction
workspace, allocator metadata and partial boundary pages are untouched.
Allocatable ownership, capacity accounting and old/new overlap remain unchanged.
Advice failure is nonfatal, and other platforms use ordinary deallocation.
Complex LU and AUTO execution do not issue this advice. This is an allocator
retention tradeoff with syscall/page-fault cost, not a hard RSS limit or a
performance-admission claim. Automatic low-level scope deallocation does not
issue the extra hint; explicit factor free and public backend finalization do.
The [two-workload current-source checkpoint](../bench/results/2026-10-16-low-dead-pages/report.md)
observes 12–14% lower RSS at about 19% higher completed-process CPU. This is a
LOW resident-RAM/CPU tradeoff; it strengthens neither the payload cap nor any
external-solver frontier claim. AUTO and complex LU retain ordinary free.
HPD already uses compact exact structural storage. Separate contexts own their
mutable work independently. Strict-budget and LOW_MEMORY execution use natural
ordering, or the audited bounded AMD path for general LU; these matrix classes
currently execute on CPU. LU, HPD and
real/complex LDL support `solver%refine=.true.` with
retained original-matrix residuals, a 32-epsilon normwise backward-error target
and at most three corrections; corrections also continue while the componentwise
backward error exceeds eps and still halves (LAPACK xGERFS style). Statistics
report correction steps and the final normwise and componentwise backward errors.
See [native refinement](native-refinement.md) for failure,
allocation and refactor-cap semantics. Production blocking/parallel
factorization and conditioning diagnostics remain planned.

`pivot_tolerance` selects LU's finite threshold in [0,1]; a negative value uses
the native default 1.0 (maximum-candidate row pivoting after scaling). A numerical
refactor always selects pivots anew. No derivative is taken through pivot choices.

Opt-in [native structural row matching](native-matching.md) supplies a nonzero
preferred diagonal for real/complex GENERAL LU without changing the original
matrix. It recomputes matching on numeric refactors and charges its graph/maps
to the host payload cap. The default policy is unchanged.

## Sylvester inertia

`sparse_inertia(solver, n_negative, n_zero, n_positive, reliability, status)`
and the `inertia_*` fields of `sparse_factor_stats_t` report the inertia of
the last native symmetric factorization from its actual pivots: 1x1 pivot
signs, and for each 2x2 block the signs of both eigenvalues (formed after
scaling by the largest block entry, the smaller one as `det/larger`, so
neither cancels). By Sylvester's law this is the inertia of the matrix that
was exactly factored. For a pencil `K - sigma*M` with SPD `M`, `n_negative`
counts generalized eigenvalues below `sigma`. Same-pattern `sparse_refactor`
over many shifts reuses the pattern and initial order and repivots afresh.

`reliability` (`FORTSPARSE_INERTIA_*`) is `RELIABLE` when every nonzero pivot
eigenvalue exceeds its local threshold `n*eps*g_k`, where `g_k` is the largest
magnitude of an input entry or of any Schur-complement entry formed in the
pivot's row(s) (a growth-scaled row norm; exactly decoupled small-scale rows,
such as zero stiffness rows of a shifted pencil, keep small thresholds).
Otherwise it is `UNCERTAIN` (counts still returned; `inertia_uncertain`
eigenvalues lie at or below their thresholds). `inertia_min_pivot` is the
smallest pivot eigenvalue magnitude and `inertia_tolerance` its threshold. Exactly zero pivot columns
are eliminated as zero 1x1 pivots: the factorization returns `SINGULAR`,
solves are refused, and inertia is still returned with `n_zero > 0` and
`reliability = SINGULAR`. The threshold is a backward-error heuristic, not a
rigorous enclosure: `RELIABLE` counts belong to the factored floating matrix
and agree with `A` only when its eigenvalues are separated from zero by more
than the backward error.

Supported families: real LDLT, complex Hermitian LDLH (real pivots) and
Cholesky (real SPD and complex HPD report `(0, 0, n)` on success). LU
(including the supernodal LU, whose static pivot perturbation can change the
inertia), complex transpose-symmetric LDLT and external backends return
`BACKEND_UNAVAILABLE` with `reliability = UNAVAILABLE`. The inertia path never
perturbs pivots.

## Optional general LU ordering

For real and complex general matrices, production native LU accepts
`solver%options%ordering=FORTSPARSE_ORDER_AMD` when SuiteSparse AMD is built.
AMD accepts the canonical unsymmetric CSC graph and orders its structural
`A+A^T`; numerical factorization borrows the original CSC values and applies the
column permutation through its retained sparse factors.
`FORTSPARSE_ORDER_ND` alternatively uses optional METIS nested dissection on the
undirected graph union of `A` and its transpose. Self edges and repeated mirrors
are removed before ordering. With structural matching, either ordering consumes
the matched `P*A` graph while numerical values stay in the original CSC.
DEFAULT remains natural ordering. NATIVE_REFERENCE remains natural and rejects
explicit graph orderings. Missing graph libraries return BACKEND_UNAVAILABLE.
`native_ordering_available(ordering)` queries the built native graph capability
without allocating or factoring. Matrix-class, device and memory-policy support
are still checked by factorization; this query is not a strict-cap guarantee.

Value-only refactors reuse the retained column order when the canonical pattern
and requested ordering strategy match; numerical row pivots are selected anew.
A changed pattern recomputes the order, and a changed ordering strategy invalidates
pattern-reuse reporting. Normal, transpose and adjoint solves, refinement and
implicit solve derivatives act on the original matrix.

Uncapped AUTO-memory AMD uses the existing opaque `amd_order` path: actual
ordering and retained owned payload are reported, but complete host peak/factor
memory remains unknown (`-1`) and `strict_memory_supported` is false. Positive
host caps and LOW_MEMORY instead use the audited AMD 3.3.4 `amd_2` caller-workspace
path, including old-factor, user-vector, refinement and permutation-input overlap
in cap admission and phase peak statistics. An unavailable audited ABI fails
explicitly. See [bounded AMD](native-amd-bounded.md) for exact workspace counts,
release trust, reuse and verification limits. Use measured process RSS for process
memory comparisons; these changes make no performance-admission claim.

METIS owns opaque internal workspace: ND always reports unknown complete peak
(`-1`) and `strict_memory_supported=.false.`. Positive host caps or LOW_MEMORY
with ND fail before calling METIS; use natural or bounded AMD for those contracts.
This does not extend ND to the complex structured reference families.

## Opt-in supernodal real LU

`solver%options%lu_method=FORTSPARSE_LU_SUPERNODAL` selects the static-structure
supernodal LU ([design](lu-supernodal.md)) for real GENERAL matrices on the
explicit `FORTSPARSE_BACKEND_NATIVE` CPU backend. The default stays
`FORTSPARSE_LU_SCALAR`; `FORTSPARSE_LU_AUTO` behaves as SCALAR until its
thresholds are calibrated. The family reports `native_factor_family='LU_SUPERNODAL'`.

- Analysis: maximum-product matching with dual scaling, then AMD (bounded
  caller-workspace AMD when available) or METIS ND on the matched symmetrized
  pattern, an etree postorder, and the shared supernode partition and schedule.
  DEFAULT ordering selects METIS ND when no host cap is set, LOW_MEMORY is not
  requested and METIS is built; otherwise AMD. An explicit AMD or ND request
  wins. `ordering_used` reports the ordering actually used. `structural_matching` does not apply; the matching is always
  weighted (`matching_used=FORTSPARSE_MATCHING_WEIGHTED`).
- Factor: serial left-looking panels, partial pivoting restricted to each
  supernode diagonal block, and pivots below `pivot_perturbation*||A~||_inf`
  (default 1e-13) replaced by `sign(p)` times that threshold. Statistics report
  `perturbed_pivots`, `block_pivot_swaps`, `pivot_growth`, `min_pivot` and
  `scaling_log10_range`.
- Gates: a non-finite factor entry, more than `max(0,n/1000)` perturbed pivots
  or growth above `1/(64*sqrt(eps))` fall back to the scalar LU (threshold 0.1,
  structural matching, AMD unless ND was requested). With `verify_factor`
  (default true) and always after perturbations, one solve against the original
  matrix is refined at most three times to the `32*eps` backward-error gate;
  failure falls back with reason REFINEMENT. A perturbed factor retains the
  original matrix for refinement even when `refine` is false; if that does not
  fit the cap it falls back instead. Structurally singular matrices go to the
  scalar LU (reason SINGULAR), which reports the singularity. The fallback runs
  its own preflight under the same cap: a refusal returns RESOURCE_LIMIT and keeps
  `lu_fallback_reason`. A successful fallback returns OK with `lu_fallback=.true.`,
  family `LU_REAL` and the supernodal attempt's pivot counters. A refined user solve
  that misses the gate refactors once with the fallback when no host cap is set.
- Refactor: a same-pattern refactor reuses matching, scaling, order, partition
  and panels (`reused=.true.`). An entry outside the static structure or a
  failed gate repeats the full analysis. NORMAL, TRANSPOSE/ADJOINT and block solves,
  refinement, implicit JVP/VJP and context copies use the retained factors.
- Unsupported: complex values (BACKEND_UNAVAILABLE), structured symmetries,
  NATIVE_REFERENCE, devices and OpenMP execution (one active thread). ND is
  refused under a host cap or LOW_MEMORY. LOW_MEMORY currently changes only
  the matching cost cache and the AMD workspace.

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
