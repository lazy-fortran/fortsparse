# Native CPU matrix classes

Explicit `FORTSPARSE_BACKEND_NATIVE` now dispatches real SPD matrices to the
existing Cholesky engine and other classes to the following new CPU references:

| Value/class | Native factors | Current execution |
|---|---|---|
| Real general | Row-scaled threshold-pivoted sparse LU | Default AMD (see below), natural or ND order, retained CSC L/U, serial factorization |
| Complex general | Row-scaled threshold-pivoted sparse LU | Default AMD (see below), natural or ND order, retained CSC L/U, serial factorization |
| Complex Hermitian positive definite (`SYM_SPD`) | Supernodal Cholesky of the real embedding | Uncapped CPU; under a host cap, LOW_MEMORY or NATIVE_REFERENCE the simplicial complex Cholesky |
| Real symmetric indefinite (`SYM_INDEFINITE` or real `SYM_HERMITIAN`) | [Multifrontal LDLT](#multifrontal-ldlt) | Uncapped CPU, serial; under a host cap, LOW_MEMORY or NATIVE_REFERENCE the scalar Bunch--Kaufman LDLT |
| Complex Hermitian indefinite (`SYM_HERMITIAN`) | Multifrontal LDLT of the real embedding | Uncapped CPU; under a host cap, LOW_MEMORY or NATIVE_REFERENCE the scalar Bunch--Kaufman LDLH |
| Complex transpose-symmetric indefinite (`SYM_INDEFINITE`) | Sparse complex Bunch--Kaufman LDLT | Default AMD; restricted 1x1/2x2 pivots with delays, serial reference |

The symmetric and Hermitian CSC factors, real and complex, accept explicit AMD
or ND ordering. The real embedding of a complex Hermitian `A = Ar + i Ai`
is `E = [[Ar, -Ai], [Ai, Ar]]` with interleaved unknowns: every eigenvalue of
`A` appears twice in `E`, so the reported inertia is half that of `E`, and the
statistics name the families `LDLH_EMBEDDED` and `HPD_EMBEDDED`.

## Default ordering

`FORTSPARSE_ORDER_DEFAULT` selects AMD when the matrix has at least 4096
columns and an entry more than 12 positions from the diagonal, unless the
pattern fills its own envelope (every column is one run of rows around the
diagonal and the runs never move back up, as in dense bands and
block-tridiagonal matrices with dense blocks), where natural order creates no
fill. Otherwise it keeps natural order. This holds for real SPD and complex HPD Cholesky, real and complex
general LU, and the scalar symmetric indefinite LDL^T/LDL^H factors. The
multifrontal LDL^T orders DEFAULT as AMD.
`FORTSPARSE_ORDER_NATURAL` keeps the input order; `ordering_used` reports
natural order as 0. AMD needs no external library: the project-owned port of
SuiteSparse AMD (`fortsparse_amd.f90`, the algorithm of Amestoy, Davis and
Duff, 1996) orders as `amd_order` does, also under `fpm`. General LU takes the
accounted caller-workspace AMD for the default, so exact host peaks and caps
keep working; so do the HPD and indefinite factors and real SPD under a host
cap, LOW_MEMORY or refinement, where `sparse_analyze` plans the AMD layout that
the factorization will use. Only NATIVE_REFERENCE keeps natural order by default.

`FORTSPARSE_ORDER_ND` needs no external library either: `fortsparse_nd` is a
Fortran port of METIS 5.2.1 `METIS_NodeND` with default options and seed 0
(multilevel node bisection, minimum-degree leaves). Its random numbers, sorts
and priority queue are those of METIS's GKlib, with glibc's `rand`, so it
returns the orderings of METIS on Linux on every platform;
`test_fortsparse_nd` compares them with orderings recorded from METIS. Its run
time equals METIS's within a few percent on grid and random graphs of
10^3 to 2*10^5 vertices.

### Multifrontal LDLT

The real symmetric indefinite factor `P^T A P = L D L^T` is multifrontal
(Duff and Reid 1983, as MA57 and MUMPS). Each supernode of the elimination
tree owns a dense front: its fully summed variables, which are the variables
its children delayed plus its own columns, followed by the rows of its
structure. Children's contribution blocks are added into the parent front
(extend-add). A fully summed variable `k` is accepted as a 1x1 pivot if
`|a_kk| >= u max_i |a_ik|` over all remaining rows of the front, or as a 2x2
pivot with the largest fully summed entry `r` of its column if
`u |P^-1| (max_i |a_ik|, max_i |a_ir|)^T <= (1, 1)^T`, with `u = 0.01`. The 2x2
block is scaled by its largest entry before the test and the inverse, so no
product over- or underflows. Variables that fail stay in the front, are
updated, and move to the parent front as delayed pivots. At a root, where
nothing can be delayed, the unrestricted Bunch--Kaufman choice is taken, so a
pivot always exists. An exactly zero column is a zero 1x1 pivot.

Pivots are sought in windows of 32 fully summed columns that are updated by
rank-1/2 steps. A candidate is tested once per window. The rest of the front
is updated by BLAS-3 after each window (fortnum `dgemm`). The front is
assembled and factored in place at the end of the factor storage, so only the
contribution block is copied, to a stack.

`DEFAULT` and `AMD` order with AMD, `ND` with nested dissection, and `NATURAL`
keeps the input order. A variable whose
diagonal fails the pivot test is first matched with its largest unmatched
neighbour, and the pair is one vertex of the ordered graph, the member with
the larger diagonal first (compressed ordering, Duff and Pralet 2005). A
saddle-point constraint row thus meets its partner in the same front instead
of being delayed. The elimination tree is postordered and its fundamental
supernodes are amalgamated with the CHOLMOD rules, keeping matched pairs
together. `sparse_refactor` with the same pattern and ordering option reuses
the ordering and the symbolic structure and chooses its pivots afresh.

On a 2947-unknown free-boundary MHD saddle-point matrix with 111 746 entries
and 1013 zero diagonals (one thread), the factor holds 2.5e5 entries with
about 180 delayed pivots, as MUMPS 5.8 with its compressed AMD ordering does
(2.3e5 entries, 171 delays):

| Solver | Analysis + factorization | Refactorization |
|---|---|---|
| Multifrontal LDLT | 8.4-9.0 ms | 3.1-3.4 ms |
| MUMPS LDLT | 9.2-10.1 ms | 3.7-3.9 ms |
| Scalar Bunch--Kaufman LDLT | 1060 ms | - |

The inertia and its reliability follow [Sylvester inertia](#sylvester-inertia),
with `g_k` the largest input entry of the pivot row plus `sum_j l_kj^2 |d_j|`
over the earlier pivots that updated the row, or the largest entry of the
pivot column if larger. A 2x2 pivot `D_j = [d1 e; e d2]` contributes
`|l_j|^T (|D_j| + |e| I) |l_j|`, a positive semidefinite majorant of `|D_j|`.
The statistics report `native_factor_family = 'LDLT_MULTIFRONT'`,
`factor_nnz` (stored entries, with amalgamation zeros), `symbolic_nnz`
(`nnz(L)` of the order), `block_pivot_swaps` (2x2 pivots) and the retained
payload as `host_current_bytes`. The AMD and BLAS workspaces are opaque, so
`host_peak_bytes` is unknown and no strict cap applies.

### Scalar LDLT

Under a host cap, LOW_MEMORY or NATIVE_REFERENCE, real and complex symmetric
indefinite matrices use a sparse Bunch--Kaufman elimination whose memory is
fully accounted. Within a fill-reducing order the next vertex `k` is accepted
as a 1x1 pivot if `|a_kk| >= u max_i |a_ik|`, or else as a 2x2 pivot with its
first remaining neighbour `p` in the order (its elimination-tree parent) if
every multiplier stays below `1/u`. Otherwise `k` is delayed to the position
right after `p`. A vertex that comes up again before any elimination took
place takes the ordinary unrestricted Bunch--Kaufman pivot, so a pivot always
exists. The threshold is `u = 0.1`. Natural order keeps unrestricted
Bunch--Kaufman with position swaps; see
[the measurements](../bench/results/2026-10-30-symmetric-ordering/README.md).
OpenMP workers merge disjoint Schur rows (`ldlt_active_threads`); the merge
workspace costs `4*n + 12*n*threads` bytes and is charged, and LOW_MEMORY
retains entry-wise updates. The initial symmetric order is retained for
same-pattern refactors; each refactor chooses its pivots afresh.
explicit AMD is unavailable for NATIVE_REFERENCE. The common caller-workspace
AMD path supports host caps and LOW_MEMORY for the real and complex LDL and
the HPD factors. Its permutation input, retained initial order (`4*n` bytes),
delay marks (`4*n` bytes during factorization), graph workspace, factor
growth, retained factors and optional original-matrix refinement storage are
all charged; ordered HPD factors retain the permuted pattern plus a `4*n`
order and a `4*nnz` source map. Opaque uncapped AMD reports complete peak
memory as unknown. Ordered pre-factor planning remains unavailable for these
families.

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
a no-op when the compiler passes `test_fortsparse_compiler_self_assignment`.
The live-backend finalization defect was reproduced with GNU Fortran 11 and
reported for GNU Fortran 12.3; GNU Fortran 16.2 passed the probe.
CMake rejects an affected compiler during
configuration. The fpm build provides the same test, and defined assignment
checks the compiler before returning from an assignment with aliased factors,
stopping before the compiler can release the user's factors. An explicit
`call sparse_copy(a, a, status)` avoids the assignment temporary.
The lifetime probe does not check leaks from temporary backend-array copies;
prefer the explicit copy when copying an object onto itself.
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
ordering or the audited bounded AMD path; these matrix classes
currently execute on CPU. LU, HPD and
real/complex LDL support `solver%refine=.true.` with
retained original-matrix residuals, a 32-epsilon normwise backward-error target
and at most three corrections; corrections also continue while the componentwise
backward error exceeds eps and still halves (LAPACK xGERFS style). Statistics
report correction steps and the final normwise and componentwise backward errors.
See [native refinement](native-refinement.md) for failure,
allocation and refactor-cap semantics. Production blocking/parallel
factorization and conditioning diagnostics remain planned.

`pivot_tolerance` selects LU's finite threshold in [0,1]: the preferred (diagonal
or matched) row is kept while its row-scaled magnitude is at least the threshold
times the column maximum. A negative value selects the native default: 1.0
(maximum-candidate row pivoting after scaling) without refinement and 0.01 with
`refine=.true.`. Statistics report the resolved `effective_pivot_tolerance`.
A numerical refactor always selects pivots anew. No derivative is taken through
pivot choices.

The refined default follows the KAMEL QL-Balance time-step system (n=43804,
natural order, entries 5.8e-41..6.5e23, solution components 1e-10..1e14). Rows of one field
couple to tiny-solution fields with coefficients up to 1e24; these dominate the
row maxima but not `|A||x|`, so maximum pivoting selects off-diagonal rows. With
threshold 1.0 the unrefined componentwise backward error is 1.1e-4 (transpose
3.7e-2); refinement reaches 1.4e-15 at the three-correction cap and stagnates near
4e-7 for the transpose. With 0.01 (factor size identical for 0.001 to 0.02):
2.4e-16 unrefined and 1.8e-16 (transpose 1.7e-16) after one correction, with 15%
fewer factor entries; forward agreement with UMFPACK stays at 1.6e-10.
Row-sum scaling selects the same pivots as row-max scaling; Ruiz row/column
equilibration made the unrefined error worse (4.9e-2). Unrefined thresholds
below 1 can lose normwise accuracy (with AMD and 0.01, `gre_1107` moves from
3.6e-16 to 1.7e-14 and `lns_3937` from 8e-17 to 5.7e-14), so the unrefined default
stays 1.0. On 18 Harwell-Boeing
general matrices with natural and AMD orderings, refined solves reached
componentwise errors of at most 4.5e-16 with either threshold, except `fs_760_3`
with AMD, where both stop at 1.4e-12 to 2e-12 (normal solve). Factor entries
with 0.01 changed by -36% to +31%. These are single-run accuracy observations
on one Mac build, not performance admission.

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
`solver%options%ordering=FORTSPARSE_ORDER_AMD` in every build.
AMD accepts the canonical unsymmetric CSC graph and orders its structural
`A+A^T`; numerical factorization borrows the original CSC values and applies the
column permutation through its retained sparse factors.
`FORTSPARSE_ORDER_ND` alternatively uses nested dissection on the
undirected graph union of `A` and its transpose. Self edges and repeated mirrors
are removed before ordering. With structural matching, either ordering consumes
the matched `P*A` graph while numerical values stay in the original CSC.
DEFAULT orders as described under [default ordering](#default-ordering).
NATIVE_REFERENCE remains natural and rejects explicit graph orderings.
`native_ordering_available(ordering)` queries the built native graph capability
without allocating or factoring. Matrix-class, device and memory-policy support
are still checked by factorization; this query is not a strict-cap guarantee.

Value-only refactors reuse the retained column order when the canonical pattern
and requested ordering strategy match; numerical row pivots are selected anew.
A changed pattern recomputes the order, and a changed ordering strategy invalidates
pattern-reuse reporting. Normal, transpose and adjoint solves, refinement and
implicit solve derivatives act on the original matrix.

Explicit uncapped AUTO-memory AMD uses the existing opaque `amd_order` path: actual
ordering and retained owned payload are reported, but complete host peak/factor
memory remains unknown (`-1`) and `strict_memory_supported` is false. Positive
host caps, LOW_MEMORY and the implicit default instead use the caller-workspace
path (audited AMD 3.3.4 `amd_2` when built, otherwise the project-owned core),
including old-factor, user-vector, refinement and permutation-input overlap
in cap admission and phase peak statistics. See [bounded AMD](native-amd-bounded.md) for exact workspace counts,
release trust, reuse and verification limits. Use measured process RSS for process
memory comparisons; these changes make no performance-admission claim.

The nested-dissection workspace is not budgeted: ND always reports unknown
complete peak (`-1`) and `strict_memory_supported=.false.`. Positive host caps or
LOW_MEMORY with ND fail before ordering; use natural or bounded AMD for those
contracts.
The same applies to the complex HPD and complex LDL factors.

## Opt-in supernodal real and complex LU

`solver%options%lu_method=FORTSPARSE_LU_SUPERNODAL` selects the static-structure
supernodal LU ([design](lu-supernodal.md)) for real or complex GENERAL matrices on the
explicit `FORTSPARSE_BACKEND_NATIVE` CPU backend. The default stays
`FORTSPARSE_LU_SCALAR`; `FORTSPARSE_LU_AUTO` behaves as SCALAR until its
thresholds are calibrated. The family reports `native_factor_family='LU_SUPERNODAL'` for real factors and
`'LU_SUPERNODAL_Z'` for complex factors.

- Analysis: maximum-product matching with dual scaling, then AMD (bounded
  caller-workspace AMD when available) or ND on the matched symmetrized
  pattern, an etree postorder, and the shared supernode partition and schedule.
  DEFAULT ordering selects ND when no host cap is set and LOW_MEMORY is not
  requested; otherwise AMD. An explicit AMD or ND request
  wins. `ordering_used` reports the ordering actually used. `structural_matching` does not apply; the matching is always
  weighted (`matching_used=FORTSPARSE_MATCHING_WEIGHTED`).
- Factor: left-looking panels, partial pivoting restricted to each
  supernode diagonal block, and pivots below `pivot_perturbation*||A~||_inf`
  (default 1e-13) replaced by `sign(p)` times that threshold. Statistics report
  `perturbed_pivots`, `block_pivot_swaps`, `pivot_growth`, `min_pivot` and
  `scaling_log10_range`.
- Threads: `nthreads` (0 = OpenMP default) runs one OpenMP task per
  elimination-subtree group, depending on its predecessor groups, as the SPD
  supernodal factor does. The team is one below 1e6 flops, at most a quarter
  of the groups, halved under a host cap until the per-thread row map, update
  scratch and statistics fit (`budget_reduced_threads`), and one under
  LOW_MEMORY. `active_threads` reports the actual team. Each destination applies
  its updates in a fixed order, so pivot counts and gate decisions do not depend
  on the team. BLAS calls inside tasks use the BLAS library's own threading,
  as for SPD: with an OpenMP BLAS from another runtime (Homebrew OpenBLAS uses
  libomp) set `OMP_NUM_THREADS=1`/`OPENBLAS_NUM_THREADS=1` to avoid nested
  oversubscription. The 64-wide panels at the top of the elimination tree form
  a serial chain: on KIN6D the tree critical path bounds subtree-task speedup
  at about 2.4 (ND), and 1.5-2.0 was measured
  ([stage 7 notes](lu-supernodal.md#9-openmp-and-generated-kernels)).
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
  family `LU_REAL` or `LU_COMPLEX` and the supernodal attempt's pivot counters. A refined user solve
  that misses the gate refactors once with the fallback when no host cap is set.
- Refactor: a same-pattern refactor reuses matching, scaling, order, partition
  and panels (`reused=.true.`). An entry outside the static structure or a
  failed gate repeats the full analysis. NORMAL, TRANSPOSE/ADJOINT and block solves,
  refinement, implicit JVP/VJP and context copies use the retained factors.
- Unsupported: structured symmetries,
  NATIVE_REFERENCE and devices. ND is
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
