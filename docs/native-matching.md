# Opt-in structural matching for native general LU

Set `solver%options%structural_matching=.true.` for native CPU GENERAL real or
complex LU. The default is false. Other matrix classes, accelerators and
selected external backends return `FORTSPARSE_BACKEND_UNAVAILABLE`; AUTO does
not silently ignore this option if it selects an external backend. No new
matching dependency or C ABI option is introduced.

The native Fortran graph analysis computes a deterministic maximum-cardinality
bipartite matching between CSC columns and original rows, excluding numerical
zeros. A perfect match supplies an original-column to preferred-original-row
bijection. A Hall-deficient graph returns `FORTSPARSE_SINGULAR`; nonfinite or
noncanonical input returns `FORTSPARSE_INVALID_MATRIX`. A perfect graph match
alone does not prove numerical nonsingularity or good conditioning.

The algorithm uses breadth-first shortest augmenting-path layers and an
iterative depth-first search (Hopcroft--Karp), with no recursion or dense matrix.
It visits roots in ascending original-column order and original rows in stored
CSC order. This differs from the external SciPy diagnostic's row-oriented CSR
traversal; maximum matchings are nonunique. Neither numerical weighting nor
numerical equilibration is added by matching. Existing row scaling and pivot
thresholds are preserved. No pivot is accepted merely because it was matched:
the same numerical threshold still applies, and singleton pivots retain priority.

With AMD, analysis sorts matched row labels for the graph of `P*A`, then orders
its structural symmetrization. Numerical values and caller CSC arrays remain
unchanged. LU receives the column order and borrowed original-row preferences;
normal, transpose, adjoint and implicit derivative solves still act on the
original equation. Without AMD, the natural column order uses the same matched
row preference. Matching is a discrete solver choice, not a differentiable
operator; solve JVP/VJP continue to differentiate `A*x=b`.

Every factor/refactor recomputes matching and its ordering, including same stored
patterns: changing a stored value to zero can invalidate a previous match.
Pattern/order reuse is conservatively reported false for these requests and
when switching from matching back to an ordinary ordering. Factors remain
retained for subsequent RHS solves. Failure makes the public factor unusable
until a successful factorization; transient matching arrays are released.

For default four-byte indices, analysis owns two `n`-entry maps and five
`n`-entry work vectors: `28*n` bytes. AMD also needs `4*nnz` bytes for sorted
matched graph labels. These arrays coexist with old factors, retained originals,
user vectors and AMD's separately charged workspace during analysis. Only the
`4*n` preferred-row map remains during numerical factorization; it is charged
against the factor budget, then released. Actual index storage size determines
the implementation's byte count. There is no persistent numerical CSC copy.
LOW_MEMORY uses the same graph storage and existing compact LU growth policy.
Positive caps and LOW with AMD require its existing audited caller-workspace ABI;
uncapped opaque AMD retains its existing unknown complete-peak status. Existing
refinement admission conservatively reserves the future checker quota even
before those arrays are allocated. Owned array payload is distinct from RSS.

Independent checks exhaust all 512 three-by-three graphs against Hall's subset
criterion, exercise long augmentations and Hall failures with no empty rows or
columns, and verify real/complex original-matrix N/T/H, block RHS, refactors,
implicit JVP/VJP and failure recovery. Linux allocation interception checks exact
retained/peak payload, cap refusal before allocation, cleanup and recovery after
each reached allocation failure, including bounded AMD and complex checker
storage. These checks are numerical/structural verification, not certificates.

The concrete motivating input is the public manufactured KIN6D first coupled
Jacobian (`n=5525`, `nnz=203054`); the retained external diagnostic recipe and
provenance are under `sparse-dev/benchmarks/2026-10-16-kin6d-native-migration/`.
The original larger scientific convergence test remains a separate qualification
gate. A bounded successful factor does not establish performance admission or
predict the largest mesh.

Algorithm reference: Hopcroft and Karp, *An n^(5/2) Algorithm for Maximum
Matchings in Bipartite Graphs*, SIAM Journal on Computing 2(4), 1973,
[doi:10.1137/0202019](https://doi.org/10.1137/0202019). The external diagnostic
uses [SciPy's maximum-bipartite-matching contract](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csgraph.maximum_bipartite_matching.html)
and separately retained SciPy 1.18.1 traversal source. Production matching is
project-owned Fortran; SciPy is an external diagnostic oracle only.
