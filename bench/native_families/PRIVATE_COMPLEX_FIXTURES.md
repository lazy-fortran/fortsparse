# Private complex consumer adapter

`complex-mm` adds one concrete benchmark-only interchange case to the three
existing executables. It uses general complex LU, the original full matrix,
and `rhs2`: one physical RHS followed by one independently manufactured RHS.
It does not use NEO-2's real three-physical/one-known convention. No solver
algorithm, product API, frozen synthetic fixture or accuracy threshold changes.

The common Fortran reader accepts only canonical complex-general MatrixMarket:
coordinate matrix entries sorted by column and unique row, with separate real
and imaginary values; array data column-major with real/imaginary pairs.
`matrix.mtx`, `rhs.mtx` and `centers.mtx` are authoritative. Arrays have two
columns in the declared physical/known order. Duplicate/surplus/truncated,
nonfinite, malformed and unsupported input is rejected. No symmetry expansion,
coalescing, zero removal, scaling or other numeric transformation is hidden.
Any extraction canonicalization belongs in the hashed fixture provenance.
Candidate export preserves both binary64 components and refuses overwrite.

The existing original-full-CSC cheap residual and center checks remain in each
child, along with input/export costs. Public tiny tests independently generate
N/T/H actions for a non-Hermitian complex matrix. The actual saved MEPHIT
consumer collection declares N only; no transpose RHS is inferred from a
physical normal RHS. Direct UMFPACK remains in its separate GPL executable,
and Intel PARDISO remains separately linked from FortSparse/OpenBLAS.

## Independent audit and eligibility

`complex_fixture_audit.py` computes Gaussian-rational residual components from
stored binary64 matrix, RHS and candidate values. The complex absolute value
is Euclidean modulus, never an undeclared real/imaginary L1 norm. For a positive
rational squared modulus q, scale its root to at least 128 bits, take the integer
floor square root, and bracket with adjacent rational endpoints. Exact squares
and zero remain exact. Monotonic addition, maximum and positive division give
conservative lower/upper norm and error-ratio bounds.

Acceptance requires the ratio **upper** bound to satisfy each frozen threshold.
A lower bound above the threshold fails; an enclosure overlapping the threshold
is `audit_indeterminate`, never success. Normwise backward error uses
`||r||inf/(||A||inf*||x||inf+||b||inf)`. Row-scaled error uses
`max_i |r_i|/(sum_j |A_ij|*||x||inf+|b_i|)`, including the global vector norm.
Physical and known-vector errors both divide by the respective center infinity
norm. Zero-over-zero is zero; positive-over-zero fails. These are discrete
residual/empirical center checks, not a forward solution enclosure or a
continuous scientific certificate. Rounded manufactured RHS remain rounded
fixed inputs, rather than becoming an asserted exact solution identity.

The independently frozen policy sets backward gates 1e-12, physical/known
center gates 1e-10, 128 relative enclosure bits, an audit budget of 0.5 s, and
combined child/audit eligibility below 0.9 s. Audit deadlines include candidate
reading and are checked during exact arithmetic. The runner passes only the
remaining combined budget to the audit. Timeout, ambiguous enclosure or budget
exhaustion stays explicit and cannot qualify a case. Extraction, hashing and
one-time source parsing are separate preparation costs. The full 63988 system
may be ineligible; complete 606/814 surface systems can be selected by external
feasibility only, before native execution. No native result may select a case
or relax the frozen policy.

## Explicit private workflow

`run_private_complex.py` reuses the existing launch/privacy/qualification
boundary with a concrete complex fixture object. Staging filenames, driver
case, workload and audit module come from that object; the existing real
workflow keeps its original semantics. An explicit collection member ID,
collection schema, frozen policy hash, file sizes/hashes, roles and dimensions
are checked. Qualification receipts bind the member ID as well as collection
and policy hashes, preventing cross-surface reuse.

```sh
python bench/native_families/run_private_complex.py \
  --fixture "$private_root/manifest.json" \
  --fixture-id mephit-cylinder-surface64-v1 \
  --policy-sha256 221ee99ec8b4cc3ebb260beadb228bfcf9c8027e51ed9b761e179ee46fc5f3e2 \
  --build "$benchmark_build" --output "$new_private_output" --mode prepare
```

The local archive is
`/home/ert/data/sparse-dev-private/mephit-cylinder-2026-10-16`. IDs are
`mephit-cylinder-full-v1`, `mephit-cylinder-surface64-v1` and
`mephit-cylinder-surface154-v1`. Numeric inputs, centers, exports and audit
receipts remain under that archive and outside Git. Preparation runs no solver.

External qualification needs an explicit controller-coordinated window. Native
mode rechecks successful qualification for all six external configurations,
including both direct UMFPACK control variants, against this exact member,
policy, binary and candidate hashes; it independently re-audits every candidate.
Missing, stale, failed, timed-out or ineligible external evidence blocks native
execution. The supplied adapter has no qualification or performance claim.
The historical real/private fixed-scope report is not a report for this complex
collection; any publication needs its own declared scope and evidence.

The historical physical center is a saved UMFPACK result. The physical RHS was
reconstructed from retained physical inputs with a declared grouping; original
pre-solve RHS memory bits and producer executable identity were not recorded.
Those provenance limits remain explicit and are not cured by solver agreement.
