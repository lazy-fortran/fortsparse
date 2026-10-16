# Private real consumer fixtures

`run_private_real.py` adapts a hashed private fixture manifest to the existing
native-family executables. The only new driver case is `real-mm`, with normal
`rhs4`: the original three physical RHS columns followed by a manufactured
known-vector RHS. Duplicate physical columns remain present. No application
physics or private numeric fixture is copied into this repository. The frozen
synthetic fixtures and `contract.json` remain unchanged.

The benchmark-only reader accepts square, real-general MatrixMarket coordinate
data already ordered by column and strictly increasing row, without duplicates.
It uses O(n+nnz) storage, validates finite values and indices, and rejects
truncation, surplus records, unsupported headers and malformed records. Array
files contain one scalar per record in column-major order. There is no implicit
symmetrization, reordering, duplicate summation or changed numeric value.

All three existing executables retain their timed original-full-CSC residual
and known/reference-center checks. The common I/O module adds a candidate
MatrixMarket export, with its cost included in each child equally. Real data
uses GENERAL symmetry; Intel PARDISO retains true CSR conversion. Direct
UMFPACK stays in its separate GPL executable; neither the shared I/O module
nor libfortsparse links UMFPACK. Intel remains separately linked from the
FortSparse/OpenBLAS executable.

The child's existing CSV forward diagnostic uses the historical aggregate
normalization. It is a preliminary check, not the private accuracy-policy
decision. `real_fixture_audit.py` makes that independent decision, using exact
rational arithmetic on stored binary64 A, b and candidate x. It applies the
unchanged policy's global and row-scaled residual thresholds per RHS, the
manufactured known-vector comparison, and physical comparisons against BOTH
retained independent centers. Manufactured RHS rounding remains explicit in
the private provenance; a known vector is not falsely declared an exact
solution to its rounded binary64 RHS. Exact residuals do not establish a
forward-error enclosure or continuous scientific accuracy.

The audit is outside the timed child. Its duration, child-plus-audit duration,
and complete adapter duration are recorded separately. Child-only time is not
a certificate-inclusive workflow. The exact audit is a benchmark oracle, not
a new production checker or mandatory application Python layer.

## Local recipe

Build the three existing targets and `test_native_family_real_io` with the
native CMake settings in the family README. For sequential MUMPS, set both
`METIS_INCLUDE_DIR` and `METIS_LIBRARY` to its matching local installation:
the imported MUMPS target can already own METIS before FortSparse's own cache
variables are consulted. Use the direct-UMFPACK README's explicit cache entries
for the local OpenBLAS provision; keep Intel oneMKL in its separate executable.

Affected checks use public tiny data only:

```sh
ctest --test-dir BUILD -j 2 --output-on-failure \
  -R '^native_family_(real_|fixtures$|selector$)'
```

The private NEO-2 fixture's manifest is local at
`/home/ert/data/sparse-dev-private/neo2-lorentz-tag6-2026-10-16/fixture.json`.
Its frozen policy SHA256 is
`c03da85fc3dfae66b8c09874a69faf78d10de0d9b368b18d6b91813046a42d54`.
Numeric outputs must stay beneath that private archive and outside Git.
Output directories must be new; existing evidence is never overwritten.

```sh
python bench/native_families/run_private_real.py \
  --fixture "$fixture_manifest" --policy-sha256 "$frozen_policy_sha" \
  --build "$benchmark_build" --output "$private_output" --mode prepare
```

`prepare` verifies hashes and constructs private input symlinks without running
a solver. During a controller-coordinated window, `--mode qualify --window-note
"WINDOW"` runs one attempt per existing external peer: SuperLU, UMFPACK IPC,
direct UMFPACK, sequential MUMPS and Intel PARDISO. `--backends umfpack` permits
an explicit qualification subset; it does not reduce native prerequisites.
Unavailable executables, failures, timeouts and failed accuracy gates remain
in the receipt. Each child is bounded below one second, with fixed physical
core affinity and explicit thread settings. This adapter does not launch a
five-trial timing or scaling campaign.

Native mode requires `--backend native` or `native-amd`, a coordinated window
note, and `--qualification RECEIPT [RECEIPT ...]` covering all five external
peers exactly once. It rechecks fixture/policy identities, comparator binary
and candidate hashes, and independently re-audits every external candidate
before any native execution. A stale binary, partial coverage or merely
asserted success is rejected. The independent prerequisite audits are setup
and have their own recorded duration. Root integration creates new binaries,
so it must produce fresh qualification receipts for those binaries.

An explicit `--mode measure --backend native-amd --qualification RECEIPT`
launches matched measurements after rechecking those same prerequisites. It
rotates all five external peers and the selected native implementation across
five fresh trials at thread budgets 1/2/4/8/16. `--trials` and `--thread-sweep`
permit explicitly scoped follow-ups. Every trial uses the unchanged four RHS,
has a unique private candidate export and receives its own independent exact
audit. Failures and timeouts stay in flushed `attempts.jsonl`; later trials
still run. The final receipt retains all attempts, raw hashes, source identity,
requested budgets and separate child/audit durations. This mode measures cold
factor-plus-RHS4 workflows only; it does not establish warm/reuse, transpose or
whole-application admission.

The five-peer comparison scope follows the existing family contract. It is
explicitly broader than the first four-peer private extraction qualification;
the frozen numeric accuracy policy itself has not changed. IPC parent and
helper RSS remain separate and their aggregate memory is unknown.

## Verification receipt

`private_real_qualification.json` retains only hashes and aggregate diagnostics,
not physical inputs, matrix values, RHS or solution vectors. On 2026-10-16,
all five external peers passed the fixed policy through the actual shared
drivers, after 16 distinct affected CTests passed with `-fcheck=all`. The
original four-peer private receipt was retained and IPC was qualified in a
separate fresh receipt. Each child took less than 0.11 seconds in this one-shot
check; the exact audit took approximately 0.53–0.55 seconds. These receipts
establish bounded adapter correctness, not performance ranking, calibrated
baseline, memory admission or frontier membership. No native execution on
the private fixture accompanied this change.
