# Native structural matching: bounded coupled-Jacobian evidence

Source `8914ede`; this is a correctness/graph-fill diagnostic, with one trial per
configuration. It makes no performance/Pareto or full scientific-case claim.
The original KIN6D radial-level 8/16/32 convergence test remains a separate gate.
No large-mesh factorization was run for this checkpoint.

The public manufactured first Jacobian from KIN6D baseline
`73891c54a71c3a080826043fd02924779a6e79a3` has `n=5525`, `nnz=203054`.
Its SHA256 is `707d5e280ed1f345f12581a1086be4bfe42ff00e5a1291e87b3176b4b5288078`.
The original matrix and all physics/constraints remain unchanged. There are
737 zero/missing original diagonals, 640 numeric row singletons, maximum row
and column degrees 4148/3411, and 95.77% transpose-present nonzero pattern.

| Preparation, pivot threshold | Factor entries | Retained / peak owned bytes |
|---|---:|---:|
| AMD, 1.0 | 10,179,016 | 138,769,328 / 189,100,976 |
| COLAMD, 1.0 | 7,210,663 | 152,138,672 / 202,470,320 |
| External SciPy structural matching + AMD, 0.1 | 4,895,017 | 76,641,200 / 101,807,024 |
| External SciPy structural matching + COLAMD, 0.1 | 5,324,450 | 76,641,200 / 101,807,024 |
| Native structural matching + AMD, 0.1 | 4,321,295 | 76,641,200 / 101,851,224 |

The last row uses the public native API and a 512 MiB owned-array cap. Its
matching/ordering analysis took 0.00718 s, numerical factorization 0.57895 s,
and the enclosing sparse-factor call/statistics sample 0.58617 s. Independent
quadruple-precision accumulation of original COO entries gave normwise backward
error `2.387e-15`, componentwise backward error `2.152e-14`, and relative known-vector
difference `1.261e-9`. The manufactured RHS is rounded to binary64; no rigorous
forward-error certificate follows from these numbers. Refinement is disabled.

The native matcher returns a perfect nonzero original-row bijection, independently
checked against the original entries. It differs from SciPy's valid matching at
2210 columns: native traversal starts from original CSC columns, while SciPy
1.18.1 starts from original CSR rows. Their DFS schedules also differ. Equality
of maximum matchings is neither required nor claimed. The exact permutations,
SciPy source, installed extension hash and both traversal descriptions are retained.

The first four rows measure direct LU factor-owned payload. The last measures
public-backend payload, including live preferred-row and ordering input arrays
at the numerical peak, and all charged matching/AMD analysis arrays at their
actual phase peaks. These scopes are distinct; none is process RSS. COLAMD's
lower factor-entry count in the unmatched diagnostic coincided with HIGHER
allocated capacity; it is not a payload-memory win. Complete process CPU, RSS,
controlled repetitions and total-workflow comparisons are unavailable here.

Six initial unmatched probes at 1.0/0.1 and both matched probes at 1.0 timed out
at a one-second numerical deadline. Their missing counters remain blank in
[the compact evidence](2026-10-16-structural-matching.csv). Two explicitly
extended unmatched AMD/COLAMD probes completed within a 15-second limit; these
are diagnostic counters excluded from ordinary sub-second benchmark admission.
Preprocessing/input parsing and known-RHS construction precede the READY/start
handshake. The native producer also has a separate matching-only preflight to
retain its permutation, so no cold-cache claim follows. Host load was shared;
these unreplicated times do not establish a ranking or robust range.

## Verification and reproduction

Eight affected bounds-checked Debug CTests passed, including the new graph/API/
allocation oracles, generated row hints, existing general LU, bounded AMD and
original NEO2 spline regression. Focused fpm graph/API tests also passed without
AMD. Default matching-disabled policies remain unchanged. Source details and
workspace/refactor semantics are in [native matching](../../docs/native-matching.md).

Retained public input/capture patch and all raw failures, Fortran diagnostic
sources, C ordering shim, runners, permutations, exact commands and compiler/
affinity/thread identities are under:
`/home/ert/code/sparse-dev/benchmarks/2026-10-16-kin6d-native-migration/nonnested-csc-probe/`.
`native-matching/provenance.json` pins source `8914ede`, the Release static library,
probe ELF and original input. Probe ELF SHA256:
`c43737db2762edcf0335fbd65b922f33c496c1ba4f0db133560a963afeb63440`.
The native numerical deadline is one second, OMP/OpenBLAS budgets are one, and
inherited affinity is recorded. Configure native-only Release with AMD on,
METIS/external factorizations off; build with at most two cores. Compile retained
`probe_native.f90` against that library and its module directory, then execute
`native-matching/run_native.py` in that directory only in a controller-approved
window. No Python is linked to or required by production matching.
