# Native LU work diagnostic, 2026-10-16

The unchanged KIN6D `nonnested_coupled_solovev_p2` test (`solovev 2 1e6`,
8/16/32 radial levels) timed out at **900.64 s**. Native structural matching,
matched-graph AMD and pivot threshold 0.1 pass the first two levels. The largest
level remains unqualified. This attempt compiled FortSparse aca8c86; ee1cdc9
only wraps a retained-byte sum for older Fortran compilers. Current consumer
checkpoint d7dbf24 pins ee1cdc9; eleven other affected checks pass. No completed
full-case CPU/RSS delta is available. Earlier six-case measurements remain
explicitly at the historical dfa0b79 checkpoint in its draft PR.

A one-second sample during that qualification run captured 146 cycles samples
with none lost: 93.38% in `spd_scatter_update`, 5.51% in LU DFS `visit`, 1.11%
in `factor_raw`. This sampled window alone does not identify whole-run fill,
pivot instability or a root cause.

One separate diagnostic replay uses the original first 19,989-row Jacobian
(801,006 entries), SHA256
`13ffe01f42f89dcefd7a692a9d0a725e0983b84fbe7f2d4445f3ab5cb447b1e5`.
It links ee1cdc9 and wraps the generated real scatter entry point, forwarding
every call unchanged. It requests matching + AMD + threshold 0.1, one OpenMP/
BLAS thread, a **2 GiB owned-payload cap**, and a 15-second diagnostic deadline.
The safe pre-termination counter snapshot records:

| Prefix counter | Count |
|---|---:|
| Scatter calls | 31,484,894 |
| Scatter scalar entries | 21,691,183,792 |
| Calls with zero multiplier | 29 |
| Entries in zero-multiplier calls | 21,732 |

The factor did not return before termination. Final factor nnz, allocated
capacity, factor memory peaks and residuals are unavailable. Zero multipliers
account for about one part per million of observed entries and cannot explain
the dominant prefix work. Do not infer a completed-factor cost or optimize zero
skips from this result. Next distinguish fill/pivot departures from repeated
updates and per-entry cost using factor-ordinal and capacity counters.

These are instrumented, deadline-censored diagnostics on a shared host;
the replay overlapped the full qualification run. They are excluded from
sub-second timing/admission evidence. Source, binary/library hashes, wrapper,
safe signal/pipe protocol, raw counters and commands are retained under
`sparse-dev/benchmarks/2026-10-16-kin6d-native-migration/`, in
`qualification-matching/` and `nonnested-csc-probe/19989-diagnostic/`.
