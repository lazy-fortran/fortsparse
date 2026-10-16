# Private real NEO-2 aggregate metrics

This checkpoint retains benchmark metrics, dimensionless aggregate error
diagnostics, statuses and provenance hashes. Matrix, RHS, candidate exports and
per-RHS exact audit records remain in the private archive. The frozen fixture
and policy identities are recorded in `recipe.json`; it deliberately omits the
private fixture manifest and numeric exports.

The 150 trials cover six backends, five requested thread budgets and five fresh
processes, with all four RHS. Every wrapper/child/exact-audit status is retained
in `attempt-statuses.csv`. All passed. `summary.csv` has 30 successful-trial
configuration groups. Requested-budget wall speedup is each backend's own
one-thread median divided by its requested-budget median; requested-budget
efficiency divides this ratio by the budget. These exclude exact auditing and
are not active-team efficiencies: native factors remain serial and external
active teams are unknown. `audits.csv` keeps exact benchmark-audit costs and scalar
error diagnostics separately, including available data from failed audits when
present. Empty summaries with count zero mean unavailable.

`fronts.csv` compares each native budget with every external budget. Three of
five native configurations are robustly nondominated for child-wall/RSS and
completed-child-CPU/RSS; two remain indeterminate. This is one frozen physical
fixture, not broad production admission. UMFPACK IPC aggregate memory is
unknown; the report never sums parent/helper RSS or plots its parent as an
aggregate point.

Child wall includes input/reference I/O, original-A checks, numeric candidate
export and teardown. Timed factor/solve calls are separate. Completed child CPU
uses wait4 and includes reaped helpers once. Exact rational residual auditing
is separately timed benchmark verification, not a forward-error certificate.
Audit CPU and peak memory are unavailable. Child+audit wall is a partial sum
and excludes prerequisite audits, hashing and adapter orchestration.

Reproduce only aggregate reports/plots, without rerunning solvers:

```sh
python3 bench/native_families/report_private_real.py \
  /home/ert/data/sparse-dev-private/neo2-lorentz-tag6-2026-10-16/scaling-c7207e7 \
  --output /tmp/fortsparse-private-neo2-aggregate --plots
```

PNGs remain outside Git. `scaling-and-phases.png` separates child wall, completed
CPU, factor and solve calls; `audit-costs.png` shows the separate oracle cost;
`time-rss.png` displays measured single-process peaks and omits IPC's unknown
aggregate point. All markers show medians with full five-trial ranges.
