# Pruned LU drift checkpoint, source 7e1f123

All 1,200 trials passed. All 80 native configuration/metric assessments remain
**indeterminate** against the five external peers and every requested thread
budget. This cohort establishes no robust observed-range frontier membership;
indeterminate is distinct from a proved external dominator or numerical failure.

This checkpoint contains only the pruned drift campaign frozen at
`7e1f123ea21fd8412480ec6a6ddf0696893acb2a`. Its frozen ELF identities and dynamic
library resolution are in `provenance.json`, with raw path/hash receipts in
`recipe.json` and `reproduction.json`. Historical c6fa39b chain measurements
are not joined. The historical external-only calibration is retained solely as
the pre-tuning size gate: real drift n=16,384 and complex drift n=4,096.

Coverage is **historical-five**: SuperLU, UMFPACK IPC, standalone default-control
UMFPACK, MUMPS and Intel oneMKL PARDISO. It does not include the subsequently
added direct UMFPACK no-refinement comparator and claims no contract-v4-six
coverage. There are five fresh trials for each of 240 configurations spanning
training/heldout coefficients, normal solves, cold RHS8 and prime+refactor,
with requested physical thread budgets 1/2/4/8/16.

`summary.csv` retains successful-trial counts and min/median/max values for
complete subprocess wall, completed wait4 CPU, separate parent/helper CPU/RSS,
solver phases, fill, attributed memory diagnostics and numerical errors.
`fronts.csv` retains all 80 assessments, each against all 25 peer/thread
configurations. Maximum observed backward error is 1.051e-15; maximum
known-solution relative error is 2.679e-15. These are numerical diagnostics,
not a forward-error certificate.

`requested-scaling.csv` reports each backend's own one-thread wall median
ratio t1/tP and requested-budget efficiency (t1/tP)/P. Native factorization
reports one active thread at every budget; external actual teams are unknown.
The ratios therefore describe requested budgets, not measured active-team
parallel efficiencies. Setup, solves, checking and teardown remain in process
wall/CPU; factor, solve and refactor calls are also retained separately.

UMFPACK IPC aggregate memory remains unknown. Parent VmHWM alone can safely
exclude a possible dominator when native RSS is strictly smaller than that
lower bound. It cannot establish IPC aggregate peak RSS. Helper lifecycle
ru_maxrss may retain inherited highwater; neither sum nor maximum is an
aggregate measurement. Time/RSS figures omit an IPC aggregate point explicitly.

Reproduce the report and PNGs without executing a solver:

```sh
python3 bench/results/native-families-2026-10-16/pruned-five-peer-drifts/reproduce.py \
  --output /tmp/fortsparse-pruned-drift-report --plots
```

The recipe verifies pinned raw, source/binary provenance and calibration
identities, then uses the existing historical-five report path. PNGs and raw
attempts remain outside Git. `publication.json` records hashes and ephemeral
slopbox URLs for the eight inspected plots; those URLs expire after three days.
