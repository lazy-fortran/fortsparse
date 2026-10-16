# LOW growth A/B: completion, traffic and RSS

The ordinary frozen comparison has 240 interleaved CPU0 attempts: baseline
7e1f123 completes 100 and times out on all 20 real-LOW attempts; candidate
2f0cb01 completes all 120. Its 48 version/configuration groups each contain
five trials across real 16k/complex4k, training/heldout, cold RHS8 and
prime+refactor. `summary.csv` is regenerated directly from immutable raw
receipts, independently of the previously corrected archived summary filter.
No solver, build or new timing invocation is part of this report recipe.

The three policies are AUTO without a cap, AUTO with a 67,108,864-byte cap
(64 MiB), and LOW without a cap. **The separate [positive-cap LOW checkpoint](../low-positive-cap/README.md)
contains 40/40 passes at source cac6656**, with real 16 MiB/complex 5 MiB caps.
It is not merged with this 240-trial uncapped-policy A/B. Dense80 instrumented
retry tests remain distinct from that ordinary positive-cap campaign.
Later production integrations f4f329d/9354c60 are distinct from these measured
snapshot identities. The manifest retains executable/library hashes and the
candidate's additional METIS dependency.

For training complex LOW, median whole-process wall falls from about 272 to 25 ms
for RHS8 and 259 to 30 ms for refactor. Candidate real LOW completes around 74 ms
and 90 ms respectively; the killed baselines provide censoring evidence, not
finite completed-time speedup estimates. Full ranges, completed wait4 CPU,
factor/solve/refactor phases, fill and numerical diagnostics remain in the CSV.
There is no external solver or parallel scaling comparison in this native-only
policy A/B, and no frontier admission follows.

Candidate real training RHS8 LOW retains 10,611,660 bytes versus AUTO 13,891,596,
but median process RSS is higher: 33,440 KiB (32.66 MiB) versus 29,948 KiB
(29.25 MiB). Retained live payload is not process peak RSS. Neither a lower RSS
claim nor an allocator-cause claim is made; allocator diagnosis is separate.

`diagnostic-counters.csv` summarizes the independently instrumented real LOW
replays. Growth events fall from 32,767 to 317 and cumulative factor-prefix
copying from 54,871,019,988 to 148,759,548 bytes. The 54.87 decimal GB is cumulative
traffic, **not peak RAM**. Diagnostic phase-maximum owned payload increases
slightly, 14,899,228 to15,187,408 bytes. Instrumented wall times are excluded
from ordinary runtime comparisons. The older source 3b9384b sampling profile
has 740 samples, no losses, 64.99% own cycles in reserve and 5.23% in visit; it is
identified separately from the ordinary 7e1f123 baseline.

The dense80 diagnostic records one exact retry: success at 144,812 bytes and
failure at 144,811, with required 38,880 versus available 38,879 bytes. Peak is
the maximum across sequential attempts, never their sum. Detailed source,
counter, retry and profile documents stay in the retained diagnostic archive;
`recipe.json` pins their paths/hashes.

Regenerate aggregate reports/plots without executing a solver:

```sh
python3 bench/results/native-families-2026-10-16/low-growth-ab/reproduce.py \
  --output /tmp/fortsparse-low-growth-report --plots
```

The recipe verifies raw, ordinary manifest and diagnostic-manifest identities,
then checks primary diagnostic files against the manifest. It keeps all 20
watchdog failures visible and leaves their completed CPU/memory unavailable.
The three inspected PNGs remain outside Git; `publication.json` retains their
hashes and ephemeral slopbox URLs, which expire after three days.
