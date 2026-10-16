# Rejected compact LU visited marks

The byte-mark candidate `3033a3a` was not integrated. On 2026-10-16, a quiet,
interleaved comparison against frozen indexed baseline `c6fa39b` ran 80 fresh
processes: real drift n=16384 and complex drift n=4096, training/held-out
coefficients, RHS8/refactor workflows, five trials of each binary, one pinned
physical core. Both used the same system runtime libraries through an explicit
`LD_LIBRARY_PATH=/usr/lib`. No external solver was remeasured in this experiment.

All 80 accuracy/status checks passed. All 40 matched pairs had identical factor
fill and exactly `3*n` fewer attributed owned bytes in the candidate. Nevertheless,
real training RHS8 median total time was 60.776 versus 60.661 ms, and factor time
was 46.482 versus 46.489 ms (baseline versus candidate). Held-out real factor time
was 46.384 versus 46.656 ms. Real refactor phase medians were effectively unchanged.
Complex n=4096 results showed only small, inconsistent changes. RSS did not show a
reproducible decrease: real training RHS8 medians were 29420 versus 29572 KiB.
The reduced marker payload did not establish an end-to-end time or RSS benefit.

`summary.csv` retains medians and observed ranges, including complete-process
CPU and RSS. `recipe.json` records sources, binary/library hashes, exact case
product, interleaving, environment, affinity, timeout, accuracy and raw receipts.
Ranges from five trials are descriptive, not confidence intervals. Cross-campaign
absolute timings differ from earlier checkpoints; these A/B measurements do not
establish frontier admission. The rejected code remains on its experiment branch
for recoverable cleanup; the production implementation keeps generation stamps.
