# Original KIN6D coupled-case comparison attempt

One unchanged full-case pair on the original Solovev P2 levels 8/16/32:
`test_nonnested_coupled solovev 2 1e6`. Baseline KIN6D `73891c54` passes;
native `fc8b36a`, using FortSparse `886f7dc`, reaches the 900-second watchdog.
The test drivers are byte-identical. No scientific tolerance or mesh changed.

| Outcome | Wall s | CPU s, user + system | Process peak RSS KiB |
|---|---:|---:|---:|
| Original baseline: PASS | 6.049655176 | 11.012585 | 642452 |
| Native: TIMEOUT, failed attempt | 900.025097186 | 895.873448 | 7398032 |

No successful native CPU/RSS delta is available. Its values describe consumed
resources before termination, not a completed solve or a feasible Pareto point.
Buffered native stdout was empty, so its last completed scientific stage is
unknown. No retry or deadline increase followed. The child was reaped and no
observer, watchdog or target process remained.

Both runs requested four OpenMP threads and one BLAS thread on physical CPUs
0–3, with passive OpenMP waiting. Actual solver teams were not measured. This
is one pair on a shared AMD Ryzen 9 5950X host; unrelated services remained.
Startup is included, compilation and profiling are excluded. CPU comes from
one `wait4` user/system record, without adding helper CPU again. RSS is the
individual process peak, not a concurrent workflow sum.

The small C spawn/wait observer was independently qualified against post-exec
memory and affinity readback under a larger Python parent. Source, executable,
cache and all resolved-library hashes were unchanged after the run. Numeric
data, flags and source/driver/ELF identities are retained in `plot-data.json`;
the external campaign retains raw records and the exact observer recipe.

A separate native registered correctness run passed in 797.81 elapsed seconds
under different load/affinity. That observation is preserved and is not
substituted for the failed measured attempt. The subsequent KIN6D `9eb9552`
change only classifies this long test as `slow`. This application case is
outside the sub-second solver benchmark corpus and its SPD frontier claim.

Generate the inspected figure with optional Matplotlib:

```sh
python3 plot_timeout.py --input plot-data.json --output /tmp/kin6d-nn-timeout.png
```

Keep generated PNGs outside source history and publish them through slopbox.
No repeated timing, scaling, GPU or full-case performance admission is claimed.
