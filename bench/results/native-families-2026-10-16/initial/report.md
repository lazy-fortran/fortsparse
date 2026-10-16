# Native sparse-family initial baseline

Source revision: `42b8acaf5ac7eb50a648229a473e6a580ec162c8`. Raw attempts SHA256: `5eeda7dc902e8cce880a6723d58e33d8b3648ff714c0094fde5d85e30e5c6852`.

2500 attempts: 2350 ok, 150 timeout.
CPU interpretation: **parent main-entry CPU_TIME interval; helper reaped lifetime CPU separately; complete workflow CPU unavailable**. Prior raw fields are preserved.

The rows below show training, normal operation, RHS8, one requested thread. Wall is complete subprocess elapsed time; factor is a timed library call. Ranges span all successful trials. Timeout rows have no factor/RSS sample.

| Family | Backend | OK/all | Wall median [min,max] ms | Factor median ms | Parent peak MiB |
|---|---|---:|---:|---:|---:|
| neo2-mini5 | intel-pardiso | 5/5 | 9.874 [7.889,15.219] | 2.023 | 19.125 |
| neo2-mini5 | mumps | 5/5 | 5.200 [4.493,8.258] | 0.501 | 11.398 |
| neo2-mini5 | native | 5/5 | 5.733 [5.392,9.118] | 0.038 | 6.418 |
| neo2-mini5 | superlu | 5/5 | 5.490 [4.954,5.748] | 0.095 | 9.137 |
| neo2-mini5 | umfpack | 5/5 | 6.931 [6.212,8.964] | 2.128 | 6.301 |
| drift-real | intel-pardiso | 5/5 | 101.451 [93.978,107.533] | 46.585 | 47.711 |
| drift-real | mumps | 5/5 | 46.679 [44.804,51.034] | 22.079 | 36.121 |
| drift-real | native | 5/5 | 819.285 [771.962,841.266] | 726.998 | 92.227 |
| drift-real | superlu | 5/5 | 66.832 [58.319,72.936] | 36.367 | 37.586 |
| drift-real | umfpack | 5/5 | 67.482 [62.868,72.048] | 38.467 | 18.695 |
| drift-complex | intel-pardiso | 5/5 | 109.651 [108.333,118.551] | 52.034 | 56.664 |
| drift-complex | mumps | 5/5 | 58.255 [57.246,60.550] | 31.282 | 46.238 |
| drift-complex | native | 0/5 | timeout | unavailable | unavailable |
| drift-complex | superlu | 5/5 | 85.692 [84.541,88.939] | 50.244 | 57.344 |
| drift-complex | umfpack | 5/5 | 85.705 [84.660,100.966] | 52.277 | 18.438 |
| hpd-chain | intel-pardiso | 5/5 | 50.951 [49.525,55.225] | 14.060 | 39.480 |
| hpd-chain | mumps | 5/5 | 32.723 [31.857,33.365] | 12.992 | 27.586 |
| hpd-chain | native | 5/5 | 21.739 [21.100,26.640] | 1.953 | 15.770 |
| hpd-chain | superlu | 5/5 | 23.577 [23.278,26.820] | 7.046 | 24.645 |
| hpd-chain | umfpack | 5/5 | 30.608 [30.073,31.299] | 12.218 | 16.336 |
| ldlt-chain | intel-pardiso | 5/5 | 42.782 [40.481,45.857] | 14.181 | 34.723 |
| ldlt-chain | mumps | 5/5 | 43.917 [43.401,46.631] | 23.568 | 26.223 |
| ldlt-chain | native | 5/5 | 23.860 [22.373,32.267] | 5.164 | 20.004 |
| ldlt-chain | superlu | 5/5 | 23.754 [21.816,31.151] | 6.433 | 22.746 |
| ldlt-chain | umfpack | 5/5 | 31.127 [29.202,31.561] | 13.419 | 16.852 |

Timeout counts (all thread budgets, variants and trials):

- drift-complex, refactor: 50.
- drift-complex, rhs8: 50.
- drift-real, refactor: 50.

Interpretation and limits:

- Natural-order native LU is a clear bottleneck on the retained 2D drift cases. Complex drift times out for both workloads; real drift times out for refactor.
- The HPD chain has a small native factor-call cost; total process cost is much larger. The chain and mini5 fixtures do not establish broad production leadership.
- Native records report one active factor thread at every requested budget. External active teams are unknown (-1), not evidence that they used the requested team.
- Parent VmHWM and helper RUSAGE_CHILDREN ru_maxrss remain separate. Helper RSS can include inherited/pre-exec highwater. Never sum them or use their maximum as a validated aggregate lower bound. UMFPACK aggregate memory is unknown; its parent VmHWM alone is a safe lower bound.
- CPU in the initial receipts is the parent main-entry interval plus separately recorded helper lifetime CPU. These do not measure complete workflow CPU. New wait4 instrumentation is a separate measurement checkpoint.
- Successful-only timing ranges do not replace timeout evidence. A timeout marks a censored observation at the configured deadline, not a measured factor time.
- No global Pareto/admission claim follows: this is a bounded fixture corpus, several timings are censored, true helper aggregate memory is unknown, and the shared host has background load.
- Observed pre-launch one-minute host load: 4.64–11.23. Core affinity reduces placement variation; it does not eliminate shared-host interference.
