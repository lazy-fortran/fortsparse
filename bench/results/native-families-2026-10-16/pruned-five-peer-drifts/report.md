# Five-peer completed CPU/RSS checkpoint

1200 completed attempts; 240 five-trial configuration groups. Source `7e1f123ea21fd8412480ec6a6ddf0696893acb2a`. Native backend labels: native-amd. All process CPU uses completed wait4 user+system counters including reaped helpers.

Five external peers: SuperLU, UMFPACK IPC, standalone direct UMFPACK, MUMPS and Intel oneMKL PARDISO. Requested physical thread budgets are 1/2/4/8/16. Native active factor teams reported: [1]. External active teams are retained separately in the numeric CSV; negative values mean unknown.

External-only calibration retained 2,052 attempts: 2052 ok. Rejected observations under the strict <0.5s process-wall/status gate: 1. Rejected observations and unmeasured candidate tails are retained; their success is never inferred. Earlier differently sized diagnostics remain separate.

The following counts compare each of five native requested-thread configurations with all 25 external peer/thread configurations at the same family, actual size, workload, variant and operation:

| Family | Workload | Variant | Wall: robust / indeterminate / dominated | CPU: robust / indeterminate / dominated |
|---|---|---|---:|---:|
| drift-complex | refactor | heldout | 0 / 5 / 0 | 0 / 5 / 0 |
| drift-complex | refactor | training | 0 / 5 / 0 | 0 / 5 / 0 |
| drift-complex | rhs8 | heldout | 0 / 5 / 0 | 0 / 5 / 0 |
| drift-complex | rhs8 | training | 0 / 5 / 0 | 0 / 5 / 0 |
| drift-real | refactor | heldout | 0 / 5 / 0 | 0 / 5 / 0 |
| drift-real | refactor | training | 0 / 5 / 0 | 0 / 5 / 0 |
| drift-real | rhs8 | heldout | 0 / 5 / 0 | 0 / 5 / 0 |
| drift-real | rhs8 | training | 0 / 5 / 0 | 0 / 5 / 0 |

Actual measured sizes: drift-complex n=4096, drift-real n=16384.

Assessment rule and limits:

- Robust observed-range nondominance excludes every external configuration by native maximum time < external minimum time OR native maximum process RSS < external minimum process RSS. Strict inequality treats observed overlap as indeterminate. These trial ranges are descriptive, not confidence intervals.
- IPC aggregate RSS is unknown. Only its parent VmHWM is used as a safe lower bound for excluding a potential dominator. Helper ru_maxrss can retain inherited highwater; no sum or maximum of parent/helper is an aggregate peak. IPC cannot be declared a definite dominator with these memory diagnostics.
- A definite dominator has strictly lower maximum time and process RSS than the native minimum values. Such a result requires known single-process RSS.
- Nondominance is a time/memory tradeoff, distinct from matching or beating the same best comparator on both axes. This report does not promote a library default or claim broad production/frontier completion. It covers 2 fixture families, operations N, workloads refactor,rhs8 only.
- CPU and elapsed front assessments are separate. The refactor workflow includes prime factorization. The numeric CSV separates timed analysis, factor, numeric, solve and refactor calls from complete subprocess wall and CPU.
- Native labels with unavailable attributed peak bytes: ['native-amd']. Process VmHWM is valid RSS evidence and must not be called attributed factor memory. Unknown summaries have count zero and must not be plotted as zero.
- Five-peer selection differs from the older four-peer diagnostic. No cross-checkpoint speedup at different sizes or host loads is inferred.

Calibration rejections:

- drift-complex, level 128, umfpack, 4 threads, trial 0: status ok, process wall 0.747116s.
- Pre-launch one-minute host load: 0.417–2.927. Affinity and controller quiescence retain ordinary resident services; this is not an exclusive host.

Reproduction and provenance:

- `/home/ert/code/sparse-dev/benchmarks/2026-10-16-native-families/pruned-five-peer-drifts`: raw SHA256 `bbb373762d0b1347c5ee1c97c3c8de703dff4f9e1add4b11427c774d9f4acda1`.
- Calibration selection SHA256 `cb39c751f3832a279e3e378ef33b47fc3cf9c714ee5128e3aedf380e790a42dc`; all rejected observations remain in its sibling `attempts.jsonl`.
