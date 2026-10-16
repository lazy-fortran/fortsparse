# Five-peer completed CPU/RSS checkpoint

2400 completed attempts; 480 five-trial configuration groups. Source `c6fa39b692c4d547a6032e03aa5e9343ab4734e3`. Native backend labels: native, native-amd. All process CPU uses completed wait4 user+system counters including reaped helpers.

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
| drift-real | rhs8 | heldout | 0 / 4 / 1 | 0 / 4 / 1 |
| drift-real | rhs8 | training | 0 / 4 / 1 | 0 / 4 / 1 |
| hpd-chain | refactor | heldout | 5 / 0 / 0 | 5 / 0 / 0 |
| hpd-chain | refactor | training | 5 / 0 / 0 | 5 / 0 / 0 |
| hpd-chain | rhs8 | heldout | 5 / 0 / 0 | 5 / 0 / 0 |
| hpd-chain | rhs8 | training | 5 / 0 / 0 | 5 / 0 / 0 |
| ldlt-chain | refactor | heldout | 5 / 0 / 0 | 5 / 0 / 0 |
| ldlt-chain | refactor | training | 5 / 0 / 0 | 5 / 0 / 0 |
| ldlt-chain | rhs8 | heldout | 4 / 1 / 0 | 4 / 1 / 0 |
| ldlt-chain | rhs8 | training | 4 / 1 / 0 | 4 / 1 / 0 |

Actual measured sizes: drift-complex n=4096, drift-real n=16384, hpd-chain n=16384, ldlt-chain n=16384.

Assessment rule and limits:

- Robust observed-range nondominance excludes every external configuration by native maximum time < external minimum time OR native maximum process RSS < external minimum process RSS. Strict inequality treats observed overlap as indeterminate. These trial ranges are descriptive, not confidence intervals.
- IPC aggregate RSS is unknown. Only its parent VmHWM is used as a safe lower bound for excluding a potential dominator. Helper ru_maxrss can retain inherited highwater; no sum or maximum of parent/helper is an aggregate peak. IPC cannot be declared a definite dominator with these memory diagnostics.
- A definite dominator has strictly lower maximum time and process RSS than the native minimum values. Such a result requires known single-process RSS.
- Nondominance is a time/memory tradeoff, distinct from matching or beating the same best comparator on both axes. This report does not promote a library default or claim broad production/frontier completion. It covers 4 fixture families, operations N, workloads refactor,rhs8 only.
- CPU and elapsed front assessments are separate. The refactor workflow includes prime factorization. The numeric CSV separates timed analysis, factor, numeric, solve and refactor calls from complete subprocess wall and CPU.
- Native labels with unavailable attributed peak bytes: ['native-amd']. Process VmHWM is valid RSS evidence and must not be called attributed factor memory. Unknown summaries have count zero and must not be plotted as zero.
- Five-peer selection differs from the older four-peer diagnostic. No cross-checkpoint speedup at different sizes or host loads is inferred.

Calibration rejections:

- drift-complex, level 128, umfpack, 4 threads, trial 0: status ok, process wall 0.747116s.
- Pre-launch one-minute host load: 1.501–2.687. Affinity and controller quiescence retain ordinary resident services; this is not an exclusive host.

Reproduction and provenance:

- `/home/ert/code/sparse-dev/benchmarks/2026-10-16-native-families/openblas-drifts`: raw SHA256 `48fba272d2509e4fde40b5ce78fd4f6855f022ac61a237038fea9adcee911fe4`.
- `/home/ert/code/sparse-dev/benchmarks/2026-10-16-native-families/openblas-chains`: raw SHA256 `d3d4bf3f169a5dd38a946898f20a7d04cde2bddaef2a147bc377b334d3ea6572`.
- Calibration selection SHA256 `cb39c751f3832a279e3e378ef33b47fc3cf9c714ee5128e3aedf380e790a42dc`; all rejected observations remain in its sibling `attempts.jsonl`.
