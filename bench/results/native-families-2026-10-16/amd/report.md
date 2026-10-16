# Native sparse-family completed checkpoint

Source revision: `52dde039f95c189532533d99440d0fee3055f554`. Raw attempts SHA256: `50586c748a7cd3a429259ad29c5a2a1a928c10738a2934f3910fa1c503c636cd`.

1000 attempts: 1000 ok.
CPU interpretation: **completed workflow wait4 user+system CPU including reaped descendants; parent cumulative pre-CSV and main interval separate**. Prior raw fields are preserved.

The rows below show training, normal operation, RHS8, one requested thread. Wall is complete subprocess elapsed time; factor is a timed library call. Ranges span all successful trials. Timeout rows have no factor/RSS sample.

| Family | Backend | OK/all | Wall median [min,max] ms | Factor median ms | Workflow CPU median ms | Parent peak MiB |
|---|---|---:|---:|---:|---:|---:|
| drift-real | intel-pardiso | 5/5 | 88.043 [87.626,88.181] | 41.947 | 86.173 | 47.719 |
| drift-real | mumps | 5/5 | 41.959 [41.147,42.305] | 20.758 | 40.418 | 36.258 |
| drift-real | native-amd | 5/5 | 111.966 [110.234,114.368] | 86.307 | 109.982 | 28.844 |
| drift-real | superlu | 5/5 | 56.962 [56.213,60.084] | 31.038 | 54.870 | 40.250 |
| drift-real | umfpack | 5/5 | 58.561 [57.955,58.606] | 31.777 | 56.788 | 18.617 |
| drift-complex | intel-pardiso | 5/5 | 104.873 [104.258,105.708] | 50.124 | 102.765 | 56.676 |
| drift-complex | mumps | 5/5 | 55.820 [54.409,56.402] | 29.871 | 54.016 | 46.195 |
| drift-complex | native-amd | 5/5 | 242.407 [241.519,251.757] | 191.285 | 239.547 | 34.785 |
| drift-complex | superlu | 5/5 | 79.261 [74.990,80.014] | 46.600 | 77.367 | 57.250 |
| drift-complex | umfpack | 5/5 | 81.666 [80.079,87.879] | 48.408 | 79.850 | 18.410 |

Timeout counts (all thread budgets, variants and trials):

- None.

Interpretation and limits:

- AMD is an explicit native ordering option; default ordering is unchanged. Native AMD still takes more total wall time than the best declared external comparator. This checkpoint does not admit it to the performance frontier.
- UMFPACK here uses the helper/IPC path. A standalone direct UMFPACK comparator is still required before judging best-library performance.
- Completed workflow CPU uses wait4 user+system usage, including reaped helper CPU and termination. Helper CPU must not be added again. Parent pre-CSV cumulative CPU and main interval remain separate in the CSV.
- Native AMD attributed peak bytes are unavailable (-1), never zero. Parent VmHWM is a process RSS measurement, not attributed factor memory.
- Native records report one active factor thread at every requested budget. External active teams are unknown (-1), not evidence that they used the requested team.
- Parent VmHWM and helper RUSAGE_CHILDREN ru_maxrss remain separate. Helper RSS can include inherited/pre-exec highwater. Never sum them or use their maximum as a validated aggregate lower bound. UMFPACK aggregate memory is unknown; its parent VmHWM alone is a safe lower bound.
- Successful-only timing ranges do not replace timeout evidence. A timeout marks a censored observation at the configured deadline, not a measured factor time.
- No global Pareto/admission claim follows: this is a bounded fixture corpus, the initial checkpoint retains 150 censored timings, true helper aggregate memory is unknown, and the shared host has background load.
- Observed pre-launch one-minute host load: 0.37–2.87. Core affinity reduces placement variation; it does not eliminate shared-host interference.
