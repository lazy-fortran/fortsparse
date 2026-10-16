# Private real NEO-2 aggregate benchmark

Fixture `neo2-lorentz-tag6-4856-v1`, source `c7207e78fddc35a0fb84fe50dd8ec7fe3b5fe4d6`. 150 attempts: 150 pass.

Raw private attempts SHA256 `60583a9d47f7e6b0c493089ff372505be269a7a28dcac3d1f7b88ceea64effb2`. Frozen accuracy-policy SHA256 `c03da85fc3dfae66b8c09874a69faf78d10de0d9b368b18d6b91813046a42d54`; manifest SHA256 `08741c8f0261d20ff2d04b5a2c87706c9f6c9aebe5d6047304656a5ffdfa4516`.

One real NEO-2 matrix, n=4,856, nnz=28,898. Four RHS contain three physical columns (one duplicate) and one manufactured known-solution column. Five fresh trials for each backend and requested budget 1/2/4/8/16; all five external peers qualified with the same frozen executable set before native measurement. No matrix, RHS, candidate solution or per-RHS numeric audit is published.

One requested thread: successful child-process medians [min,max]; timed factor/solve medians are separate library calls. RSS is process VmHWM.

| Backend | Child wall ms | Completed child CPU ms | Factor ms | Four-RHS solve ms | Parent peak MiB |
|---|---:|---:|---:|---:|---:|
| intel-pardiso | 100.458 [99.965,101.877] | 98.759 [98.255,100.140] | 15.245 | 4.955 | 27.453 [27.395,27.574] |
| mumps | 87.796 [86.191,88.532] | 86.173 [84.577,86.853] | 8.231 | 1.034 | 18.488 [18.379,19.277] |
| native-amd | 88.360 [87.927,89.957] | 86.716 [86.290,88.428] | 10.082 | 0.727 | 11.176 [11.113,11.438] |
| superlu | 89.387 [89.007,95.771] | 87.748 [87.415,91.803] | 9.850 | 1.483 | 16.305 [15.520,18.402] |
| umfpack | 92.769 [92.421,94.589] | 91.152 [90.645,92.221] | 13.290 | 1.283 | 8.961 [8.832,9.133] |
| umfpack-direct | 93.751 [92.998,101.419] | 92.090 [91.324,99.608] | 10.075 | 3.672 | 16.164 [16.090,16.199] |

Observed-range child-time/RSS assessments against every external peer/thread configuration:

- process_wall_s: 3 robust_nondominated, 2 indeterminate.
- process_tree_cpu_s: 3 robust_nondominated, 2 indeterminate.

Exact benchmark audit and timing limits:

- intel-pardiso: per-budget median audit+reader wall 0.539–0.552s; maximum exact-audit normwise backward diagnostic 9.617e-17, row-scaled 2.401e-16.
- mumps: per-budget median audit+reader wall 0.543–0.558s; maximum exact-audit normwise backward diagnostic 1.200e-16, row-scaled 5.264e-16.
- native-amd: per-budget median audit+reader wall 0.544–0.554s; maximum exact-audit normwise backward diagnostic 1.430e-16, row-scaled 5.188e-16.
- superlu: per-budget median audit+reader wall 0.545–0.555s; maximum exact-audit normwise backward diagnostic 1.527e-16, row-scaled 2.464e-13.
- umfpack: per-budget median audit+reader wall 0.543–0.555s; maximum exact-audit normwise backward diagnostic 2.428e-16, row-scaled 4.693e-16.
- umfpack-direct: per-budget median audit+reader wall 0.544–0.551s; maximum exact-audit normwise backward diagnostic 4.722e-17, row-scaled 1.542e-16.

- Completed child CPU uses wait4 user+system usage including reaped helpers; helper CPU must not be added again. Child wall includes input/reference I/O, original-A checks, candidate export and teardown. Factor/solve calls occupy only part of it; the remainder is not an isolated I/O measurement.
- Exact rational residual auditing runs separately after each child. It is a benchmark correctness oracle and empirical forward regression, not a rigorous forward-error or continuous certificate. Child timings exclude this audit. Audit CPU and peak memory were not measured. Child RSS is not audit-inclusive application memory. The child+audit wall sum excludes prerequisite audits, hashing, dependency inspection and other orchestration; it is not whole-application elapsed time.
- The CSV derives requested-budget wall speedup from each backend’s own one-thread median divided by its requested-budget median; efficiency divides that speedup by the requested budget. These exclude exact-audit costs. Native factor teams are one and external active teams are unknown in this checkpoint: these are requested-budget ratios, not measured active-team efficiencies.
- Observed-range nondominance uses strict native max time < external min time OR native max RSS < external min RSS against all 25 external configurations. IPC aggregate RSS is unknown; only parent VmHWM is a safe lower bound for excluding a potential dominator. Helper lifecycle ru_maxrss is separate and may include inherited highwater. IPC has no aggregate point in the figures.
- Nondominance is distinct from matching/beating the same comparator on both axes. No broad matrix-class admission, certificate-inclusive timing or whole application result follows from this one frozen private fixture.
- Every wrapper, child and exact-audit status is preserved in aggregate-safe `attempt-statuses.csv`; failed trial counts are retained. Empty numeric summaries with count zero mean unavailable, never zero-valued performance.
- Pre-launch one-minute host load 1.377–1.846; controller quiescence and fixed affinity do not make the host exclusive.
