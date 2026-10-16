# Interleaved scalar versus generated indexed LU

Raw attempts: `/home/ert/code/sparse-dev/benchmarks/2026-10-16-native-families/indexed-ab/attempts.jsonl`; SHA256 `a46fcbcfc74adad1dfd8aefda7de83ce82a1fed88c4f974245518ff7da7fcf7e`.

80/80 attempts completed successfully. Sixteen five-trial revision/configuration groups form eight paired comparisons. Each trial alternates execution order. Both binaries use explicit AMD and one requested/actual native factor thread on CPU 0; N solves, training and heldout coefficients, n=16,384. Refactor workflows include the prime factorization.

Scalar source: `52dde039f95c189532533d99440d0fee3055f554`. Indexed source: `b62e1de981dd848df3d68a63d6ef43e17433a7f6`. Binary hashes and paths are retained in `provenance.json` and `recipe.json`.

Complete subprocess wall milliseconds, median [minimum,maximum]:

| Family | Workload | Variant | Scalar | Indexed | Paired speedup |
|---|---|---|---:|---:|---:|
| drift-complex | refactor | heldout | 388.197 [383.377,392.743] | 159.364 [154.992,160.252] | 2.441 [2.417,2.505]× |
| drift-complex | refactor | training | 386.393 [384.685,388.516] | 157.483 [156.168,158.557] | 2.454 [2.430,2.488]× |
| drift-complex | rhs8 | heldout | 244.288 [239.874,245.721] | 110.860 [109.235,112.056] | 2.196 [2.180,2.217]× |
| drift-complex | rhs8 | training | 246.124 [243.713,257.114] | 111.983 [111.588,117.294] | 2.191 [2.099,2.230]× |
| drift-real | refactor | heldout | 174.018 [172.055,174.840] | 94.267 [93.474,96.435] | 1.836 [1.806,1.862]× |
| drift-real | refactor | training | 172.079 [170.310,180.578] | 94.469 [94.048,102.543] | 1.822 [1.761,1.862]× |
| drift-real | rhs8 | heldout | 110.830 [110.244,114.730] | 66.870 [66.315,69.921] | 1.665 [1.585,1.695]× |
| drift-real | rhs8 | training | 109.571 [108.968,113.236] | 66.164 [65.275,66.789] | 1.648 [1.636,1.735]× |

Completed workflow CPU milliseconds, median [minimum,maximum]:

| Family | Workload | Variant | Scalar | Indexed |
|---|---|---|---:|---:|
| drift-complex | refactor | heldout | 385.237 [380.128,389.603] | 157.179 [153.047,158.280] |
| drift-complex | refactor | training | 382.782 [381.661,385.428] | 155.618 [154.403,156.523] |
| drift-complex | rhs8 | heldout | 241.708 [237.681,243.413] | 109.124 [107.632,110.390] |
| drift-complex | rhs8 | training | 243.623 [241.549,254.384] | 110.348 [109.894,115.608] |
| drift-real | refactor | heldout | 171.969 [169.638,172.942] | 92.766 [91.974,94.833] |
| drift-real | refactor | training | 169.916 [168.391,178.467] | 92.731 [92.558,98.624] |
| drift-real | rhs8 | heldout | 108.816 [108.358,110.005] | 65.129 [64.830,68.484] |
| drift-real | rhs8 | training | 107.577 [107.377,111.623] | 64.647 [63.797,65.232] |

Scope and limits:

- This is a native scalar/indexed A/B diagnostic, not an external solver comparison or frontier admission. Sizes retain the earlier four-external calibration; this run does not recalibrate them or add direct UMFPACK evidence.
- Every indexed wall and CPU trial is below every scalar trial in its matched configuration. The retained paired ratios use matching trial IDs.
- Both versions report exactly 762,644 factor entries throughout. Process RSS ranges overlap in every matched group; no RSS reduction is established. Attributed peak bytes are unavailable in both, not zero.
- Workflow CPU is completed wait4 user+system time including reaped descendants. Helper CPU is not added again. Parent pre-CSV and main-entry CPU remain distinct. These native runs have no helper CPU/RSS.
- `summary.csv` separates factor, numeric, analysis, solve, refactor, prime-factor, main-workflow and complete subprocess phases; unknown metrics have count zero. `paired.csv` retains all eight wall/CPU ratio ranges. Both versions pass the original-A known-solution/backward-error oracle.
- Observed pre-launch one-minute load: 0.547–0.716. Core affinity and interleaving do not eliminate shared-host variability.
