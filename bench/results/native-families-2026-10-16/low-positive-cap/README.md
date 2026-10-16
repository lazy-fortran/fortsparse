# Positive-cap LOW: separate ordinary checkpoint

Source `cac66561dce6178e04804da4af9c95c6fd6b19a3`; 40/40 ordinary native LOW trials passed, with five trials per family/workload/variant. CPU0, one requested and actual native factor thread. Real 16k uses 16 MiB (16,777,216 bytes); complex 4k uses 5 MiB (5,242,880 bytes). Every reported owned-payload peak is at or below its requested cap.

| Family | Workload | Variant | Wall median [min,max] ms | Completed CPU median ms | Max owned payload bytes | Cap bytes |
|---|---|---|---:|---:|---:|---:|
| drift-complex | refactor | heldout | 29.058 [28.486,31.767] | 28.056 | 3,031,240 | 5,242,880 |
| drift-complex | refactor | training | 28.627 [28.511,31.482] | 27.747 | 3,031,240 | 5,242,880 |
| drift-complex | rhs8 | heldout | 24.444 [24.272,26.235] | 23.465 | 4,327,080 | 5,242,880 |
| drift-complex | rhs8 | training | 25.376 [24.746,26.212] | 23.982 | 4,327,080 | 5,242,880 |
| drift-real | refactor | heldout | 87.371 [86.034,90.405] | 85.652 | 10,677,196 | 16,777,216 |
| drift-real | refactor | training | 87.147 [85.519,104.930] | 85.421 | 10,677,196 | 16,777,216 |
| drift-real | rhs8 | heldout | 72.151 [70.715,78.084] | 70.517 | 15,187,408 | 16,777,216 |
| drift-real | rhs8 | training | 72.108 [69.944,75.078] | 70.546 | 15,187,408 | 16,777,216 |

This is a separate source/binary checkpoint from the 240-trial 7e1f123/2f0cb01 uncapped-policy A/B. Its original 100 baseline passes and 20 timeouts are unchanged. The added LOW cap CLI is represented by the source cac6656 binary identity; no capped-LOW result is attributed to either older snapshot. This is native-only cap feasibility evidence, not an external-frontier comparison or matched cap-overhead A/B.

The caps cover solver-owned payload, not whole-process RSS. Process VmHWM, exact owned-current/peak diagnostics, phase times, completed wait4 CPU and numerical errors remain separate in `summary.csv`. All 40 raw policy/limit/peak records were checked directly. Dense80 instrumented retry diagnostics remain distinct from these ordinary measurements.

The controller reports this campaign began 13.289 s after a separate 80-trial A/B completed, without overlap. This is retained as a controller scheduling qualification, not reconstructed from benchmark elapsed times.

`recipe.json` pins raw/manifest, ELF and runtime-manifest hashes and the externally retained acquisition `run.py`. The acquisition recipe was not executed for reporting. Raw receipts and binaries stay outside Git. Existing LOW A/B PNGs describe only the earlier 240-trial campaign; they have not been relabelled as positive-cap data.
