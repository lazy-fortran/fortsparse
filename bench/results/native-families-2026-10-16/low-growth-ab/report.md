# LOW growth: ordinary native A/B and separate diagnostics

Baseline `7e1f123ea21fd8412480ec6a6ddf0696893acb2a`, candidate `2f0cb01f902d5b38768402eddefbcea3de225692`. 240 interleaved ordinary attempts on CPU0: baseline 100 passed/20 timed out; candidate 120 passed. Each version/policy/case/workload/variant has five trials. Raw SHA256 `7292b8d45799c621f12e06d9a066450625aad6e73e94407948bd0181c2495832`.

Policies are AUTO without a cap, AUTO with a 64 MiB cap, and LOW without a cap. The optional positive LOW-cap campaign is not present here. Native factorization is serial. No external-frontier or thread-scaling claim follows.

LOW complete process wall milliseconds, median [minimum,maximum]:

| Family | Workload | Variant | Baseline LOW | Candidate LOW |
|---|---|---|---:|---:|
| drift-real | rhs8 | training | 5 timeouts ≥900 ms | 74.274 [73.589,75.436] |
| drift-real | rhs8 | heldout | 5 timeouts ≥900 ms | 73.841 [71.750,77.543] |
| drift-real | refactor | training | 5 timeouts ≥900 ms | 89.606 [87.306,91.499] |
| drift-real | refactor | heldout | 5 timeouts ≥900 ms | 90.555 [86.179,91.916] |
| drift-complex | rhs8 | training | 271.871 [269.982,283.114] | 24.991 [24.460,25.990] |
| drift-complex | rhs8 | heldout | 272.781 [269.888,279.471] | 25.699 [24.853,26.679] |
| drift-complex | refactor | training | 259.390 [258.442,264.513] | 29.799 [28.819,30.037] |
| drift-complex | refactor | heldout | 264.588 [257.480,278.118] | 30.557 [29.248,32.846] |

Memory and scope:

- Candidate real training RHS8 LOW: retained 10,611,660 bytes, median 33,440 KiB (32.66 MiB) RSS. AUTO: retained 13,891,596 bytes, median 29,948 KiB (29.25 MiB) RSS. Lower retained owned payload does not establish lower RSS; this LOW RSS is higher. Allocator behavior is investigated separately, not attributed here.
- `summary.csv` retains all 48 version/configuration groups, timeouts, min/median/max wall and completed wait4 CPU, factor/solve/refactor, payload, process RSS and numerical errors. Timeout CPU is incomplete and remains unavailable. No finite completed-time speedup is assigned to the real baseline LOW timeouts.
- Aggregation reads only raw receipts, never the archived summary.json. Its corrected filter required no new solver invocations. The candidate ELF also links METIS while the baseline does not; library identities are retained. Small AUTO timing differences are descriptive, not an isolated algorithm gain.
- The source and binary identities are the measured snapshots above. Later production integrations f4f329d/9354c60 do not change this campaign’s identity. AUTO+capped64MiB is not LOW+capped64MiB.

Separate instrumented profile/counter evidence:

| Quantity | Exact-growth diagnostic | Small-slack diagnostic |
|---|---:|---:|
| growth_events | 32,767 | 317 |
| cumulative_factor_prefix_copied_bytes | 54,871,019,988 | 148,759,548 |
| cumulative_factor_capacity_requested_bytes | 54,879,975,108 | 158,142,540 |
| cumulative_process_allocation_calls | 65,604 | 704 |
| cumulative_process_allocation_requested_bytes | 54,898,845,950 | 177,013,382 |
| retained_array_payload_bytes | 10,263,804 | 10,611,660 |
| phase_maximum_array_payload_bytes | 14,899,228 | 15,187,408 |
| factor_nnz | 762,644 | 762,644 |

- The 54,871,019,988 copied bytes (~54.87 decimal GB) are cumulative traffic, not peak memory. Diagnostic phase-maximum payload rises from 14,899,228 to 15,187,408 bytes; retained payload also rises slightly versus old exact LOW.
- The historical perf sample uses source 3b9384b: 740 samples, zero lost, reserve 64.99% own cycles and visit 5.23%. It is separate from the ordinary 7e1f123 A/B snapshot. Instrumented diagnostic runtimes and allocations are not ordinary performance measurements. Whole-process allocation counters also include runtime work.
- Dense80 diagnostic evidence retains one exact retry after the hinted growth attempt: success with peak 144,812 bytes at that cap, failure one byte below, required 38,880 versus available 38,879. Actual peak is max across sequential attempts, never their sum. These diagnostics are not the pending ordinary positive LOW-cap campaign. Source-only final cleanup is distinguished in the diagnostic archive from the final replay’s instrumented algorithm.
