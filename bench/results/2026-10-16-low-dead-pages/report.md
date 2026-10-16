# Real LOW dead-page release: current-source checkpoint

The real-only 1 MiB page-release port reduces observed process RSS at a CPU/time cost on two frozen LOW workloads. This is a bounded native-before/after comparison, not a frontier or production-class admission. All 20 ordinary children passed; no failures, retries or size changes occurred. Two separate instrumented children prove the actual helper calls and are excluded from every timing, CPU and RSS aggregate.

| Frozen drift128 workload | Whole-process wall ms, baseline → port | Completed-process CPU ms | Peak RSS MiB |
|---|---:|---:|---:|
| rhs8 | 105.922 → 121.640 (+14.84%) | 101.579 → 120.586 (+18.71%) | 32.133 → 27.625 (-14.03%) |
| refactor | 114.605 → 137.044 (+19.58%) | 113.755 → 136.113 (+19.65%) | 26.344 → 23.180 (-12.01%) |

Each point has five fresh-process trials. `summary.csv` retains every observed minimum/maximum and factor/solve/refactor phases separately. Port RSS ranges lie strictly below baseline ranges in both groups; its whole wall and CPU ranges lie strictly above baseline. Refactor workflow CPU includes initial factorization, retained changed-value refactor, solve, checking and exit; it is not the refactor phase CPU.

All ten paired trials have identical factor entries (762644), backward/known-x forward errors and reuse flags. Reported final current payload is 10742732 B for both builds. Reported operation peak is 15318480 B for RHS8 and 10808268 B after refactor, unchanged by the port. The refactor statistic describes the latest factor operation; it is not the maximum payload of the complete workflow. These uncapped controls request host cap zero and LOW policy one; existing cap/fault correctness tests separately verify the unchanged quota semantics.

The independent 2×2/3×3 and large real/complex lifetime oracles, page guards, denied-advice fallback and all 25 discovered allocation-fault positions per scalar type pass. Copied-factor checks concern allocatable factor storage. The public `sparse_vector` pointer-copy ownership defect remains separately tracked; this is not whole-context copy coverage.

Both instrumented fixture invocations executed 634 helper entries, 50 eligible payload entries and 50 actual helper-origin MADV_DONTNEED calls, with 88064000 cumulative advised bytes each. There were zero syscall failures, invalid interior spans or unrelated wrapped calls. Eligibility checks the original payload size: boundary trimming can make an advised span smaller than 1 MiB. The cumulative byte count is not peak or retained storage. Wrappers surround the helper entry and madvise; each accepted syscall is proven strictly inside that helper payload. The diagnostic ELF links the exact candidate core archive/driver objects plus the counter object.

Baseline source is `8d7f634`; port source is `fda5cbd`. The public patch integration is `4fd03efc9397be1c80cf6f2cebab49ef614f45a8`; canonical main was not retimed under that identity. Both are native-only Release builds with AMD 3.3.4, CUDA/SuperLU/UMFPACK/MUMPS/METIS disabled and matched `-O3 -march=znver3 -DNDEBUG` flags. Resolved shared-library paths and SHA256 values match. Fixture/driver/records/process-metrics sources are identical. The original independent Fortran fixture generates n=16384, nnz=81408; no Python numerical kernel or new permutation is used.

The predeclared groups are training drift-real level128, normal operation, explicit AMD+LOW, RHS8 and refactor. Five alternating pairs per group run on CPU30 with one OpenMP and BLAS thread, a 0.9 s wall watchdog and no task-owned concurrent jobs; normal resident services remain. The two diagnostics use a separate 3 s wall watchdog, not a CPU budget. Existing external-only fixture selection is unchanged. No external sparse-factor peer is compared in this checkpoint.

The retained qualified C posix_spawn/wait4 observer reports complete child user+system CPU through exit and completed-process RSS. Its exact ELF SHA matches the prior heavy-Python-parent RSS/watchdog qualification. The observer verifies CPU30 affinity. Driver post-exec VmHWM equals observer peak RSS in all 20 ordinary children. Neither quantity represents aggregate multi-process memory; these native-only children have no sparse helper.

Raw evidence and reproduction recipes are retained outside Git at:

`/home/ert/code/sparse-dev/benchmarks/2026-10-16-native-families/low-real-deadpage-current-ab`

`trials.csv` retains all20ordinaryscalar records, `diagnostics.csv` retains both excluded call-count observations, and `provenance.json` pins the source/ELF/library/raw identities. See [recipe](recipe.md) for portable table regeneration and the bounded measurement controls.

Restore the two source pins in isolated worktrees before rebuilding. `freeze.py` retains the exact configuration/link recipe and verifies flags, input-source identities, shared-library closure and observer qualification; both CMake caches are archived. Run `run.py --go` only in a fresh archive after acquiring the measurement window. It preserves every attempt and keeps the two diagnostic receipts separate. Original raw files are intentionally not overwritten. No generated plots, ELFs or large data are committed.

| Retained artifact | SHA256 |
|---|---|
| `manifest.json` | `f4aefe1db3e2a3564eafcbef239b01b08114f49d7119d962f8594d302ddd48ca` |
| `attempts.jsonl` | `0828d17ede0c949def0e942f3f83df2c04015d6df925ea2ab7c8eeb30652b4b6` |
| `diagnostics.jsonl` | `a94270d37977f050b2e667d75f6683064961e5cd6d51ef87ff8f10e7b4766741` |
| `run.py` | `85557b4af4054b6b4598820ca79f4d5b2c2be3640e8dd191aa6825c8428f50b9` |
| `freeze.py` | `e1102c91032472a113030d79f08823f19b9291bef7430240d5a2877665caafd5` |

Frozen baseline ELF: `a11319dd732691cdb21da33a68d89bafe217f070a0972d5bd3f9014ed843c6cc`.
Frozen port ELF: `f4798066cb09414460380f4e53978adf35f61b502fc5a91742e57a1d5b57b4ef`.
Instrumented ELF: `7f93db0554a43f2e2bbd3946153269f9e73d603415efa5971b2aba0867f3b821`.
