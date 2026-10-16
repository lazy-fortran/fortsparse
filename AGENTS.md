# FortSparse agent rules

## Active delivery

Implement native single-node coverage and bounded-memory execution on `main`,
using `ROADMAP.md` as the authoritative work order. Prioritize feature completeness
and correctness before further performance tuning. Deliver native-only consumer
branches and draft PRs with affected-test evidence and explicit unsupported paths. Read README.md, ROADMAP.md
and the affected interface before changing behavior. Keep external backends
available. Reach a usable GPU SPD checkpoint, then prioritize CPU matrix classes;
full GPU frontier admission remains required later. MPI, multi-GPU execution
and disk-backed factors are outside this plan.

After functional consumer coverage and correctness, resume profiling and optimizing
unresolved cases until native execution is on the measured time/memory Pareto front. A working factorization, a microkernel
win or a fallback to another library is not completion. Do not claim a novel
algorithm merely because the implementation and generated schedules are new.

Commit and push complete, locally verified increments to `main`
regularly. Publish performance numbers after implementation milestones and
material optimizations. Preserve unrelated changes, stage explicit paths and
keep worker integration under one controller. Remove task worktrees only after
their work is integrated and their owned processes have stopped.

Every consumer draft PR includes matched baseline/native CPU time and memory
numbers for the exercised tests and application runs, with absolute and relative
changes. Record source/dependency pins, compiler flags, input identity, affinity,
thread counts, repetitions and correctness status. New tests without an equivalent
baseline and blocked comparisons stay explicitly unavailable. Long application
correctness runs are separate from the sub-second solver benchmark corpus; label
unreplicated observations and shared-host limits, and do not infer Pareto admission.
Use spawn/exec-aware observers: Python fork and observed Python spawn paths can
impose an RSS high-water floor on a tiny executable. Qualify the observer with
a small C spawn/wait helper or actual post-exec memory evidence. Retain per-process evidence and
never present separate process peaks as aggregate workflow memory.

## Mathematics and implementation

- Retain project-owned continuous mathematics in FortSym generators. Generate
  primal/JVP/VJP together; keep emitted CPU/CUDA sources reproducible. Change
  definitions or schedules rather than hand-editing generated output.
- Keep graph analysis, indexed operator structure, hardware schedules and
  backend orchestration separate. Use BLAS/LAPACK for suitable contractions.
- Differentiate the solve equation, not factorization iterations. Preserve
  sparse cotangents, factor reuse and independent residual evaluation.
- Use explicit, reentrant state and bounded reusable workspaces. Numerical
  hot paths allocate nothing. Include padding, task metadata, pinned storage,
  allocator pools and CPU/GPU copies in memory accounting.
- Enforce explicit solver-owned host/device allocation limits before allocation,
  including analysis, transient growth, reuse overlap and derivative work.
  LOW_MEMORY trades concurrency/workspace for RAM without changing accuracy.
  Mark predictions as layout sizes, bounds or estimates; never equate an
  estimate-only vendor gate with a hard cap. Account process overhead separately.
- OpenMP and optional native CUDA are production execution paths. Control
  task/BLAS thread budgets, NUMA placement and device residency. Record actual
  active threads as well as requested threads. CPU builds require no CUDA.
- Inspect licensing and matching header/library ABIs before adding wrappers.
  Preserve the MIT core and separate GPL UMFPACK helper. Fix generic numerical
  or tooling defects in their owning library when they block this work.

## Measurement contract

- Freeze correctness, representative performance and held-out fixtures before
  tuning. Select scaling sizes using external baselines alone. Never shrink a
  failing native case, remove a difficult case or weaken accuracy to manufacture
  a win. Retain failures, timeouts and unavailable comparisons.
- Ordinary tests and benchmark invocations must finish below one second on
  the declared reference machine. Use bounded replay for useful sampling of
  short kernels. Identify profiler overhead separately; instrumented runs are
  diagnostic evidence, never reported as ordinary performance measurements.
- Compare SuperLU, UMFPACK and explicitly named Intel/Panua PARDISO products;
  include MUMPS and other relevant available solvers. Keep incompatible BLAS
  runtimes in separate executables and verify dependency resolution.
- Measure cold numerical-runtime and warm-runtime workflows, analysis,
  factorization, refactorization, one/many RHS, derivatives, checking and GPU
  transfers. State reuse counts. Do not confuse a fresh process with a cold
  filesystem cache or omit setup when reporting total cost.
- Report wall time AND process CPU time (user plus system), including helper
  CPU separately. Full-lifetime CPU needs completed-process user/system counters.
  An in-process cumulative sample ends at that sample; a timer started in `main`
  also excludes loader/startup. Label those intervals and preserve old receipts.
  On Linux, `wait4` includes reaped-child CPU: never add helper CPU twice.
  Sweep 1, 2, 4, 8 and 16 physical cores where available, with
  fixed affinity, explicit thread settings and controlled background load.
  Report speedup, efficiency, actual teams and serial/overhead-limited cases.
- Report host/device memory and factor fill. Distinguish exact process peak
  RSS, sampled/PSS measurements and attributed storage/workspace bytes.
  Report UMFPACK parent/helper memory separately; summing RSS double-counts
  shared mappings. Unknown or incomparable memory is not evidence of a win.
  A helper's own post-exec Linux VmHWM is a safe individual-process lower
  bound for the process-tree resident union peak. Inherited child ru_maxrss
  is unsuitable. Validate the requested diagnostic and each lifetime; malformed,
  missing or ambiguous receipts stay unavailable. Never call a sum or maximum
  of separate process peaks the aggregate peak.
- Measure repeated factor/free lifetimes as well as retained-factor reuse. Lower
  live allocation bytes do not establish lower RSS: inspect allocator caching
  and validate actual page release with context-local workspace ownership.
  Distinguish cumulative allocation/copy traffic from peak live storage.
  A lower payload cap or fewer live bytes does not guarantee a smaller process.
- Use five fresh-process timing trials by default; retain raw numeric evidence,
  medians/ranges, commands, source/binary identities, hardware and tool/library
  versions. Keep large traces and datasets outside source history, retaining
  compact results and reproduction recipes in the repository.
- Compare each declared case against the best measured external time/memory
  frontier at matched accuracy, using each solver's tuned workflow. Retain
  default and accuracy-qualified tuned controls as distinct configurations
  (including UMFPACK refinement steps); record runtime controls and ABI. Provide
  common-ordering diagnostics where supported and label their limits.
  Match-or-beat requires time and memory against the SAME frontier comparator;
  overlap within observed variability is indeterminate. Unknown-memory solvers
  only obstruct frontier membership when their time could actually dominate.
  Aggregate speedups must not conceal regressions or failures.
- Distinguish Pareto membership from the stronger match/beat gate. A native
  time/memory tradeoff can be nondominated without beating an existing point
  on both axes. Keep the strict same-comparator result unchanged and report
  nondominance separately. Robust observed-range membership must exclude EVERY
  measured external configuration as a possible dominator by strict time OR
  memory separation. Unknown aggregate memory, missing comparators and range
  overlaps remain explicit; median-front membership is a weaker claim.
- GPU admission compares native CUDA with actual GPU-library execution at the
  same FP64 accuracy, device and host-thread budget. Follow
  `bench/GPU_BENCHMARKS.md`; CPU-versus-GPU plots cannot establish GPU-library
  competitiveness. Record CPU fallback, unavailable offload and unknown memory.
- Generic emission defects are repaired in the owning generator (FortGen,
  through FortSym where appropriate), then regenerated in FortSparse. Retain
  compiled-kernel evidence as well as independent value/JVP/VJP checks.

## Profile-driven optimization loop

1. Reproduce a losing case on the current revision and identify whether setup,
   ordering, allocation, numerical work, solve, checking or transfers dominate.
   Profile cold and reused paths separately.
2. Use CPU sampling/counters (`perf`) and allocation/call tracing when useful.
   Inspect call stacks, cycles, instructions, cache misses, bandwidth, branch
   behavior, synchronization, task granularity and allocator traffic. If hardware
   counters are unavailable, record that limit and use supported user-space
   profiling; do not silently change system security settings.
3. On GPUs use Nsight Systems/Compute or supported equivalents for launch,
   copy, synchronization and kernel timelines. Inspect FP64 throughput, memory
   traffic, occupancy, register pressure and spills. Count transfers and live
   storage, including retained workspaces and vendor allocations where measurable.
4. Estimate algorithmic complexity, unavoidable bytes, arithmetic intensity and
   the critical path. Optimize representation/algorithm, data movement and
   parallel work before instruction details. Inspect vectorization reports,
   assembly and GPU instructions for the measured hot kernels.
   For sparse fill regressions, inspect input/factor nnz, allocated capacity,
   zero diagonals, pivot departures and singleton removal. Compare matching and
   ordering before tuning scatter instructions. Fewer factor entries need not
   mean lower payload when geometric capacities change.
   Count cumulative scatter entries and zero-multiplier work: final factor nnz
   alone does not describe elimination cost. Identify mesh, nonlinear iteration
   and factor ordinal, including preconditioner replacements. A sampled hotspot
   cannot distinguish excessive fill from repeated updates or expensive entries.
   For masked panels, retain active-lane histograms weighted by row entries
   and arithmetic, not only call counts. Check packed instructions in the
   compiled release kernel; a SIMD directive does not establish vectorization.
   Distinguish inlined DWARF frames from actual calls: inspect disassembly
   before attributing a sampled symbolic leaf to per-entry call overhead.
5. Make one bounded change addressing measured evidence. Generate a bounded
   ahead-of-time variant set for tiles/layouts/unrolling/fusion where justified.
   Fewer symbolic operations do not imply lower elapsed time. Avoid unbounded
   autotuning, speculative frameworks and unrelated cleanup.
6. Recheck independent values, derivative identities and error semantics, then
   measure microkernel, factorization, complete solve and consumer impact.
   Reject changes that merely move time or memory outside the measurement.
7. Retain a specialization only for reproducible end-to-end gains without
   violating the case's memory/accuracy contract. Update evidence, commit and
   push the verified checkpoint, then profile the remaining losses. Stop
   optimizing a case once the declared frontier criterion is satisfied unless
   another unresolved case or a regression justifies more work.

Keep the original matrix fixed across ordering/matching experiments. Record
permutations and RHS/solution transformations, and check residuals against the
original matrix with an independent oracle. Capture manufactured inputs with
source pins, exact export patches and hashes outside Git. Never weaken a
scientific gate to obtain a smaller benchmark. Preserve failed exporters and
deadline-censored attempts; a timed-out factor has no completed fill/solve
result. Longer profiling diagnostics are explicitly excluded from the
sub-second benchmark admission corpus.

For accuracy regressions, capture the original matrix, true RHS and rounded
factors before changing arithmetic. Isolate factorization, triangular solves,
scaling, residual evaluation and reconstruction. Compare the requested
observable as well as normwise errors: cancellation can amplify tiny solution
differences. Separate solve error from the rounding of the final projection.
A higher-precision residual-corrected reference is numerical evidence, not a
rigorous enclosure. Preserve existing gates and report their failures; do not
equate agreement with a historical floating result with mathematical truth.

Check the surrounding algorithm before blaming the sparse solver. Record
whether outer iterations actually converge. Verify projected eigensystems with
independent known eigenvalues and true eigenvector residuals, including the
final projected column and happy breakdown. Apply shared correctness repairs
identically to baseline and candidate, retaining the original failed evidence.
For affine operator callbacks, compare the offset scale with the normalized
probe before subtracting F(0). Verify any scaling repair on independent spectra
and original-operator residuals, including large offsets and finite-range
failures; preserve the callback contract and physical case controls.

## Verification

Use the Fortran skill, the project's native CMake/CTest and fpm authorities,
and one controller-owned resident fo Gremlin lane where supported. Run affected
named tests with bounded parallelism, never the full suite as an edit ritual.
Avoid profiling/benchmarking while builds or other heavy campaigns compete.

Use independent known solutions, small dense references, structural/adjoint
identities and external libraries as behavioral oracles. Cover permutations,
disconnected blocks, failed pivots, invalid inputs, block RHS, reuse, derivatives
and CPU/GPU agreement. GPU claims need actual-device evidence. Keep backward
error, forward bounds and rigorous certificates distinct; return or document
unavailable/failed status for unsupported claims. Report verification limits.
