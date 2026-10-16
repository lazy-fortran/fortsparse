# FortSparse roadmap

Active plan, 2026-10-16. Checked items describe the implemented baseline;
unchecked items are planned capabilities, not implementation claims.

## Goal and scope

Develop native sparse solvers alongside the retained external backends. Admit
native paths case by case when independent correctness, accuracy, elapsed-time
and memory evidence justify replacing an existing backend for that workload.
Keep external solvers available for capabilities or performance not yet covered.

The target is the fastest validated total solve at the requested accuracy and
memory budget on a declared corpus and hardware. Time and minimum memory can
conflict: retain the measured Pareto frontier and explicit per-case limits.
Report native performance separately from performance obtained by fallback.

Execution is **single node, single process**. OpenMP, NUMA-aware CPU execution
and native CUDA are part of the first production slice. MPI, distributed
ordering and multi-node execution are out of scope. Multiple GPUs within one
node are a later, evidence-driven milestone. Comparison backends use their
sequential or shared-memory configurations.

Use the smallest complete implementation for each admitted case. Keep state
explicit and reentrant, and expose native and external paths through the same
status-bearing solve contracts.

## Implemented baseline

- [x] Real/complex CSC storage and construction, sparse products, status API,
  backend dispatch and solver lifetime rules.
- [x] In-process SuperLU and a separate GPL UMFPACK helper.
- [x] Optional MUMPS backend, retained factors, block solves, refinement,
  pattern reuse and solver diagnostics.
- [x] Real/complex equation-level solve JVP/VJP products with independent
  numerical and adjoint checks.
- [x] Backend selection and BLAS discovery, including the fetched-SuperLU
  discovery fix verified through KIN6D.

- [x] Native real SPD simplicial reference and blocked factorization, optional
  AMD/METIS ordering, retained analysis, block RHS and implicit derivatives.
- [x] Retained FortSym primal/JVP/VJP leaves with generated CPU/CUDA schedules;
  actual-device correctness checks and context reuse.
- [x] OpenMP destination tasks grouped by complete elimination subtrees,
  compact owning factor storage and explicit actual-team statistics.

These implementation checks do not establish performance admission. The frozen
corpus, external baselines and current per-case results live in `bench/`.

The [2026-10-16 CPU MVP report](bench/results/2026-10-16-mvp.md) now establishes
robust observed-range nondominance for all 115 frozen case/workload groups on
the declared host, separately for elapsed/RSS and whole-workflow CPU/RSS.
The stronger same-comparator gate passes 109 groups; four remain indeterminate
and two are nondominated tradeoffs. GPU correctness is verified, but its
performance is not admitted. These are case-specific measured configurations,
not universal or AUTO-selection claims.

## One mathematical source, several execution schedules

Follow the KIN6D separation between mathematics, indexed operator structure
and hardware scheduling:

| Layer | Owner and retained source |
|-------|---------------------------|
| Mathematical definitions | FortSym definitions of project-owned numerical leaves: panels, block pivots, Schur updates, scaling, residuals and their derivatives. |
| Indexed operators | Contraction, reduction, gather/scatter domains and block shapes; do not expand an entire sparse factorization into scalar expressions. |
| Execution schedules | Bounded generated CPU SIMD/OpenMP and CUDA variants, including appropriate BLAS/LAPACK lowering for large contractions. |
| Sparse orchestration | FortSparse Fortran code for ordering, elimination trees, pivot decisions, task dependencies, memory management and backend selection. |

Sparse symbolic analysis of graph structure and fill is distinct from FortSym
algebra. Keep both explicit. FortSparse owns its generators; add missing
generic symbolic/compiler capabilities to FortSym only for a concrete need.

Retain definitions, assumptions and generation recipes. Generate primal,
JVP and VJP from the same mathematical source, with enclosure evaluation where
supported; add HVP only for a real consumer. Numerical probes complement
algebraic identities and independent references; they do not constitute proof.
Change definitions or schedules rather than hand-editing emitted kernels.

## Matrix, derivative and error contracts

Start with real SPD matrices from an explicitly justified consumer, such as a
fixed-boundary linear Grad–Shafranov discretization with Dirichlet conditions.
A nonlinear or free-boundary equilibrium is not automatically SPD. Specify
storage, permutations, duplicates, index range and singular/nonpositive-pivot
behavior before admitting the path. Distinguish complex Hermitian and complex
transpose-symmetric systems; symmetry alone never implies positive definiteness.

Differentiate the implicit equation, not factorization iterations:
`A dx = db - dA x`. The adjoint uses the appropriate transpose or conjugate
transpose and a declared real pairing for complex variables. Preserve sparse
matrix cotangents, factor reuse and block RHS support. Native derivative paths
must support the required transpose solves without a dense global Jacobian.

Check residuals against the original matrix independently of the fast path.
Define rounding, overflow, scaling and refinement behavior. Distinguish a
backward-error estimate, a conditioning-dependent forward bound, a validated
discrete enclosure and any downstream continuous-model certificate. Return
explicit unavailable/failed status when a claim is unsupported. These contracts
belong in the first slice. FortNum owns reusable verified arithmetic.

Mixed precision, fast math and later compression require an explicit accuracy
budget and independent checking. A successful factorization or small computed
residual alone is not a rigorous forward-error certificate.

## Execution and memory design

Use compact elimination/supernode metadata, explicit factor/update maps and
bounded reusable workspaces. Each live block has an authoritative owner.
CPU/GPU layouts may differ where measurements justify them; account for all
copies, pinned buffers, allocator pools and peak live update storage.

Use a persistent OpenMP task region with dependencies, coarsen small subtrees,
place work and storage with NUMA locality, and budget BLAS threads together
with task parallelism. Avoid accidental nested oversubscription.

Keep CUDA optional behind ISO_C_BINDING and a thin C ABI. Use reusable device
storage, streams/events and explicit asynchronous dependencies. A CPU build
must not require CUDA. Choose CPU/GPU placement from measured front shapes,
available parallelism, transfers and memory capacity. Verify factorization,
refactorization, RHS batching and transpose solves on an actual device.

Allocate during setup or bounded workspace growth, never inside hot numerical
kernels. SPD analysis sizes the required structures. Later indefinite paths
must account explicitly for delayed pivots and factor growth: enforce a memory
budget, repair task dependencies as needed, and return a clear failure/fallback
status. Reuse keys must cover every structural assumption on which reuse rests.

## Benchmark and admission rules

Freeze a small correctness corpus, a representative performance corpus and
held-out cases before optimization. Include:

- SPD meshes with anisotropy/scaling and varied ordering/fill; low-fill bands
  and fronts large enough to exercise dense kernels.
- Indefinite and adversarial block-pivot cases, complex and unsymmetric cases
  as those matrix classes enter the roadmap.
- Repeated patterns with changing values; one/many RHS; primal, tangent and
  adjoint workloads; cold setup and realistic reuse counts.
- Singular/invalid inputs, ill-conditioned cases, memory-budget failures and
  explicit error/certificate status checks.

Use named consumers such as KIN6D equilibrium/MHD, Fourier-coupled systems and
NEO-2 banded workloads when fixtures are available, plus public sparse matrices.
Check each fixture's mathematical class rather than inferring it from its name.
Compare with appropriate available CHOLMOD, SuperLU/SuperLU_MT, MUMPS, UMFPACK,
KLU, PARDISO, SPRAL and cuDSS configurations subject to access and licensing.
Do not implement every adapter before the first useful comparison.

Record ordering, conversion, analysis, factorization, solve, refinement,
transfers, derivative and verification costs, including amortized reuse.
Measure peak host/device memory, fill, allocations/synchronizations and achieved
accuracy. Use both common-ordering diagnostics and each solver's tuned whole
workflow. Retain raw trials, hardware, compiler/library revisions and settings.

Compare against the best known fixed-input baseline, not merely the previous
commit. Predeclare per-case accuracy, time and memory gates with measured noise
allowances. Aggregate speedups must not hide a failed case or a memory
regression. Record fallback explicitly; universal leadership is not a milestone.

Run focused correctness tests and dedicated performance fixtures. Use affected
named CTest tests with bounded parallelism and the supported fpm/fo path where
applicable; do not make a full-suite run a development ritual. GPU claims need
real GPU evidence; an unavailable runner is a reported validation limit.

## Microoptimization loop

Optimize in this order: mathematical representation and complexity, algorithm,
data movement, parallel work, schedule, then instructions.

1. Collect actual panel/front/RHS shapes and identify the measured bottleneck.
2. Estimate unavoidable bytes, arithmetic intensity, launch costs and critical
   path; calibrate against measured bandwidth and compute throughput.
3. Generate a bounded ahead-of-time variant set for tile size, layout, unroll,
   fusion and reduction strategy. Avoid an unbounded tuning framework.
4. Inspect vectorization reports and assembly; on GPUs inspect instructions,
   occupancy, registers and spills. Fewer symbolic operations need not be faster.
5. Recheck independent value/derivative oracles and arithmetic budgets, then
   measure microkernel, factorization, full solve and consumer workloads.
6. Retain variants only for reproducible end-to-end benefits, with hardware and
   selection provenance. Bound compile time and generated code size.

Keep a simple reference path. Put reusable optimizations in the owning library
and change generated schedules rather than maintaining divergent handwritten
mathematical implementations.

## Milestones and exit evidence

- [x] **M5.0 — Baselines and contracts.** Freeze the first real SPD fixtures,
  matrix/derivative/error contracts, raw external-solver baselines and memory
  accounting. Name the consumer and independent correctness oracle for each.
- [ ] **M5.1 — Generated numerical leaves.** Retain FortSym panel/update and
  residual definitions; generate CPU and CUDA schedules with reusable workspace.
  Exit with reproducible generation, independent value/JVP/VJP checks, actual
  GPU results and matched microbenchmarks for the measured shape distribution.
- [x] **M5.2 — Native SPD reference.** Implement a simple simplicial real SPD
  analysis/factor/solve path, permutation handling and nonpositive-pivot
  rejection. Verify factor reuse, block RHS, derivatives and error statuses
  against independent small-matrix and retained-backend oracles.
- [x] **M5.3 — Native production SPD, CPU MVP.** Add blocked/supernodal factorization,
  generated kernels, OpenMP task scheduling and CPU/GPU placement. Admit cases
  only after full-workflow time, accuracy and peak-memory gates pass, including
  setup, transfers, reuse, requested derivatives and checking.
  CPU admission is established on the declared corpus. The optional native
  GPU engine is verified and measured but remains outside performance admission.
- [ ] **M5.4 — Additional matrix classes.** Admit complex Hermitian positive
  definite and robust symmetric-indefinite paths as separate slices. The latter
  needs stable 1x1/2x2 pivots, delayed-pivot memory handling and adversarial
  stability tests. Add unsymmetric LU only for a concrete consumer with its own
  stability, transpose-solve, derivative and performance gates.
- [ ] **M5.G — Matched GPU admission.** Compare native CUDA with cuDSS and
  standalone CHOLMOD CUDA on the frozen FP64 SPD workflows in
  `bench/GPU_BENCHMARKS.md`. Establish actual offload, complete host/device
  memory accounting and external-only size calibration before native tuning.
  Profile losing cases; improve generic emission upstream and regenerate the
  consumer. Admit host and device-resident workflows separately, preserving
  the CPU frontier and all failed or unavailable comparisons.
- [ ] **M5.5 — Measured extensions within one node.** Consider multiple GPUs
  and approximate/compressed factors or preconditioners only for workloads
  that justify them, with explicit memory and accuracy budgets.

Develop M5.1 and the SPD reference together so a CPU-only prototype does not
dictate the production architecture. A slow reference is a verifier, not a
performance admission. The real-SPD CPU MVP exit is established on the frozen
corpus, including HB/bcsstk09. Further work is driven by a concrete consumer:
many-RHS batching and stronger thread scaling can improve existing tradeoffs;
GPU work needs a justified crossover case. Wider standalone microbenchmark
coverage in M5.1 remains open. Explicit AMD/ND and thread configurations remain
independent measurements; a case-specific tuned win does not establish that
AUTO or the default ordering always wins.

## Ownership and retained work

FortSparse owns sparse formats, graph analysis, factorization, sparse products,
solve diagnostics and their execution. FortSym owns generic symbolic generation;
FortNum owns generic dense/iterative numerics and verified arithmetic. FortFEM,
KIN6D and FortML retain their mathematical assembly and application semantics;
FortOpt owns optimization. Fix generic missing capabilities upstream rather
than duplicating them in a consumer.

CSR/ELL views, batched sparse products and explicit-pattern compact-support
assembly remain consumer-driven work. They must not force a change to the CSC
direct-solve API. Kernel formulas stay in FortML; PyTorch/KeOps comparisons
belong to matched FortML workloads. Stabilize/version public APIs for actual
consumers. Enzyme may be a derivative oracle, never a correctness dependency.

## Provenance and licensing

Retain source revisions, mathematical definitions, generator inputs, license
checks, numeric fixtures and benchmark recipes. Preserve the MIT core and the
existing separate GPL UMFPACK helper boundary; inspect each added dependency's
terms before choosing an integration boundary. Public algorithm descriptions
do not grant permission to copy incompatible implementations.

The local research collection is `~/code/sparse-dev/`: repositories, downloaded
papers, library cards, availability records and a report. It is supporting
provenance, not a build dependency. Keep generated PDFs, plots and large
benchmark data out of source history. Starting references:

- [KIN6D symbolic-generation policy](https://github.com/itpplasma/kin6d/blob/73891c54a71c3a080826043fd02924779a6e79a3/docs/tooling.md).
- [KIN6D performance policy](https://github.com/itpplasma/kin6d/blob/73891c54a71c3a080826043fd02924779a6e79a3/docs/performance.md).
- [Olaf Schenk's ETH thesis record](https://doi.org/10.3929/ethz-a-003876213),
  *Scalable Parallel Sparse LU Factorization Methods on Shared Memory
  Multiprocessors* (2000); retained PDF and metadata in the local collection.
- [SPRAL documentation](https://ralna.github.io/spral/), for a public
  shared-memory solver reference alongside the other retained solver records.
