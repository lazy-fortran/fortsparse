# FortSparse roadmap

Active plan, 2026-10-16. Checked items describe the implemented baseline;
unchecked items are planned capabilities, not implementation claims.

## Goal and scope

Deliver native sparse factorization and solves for the audited single-node
requirements of ITPplasma codes. Retain external backends as explicit alternatives
and independent comparison oracles. BLAS/LAPACK and graph-ordering dependencies
are permitted; an external sparse factorization is not native execution.

The release requirement is measured time/memory Pareto-front membership for
every declared representative workload at matched accuracy and hardware.
Pursue the fastest and the smallest-memory configurations, reporting their
tradeoffs: one configuration need not minimize both. Retain the stronger
same-comparator match/beat result separately. Missing comparator access or
unknown memory limits the claim; fallback never establishes native admission.

Predictable memory, user-supplied host/device limits and low-memory execution
are first-class requirements, especially for NEO-2 and cluster-node sizing.
A limit covers solver-owned allocations, not the entire process. Publish
process/runtime overhead and explicit headroom separately.

Execution is **single node, single process**, with OpenMP, NUMA-aware CPU
execution and optional CUDA. MPI, distributed execution, multiple GPUs and
out-of-core factors are outside this plan. Reach a usable GPU SPD checkpoint,
then prioritize CPU matrix-class coverage. Full GPU Pareto admission remains
a requirement later; it does not block CPU development.

Use `main` for coherent, locally verified integration and regular publication.
Use explicit, reentrant state and preserve existing calling conventions through
defaults. Update checked milestones only when their executable evidence exists.

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
kernels. Use checked 64-bit byte arithmetic. Add `host_memory_limit_bytes` and
`device_memory_limit_bytes` (zero means unlimited), with AUTO and LOW_MEMORY
policies, independently of UMFPACK's existing estimate-only memory gate.

A reusable analysis/planning API reports factor storage, persistent metadata,
peak workspace, supported RHS batch size and prediction quality (exact layout,
conservative bound or estimate). Count ordering and analysis scratch, padding,
task metadata, pinned buffers, pools, all copies and old/new overlap during
refactorization. Reserve before allocating, including growth and derivative work.
Opaque allocations cannot be advertised as strictly bounded: use an accountable
path or report the strict-budget operation unsupported.

LOW_MEMORY reduces concurrency, scratch, supernode padding and RHS batches;
it never silently changes solver, drops accuracy or spills to disk. If the
factors cannot fit, return budget-exceeded status with phase and required versus
available bytes. Enforce numerical-pivot growth limits even when pre-factor
fill is only an estimate. Preserve safe cleanup and recovery after exhaustion.
Report current/peak attributed storage separately from process RSS and runtime
cost; provide cluster-sizing guidance with explicit headroom and uncertainty.
Reuse keys cover every structural and scheduling assumption.

The first implemented hard-cap slice is natural-order CPU real-SPD simplicial
execution, with exact payload planning, capped graph scratch and retained
workspace for solves/JVP. See [memory contracts](docs/memory.md). Production
supernode/ordering and CUDA caps remain unimplemented; their strict requests
fail explicitly. This does not complete M6.1 or confer performance admission.

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

- [x] **M5.0--M5.3: real-SPD CPU MVP.** Native simplicial and production
  supernodal Cholesky, generated leaves, reuse, block RHS, implicit derivatives
  and observed-range CPU frontier evidence on the frozen 115-group corpus.
  Wider standalone kernel coverage and stronger match/beat outcomes remain open.
- [ ] **M6.0: consumer coverage.** Audit actual call paths, matrix classes,
  storage/index limits, transpose operations, RHS/reuse, accuracy, threading and
  memory needs across ITPplasma codes. Freeze consumer and held-out fixtures,
  including NEO-2; distinguish unavailable evidence from verified coverage.
- [ ] **M6.G: usable GPU SPD checkpoint.** Device triangular solves, bounded
  RHS batches, retained factors across solves/refactors, no mandatory factor
  downloads; independent primal/JVP/VJP and actual-device tests, complete memory
  attribution and matched GPU-library diagnostics. Full GPU frontier is deferred
  until CPU coverage, not removed. Current CUDA factorization is implemented;
  host triangular solves and incomplete peak-VRAM evidence remain limitations.
- [ ] **M6.1: bounded-memory native execution.** Planning, caps and RAM-only
  LOW_MEMORY for existing SPD paths first; require the same contract for every
  new family. Test analysis/numeric exhaustion, reuse overlap and recovery.
- [ ] **M6.2: real unsymmetric LU.** Sparse pivoted reference, then blocked/
  supernodal production execution; scaling, threshold pivoting, normal/transpose
  solves, refinement, diagnostics, bounded growth and consumer frontier gates.
- [ ] **M6.3: complex unsymmetric LU.** Add conjugate-transpose solves and
  explicit real-pairing derivatives; retain separate complex-class evidence.
- [ ] **M6.4: complex Hermitian positive-definite Cholesky.** Generated complex
  leaves, correct conjugation, bounded-memory CPU schedules and admission.
- [ ] **M6.5: indefinite native factors.** Separately admit pivoted real
  symmetric, complex Hermitian and complex transpose-symmetric factorizations;
  stable 1x1/2x2 pivots, delayed-pivot dependencies, growth caps and adversarial
  stability tests. Do not substitute unpivoted LDL for robust indefinite solves.
- [ ] **M6.6: consumer integration and CPU independence.** Consumer-required
  sparse preconditioners (Krylov remains in FortNum), native configuration tests
  without external sparse solvers, representative/held-out admission evidence,
  cluster-memory predictions and migration documentation. Expand index width
  when the consumer inventory demonstrates a requirement.
- [ ] **M6.7: GPU admission and coverage.** Profile unresolved frozen SPD cases
  to the GPU-library Pareto frontier, then add consumer-required GPU classes
  with the same correctness, derivative, memory and performance gates.

For each family implement a small independent reference, verify the public
lifecycle and derivative contracts, then profile production execution until the
declared gates pass. Maintain the existing CPU corpus throughout. Native default
selection is promoted only for admitted classes/configurations; correctness or
an optimized microkernel alone is not performance admission.

Normal, transpose and conjugate-transpose solves reuse the same factors; public
APIs retain normal solves by default. Diagnostics distinguish invalid matrices,
nonpositive pivots, singularity, allocation failure, budget exhaustion and
unsupported claims. Never relabel an estimate as a rigorous error certificate.

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
