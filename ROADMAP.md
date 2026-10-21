# FortSparse roadmap

Active plan and controller handoff, 2026-10-16. This is the single authoritative
plan; do not create a competing PLAN.md. Checked items describe the implemented
baseline; unchecked items are planned capabilities, not implementation claims.

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

Current delivery priority is **feature completeness and correctness before
further performance tuning**. Maintain the seven existing native-only consumer
draft PRs; do not merge them during this handoff. SPEC, JOREK and M3D-C1 are
deferred forks and have no current priority. Track actual paths and merge blockers
in [native consumer migrations](docs/native-consumers.md). Pin verified library revisions, preserve each
consumer's mathematical/storage/lifetime contracts, and run affected independent
checks. Each draft PR reports matched baseline/native CPU-time and memory changes
for the exercised tests and runs, including failures and unavailable comparisons.
Document unsupported paths explicitly; a draft PR is not validation or
release admission. The Pareto requirement above remains a later gate.

Use `main` for coherent, locally verified integration and regular publication.
Use explicit, reentrant state and preserve existing calling conventions through
defaults. Update checked milestones only when their executable evidence exists.

## Current integration checkpoint

### October 21 continuation checkpoint

The user requested that agents conclude and preserve current evidence. No
FortSparse source job remains active in this lane. Recheck forge state before
promotion; required checks execute on owned hardware in the normal merge queue,
with PR CI skipped. Preserve the original one-second native test deadlines.

[PR #26](https://github.com/lazy-fortran/fortsparse/pull/26) merged as
6d17059c7741b280d04998f479a160c341506629. It batches uncapped real-LDLT scalar
recovery contractions by disjoint rows while retaining canonical contraction
order. Strict-cap and LOW workspace paths remain excluded. Exact candidate
38ea54ca95f481a0b9cb7daf2cd024ed8962b3bf passed six native tests at each of
1/2/4/8/16 threads (30/30, cluster job 3321634) with unchanged deadlines.

A matched [GLISS](https://github.com/itpplasma/GLISS) experiment, job3321625,
retains 150 fresh processes: five trials at two public P4 configurations,
1/2/4/8/16 cores, block control, sparse baseline and repaired sparse candidate.
At one thread the sparse repair improves solve medians by 17–21% and total wall
by 6–7%, with similar CPU cost. Sparse process memory is about 125 MiB versus
about 65 MiB for the faster block control. More threads slow all solve paths.
This is a bounded consumer improvement, not strong-scaling or Pareto admission.

[PR #27](https://github.com/lazy-fortran/fortsparse/pull/27) merged as
cc706d011f925d539308da9b81368abd4b42dfdb after strict owned CI #24 merged as
0a625bad57e3829e47ea040dd3a77a1e57440594. It adds an explicit LP64 DGESV
interface to the existing mixed-LU test; production solver bytes are unchanged.
Official fpm 0.13.0/GNU 12.2 reproduced the missing-interface failure at 6d17059
(job 21817508), then passed the affected mixed-LU/precision controls at exact
8db2fe7eecd93999c191cdbc7b0e38c328d81247 (job 21817529). The separate Fo missing
OpenMP link propagation remains with the tools owner; do not add a workaround
or claim that adapter defect repaired by this source change.

The large private [KIN6D](https://github.com/itpplasma/kin6d) rest-block factor/count remains unqualified. Three
six-hour pilots timed out; bounded extracted-fixture replay also remained
incomplete. Profiles locate scalar recovery after blocked-path fallback.
Panel-cap experiments and an explicit-zero hypothesis did not establish a
qualified repair and were not retained. No private matrix is published, no
large-C inertia/solve claim is made, and no six-hour campaign should restart
until the actual fixture admits a completed factor/count/solve under its
unchanged numerical gates. Root/controller review owns any future promotion.

Local canonical receipts, reports, scripts and two inspected comparison/scaling
plots are indexed by ~/data/handoffs/20261021_completion/sparse_gliss.md and
sparse_gliss/report.md in that directory. Raw application inputs remain outside
Git with checksums. The full sparse_gliss manifest was verified for 1618 files.
[GLISS PR #59](https://github.com/itpplasma/GLISS/pull/59) retains the qualified 6d17059 production pin; do not repin it to the
later test-only #27 change. Its normal promotion remains held until its strict
CI #58 actually merges and owned required gates pass.

The October 16 checkpoint below is historical and must not override this
revision-bound continuation record.

FortSparse had **no open PRs** at the 2026-10-16 handoff check. Published `main`
was `5ad34ddaba1ac51de0bff5f0bb4dc5c060bd33ff` before this planning update.
There is therefore no outstanding FortSparse PR to merge. Consumer PRs remain
drafts, with their existing pins, tests and unresolved gates recorded in
[native consumer migrations](docs/native-consumers.md). Recheck remote state
when taking over; a local worker branch or proposal is not an open PR or an
integrated implementation. No workers are running at this checkpoint.

The real-SPD CPU MVP has the bounded corpus admission described below. Other
families have useful native references but incomplete production admission.
Do not use an overall completion percentage as an acceptance criterion.

## Priority work order for the next controller

These priorities order delivery; the M6 identifiers below retain their existing
meaning. Complete one bounded change and its affected checks before expanding
scope. A documented application-science blocker need not stop independent
FortSparse correctness or coverage work.

1. **Close CPU correctness and ownership gaps.** Inventory each real consumer's
   missing operation against [native CPU support](docs/native-cpu.md). Repair
   reproducible original-matrix, transpose/adjoint, pivot, reuse, derivative and
   exhaustion defects before optimizing. Public solver contexts now have an
   explicit copy contract (`sparse_copy` and defined assignment): native CPU
   factors are deep-copied under the destination host cap; external and CUDA
   factors are not copied and leave an unfactored destination. gfortran 16
   does not invoke the defined assignment for a `sparse_solver_t` component of
   an enclosing derived type, so such components must be copied explicitly.
   Retain independent known-solution oracles and explicit unsupported paths.
2. **Bound consumer qualification work.** Apply the NEO-2 stop rule below before
   further golden investigation. Diagnose FortFEM's two failing CI checks and
   KAMEL's golden/Linux/macOS failures, distinguishing baseline failures from
   migration regressions. NEO-EQ's repaired frozen pair agrees to about 3.3e-7
   in the reported field diagnostic, but both outer iterations diverge. Treat
   physical convergence as a separate application issue: no unbounded mesh,
   basis, damping or case search, and no invented peer tolerance. MEPHIT's
   cylinder is qualified; full optional MFEM physics remains unverified.
   NEO-RT's declared background approximation remains explicit. Keep all seven
   consumer PRs unmerged until separately authorized.
3. **Finish predictable memory for required CPU paths.** Extend analysis and
   pre-factor plans beyond compact natural/AMD SPD, and enforce owned-payload
   limits through conversion, ordering, replacement overlap, factor growth,
   refinement and derivatives. Add RAM-only LOW modes where supported; reject
   unsupported strict requests. In particular, ND currently has no strict host
   cap. Record exact sizes, conservative bounds and estimates distinctly, with
   separate process overhead/headroom for NEO-2 and node sizing. The real LOW
   page-advice change is already on main: measured RSS falls 12--14% at about
   19% more CPU on its two frozen workloads. Preserve that tradeoff; do not
   generalize it to AUTO, complex factors or a stronger payload bound.
4. **Complete required CPU families, then production schedules.** Fill only
   demonstrated consumer gaps in complex LU, HPD and pivoted indefinite paths,
   their ordering, retained operations and planning. References already exist;
   avoid rebuilding them. After the affected correctness contracts pass, add
   generated blocked execution and OpenMP scaling. General real LU is the most
   urgent practical performance gap: the controlled full KIN6D case times out
   at 900 seconds against a roughly 6-second baseline and reaches about 7 GiB
   RSS. Retain this failure and the original physical inputs. Implement the
   static-structure supernodal LU below before repeating the full consumer run.
5. **Repair blocking tooling in its owner; qualify pending generated work.**
   Fo's newer CMake File API route omits the scanner fallback used by its older
   route. The tiny BVLS reproducer retains a readable source on two failed
   resident generations; diagnostic regex fallback succeeds. Implement and
   independently test the bounded upstream integration repair, preserving
   preprocessing/active-branch semantics and parser diagnostics; then rerun
   the original focused FortNum gate. Do not alter scientific source, replace
   shared installed tools or label a diagnostic workaround a passing gate.
   The unpublished FortNum bicubic candidate `97d5a3b` has native/fpm,
   symbolic and full regeneration evidence, but its resident Fo gate remains
   blocked. After that gate passes, trial its generated Horner value/partial
   jet in an isolated NEO-2 adapter for NEO-RT, preserving the existing primal
   interpolation and physical checks. Measure the complete consumer workload
   before upstream promotion and public pin updates. Do not assume fewer
   operations make the analytical JVP/VJP schedules fastest.
6. **Resume CPU Pareto admission, then extend GPU admission.** After functional
   correctness, compare each declared CPU workload with the best available
   external frontier at matched accuracy, with CPU time, RSS, payload and
   1/2/4/8/16-core scaling. Then complete usable GPU SPD memory evidence and
   matched actual-GPU-library admission, followed by consumer-required GPU
   families. GPU host caps, LOW and complete peak-VRAM admission remain open.
   CPU-versus-GPU comparisons do not establish GPU-library competitiveness.
   Keep the full CPU/GPU Pareto requirement; no MPI or multi-GPU expansion.
7. **Publish final evidence and close the handoff.** Refresh affected consumer
   pins only after upstream publication and recheck the exact public closure.
   Put current CPU/RSS numbers or explicit unavailable results in each draft;
   preserve older tables at their measured source identities. Push verified
   FortSparse increments regularly. Publish PNG Pareto/scaling plots through
   ephemeral slopbox, retaining recipes, numeric evidence and provenance in Git.
   Recoverably retire completed worktrees/branches only after integration,
   stopped owned jobs and an audit of all remaining generator dependencies.
   Preserve negative experiments and useful unpublished work.

Tooling repair in item 5 may run earlier when it blocks a required check in
items 1--4. Do not start unrelated upstream repairs or performance campaigns.

## NEO-2 accuracy investigation and stop rule

The Lorentz golden is highly sensitive. Cancellation in downstream derivative
projections can amplify tiny solution differences; historical floating-point
tables are comparison data, not an exact mathematical oracle. Generated
compensated solves passed the strict final transport comparison at an earlier
solver-only checkpoint. This is not a full golden pass for the combined branch.

First reuse the retained matrix/RHS/factor captures, independent references and
failed-gate receipts. For any genuinely new discrepancy, isolate assembly,
factorization, triangular solve, residual evaluation and observable projection.
Check the original equations, an independent reference and the reported physical
observable. Include conditioning, representable-input rounding and projection
cancellation in the explanation. Higher-precision evidence is not a rigorous
certificate, and a small backward error alone does not bound forward error.

**Stop solver-side tuning** when independent checks support the solver's stated
accuracy contract, remaining discrepancies are explained within a justified
conditioning/rounding analysis, and there is no reproducible FortSparse defect.
Write the evidence and remaining limitation once, then move to the next CPU
priority. Do not keep changing kernels to reproduce cancellation-amplified
intermediate tables. If that explanation cannot be established, record the
specific unresolved question and the next discriminating check; do not launch
an open-ended sequence of full golden runs.

Preserve all existing gates and their failures. Classify provenance/metadata
differences separately against their actual schema contract. If scientific
acceptance needs an observable accuracy budget that has not been supplied,
leave the application gate open for its scientific owner. Do not invent a
tolerance, silently refresh goldens, discard intermediates, or promote the
consumer PR. Repeat an expensive run only after a concrete fix or new check
can resolve the recorded question. A new independently reproduced solver defect
reopens solver work; unexplained historical-bit disagreement alone does not.

## General LU: static-structure supernodal design

Decided 2026-10-16: close the large-front LU gap with a PARDISO-style
static-structure supernodal LU, specified in [the LU design](docs/lu-supernodal.md).
Maximum-product matching with dual scaling, ND/AMD on the matched A+Aᵀ, a static
symmetric symbolic structure with exact payload planning, left-looking BLAS-3
updates reusing the SPD schedule, pivoting restricted to supernode diagonal
blocks with counted perturbation, and an independent original-matrix gate that
falls back, reported, to the existing threshold-pivoting scalar LU.

The [static-structure diagnostic](bench/results/2026-10-16-lu-static-structure/README.md)
on both frozen KIN6D Jacobians gives static L+U of 0.97 M and 4.94 M entries
against 2.0--57.5 M for pivoted scalar LU and 2.2/15.3 M for SuperLU, with zero
perturbations and componentwise backward error near 5e-16 after one refinement
step in dense emulation. Structural (unweighted) matching fails this scheme. These
are structure/accuracy findings on one consumer, not performance admission.

The width-eight prefix source-block experiment is superseded and will not be
implemented; its 60 negative observations remain in the
[panel evidence](bench/results/2026-10-16-lu-panels/README.md). The public
default stays scalar until the staged implementation passes its correctness
tests and the frozen-corpus frontier comparison, which adds ParU as a peer.

## Local handoff evidence

These paths retain supporting evidence and candidate work on the current host;
they are not build dependencies. Published contracts and benchmark summaries
remain in this repository. Inspect worktree status before resuming any candidate.

| Work | Retained location |
|---|---|
| Exact source-block proposal and structural diagnostic | `~/code/sparse-dev/benchmarks/2026-10-16-kin6d-native-migration/nonnested-csc-probe/lu-source-structure/` |
| FortNum bicubic candidate and verification receipt | `~/code/worktrees/fortnum-bicubic-horner` at `97d5a3b`; `~/code/sparse-dev/migrations/neort-native/restoration-2026-10-16/bicubic-horner/receipt.json` |
| Fo two-generation BVLS reproducer and repair proposal | `/var/tmp/fo-captured-bvls-audit/REPORT.md`, `project/` and `receipts/` alongside it |
| Recoverable consumer-worker cleanup and dependency correction | `~/code/sparse-dev/cleanup/2026-10-16-completed-consumer-audit/` |

Protect `/var/tmp/neo2-curl-gen-build` and the active generator source closure:
the shared FortSym build is still required despite its historical consumer
name. Its temporary quarantine was reversed. A clean tree or an absence of
running processes alone does not establish that a build directory is unused;
audit manual generation recipes and obtain the remaining owners' release.

## Implemented baseline

- [x] Real/complex CSC storage and construction, sparse products, status API,
  backend dispatch and solver lifetime rules.
- [x] In-process SuperLU and a separate GPL UMFPACK helper.
- [x] Optional MUMPS backend, retained factors, block solves, refinement,
  pattern reuse and solver diagnostics.
- [x] Retained normal/transpose/adjoint solves across SuperLU, UMFPACK IPC
  and MUMPS, with operation-aware refinement and helper capability negotiation.
- [x] Real/complex equation-level solve JVP/VJP products with independent
  numerical and adjoint checks.
- [x] Backend selection and BLAS discovery, including the fetched-SuperLU
  discovery fix verified through KIN6D.

- [x] Native real SPD simplicial reference and blocked factorization, optional
  AMD/METIS ordering, retained analysis, block RHS and implicit derivatives.
- [x] Native CPU sparse real/complex pivoted LU references, complex HPD
  simplicial Cholesky and real/complex symmetric or Hermitian pivoted LDL
  references, with original-matrix refinement, retained
  transpose/adjoint solves and allocation caps. Production admission remains open.
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

Implemented host caps cover natural-order CPU SPD compact and blocked execution,
and the CPU factor references. General LU and real SPD also support explicit
AMD with audited caller-owned ordering workspace, including LOW_MEMORY. General
LU accounts for reuse overlap; ordered SPD accounts for analysis, panels and
original-matrix refinement.
Bounded panels count descriptors, graph/task
scratch and requested-team workspace, reduce teams to fit and avoid opaque BLAS
workspace. Compact natural SPD and explicit-AMD LOW_MEMORY SPD have exact
payload planning, including ordering and symbolic-count phase peaks. Blocked
planning, SPD ND host caps and wider production planning remain unimplemented.
CUDA real SPD now enforces a separate explicit device-allocation payload cap,
including plan-replacement overlap and failure diagnostics. CUDA host caps and
LOW_MEMORY remain unsupported; opaque runtime/driver usage stays outside the
owned payload scope. Strict unsupported requests fail explicitly. See
[memory contracts](docs/memory.md). This does not complete M6.1
or confer performance admission.

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
  device triangular solves and RHS tiling are implemented and verified. Complete
  peak-VRAM benchmark admission remains open.
- [ ] **M6.1: bounded-memory native execution.** Planning, caps and RAM-only
  LOW_MEMORY for existing SPD paths first; require the same contract for every
  new family. Test analysis/numeric exhaustion, reuse overlap and recovery.
- [ ] **M6.2: real unsymmetric LU.** Sparse pivoted reference, then blocked/
  supernodal production execution; scaling, threshold pivoting, normal/transpose
  solves, refinement, diagnostics, bounded growth and consumer frontier gates.
  The KIN6D coupled gauge/voltage Jacobian exposes a practical ordering/fill
  gap. Opt-in structural matching with matched-graph AMD and generated pivot-row
  preferences is implemented with original-matrix and cap-accounting checks.
  Optional METIS ND now supports general real/complex LU without strict caps.
  The matched-AMD full consumer case times out at 900 seconds on its largest
  mesh. A frozen 19,989-row replay completes with unmatched ND in 5.49 seconds
  of factor time and 608 MB retained payload; this shared-host diagnostic is
  not performance admission. An unchanged full consumer run passes with
  unmatched ND in 797.81 shared-host elapsed seconds, but a subsequent controlled
  CPU/RSS attempt times out at 900 seconds while its original baseline passes.
  Large general LU remains a performance gap; retain the
  [failed attempt and reproduction data](bench/results/2026-10-16-kin6d-nonnested/README.md).
  A generated real width-eight rolling target-panel prototype is available
  only through low-level opt-in. Independent closure/pivot/derivative and exact
  payload/fault-recovery tests pass; all three bounded A/B cohorts pass their
  residual gates but do not establish a reproducible workflow gain over scalar
  execution. Full-active and contiguous-range generated SIMD are implemented;
  the static-structure supernodal design above supersedes further panel work.
  No public default or family admission follows from this prototype.
  Generated compensated triangular solves improve the NEO-2 final
  transport comparison while retaining allocation-free solves and capped
  correction storage; intermediate application checks remain open.
- [ ] **M6.3: complex unsymmetric LU.** The serial native reference already
  supports normal/transpose/adjoint retained solves, real-pairing implicit
  derivatives, refinement, AMD/ND and bounded natural/AMD execution. Generated
  compensated triangular accumulation is implemented. Production blocking,
  parallel factorization and separate complex consumer frontier admission remain
  open; the real panel prototype does not extend to complex factors.
- [ ] **M6.4: complex Hermitian positive-definite Cholesky.** Generated complex
  leaves, conjugation checks, retained operations/refinement and bounded-memory
  natural simplicial execution are implemented. Structured complex ordering,
  production blocking/parallel schedules, pre-factor plans and admission remain
  open.
- [ ] **M6.5: indefinite native factors.** Sparse Bunch--Kaufman references for
  real symmetric, complex Hermitian and complex transpose-symmetric matrices
  already implement unrestricted 1x1/2x2 pivots, retained operations/refinement
  and growth caps. Real symmetric AMD is implemented; complex structured
  ordering, parallel blocked schedules, pre-factor plans and representative
  adversarial/consumer performance admission remain open. Do not substitute
  unpivoted LDL for robust indefinite solves. Exact current class contracts are
  retained in [native CPU support](docs/native-cpu.md).
- [ ] **M6.6: consumer integration and CPU independence.** Consumer-required
  sparse preconditioners (Krylov remains in FortNum), native configuration tests
  without external sparse solvers, representative/held-out admission evidence,
  cluster-memory predictions and migration documentation. The whole-context
  copy contract is implemented and tested for every native CPU family and the
  SuperLU/UMFPACK-IPC handles; MUMPS and CUDA copies are compile-checked
  only on hosts without those builds. Expand index width
  when the consumer inventory demonstrates a requirement.
  Native real CSR [ILUTP](docs/ilutp.md) with explicit fill caps, column pivots,
  bounded payload and allocation-free retained normal apply is implemented.
  First-party draft migrations and measured CPU/RSS changes are tracked in
  [native consumers](docs/native-consumers.md). Unqualified scientific or runtime
  paths remain explicit. The deferred SPEC fork checkpoint does not take current
  priority; representative preconditioner performance admission remains pending.
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
Before frontier admission or consumer PRs use the supernodal LU, scan its new
sources and generated leaves for code similarity against the retained solver
snapshots (SuperLU, SuiteSparse, MUMPS, PaStiX, STRUMPACK, SPRAL), review every
hit and record the result with the algorithm citations in the LU design.

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
