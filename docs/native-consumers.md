# Native-only consumer migrations

Status on 2026-10-16. Completeness and unchanged scientific correctness gates
precede performance tuning. These are draft migrations, not release admission.
No consumer PR is merged as part of this work. Dense BLAS/LAPACK and optional
AMD/METIS graph ordering are allowed; external sparse factorizations are not.

| Consumer | State | Evidence or remaining gate |
|---|---|---|
| FortFEM | [Draft PR 82](https://github.com/lazy-fortran/fortfem/pull/82) | Native CMake/fpm checks pass. Preserves KIN6D's newer scientific modules; exact Hermitian Gram assembly fixed and independently checked upstream. Nine-repeat adapter CPU/RSS table included; that table predates the Gram fix. |
| NEO-2 | [Draft PR 212](https://github.com/itpplasma/NEO-2/pull/212) | Combined e9e2877d builds PAR/QL and passes 21 affected checks, including int64 payload limits/LOW policy, terminal rank failure, complete Arnoldi projection and generated full covariant curl from consistent interpolants. Generated compensated LU solves passed the original strict final Lorentz transport comparison at an earlier solver-only checkpoint; intermediate and metadata golden differences remain under investigation. No full application pass is claimed for the combined branch. Current three-pair provider-isolation spline measurements show CPU -16.18% and peak process RSS -28.80%; the provider-independent Arnoldi control is labeled separately. Historical Lorentz CPU +7.07%, largest-process RSS +3.71% predates compensation. Not ready to merge. |
| MEPHIT | [Draft PR 71](https://github.com/itpplasma/MEPHIT/pull/71) | Native FDM, both optional MFEM sites and shared FreeFem module. Stock 10/10 and minimal-runtime/MFEM 11/11 focused checks pass; complete cylinder passes unchanged scientific gates. Corrected stock-native CPU +20.89%; component RSS and optional-runtime results included. Failed initial IPC run and its bounded protocol repair retained. |
| KAMEL | [Draft PR 357](https://github.com/itpplasma/KAMEL/pull/357) | Actual KIM/QL builds and 11 affected checks pass. Matched five-case CPU/RSS table included; long periodic/global CPU -0.7%, RSS +2.4%. Coarse diagnostics remain explicitly qualified; full golden/cross-platform admission pending. |
| KIN6D | [Draft PR 9](https://github.com/itpplasma/kin6d/pull/9) | Native-only adapters, updated FortFEM and explicit assembly-roundoff normalization. Six matched historical cases at dfa0b79: 44/44 measured processes pass; CPU/RSS tables included. Vertical Newton CPU +14.9%, RSS -6.2%; direct/Newton CPU +1.0%, RSS -12.2%. Current fc8b36a selects optional ND without matching for the coupled Jacobian. The original 8/16/32-level nonnested test now passes unchanged in 797.81 shared-host elapsed seconds; eleven other affected checks pass. Previous matching+AMD timeout remains recorded. A subsequent controlled pair has a passing baseline and a native 900-second timeout; [observed failed-attempt costs](../bench/results/2026-10-16-kin6d-nonnested/README.md) are retained, with no successful full-case delta. |
| NEO-EQ | [Draft PR 9](https://github.com/itpplasma/NEO-EQ/pull/9) | Native pressure/current and shared FreeFem Ampere paths build; focused checks and component CPU/RSS tables included. Both legacy and native AMD complete the frozen public KiLCA solver phase after shared mesh and complete-Arnoldi repairs, but both outer iterations remain nonconverged and the final peer fields disagree. Independent captured-system checks agree with SuperLU; these do not establish coupled scientific equivalence. A legacy-only larger-basis diagnostic also fails to establish convergence. The earlier natural-order native attempt and unchanged legacy postprocessing time out. The scaled affine-action repair passes native 11/11 and legacy 8/8 affected checks; the unchanged physical pair completes with field differences about 3.3e-7, but both steady iterations still diverge and no peer-accuracy pass is assigned. No tolerance relaxation or successful full-pipeline comparison claimed. FreeFem runtime migration is not an isolated solver comparison. |
| NEO-RT | [Draft PR 198](https://github.com/itpplasma/NEO-RT/pull/198) | Optional native application builds after explicit NEO-2 field-provider lifetime and signed-curl repairs. Default 5/5 and optional 6/6 checks pass, including unchanged nonzero ripple-plateau transport gates. The complete-curl dependency refresh at 08c7c61 passes both configurations and unchanged transport gates. Twelve fresh optional measurements pass; ripple CPU is 12.237437 s and process RSS 47728 KiB. Earlier default and optional tables remain explicitly historical; the current optional cost is higher and its original baseline was unbuildable. The original optional application was unbuildable, so its baseline delta remains unavailable; new optional absolute costs are reported separately. |

GLISS (itpplasma/GLISS) hand-writes a block-tridiagonal symmetric-indefinite
LDLᵀ for shifted pencils K − σM with Sylvester inertia and a reliability flag.
Native LDLᵀ now reports exact inertia with reliability status, and its counts
agree with GLISS on real GLISS pencils. No migration is planned: on
those pencils the serial scalar native LDLᵀ is 14–40× slower per shift
(unqualified Mac runs) than GLISS's dense-block kernel. Reconsider only after a
supernodal pivoted LDLᵀ is within about 2× per shift. Static pivot
perturbation (supernodal LU) is unsuitable for inertia.

CHEASE (itpplasma/chease) factors its Grad–Shafranov bicubic-Hermite system
with a hand-written banded LDLᵀ (`ALDLT`, `src-f90/matrix.f90`), refactored in
the ballooning loop. It is the strongest remaining candidate: a sparse
ordering could change its cost class. Its owners are actively tuning that
solver under EPFL parity constraints, so no FortSparse migration is planned.
ARES and tiago hand-write block-tridiagonal (Thomas-type) solvers in the right
cost class; they need no migration. A 2026-10-16 audit of all itpplasma and
lazy-fortran repositories found no other meaningful candidate.

The audited default paths of libneo, SIMPLE, GORILLA/GORILLA_APPLETS, ARES,
NEO, tiago, KiLCA/KiLCA-FLR2 and REALSPACE_WITH_COLLISIONS have no external
sparse factorization to replace. Avoid adding an unused dependency or a nominal
migration PR. Re-audit if their actual build graph changes. Workflow repositories
such as rmp_torque are not solver owners.

SPEC, JOREK and M3D-C1 are deferred forks. Retain their isolated checkpoints,
but allocate no current migration priority to them. The old magdif Python/paper
prototype is likewise deferred; this does not defer the separate NEO-EQ code.

Each genuine migration must retain storage conventions, matrix class, index
normalization, factor/reuse/free semantics, original-matrix error checks and
consumer tolerances. An external fallback or a smaller residual alone does not
establish scientific equivalence. Keep failures visible in draft PRs, repair
generic defects upstream, then pin and recheck the original consumer.

Every draft includes baseline/native CPU seconds, individual-process peak RSS,
absolute/relative deltas, exact tested commands, repeats, thread settings and
correctness status. Short process timings include startup. Failed or unbuilt
paths have no successful comparison. Multi-process component peaks are not an
aggregate concurrent peak, and waited-descendant CPU is counted only once.
These shared-host migration observations do not establish a Pareto-front claim.

The FortSparse public C ABI and native-default fpm configuration support these
migrations. C ABI input conversion and solve staging now expose strict owned
host-payload caps. Structured RAW input conversion supports the CPU matrix
classes, including old/new factor overlap. Explicit-AMD real SPD supports exact
compact LOW planning and strict execution; real LDLT and general LU support
bounded AMD execution. Optional general-LU ND is available without a strict
budget. CUDA real SPD enforces explicit device-payload caps; GPU host caps,
GPU LOW, ND host caps and broader production planning remain open. A solver
payload cap does not bound consumer assembly, opaque runtime memory or
whole-process RSS. These limits
and the later CPU/GPU Pareto requirement remain active roadmap work.
