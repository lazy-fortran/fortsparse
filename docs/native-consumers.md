# Native-only consumer migrations

Status on 2026-10-16. Completeness and unchanged scientific correctness gates
precede performance tuning. These are draft migrations, not release admission.
No consumer PR is merged as part of this work. Dense BLAS/LAPACK and optional
AMD/METIS graph ordering are allowed; external sparse factorizations are not.

| Consumer | State | Evidence or remaining gate |
|---|---|---|
| FortFEM | [Draft PR 82](https://github.com/lazy-fortran/fortfem/pull/82) | Native CMake/fpm checks pass. Preserves KIN6D's newer scientific modules; exact Hermitian Gram assembly fixed and independently checked upstream. Nine-repeat adapter CPU/RSS table included; that table predates the Gram fix. |
| NEO-2 | [Draft PR 212](https://github.com/itpplasma/NEO-2/pull/212) | PAR/QL builds and all four focused checks pass, including the unchanged spline gate after upstream row-singleton preprocessing. Lorentz stages complete; strict application golden still fails. Fresh legacy repeats have identical scientific datasets but differing metadata. Draft reports native CPU +7.07%, largest-process RSS +3.71% for the matched Lorentz run. Not ready to merge. |
| MEPHIT | [Draft PR 71](https://github.com/itpplasma/MEPHIT/pull/71) | Native FDM, both optional MFEM sites and shared FreeFem module. Stock 10/10 and minimal-runtime/MFEM 11/11 focused checks pass; complete cylinder passes unchanged scientific gates. Corrected stock-native CPU +20.89%; component RSS and optional-runtime results included. Failed initial IPC run and its bounded protocol repair retained. |
| KAMEL | [Draft PR 357](https://github.com/itpplasma/KAMEL/pull/357) | Actual KIM/QL builds and 11 affected checks pass. Matched five-case CPU/RSS table included; long periodic/global CPU -0.7%, RSS +2.4%. Coarse diagnostics remain explicitly qualified; full golden/cross-platform admission pending. |
| KIN6D | [Draft PR 9](https://github.com/itpplasma/kin6d/pull/9) | Native-only adapters, updated FortFEM and explicit assembly-roundoff normalization. Six matched historical cases at dfa0b79: 44/44 measured processes pass; CPU/RSS tables included. Vertical Newton CPU +14.9%, RSS -6.2%; direct/Newton CPU +1.0%, RSS -12.2%. Current d7dbf24 opts into native matching for the coupled Jacobian; eleven other affected checks pass. Original 8/16/32-level nonnested test times out at 900.64 seconds on the largest level after successful first two levels. No successful full-case CPU/RSS delta claimed. |
| NEO-EQ | [Draft PR 9](https://github.com/itpplasma/NEO-EQ/pull/9) | Native pressure/current and shared FreeFem Ampere paths build; three affected checks and fresh component CI pass. Three-repeat CPU/RSS tables included. Full physical run blocked by stale minimal input and missing current device field data; no tolerance/control relaxation. FreeFem runtime migration is not an isolated solver comparison. |
| NEO-RT | [Draft PR 198](https://github.com/itpplasma/NEO-RT/pull/198) | Optional coupling pins native NEO-2; default standalone remains independent and three focused checks pass with CPU/RSS tables. Optional application compile is blocked by an existing incompatible field-adapter interface; optional-native runtime deltas unavailable. |

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
migrations. C ABI input conversion and solve staging do not yet expose a strict
user memory cap. Explicit-AMD real SPD and pivoted LDLT CPU factors support
payload caps; ordered planning and strict ND/GPU caps remain unsupported. These
limits and the later CPU/GPU Pareto requirement remain active roadmap work.
