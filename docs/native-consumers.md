# Native-only consumer migrations

Status on 2026-10-16. Completeness and unchanged scientific correctness gates
precede performance tuning. These are draft migrations, not release admission.
No consumer PR is merged as part of this work. Dense BLAS/LAPACK and optional
AMD/METIS graph ordering are allowed; external sparse factorizations are not.

| Consumer | State | Evidence or remaining gate |
|---|---|---|
| FortFEM | [Draft PR 82](https://github.com/lazy-fortran/fortfem/pull/82) | Native factor/solve, derivative, C API and actual fpm checks pass. Preserve KIN6D's newer scientific dependency sources when selecting its pin. |
| NEO-2 | [Draft PR 212](https://github.com/itpplasma/NEO-2/pull/212) | PAR/QL builds, adapter/lifecycle and Arnoldi checks pass; external sparse linkage absent. Unchanged spline node-reconstruction regression fails. Lorentz output differences remain under investigation; the retained old run also fails the strict golden comparator, so not every discrepancy is attributable to migration. Not ready to merge. |
| MEPHIT | In progress | Native Fortran adapters and default binary checks pass; optional MFEM sites and the FreeFem Maxwell daemon both require native migration and independent application checks. Binary linkage alone is insufficient. |
| KAMEL | In progress | Native real/complex compatibility adapter tests pass; actual KIM/QL consumer builds and checks pending. |
| KIN6D | In progress | Explicit native adapters and dependency configuration; scientific dependency pin and affected physics/derivative checks pending. |
| NEO-EQ | In progress | Actual pressure/current complex sparse paths identified; native adapter and application build checks pending. |
| NEO-RT | Default path independent | Optional non-standalone path inherits NEO-2; update its dependency after a qualified NEO-2 revision exists. |

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

The FortSparse public C ABI and native-default fpm configuration support these
migrations. C ABI input conversion and solve staging do not yet expose a strict
user memory cap. Explicit-AMD real SPD CPU factors support payload caps and
refinement; ordered planning and strict ND/GPU caps remain unsupported. These
limits and the later CPU/GPU Pareto requirement remain active roadmap work.
