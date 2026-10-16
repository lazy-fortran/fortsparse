# Single-node consumer follow-up

Audit date: 2026-10-16. This extends [consumer coverage](consumer-coverage.md)
and [its fixture manifest](consumer-fixtures.json). Source inspection and exact
five-by-five arithmetic only: no application build, production run or solver
benchmark. The source files read here match the pinned repository commits;
unrelated `.gitignore`/`.sloptools` changes were preserved. References below use
NEO-2 `2b8d5c0e`, MEPHIT `5bab300c`, KIN6D `73891c54`; full revisions are in the
manifest. FortFEM's adapter is a separate dependency, not duplicated KIN6D code.

## What single-node admission means

The native plan covers local CPU/GPU factorization and solves, with no native
MPI wrapper, distributed factors, multi-GPU or disk-backed factors. An existing
MPI application can still call a local native factorization for a complete
matrix held by each rank. That is distinct from gathering/distributing a global
matrix or replacing a distributed factorization collectively.

| Existing execution path | Single-node/local interpretation | Excluded or unverified requirement |
|---|---|---|
| NEO-2 PAR and QL complete CSC per solver call | Local algebra can use a native backend without MPI in FortSparse. Application orchestration stays outside it. | No MPI application integration executed. QL multispecies uses rank/species collectives; one rank is not certified equivalent to a multispecies run. |
| MEPHIT FDM and optional serial MFEM matrix | Complete matrix and RHS are passed locally. No distributed sparse contract appears at these call sites. | Full application/threading/runtime not executed. |
| KIN6D GS, ideal, Ampere and coupled operators | Local CSC/CSR/real embedding and local band solves; native local replacements fit the algebra. | Existing MPI-MUMPS runtime initialization is backend-specific, not a requirement to add native MPI. |
| JOREK MUMPS adapter | The inspected adapter accepts an application communicator; a one-rank communicator is not rejected by its interface. This is a possible one-rank baseline, not an executed configuration. | Distributed COO ownership, collective factorization and MPI broadcast semantics remain outside the native plan even if all ranks run on one physical node. |
| M3D-C1 PETSc MPIAIJ | One-rank/local algebra may provide an extraction/reference path; a working one-rank application configuration was not established. | PETSc distributed ownership, KSP/preconditioner composition and MPI execution are not replaced by the local sparse API. |
| GPEC, MARS-Q, SPEC, STELLOPT structured paths | Previously inspected band/dense/local preconditioner calls are single-process algebra, even inside an MPI application. | Native general sparse admission does not automatically establish band/block solver or application parity. |

NEO-2 `COMMON/collision_operator_mems.f90:52` checks process count against
species count for multispecies work. A mismatch with more than one process
stops; one process is allowed only with an explicit testing warning. QL axi
`ripple_solver_axi_test.f90:269` sets species from rank and uses allgather at
`:2510`. Therefore a one-rank multispecies run must not silently become the
reference for the intended species-coupled workload.

KIN6D `src/fem2d/gs_direct.f90:80` initializes/queries MPI only under
`KIN6D_MUMPS_MPI`, and requires serialized thread support there. Sequential
MUMPS and a future local native solver do not need that runtime boundary.
`gs_direct_available()` currently returns true (`:32`); fallback source
branches are not evidence that the current configuration exercises them.

## NEO-2 current routing and missing native seam

All inspected PAR/QL calls enter the generic `COMMON/sparse_mod.f90` API.
This revision has no FortSparse route. `COMMON/CMakeLists.txt:80` compiles
`umf4_f77wrapper_ver_4_5.c` with DLONG and the complex wrapper with ZLONG;
`:108`/`:110` links legacy SuiteSparse or SuiteSparse::umfpack. The fetched
SuiteSparse version is v7.6.1 (`cmake/suitesparse.cmake:7`). This is the actual
integration route, regardless of older experimental branches.

| QL family | Verified matrix and solve contract | Consequence |
|---|---|---|
| Axi `NEO-2-QL/ripple_solver_axi_test.f90:217` | Real double amat_sp, assembled/compressed to general 1-based CSC (`:2431`), no symmetry guarantee. iopt=1 factor at `:2451`, iopt=2 three source RHS at `:2475`, then correction/propagator RHS, iopt=3 free at `:2937`. | Requires real general LU and retained one/many-RHS solves. Rank/species exchange is outside the factorization. |
| Arnoldi `ripple_solver_ArnoldiOrder2_test.f90:3489` | `problem_type` selects a real factor of DBLE(amat_sp), otherwise a complex factor. In the real branch a complex RHS is split into real/imaginary solves with the same real factor (`:4558`); complex branch directly solves at `:4584`. | Preserve real and complex general LU, exact branch semantics and factor lifetime. Do not reinterpret the real branch as Hermitian Cholesky or factor once per RHS component. |

For eventual integration, the high-level real/complex vector and multiple-RHS
`sparse_solve` overloads are the seam, before UMFPACK's index/value conversion.
They normalize full column indices to CSC and interpret iopt 0/1/2/3 as
factor-solve-free/factor/solve/free. Dispatching only a one-shot solve would lose
retained-factor behavior. Method 2 requests refinement; default method 3 omits
it. The current single module-global factor state is not a reentrant API;
introducing native handles must make ownership explicit and preserve the
caller lifecycle. No wrapper implementation is changed by this audit.

A native route must also preserve default-integer caller validation before
conversion to kind-8 external indices, input overwrite behavior, real versus
complex RHS products and status propagation. Integrating at the UMFPACK C
wrappers would preserve avoidable split-complex copies and bind native policy
to an external ABI. This is a recommendation from the traced layering, not an
implemented route or a claim about future performance.

## MEPHIT remaining sparse paths

`src/mephit_iter.F90:941` is the only production Fortran `sparse_solve` call
found outside `sparse_mod`: complex general FDM, one RHS, iopt omitted. Shared
`sparse_mod.f90:361` normalizes column indices and uses the complex UMFPACK
wrapper; method 2/3 and factor lifecycle mirror the NEO-2 design. Real and
multiple-RHS overloads exist in that module, but their existence does not
establish additional current physics callers.

The optional MFEM MDE path (`src/mephit_fem.cpp:108`) constructs a real block
matrix from K/G/G_minus/K (`src/magnetic_differential_equation.cpp:116`), then
uses UMFPackSolver with METIS and one Mult. This represents coupled real/imaginary
unknowns; positive definiteness or symmetry is not established by the assembly
inspection. The utility `src/miscellaneous.cpp:93` copies a general SparseMatrix,
adds penalty diagonals, finalizes and factors it through another UMFPackSolver.
No in-repository caller of that utility was found. It is a compiled capability,
not a demonstrated additional production workload. The CMake link at `:150`
retains SuiteSparse; no FortSparse native route is present.

## Remaining KIN6D families

These are matrix contracts and retained-factor requirements, not new independent
physics implementations. Default indices are 1-based Fortran integers in CSC;
CSR arrays can be reused as CSC only where the operator's symmetry justifies it.

| Family and source | Class and current realization | Reuse, derivative or accuracy consequence |
|---|---|---|
| Nonnested coupled `src/mhd/nonnested_coupled.f90:1035` | Explicitly nonsymmetric real Jacobian, reduced fixed/tied DOFs, general LU. | Old factor is right preconditioner for FGMRES (`:929`); current residual tolerance 1e-10, otherwise refactor current Jacobian. Cannot replace by SPD factorization. |
| Nonnested equilibrium `src/mhd/equilibrium_nonnested.f90:437` | Masked real SPD factor; complex RHS split into two real solves (`:484`). | Reuses real factor; public complex state does not imply complex matrix factorization. Fixed-data derivatives reuse symmetric solve. |
| Per-n Ampere `src/mhd/mhd_ampere_n.f90:250` | Hermitian positive-definite H=A+iB embedded in doubled real SPD [[A,-B],[B,A]]. Compatibility branch directly factors complex CSC (`:242`). | Real SPD native path can cover current embedding; direct Hermitian native path could avoid doubling but needs independent evidence. Preserve block ordering/conjugation. |
| Equilibrium Riesz `src/equilibrium/equilibrium_residual.f90:352` | Explicit gs_csr_spd_matrix, real SPD, two RHS. Compatibility general SuperLU path at `:366`. | Factor once for both RHS. Generic fallback hint does not change the proven SPD construction. |
| GS majorant `src/equilibrium/gs_majorant.f90:329` | Real energy SPD and recovery systems via GS direct; FortFEM general adapter fallback (`:337`). | Repeated energy inverse/recovery/refinement RHS (`:466`, `:1303`). Retain residual/refinement semantics and caller certificate separation. |
| GS free boundary `src/equilibrium/gs_free_boundary.f90:1390` | Real symmetric-indefinite path only after asymmetry check; otherwise general LU through FortFEM (`:1428`). | Two-column Schur correction block; separate transpose factor for nonsymmetric adjoint at `:1581`. Symmetry hint must not erase general path. |
| Ideal free boundary `src/mhd/ideal_free_boundary.f90:1647` | Symmetric KKT matrices can be indefinite; general multiplier C^T block (`:1933`) and off-diagonal block subproblems (`:2440`). | Old factor preconditions FGMRES; fixed-pattern numeric refactors (`:1671`); general projection factor retained for repeated RHS (`:1950`). Requires both indefinite and general support. |
| Resolved strip/layer `src/mhd/response_strip.f90:259`, `mhd_layer_n.f90:234` | Complex general band ZGBSV, not a sparse backend caller. | One RHS; preserve bandwidth and complex signs. Default layer grid is 2000 then 4000 for extrapolation (`mhd_layer_n.f90:133`), not a presumed tiny fixture. |
| Strip preconditioner `src/mhd/response_strip_preconditioner.f90:170` | Retained complex band ZGBTRF per Fourier mode. | ZGBTRS chooses N or conjugate-transpose C (`:218`, `:242`), alongside FFT/rotation/scaling. Generic sparse coverage does not establish this entire adjoint composition. |
| CSC-only users | nonnested_mapped, nonnested_transport, mhd_ideal_n and quantum_collision_kernel use sparse storage/products in inspected imports/calls. | They do not introduce another direct-factor matrix family. Sparse products remain a separate capability. |

This narrows the former discovery-only gaps. It does not establish global
application coverage of every optional backend or configuration, a certification
result for a new solver, or single-node performance/memory admission.

## Numeric fixture now retained

[NEO-2 mini5](../bench/fixtures/neo2-mini5/README.md) now contains exact integer
MatrixMarket data, its all-ones RHS, a rational solution and nonzero exact
determinant. The independent oracle uses integer signed-permutation determinants
and Cramer's rule, followed by exact rational Ax=b. No solver implementation
was copied; the original MIT license and source revision/hash are retained.
The fixture is a small general nonsymmetric wrapper case, not a kinetic
performance/fill surrogate. Other existing test entry points remain in the
manifest without fabricated application dumps or runtime claims.

A later [private NEO-2 capture and benchmark](../bench/results/private-neo2-2026-10-16/report.md)
adds a real general system with 4,856 unknowns and 28,898 nonzeros, taken before
RHS overwrite from an unchanged NEO-2 binary at revision `7eeee3cd`. Three
physical RHS columns, including their original duplicate, and a manufactured
fourth RHS pass the frozen policy for native LU and all five external backends
across the requested thread budgets. The report distinguishes child timings,
independent exact residual audits and empirical forward comparisons. It establishes one local algebraic
fixture, not whole-application or QL/multispecies admission. Its capture hashes,
private-data boundary and adapter recipe are retained in
[the fixture documentation](../bench/native_families/PRIVATE_REAL_FIXTURES.md).
