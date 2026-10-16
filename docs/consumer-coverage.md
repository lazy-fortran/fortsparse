# ITPplasma sparse consumer coverage

Source audit: 2026-10-16. This records inspected call contracts, not completed
application validation. [The fixture manifest](consumer-fixtures.json) pins each
local repository revision and distinguishes inspected paths from discovery-only
matches. No application was run, no matrix was exported, and no cluster host was
accessed. The local checkouts were read without modification.

Subsequent work captured a private real NEO-2 system from an unchanged binary
at a newer pinned revision. The [aggregate benchmark](../bench/results/private-neo2-2026-10-16/report.md)
records its 4,856 unknowns, four RHS, frozen accuracy policy and five-peer
comparison. It supplements this source audit; it does not validate a complete
NEO-2 application run or the other consumer families below. Numeric inputs and
solutions remain private.

## Priority and admission boundary

1. **Real general sparse LU, then complex general sparse LU**, with threshold
   pivoting, retained factors and many RHS, directly address NEO-2 and MEPHIT.
   NEO-2 makes memory a first-class acceptance criterion. SPD-only code does not
   cover these consumers; the inspected assembly does not establish symmetry.
2. **Real SPD simplicial Cholesky** is a smaller independent native slice for
   KIN6D's GS path. It is useful but cannot stand in for ecosystem LU coverage.
3. **Transpose/conjugate-transpose and symbolic/numeric reuse** must be explicit
   contracts. GPEC already exercises conjugate-transpose band solves. KIN6D's
   indefinite Newton systems retain symbolic analysis across value changes.
4. Keep mature **pivoted indefinite, large parallel and distributed** backends.
   JOREK and M3D-C1 require communicator/ownership semantics beyond the current
   local factor-once API. SPEC's ILUTP and GMRES are a separate sparse
   preconditioning/iterative need, not evidence for a Cholesky replacement.

Native sparse code may call BLAS/LAPACK for dense blocks. Exploiting a proven
band structure can be preferable to general sparse LU: NEO-2's spline path and
STELLOPT/THRIFT already do this. No timing, fill or memory superiority follows
from implementation language or a shared API.

## Inspected contracts

References below are repository-relative `file:line` at the manifest revision.
Default Fortran `INTEGER` width is compiler-dependent; “kind 8” is the actual
source spelling, not an assertion that all builds use the same ABI. No sparse
index-width or allocation-limit capability should be inferred from a solver name.

| Consumer and inspected path | Matrix, storage and indices | RHS, lifetime and transpose | Accuracy and memory evidence |
|---|---|---|---|
| NEO-2 PAR, `NEO-2-PAR/ripple_solver.f90:2251` | Real double general system; COO assembly sorted/compressed to 1-based CSC; default-integer caller arrays. No SPD proof found. | `:2271` factors, `:2283` solves 3 source columns, `:2302` repeatedly solves integral corrections; `:2762` frees. | Relative integral-iteration stopping at `:2309` is not a linear forward-error bound. Full LU fill and peak RSS unknown. |
| NEO-2 QL, `NEO-2-QL/ripple_solver_ArnoldiOrder2_test.f90:3493` | Calls both real `DBLE(amat_sp)` and complex `amat_sp` branches. General CSC contract; branch conditions and entire QL solver not audited. | Factor/solve stages at `:3493`, `:3698`; additional solves discovered at `:4561`, `:4584`. No transpose call established. | Real conversion can create a temporary. Large QL production sizes and branch-specific error requirements unknown. |
| NEO-2 shared wrapper, `COMMON/sparse_mod.f90:651` | Accepts full column indices or CSC; converts 1-based default integers to zero-based `INTEGER(kind=8)` Ai/Ap at `:1120`. Real and split-complex UMFPACK paths. | iopt 1 factor, 2 solve, 3 free; module-global factor handle. `sys=0` at `:31`. Multiple RHS wrapper loops columns at `:1321`. Symbolic and numeric retained until free. | Method 2 requests refinement, default method 3 omits it (`:1148`); no certificate. Copies/workspaces allocated per call, including reuse calls. |
| NEO-2 spline, `COMMON/spline_cof.f90:175` | Common switches `(2,4)` use permuted real band matrix, lower/upper bandwidth 17/14. Other switches use sparse assembly. | `:679` calls DGBSV once; off-band insertion rejected at `:670`. | Eight-node reconstruction regression exists; this does not sample kinetic-LU fill. |
| MEPHIT FDM, `src/mephit_iter.F90:907` | Complex double matrix of dimension `mesh%npoint-1`, COO/full column indices passed to `sparse_mod`; default integers. No Hermitian/SPD claim established. | One RHS overwritten, no iopt supplied at `:941`, so wrapper factor/solve/free per call. No transpose found in inspected path. | Explicit residual at `:945` and relative-error logging at `:957`; no enclosing bound or rejection threshold established there. Matrix values saved by FDM I/O, smallest portable dataset unknown. |
| MEPHIT MFEM branch, `src/mephit_fem.cpp:108` | MFEM real SparseMatrix for MDE, assembled externally; CSR representation managed by MFEM, internal index ABI not audited. | UMFPackSolver with METIS ordering; SetOperator and one Mult at `:113`. | Real assembly uses complex input fields; algebraic class not proven. Framework lifetime/copy cost and peak unknown. |
| KIN6D GS, `src/fem2d/gs_direct.f90:36` | Full real symmetric CSR, interpreted as CSC by symmetry; default SPD, optional indefinite. Default integer indices. | Cached factor and repeated solves. Fixed thread count at factor/solve. SPD failure is recorded without retry (`:40`). | Existing element-action and iterative oracles; no native sparse certificate follows from factor success. Statistics requested at `:72`. |
| KIN6D ideal map, `src/mhd/ideal_map.f90:1583` | Damped real symmetric system explicitly may be indefinite. Full sparse input; ordering/pivoting backend-dependent. | `:1605` value refactor with fixed pattern; symbolic reuse through MUMPS. | Indefinite pivoting required. Whole-application accuracy, peak and asymptotic fill not measured by this audit. |
| GPEC singular-layer path, `dcon/sing.f:949` | Complex general band LU; LAPACK storage `3*mband+1`, default-integer pivots. | `:998` normal and `:1005` conjugate-transpose ZGBTRS, multiple RHS `mpert`; factors reused. | No sparse migration need established: preserve band and conjugation semantics. Tolerance, largest bandwidth and memory budget unknown. |
| MARS-Q block path, `MarsQ_2FK/pams.f:2225` | Complex dense block LU via ZGETRF; no general CSC/COO path found by discovery scan. | `:2233` two normal ZGETRS calls with N RHS; factors reused. | Block recurrence outside these routines not audited; no sparse-replacement case established. |
| JOREK MUMPS adapter, `solvers/mod_mumps.f90:28` | DMUMPS real general (`SYM=0`); central or distributed assembled COO via irn/jcn at `:72`, `:90`. `datatypes/mod_integer_types.f90:8` selects int32/int64 with INTSIZE64; MUMPS ABI consistency not checked. | Separate analysis JOB=1, numeric JOB=2, solve JOB=3; vector and multiple-RHS entry points (`:169`); communicator and broadcast ownership explicit. | BLR accuracy parameter epsilon_BLR at `:85`; ICNTL workspace percentages are relaxation choices, not measured peaks. Other solver adapters remain uninspected. |
| M3D-C1 PETSc path, `unstructured/petsc_matrix.f90:106` | Distributed MPIAIJ, local/global DOF counts, real/complex conditional insertion; zero-based MatSetValues indices (`:298`). PETSc index width build-dependent. | In-place KSPSolve at `:325`; solver/context retained. SCOREC path also uses PetscInt and MPIAIJ (`m3dc1_scorec/src/m3dc1_matrix.cc:677`). | Block preallocation at `:139`; SCOREC has multiple Krylov/preconditioner policies. No unique LU or tolerance policy established for the application. |
| SPEC Beltrami, `src/mp00ac.f90:377` | Direct path is dense real LU with condition estimate/refinement; GMRES path builds sparse ILUTP from matrixS/jdMAS/idMAS (`:437`). | Dense primal factor reused for 2 derivative RHS (`:397`); sparse preconditioner reused by GMRES. | DGECON, DGERFS and `rcond>=machprec`; epsILU controls approximation. CSR-like preconditioner inputs inferred from ILUTP call, row-index convention not independently audited. |
| STELLOPT THRIFT, `THRIFT/Sources/thrift_plasma_solver_mod.f90:1082` | Routine named “sparse” actually calls dense DGESV; `:1116` permutes species-major dense input to band with kl=ku=num_species and calls DGBSV (`:1156`). Density uses DGTSV (`:1068`). | One RHS, no persistent factor reuse established. | NaN/error checks; band result is a structural optimization, not a general sparse consumer. Original dense input remains resident. |

The [single-node follow-up](consumer-single-node.md) traces QL axi/Arnoldi,
MEPHIT and remaining KIN6D direct/band families. Its newer classifications
supersede the corresponding discovery-only gaps below; other limits remain.

## Coverage status and remaining inspection

“Scoped no need” means no general sparse direct requirement was established in
the inspected production family, not proof that every file, submodule or optional
build lacks one. Discovery searches included common sparse-library names and
factor/solve APIs across Fortran/C/C++ source; false positives and bundled
libraries were not counted as application use.

| Repository | Status | Evidence and limit |
|---|---|---|
| NEO-2 | Inspected primary PAR/shared wrapper; partial QL | Other QL axi branch discovered, not traced. No FortSparse call in inspected current wrapper; earlier branch history is not current behavior. |
| MEPHIT | Partial consumer coverage | FDM and MFEM calls inspected; other physics and sparse_mod variants not fully traced. |
| KIN6D | Partial consumer coverage | GS and ideal_map inspected; additional sparse callers in response_strip, nonnested_coupled, mhd_layer_n, mhd_ampere_n, ideal_free_boundary, equilibrium_nonnested and equilibrium_residual need family-specific audit. |
| JOREK | Partial consumer coverage | MUMPS adapter inspected; solve_mumps_all, legacy mumps_module, PaStiX, STRUMPACK, WSMP and application assembly/call graph uninspected. |
| M3D-C1 | Partial consumer coverage | Fortran PETSc lifecycle and selected SCOREC allocation/solver calls inspected; other matrix backends and options not exhaustively traced. |
| GPEC | Partial structured coverage | Complex band solves inspected in dcon/sing; band calls also discovered in equil/spline, equil/bicube, dcon/fourfit, rmatch/deltac, rdcon/gal; those families uninspected. |
| MARS-Q | Partial structured coverage | Dense block factor/solve inspected; full recurrence and all versions uninspected. |
| SPEC | Partial direct/preconditioner coverage | mp00ac inspected; other matrix families uninspected. |
| STELLOPT | Partial structured coverage | THRIFT inspected; VMEC2000 blocktridiagonalsolver_bst imports dense BLAS/LAPACK timing and parallel block machinery; its numerical recurrence remains uninspected. Other STELLOPT applications uninspected. |
| VMEC standalone | Uninspected/not located | No independent `/home/ert/code/VMEC` checkout. STELLOPT/VMEC2000 is the available scoped pointer, not a claim about all VMEC implementations. |
| libneo | Scoped no general sparse need | `src/magfie/analytical_gs_solver.f90:396` uses 7x7 dense DGESV; `src/interpolate/acc_spline_build_order5.f90:98` implements Thomas/cyclic tridiagonal solves. Other optional families not semantically audited. |
| SIMPLE | Scoped no general sparse need | `src/orbit_symplectic.f90:738` and diagnostic Newton paths use small dense DGESV; `src/orbit_rk_core.f90:169` supplies device dense LU. Orbit physics not a sparse-LU benchmark. |
| GORILLA | Scoped no general sparse need | `SRC/differentiate.f90:53` calls DGESV(3,3) for local differentiation. Remaining source only discovery-scanned. |
| NEO-RT | Scoped no general sparse need | `POTATO/SRC/equimoments.f90:158` calls dense DGESV with two RHS. sample_matrix files are sampling objects, not evidence of a sparse linear solver. Remaining families uninspected. |

## NEO-2 memory and cluster sizing

The PAR dimension is computed at `ripple_solver.f90:853`; COO capacity is
counted before allocating irow, icol, amat_sp, ipcol and bvec_sp (`:1884`).
Compression reports actual n and nnz at `:2253`. Source and flux arrays alone
hold six real vectors (`:1157`), with other state and auxiliary vectors present.

For real double input with caller integer width I bytes, the assembled COO
row/column/value arrays occupy `(2*I+8)*nnz` bytes before spare capacity and
compression scratch. CSC pointers add `I*(n+1)`. The real wrapper allocates
kind-8 Ai/Ap plus a real solution vector on every call (`sparse_mod.f90:1103`),
adding `8*(nnz+n+1)+8*n` bytes on the conventional kind-8 ABI. These are partial
accounting terms, not a peak-memory prediction. Caller capacities can exceed
compressed nnz, and factors/workspace/state dominate some regimes.

The complex wrapper additionally splits values into valx/valz and has xx/xz,
bx/bz (`:1191`): `16*nnz+32*n` bytes at 64-bit reals, in addition to kind-8 CSC
copies and the caller's original complex values. Reuse calls still allocate
these conversion arrays. QL `DBLE(amat_sp)` branches can add another value
conversion. Symbolic analysis and numeric handles survive until explicit free;
multiple simultaneous handles and wrapper reentrancy are not supported by the
inspected module-global design.

The complex diagnostic prints symbolic upper estimates and numeric LU size,
memory and fill (`:1230`, `:1253`) when sparse_talk is enabled. These diagnostic
formulas should be checked against the actual linked UMFPACK version before
using them as an authoritative byte count. The real path does not print the
same detailed block. `DOC/UserDoc/Input.tex:51` states at-least-cubic growth with
Laguerre order; this is application guidance, not a fitted scaling law or a
proved total-memory model.

Cluster templates request `21*1024`, `40*1024`, or `50*1024` memory units in
`DOC/condor_submit_example:17`, `condor_submit_all_reconstructions_example:17`,
and `condor_submit_failed.template:12`. They are examples, not measured minimum
RAM. Some request multiple processes while setting OMP_NUM_THREADS=1. The
communicator/process replication pattern must be measured before converting
per-solver bytes into per-job RAM. No protected host was contacted.

Before native LU promotion, collect n, allocated capacity, compressed nnz,
index widths, L/U nnz, symbolic and numerical workspace estimates, peak RSS,
RHS count, conversion/setup time and solver residual for an owner-approved
small PAR/QL case and at least one representative larger case. Record MPI ranks,
threads, BLAS, accuracy and any scheduler memory headroom separately. Never
extrapolate cluster capacity from nnz alone or run a production case merely to
populate this audit.

## Reproducible small fixtures

The manifest lists existing constructors/tests, not new application dumps.
NEO-2's [5x5 mini numeric fixture](../bench/fixtures/neo2-mini5/README.md)
and eight-node spline regression are available without
a kinetic run. Neither substitutes for a production kinetic fill distribution.
KIN6D's manufactured GS regression constructs a mesh internally, but its full
multi-mesh nonlinear runtime was not measured here. MEPHIT/JOREK/M3D-C1/GPEC
have no inspected portable tiny application matrix in this audit; fixture
admission there remains blocked on source-specific extraction and independent
oracle selection. Do not label invented Poisson matrices as those applications.
