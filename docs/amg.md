# Smoothed aggregation and CSR smoothers

`fortsparse` exports `csr_t`, `amg_params_t`, `amg_handle_t`, `amg_setup`,
`amg_setup_csr`, `amg_apply`, and `amg_destroy`. These APIs are independent of
external sparse backends. Matrices use owned, one-based FP64 CSR with strictly
increasing columns in each row, no duplicates, a stored positive diagonal and
both halves of the symmetric matrix. Setup checks storage, finite values,
symmetry and positive diagonals. The caller must supply an SPD operator:
coarse LDLT inertia is checked, but this is not a proof that the fine operator
is SPD. Singular Neumann problems need constraints before setup.

```fortran
use fortsparse, only: dp, csr_t, amg_params_t, amg_handle_t, &
    fortsparse_status_t, amg_setup, amg_apply, amg_destroy
! Fill A%nrow, A%ncol, A%row_ptr, A%col and A%val before setup.
type(csr_t) :: A
type(amg_params_t) :: params
type(amg_handle_t) :: handle
type(fortsparse_status_t) :: info
! r and z have length A%nrow, and are separate caller-owned arrays.
! Optional: params%near_nullspace(n,k) = caller's k candidate modes.
call amg_setup(A, params, handle, info)
! Check info%code == 0 before applying the preconditioner.
call amg_apply(handle, r, z, info)
call amg_destroy(handle)
```

For a consumer such as [FortFEM](https://github.com/lazy-fortran/fortfem),
`amg_setup_csr(n, row_ptr, columns, values, params, handle, info)` accepts raw
arrays. `amg_apply(handle, r, z)` is the procedure-pointer-free PCG hook, with
an optional fourth status argument. It starts every application from zero,
returns a fixed linear correction and allocates no numerical workspace.
Check the optional status when errors must be reported. An invalid handle or
shape returns zero output and a failing status. A handle owns mutable workspace,
so concurrent calls need separate handles. Repeated setup replaces the old
hierarchy, and destroy is idempotent. There is no global registry or hidden
external-backend fallback.

## Setup and cycles

An edge is strong when
`abs(A(i,j))/sqrt(A(i,i))/sqrt(A(j,j)) >= strength_threshold`, excluding
self-edges and stored zeros. The default threshold is 0.02. A greedy seed pass
forms disjoint strong-neighbor aggregates, and a second pass attaches leftover
vertices. Graph selection is ordered and serial. Strength construction, local
QR, sparse transpose/products and row kernels use OpenMP.

Each aggregate uses twice-reorthogonalized modified Gram-Schmidt on the caller's
candidate modes. Dependent local columns are dropped at relative norm 1e-12.
A zero local candidate space receives a normalized constant. Local ranks may
differ, so even redundant candidate vectors are supported. The tentative
prolongator T and propagated candidates Bc satisfy T Bc = B up to the local
rank tolerance. When no modes are supplied, B is the constant vector.

The prolongator is `P = (I - omega D^-1 A) T`, where
`omega = prolongator_damping / max_i(sum_j(abs(A(i,j)))/A(i,i))` and the default
`prolongator_damping` is 4/3. Restriction is P transpose. Two sparse products
form RAP, whose paired entries are averaged to remove asymmetric roundoff.
No interpolation truncation or anisotropy-dependent parameter fitting is used.
Coarsening stops at 64 rows, 20 levels, or when the tentative space does not
reduce the dimension. The existing `native_ldlt_t` factors the last operator.
Noncoarsening problems can therefore make this coarse solve expensive.

Both `cycle=1` (V, default) and `cycle=2` (W) use the same fixed polynomial
before and after coarse correction. These symmetric cycles are intended for
PCG. Defaults are two Chebyshev steps per smoothing application and a spectral
interval ratio of 20. All hierarchy matrices and vector workspaces are retained.
There is no sparse fill truncation or implemented strict hierarchy allocation
cap. `memory_limit_bytes /= 0` fails with `FORTSPARSE_RESOURCE_LIMIT` before
hierarchy allocation, including the raw-array bridge. This is an explicit
unsupported path, not a claimed memory bound.

## Standalone smoothers

`csr_smoother_t` is set up with
`smoother_setup(A, method, handle, info, degree, ratio)` and applied with
`smoother_apply(A, handle, b, x, info)`. Apply updates the supplied iterate x.
Use `SMOOTHER_L1_JACOBI` or `SMOOTHER_CHEBYSHEV`, and release workspace with
`smoother_destroy(handle)`. The matrix must remain unchanged between setup and
apply. These routines assume SPD and validate canonical storage and positive
diagonal entries during setup.

l1-Jacobi scales the residual by the inverse absolute row sum, with unit
relaxation. Chebyshev uses diagonal scaling and the three-term recurrence on
`[upper/ratio, upper]`. Twelve power steps estimate the largest eigenvalue of
the symmetrically scaled operator. The reported `eigenvalue_estimate()` is a
diagnostic, whereas `eigenvalue_bound()` supplies a conservative Gershgorin
upper bound for the polynomial. Using the bound avoids turning an
underestimated Ritz value into an unstable polynomial. Both smoothers own
reusable residual and direction vectors and thread their row updates.

## Verification and measurements

`test_fortsparse_csr_smoothers` checks tentative candidate reproduction with a
redundant mode, sparse transpose/products and RAP against dense operations,
known spectral bounds, and standalone residual contraction.
`test_fortsparse_amg` checks 2D/3D P1-like and Q1 operators, V/W symmetry and
positivity, multi-mode/raw-array use, setup failure, lifetime, manufactured
solutions and LDLT agreement to 1e-10. On Linux, `test_amg_allocations` observes
five repeated applications of each cycle through the existing malloc/calloc/
realloc wrapper after OpenMP warmup. The ordinary tests have a one-second
CTest budget.

The synthetic P1-like generator uses Dirichlet 5/7-point graphs. Q1 uses exact
tensor products of 1D stiffness and consistent mass matrices, with boundary
degrees of freedom eliminated. Benchmarks use these published formulas without
machine or experiment data. `bench_amg SIDE DIM FAMILY EPS THETA SOLVER MAXITER`
compares AMG-PCG, ILUTP-PCG and native LDLT. Family 1 is P1-like, family 2 is Q1.
ILUTP uses fill 12/12, drop tolerance 1e-3 and pivot tolerance zero. It is
not generally symmetric, so its use in PCG is an empirical requested comparison,
not a guarantee of PCG applicability. All iterative runs use the same true
relative residual target of 1e-12 and report convergence failures.

Remote convergence and timing tables are in [AMG results](../bench/AMG_RESULTS.md).
