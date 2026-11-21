# Preconditioned conjugate gradients

`pcg_setup` (or `pcg_setup_csr` on raw arrays) and `pcg_solve` solve
`A x = b` for a real symmetric positive definite CSR matrix, and for a
semidefinite one when `b` lies in the range of `A`; the kernel component of
`x` is then the one of the initial guess. The matrix uses the CSR convention
of [AMG](amg.md): one-based, both triangles, strictly increasing columns.

```fortran
use fortsparse, only: dp, pcg_params_t, pcg_handle_t, fortsparse_status_t, &
    pcg_setup_csr, pcg_solve, pcg_destroy, PCG_IC0
type(pcg_params_t) :: params
type(pcg_handle_t) :: handle
type(fortsparse_status_t) :: info
integer :: iterations
real(dp) :: residual
params%preconditioner = PCG_IC0   ! PCG_NONE, PCG_JACOBI (default), PCG_IC0, PCG_AMG
params%rtol = 1.0e-10_dp          ! stop at ||b - A x|| <= rtol ||b||
params%max_iterations = 1000
call pcg_setup_csr(n, row_ptr, col, val, params, handle, info)
x = 0                              ! initial guess
call pcg_solve(handle, params, b, x, info, iterations, residual)
call pcg_destroy(handle)
```

`info%code` is `FORTSPARSE_OK` when the tolerance is met,
`FORTSPARSE_NOT_CONVERGED` at the iteration limit, and
`FORTSPARSE_NOT_POSITIVE_DEFINITE` when `p^T A p` or `r^T z` is not positive,
as for an inconsistent semidefinite system. IC(0) factors on the pattern of
the lower triangle; a pivot that is not positive restarts it with the
diagonal scaled by `1 + s`, `s` doubling from `1e-3` (Manteuffel 1980), which
semidefinite matrices need. `PCG_AMG` uses the smoothed aggregation of
`params%amg` and needs an SPD matrix.

Iterations to `rtol = 1e-8` on FortFEM's spline curl-curl matrices
(p = 2, two Fourier modes), from zero, Jacobi and IC(0):

| case | n | A + 1e-2 M1 | A, b = A x |
|---|---|---|---|
| torus, level 1 / 2 / 3 | 740 / 1980 / 6260 | 593 / 963 / 1528, 249 / 394 / 624 | 59 / 81 / 116, 21 / 28 / 43 |
| ring torus, level 1 / 2 / 3 | 660 / 1700 / 5220 | 653 / 1093 / 1960, 154 / 369 / 728 | 107 / 160 / 230, 21 / 42 / 71 |

The counts grow with refinement; H(curl) problems need an auxiliary-space
preconditioner.

## AMS for H(curl)

`pcg_setup_ams(A, G, Pi, params, handle, info, mass)` preconditions an
H(curl) matrix `A = C^T M2 C + beta M` (curl-curl, with or without a mass
term) by the auxiliary-space Maxwell solver (Hiptmair and Xu 2007, Kolev and
Vassilevski 2009): l1-Jacobi smoothing on `A`, and corrections in the
gradient space `G` (n1 x n0) and in the nodal vector space `Pi` (n1 x 3 n0,
x, y, z column blocks), each solved by smoothed aggregation AMG. Pass the
1-form mass matrix `M` as `mass` and `beta` as `params%ams%beta`: each
Cartesian component of the nodal space is then solved with the vector
Laplacian of Hiptmair and Xu, `G^T M G + beta Pi_d^T M Pi_d`. Without
`mass`, the component blocks `Pi_d^T A Pi_d` are used; for curl-curl they
lack a derivative, and AMG converges slowly on them. `params%ams` also
selects the multiplicative cycle (`AMS_MULTIPLICATIVE`, default) or the
additive one, and the smoothing sweeps. When `A G = 0` the gradient
correction is left out (`ams_uses_gradient`). The auxiliary operators get
`aux_shift` times their diagonal added, so that a semidefinite `A` gives
definite AMG problems. `ams_setup_csr` and `ams_apply` give the
preconditioner alone.

Iterations to `1e-8` on the curl-curl matrices above, multiplicative cycle,
two smoothing sweeps, with `mass` and `beta` (IC(0) in parentheses):

| case, level 1 / 2 / 3 | A + 1e-2 M1 | A, b = A x |
|---|---|---|
| torus | 19 / 22 / 26 (249 / 394 / 644) | 16 / 17 / 21 (21 / 28 / 43) |
| ring torus | 36 / 41 / 47 (154 / 369 / 727) | 24 / 27 / 30 (21 / 42 / 71) |

Without `mass` the torus needs 30 / 24 / 35 and 19 / 23 / 35, the ring torus
64 / 63 / 76 and 40 / 44 / 54.
