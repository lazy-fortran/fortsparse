# User guide

Build with the [CMake or fpm commands](building.md), then factor a matrix,
solve one or more right-hand sides and release its factors. This complete
program solves a symmetric positive-definite system using the native CPU
backend:

```fortran
program first_solve
    use fortsparse, only: dp, csc_t, sparse_solver_t, fortsparse_status_t, &
        csc_from_triplet, sparse_factor, sparse_solve, sparse_free, status_ok, &
        FORTSPARSE_BACKEND_NATIVE, FORTSPARSE_SYM_SPD
    implicit none
    type(csc_t) :: matrix
    type(sparse_solver_t) :: solver
    type(fortsparse_status_t) :: status
    real(dp) :: solution(3)

    call csc_from_triplet(3, 3, [1, 1, 2, 2, 2, 3, 3], &
        [1, 2, 1, 2, 3, 2, 3], &
        [4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], &
        matrix, status)
    if (.not. status_ok(status)) error stop "matrix construction failed"
    solver%backend_id = FORTSPARSE_BACKEND_NATIVE
    solver%options%symmetry = FORTSPARSE_SYM_SPD
    call sparse_factor(solver, matrix, status)
    if (.not. status_ok(status)) error stop "factorization failed"
    call sparse_solve(solver, [6.0_dp, 10.0_dp, 8.0_dp], solution, status)
    if (.not. status_ok(status)) error stop "solve failed"
    print *, solution
    call sparse_free(solver)
end program first_solve
```

The solution is `[1, 2, 3]`. Triplet row and column indices are one-based.
`csc_from_triplet` sorts entries and combines duplicates into compressed sparse
column storage. Supply the full symmetric matrix for the native
positive-definite path, not just one triangle.

Check `status` after every call that returns it. `status%message` explains
failures, including invalid matrices, singular or nonpositive pivots,
unavailable backends and resource limits. Select a backend before factoring;
automatic selection depends
on the configured build. An explicit unavailable backend reports an error.

A solver retains its factors until `sparse_free` or `sparse_destroy`. Repeated
`sparse_solve` calls reuse them; `sparse_refactor` reuses symbolic work only
when its same-pattern contract holds. Use the supported copy operations rather
than copying private backend storage; see the ownership restrictions in
[native matrix classes](native-cpu.md).

For derivative actions, the tangent call solves the linearized equation and
the adjoint call uses transpose or conjugate-transpose factors. They compute
derivatives of the specified solve equation; they do not differentiate an
application's matrix assembly automatically. See [the derivative contract and
examples](backends.md#implicit-derivatives).
