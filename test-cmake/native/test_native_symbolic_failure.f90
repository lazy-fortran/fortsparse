program test_native_symbolic_failure
    use, intrinsic :: iso_c_binding, only: c_int
    use fortsparse
    implicit none
    type(csc_t) :: matrix
    type(sparse_solver_t) :: solver
    type(fortsparse_status_t) :: status
    real(dp) :: result(3)
    interface
        subroutine fail_allocation(enabled) bind(C, name='fsparse_test_symbolic_fail')
            import :: c_int
            integer(c_int), value :: enabled
        end subroutine fail_allocation
    end interface

    call csc_from_triplet(3, 3, [1,2,3], [1,2,3], [2.0_dp,3.0_dp,4.0_dp], matrix, status)
    if (.not. status_ok(status)) error stop 'construct'
    solver%backend_id = FORTSPARSE_BACKEND_NATIVE
    solver%options%symmetry = FORTSPARSE_SYM_SPD
    call fail_allocation(1_c_int)
    call sparse_factor(solver, matrix, status)
    if (status%code /= FORTSPARSE_OUT_OF_MEMORY) error stop 'allocation failure status'
    call sparse_solve(solver, [2.0_dp,6.0_dp,12.0_dp], result, status)
    if (status%code /= FORTSPARSE_NOT_FACTORED) error stop 'failed factor invalidated'
    call sparse_destroy(solver)
    call sparse_destroy(solver)
    call fail_allocation(0_c_int)
    call sparse_factor(solver, matrix, status)
    if (.not. status_ok(status)) error stop 'recover after failed allocation'
    call sparse_solve(solver, [2.0_dp,6.0_dp,12.0_dp], result, status)
    if (.not. status_ok(status)) error stop 'recovered solve status'
    if (maxval(abs(result - [1.0_dp,2.0_dp,3.0_dp])) > 1.0e-12_dp) error stop 'recovered values'
    call sparse_destroy(solver)
end program test_native_symbolic_failure
