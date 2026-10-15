program test_fortsparse_mumps_mpi_lifecycle
    use fortsparse, only: dp, csc_t, sparse_solver_t, sparse_factor, sparse_solve, &
        sparse_destroy, sparse_refactor, fortsparse_status_t, status_ok, &
        FORTSPARSE_BACKEND_MUMPS, FORTSPARSE_BACKEND_UNAVAILABLE
    use fs_test_mpi, only: test_mpi_init, test_mpi_finalize
    implicit none

    type(csc_t) :: A
    type(sparse_solver_t) :: solver, retained, retained_factor, retained_free
    type(fortsparse_status_t) :: status
    real(dp) :: b(2), x(2)
    logical :: reused

    A%nrow = 2
    A%ncol = 2
    A%nnz = 2
    A%col_ptr = [1, 2, 3]
    A%row_idx = [1, 2]
    A%val = [2.0_dp, 3.0_dp]
    b = [4.0_dp, 9.0_dp]
    solver%backend_id = FORTSPARSE_BACKEND_MUMPS
    call sparse_factor(solver, A, status)
    if (status%code /= FORTSPARSE_BACKEND_UNAVAILABLE) &
        error stop "factor before MPI_Init must return unavailable"
    call test_mpi_init()
    call sparse_factor(solver, A, status)
    if (.not. status_ok(status)) error stop "factor after MPI_Init failed"
    call sparse_solve(solver, b, x, status)
    if (.not. status_ok(status)) error stop "solve after MPI_Init failed"
    if (maxval(abs(x - [2.0_dp, 3.0_dp])) > 1.0e-14_dp) &
        error stop "incorrect MPI MUMPS solution"
    retained%backend_id = FORTSPARSE_BACKEND_MUMPS
    retained_factor%backend_id = FORTSPARSE_BACKEND_MUMPS
    retained_free%backend_id = FORTSPARSE_BACKEND_MUMPS
    call sparse_factor(retained, A, status)
    if (.not. status_ok(status)) error stop "retained factor failed"
    call sparse_factor(retained_factor, A, status)
    if (.not. status_ok(status)) error stop "replacement factor failed"
    call sparse_factor(retained_free, A, status)
    if (.not. status_ok(status)) error stop "late free factor failed"
    call sparse_destroy(solver)
    call test_mpi_finalize()
    call sparse_solve(retained, b, x, status)
    if (status%code /= FORTSPARSE_BACKEND_UNAVAILABLE) &
        error stop "retained solve after MPI_Finalize must return unavailable"
    call sparse_refactor(retained, A, status, reused)
    if (status%code /= FORTSPARSE_BACKEND_UNAVAILABLE .or. reused) &
        error stop "retained refactor after MPI_Finalize must return unavailable"
    call sparse_factor(retained_factor, A, status)
    if (status%code /= FORTSPARSE_BACKEND_UNAVAILABLE) &
        error stop "retained factor after MPI_Finalize must return unavailable"
    call sparse_destroy(retained_free)
    call sparse_factor(solver, A, status)
    if (status%code /= FORTSPARSE_BACKEND_UNAVAILABLE) &
        error stop "factor after MPI_Finalize must return unavailable"
end program test_fortsparse_mumps_mpi_lifecycle
