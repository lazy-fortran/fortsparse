program test_fortsparse_mumps_mpi_ranks
    ! Ranks factor different systems at different times. MUMPS must run on a
    ! per-process communicator; a shared one would hang or mix the systems.
    use mpi, only: MPI_COMM_WORLD, MPI_Comm_rank, MPI_Barrier
    use fortsparse, only: dp, sparse_solver_t, sparse_factor, sparse_solve, &
        sparse_destroy, fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_MUMPS
    use fs_test_mpi, only: test_mpi_init, test_mpi_finalize
    implicit none

    type(sparse_solver_t) :: solver
    type(fortsparse_status_t) :: status
    real(dp) :: b(2), x(2), d(2)
    integer :: rank, ierr, rep

    call test_mpi_init()
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    solver%backend_id = FORTSPARSE_BACKEND_MUMPS
    ! Rank r solves diag(2+r, 3+r) x = d*[1, 2], and does it r+1 times.
    d = [2.0_dp + rank, 3.0_dp + rank]
    do rep = 1, rank + 1
        call sparse_factor(solver, 2, 2, 2, [1, 2, 3], [1, 2], d, status)
        if (.not. status_ok(status)) error stop "FAIL: factor"
        b = d*[1.0_dp, 2.0_dp]
        call sparse_solve(solver, b, x, status)
        if (.not. status_ok(status)) error stop "FAIL: solve"
        if (maxval(abs(x - [1.0_dp, 2.0_dp])) > 1.0e-14_dp) &
            error stop "FAIL: wrong per-rank solution"
    end do
    call sparse_destroy(solver)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call test_mpi_finalize()
    if (rank == 0) print *, "PASS"
end program test_fortsparse_mumps_mpi_ranks
