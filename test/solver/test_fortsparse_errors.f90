program test_fortsparse_errors
    ! Error paths on the solver: a singular matrix reports FORTSPARSE_SINGULAR,
    ! a solve before factor reports FORTSPARSE_NOT_FACTORED, an unknown backend
    ! id reports FORTSPARSE_BACKEND_UNAVAILABLE, and the UMFPACK_IPC backend
    ! reports the same when no helper is present (the build under fpm/fo ships
    ! no helper). None of these crash.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: dp, csc_t, csc_from_triplet, sparse_solver_t, &
        sparse_factor, sparse_solve, sparse_free, fortsparse_status_t, &
        FORTSPARSE_SINGULAR, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_BACKEND_UMFPACK_IPC
    implicit none

    integer :: nfail
    nfail = 0

    call test_singular(nfail)
    call test_not_factored(nfail)
    call test_unknown_backend(nfail)
    call test_umfpack_ipc_no_helper(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    ! A matrix with a structurally empty column is singular.
    subroutine test_singular(nfail)
        integer, intent(inout) :: nfail

        type(csc_t)               :: A
        type(sparse_solver_t)     :: solver
        type(fortsparse_status_t) :: status
        integer                   :: rows(2), cols(2)
        real(dp)                  :: vals(2)

        ! Column 2 has no entries -> singular.
        rows = [1, 3]
        cols = [1, 3]
        vals = [1.0_dp, 1.0_dp]
        call csc_from_triplet(3, 3, rows, cols, vals, A, status)
        call sparse_factor(solver, A, status)
        call check_code("singular", status%code, FORTSPARSE_SINGULAR, nfail)
        call sparse_free(solver)
    end subroutine test_singular

    ! Solving before factoring must report not-factored.
    subroutine test_not_factored(nfail)
        integer, intent(inout) :: nfail

        type(sparse_solver_t)     :: solver
        type(fortsparse_status_t) :: status
        real(dp)                  :: b(2), x(2)

        b = [1.0_dp, 2.0_dp]
        call sparse_solve(solver, b, x, status)
        call check_code("not_factored", status%code, FORTSPARSE_NOT_FACTORED, &
            nfail)
    end subroutine test_not_factored

    ! An out-of-range backend id reports unavailable, not a crash.
    subroutine test_unknown_backend(nfail)
        integer, intent(inout) :: nfail

        type(csc_t)               :: A
        type(sparse_solver_t)     :: solver
        type(fortsparse_status_t) :: status
        integer                   :: rows(1), cols(1)
        real(dp)                  :: vals(1)

        rows = [1]; cols = [1]; vals = [1.0_dp]
        call csc_from_triplet(1, 1, rows, cols, vals, A, status)
        solver%backend_id = 99
        call sparse_factor(solver, A, status)
        call check_code("unknown_backend", status%code, &
            FORTSPARSE_BACKEND_UNAVAILABLE, nfail)
        call sparse_free(solver)
    end subroutine test_unknown_backend

    ! Selecting the UMFPACK_IPC backend without a helper reports unavailable.
    subroutine test_umfpack_ipc_no_helper(nfail)
        integer, intent(inout) :: nfail

        type(csc_t)               :: A
        type(sparse_solver_t)     :: solver
        type(fortsparse_status_t) :: status
        integer                   :: rows(1), cols(1)
        real(dp)                  :: vals(1)

        rows = [1]; cols = [1]; vals = [1.0_dp]
        call csc_from_triplet(1, 1, rows, cols, vals, A, status)
        solver%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC
        call sparse_factor(solver, A, status)
        call check_code("umfpack_ipc_no_helper", status%code, &
            FORTSPARSE_BACKEND_UNAVAILABLE, nfail)
        call sparse_free(solver)
    end subroutine test_umfpack_ipc_no_helper

    subroutine check_code(label, got, want, nfail)
        character(*), intent(in)    :: label
        integer,      intent(in)    :: got, want
        integer,      intent(inout) :: nfail
        if (got /= want) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,i0,a,i0)") "FAIL [", label, "] got ", &
                got, " want ", want
        end if
    end subroutine check_code

end program test_fortsparse_errors
