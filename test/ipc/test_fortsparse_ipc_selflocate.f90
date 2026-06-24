program test_fortsparse_ipc_selflocate
    ! Exercise exe-directory discovery of the GPL helper. This executable is
    ! built into the same directory as fortsparse_umfpack_helper, and ctest runs
    ! it WITHOUT FORTSPARSE_UMFPACK_HELPER in the environment. The out-of-process
    ! UMFPACK backend must still resolve the helper sitting next to the program,
    ! with no PATH entry and no environment variable. Solves a known real system
    ! and checks the solution.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: dp, csc_t, csc_from_triplet, &
        sparse_solver_t, sparse_factor, sparse_solve, sparse_free, &
        fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_UMFPACK_IPC
    implicit none

    integer :: nfail
    nfail = 0

    call run_real(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    ! Real system A x = b solved through the helper found next to this program.
    subroutine run_real(nfail)
        integer, intent(inout) :: nfail

        type(csc_t)               :: A
        type(sparse_solver_t)     :: solver
        type(fortsparse_status_t) :: status
        integer                   :: rows(7), cols(7)
        real(dp)                  :: vals(7)
        real(dp)                  :: b(3), x(3), xe(3)

        rows = [1, 1, 2, 2, 2, 3, 3]
        cols = [1, 2, 1, 2, 3, 2, 3]
        vals = [4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]
        call csc_from_triplet(3, 3, rows, cols, vals, A, status)
        call check_true("real_build", status_ok(status), nfail)

        solver%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC
        xe = [1.0_dp, 2.0_dp, 3.0_dp]
        b = [6.0_dp, 10.0_dp, 8.0_dp]
        call sparse_factor(solver, A, status)
        call check_true("real_factor", status_ok(status), nfail)

        call sparse_solve(solver, b, x, status)
        call check_true("real_solve", status_ok(status), nfail)
        call check_err("real_x", x, xe, nfail)
        call sparse_free(solver)
    end subroutine run_real

    subroutine check_err(label, got, want, nfail)
        character(*), intent(in)    :: label
        real(dp),     intent(in)    :: got(:), want(:)
        integer,      intent(inout) :: nfail
        if (maxval(abs(got - want)) > 1.0e-9_dp) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,es13.6)") "FAIL [", label, &
                "] max|err| ", maxval(abs(got - want))
        end if
    end subroutine check_err

    subroutine check_true(label, cond, nfail)
        character(*), intent(in)    :: label
        logical,      intent(in)    :: cond
        integer,      intent(inout) :: nfail
        if (.not. cond) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a)") "FAIL [", label, "] expected .true."
        end if
    end subroutine check_true

end program test_fortsparse_ipc_selflocate
