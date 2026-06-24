program test_fortsparse_umfpack_ipc
    ! Round-trip the out-of-process UMFPACK backend through the GPL helper. The
    ! ctest harness sets FORTSPARSE_UMFPACK_HELPER to the built helper. Selects
    ! FORTSPARSE_BACKEND_UMFPACK_IPC and solves a known real and complex system
    ! over shared memory, checking the solution and residual.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: dp, csc_t, csc_z_t, csc_from_triplet, csc_matvec, &
        sparse_solver_t, sparse_factor, sparse_solve, sparse_free, &
        fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_UMFPACK_IPC
    implicit none

    integer :: nfail
    nfail = 0

    if (.not. helper_present()) then
        ! No GPL helper in this build (the MIT-only fpm/fo path). The
        ! out-of-process backend is unavailable by design, so there is nothing
        ! to round-trip. ctest builds the helper and sets the variable.
        write (*, "(a)") "SKIP: FORTSPARSE_UMFPACK_HELPER not set"
        write (*, "(a)") "PASS"
        stop 0
    end if

    call run_real(nfail)
    call run_complex(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    ! True when the ctest harness points at the built GPL helper binary.
    logical function helper_present()
        character(4096) :: buf
        integer         :: n

        call get_environment_variable("FORTSPARSE_UMFPACK_HELPER", buf, n)
        helper_present = (n > 0)
    end function helper_present

    ! Real system A x = b solved through the helper, second RHS reuses factors.
    subroutine run_real(nfail)
        integer, intent(inout) :: nfail

        type(csc_t)               :: A
        type(sparse_solver_t)     :: solver
        type(fortsparse_status_t) :: status
        integer                   :: rows(7), cols(7)
        real(dp)                  :: vals(7)
        real(dp)                  :: b(3), x(3), x2(3), xe(3), xe2(3)

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

        xe2 = [0.0_dp, 1.0_dp, 0.0_dp]
        call sparse_solve(solver, csc_matvec(A, xe2), x2, status)
        call check_true("real_solve2", status_ok(status), nfail)
        call check_err("real_x2", x2, xe2, nfail)
        call sparse_free(solver)
    end subroutine run_real

    ! Complex system A x = b solved through the helper.
    subroutine run_complex(nfail)
        integer, intent(inout) :: nfail

        type(csc_z_t)             :: A
        type(sparse_solver_t)     :: solver
        type(fortsparse_status_t) :: status
        integer                   :: rows(5), cols(5)
        complex(dp)               :: vals(5), b(3), x(3), xe(3)

        rows = [1, 2, 2, 3, 3]
        cols = [1, 1, 2, 2, 3]
        vals = [cmplx(2.0_dp, 1.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
            cmplx(3.0_dp, -1.0_dp, dp), cmplx(1.0_dp, 1.0_dp, dp), &
            cmplx(2.0_dp, 0.0_dp, dp)]
        call csc_from_triplet(3, 3, rows, cols, vals, A, status)
        call check_true("cplx_build", status_ok(status), nfail)

        solver%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC
        xe = [cmplx(1.0_dp, -1.0_dp, dp), cmplx(0.0_dp, 2.0_dp, dp), &
            cmplx(3.0_dp, 0.0_dp, dp)]
        b = csc_matvec(A, xe)
        call sparse_factor(solver, A, status)
        call check_true("cplx_factor", status_ok(status), nfail)

        call sparse_solve(solver, b, x, status)
        call check_true("cplx_solve", status_ok(status), nfail)
        call check_zerr("cplx_x", x, xe, nfail)
        call sparse_free(solver)
    end subroutine run_complex

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

    subroutine check_zerr(label, got, want, nfail)
        character(*), intent(in)    :: label
        complex(dp),  intent(in)    :: got(:), want(:)
        integer,      intent(inout) :: nfail
        if (maxval(abs(got - want)) > 1.0e-9_dp) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,es13.6)") "FAIL [", label, &
                "] max|err| ", maxval(abs(got - want))
        end if
    end subroutine check_zerr

    subroutine check_true(label, cond, nfail)
        character(*), intent(in)    :: label
        logical,      intent(in)    :: cond
        integer,      intent(inout) :: nfail
        if (.not. cond) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a)") "FAIL [", label, "] expected .true."
        end if
    end subroutine check_true

end program test_fortsparse_umfpack_ipc
