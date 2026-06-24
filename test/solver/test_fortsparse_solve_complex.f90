program test_fortsparse_solve_complex
    ! Complex solve on the default (SuperLU) backend: factor a small complex
    ! system, solve, and check the solution and residual.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: dp, csc_z_t, csc_from_triplet, csc_matvec, &
        sparse_solver_t, sparse_factor, sparse_solve, sparse_free, &
        sparse_solve_once, fortsparse_status_t, status_ok
    implicit none

    integer :: nfail
    nfail = 0

    call run(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    subroutine run(nfail)
        integer, intent(inout) :: nfail

        type(csc_z_t)             :: A
        type(sparse_solver_t)     :: solver
        type(fortsparse_status_t) :: status
        integer                   :: rows(4), cols(4)
        complex(dp)               :: vals(4)
        complex(dp)               :: b(2), x(2), xe(2)
        complex(dp), parameter    :: ii = (0.0_dp, 1.0_dp)

        ! A = [[2, 1-i],[1+i, 3]]
        rows = [1, 1, 2, 2]
        cols = [1, 2, 1, 2]
        vals = [(2.0_dp, 0.0_dp), (1.0_dp, -1.0_dp), &
            (1.0_dp, 1.0_dp), (3.0_dp, 0.0_dp)]
        call csc_from_triplet(2, 2, rows, cols, vals, A, status)
        call check_true("build_ok", status_ok(status), nfail)

        xe = [(1.0_dp, 0.0_dp), ii]
        b = [(3.0_dp, 1.0_dp), (1.0_dp, 4.0_dp)]

        call sparse_factor(solver, A, status)
        call check_true("factor_ok", status_ok(status), nfail)
        call sparse_solve(solver, b, x, status)
        call check_true("solve_ok", status_ok(status), nfail)
        call check_err("solve_x", x, xe, nfail)
        call check_residual("solve_res", A, x, b, nfail)
        call sparse_free(solver)

        call sparse_solve_once(A, b, x, status)
        call check_true("once_ok", status_ok(status), nfail)
        call check_err("once_x", x, xe, nfail)
    end subroutine run

    subroutine check_err(label, got, want, nfail)
        character(*), intent(in)    :: label
        complex(dp),  intent(in)    :: got(:), want(:)
        integer,      intent(inout) :: nfail
        if (maxval(abs(got - want)) > 1.0e-9_dp) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,es13.6)") "FAIL [", label, &
                "] max|err| ", maxval(abs(got - want))
        end if
    end subroutine check_err

    subroutine check_residual(label, A, x, b, nfail)
        character(*),  intent(in)    :: label
        type(csc_z_t), intent(in)    :: A
        complex(dp),   intent(in)    :: x(:), b(:)
        integer,       intent(inout) :: nfail
        real(dp)                     :: res

        res = maxval(abs(csc_matvec(A, x) - b))
        if (res > 1.0e-9_dp) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,es13.6)") "FAIL [", label, &
                "] residual ", res
        end if
    end subroutine check_residual

    subroutine check_true(label, cond, nfail)
        character(*), intent(in)    :: label
        logical,      intent(in)    :: cond
        integer,      intent(inout) :: nfail
        if (.not. cond) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a)") "FAIL [", label, "] expected .true."
        end if
    end subroutine check_true

end program test_fortsparse_solve_complex
