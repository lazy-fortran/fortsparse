program test_fortsparse_mumps_default
    ! A sequential MUMPS build with ZMUMPS makes MUMPS the default backend: a
    ! solver that names no backend factors real and complex systems through
    ! mumps_backend_t, and sparse_solve_once solves both correctly.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: dp, csc_t, csc_z_t, sparse_solver_t, sparse_factor, &
        sparse_solve, sparse_destroy, sparse_solve_once, fortsparse_status_t, &
        status_ok, FORTSPARSE_BACKEND_MUMPS, FORTSPARSE_BACKEND_DEFAULT
    use fortsparse_mumps, only: mumps_backend_t
    implicit none

    type(csc_t) :: A
    type(csc_z_t) :: Az
    type(sparse_solver_t) :: solver, zsolver
    type(fortsparse_status_t) :: status
    real(dp) :: b(3), x(3)
    complex(dp) :: bz(3), xz(3), xz_exact(3)
    integer :: nfail

    nfail = 0
    ! [4 1 0; 1 3 1; 0 1 2] with solution (1, 2, 3).
    A%nrow = 3
    A%ncol = 3
    A%nnz = 7
    A%col_ptr = [1, 3, 6, 8]
    A%row_idx = [1, 2, 1, 2, 3, 2, 3]
    A%val = [4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]
    b = [6.0_dp, 10.0_dp, 8.0_dp]
    Az%nrow = 3
    Az%ncol = 3
    Az%nnz = 7
    Az%col_ptr = A%col_ptr
    Az%row_idx = A%row_idx
    Az%val = cmplx(A%val, 0.5_dp, dp)
    xz_exact = [(1.0_dp, -1.0_dp), (2.0_dp, 0.5_dp), (0.0_dp, 3.0_dp)]
    bz = matvec_z(Az, xz_exact)

    call check("default_is_mumps", FORTSPARSE_BACKEND_DEFAULT == &
        FORTSPARSE_BACKEND_MUMPS)
    call check("solver_default_is_mumps", &
        solver%backend_id == FORTSPARSE_BACKEND_MUMPS)
    call sparse_factor(solver, A, status)
    call check("real_factor", status_ok(status))
    call check("real_backend_mumps", is_mumps(solver))
    call sparse_solve(solver, b, x, status)
    call check("real_solve", status_ok(status))
    call check("real_solution", maxval(abs(x - [1.0_dp, 2.0_dp, 3.0_dp])) < 1.0e-12_dp)
    call sparse_destroy(solver)

    call sparse_factor(zsolver, Az, status)
    call check("complex_factor", status_ok(status))
    call check("complex_backend_mumps", is_mumps(zsolver))
    call sparse_destroy(zsolver)

    x = 0.0_dp
    call sparse_solve_once(A, b, x, status)
    call check("real_once", status_ok(status))
    call check("real_once_solution", &
        maxval(abs(x - [1.0_dp, 2.0_dp, 3.0_dp])) < 1.0e-12_dp)
    call sparse_solve_once(Az, bz, xz, status)
    call check("complex_once", status_ok(status))
    call check("complex_once_solution", maxval(abs(xz - xz_exact)) < 1.0e-12_dp)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"

contains

    logical function is_mumps(s)
        type(sparse_solver_t), intent(in) :: s

        is_mumps = .false.
        if (.not. allocated(s%backend)) return
        select type (backend => s%backend)
        type is (mumps_backend_t)
            is_mumps = .true.
        end select
    end function is_mumps

    function matvec_z(M, v) result(y)
        type(csc_z_t), intent(in) :: M
        complex(dp),   intent(in) :: v(:)
        complex(dp) :: y(M%nrow)
        integer :: j, k

        y = (0.0_dp, 0.0_dp)
        do j = 1, M%ncol
            do k = M%col_ptr(j), M%col_ptr(j + 1) - 1
                y(M%row_idx(k)) = y(M%row_idx(k)) + M%val(k)*v(j)
            end do
        end do
    end function matvec_z

    subroutine check(name, ok)
        character(*), intent(in) :: name
        logical,      intent(in) :: ok

        if (.not. ok) then
            write (error_unit, "(a,a)") "FAIL: ", name
            nfail = nfail + 1
        end if
    end subroutine check
end program test_fortsparse_mumps_default
