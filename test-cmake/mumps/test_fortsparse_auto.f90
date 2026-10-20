program test_fortsparse_auto
    ! FORTSPARSE_BACKEND_AUTO is the default backend. Below the measured
    ! crossover it factors with SuperLU, above it with MUMPS; on both sides,
    ! real and complex, the solution is correct, a same-pattern refactor keeps
    ! the chosen backend (and the MUMPS analysis), and sparse_solve_once works.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: dp, csc_t, csc_z_t, sparse_solver_t, sparse_factor, &
        sparse_refactor, sparse_solve, sparse_destroy, sparse_solve_once, &
        fortsparse_status_t, status_ok, sparse_auto_backend, &
        FORTSPARSE_BACKEND_AUTO, FORTSPARSE_BACKEND_DEFAULT, &
        FORTSPARSE_BACKEND_MUMPS, FORTSPARSE_BACKEND_SUPERLU
    use fortsparse_mumps, only: mumps_backend_t
    use fortsparse_superlu, only: superlu_backend_t
    use fs_test_mpi, only: test_mpi_init, test_mpi_finalize, test_mpi_parallel
    use fs_test_matrices, only: laplace2d
    implicit none

    type(csc_t) :: A
    type(csc_z_t) :: Az
    integer :: nfail, m, kind

    nfail = 0
    call check("default_is_auto", FORTSPARSE_BACKEND_DEFAULT == &
        FORTSPARSE_BACKEND_AUTO)

    ! Small side: a 4 x 4 grid.
    call laplace2d(4, A)
    call check("small_rule_superlu", sparse_auto_backend(A%nnz, &
        .false.) == FORTSPARSE_BACKEND_SUPERLU)
    call check_real(A, FORTSPARSE_BACKEND_SUPERLU, "small_real")
    call complexify(A, Az)
    call check_complex(Az, FORTSPARSE_BACKEND_SUPERLU, "small_complex")

    ! A bounded grid above the measured crossover. MPI AUTO must recognize
    ! the application's runtime instead of silently keeping SuperLU.
    m = 96
    if (test_mpi_parallel) call check("before_init_rule", &
        sparse_auto_backend(5*m*m - 4*m, .false.) == FORTSPARSE_BACKEND_SUPERLU)
    call test_mpi_init()
    kind = sparse_auto_backend(5*m*m - 4*m, .false.)
    call check("large_rule_mumps", kind == FORTSPARSE_BACKEND_MUMPS)
    call laplace2d(m, A)
    call check_real(A, FORTSPARSE_BACKEND_MUMPS, "large_real")
    call complexify(A, Az)
    call check("large_complex_rule_mumps", sparse_auto_backend( &
        Az%nnz, .true.) == FORTSPARSE_BACKEND_MUMPS)
    call check_complex(Az, FORTSPARSE_BACKEND_MUMPS, "large_complex")
    call test_mpi_finalize()
    if (test_mpi_parallel) call check("after_finalize_rule", &
        sparse_auto_backend(A%nnz, .false.) == FORTSPARSE_BACKEND_SUPERLU)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"

contains

    subroutine check_real(M, expect, tag)
        type(csc_t), intent(in) :: M
        integer, intent(in) :: expect
        character(*), intent(in) :: tag
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(csc_t) :: M2
        real(dp), allocatable :: xe(:), b(:), x(:)
        logical :: reused
        integer :: i

        xe = [(sin(0.3_dp*i), i=1, M%ncol)]
        b = matvec_r(M, xe)
        allocate (x(M%ncol))
        call check(tag//"_solver_default", solver%backend_id == &
            FORTSPARSE_BACKEND_AUTO)
        call sparse_factor(solver, M, status)
        if (.not. status_ok(status)) write (error_unit, '(a)') status%msg
        call check(tag//"_factor", status_ok(status))
        call check(tag//"_backend", backend_of(solver) == expect)
        call sparse_solve(solver, b, x, status)
        call check(tag//"_solve", status_ok(status))
        call check(tag//"_solution", maxval(abs(x - xe)) < 1.0e-9_dp)
        M2 = M
        M2%val = 2.0_dp*M%val
        call sparse_refactor(solver, M2, status, reused)
        call check(tag//"_refactor", status_ok(status))
        call check(tag//"_refactor_backend", backend_of(solver) == expect)
        if (expect == FORTSPARSE_BACKEND_MUMPS) &
            call check(tag//"_refactor_reused", reused)
        call sparse_solve(solver, 2.0_dp*b, x, status)
        call check(tag//"_refactor_solution", status_ok(status) .and. &
            maxval(abs(x - xe)) < 1.0e-9_dp)
        call sparse_destroy(solver)
        x = 0.0_dp
        call sparse_solve_once(M, b, x, status)
        call check(tag//"_once", status_ok(status) .and. &
            maxval(abs(x - xe)) < 1.0e-9_dp)
    end subroutine check_real

    subroutine check_complex(M, expect, tag)
        type(csc_z_t), intent(in) :: M
        integer, intent(in) :: expect
        character(*), intent(in) :: tag
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        complex(dp), allocatable :: xe(:), b(:), x(:)
        integer :: i

        xe = [(cmplx(sin(0.3_dp*i), cos(0.2_dp*i), dp), i=1, M%ncol)]
        b = matvec_z(M, xe)
        allocate (x(M%ncol))
        call sparse_factor(solver, M, status)
        call check(tag//"_factor", status_ok(status))
        call check(tag//"_backend", backend_of(solver) == expect)
        call sparse_solve(solver, b, x, status)
        call check(tag//"_solution", status_ok(status) .and. &
            maxval(abs(x - xe)) < 1.0e-9_dp)
        call sparse_destroy(solver)
        x = 0.0_dp
        call sparse_solve_once(M, b, x, status)
        call check(tag//"_once", status_ok(status) .and. &
            maxval(abs(x - xe)) < 1.0e-9_dp)
    end subroutine check_complex

    ! A + 0.5 i diag: nonsingular complex symmetric.
    subroutine complexify(M, Z)
        type(csc_t), intent(in) :: M
        type(csc_z_t), intent(out) :: Z
        integer :: j, k

        Z%nrow = M%nrow
        Z%ncol = M%ncol
        Z%nnz = M%nnz
        Z%col_ptr = M%col_ptr
        Z%row_idx = M%row_idx
        Z%val = cmplx(M%val, 0.0_dp, dp)
        do j = 1, M%ncol
            do k = M%col_ptr(j), M%col_ptr(j + 1) - 1
                if (M%row_idx(k) == j) Z%val(k) = Z%val(k) + (0.0_dp, 0.5_dp)
            end do
        end do
    end subroutine complexify

    integer function backend_of(s)
        type(sparse_solver_t), intent(in) :: s

        backend_of = -1
        if (.not. allocated(s%backend)) return
        select type (backend => s%backend)
        type is (mumps_backend_t)
            backend_of = FORTSPARSE_BACKEND_MUMPS
        type is (superlu_backend_t)
            backend_of = FORTSPARSE_BACKEND_SUPERLU
        end select
    end function backend_of

    function matvec_r(M, v) result(y)
        type(csc_t), intent(in) :: M
        real(dp), intent(in) :: v(:)
        real(dp) :: y(M%nrow)
        integer :: j, k

        y = 0.0_dp
        do j = 1, M%ncol
            do k = M%col_ptr(j), M%col_ptr(j + 1) - 1
                y(M%row_idx(k)) = y(M%row_idx(k)) + M%val(k)*v(j)
            end do
        end do
    end function matvec_r

    function matvec_z(M, v) result(y)
        type(csc_z_t), intent(in) :: M
        complex(dp), intent(in) :: v(:)
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
        logical, intent(in) :: ok

        if (.not. ok) then
            write (error_unit, "(a,a)") "FAIL: ", name
            nfail = nfail + 1
        end if
    end subroutine check
end program test_fortsparse_auto
