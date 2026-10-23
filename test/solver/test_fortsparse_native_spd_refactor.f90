program test_fortsparse_native_spd_refactor
    ! Same-pattern refactor of the native SPD backend: the pattern checks of the
    ! analysis are not repeated, but non-finite and asymmetric values are still
    ! rejected, a non-positive-definite matrix keeps the analysis (the next
    ! refactor reuses it) and a changed pattern is analysed again.
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortsparse
    implicit none
    integer, parameter :: n = 24
    integer :: ordering
    do ordering = FORTSPARSE_ORDER_DEFAULT, FORTSPARSE_ORDER_ND
        call refactor_checks(ordering)
    end do
    print '(a)', 'Native SPD same-pattern refactor checks pass'
contains
    subroutine check(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) error stop message
    end subroutine check

    subroutine laplacian(shift, a)
        ! Pentadiagonal SPD for shift > 0, indefinite for a negative shift.
        real(dp), intent(in) :: shift
        type(csc_t), intent(out) :: a
        type(fortsparse_status_t) :: status
        real(dp) :: dense(n, n)
        integer :: i, j, count, rows(5*n), cols(5*n)
        real(dp) :: vals(5*n)
        dense = 0
        do i = 1, n
            dense(i, i) = 4 + shift
            if (i < n) dense(i, i + 1) = -1
            if (i < n) dense(i + 1, i) = -1
            if (i < n - 2) dense(i, i + 3) = -0.5_dp
            if (i < n - 2) dense(i + 3, i) = -0.5_dp
        end do
        count = 0
        do j = 1, n
            do i = 1, n
                if (dense(i, j) == 0) cycle
                count = count + 1
                rows(count) = i; cols(count) = j; vals(count) = dense(i, j)
            end do
        end do
        call csc_from_triplet(n, n, rows(:count), cols(:count), vals(:count), a, status)
        call check(status_ok(status), 'fixture')
    end subroutine

    subroutine refactor_checks(ordering)
        integer, intent(in) :: ordering
        type(sparse_solver_t) :: solver
        type(csc_t) :: a, b, bad
        type(fortsparse_status_t) :: status
        real(dp) :: x(n), rhs(n), exact(n)
        integer :: i
        logical :: reused
        do i = 1, n
            exact(i) = 1 + 0.1_dp*i
        end do
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry = FORTSPARSE_SYM_SPD
        solver%options%ordering = ordering
        call laplacian(0.0_dp, a)
        call sparse_factor(solver, a, status)
        call check(status_ok(status), 'factor')
        ! Same pattern, new values: reused, and the solve is exact.
        call laplacian(1.0_dp, b)
        call sparse_refactor(solver, b, status, reused)
        call check(status_ok(status) .and. reused, 'refactor reuses the analysis')
        call matvec(b, exact, rhs)
        call sparse_solve(solver, rhs, x, status)
        call check(status_ok(status), 'solve')
        call check(maxval(abs(x - exact)) < 1.0e-12_dp, 'refactor solution')
        ! Non-finite and asymmetric values are rejected on the refactor path.
        bad = b
        bad%val(3) = ieee_value(1.0_dp, ieee_quiet_nan)
        call sparse_refactor(solver, bad, status, reused)
        call check(.not. status_ok(status), 'NaN value rejected')
        bad = b
        bad%val(3) = bad%val(3)*(1 + 1.0e-12_dp)
        call sparse_refactor(solver, bad, status, reused)
        call check(.not. status_ok(status), 'asymmetric value rejected')
        ! Not positive definite: fails, and keeps the analysis.
        call laplacian(-3.0_dp, bad)
        call sparse_refactor(solver, bad, status, reused)
        call check(status%code == FORTSPARSE_NOT_POSITIVE_DEFINITE, 'indefinite matrix')
        call sparse_refactor(solver, b, status, reused)
        call check(status_ok(status) .and. reused, 'analysis kept after indefinite')
        call sparse_solve(solver, rhs, x, status)
        call check(status_ok(status) .and. maxval(abs(x - exact)) < 1.0e-12_dp, &
            'solution after indefinite')
        ! A changed pattern is analysed again and validated in full.
        bad = b
        bad%row_idx(2) = bad%row_idx(2) + 1
        call sparse_refactor(solver, bad, status, reused)
        call check(.not. status_ok(status), 'unsorted pattern rejected')
        call sparse_refactor(solver, a, status, reused)
        call check(status_ok(status), 'refactor after rejection')
    end subroutine

    subroutine matvec(a, v, y)
        type(csc_t), intent(in) :: a
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: y(:)
        integer :: j, p
        y = 0
        do j = 1, a%ncol
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                y(a%row_idx(p)) = y(a%row_idx(p)) + a%val(p)*v(j)
            end do
        end do
    end subroutine
end program
