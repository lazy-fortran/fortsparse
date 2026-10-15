program test_fortsparse_superlu_refactor
    ! SuperLU refactorization keeps the column ordering of the retained real
    ! factorization when the pattern and symmetry hint are unchanged. On new
    ! values it must report reuse and give exactly the solution of a fresh
    ! factorization of the same matrix; a changed pattern or hint falls back to
    ! a full factorization (no reuse) that still solves the new system.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: dp, csc_t, csc_from_triplet, csc_matvec, &
        sparse_solver_t, sparse_factor, sparse_refactor, sparse_solve, &
        sparse_free, fortsparse_status_t, status_ok, FORTSPARSE_SYM_GENERAL, &
        FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_BACKEND_SUPERLU
    implicit none

    integer, parameter :: m = 10, n = m*m
    integer :: nfail, i
    type(csc_t) :: A1, A2, A3
    type(sparse_solver_t) :: solver, fresh
    type(fortsparse_status_t) :: status
    real(dp) :: xe(n), b(n), x(n), xf(n)
    logical :: reused

    nfail = 0
    xe = [(sin(0.37_dp*i) + 0.5_dp, i = 1, n)]
    call laplacian(1.3_dp, 0.0_dp, .false., A1)
    call laplacian(2.7_dp, 0.4_dp, .false., A2)
    call laplacian(2.7_dp, 0.4_dp, .true., A3)

    solver%backend_id = FORTSPARSE_BACKEND_SUPERLU
    solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
    call sparse_factor(solver, A1, status)
    call check('factor', status_ok(status), nfail)

    ! Same pattern, new values: reused and identical to a fresh factor.
    call sparse_refactor(solver, A2, status, reused)
    call check('refactor_status', status_ok(status), nfail)
    call check('refactor_reused', reused, nfail)
    b = csc_matvec(A2, xe)
    call sparse_solve(solver, b, x, status)
    call check('refactor_solve', status_ok(status), nfail)
    call check('refactor_error', norm2(x - xe) <= 1.0e-10_dp*norm2(xe), nfail)
    fresh%backend_id = FORTSPARSE_BACKEND_SUPERLU
    fresh%options%symmetry = FORTSPARSE_SYM_INDEFINITE
    call sparse_factor(fresh, A2, status)
    call sparse_solve(fresh, b, xf, status)
    call check('refactor_equals_fresh', all(x == xf), nfail)
    call sparse_free(fresh)

    ! Changed pattern: full factor.
    call sparse_refactor(solver, A3, status, reused)
    call check('pattern_status', status_ok(status), nfail)
    call check('pattern_not_reused', .not. reused, nfail)
    b = csc_matvec(A3, xe)
    call sparse_solve(solver, b, x, status)
    call check('pattern_error', norm2(x - xe) <= 1.0e-10_dp*norm2(xe), nfail)

    ! Changed symmetry hint: full factor.
    solver%options%symmetry = FORTSPARSE_SYM_GENERAL
    call sparse_refactor(solver, A3, status, reused)
    call check('hint_status', status_ok(status), nfail)
    call check('hint_not_reused', .not. reused, nfail)
    call sparse_solve(solver, b, x, status)
    call check('hint_error', norm2(x - xe) <= 1.0e-10_dp*norm2(xe), nfail)
    call sparse_free(solver)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"

contains

    subroutine laplacian(shift, skew, diagonal_links, A)
        ! 5-point Laplacian on the m x m grid minus shift (indefinite for
        ! shift > lowest eigenvalue), plus skew on the x-links (unsymmetric
        ! values, symmetric pattern); diagonal_links adds (p, p + m + 1).
        real(dp), intent(in) :: shift, skew
        logical, intent(in) :: diagonal_links
        type(csc_t), intent(out) :: A
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        integer :: i, j, p, k
        allocate (rows(7*n), cols(7*n), vals(7*n))
        k = 0
        do j = 1, m
            do i = 1, m
                p = i + m*(j - 1)
                call add(rows, cols, vals, k, p, p, 4 - shift)
                if (i > 1) call add(rows, cols, vals, k, p, p - 1, -1 - skew)
                if (i < m) call add(rows, cols, vals, k, p, p + 1, -1 + skew)
                if (j > 1) call add(rows, cols, vals, k, p, p - m, -1.0_dp)
                if (j < m) call add(rows, cols, vals, k, p, p + m, -1.0_dp)
                if (diagonal_links .and. i < m .and. j < m) then
                    call add(rows, cols, vals, k, p, p + m + 1, -0.1_dp)
                    call add(rows, cols, vals, k, p + m + 1, p, -0.1_dp)
                end if
            end do
        end do
        call csc_from_triplet(n, n, rows(:k), cols(:k), vals(:k), A, status)
        call check('build', status_ok(status), nfail)
    end subroutine laplacian

    subroutine add(rows, cols, vals, k, r, c, v)
        integer, intent(inout) :: rows(:), cols(:), k
        real(dp), intent(inout) :: vals(:)
        integer, intent(in) :: r, c
        real(dp), intent(in) :: v
        k = k + 1
        rows(k) = r; cols(k) = c; vals(k) = v
    end subroutine add

    subroutine check(label, cond, nfail)
        character(*), intent(in) :: label
        logical, intent(in) :: cond
        integer, intent(inout) :: nfail
        if (.not. cond) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a)") "FAIL [", label, "]"
        end if
    end subroutine check

end program test_fortsparse_superlu_refactor
