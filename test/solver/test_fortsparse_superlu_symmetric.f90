program test_fortsparse_superlu_symmetric
    ! SuperLU with a symmetric hint (INDEFINITE, SPD, HERMITIAN) uses its
    ! symmetric-mode ordering and threshold pivoting. Saddle-point check: the
    ! 5-point Laplacian on an m x m grid, shifted below its (2,2) eigenvalue
    ! (indefinite), bordered by a mean-value Lagrange row and column with a
    ! zero diagonal, so a diagonal pivot is not always admissible. The
    ! manufactured solution must be recovered with every hint, real and
    ! complex (complex symmetric for INDEFINITE, Hermitian for HERMITIAN).
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: dp, csc_t, csc_z_t, csc_from_triplet, csc_matvec, &
        sparse_solver_t, sparse_factor, sparse_solve, sparse_free, &
        fortsparse_status_t, status_ok, FORTSPARSE_SYM_GENERAL, &
        FORTSPARSE_SYM_SPD, FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_SYM_HERMITIAN
    implicit none

    integer, parameter :: m = 12, n = m*m + 1
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: nfail
    nfail = 0

    call run_real(FORTSPARSE_SYM_GENERAL, 'real_general', nfail)
    call run_real(FORTSPARSE_SYM_INDEFINITE, 'real_indefinite', nfail)
    call run_spd(nfail)
    call run_complex(FORTSPARSE_SYM_GENERAL, 'complex_general', nfail)
    call run_complex(FORTSPARSE_SYM_INDEFINITE, 'complex_symmetric', nfail)
    call run_complex(FORTSPARSE_SYM_HERMITIAN, 'complex_hermitian', nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    subroutine run_real(symmetry, label, nfail)
        integer, intent(in) :: symmetry
        character(*), intent(in) :: label
        integer, intent(inout) :: nfail
        type(csc_t) :: A
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        real(dp) :: xe(n), b(n), x(n)
        integer :: i
        call saddle(rows, cols, vals, .true.)
        call csc_from_triplet(n, n, rows, cols, vals, A, status)
        call check(label//'_build', status_ok(status), nfail)
        xe = [(sin(0.37_dp*i) + 0.5_dp, i = 1, n)]
        b = csc_matvec(A, xe)
        solver%options%symmetry = symmetry
        call sparse_factor(solver, A, status)
        call check(label//'_factor', status_ok(status), nfail)
        call sparse_solve(solver, b, x, status)
        call check(label//'_solve', status_ok(status), nfail)
        call check(label//'_error', norm2(x - xe) <= 1.0e-10_dp*norm2(xe), nfail)
        call sparse_free(solver)
    end subroutine run_real

    subroutine run_spd(nfail)
        ! The unshifted Laplacian block alone (SPD).
        integer, intent(inout) :: nfail
        type(csc_t) :: A
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        real(dp) :: xe(n - 1), b(n - 1), x(n - 1)
        integer :: i
        call saddle(rows, cols, vals, .false.)
        call csc_from_triplet(n - 1, n - 1, rows, cols, vals, A, status)
        call check('spd_build', status_ok(status), nfail)
        xe = [(cos(0.21_dp*i), i = 1, n - 1)]
        b = csc_matvec(A, xe)
        solver%options%symmetry = FORTSPARSE_SYM_SPD
        call sparse_factor(solver, A, status)
        call check('spd_factor', status_ok(status), nfail)
        call sparse_solve(solver, b, x, status)
        call check('spd_solve', status_ok(status), nfail)
        call check('spd_error', norm2(x - xe) <= 1.0e-11_dp*norm2(xe), nfail)
        call sparse_free(solver)
    end subroutine run_spd

    subroutine run_complex(symmetry, label, nfail)
        ! Real saddle matrix plus i*0.3 on the Laplacian diagonal (complex
        ! symmetric) or a Hermitian skew coupling i*0.2 on the x-neighbours.
        integer, intent(in) :: symmetry
        character(*), intent(in) :: label
        integer, intent(inout) :: nfail
        type(csc_z_t) :: A
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        complex(dp), allocatable :: zvals(:)
        complex(dp) :: xe(n), b(n), x(n)
        integer :: i
        call saddle(rows, cols, vals, .true.)
        zvals = cmplx(vals, 0.0_dp, dp)
        do i = 1, size(rows)
            if (rows(i) == cols(i) .and. rows(i) < n) then
                if (symmetry /= FORTSPARSE_SYM_HERMITIAN) &
                    zvals(i) = zvals(i) + cmplx(0.0_dp, 0.3_dp, dp)
            else if (rows(i) == cols(i) + 1 .and. rows(i) < n) then
                if (symmetry == FORTSPARSE_SYM_HERMITIAN) &
                    zvals(i) = zvals(i) + cmplx(0.0_dp, 0.2_dp, dp)
            else if (cols(i) == rows(i) + 1 .and. cols(i) < n) then
                if (symmetry == FORTSPARSE_SYM_HERMITIAN) &
                    zvals(i) = zvals(i) - cmplx(0.0_dp, 0.2_dp, dp)
            end if
        end do
        call csc_from_triplet(n, n, rows, cols, zvals, A, status)
        call check(label//'_build', status_ok(status), nfail)
        xe = [(cmplx(sin(0.37_dp*i), cos(0.11_dp*i), dp), i = 1, n)]
        b = csc_matvec(A, xe)
        solver%options%symmetry = symmetry
        call sparse_factor(solver, A, status)
        call check(label//'_factor', status_ok(status), nfail)
        call sparse_solve(solver, b, x, status)
        call check(label//'_solve', status_ok(status), nfail)
        call check(label//'_error', norm2(abs(x - xe)) <= 1.0e-10_dp*norm2(abs(xe)), &
            nfail)
        call sparse_free(solver)
    end subroutine run_complex

    subroutine saddle(rows, cols, vals, bordered)
        ! Triplets of the 5-point Laplacian (h = 1) on the m x m grid; with
        ! bordered, shifted between its (1,2) and (2,2) eigenvalues and
        ! bordered by the all-ones Lagrange row/column (zero corner).
        integer, allocatable, intent(out) :: rows(:), cols(:)
        real(dp), allocatable, intent(out) :: vals(:)
        logical, intent(in) :: bordered
        real(dp) :: shift, lam(3)
        integer :: i, j, p, k
        allocate (rows(7*m*m + 1), cols(7*m*m + 1), vals(7*m*m + 1))
        shift = 0
        if (bordered) then
            ! Eigenvalues 4 - 2 cos(a pi/(m+1)) - 2 cos(b pi/(m+1)), (a, b) =
            ! (1, 1), (1, 2), (2, 2).
            lam(1) = 4 - 4*cos(pi/(m + 1))
            lam(2) = 4 - 2*cos(pi/(m + 1)) - 2*cos(2*pi/(m + 1))
            lam(3) = 4 - 4*cos(2*pi/(m + 1))
            shift = (lam(2) + lam(3))/2
        end if
        k = 0
        do j = 1, m
            do i = 1, m
                p = i + m*(j - 1)
                call add(rows, cols, vals, k, p, p, 4 - shift)
                if (i > 1) call add(rows, cols, vals, k, p, p - 1, -1.0_dp)
                if (i < m) call add(rows, cols, vals, k, p, p + 1, -1.0_dp)
                if (j > 1) call add(rows, cols, vals, k, p, p - m, -1.0_dp)
                if (j < m) call add(rows, cols, vals, k, p, p + m, -1.0_dp)
                if (bordered) then
                    call add(rows, cols, vals, k, p, n, 1.0_dp)
                    call add(rows, cols, vals, k, n, p, 1.0_dp)
                end if
            end do
        end do
        if (bordered) call add(rows, cols, vals, k, n, n, 0.0_dp)
        rows = rows(:k); cols = cols(:k); vals = vals(:k)
    end subroutine saddle

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

end program test_fortsparse_superlu_symmetric
