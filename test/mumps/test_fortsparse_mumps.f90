program test_fortsparse_mumps
    ! MUMPS backend: SPD residuals on a 2D Laplacian and a P3 Poisson matrix,
    ! factor reuse over many right-hand sides, block solve equal to repeated
    ! single solves, thread invariance, and the symmetric-indefinite and
    ! general paths.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: dp, i8, csc_t, csc_matvec, sparse_solver_t, &
        sparse_factor, sparse_solve, sparse_solve_block, &
        sparse_destroy, sparse_factor_stats, sparse_factor_stats_t, &
        fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_MUMPS, &
        FORTSPARSE_SYM_SPD, FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_SYM_GENERAL, &
        FORTSPARSE_SINGULAR
    use fs_test_matrices, only: laplace2d, p3_poisson, inf_residual
    implicit none

    integer :: nfail
    real(dp), parameter :: pi = acos(-1.0_dp)

    nfail = 0
    call test_laplacian(nfail)
    call test_p3_energy(nfail)
    call test_p3_reuse(nfail)
    call test_block_equals_single(nfail)
    call test_thread_invariance(nfail)
    call test_indefinite_and_general(nfail)
    call test_not_positive_definite(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"

contains

    subroutine check_true(name, ok, nfail)
        character(*), intent(in)    :: name
        logical,      intent(in)    :: ok
        integer,      intent(inout) :: nfail

        if (.not. ok) then
            write (error_unit, "(a,a)") "FAIL: ", name
            nfail = nfail + 1
        end if
    end subroutine check_true

    function pattern(n, k) result(x)
        integer, intent(in) :: n, k
        real(dp) :: x(n)
        integer :: i

        do i = 1, n
            x(i) = sin(1.234_dp*i + 0.7_dp*k) + 0.5_dp*cos(0.37_dp*i*k)
        end do
    end function pattern

    subroutine spd_solver(solver, nthreads)
        type(sparse_solver_t), intent(out) :: solver
        integer, intent(in) :: nthreads

        solver%backend_id = FORTSPARSE_BACKEND_MUMPS
        solver%refine = .true.
        solver%options%symmetry = FORTSPARSE_SYM_SPD
        solver%options%nthreads = nthreads
    end subroutine spd_solver

    subroutine test_laplacian(nfail)
        integer, intent(inout) :: nfail
        type(csc_t) :: A
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        real(dp), allocatable :: b(:), x(:), xt(:)
        real(dp) :: res, err

        call laplace2d(150, A)
        call spd_solver(solver, 2)
        call sparse_factor(solver, A, status)
        call check_true("laplace factor", status_ok(status), nfail)
        if (.not. status_ok(status)) return
        xt = pattern(A%ncol, 1)
        b = csc_matvec(A, xt)
        allocate (x(A%ncol))
        call sparse_solve(solver, b, x, status)
        call check_true("laplace solve", status_ok(status), nfail)
        res = inf_residual(A, x, b)
        err = maxval(abs(x - xt))/maxval(abs(xt))
        write (*, "(a,es10.2,a,es10.2)") "laplace2d 150^2: residual", res, &
            " error", err
        call check_true("laplace residual", res <= 1.0e-14_dp, nfail)
        call check_true("laplace error", err <= 1.0e-9_dp, nfail)
        call sparse_destroy(solver)
    end subroutine test_laplacian

    ! The P3 interpolant of sin(pi x) sin(pi y) must have stiffness energy
    ! pi^2/2 up to the interpolation error, which validates the global dof
    ! numbering and element matrices independently of the solver.
    subroutine test_p3_energy(nfail)
        integer, intent(inout) :: nfail
        type(csc_t) :: A
        real(dp), allocatable :: xy(:, :), u(:), au(:)
        real(dp) :: energy, exact
        integer :: k

        call p3_poisson(16, A, xy)
        allocate (u(A%ncol))
        do k = 1, A%ncol
            if (xy(1, k) <= 1.0e-12_dp .or. xy(1, k) >= 1.0_dp - 1.0e-12_dp .or. &
                xy(2, k) <= 1.0e-12_dp .or. xy(2, k) >= 1.0_dp - 1.0e-12_dp) then
                u(k) = 0.0_dp
            else
                u(k) = sin(pi*xy(1, k))*sin(pi*xy(2, k))
            end if
        end do
        au = csc_matvec(A, u)
        energy = dot_product(u, au)
        exact = 0.5_dp*pi**2
        write (*, "(a,es12.5,a,es12.5)") "p3 energy", energy, " exact", exact
        call check_true("p3 energy", abs(energy - exact) <= 1.0e-5_dp*exact, nfail)
    end subroutine test_p3_energy

    subroutine test_p3_reuse(nfail)
        integer, intent(inout) :: nfail
        integer, parameter :: nrhs = 25
        type(csc_t) :: A
        type(sparse_solver_t) :: solver
        type(sparse_factor_stats_t) :: stats0, stats1
        type(fortsparse_status_t) :: status
        real(dp), allocatable :: b(:), x(:)
        real(dp) :: worst, res
        integer :: k

        call p3_poisson(40, A)
        call spd_solver(solver, 4)
        call sparse_factor(solver, A, status)
        call check_true("p3 factor", status_ok(status), nfail)
        if (.not. status_ok(status)) return
        call sparse_factor_stats(solver, stats0, status)
        call check_true("p3 stats", status_ok(status) .and. &
            stats0%factor_nnz > int(A%nnz, i8)/2, nfail)
        allocate (x(A%ncol))
        worst = 0.0_dp
        do k = 1, nrhs
            b = pattern(A%ncol, k)
            call sparse_solve(solver, b, x, status)
            call check_true("p3 solve", status_ok(status), nfail)
            res = inf_residual(A, x, b)
            worst = max(worst, res)
        end do
        call sparse_factor_stats(solver, stats1, status)
        write (*, "(a,i0,a,es10.2,a,i0,a,f8.1,a)") "p3 n=", A%ncol, &
            " worst residual over 25 RHS", worst, " factor nnz", &
            stats1%factor_nnz, " mem", stats1%factor_mem_mb, " MB"
        call check_true("p3 residual", worst <= 1.0e-14_dp, nfail)
        call check_true("p3 factor reused", stats1%factor_nnz == stats0%factor_nnz, &
            nfail)
        call sparse_destroy(solver)
    end subroutine test_p3_reuse

    subroutine test_block_equals_single(nfail)
        integer, intent(inout) :: nfail
        integer, parameter :: nrhs = 9
        type(csc_t) :: A
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        real(dp), allocatable :: B(:, :), X(:, :), x1(:)
        real(dp) :: diff
        integer :: k

        call p3_poisson(30, A)
        call spd_solver(solver, 4)
        call sparse_factor(solver, A, status)
        call check_true("block factor", status_ok(status), nfail)
        if (.not. status_ok(status)) return
        allocate (B(A%ncol, nrhs), X(A%ncol, nrhs), x1(A%ncol))
        do k = 1, nrhs
            B(:, k) = pattern(A%ncol, k)
        end do
        call sparse_solve_block(solver, B, X, status)
        call check_true("block solve", status_ok(status), nfail)
        diff = 0.0_dp
        do k = 1, nrhs
            call sparse_solve(solver, B(:, k), x1, status)
            diff = max(diff, maxval(abs(X(:, k) - x1))/maxval(abs(x1)))
        end do
        write (*, "(a,es10.2)") "block vs single max rel diff", diff
        call check_true("block equals single", diff <= 1.0e-13_dp, nfail)
        call sparse_destroy(solver)
    end subroutine test_block_equals_single

    subroutine test_thread_invariance(nfail)
        integer, intent(inout) :: nfail
        integer, parameter :: counts(3) = [1, 2, 8]
        type(csc_t) :: A
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        real(dp), allocatable :: b(:), x(:, :)
        integer :: i
        real(dp) :: d

        call p3_poisson(60, A)
        b = pattern(A%ncol, 3)
        allocate (x(A%ncol, size(counts)))
        do i = 1, size(counts)
            call spd_solver(solver, counts(i))
            call sparse_factor(solver, A, status)
            call check_true("thread factor", status_ok(status), nfail)
            if (.not. status_ok(status)) return
            call sparse_solve(solver, b, x(:, i), status)
            call check_true("thread solve", status_ok(status), nfail)
            call sparse_destroy(solver)
        end do
        do i = 2, size(counts)
            d = norm2(x(:, i) - x(:, 1))/norm2(x(:, 1))
            write (*, "(a,i0,a,es10.2)") "threads ", counts(i), &
                " vs 1: relative difference", d
            call check_true("thread invariance", d <= 1.0e-13_dp, nfail)
        end do
    end subroutine test_thread_invariance

    ! Symmetric indefinite (shifted Laplacian, LDL^T) and a nonsymmetric
    ! (scaled subdiagonal, LU) system.
    subroutine test_indefinite_and_general(nfail)
        integer, intent(inout) :: nfail
        type(csc_t) :: A, G
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        real(dp), allocatable :: b(:), x(:)
        integer :: j, e

        call laplace2d(60, A)
        G = A
        do j = 1, A%ncol
            do e = A%col_ptr(j), A%col_ptr(j + 1) - 1
                if (A%row_idx(e) == j) A%val(e) = A%val(e) - 4.05_dp
                if (G%row_idx(e) == j + 1) G%val(e) = 1.3_dp*G%val(e)
            end do
        end do
        solver%backend_id = FORTSPARSE_BACKEND_MUMPS
        solver%refine = .true.
        solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        call sparse_factor(solver, A, status)
        call check_true("indefinite factor", status_ok(status), nfail)
        b = pattern(A%ncol, 2)
        allocate (x(A%ncol))
        if (status_ok(status)) then
            call sparse_solve(solver, b, x, status)
            call check_true("indefinite residual", &
                status_ok(status) .and. inf_residual(A, x, b) <= 1.0e-13_dp, &
                nfail)
        end if
        call sparse_destroy(solver)

        solver%backend_id = FORTSPARSE_BACKEND_MUMPS
        solver%refine = .true.
        solver%options%symmetry = FORTSPARSE_SYM_GENERAL
        call sparse_factor(solver, G, status)
        call check_true("general factor", status_ok(status), nfail)
        if (status_ok(status)) then
            call sparse_solve(solver, b, x, status)
            call check_true("general residual", &
                status_ok(status) .and. inf_residual(G, x, b) <= 1.0e-13_dp, &
                nfail)
        end if
        call sparse_destroy(solver)
    end subroutine test_indefinite_and_general

    ! An indefinite matrix declared SPD must fail with a singular status
    ! instead of returning a wrong solution.
    subroutine test_not_positive_definite(nfail)
        integer, intent(inout) :: nfail
        type(csc_t) :: A
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        integer :: j, e

        call laplace2d(20, A)
        do j = 1, A%ncol
            do e = A%col_ptr(j), A%col_ptr(j + 1) - 1
                if (A%row_idx(e) == j) A%val(e) = A%val(e) - 4.05_dp
            end do
        end do
        call spd_solver(solver, 1)
        call sparse_factor(solver, A, status)
        call check_true("indefinite declared spd fails", &
            status%code == FORTSPARSE_SINGULAR, nfail)
        call sparse_destroy(solver)
    end subroutine test_not_positive_definite

end program test_fortsparse_mumps
