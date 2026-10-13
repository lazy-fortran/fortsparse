program test_fortsparse_mumps_refactor
    ! sparse_refactor on the real MUMPS backend: new values on an unchanged
    ! pattern reuse the symbolic analysis and match a fresh factorization to
    ! roundoff; a changed pattern, size or symmetry option falls back to a
    ! full analysis; a numerical failure leaves the solver unfactored and the
    ! next refactor recovers.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: dp, csc_t, csc_from_triplet, sparse_solver_t, &
        sparse_factor, sparse_refactor, sparse_solve, sparse_destroy, &
        sparse_factor_stats, sparse_factor_stats_t, &
        fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_MUMPS, &
        FORTSPARSE_SYM_SPD, FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_SYM_GENERAL, &
        FORTSPARSE_SINGULAR
    use fs_test_matrices, only: laplace2d
    use fs_test_mpi, only: test_mpi_init, test_mpi_finalize
    implicit none

    integer :: nfail

    call test_mpi_init()
    nfail = 0
    call test_values(FORTSPARSE_SYM_SPD, .false., "spd", nfail)
    call test_values(FORTSPARSE_SYM_SPD, .true., "spd refined", nfail)
    call test_values(FORTSPARSE_SYM_INDEFINITE, .false., "indefinite", nfail)
    call test_values(FORTSPARSE_SYM_GENERAL, .false., "general", nfail)
    call test_pattern_change(nfail)
    call test_options_change(nfail)
    call test_failure_recovery(nfail)
    call test_unfactored(nfail)
    call test_raw_ownership(nfail)

    call test_mpi_finalize()
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

    subroutine make_solver(solver, sym, refine)
        type(sparse_solver_t), intent(out) :: solver
        integer, intent(in) :: sym
        logical, intent(in) :: refine

        solver%backend_id = FORTSPARSE_BACKEND_MUMPS
        solver%refine = refine
        solver%options%symmetry = sym
        solver%options%nthreads = 2
    end subroutine make_solver

    ! Same pattern, new symmetric values: off-diagonals scaled, diagonal
    ! shifted; keeps the Laplacian positive definite.
    subroutine change_values(A, B)
        type(csc_t), intent(in)  :: A
        type(csc_t), intent(out) :: B

        integer :: j, e

        B = A
        do j = 1, A%ncol
            do e = A%col_ptr(j), A%col_ptr(j + 1) - 1
                if (A%row_idx(e) == j) then
                    B%val(e) = 1.1_dp*A%val(e) + 0.01_dp*real(mod(j, 7), dp)
                else
                    B%val(e) = 0.8_dp*A%val(e)
                end if
            end do
        end do
    end subroutine change_values

    function rhs(n) result(b)
        integer, intent(in) :: n
        real(dp) :: b(n)
        integer :: i

        do i = 1, n
            b(i) = sin(0.37_dp*real(i, dp)) + 0.5_dp
        end do
    end function rhs

    real(dp) function reldiff(x, y)
        real(dp), intent(in) :: x(:), y(:)

        reldiff = maxval(abs(x - y))/maxval(abs(y))
    end function reldiff

    subroutine test_values(sym, refine, name, nfail)
        integer,      intent(in)    :: sym
        logical,      intent(in)    :: refine
        character(*), intent(in)    :: name
        integer,      intent(inout) :: nfail

        type(csc_t) :: A, B, Bn
        type(sparse_solver_t) :: s_re, s_fresh
        type(fortsparse_status_t) :: st
        type(sparse_factor_stats_t) :: stats
        real(dp), allocatable :: b0(:), x1(:), x2(:)
        logical :: reused
        integer :: k

        call laplace2d(24, A)
        allocate (b0(A%ncol), x1(A%ncol), x2(A%ncol))
        b0 = rhs(A%ncol)
        call make_solver(s_re, sym, refine)
        call sparse_factor(s_re, A, st)
        call check_true(name//" initial factor", status_ok(st), nfail)
        ! Three successive value changes on the one analysis.
        B = A
        do k = 1, 3
            call change_values(B, Bn)
            B = Bn
            call sparse_refactor(s_re, B, st, reused)
            call check_true(name//" refactor status", status_ok(st), nfail)
            call check_true(name//" analysis reused", reused, nfail)
            call sparse_solve(s_re, b0, x1, st)
            call check_true(name//" solve after refactor", status_ok(st), nfail)
            call make_solver(s_fresh, sym, refine)
            call sparse_factor(s_fresh, B, st)
            call sparse_solve(s_fresh, b0, x2, st)
            call check_true(name//" refactor equals fresh factor", &
                reldiff(x1, x2) < 1.0e-12_dp, nfail)
            call sparse_destroy(s_fresh)
        end do
        call sparse_factor_stats(s_re, stats, st)
        call check_true(name//" refactor reports no analysis time", &
            status_ok(st) .and. stats%analysis_s == 0.0_dp, nfail)
        call sparse_destroy(s_re)
    end subroutine test_values

    ! Same n, different pattern (tridiagonal vs 5-point) and different n.
    subroutine test_pattern_change(nfail)
        integer, intent(inout) :: nfail

        type(csc_t) :: A, T, L
        type(sparse_solver_t) :: s
        type(fortsparse_status_t) :: st
        real(dp), allocatable :: b0(:), x(:), xr(:)
        integer, allocatable :: r(:), c(:)
        real(dp), allocatable :: v(:)
        logical :: reused
        integer :: n, i, k

        call laplace2d(20, A)
        n = A%ncol
        allocate (r(3*n), c(3*n), v(3*n))
        k = 0
        do i = 1, n
            k = k + 1; r(k) = i; c(k) = i; v(k) = 3.0_dp
            if (i > 1) then
                k = k + 1; r(k) = i; c(k) = i - 1; v(k) = -1.0_dp
            end if
            if (i < n) then
                k = k + 1; r(k) = i; c(k) = i + 1; v(k) = -1.0_dp
            end if
        end do
        call csc_from_triplet(n, n, r(1:k), c(1:k), v(1:k), T, st)
        b0 = rhs(n)
        allocate (x(n), xr(n))
        call make_solver(s, FORTSPARSE_SYM_SPD, .false.)
        call sparse_factor(s, A, st)
        call sparse_refactor(s, T, st, reused)
        call check_true("pattern change status", status_ok(st), nfail)
        call check_true("pattern change detected", .not. reused, nfail)
        call sparse_solve(s, b0, x, st)
        call check_true("pattern change solves", status_ok(st), nfail)
        call sparse_destroy(s)
        call make_solver(s, FORTSPARSE_SYM_SPD, .false.)
        call sparse_factor(s, T, st)
        call sparse_solve(s, b0, xr, st)
        call check_true("pattern change equals fresh factor", &
            reldiff(x, xr) < 1.0e-12_dp, nfail)
        ! The new pattern is now the reference: a value change reuses again.
        call sparse_refactor(s, T, st, reused)
        call check_true("same pattern after fallback reused", reused, nfail)
        ! A different size is also a pattern change.
        call laplace2d(12, L)
        call sparse_refactor(s, L, st, reused)
        call check_true("size change detected", status_ok(st) .and. .not. reused, &
            nfail)
        deallocate (x)
        allocate (x(L%ncol))
        call sparse_solve(s, rhs(L%ncol), x, st)
        call check_true("size change solves", status_ok(st), nfail)
        call sparse_destroy(s)
    end subroutine test_pattern_change

    ! A different symmetry option on the same pattern cannot reuse the
    ! analysis either.
    subroutine test_options_change(nfail)
        integer, intent(inout) :: nfail

        type(csc_t) :: A
        type(sparse_solver_t) :: s
        type(fortsparse_status_t) :: st
        logical :: reused

        call laplace2d(16, A)
        call make_solver(s, FORTSPARSE_SYM_SPD, .false.)
        call sparse_factor(s, A, st)
        s%options%symmetry = FORTSPARSE_SYM_GENERAL
        call sparse_refactor(s, A, st, reused)
        call check_true("option change status", status_ok(st), nfail)
        call check_true("option change detected", .not. reused, nfail)
        call sparse_destroy(s)
    end subroutine test_options_change

    ! SPD refactor of a negative definite matrix fails as singular, leaves the
    ! solver unfactored, and the next refactor with good values recovers.
    subroutine test_failure_recovery(nfail)
        integer, intent(inout) :: nfail

        type(csc_t) :: A, Bad
        type(sparse_solver_t) :: s
        type(fortsparse_status_t) :: st
        real(dp), allocatable :: b0(:), x(:)
        logical :: reused

        call laplace2d(16, A)
        b0 = rhs(A%ncol)
        allocate (x(A%ncol))
        Bad = A
        Bad%val = -A%val
        call make_solver(s, FORTSPARSE_SYM_SPD, .false.)
        call sparse_factor(s, A, st)
        call sparse_refactor(s, Bad, st, reused)
        call check_true("indefinite refactor reported singular", &
            st%code == FORTSPARSE_SINGULAR, nfail)
        call sparse_solve(s, b0, x, st)
        call check_true("failed refactor leaves solver unfactored", &
            .not. status_ok(st), nfail)
        call sparse_refactor(s, A, st, reused)
        call check_true("recovery status", status_ok(st), nfail)
        call check_true("recovery is a full analysis", .not. reused, nfail)
        call sparse_solve(s, b0, x, st)
        call check_true("recovery solves", status_ok(st), nfail)
        call sparse_destroy(s)
    end subroutine test_failure_recovery

    subroutine test_unfactored(nfail)
        integer, intent(inout) :: nfail

        type(csc_t) :: A
        type(sparse_solver_t) :: s
        type(fortsparse_status_t) :: st
        real(dp), allocatable :: x(:)
        logical :: reused

        call laplace2d(10, A)
        allocate (x(A%ncol))
        call make_solver(s, FORTSPARSE_SYM_SPD, .false.)
        call sparse_refactor(s, A, st, reused)
        call check_true("refactor on fresh solver factors", &
            status_ok(st) .and. .not. reused, nfail)
        call sparse_solve(s, rhs(A%ncol), x, st)
        call check_true("fresh-solver refactor solves", status_ok(st), nfail)
        call sparse_destroy(s)
    end subroutine test_unfactored

    ! Raw arrays may have extra capacity or strides, and remain caller-owned.
    subroutine test_raw_ownership(nfail)
        integer, intent(inout) :: nfail

        integer, parameter :: n = 3, nz = 7
        integer, parameter :: cp(n + 1) = [1, 3, 6, 8]
        integer, parameter :: ri(nz) = [1, 2, 1, 2, 3, 2, 3]
        real(dp), parameter :: values(nz) = [4.0_dp, -2.0_dp, -1.0_dp, &
            5.0_dp, -0.5_dp, -3.0_dp, 6.0_dp]
        real(dp), parameter :: dense(n, n) = reshape([4.0_dp, -2.0_dp, 0.0_dp, &
            -1.0_dp, 5.0_dp, -0.5_dp, 0.0_dp, -3.0_dp, 6.0_dp], [n, n])
        real(dp), parameter :: exact(n) = [1.0_dp, -2.0_dp, 0.5_dp]
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        integer, allocatable :: ptr(:), idx(:)
        real(dp), allocatable :: val(:)
        real(dp) :: b(n), x(n), changed(nz)
        integer :: iref, stride, phase
        logical :: reused

        changed = 1.5_dp*values
        changed([1, 4, 7]) = changed([1, 4, 7]) + 1.0_dp
        do iref = 0, 1
            do stride = 1, 2
                call make_solver(solver, FORTSPARSE_SYM_GENERAL, iref == 1)
                ! The initially unrefined case switches off -> on -> off.
                do phase = 1, 3
                    allocate (ptr(stride*(n + 2)), idx(stride*(nz + 1)), &
                        val(stride*(nz + 1)))
                    ptr = -1; idx = -1; val = -999.0_dp
                    ptr(1:stride*n + 1:stride) = cp
                    idx(1:stride*(nz - 1) + 1:stride) = ri
                    if (phase == 1) then
                        val(1:stride*(nz - 1) + 1:stride) = values
                        b = matmul(dense, exact)
                        call sparse_factor(solver, n, n, nz, ptr(1::stride), &
                            idx(1::stride), val(1::stride), status)
                    else
                        solver%refine = (phase == 2)
                        val(1:stride*(nz - 1) + 1:stride) = changed
                        b = 1.5_dp*matmul(dense, exact) + exact
                        call sparse_refactor(solver, n, n, nz, ptr(1::stride), &
                            idx(1::stride), val(1::stride), status, reused)
                        call check_true("raw refactor reuses saved pattern", &
                            reused, nfail)
                    end if
                    call check_true("raw factor/refactor status", &
                        status_ok(status), nfail)
                    ! Both factors and refinement data must outlive the input.
                    ptr = 0; idx = 0; val = 0.0_dp
                    deallocate (ptr, idx, val)
                    if (.not. status_ok(status)) exit
                    call sparse_solve(solver, b, x, status)
                    call check_true("raw solve after input release", &
                        status_ok(status), nfail)
                    if (status_ok(status)) &
                        call check_true("raw nonsymmetric known solution", &
                            reldiff(x, exact) < 1.0e-12_dp, nfail)
                end do
                call sparse_destroy(solver)
            end do
        end do
    end subroutine test_raw_ownership

end program test_fortsparse_mumps_refactor
