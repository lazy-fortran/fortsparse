program test_fortsparse_ldlt_supernodal
    ! Supernodal LDL^T (options%ldlt_supernodal, FORTSPARSE_SYM_INDEFINITE).
    ! Oracles: the analytic spectrum of the shifted 5-point Dirichlet
    ! Laplacian (inertia), backward errors of solves against the input matrix, the
    ! unrestricted scalar Bunch--Kaufman reference on a dense indefinite
    ! matrix with zero diagonal (2x2 pivots inside one supernode), and team
    ! sizes 1 and 4 (same inertia, small residuals).
    use fortsparse
    implicit none
    integer, parameter :: grid = 24, n = grid*grid, dense = 72
    real(dp), parameter :: pi = acos(-1.0_dp)
    type(csc_t) :: a
    type(sparse_solver_t) :: solver, reference
    type(fortsparse_status_t) :: status
    type(sparse_factor_stats_t) :: stats
    real(dp) :: shift
    real(dp), parameter :: shifts(3) = [0.0_dp, 1.37_dp, 4.11_dp]
    integer(i8) :: nneg, nzero, npos
    integer :: failures, reliability, threads, k
    logical :: reused
    failures = 0

    ! Shifted Laplacian: inertia from the analytic eigenvalues.
    do threads = 1, 4, 3
        do k = 1, 3
            shift = shifts(k)
            call laplacian(grid, shift, a)
            call configure(solver, threads)
            if (k == 1) then
                call sparse_factor(solver, a, status)
            else
                call sparse_refactor(solver, a, status, reused)
                call expect(reused, 'same-pattern refactor reuses the analysis')
            end if
            call expect(status_ok(status), 'factor of the shifted Laplacian')
            if (.not. status_ok(status)) cycle
            call sparse_factor_stats(solver, stats, status)
            call expect(stats%native_factor_family == 'LDLT_SUPERNODAL', 'factor family')
            call sparse_inertia(solver, nneg, nzero, npos, reliability, status)
            call expect(status_ok(status), 'inertia available')
            call expect(nneg == below(grid, shift) .and. nzero == 0 .and. &
                npos == n - below(grid, shift), 'inertia matches the spectrum')
            call expect(reliability == FORTSPARSE_INERTIA_RELIABLE, 'inertia reliable')
            call residual_check(solver, a, 1.0e-12_dp)
        end do
        call sparse_destroy(solver)
    end do

    ! Dense indefinite matrix with zero diagonal: one supernode, 2x2 pivots.
    call dense_zero_diagonal(dense, a)
    call configure(solver, 1)
    call sparse_factor(solver, a, status)
    call expect(status_ok(status), 'factor of the zero-diagonal dense matrix')
    if (status_ok(status)) then
        call residual_check(solver, a, 1.0e-12_dp)
        reference%backend_id = FORTSPARSE_BACKEND_NATIVE
        reference%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        call sparse_factor(reference, a, status)
        call expect(status_ok(status), 'scalar Bunch--Kaufman reference')
        call sparse_inertia(solver, nneg, nzero, npos, reliability, status)
        block
            integer(i8) :: rneg, rzero, rpos
            integer :: rrel
            call sparse_inertia(reference, rneg, rzero, rpos, rrel, status)
            call expect(rneg == nneg .and. rzero == nzero .and. rpos == npos, &
                'inertia agrees with the scalar reference')
        end block
        call sparse_destroy(reference)
    end if
    call sparse_destroy(solver)

    ! Explicitly unsupported: a strict host cap.
    call laplacian(grid, 1.37_dp, a)
    call configure(solver, 1)
    solver%options%host_memory_limit_bytes = 1000000000_i8
    call sparse_factor(solver, a, status)
    call expect(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, 'capped LDL^T unavailable')
    call sparse_destroy(solver)

    if (failures /= 0) then
        print '(a,i0)', 'supernodal LDL^T failures: ', failures
        error stop 1
    end if
    print '(a)', 'supernodal LDL^T: ok'
contains
    subroutine configure(s, threads)
        type(sparse_solver_t), intent(inout) :: s
        integer, intent(in) :: threads
        s%backend_id = FORTSPARSE_BACKEND_NATIVE
        s%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        s%options%ldlt_supernodal = .true.
        s%options%ordering = FORTSPARSE_ORDER_AMD
        s%options%nthreads = threads
    end subroutine configure

    subroutine expect(condition, what)
        logical, intent(in) :: condition
        character(*), intent(in) :: what
        if (condition) return
        failures = failures + 1
        print '(2a)', 'FAIL: ', what
    end subroutine expect

    subroutine residual_check(s, m, tolerance)
        ! Normwise backward error |b - A x|/(|A| |x| + |b|) (infinity norms).
        type(sparse_solver_t), intent(inout) :: s
        type(csc_t), intent(in) :: m
        real(dp), intent(in) :: tolerance
        real(dp), allocatable :: xs(:), bs(:), r(:), rowsum(:)
        real(dp) :: berr
        integer :: i, j, p
        allocate (xs(m%ncol), bs(m%ncol), r(m%ncol), rowsum(m%ncol))
        bs = [(sin(0.37_dp*i) + 0.1_dp, i=1, m%ncol)]
        call sparse_solve(s, bs, xs, status)
        call expect(status_ok(status), 'solve')
        r = bs
        rowsum = 0
        do j = 1, m%ncol
            do p = m%col_ptr(j), m%col_ptr(j + 1) - 1
                r(m%row_idx(p)) = r(m%row_idx(p)) - m%val(p)*xs(j)
                rowsum(m%row_idx(p)) = rowsum(m%row_idx(p)) + abs(m%val(p))
            end do
        end do
        berr = maxval(abs(r))/(maxval(rowsum)*maxval(abs(xs)) + maxval(abs(bs)))
        if (berr > tolerance) print '(a,es10.2)', 'backward error ', berr
        call expect(berr <= tolerance, 'backward error')
    end subroutine residual_check

    integer(i8) function below(g, sigma) result(count)
        ! Eigenvalues 4 - 2 cos(i pi/(g+1)) - 2 cos(j pi/(g+1)) below sigma.
        integer, intent(in) :: g
        real(dp), intent(in) :: sigma
        integer :: i, j
        count = 0
        do j = 1, g
            do i = 1, g
                if (4 - 2*cos(i*pi/(g + 1)) - 2*cos(j*pi/(g + 1)) < sigma) &
                    count = count + 1
            end do
        end do
    end function below

    subroutine laplacian(g, sigma, m)
        integer, intent(in) :: g
        real(dp), intent(in) :: sigma
        type(csc_t), intent(out) :: m
        integer :: rows(5*g*g), cols(5*g*g), i, j, k, c
        real(dp) :: vals(5*g*g)
        k = 0
        do j = 1, g
            do i = 1, g
                c = i + (j - 1)*g
                k = k + 1; rows(k) = c; cols(k) = c; vals(k) = 4 - sigma
                if (i > 1) then
                    k = k + 1; rows(k) = c - 1; cols(k) = c; vals(k) = -1
                end if
                if (i < g) then
                    k = k + 1; rows(k) = c + 1; cols(k) = c; vals(k) = -1
                end if
                if (j > 1) then
                    k = k + 1; rows(k) = c - g; cols(k) = c; vals(k) = -1
                end if
                if (j < g) then
                    k = k + 1; rows(k) = c + g; cols(k) = c; vals(k) = -1
                end if
            end do
        end do
        call csc_from_triplet(g*g, g*g, rows(:k), cols(:k), vals(:k), m, status)
    end subroutine laplacian

    subroutine dense_zero_diagonal(nd, m)
        integer, intent(in) :: nd
        type(csc_t), intent(out) :: m
        integer :: rows(nd*nd), cols(nd*nd), i, j, k
        real(dp) :: vals(nd*nd)
        k = 0
        do j = 1, nd
            do i = 1, nd
                k = k + 1; rows(k) = i; cols(k) = j
                vals(k) = 0
                if (i /= j) vals(k) = cos(0.3_dp*(i + j)) + 1.0_dp/(1 + abs(i - j))
            end do
        end do
        call csc_from_triplet(nd, nd, rows(:k), cols(:k), vals(:k), m, status)
    end subroutine dense_zero_diagonal
end program test_fortsparse_ldlt_supernodal
