program test_fortsparse_ldlt_supernodal
    ! Supernodal LDL^T (options%ldlt_supernodal, FORTSPARSE_SYM_INDEFINITE).
    ! Oracles: the analytic spectrum of the shifted 5-point Dirichlet
    ! Laplacian (inertia), backward errors of solves against the input matrix, the
    ! unrestricted scalar Bunch--Kaufman reference on a dense indefinite
    ! matrix with zero diagonal (2x2 pivots inside one supernode), team
    ! sizes 1 and 4 (same inertia, small residuals), and an exactly singular
    ! bordered matrix whose count must stay uncertified under row scaling.
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
    integer, parameter :: powers(4) = [0, 0, -30, 30]
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

    ! Rows scaled over ten decades, S A S with powers of two: a congruence,
    ! so the inertia is that of the shifted Laplacian, and every pivot stays
    ! far above the rounding threshold of its own row.
    call laplacian(grid, 1.37_dp, a)
    call scale_rows(grid, a)
    call configure(solver, 1)
    call sparse_factor(solver, a, status)
    call expect(status_ok(status), 'factor of the row-scaled Laplacian')
    if (status_ok(status)) then
        call sparse_inertia(solver, nneg, nzero, npos, reliability, status)
        call expect(nneg == below(grid, 1.37_dp) .and. nzero == 0, &
            'row-scaled inertia matches the spectrum')
        call expect(reliability == FORTSPARSE_INERTIA_RELIABLE, &
            'row-scaled inertia reliable')
    end if
    call sparse_destroy(solver)

    ! Negative control: the dyadic shift 11/8 and an extra row and column
    ! c = A w, d = w^T A w with integer w, all exact, then the same power-of-
    ! two scalings. The matrix is exactly singular, so one pivot is rounding
    ! of its own row, whether that row is tiny or huge against the rest. The
    ! count must not be certified.
    do k = 1, size(powers)
        do threads = 1, 4, 3
            call laplacian(grid, 1.375_dp, a)
            call add_dependent_row(a)
            if (k > 1) call scale_rows(grid, a)
            call scale_last(a, powers(k))
            call configure(solver, threads)
            call sparse_factor(solver, a, status)
            if (status_ok(status)) then
                call sparse_inertia(solver, nneg, nzero, npos, reliability, status)
                call expect(reliability /= FORTSPARSE_INERTIA_RELIABLE, &
                    'rounding-level pivot is not certified')
            else
                call expect(status%code == FORTSPARSE_SINGULAR, &
                    'rounding-level pivot fails as singular only')
            end if
            call sparse_destroy(solver)
        end do
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

    subroutine scale_rows(g, m)
        ! m := S m S, S = 2**k with k rising from -17 to 17 over the grid rows.
        integer, intent(in) :: g
        type(csc_t), intent(inout) :: m
        integer :: j, p, k(m%ncol)
        k = [(nint(-17 + 34.0_dp*((j - 1)/g)/(g - 1)), j=1, m%ncol)]
        do j = 1, m%ncol
            do p = m%col_ptr(j), m%col_ptr(j + 1) - 1
                m%val(p) = scale(m%val(p), k(m%row_idx(p)) + k(j))
            end do
        end do
    end subroutine scale_rows

    subroutine add_dependent_row(m)
        ! m := [m, c; c^T, d] with c = m w and d = w^T m w for the integer
        ! vector w = 1 + mod(i, 3): exact for dyadic entries of m.
        type(csc_t), intent(inout) :: m
        integer :: rows(m%nnz + 2*m%nrow + 1), cols(m%nnz + 2*m%nrow + 1), i, j, p, k, extra
        real(dp) :: vals(m%nnz + 2*m%nrow + 1), c(m%nrow), d
        extra = m%ncol + 1
        k = 0
        c = 0
        do j = 1, m%ncol
            do p = m%col_ptr(j), m%col_ptr(j + 1) - 1
                k = k + 1; rows(k) = m%row_idx(p); cols(k) = j; vals(k) = m%val(p)
                c(m%row_idx(p)) = c(m%row_idx(p)) + m%val(p)*(1 + mod(j, 3))
            end do
        end do
        d = 0
        do i = 1, m%nrow
            d = d + c(i)*(1 + mod(i, 3))
            k = k + 1; rows(k) = i; cols(k) = extra; vals(k) = c(i)
            k = k + 1; rows(k) = extra; cols(k) = i; vals(k) = c(i)
        end do
        k = k + 1; rows(k) = extra; cols(k) = extra; vals(k) = d
        call csc_from_triplet(extra, extra, rows(:k), cols(:k), vals(:k), m, status)
    end subroutine add_dependent_row

    subroutine scale_last(m, power)
        ! Last row and column of m times 2**power.
        type(csc_t), intent(inout) :: m
        integer, intent(in) :: power
        integer :: j, p
        do j = 1, m%ncol
            do p = m%col_ptr(j), m%col_ptr(j + 1) - 1
                if (m%row_idx(p) == m%nrow) m%val(p) = scale(m%val(p), power)
                if (j == m%ncol) m%val(p) = scale(m%val(p), power)
            end do
        end do
    end subroutine scale_last

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
