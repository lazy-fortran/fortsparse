program test_fortsparse_ldlt_multifrontal
    ! Multifrontal LDL^T (native CPU, FORTSPARSE_SYM_INDEFINITE). Oracles: the
    ! analytic spectrum of the shifted 5-point Dirichlet Laplacian (inertia),
    ! backward errors of solves against the input matrix, the scalar
    ! Bunch--Kaufman reference (NATIVE_REFERENCE) on a dense matrix with zero
    ! diagonal and on block-tridiagonal matrices, Sylvester's law on
    ! saddle-point matrices [[H, B^T], [B, 0]] whose zero block forces 2x2 and
    ! delayed pivots, and an exactly singular bordered matrix whose count
    ! must stay uncertified under row scaling.
    use fortsparse
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    interface
        integer(c_int) function amd_order(n, nnz, ptr, row, perm, inverse) &
                bind(C, name='fsparse_native_order_amd')
            import :: c_int
            integer(c_int), value :: n, nnz
            integer(c_int), intent(in) :: ptr(*), row(*)
            integer(c_int), intent(out) :: perm(*), inverse(*)
        end function amd_order
        ! The same ordering on 64-bit indices.
        integer(c_int) function amd_order_wide(n, nnz, ptr, row, perm, inverse) &
                bind(C, name='fsparse_native_order_amd_general_wide')
            import :: c_int
            integer(c_int), value :: n, nnz
            integer(c_int), intent(in) :: ptr(*), row(*)
            integer(c_int), intent(out) :: perm(*), inverse(*)
        end function amd_order_wide
    end interface
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
    integer :: order(n), inverse(n), code
    integer(i8) :: fill32, fill64
    logical :: reused
    failures = 0

    ! Shifted Laplacian: inertia from the analytic eigenvalues.
    do threads = 1, 1
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
            call expect(stats%native_factor_family == 'LDLT_MULTIFRONT', 'factor family')
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

    ! The same with the AMD ordering applied explicitly (natural order of
    ! P^T A P), on 32-bit and on 64-bit indices. Inertia and its reliability
    ! must not depend on the ordering, and both index widths order alike.
    call laplacian(grid, 1.37_dp, a)
    call scale_rows(grid, a)
    code = amd_order(n, int(a%nnz), a%col_ptr, a%row_idx, order, inverse)
    call expect(code == 0, 'AMD order')
    if (code == 0) fill32 = ordered_inertia(a, order, 'AMD')
    code = amd_order_wide(n, int(a%nnz), a%col_ptr, a%row_idx, order, inverse)
    call expect(code == 0, '64-bit AMD order')
    if (code == 0) then
        fill64 = ordered_inertia(a, order, '64-bit AMD')
        call expect(fill32 == fill64, '32- and 64-bit AMD fill agree')
    end if

    ! Negative control: the dyadic shift 11/8 and an extra row and column
    ! c = A w, d = w^T A w with integer w, all exact, then the same power-of-
    ! two scalings. The matrix is exactly singular, so one pivot is rounding
    ! of its own row, whether that row is tiny or huge against the rest. The
    ! count must not be certified.
    do k = 1, size(powers)
        do threads = 1, 1
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

    ! Dense indefinite matrix with zero diagonal: one front, 2x2 pivots.
    call dense_zero_diagonal(dense, a)
    call configure(solver, 1)
    call sparse_factor(solver, a, status)
    call expect(status_ok(status), 'factor of the zero-diagonal dense matrix')
    if (status_ok(status)) then
        call residual_check(solver, a, 1.0e-12_dp)
        reference%backend_id = FORTSPARSE_BACKEND_NATIVE_REFERENCE
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

    ! Block-tridiagonal dense blocks: a chain of fronts.
    call block_tridiagonal(8, 40, a)
    reference%backend_id = FORTSPARSE_BACKEND_NATIVE_REFERENCE
    reference%options%symmetry = FORTSPARSE_SYM_INDEFINITE
    call sparse_factor(reference, a, status)
    call expect(status_ok(status), 'scalar reference of the block-tridiagonal matrix')
    do k = 1, 2
        call configure(solver, 1)
        solver%options%ordering = merge(FORTSPARSE_ORDER_NATURAL, &
            FORTSPARSE_ORDER_DEFAULT, k == 1)
        call sparse_factor(solver, a, status)
        call expect(status_ok(status), 'factor of the block-tridiagonal matrix')
        if (.not. status_ok(status)) cycle
        call residual_check(solver, a, 1.0e-12_dp)
        call sparse_inertia(solver, nneg, nzero, npos, reliability, status)
        block
            integer(i8) :: rneg, rzero, rpos
            integer :: rrel
            call sparse_inertia(reference, rneg, rzero, rpos, rrel, status)
            call expect(rneg == nneg .and. rzero == nzero .and. rpos == npos, &
                'block-tridiagonal inertia agrees with the scalar reference')
        end block
        call expect(reliability == FORTSPARSE_INERTIA_RELIABLE, &
            'block-tridiagonal inertia reliable')
        call sparse_destroy(solver)
    end do
    call sparse_destroy(reference)

    ! Saddle points: Sylvester's law gives the inertia of [[H, B^T], [B, 0]]
    ! with H positive definite and B of full row rank as (nh, mb, 0) for
    ! (+, -, 0). Every constraint row has a zero diagonal. In the natural
    ! order with the constraints first, each of them fails as a 1x1 pivot
    ! and is delayed or paired; AMD orders it next to a matched neighbour.
    do k = 1, 4
        call saddle_point(grid, 7*k, a, first=k > 2)
        call configure(solver, 1)
        if (k > 2) solver%options%ordering = FORTSPARSE_ORDER_NATURAL
        call sparse_factor(solver, a, status)
        call expect(status_ok(status), 'factor of the saddle-point matrix')
        if (.not. status_ok(status)) cycle
        call sparse_inertia(solver, nneg, nzero, npos, reliability, status)
        call expect(nneg == count_constraints(grid, 7*k) .and. &
            nzero == 0 .and. npos == n, 'saddle-point inertia (n, m, 0)')
        call expect(reliability == FORTSPARSE_INERTIA_RELIABLE, &
            'saddle-point inertia reliable')
        call residual_check(solver, a, 1.0e-12_dp)
        ! new values on the same pattern: refactor reuses the analysis
        call saddle_point(grid, 7*k, a, first=k > 2, weight=3.0_dp)
        call sparse_refactor(solver, a, status, reused)
        call expect(status_ok(status) .and. reused, 'saddle-point refactor')
        call residual_check(solver, a, 1.0e-12_dp)
        call sparse_destroy(solver)
    end do

    ! A strict host cap selects the scalar factor.
    call laplacian(grid, 1.37_dp, a)
    call configure(solver, 1)
    solver%options%host_memory_limit_bytes = 1000000000_i8
    call sparse_factor(solver, a, status)
    call expect(status_ok(status), 'capped LDL^T factors')
    call sparse_factor_stats(solver, stats, status)
    call expect(stats%native_factor_family == 'LDLT_REAL', 'capped LDL^T is the scalar family')
    call sparse_destroy(solver)

    if (failures /= 0) then
        print '(a,i0)', 'multifrontal LDL^T failures: ', failures
        error stop 1
    end if
    print '(a)', 'multifrontal LDL^T: ok'
contains
    integer(i8) function ordered_inertia(m, perm, what) result(fill)
        ! Factor P^T m P (perm(new) = old) in its natural order and check its
        ! inertia; returns nnz(L).
        type(csc_t), intent(in) :: m
        integer, intent(in) :: perm(:)
        character(*), intent(in) :: what
        type(csc_t) :: b
        type(sparse_solver_t) :: s
        integer :: rows(m%nnz), cols(m%nnz), position(m%ncol), j, p
        real(dp) :: vals(m%nnz)
        position(perm) = [(j, j=1, m%ncol)]
        do j = 1, m%ncol
            do p = m%col_ptr(j), m%col_ptr(j + 1) - 1
                rows(p) = position(m%row_idx(p))
                cols(p) = position(j)
                vals(p) = m%val(p)
            end do
        end do
        call csc_from_triplet(m%nrow, m%ncol, rows, cols, vals, b, status)
        call configure(s, 1)
        s%options%ordering = FORTSPARSE_ORDER_NATURAL
        fill = -1
        call sparse_factor(s, b, status)
        call expect(status_ok(status), 'factor in the '//what//' order')
        if (status_ok(status)) then
            call sparse_factor_stats(s, stats, status)
            fill = stats%symbolic_nnz
            call sparse_inertia(s, nneg, nzero, npos, reliability, status)
            call expect(nneg == below(grid, 1.37_dp) .and. nzero == 0, &
                what//': row-scaled inertia matches the spectrum')
            call expect(reliability == FORTSPARSE_INERTIA_RELIABLE, &
                what//': row-scaled inertia reliable')
        end if
        call sparse_destroy(s)
    end function ordered_inertia

    subroutine configure(s, threads)
        type(sparse_solver_t), intent(inout) :: s
        integer, intent(in) :: threads
        s%backend_id = FORTSPARSE_BACKEND_NATIVE
        s%options%symmetry = FORTSPARSE_SYM_INDEFINITE
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

    subroutine block_tridiagonal(nb, w, m)
        ! Dense w by w blocks on the three block diagonals; indefinite.
        integer, intent(in) :: nb, w
        type(csc_t), intent(out) :: m
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        integer :: i, j, k
        allocate (rows(3*w*w*nb), cols(3*w*w*nb), vals(3*w*w*nb))
        k = 0
        do j = 1, nb*w
            do i = max(1, ((j - 1)/w - 1)*w + 1), min(nb*w, ((j - 1)/w + 2)*w)
                k = k + 1; rows(k) = i; cols(k) = j
                vals(k) = cos(0.3_dp*(i + j)) + 1.0_dp/(1 + abs(i - j))
                if (i == j) vals(k) = vals(k) + 3*cos(0.7_dp*i)
            end do
        end do
        call csc_from_triplet(nb*w, nb*w, rows(:k), cols(:k), vals(:k), m, status)
    end subroutine block_tridiagonal

    integer function count_constraints(g, stride) result(mb)
        integer, intent(in) :: g, stride
        mb = (g*g - 1)/stride + 1
    end function count_constraints

    subroutine saddle_point(g, stride, m, first, weight)
        ! H = shifted Laplacian (positive definite, sigma = -0.5) on g*g
        ! unknowns, B with one row per stride-th unknown i: B(r, i) = 1,
        ! B(r, i+1) = w (w = weight, default 0.5), so B has full row rank.
        ! With first, the constraint rows are numbered before H.
        integer, intent(in) :: g, stride
        type(csc_t), intent(out) :: m
        logical, intent(in) :: first
        real(dp), intent(in), optional :: weight
        type(csc_t) :: h
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        integer :: nh, mb, j, p, k, r, i
        real(dp) :: w
        w = 0.5_dp
        if (present(weight)) w = weight
        call laplacian(g, -0.5_dp, h)
        nh = g*g
        mb = count_constraints(g, stride)
        allocate (rows(h%nnz + 4*mb), cols(h%nnz + 4*mb), vals(h%nnz + 4*mb))
        k = 0
        do j = 1, nh
            do p = h%col_ptr(j), h%col_ptr(j + 1) - 1
                k = k + 1; rows(k) = h%row_idx(p); cols(k) = j; vals(k) = h%val(p)
            end do
        end do
        do r = 1, mb
            i = 1 + (r - 1)*stride
            k = k + 1; rows(k) = nh + r; cols(k) = i; vals(k) = 1
            k = k + 1; rows(k) = i; cols(k) = nh + r; vals(k) = 1
            if (i < nh) then
                k = k + 1; rows(k) = nh + r; cols(k) = i + 1; vals(k) = w
                k = k + 1; rows(k) = i + 1; cols(k) = nh + r; vals(k) = w
            end if
        end do
        if (first) then
            ! constraints first: unknown v of [H, B^T; B, 0] becomes v + mb
            rows(:k) = merge(rows(:k) - nh, rows(:k) + mb, rows(:k) > nh)
            cols(:k) = merge(cols(:k) - nh, cols(:k) + mb, cols(:k) > nh)
        end if
        call csc_from_triplet(nh + mb, nh + mb, rows(:k), cols(:k), vals(:k), m, status)
    end subroutine saddle_point

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
end program test_fortsparse_ldlt_multifrontal
