program test_fortsparse_native_ldlt
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csc, only: csc_t, csc_from_triplet, csc_matvec
    use fortsparse_status, only: fortsparse_status_t, FORTSPARSE_OK, &
                                 FORTSPARSE_SINGULAR, FORTSPARSE_INVALID_MATRIX, &
                                 FORTSPARSE_NOT_FACTORED, FORTSPARSE_RESOURCE_LIMIT
    use fortsparse_native_ldlt, only: native_ldlt_t
    use fortsparse_ldlt_leaves, only: ldlt_pair_jvp, ldlt_pair_vjp, &
                                      ldlt_update_jvp, ldlt_update_vjp
    implicit none
    type(native_ldlt_t) :: factor, compact
    type(csc_t) :: a
    type(fortsparse_status_t) :: status
    integer :: i, j, n, trial, permutation(2), before(2)
    integer, allocatable :: rows(:), columns(:)
    integer(i8) :: peak, low_peak, retained
    real(dp) :: pair(2, 2), rhs2(2), result2(2), three(3, 3)
    real(dp), allocatable :: dense(:, :), expected(:), b(:), x(:), values(:)

    a%nrow = 2
    a%ncol = 2
    call factor%factor(a, status)
    if (status%code /= FORTSPARSE_INVALID_MATRIX) error stop 'unallocated CSC accepted'
    ! An unpivoted LDL cannot factor this invertible matrix.
    pair = reshape([0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp], [2, 2])
    call check_dense(pair)
    deallocate (a%val)
    allocate (a%val(0))
    call factor%factor(a, status)
    if (status%code /= FORTSPARSE_INVALID_MATRIX) error stop 'short values accepted'
    call check_dense(1.0e-300_dp*pair)
    call check_dense(1.0e300_dp*pair)
    ! A nonadjacent 2x2 pivot with a trailing row needs both maps and columns.
    three = reshape([0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 3.0_dp, 4.0_dp, &
                     2.0_dp, 4.0_dp, 0.0_dp], [3, 3])
    call check_dense(three)
    ! Symmetric 1x1 interchange, then same-pattern numerical repivoting.
    pair = reshape([0.01_dp, 1.0_dp, 1.0_dp, 10.0_dp], [2, 2])
    call from_dense(pair, a)
    call factor%factor(a, status)
    call ok('interchange factor')
    call factor%permutation(before, status)
    call ok('permutation')
    if (before(1) /= 2) error stop '1x1 symmetric pivot interchange missing'
    pair = reshape([10.0_dp, 1.0_dp, 1.0_dp, 0.01_dp], [2, 2])
    call from_dense(pair, a)
    call factor%refactor(a, status)
    call ok('same pattern refactor')
    call factor%permutation(permutation, status)
    call ok('repivot permutation')
    if (permutation(1) /= 1) error stop 'numeric refactor did not repivot'
    rhs2 = matmul(pair, [2.0_dp, -3.0_dp])
    call factor%solve(rhs2, result2, status, transpose=.true.)
    call ok('retained transpose solve')
    if (maxval(abs(result2 - [2.0_dp, -3.0_dp])) > 1.0e-12_dp) &
        error stop 'repivot solve'

    ! Deterministic mixed-sign symmetric matrices: independent known solution
    ! and dense matmul residuals, with 1x1/2x2 pivots after several steps.
    do trial = 1, 40
        n = 2 + mod(trial, 14)
        allocate (dense(n, n))
        do j = 1, n
            do i = 1, j
                dense(i, j) = sin(real(7*i + 13*j + trial, dp))
                dense(j, i) = dense(i, j)
            end do
            dense(j, j) = dense(j, j) + 0.2_dp*(-1.0_dp)**j
        end do
        call check_dense(dense)
        deallocate (dense)
    end do

    ! Arrowhead elimination creates a sparse clique and exercises row growth.
    n = 25
    allocate (dense(n, n), expected(n), b(n), x(n))
    dense = 0.0_dp
    dense(1, 1) = 100.0_dp
    do i = 2, n
        dense(i, i) = 2.0_dp*(-1.0_dp)**i
        dense(1, i) = 0.25_dp
        dense(i, 1) = 0.25_dp
    end do
    expected = [(sin(real(i, dp)), i=1, n)]
    b = matmul(dense, expected)
    call from_dense(dense, a)
    call factor%factor(a, status)
    call ok('growing sparse factor')
    peak = factor%peak_bytes()
    retained = factor%retained_bytes()
    if (peak <= retained) error stop 'factor scratch omitted from peak'
    call compact%factor(a, status, low_memory=.true.)
    call ok('low-memory factor')
    low_peak = compact%peak_bytes()
    if (low_peak > peak) error stop 'low-memory growth increased payload peak'
    call compact%solve(b, x, status)
    call ok('low-memory solve')
    if (maxval(abs(x - expected)) > 1.0e-12_dp) error stop 'low-memory accuracy'
    call factor%factor(a, status, host_limit_bytes=peak)
    call ok('exact peak cap factor')
    if (factor%peak_bytes() /= peak) error stop 'nonreproducible payload peak'
    retained = factor%retained_bytes()
    x = b
    call factor%solve_inplace(x, status)
    call ok('inplace under exact memory cap')
    if (maxval(abs(x - expected)) > 1.0e-12_dp) error stop 'capped inplace accuracy'
    if (factor%retained_bytes() /= retained) error stop 'capped inplace grew storage'
    if (factor%required_bytes() /= peak) error stop 'successful required diagnostic'
    call factor%refactor(a, status, host_limit_bytes=peak)
    call ok('cap destructive refactor')
    call factor%factor(a, status, host_limit_bytes=peak - 1_i8)
    if (status%code /= FORTSPARSE_RESOURCE_LIMIT) error stop 'growth cap ignored'
    if (factor%retained_bytes() > peak - 1_i8) error stop 'cap overshot in growth'
    if (factor%available_bytes() /= peak - 1_i8) error stop 'available cap diagnostic'
    if (factor%required_bytes() <= factor%available_bytes()) &
        error stop 'failed growth requirement not reported'
    call factor%solve(b, x, status)
    if (status%code /= FORTSPARSE_NOT_FACTORED) error stop 'failed factors usable'
    call factor%factor(a, status, host_limit_bytes=1_i8)
    if (status%code /= FORTSPARSE_RESOURCE_LIMIT) error stop 'analysis cap ignored'
    if (factor%available_bytes() /= 1_i8) error stop 'analysis cap diagnostic'
    if (factor%required_bytes() <= 1_i8) error stop 'analysis requirement diagnostic'
    call factor%factor(a, status, host_limit_bytes=2_i8**40)
    call ok('int64 cap')
    deallocate (dense, expected, b, x)

    ! Large disconnected diagonal family: O(n) owned payload and sparse factors.
    n = 3000
    allocate (rows(n), columns(n), values(n), expected(n), b(n), x(n))
    do i = 1, n
        rows(i) = i
        columns(i) = i
        values(i) = (-1.0_dp)**i
        expected(i) = sin(real(i, dp))
    end do
    call csc_from_triplet(n, n, rows, columns, values, a, status)
    call ok('large diagonal CSC')
    call factor%factor(a, status, host_limit_bytes=3000000_i8)
    call ok('large sparse factor')
    if (factor%factor_nnz() /= int(n, i8)) error stop 'diagonal fill'
    b = csc_matvec(a, expected)
    call factor%solve(b, x, status)
    call ok('large sparse solve')
    if (maxval(abs(x - expected)) > 1.0e-14_dp) error stop 'large sparse accuracy'
    deallocate (rows, columns, values, expected, b, x)

    pair = 0.0_dp
    call from_dense(pair, a)
    call factor%factor(a, status)
    if (status%code /= FORTSPARSE_SINGULAR) error stop 'singular accepted'
    pair = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [2, 2])
    call from_dense(pair, a)
    call factor%factor(a, status)
    if (status%code /= FORTSPARSE_INVALID_MATRIX) error stop 'asymmetry accepted'
    a%val(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call factor%factor(a, status)
    if (status%code /= FORTSPARSE_INVALID_MATRIX) error stop 'nonfinite accepted'
    ! Singular refactor leaves no usable numerical factors but permits recovery.
    pair = reshape([2.0_dp, 1.0_dp, 1.0_dp, -2.0_dp], [2, 2])
    call from_dense(pair, a)
    call factor%factor(a, status)
    call ok('recovery initial factor')
    a%val = 1.0_dp
    call factor%refactor(a, status)
    if (status%code /= FORTSPARSE_SINGULAR) error stop 'singular refactor accepted'
    call factor%solve([1.0_dp, 1.0_dp], result2, status)
    if (status%code /= FORTSPARSE_NOT_FACTORED) error stop 'singular factors usable'
    call from_dense(pair, a)
    call factor%refactor(a, status)
    call ok('same-pattern recovery')
    ! Pattern mismatch is rejected before dropping the prior good factors.
    a%row_idx(2) = 1
    call factor%refactor(a, status)
    if (status%code /= FORTSPARSE_INVALID_MATRIX) error stop 'bad pattern accepted'
    call factor%solve(matmul(pair, [2.0_dp, -3.0_dp]), result2, status)
    call ok('good factors survive rejected input')
    if (maxval(abs(result2 - [2.0_dp, -3.0_dp])) > 1.0e-12_dp) &
        error stop 'rejected refactor damaged good factors'
    call factor%free()
    if (factor%retained_bytes() /= 0_i8) error stop 'free retained payload'
    call leaf_derivatives()
contains
    subroutine ok(label)
        character(*), intent(in) :: label
        if (status%code == FORTSPARSE_OK) return
        print *, label, status%code, trim(status%msg)
        error stop 'unexpected LDLT status'
    end subroutine ok

    subroutine from_dense(matrix, csc)
        real(dp), intent(in) :: matrix(:, :)
        type(csc_t), intent(out) :: csc
        integer :: i, j, p, n
        integer, allocatable :: r(:), c(:)
        real(dp), allocatable :: v(:)
        n = size(matrix, 1)
        allocate (r(n*n), c(n*n), v(n*n))
        p = 0
        do j = 1, n
            do i = 1, n
                if (abs(matrix(i, j)) <= 0.0_dp) cycle
                p = p + 1
                r(p) = i
                c(p) = j
                v(p) = matrix(i, j)
            end do
        end do
        call csc_from_triplet(n, n, r(:p), c(:p), v(:p), csc, status)
        call ok('small dense-to-CSC fixture')
    end subroutine from_dense

    subroutine check_dense(matrix)
        real(dp), intent(in) :: matrix(:, :)
        integer :: i, n
        integer(i8) :: before_bytes
        real(dp), allocatable :: truth(:), rhs(:, :), solution(:, :)
        n = size(matrix, 1)
        allocate (truth(n), rhs(n, 2), solution(n, 2))
        truth = [(cos(real(i, dp)), i=1, n)]
        rhs(:, 1) = matmul(matrix, truth)
        rhs(:, 2) = -2.0_dp*rhs(:, 1)
        call from_dense(matrix, a)
        call factor%factor(a, status)
        call ok('symmetric indefinite factor')
        call factor%solve_block(rhs, solution, status)
        call ok('block solve')
        if (maxval(abs(solution(:, 1) - truth)) > 2.0e-10_dp) &
            error stop 'known-solution oracle'
        if (maxval(abs(matmul(matrix, solution) - rhs)) &
            > 1.0e-12_dp*max(1.0_dp, maxval(abs(rhs)))) error stop 'residual oracle'
        call factor%solve_block(rhs, solution, status, transpose=.true.)
        call ok('transpose block solve')
        before_bytes = factor%retained_bytes()
        ! Public solve JVP forms A*dx = -A_dot*x for an unchanged RHS.
        ! A_dot=A therefore has the independent answer dx=-x.
        solution(:, 1) = -rhs(:, 1)
        call factor%solve_inplace(solution(:, 1), status)
        call ok('inplace JVP-style solve')
        if (maxval(abs(solution(:, 1) + truth)) > 2.0e-10_dp) &
            error stop 'inplace equation-level tangent'
        solution(:, 1) = rhs(:, 1)
        call factor%solve_inplace(solution(:, 1), status, transpose=.true.)
        call ok('inplace transpose solve')
        if (maxval(abs(solution(:, 1) - truth)) > 2.0e-10_dp) &
            error stop 'inplace known-solution oracle'
        if (factor%retained_bytes() /= before_bytes) error stop 'inplace retained bytes'
        if (maxval(abs(solution(:, 2) + 2.0_dp*truth)) > 4.0e-10_dp) &
            error stop 'transpose known-solution oracle'
    end subroutine check_dense

    subroutine leaf_derivatives()
        real(dp) :: inputs(5), direction(5), z(2), dz(2), z2(2), bars(5)
        real(dp) :: matrix(2, 2), matrix_dot(2, 2), seed(2), value, tangent
        inputs = [2.0_dp, 0.5_dp, -1.0_dp, 0.75_dp, 1.25_dp]
        direction = [0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp, -0.5_dp]
        seed = [0.7_dp, -0.4_dp]
        call ldlt_pair_jvp(inputs(1), inputs(2), inputs(3), inputs(4), inputs(5), &
                 direction(1), direction(2), direction(3), direction(4), direction(5), &
                           z(1), z(2), dz(1), dz(2))
        matrix = reshape([inputs(1), inputs(2), inputs(2), inputs(3)], [2, 2])
        matrix_dot = reshape([direction(1), direction(2), direction(2), &
                              direction(3)], [2, 2])
        if (maxval(abs(matmul(matrix, z) - inputs(4:5))) > 1.0e-14_dp) &
            error stop 'pair value oracle'
        if (maxval(abs(matmul(matrix, dz) + matmul(matrix_dot, z) - direction(4:5))) &
            > 1.0e-14_dp) error stop 'pair equation-level tangent oracle'
        call ldlt_pair_vjp(inputs(1), inputs(2), inputs(3), inputs(4), inputs(5), &
            seed(1), seed(2), z2(1), z2(2), bars(1), bars(2), bars(3), bars(4), bars(5))
        if (abs(dot_product(seed, dz) - dot_product(bars, direction)) > 1.0e-14_dp) &
            error stop 'pair leaf adjoint pairing'
        call ldlt_update_jvp(inputs(1), inputs(2), inputs(3), inputs(4), inputs(5), &
                 direction(1), direction(2), direction(3), direction(4), direction(5), &
                             value, tangent)
        call ldlt_update_vjp(inputs(1), inputs(2), inputs(3), inputs(4), inputs(5), &
                            seed(1), value, bars(1), bars(2), bars(3), bars(4), bars(5))
        if (abs(seed(1)*tangent - dot_product(bars, direction)) > 1.0e-14_dp) &
            error stop 'rank-two update adjoint pairing'
    end subroutine leaf_derivatives
end program test_fortsparse_native_ldlt
