program test_fortsparse_native_ldlt_z
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csc, only: csc_z_t, csc_from_triplet, csc_matvec
    use fortsparse_status, only: fortsparse_status_t, FORTSPARSE_OK, &
                                 FORTSPARSE_SINGULAR, FORTSPARSE_INVALID_MATRIX, &
           FORTSPARSE_NOT_FACTORED, FORTSPARSE_RESOURCE_LIMIT, FORTSPARSE_INTERNAL_ERROR
    use fortsparse_native_ldlt_z, only: native_ldlt_z_t
    use fortsparse_ldlt_z_leaves, only: ldlt_z_det, ldlt_z_numerators, &
                                     ldlt_z_update, ldlt_z_update_jvp, ldlt_z_update_vjp
    implicit none
    type(native_ldlt_z_t) :: factor, compact
    type(csc_z_t) :: a
    type(fortsparse_status_t) :: status
    integer :: kind, trial, n, i, j, p(2)
    integer(i8) :: peak, retained
    logical :: h
    complex(dp) :: pair(2, 2), three(3, 3), z
    complex(dp), allocatable :: matrix(:, :), truth(:), b(:), x(:), values(:)
    integer, allocatable :: rows(:), columns(:)
    do kind = 1, 2
        h = kind == 1
        pair = (0.0_dp, 0.0_dp)
        pair(1, 2) = (0.0_dp, 1.0_dp)
        pair(2, 1) = pair(1, 2)
        if (h) pair(2, 1) = conjg(pair(1, 2))
        call check_matrix(pair, h)
        call check_matrix(1.0e-300_dp*pair, h)
        call check_matrix(1.0e300_dp*pair, h)
        three = (0.0_dp, 0.0_dp)
        three(1, 3) = (0.5_dp, 2.0_dp)
        three(2, 3) = (3.0_dp, 2.0_dp)
        three(2, 2) = (3.0_dp, 0.0_dp)
        three(3, 1:2) = three(1:2, 3)
        if (h) three(3, 1:2) = conjg(three(1:2, 3))
        call check_matrix(three, h)
        do trial = 1, 35
            n = 2 + mod(trial, 11)
            allocate (matrix(n, n))
            do j = 1, n
                do i = 1, j
                    z = cmplx(sin(real(7*i + 13*j + trial, dp)), &
                              cos(real(11*i + 3*j + trial, dp)), dp)
                    matrix(i, j) = z
                    matrix(j, i) = z
                    if (h) matrix(j, i) = conjg(z)
                end do
                if (h) matrix(j, j) = real(matrix(j, j), dp)
                matrix(j, j) = matrix(j, j) + 0.2_dp*(-1.0_dp)**j
            end do
            call check_matrix(matrix, h)
            deallocate (matrix)
        end do
        ! 1x1 interchange then same-pattern coefficient-dependent repivot.
        pair = reshape([cmplx(0.01_dp, 0.0_dp, dp), (1.0_dp, 0.25_dp), &
                        (1.0_dp, 0.25_dp), (10.0_dp, 0.0_dp)], [2, 2])
        if (h) pair(2, 1) = conjg(pair(1, 2))
        call from_dense(pair, a)
        call factor%factor(a, status, hermitian=h)
        call ok('interchange factor')
        call factor%permutation(p, status)
        call ok('interchange permutation')
        if (p(1) /= 2) error stop 'missing complex symmetric interchange'
        pair(1, 1) = (10.0_dp, 0.0_dp)
        pair(2, 2) = (0.01_dp, 0.0_dp)
        call from_dense(pair, a)
        call factor%refactor(a, status, hermitian=h)
        call ok('same pattern refactor')
        call factor%permutation(p, status)
        call ok('repivot permutation')
        if (p(1) /= 1) error stop 'complex refactor did not repivot'
        call factor%refactor(a, status, hermitian=.not. h)
        if (status%code /= FORTSPARSE_INVALID_MATRIX) error stop 'class changed'
        a%val = (1.0_dp, 0.0_dp)
        call factor%refactor(a, status)
        if (status%code /= FORTSPARSE_SINGULAR) error stop 'singular refactor'
        call from_dense(pair, a)
        call factor%refactor(a, status)
        call ok('recover retained pattern')
        call check_cap(h)
    end do
    ! The class distinction must be enforced, including real Hermitian diagonal.
    pair = reshape([(1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), &
                    (0.0_dp, 1.0_dp), (2.0_dp, 1.0_dp)], [2, 2])
    call from_dense(pair, a)
    call factor%factor(a, status, hermitian=.true.)
    if (status%code /= FORTSPARSE_INVALID_MATRIX) error stop 'non-Hermitian accepted'
    pair(2, 1) = conjg(pair(1, 2))
    pair(2, 2) = (2.0_dp, 0.0_dp)
    call from_dense(pair, a)
    call factor%factor(a, status, hermitian=.false.)
    if (status%code /= FORTSPARSE_INVALID_MATRIX) error stop 'not transpose symmetric'
    a%val(1) = cmplx(ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
    call factor%factor(a, status)
    if (status%code /= FORTSPARSE_INVALID_MATRIX) error stop 'nonfinite accepted'
    pair = (0.0_dp, 0.0_dp)
    pair(1, 1) = cmplx(0.75_dp*huge(1.0_dp), 0.75_dp*huge(1.0_dp), dp)
    pair(2, 2) = (1.0_dp, 0.0_dp)
    call from_dense(pair, a)
    call factor%factor(a, status, hermitian=.false.)
    if (status%code /= FORTSPARSE_INTERNAL_ERROR) error stop 'range guard missing'
    call leaf_checks()
    call factor%free()
    if (factor%retained_bytes() /= 0_i8) error stop 'complex free payload'
contains
    subroutine ok(label)
        character(*), intent(in) :: label
        if (status%code == FORTSPARSE_OK) return
        print *, label, status%code, trim(status%msg)
        error stop 'unexpected complex LDL status'
    end subroutine ok

    subroutine from_dense(dense, csc)
        complex(dp), intent(in) :: dense(:, :)
        type(csc_z_t), intent(out) :: csc
        integer :: n, i, j, k
        integer, allocatable :: r(:), c(:)
        complex(dp), allocatable :: v(:)
        n = size(dense, 1)
        allocate (r(n*n), c(n*n), v(n*n))
        k = 0
        do j = 1, n
            do i = 1, n
                if (abs(dense(i, j)) <= 0.0_dp) cycle
                k = k + 1
                r(k) = i
                c(k) = j
                v(k) = dense(i, j)
            end do
        end do
        call csc_from_triplet(n, n, r(:k), c(:k), v(:k), csc, status)
        call ok('dense-to-CSC fixture')
    end subroutine from_dense

    subroutine check_matrix(dense, hermitian)
        complex(dp), intent(in) :: dense(:, :)
        logical, intent(in) :: hermitian
        complex(dp), allocatable :: truth(:, :), rhs(:, :), solution(:, :), op(:, :)
        integer :: i, n, operation
        integer(i8) :: bytes
        n = size(dense, 1)
        allocate (truth(n, 2), rhs(n, 2), solution(n, 2), op(n, n))
        do i = 1, n
            truth(i, 1) = cmplx(cos(real(i, dp)), sin(real(i, dp)), dp)
            truth(i, 2) = cmplx(sin(real(2*i, dp)), cos(real(3*i, dp)), dp)
        end do
        call from_dense(dense, a)
        call factor%factor(a, status, hermitian=hermitian)
        call ok('complex factor')
        bytes = factor%retained_bytes()
        do operation = 0, 2
            op = dense
            if (operation == 1) op = transpose(dense)
            if (operation == 2) op = conjg(transpose(dense))
            rhs = matmul(op, truth)
            call factor%solve_block(rhs, solution, status, operation=operation)
            call ok('complex N/T/H block solve')
            if (maxval(abs(solution - truth)) > 3.0e-10_dp) &
                error stop 'complex known-solution oracle'
            if (maxval(abs(matmul(op, solution) - rhs)) &
                > 2.0e-12_dp*max(1.0_dp, maxval(abs(rhs)))) &
                error stop 'complex independent dense residual'
            solution(:, 1) = rhs(:, 1)
            call factor%solve_inplace(solution(:, 1), status, operation=operation)
            call ok('complex N/T/H inplace solve')
            if (maxval(abs(solution(:, 1) - truth(:, 1))) > 3.0e-10_dp) &
                error stop 'complex inplace oracle'
        end do
        if (factor%retained_bytes() /= bytes) error stop 'complex solve allocated'
        call factor%solve(rhs(:, 1), solution(:, 1), status, operation=3)
        if (status%code /= FORTSPARSE_INVALID_MATRIX) error stop 'bad operation'
    end subroutine check_matrix

    subroutine check_cap(hermitian)
        logical, intent(in) :: hermitian
        integer :: i, n
        n = 25
        allocate (matrix(n, n), truth(n), b(n), x(n))
        matrix = (0.0_dp, 0.0_dp)
        matrix(1, 1) = (100.0_dp, 0.0_dp)
        do i = 2, n
            matrix(i, i) = cmplx(2.0_dp*(-1.0_dp)**i, 0.0_dp, dp)
            matrix(1, i) = cmplx(0.2_dp, 0.1_dp*(-1.0_dp)**i, dp)
            matrix(i, 1) = matrix(1, i)
            if (hermitian) matrix(i, 1) = conjg(matrix(1, i))
        end do
        truth = [(cmplx(sin(real(i, dp)), cos(real(i, dp)), dp), i=1, n)]
        b = matmul(matrix, truth)
        call from_dense(matrix, a)
        call factor%factor(a, status, hermitian=hermitian)
        call ok('complex sparse growth')
        peak = factor%peak_bytes()
        call compact%factor(a, status, hermitian=hermitian, low_memory=.true.)
        call ok('complex low memory')
        if (compact%peak_bytes() > peak) error stop 'complex low memory peak'
        call compact%solve(b, x, status)
        call ok('complex low memory solve')
       if (maxval(abs(x - truth)) > 1.0e-12_dp) error stop 'complex low memory accuracy'
        call factor%factor(a, status, hermitian=hermitian, host_limit_bytes=peak)
        call ok('complex exact cap')
        retained = factor%retained_bytes()
        x = -b
        call factor%solve_inplace(x, status)
        call ok('complex inplace equation tangent')
        if (maxval(abs(x + truth)) > 1.0e-12_dp) error stop 'complex tangent oracle'
        if (factor%retained_bytes() /= retained) error stop 'complex inplace allocated'
        call factor%refactor(a, status, host_limit_bytes=peak)
        call ok('complex capped refactor')
        call factor%factor(a, status, hermitian=hermitian, host_limit_bytes=peak - 1_i8)
        if (status%code /= FORTSPARSE_RESOURCE_LIMIT) error stop 'complex growth cap'
        if (factor%retained_bytes() > peak - 1_i8) error stop 'complex cap overshot'
        if (factor%required_bytes() <= factor%available_bytes()) &
            error stop 'complex cap diagnostics'
        call factor%solve_inplace(x, status)
        if (status%code /= FORTSPARSE_NOT_FACTORED) error stop 'failed factors usable'
        call factor%factor(a, status, hermitian=hermitian, host_limit_bytes=1_i8)
        if (status%code /= FORTSPARSE_RESOURCE_LIMIT) error stop 'complex setup cap'
        deallocate (matrix, truth, b, x)
        n = 3000
        allocate (rows(n), columns(n), values(n), truth(n), b(n), x(n))
        do i = 1, n
            rows(i) = i
            columns(i) = i
            values(i) = cmplx((-1.0_dp)**i, 0.0_dp, dp)
            if (.not. hermitian) values(i) = values(i) + (0.0_dp, 0.5_dp)
            truth(i) = cmplx(sin(real(i, dp)), cos(real(i, dp)), dp)
        end do
        call csc_from_triplet(n, n, rows, columns, values, a, status)
        call ok('large complex diagonal CSC')
        call factor%factor(a, status, hermitian=hermitian, host_limit_bytes=4000000_i8)
        call ok('large complex diagonal factor')
        if (factor%factor_nnz() /= int(n, i8)) error stop 'large diagonal fill'
        b = csc_matvec(a, truth)
        call factor%solve(b, x, status)
        call ok('large complex diagonal solve')
        if (maxval(abs(x - truth)) > 1.0e-14_dp) error stop 'large diagonal accuracy'
        deallocate (rows, columns, values, truth, b, x)
    end subroutine check_cap

    subroutine leaf_checks()
        complex(dp) :: aa, bb, dd, cc, u, v, first, second, det, direct
        complex(dp) :: arguments(5), direction(5), tangent
        real(dp) :: parts(10), delta(10), bars(10), value(2), derivative(2), seed(2)
        integer :: i
        aa = (2.0_dp, 0.3_dp)
        bb = (0.4_dp, 0.7_dp)
        dd = (-0.3_dp, 0.2_dp)
        cc = (-1.0_dp, 0.5_dp)
        u = (0.75_dp, -0.2_dp)
        v = (1.25_dp, 0.1_dp)
        det = ldlt_z_det(aa, bb, dd, cc)
        if (abs(det - (aa*cc - bb*dd)) > 1.0e-14_dp) error stop 'determinant leaf'
        call ldlt_z_numerators(aa, bb, dd, cc, u, v, first, second)
        first = first/det
        second = second/det
        if (max(abs(first*aa + second*dd - u), abs(first*bb + second*cc - v)) &
            > 1.0e-14_dp) error stop 'block numerator row-equation oracle'
        arguments = [aa, bb, dd, cc, u]
        direction = [cc, dd, bb, aa, v]
        do i = 1, 5
            parts(2*i - 1) = real(arguments(i), dp)
            parts(2*i) = aimag(arguments(i))
            delta(2*i - 1) = real(direction(i), dp)
            delta(2*i) = aimag(direction(i))
        end do
        direct = arguments(1) - arguments(2)*arguments(3) - arguments(4)*arguments(5)
        if (abs(ldlt_z_update(aa, bb, dd, cc, u) - direct) > 1.0e-14_dp) &
            error stop 'complex rank-two update oracle'
        call ldlt_z_update_jvp(parts(1), parts(2), parts(3), parts(4), parts(5), &
                               parts(6), parts(7), parts(8), parts(9), parts(10), &
                 delta(1), delta(2), delta(3), delta(4), delta(5), delta(6), delta(7), &
        delta(8), delta(9), delta(10), value(1), value(2), derivative(1), derivative(2))
        tangent = direction(1) - direction(2)*arguments(3) - arguments(2)*direction(3) &
                  - direction(4)*arguments(5) - arguments(4)*direction(5)
        if (abs(cmplx(derivative(1), derivative(2), dp) - tangent) > 1.0e-14_dp) &
            error stop 'complex rank-two tangent oracle'
        seed = [0.7_dp, -0.4_dp]
        call ldlt_z_update_vjp(parts(1), parts(2), parts(3), parts(4), parts(5), &
                  parts(6), parts(7), parts(8), parts(9), parts(10), seed(1), seed(2), &
                      value(1), value(2), bars(1), bars(2), bars(3), bars(4), bars(5), &
                               bars(6), bars(7), bars(8), bars(9), bars(10))
       if (abs(dot_product(seed, derivative) - dot_product(bars, delta)) > 1.0e-13_dp) &
            error stop 'complex real-component adjoint pairing'
    end subroutine leaf_checks
end program test_fortsparse_native_ldlt_z
