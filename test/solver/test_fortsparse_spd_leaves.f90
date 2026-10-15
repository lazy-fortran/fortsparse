program test_fortsparse_spd_leaves
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortsparse_spd_leaves, only: spd_root, spd_divide, spd_update, &
        spd_root_jvp, spd_root_vjp, spd_divide_jvp, spd_divide_vjp, &
        spd_update_jvp, spd_update_vjp, spd_panel_factor, spd_panel_solve, &
        spd_schur_update, spd_scatter_update, spd_gather_update, spd_scatter_mapped_update
    implicit none
    real(dp), parameter :: h = 1.0e-5_dp
    real(dp) :: a, b, c, da, db, dc, value, tangent, abar, bbar, cbar, fd
    real(dp) :: l(4,4), matrix(4,4), reconstruction(4,4), panel(3,4), rhs(3,4)
    real(dp) :: left(3,2), right(4,2), schur(3,4), expected(3,4)
    real(dp) :: sparse_x(10), sparse_expected(10), sparse_values(4), dot_expected
    integer, parameter :: sparse_indices(4) = [2, 5, 7, 9]
    integer :: sparse_map(10)
    integer :: i, j, info

    ! Exact independent values and nontrivial signed derivative directions.
    call close(spd_root(9.0_dp), 3.0_dp, 1.0e-15_dp)
    call close(spd_divide(7.0_dp, 4.0_dp), 1.75_dp, 1.0e-15_dp)
    call close(spd_update(7.0_dp, 3.0_dp, -2.0_dp), 13.0_dp, 1.0e-15_dp)
    do i = 1, 12
        a = real(i, dp)/3.0_dp
        b = -0.5_dp-real(i, dp)/4.0_dp
        c = 0.25_dp*real(i, dp)
        da = -0.75_dp
        db = 0.375_dp
        dc = 0.5_dp
        call spd_root_jvp(a, da, value, tangent)
        fd = (sqrt(a+h*da)-sqrt(a-h*da))/(2*h)
        call close(tangent, fd, 2.0e-9_dp)
        call spd_root_vjp(a, 1.25_dp, value, abar)
        call close(1.25_dp*tangent, abar*da, 1.0e-14_dp)
        call spd_divide_jvp(a, b, da, db, value, tangent)
        fd = ((a+h*da)/(b+h*db)-(a-h*da)/(b-h*db))/(2*h)
        call close(tangent, fd, 2.0e-9_dp)
        call spd_divide_vjp(a, b, 1.25_dp, value, abar, bbar)
        call close(1.25_dp*tangent, abar*da+bbar*db, 1.0e-14_dp)
        call spd_update_jvp(a, b, c, da, db, dc, value, tangent)
        fd = ((a+h*da)-(b+h*db)*(c+h*dc) &
             -((a-h*da)-(b-h*db)*(c-h*dc)))/(2*h)
        call close(tangent, fd, 2.0e-9_dp)
        call spd_update_vjp(a, b, c, 1.25_dp, value, abar, bbar, cbar)
        call close(1.25_dp*tangent, abar*da+bbar*db+cbar*dc, 1.0e-14_dp)
    end do

    l = 0.0_dp
    do j = 1, 4
        l(j,j) = real(j+1, dp)
        do i = j+1, 4
            l(i,j) = real(i-j, dp)/4.0_dp
        end do
    end do
    matrix = matmul(l, transpose(l))
    call spd_panel_factor(matrix, 4, 4, info)
    if (info /= 0) error stop 'SPD panel factorization failed'
    reconstruction = 0.0_dp
    do j = 1, 4
        reconstruction(j:4,j) = matrix(j:4,j)
    end do
    if (maxval(abs(reconstruction-l)) > 2.0e-14_dp) &
        error stop 'panel factor differs from independent exact factor'
    do j = 1, 4
        do i = 1, 3
            rhs(i,j) = real(2*i-j, dp)/8.0_dp
        end do
    end do
    panel = rhs
    call spd_panel_solve(l, 4, 4, panel, 3, 3)
    if (maxval(abs(matmul(panel, transpose(l))-rhs)) > 2.0e-14_dp) &
        error stop 'triangular panel solve residual'
    left = rhs(:,1:2)
    right = transpose(rhs(1:2,:))
    schur = rhs
    expected = rhs-matmul(left, transpose(right))
    call spd_schur_update(schur, 3, left, 3, right, 4, 3, 4, 2)
    if (maxval(abs(schur-expected)) > 2.0e-14_dp) &
        error stop 'Schur update differs from independent matmul'
    sparse_values = [0.2_dp, -0.3_dp, 0.5_dp, 1.2_dp]
    do i = 1, 10
        sparse_x(i) = real(i, dp)/8.0_dp
    end do
    sparse_expected = sparse_x
    sparse_expected(sparse_indices) = sparse_expected(sparse_indices) &
        -1.25_dp*sparse_values
    call spd_scatter_update(sparse_x, sparse_indices, sparse_values, 4, 1.25_dp)
    if (maxval(abs(sparse_x-sparse_expected)) > 2.0e-14_dp) &
        error stop 'sparse scatter differs from independent vector oracle'
    value = 3.0_dp
    dot_expected = value-dot_product(sparse_values, sparse_x(sparse_indices))
    call spd_gather_update(value, sparse_indices, sparse_values, sparse_x, 4)
    call close(value, dot_expected, 2.0e-14_dp)
    call spd_scatter_update(sparse_x, sparse_indices, sparse_values, 0, 1.25_dp)
    call spd_gather_update(value, sparse_indices, sparse_values, sparse_x, 0)
    call close(value, dot_expected, 2.0e-14_dp)
    if (maxval(abs(sparse_x-sparse_expected)) > 2.0e-14_dp) &
        error stop 'zero-length sparse scatter changed its argument'
    sparse_map = 0
    sparse_map(2) = 2
    sparse_map(5) = 5
    sparse_map(9) = 9
    sparse_expected([2,5,9]) = sparse_expected([2,5,9]) &
        -0.5_dp*sparse_values([1,2,4])
    call spd_scatter_mapped_update(sparse_x, sparse_indices, sparse_map, sparse_values, 4, 0.5_dp)
    if (maxval(abs(sparse_x-sparse_expected)) > 2.0e-14_dp) &
        error stop 'mapped sparse scatter mishandled a structurally absent row'
    matrix = 0.0_dp
    call spd_panel_factor(matrix, 4, 4, info)
    if (info /= 1) error stop 'missing nonpositive pivot status'
    print '(a)', 'SPD generated values, JVP/VJP and panel schedules pass'
contains
    subroutine close(actual, expected_value, tolerance)
        real(dp), intent(in) :: actual, expected_value, tolerance
        if (.not. ieee_is_finite(actual)) error stop 'generated leaf is not finite'
        if (abs(actual-expected_value) > tolerance*max(1.0_dp, abs(expected_value))) &
            error stop 'generated leaf differs from independent oracle'
    end subroutine close
end program test_fortsparse_spd_leaves
