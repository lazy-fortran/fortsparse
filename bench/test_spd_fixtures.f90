program test_spd_fixtures
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spd_fixtures, only: fixture_t, make_fixture, fixture_product, &
        fixture_spd_update, fixture_known_rhs
    implicit none
    type(fixture_t) :: a, changed
    character(5), parameter :: families(3) = ["band ", "grid2", "grid3"]
    real(dp), allocatable :: ones(:), product(:), dense(:, :), rhs(:, :)
    real(dp), allocatable :: reference(:, :), addition(:)
    integer :: family, col, p, j, held

    do family = 1, 3
        do held = 0, 1
            call make_fixture(trim(families(family)), 1, held == 1, a)
            allocate (ones(a%n), product(a%n), dense(a%n, a%n))
            ones = 1.0_dp
            call fixture_product(a, ones, product)
            if (maxval(abs(product - 0.25_dp)) > 1.0e-14_dp) &
                error stop "mass shift / row sum"
            dense = 0.0_dp
            do col = 1, a%n
                do p = a%ptr(col), a%ptr(col + 1) - 1
                    if (p > a%ptr(col)) then
                        if (a%row(p) <= a%row(p - 1)) &
                            error stop "unsorted or duplicate row"
                    end if
                    dense(a%row(p), col) = a%val(p)
                end do
            end do
            if (maxval(abs(dense - transpose(dense))) > 1.0e-14_dp) &
                error stop "asymmetric fixture"
            call check_spd(dense)
            call fixture_spd_update(a, changed)
            allocate (rhs(a%n, 2), reference(a%n, 2), addition(a%n))
            call fixture_product(changed, ones, product)
            addition = product - 0.25_dp
            if (minval(addition) <= 0.0_dp) error stop "update is not positive"
            if (maxval(addition) - minval(addition) < 0.001_dp) &
                error stop "update is uniform"
            if (any(a%ptr /= changed%ptr) .or. any(a%row /= changed%row)) &
                error stop "update changed the graph"
            do j = 1, a%n
                dense(j, j) = dense(j, j) + addition(j)
            end do
            call fixture_known_rhs(changed, 2, 3, rhs, reference)
            if (maxval(abs(rhs - matmul(dense, reference))) > 1.0e-12_dp) &
                error stop "repeated RHS does not match independent dense product"
            call check_spd(dense)
            deallocate (ones, product, dense, rhs, reference, addition)
        end do
    end do
    print '(a)', 'PASS SPD fixture mass, symmetry and dense Cholesky oracles'

contains

    subroutine check_spd(matrix)
        real(dp), intent(in) :: matrix(:, :)
        real(dp) :: factor(size(matrix, 1), size(matrix, 2))
        integer :: j, k

        ! Independent dense Cholesky checks original and updated fixtures.
        factor = matrix
        do j = 1, size(matrix, 1)
            factor(j, j) = factor(j, j) - sum(factor(j, :j - 1)**2)
            if (.not. (factor(j, j) > 0.0_dp)) error stop "not SPD"
            factor(j, j) = sqrt(factor(j, j))
            do k = j + 1, size(matrix, 1)
                factor(k, j) = (factor(k, j) - &
                    dot_product(factor(k, :j - 1), &
                                factor(j, :j - 1)))/factor(j, j)
            end do
        end do
    end subroutine check_spd
end program
