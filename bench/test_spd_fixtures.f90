program test_spd_fixtures
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spd_fixtures, only: fixture_t, make_fixture, fixture_product
    implicit none
    type(fixture_t) :: a
    character(5), parameter :: families(3) = ["band ", "grid2", "grid3"]
    real(dp), allocatable :: ones(:), product(:), dense(:, :)
    integer :: family, col, p, j, k, held

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
            ! Independent dense Cholesky establishes positive definiteness.
            do j = 1, a%n
                dense(j, j) = dense(j, j) - sum(dense(j, :j - 1)**2)
                if (.not. (dense(j, j) > 0.0_dp)) error stop "not SPD"
                dense(j, j) = sqrt(dense(j, j))
                do k = j + 1, a%n
                    dense(k, j) = (dense(k, j) - &
                        dot_product(dense(k, :j - 1), &
                                    dense(j, :j - 1)))/dense(j, j)
                end do
            end do
            deallocate (ones, product, dense)
        end do
    end do
    print '(a)', 'PASS SPD fixture mass, symmetry and dense Cholesky oracles'
end program
