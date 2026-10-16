program test_native_family_fixtures
    use native_family_fixtures, only: dp, family_fixture_t, make_family_fixture, &
                                      known_rhs, fixture_action, fixture_backward_error
    implicit none
    character(16), parameter :: names(5) = [character(16) :: &
                 'neo2-mini5', 'drift-real', 'drift-complex', 'hpd-chain', 'ldlt-chain']
    type(family_fixture_t) :: a
    complex(dp), allocatable :: dense(:, :), x(:, :), b(:, :), y(:)
    character(1024) :: mini5_path
    integer :: family, variant, op, j, p, i
    real(dp) :: row_sum
    call get_command_argument(1, mini5_path)
    do family = 1, size(names)
        do variant = 0, 1
            call make_family_fixture(trim(names(family)), 3, variant == 1, a, &
                                     trim(mini5_path))
            allocate (dense(a%n, a%n), x(a%n, 1), b(a%n, 1), y(a%n))
            dense = cmplx(0.0_dp, 0.0_dp, dp)
            do j = 1, a%n
                do p = a%ptr(j), a%ptr(j + 1) - 1
                    i = a%row(p)
                    dense(i, j) = dense(i, j) + a%val(p)
                    if (p > a%ptr(j)) &
                        call require(a%row(p) > a%row(p - 1), 'canonical row order')
                end do
            end do
            do op = 0, 2
                call known_rhs(a, trim(names(family)), 1, variant == 1, op, b, x)
                select case (op)
                case (0)
                    y = matmul(dense, x(:, 1))
                case (1)
                    y = matmul(transpose(dense), x(:, 1))
                case (2)
                    y = matmul(conjg(transpose(dense)), x(:, 1))
                end select
               call require(maxval(abs(y - b(:, 1))) < 2.e-13_dp, 'dense action oracle')
                call require(fixture_backward_error(a, x(:, 1), b(:, 1), op) &
                             < 1.e-15_dp, 'original-matrix residual')
                if (family == 1 .and. variant == 0 .and. op == 0) &
                    call require(maxval(abs(b(:, 1) - 1.0_dp)) < 3.e-14_dp, &
                                 'actual mini5 all-ones RHS and rational oracle')
            end do
            if (a%mtype == 4) then
                call require(maxval(abs(dense - conjg(transpose(dense)))) &
                             < 1.e-15_dp, 'HPD Hermitian identity')
                do i = 1, a%n
                    row_sum = sum(abs(dense(i, :))) - abs(dense(i, i))
                    call require(real(dense(i, i), dp) > row_sum, &
                                 'strict positive Hermitian diagonal dominance')
                end do
            end if
            if (a%mtype == -2) then
                call require(maxval(abs(dense - transpose(dense))) < 1.e-15_dp, &
                             'indefinite transpose symmetry')
                do i = 1, a%n
                   call require(abs(dense(i, i)) == 0.0_dp, 'adversarial zero diagonal')
                end do
                do i = 1, a%n - 1, 2
                    call require(abs(dense(i, i + 1)) > 1.0_dp, &
                                 'nonzero paired pivots and tridiagonal determinant')
                end do
            end if
            deallocate (dense, x, b, y)
        end do
    end do
    call larger_linear_fixtures()
    print '(a)', 'native family fixture dense/structural/original-data oracles: PASS'
contains
    subroutine larger_linear_fixtures()
        type(family_fixture_t) :: large
        complex(dp), allocatable :: input(:),actual(:),expected(:)
        complex(dp) :: edge
        integer :: k,kind
        do kind = 1,3
            select case(kind)
            case(1)
                call make_family_fixture('hpd-chain',8192,.false.,large,'')
            case(2)
                call make_family_fixture('ldlt-chain',8192,.false.,large,'')
            case(3)
                call make_family_fixture('drift-real',128,.false.,large,'')
            end select
            call require(large%n == 16384,'largest frozen linear-size layout')
            if (kind == 3) then
                call require(size(large%val) == 5*large%n-4*128,'2D stencil count')
                cycle
            end if
            call require(size(large%val) == 3*large%n-2,'chain count')
            allocate(input(large%n),actual(large%n),expected(large%n))
            input = cmplx(1.0_dp,0.25_dp,dp)
            expected = cmplx(0.0_dp,0.0_dp,dp)
            if (kind == 1) expected = 3.21_dp*input
            do k = 1,large%n-1
                edge = cmplx(-0.4_dp,0.13_dp,dp)
                if (kind == 2) then
                    edge = cmplx(0.125_dp+0.21_dp/8,0.0_dp,dp)
                    if (mod(k,2) == 1) edge = cmplx(2.21_dp,0.0_dp,dp)
                end if
                expected(k) = expected(k)+edge*input(k+1)
                expected(k+1) = expected(k+1)+conjg(edge)*input(k)
            end do
            call fixture_action(large,input,actual,0)
            call require(maxval(abs(actual-expected)) < 2.e-15_dp, &
                         'independent large tridiagonal action')
            deallocate(input,actual,expected)
        end do
    end subroutine larger_linear_fixtures

    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (.not. condition) error stop label
    end subroutine
end program test_native_family_fixtures
