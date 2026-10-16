program test_fortsparse_hpd_leaves
    use fortsparse_kinds, only: dp
    use fortsparse_hpd_leaves, only: hpd_root_value, hpd_divide_value, &
                       hpd_update_value, hpd_conj_update_value, hpd_norm_update_value, &
                    hpd_root_jvp, hpd_divide_jvp, hpd_update_jvp, hpd_conj_update_jvp, &
                    hpd_norm_update_jvp, hpd_root_vjp, hpd_divide_vjp, hpd_update_vjp, &
                                     hpd_conj_update_vjp, hpd_norm_update_vjp
    implicit none
    real(dp) :: x(6), dx(6), y(2), dy(2), yp(2), ym(2), bar(6), seeds(2)
    real(dp), parameter :: step = 1.e-6_dp
    complex(dp) :: a, b, c, expected
    integer :: kind

    x = [2.3_dp, -0.7_dp, 1.2_dp, 0.4_dp, -0.3_dp, 0.8_dp]
    dx = [-0.3_dp, 0.2_dp, 0.4_dp, -0.1_dp, 0.6_dp, -0.2_dp]
    seeds = [0.6_dp, -0.4_dp]
    do kind = 1, 5
        call evaluate(kind, x, y)
        call jet(kind, x, dx, yp, dy)
        call require(maxval(abs(y - yp)) < 1.e-15_dp, 'primal/JVP consistency')
        call evaluate(kind, x + step*dx, yp)
        call evaluate(kind, x - step*dx, ym)
        call require(maxval(abs(dy - (yp - ym)/(2*step))) < 8.e-10_dp, &
                     'independent finite-difference JVP')
        call pullback(kind, x, seeds, yp, bar)
        call require(maxval(abs(y - yp)) < 1.e-15_dp, 'primal/VJP consistency')
        call require(abs(sum(seeds*dy) - sum(bar*dx)) < 3.e-15_dp, &
                     'real component pairing JVP/VJP identity')
    end do
    a = cmplx(x(1), x(2), dp)
    b = cmplx(x(3), x(4), dp)
    c = cmplx(x(5), x(6), dp)
    call evaluate(3, x, y)
    expected = a - b*c
    call require(abs(cmplx(y(1), y(2), dp) - expected) < 5.e-16_dp, &
                 'independent complex product oracle')
    call evaluate(4, x, y)
    expected = a - b*conjg(c)
    call require(abs(cmplx(y(1), y(2), dp) - expected) < 5.e-16_dp, &
                 'independent conjugated product oracle')
    call evaluate(5, x, y)
    call require(abs(y(1) - (x(1) - abs(cmplx(x(2), x(3), dp))**2)) &
                 < 1.e-15_dp, 'independent squared-norm oracle')
    print '(a)', 'HPD generated leaves values/JVP/VJP: PASS'
contains
    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (.not. condition) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
    end subroutine require

    subroutine evaluate(kind, x, y)
        integer, intent(in) :: kind
        real(dp), intent(in) :: x(6)
        real(dp), intent(out) :: y(2)
        y = 0.0_dp
        select case (kind)
        case (1)
            call hpd_root_value(x(1), y(1))
        case (2)
            call hpd_divide_value(x(1), x(2), x(3), y(1), y(2))
        case (3)
            call hpd_update_value(x(1), x(2), x(3), x(4), x(5), x(6), y(1), y(2))
        case (4)
            call hpd_conj_update_value(x(1), x(2), x(3), x(4), x(5), x(6), y(1), y(2))
        case (5)
            call hpd_norm_update_value(x(1), x(2), x(3), y(1))
        end select
    end subroutine evaluate

    subroutine jet(kind, x, dx, y, dy)
        integer, intent(in) :: kind
        real(dp), intent(in) :: x(6), dx(6)
        real(dp), intent(out) :: y(2), dy(2)
        y = 0.0_dp
        dy = 0.0_dp
        select case (kind)
        case (1)
            call hpd_root_jvp(x(1), dx(1), y(1), dy(1))
        case (2)
            call hpd_divide_jvp(x(1), x(2), x(3), dx(1), dx(2), dx(3), &
                                y(1), y(2), dy(1), dy(2))
        case (3)
            call hpd_update_jvp(x(1), x(2), x(3), x(4), x(5), x(6), &
                     dx(1), dx(2), dx(3), dx(4), dx(5), dx(6), y(1), y(2), dy(1), dy(2))
        case (4)
            call hpd_conj_update_jvp(x(1), x(2), x(3), x(4), x(5), x(6), &
                     dx(1), dx(2), dx(3), dx(4), dx(5), dx(6), y(1), y(2), dy(1), dy(2))
        case (5)
            call hpd_norm_update_jvp(x(1), x(2), x(3), dx(1), dx(2), dx(3), y(1), dy(1))
        end select
    end subroutine jet

    subroutine pullback(kind, x, s, y, bar)
        integer, intent(in) :: kind
        real(dp), intent(in) :: x(6), s(2)
        real(dp), intent(out) :: y(2), bar(6)
        y = 0.0_dp
        bar = 0.0_dp
        select case (kind)
        case (1)
            call hpd_root_vjp(x(1), s(1), y(1), bar(1))
        case (2)
            call hpd_divide_vjp(x(1), x(2), x(3), s(1), s(2), &
                                y(1), y(2), bar(1), bar(2), bar(3))
        case (3)
            call hpd_update_vjp(x(1), x(2), x(3), x(4), x(5), x(6), s(1), s(2), &
                             y(1), y(2), bar(1), bar(2), bar(3), bar(4), bar(5), bar(6))
        case (4)
            call hpd_conj_update_vjp(x(1), x(2), x(3), x(4), x(5), x(6), s(1), s(2), &
                             y(1), y(2), bar(1), bar(2), bar(3), bar(4), bar(5), bar(6))
        case (5)
          call hpd_norm_update_vjp(x(1), x(2), x(3), s(1), y(1), bar(1), bar(2), bar(3))
        end select
    end subroutine pullback
end program test_fortsparse_hpd_leaves
