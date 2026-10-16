    pure elemental function hpd_root(a) result(v)
        real(dp), intent(in) :: a
        real(dp) :: v
        call hpd_root_value(a, v)
    end function hpd_root

    pure elemental function hpd_divide(a, d) result(v)
        complex(dp), intent(in) :: a
        real(dp), intent(in) :: d
        complex(dp) :: v
        real(dp) :: vr, vi
        call hpd_divide_value(real(a, dp), aimag(a), d, vr, vi)
        v = cmplx(vr, vi, dp)
    end function hpd_divide

    pure elemental function hpd_update(a, b, c) result(v)
        complex(dp), intent(in) :: a, b, c
        complex(dp) :: v
        real(dp) :: vr, vi
        call hpd_update_value(real(a, dp), aimag(a), real(b, dp), aimag(b), &
                              real(c, dp), aimag(c), vr, vi)
        v = cmplx(vr, vi, dp)
    end function hpd_update

    pure elemental function hpd_conj_update(a, b, c) result(v)
        complex(dp), intent(in) :: a, b, c
        complex(dp) :: v
        real(dp) :: vr, vi
        call hpd_conj_update_value(real(a, dp), aimag(a), real(b, dp), aimag(b), &
                                   real(c, dp), aimag(c), vr, vi)
        v = cmplx(vr, vi, dp)
    end function hpd_conj_update

    pure elemental function hpd_norm_update(a, b) result(v)
        real(dp), intent(in) :: a
        complex(dp), intent(in) :: b
        real(dp) :: v
        call hpd_norm_update_value(a, real(b, dp), aimag(b), v)
    end function hpd_norm_update
