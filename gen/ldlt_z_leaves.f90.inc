    pure elemental function ldlt_z_det(a, b, d, c) result(value)
        complex(dp), intent(in) :: a, b, d, c
        complex(dp) :: value
        real(dp) :: vr, vi
        call ldlt_z_det_value(real(a, dp), aimag(a), real(b, dp), aimag(b), &
                              real(d, dp), aimag(d), real(c, dp), aimag(c), vr, vi)
        value = cmplx(vr, vi, dp)
    end function ldlt_z_det

    pure elemental subroutine ldlt_z_numerators(a, b, d, c, u, v, first, second)
        complex(dp), intent(in) :: a, b, d, c, u, v
        complex(dp), intent(out) :: first, second
        real(dp) :: v1, v2, v3, v4
        call ldlt_z_numerators_value(real(a, dp), aimag(a), real(b, dp), aimag(b), &
                  real(d, dp), aimag(d), real(c, dp), aimag(c), real(u, dp), aimag(u), &
                                     real(v, dp), aimag(v), v1, v2, v3, v4)
        first = cmplx(v1, v2, dp)
        second = cmplx(v3, v4, dp)
    end subroutine ldlt_z_numerators

    pure elemental function ldlt_z_update(a, b, c, d, e) result(value)
        complex(dp), intent(in) :: a, b, c, d, e
        complex(dp) :: value
        real(dp) :: vr, vi
        call ldlt_z_update_value(real(a, dp), aimag(a), real(b, dp), aimag(b), &
            real(c, dp), aimag(c), real(d, dp), aimag(d), real(e, dp), aimag(e), vr, vi)
        value = cmplx(vr, vi, dp)
    end function ldlt_z_update
