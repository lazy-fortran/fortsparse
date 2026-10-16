    ! Power-of-two normalization is an execution schedule for the generated
    ! Cartesian quotient. It avoids forming |b|**2 at the original scale.
    pure elemental function lu_z_divide(a,b) result(value)
        complex(dp), intent(in) :: a,b
        complex(dp) :: value
        real(dp) :: ar,ai,br,bi,vr,vi,am,bm
        integer :: ea,eb
        am=lu_z_norm(a)
        bm=lu_z_norm(b)
        if (am==0.0_dp .and. bm>0.0_dp) then
            value=(0.0_dp,0.0_dp)
            return
        end if
        ea=exponent(am)
        eb=exponent(bm)
        ar=scale(real(a,dp),-ea)
        ai=scale(aimag(a),-ea)
        br=scale(real(b,dp),-eb)
        bi=scale(aimag(b),-eb)
        call lu_z_quotient_value(ar,ai,br,bi,vr,vi)
        value=cmplx(scale(vr,ea-eb),scale(vi,ea-eb),dp)
    end function lu_z_divide

    pure elemental function lu_z_update(a,b,c) result(value)
        complex(dp), intent(in) :: a,b,c
        complex(dp) :: value
        real(dp) :: vr,vi
        call lu_z_update_value(real(a,dp),aimag(a),real(b,dp),aimag(b),real(c,dp),aimag(c),vr,vi)
        value=cmplx(vr,vi,dp)
    end function lu_z_update

    pure elemental real(dp) function lu_z_norm(a) result(value)
        complex(dp), intent(in) :: a
        ! Component infinity norm is overflow-safe for finite components and
        ! is the declared complex threshold-pivot norm.
        value=max(abs(real(a,dp)),abs(aimag(a)))
    end function lu_z_norm

    pure elemental logical function lu_z_finite(a) result(value)
        complex(dp), intent(in) :: a
        value=ieee_is_finite(real(a,dp)).and.ieee_is_finite(aimag(a))
    end function lu_z_finite
