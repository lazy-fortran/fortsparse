module fortsparse_native_residual
    ! Original-CSC checks shared by native CPU families; factors are untrusted.
    use fortsparse_kinds, only: dp
    use fortsparse_status, only: fortsparse_status_t, status_set, FORTSPARSE_OK, FORTSPARSE_SINGULAR
    use fortsparse_spd_leaves, only: spd_update, spd_divide
    use fortsparse_lu_leaves, only: lu_z_update, lu_z_finite
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: native_original_norms, native_residual_error
    interface native_original_norms
        module procedure norms_real, norms_complex
    end interface
    interface native_residual_error
        module procedure residual_real, residual_complex
    end interface
contains
    subroutine norms_real(n,ptr,rows,values,scratch,matrix_norm,transpose_norm,status)
        integer, intent(in) :: n,ptr(:),rows(:)
        real(dp), intent(in) :: values(:)
        real(dp), intent(out) :: scratch(:)
        real(dp), intent(out) :: matrix_norm,transpose_norm
        type(fortsparse_status_t), intent(out) :: status
        integer :: j,p,i
        real(dp) :: magnitude,column_sum
        call status_set(status,FORTSPARSE_OK,'')
        scratch=0.0_dp
        matrix_norm=0.0_dp
        transpose_norm=0.0_dp
        do j=1,n
            column_sum=0.0_dp
            do p=ptr(j),ptr(j+1)-1
                i=rows(p)
                magnitude=abs(values(p))
                column_sum=spd_update(column_sum,-1.0_dp,magnitude)
                scratch(i)=spd_update(scratch(i),-1.0_dp,magnitude)
                if (.not.ieee_is_finite(column_sum) .or. .not.ieee_is_finite(scratch(i))) then
                    call status_set(status,FORTSPARSE_SINGULAR,'native refinement: original norm exceeds numerical range')
                    return
                end if
            end do
            transpose_norm=max(transpose_norm,column_sum)
        end do
        do i=1,n
            matrix_norm=max(matrix_norm,scratch(i))
        end do
    end subroutine norms_real

    subroutine residual_real(ptr,rows,values,rhs,x,operation,an,residual,error, &
            denominator,componentwise)
        ! Optional componentwise (Oettli-Prager) error: denominator is caller
        ! scratch of length n receiving (|op(A)||x|+|b|)_i; see componentwise_error.
        integer, intent(in) :: ptr(:),rows(:),operation
        real(dp), intent(in) :: values(:),rhs(:),x(:)
        real(dp), intent(out) :: residual(:)
        real(dp), intent(in) :: an
        real(dp), intent(out) :: error
        real(dp), intent(out), optional :: denominator(:)
        real(dp), intent(out), optional :: componentwise
        real(dp) :: rn,xn,bn
        real(dp) :: coefficient
        integer :: j,p,i
        logical :: both
        both=present(denominator) .and. present(componentwise)
        residual=rhs
        if (both) then
            componentwise=-1.0_dp
            denominator=abs(rhs)
        end if
        do j=1,size(x)
            do p=ptr(j),ptr(j+1)-1
                i=rows(p)
                coefficient=values(p)
                if (operation==0) then
                    residual(i)=spd_update(residual(i),coefficient,x(j))
                    if (both) denominator(i)=spd_update(denominator(i), &
                        -abs(coefficient),abs(x(j)))
                else
                    residual(j)=spd_update(residual(j),coefficient,x(i))
                    if (both) denominator(j)=spd_update(denominator(j), &
                        -abs(coefficient),abs(x(i)))
                end if
            end do
        end do
        rn=0.0_dp
        xn=0.0_dp
        bn=0.0_dp
        do i=1,size(x)
            if (.not.ieee_is_finite(x(i)) .or. .not.ieee_is_finite(residual(i))) then
                error=-1.0_dp
                return
            end if
            rn=max(rn,abs(residual(i)))
            xn=max(xn,abs(x(i)))
            bn=max(bn,abs(rhs(i)))
        end do
        error=backward_error(rn,an,xn,bn)
        if (.not.both) return
        componentwise=0.0_dp
        do i=1,size(x)
            if (.not.ieee_is_finite(denominator(i))) then
                componentwise=-1.0_dp
                return
            end if
            componentwise=max(componentwise,componentwise_ratio(abs(residual(i)), &
                denominator(i),size(x)))
        end do
    end subroutine residual_real

    subroutine norms_complex(n,ptr,rows,values,scratch,matrix_norm,transpose_norm,status)
        integer, intent(in) :: n,ptr(:),rows(:)
        complex(dp), intent(in) :: values(:)
        complex(dp), intent(out) :: scratch(:)
        real(dp), intent(out) :: matrix_norm,transpose_norm
        type(fortsparse_status_t), intent(out) :: status
        integer :: j,p,i
        real(dp) :: magnitude,column_sum,row_sum
        call status_set(status,FORTSPARSE_OK,'')
        scratch=(0.0_dp,0.0_dp)
        matrix_norm=0.0_dp
        transpose_norm=0.0_dp
        do j=1,n
            column_sum=0.0_dp
            do p=ptr(j),ptr(j+1)-1
                i=rows(p)
                magnitude=abs(values(p))
                column_sum=spd_update(column_sum,-1.0_dp,magnitude)
                row_sum=spd_update(real(scratch(i),dp),-1.0_dp,magnitude)
                if (.not.ieee_is_finite(column_sum) .or. .not.ieee_is_finite(row_sum)) then
                    call status_set(status,FORTSPARSE_SINGULAR,'native refinement: original norm exceeds numerical range')
                    return
                end if
                scratch(i)=cmplx(row_sum,0.0_dp,dp)
            end do
            transpose_norm=max(transpose_norm,column_sum)
        end do
        do i=1,n
            matrix_norm=max(matrix_norm,real(scratch(i),dp))
        end do
    end subroutine norms_complex

    subroutine residual_complex(ptr,rows,values,rhs,x,operation,an,residual,error, &
            denominator,componentwise)
        ! Optional componentwise error as in residual_real, with moduli.
        integer, intent(in) :: ptr(:),rows(:),operation
        complex(dp), intent(in) :: values(:),rhs(:),x(:)
        complex(dp), intent(out) :: residual(:)
        real(dp), intent(in) :: an
        real(dp), intent(out) :: error
        real(dp), intent(out), optional :: denominator(:)
        real(dp), intent(out), optional :: componentwise
        real(dp) :: rn,xn,bn
        complex(dp) :: coefficient
        integer :: j,p,i
        logical :: both
        both=present(denominator) .and. present(componentwise)
        residual=rhs
        if (both) then
            componentwise=-1.0_dp
            denominator=abs(rhs)
        end if
        do j=1,size(x)
            do p=ptr(j),ptr(j+1)-1
                i=rows(p)
                coefficient=values(p)
                if (operation==2) coefficient=conjg(coefficient)
                if (operation==0) then
                    residual(i)=lu_z_update(residual(i),coefficient,x(j))
                    if (both) denominator(i)=spd_update(denominator(i), &
                        -abs(coefficient),abs(x(j)))
                else
                    residual(j)=lu_z_update(residual(j),coefficient,x(i))
                    if (both) denominator(j)=spd_update(denominator(j), &
                        -abs(coefficient),abs(x(i)))
                end if
            end do
        end do
        rn=0.0_dp
        xn=0.0_dp
        bn=0.0_dp
        do i=1,size(x)
            if (.not.lu_z_finite(x(i)) .or. .not.lu_z_finite(residual(i))) then
                error=-1.0_dp
                return
            end if
            rn=max(rn,abs(residual(i)))
            xn=max(xn,abs(x(i)))
            bn=max(bn,abs(rhs(i)))
        end do
        error=backward_error(rn,an,xn,bn)
        if (.not.both) return
        componentwise=0.0_dp
        do i=1,size(x)
            if (.not.ieee_is_finite(denominator(i))) then
                componentwise=-1.0_dp
                return
            end if
            componentwise=max(componentwise,componentwise_ratio(abs(residual(i)), &
                denominator(i),size(x)))
        end do
    end subroutine residual_complex

    ! One row of the componentwise (Oettli-Prager) backward error
    ! max_i |r_i|/(|op(A)||x|+|b|)_i with the LAPACK xGERFS safe denominator:
    ! d<=safe2 uses (|r|+safe1)/(d+safe1), safe1=(n+1)*tiny, safe2=safe1/eps.
    ! An exactly zero residual contributes zero (LAPACK reports one for the
    ! 0/0 row of a zero RHS). Callers mark a nonfinite denominator with -1.
    real(dp) function componentwise_ratio(r,d,n) result(ratio)
        real(dp), intent(in) :: r,d
        integer, intent(in) :: n
        real(dp) :: safe1,safe2
        ratio=0.0_dp
        if (r==0.0_dp) return
        safe1=real(n+1,dp)*tiny(1.0_dp)
        safe2=safe1/epsilon(1.0_dp)
        if (d>safe2) then
            ratio=spd_divide(r,d)
        else
            ratio=spd_divide(r+safe1,d+safe1)
        end if
    end function componentwise_ratio

    ! Evaluate r/(a*x+b) by power-of-two scheduling of generated real leaves.
    ! This avoids overflow/underflow of the norm product without accepting a
    ! nonfinite residual or solution. It is a backward-error estimate, no bound.
    real(dp) function backward_error(r,a,x,b) result(error)
        real(dp), intent(in) :: r,a,x,b
        real(dp) :: denominator,product,numerator
        integer :: ep,eb,e
        error=-1.0_dp
        if (.not.ieee_is_finite(r) .or. .not.ieee_is_finite(a) .or. &
            .not.ieee_is_finite(x) .or. .not.ieee_is_finite(b)) return
        if (r==0.0_dp) then
            error=0.0_dp
            return
        end if
        ep=-huge(0)/2
        eb=ep
        if (a>0.0_dp .and. x>0.0_dp) ep=exponent(a)+exponent(x)
        if (b>0.0_dp) eb=exponent(b)
        if (ep==eb .and. ep==-huge(0)/2) return
        e=max(ep,eb)
        denominator=0.0_dp
        if (a>0.0_dp .and. x>0.0_dp) then
            product=spd_update(0.0_dp,-fraction(a),fraction(x))
            denominator=spd_update(denominator,-1.0_dp,scale(product,ep-e))
        end if
        if (b>0.0_dp) denominator=spd_update(denominator,-1.0_dp,scale(fraction(b),eb-e))
        numerator=scale(fraction(r),exponent(r)-e)
        if (.not.ieee_is_finite(numerator)) return
        error=spd_divide(numerator,denominator)
        if (.not.ieee_is_finite(error)) error=-1.0_dp
    end function backward_error
end module fortsparse_native_residual
