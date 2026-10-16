    ! Ordered Kahan stages. FortSym emits each stage independently; composing
    ! the stages over exact reals would erase the rounding-error state.
    pure subroutine lu_r_compensated_add(value,correction,increment)
        real(dp), intent(inout) :: value,correction
        real(dp), intent(in) :: increment
        real(dp) :: x1,x2,vr,adjusted,total,lost
        x1=increment; x2=correction
@DIFFERENCE@
        adjusted=vr
        x1=value; x2=adjusted
@SUM@
        total=vr
        x1=total; x2=value
@DIFFERENCE@
        lost=vr
        x1=lost; x2=adjusted
@DIFFERENCE@
        correction=vr
        value=total
    end subroutine lu_r_compensated_add

    pure elemental function lu_r_corrected_value(value,correction) result(answer)
        real(dp), intent(in) :: value,correction
        real(dp) :: answer,x1,x2,vr
        x1=value; x2=correction
@DIFFERENCE@
        answer=vr
    end function lu_r_corrected_value

    pure elemental function lu_z_corrected_value(value,correction) result(answer)
        complex(dp), intent(in) :: value,correction
        complex(dp) :: answer
        answer=cmplx(lu_r_corrected_value(real(value,dp),real(correction,dp)), &
            lu_r_corrected_value(aimag(value),aimag(correction)),dp)
    end function lu_z_corrected_value

    pure subroutine lu_r_compensated_scatter(x,correction,indices,values,n,v)
        integer, intent(in) :: n,indices(*)
        real(dp), intent(inout) :: x(*),correction(*)
        real(dp), intent(in) :: values(*),v
        real(dp) :: x1,x2,vr
        integer :: i,row
        do i=1,n
            row=indices(i)
            x1=values(i); x2=v
@NEGATIVE_PRODUCT@
            call lu_r_compensated_add(x(row),correction(row),vr)
        end do
    end subroutine lu_r_compensated_scatter

    pure subroutine lu_r_compensated_gather(value,indices,values,x,n)
        integer, intent(in) :: n,indices(*)
        real(dp), intent(inout) :: value
        real(dp), intent(in) :: values(*),x(*)
        real(dp) :: correction,x1,x2,vr
        integer :: i
        correction=0.0_dp
        do i=1,n
            x1=values(i); x2=x(indices(i))
@NEGATIVE_PRODUCT@
            call lu_r_compensated_add(value,correction,vr)
        end do
        value=lu_r_corrected_value(value,correction)
    end subroutine lu_r_compensated_gather

    pure subroutine lu_z_compensated_scatter(x,correction,indices,values,n,v)
        integer, intent(in) :: n,indices(*)
        complex(dp), intent(inout) :: x(*),correction(*)
        complex(dp), intent(in) :: values(*),v
        real(dp) :: x1,x2,x3,x4,x5,x6,vr,vi,sr,si,cr,ci
        integer :: i,row
        x5=real(v,dp);x6=aimag(v)
        do i=1,n
            row=indices(i)
            x1=0.0_dp;x2=0.0_dp
            x3=real(values(i),dp);x4=aimag(values(i))
@COMPLEX_INCREMENT@
            sr=real(x(row),dp);si=aimag(x(row))
            cr=real(correction(row),dp);ci=aimag(correction(row))
            call lu_r_compensated_add(sr,cr,vr)
            call lu_r_compensated_add(si,ci,vi)
            x(row)=cmplx(sr,si,dp);correction(row)=cmplx(cr,ci,dp)
        end do
    end subroutine lu_z_compensated_scatter

    pure subroutine lu_z_compensated_gather(value,indices,values,x,n,adjoint)
        integer, intent(in) :: n,indices(*)
        complex(dp), intent(inout) :: value
        complex(dp), intent(in) :: values(*),x(*)
        logical, intent(in) :: adjoint
        complex(dp) :: coefficient
        real(dp) :: x1,x2,x3,x4,x5,x6,vr,vi,sr,si,cr,ci
        integer :: i
        sr=real(value,dp);si=aimag(value);cr=0.0_dp;ci=0.0_dp
        do i=1,n
            coefficient=values(i)
            if(adjoint) coefficient=conjg(coefficient)
            x1=0.0_dp;x2=0.0_dp
            x3=real(coefficient,dp);x4=aimag(coefficient)
            x5=real(x(indices(i)),dp);x6=aimag(x(indices(i)))
@COMPLEX_INCREMENT@
            call lu_r_compensated_add(sr,cr,vr)
            call lu_r_compensated_add(si,ci,vi)
        end do
        value=cmplx(lu_r_corrected_value(sr,cr),lu_r_corrected_value(si,ci),dp)
    end subroutine lu_z_compensated_gather
