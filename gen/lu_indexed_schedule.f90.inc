    ! Scatter rows must be valid and unique; no operand aliases x.
    pure subroutine lu_z_scatter_update(x,indices,values,n,v)
        integer, intent(in) :: n,indices(*)
        complex(dp), intent(inout) :: x(*)
        complex(dp), intent(in) :: values(*),v
        real(dp) :: x1,x2,x3,x4,x5,x6,vr,vi
        integer :: i
        x5=real(v,dp);x6=aimag(v)
        !$omp simd private(x1,x2,x3,x4,vr,vi)
        do i=1,n
            x1=real(x(indices(i)),dp);x2=aimag(x(indices(i)))
            x3=real(values(i),dp);x4=aimag(values(i))
@UPDATE_BODY@
            x(indices(i))=cmplx(vr,vi,dp)
        end do
    end subroutine lu_z_scatter_update
    pure subroutine lu_z_gather_update(value,indices,values,x,n,adjoint)
        integer, intent(in) :: n,indices(*)
        complex(dp), intent(inout) :: value
        complex(dp), intent(in) :: values(*),x(*)
        logical, intent(in) :: adjoint
        complex(dp) :: coefficient
        real(dp) :: x1,x2,x3,x4,x5,x6,vr,vi
        integer :: i
        do i=1,n
            coefficient=values(i)
            if(adjoint) coefficient=conjg(coefficient)
            x1=real(value,dp);x2=aimag(value)
            x3=real(coefficient,dp);x4=aimag(coefficient)
            x5=real(x(indices(i)),dp);x6=aimag(x(indices(i)))
@UPDATE_BODY@
            value=cmplx(vr,vi,dp)
        end do
    end subroutine lu_z_gather_update
