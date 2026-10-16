program test_native_lu_compensated
    use fortsparse_kinds, only: dp,i8
    use fortsparse_native_lu, only: lu_factor_t,LU_NORMAL,LU_TRANSPOSE,LU_ADJOINT
    use fortsparse_native_lu_z, only: lu_factor_z_t
    use fortsparse_status
#ifdef FS_HAVE_ALLOC_PROBE
    use, intrinsic :: iso_c_binding, only: c_int64_t
#endif
    implicit none
    integer,parameter :: n=66,nz=2*n-1
    integer :: ptr(n+1),rows(nz),j,k,op,reuse
    real(dp) :: values(nz),truth(n,2),rhs(n,2)
    complex(dp) :: zvalues(nz),ztruth(n,2),zrhs(n,2),coefficient
    type(lu_factor_t) :: factor
    type(lu_factor_z_t) :: zfactor
    type(fortsparse_status_t) :: status
    integer(i8) :: bytes,zbytes,fixed,zfixed
#ifdef FS_HAVE_ALLOC_PROBE
    integer(c_int64_t) :: calls,requested
    interface
        subroutine begin() bind(C,name='fsparse_refinement_alloc_begin')
        end subroutine
        subroutine finish(calls,requested) bind(C,name='fsparse_refinement_alloc_end')
            import c_int64_t
            integer(c_int64_t),intent(out) :: calls,requested
        end subroutine
    end interface
#endif
    ! Independent exact dyadic arrow system. Sixty-four unit contributions
    ! are lost by ordinary sequential subtraction from 1e16. The stored RHS
    ! is exact, and the small first solution component is the observable.
    k=0
    do j=1,n
        ptr(j)=k+1
        if(j>1)then
            k=k+1;rows(k)=1;values(k)=1.0_dp;zvalues(k)=(1.0_dp,0.5_dp)
        end if
        k=k+1;rows(k)=j;values(k)=1.0_dp;zvalues(k)=(1.0_dp,0.0_dp)
    end do
    ptr(n+1)=k+1
    truth(:,1)=1.0_dp;truth(1,1)=-real(n-2,dp);truth(n,1)=1e16_dp
    truth(:,2)=-truth(:,1)
    ztruth=cmplx(truth,0.0_dp,dp)
    ztruth(1,:)=cmplx(truth(1,:),0.5_dp*truth(1,:),dp)
    call factor%factor_raw(n,ptr,rows,values,status);call checked()
    call zfactor%factor_raw(n,ptr,rows,zvalues,status);call checked()
    bytes=factor%memory_bytes();zbytes=zfactor%memory_bytes()
    do reuse=1,2
        do op=LU_NORMAL,LU_ADJOINT
            rhs=truth;zrhs=ztruth
            if(op==LU_NORMAL)then
                rhs(1,1)=1e16_dp;rhs(1,2)=-1e16_dp
                zrhs(1,1)=(1e16_dp,5e15_dp);zrhs(1,2)=-zrhs(1,1)
            else
                coefficient=(1.0_dp,0.5_dp)
                if(op==LU_ADJOINT)coefficient=conjg(coefficient)
                do j=2,n
                    rhs(j,:)=truth(j,:)+truth(1,:)
                    zrhs(j,:)=ztruth(j,:)+coefficient*ztruth(1,:)
                end do
            end if
#ifdef FS_HAVE_ALLOC_PROBE
            call begin()
#endif
            call factor%solve(rhs,status,op);call checked()
            call zfactor%solve(zrhs,status,op);call checked()
#ifdef FS_HAVE_ALLOC_PROBE
            call finish(calls,requested)
            if(calls/=0.or.requested/=0)error stop 'LU solve allocated memory'
#endif
            if(any(abs(rhs-truth)>8*epsilon(1.0_dp)*max(1.0_dp,abs(truth)))) &
                error stop 'real compensated solve lost known small component'
            if(any(abs(zrhs-ztruth)>8*epsilon(1.0_dp)*max(1.0_dp,abs(ztruth)))) &
                error stop 'complex compensated N/T/H lost known component'
            if(factor%memory_bytes()/=bytes.or.zfactor%memory_bytes()/=zbytes) &
                error stop 'retained compensation storage changed during solve'
        end do
        call factor%factor_raw(n,ptr,rows,values,status,memory_limit_bytes=bytes);call checked()
        call zfactor%factor_raw(n,ptr,rows,zvalues,status,memory_limit_bytes=zbytes);call checked()
    end do
    fixed=int(storage_size(0)/8,i8)*(10_i8*n+3_i8+nz)+24_i8*n
    zfixed=int(storage_size(0)/8,i8)*(10_i8*n+3_i8+nz)+40_i8*n
    call factor%factor_raw(n,ptr,rows,values,status,memory_limit_bytes=fixed-1)
    if(status%code/=FORTSPARSE_RESOURCE_LIMIT.or.factor%failed_required_bytes/=fixed) &
        error stop 'real compensation omitted from fixed cap'
    call zfactor%factor_raw(n,ptr,rows,zvalues,status,memory_limit_bytes=zfixed-1)
    if(status%code/=FORTSPARSE_RESOURCE_LIMIT.or.zfactor%failed_required_bytes/=zfixed) &
        error stop 'complex compensation omitted from fixed cap'
    call factor%free();call zfactor%free()
    if(factor%memory_bytes()/=0.or.zfactor%memory_bytes()/=0)error stop 'compensation free'
    call factor%factor_raw(n,ptr,rows,values,status);call checked()
    call zfactor%factor_raw(n,ptr,rows,zvalues,status);call checked()
    call factor%free();call zfactor%free()
    print '(a)','Compensated LU exact dyadic N/T/H blocks, reuse, cap and allocation-free solves pass'
contains
    subroutine checked()
        if(status%code/=FORTSPARSE_OK)then
            print *,status%code,status%msg
            error stop 'unexpected LU status'
        end if
    end subroutine
end program
