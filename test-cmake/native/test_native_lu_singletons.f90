program test_native_lu_singletons
    use fortsparse
    use fortsparse_native_lu, only: lu_factor_t
    use fortsparse_native_lu_z, only: lu_factor_z_t
    implicit none
    integer, parameter :: n=6
    integer, parameter :: pr(n)=[3,5,1,6,2,4], pc(n)=[2,4,6,1,5,3]
    integer, parameter :: requested(n)=[3,4,1,2,5,6]
    integer, parameter :: expected_order(n)=[2,4,6,3,1,5]
    type(lu_factor_t) :: f
    type(lu_factor_z_t) :: zf
    type(csc_t) :: a
    type(csc_z_t) :: za
    type(fortsparse_status_t) :: status
    integer :: i,j,k,operation,prune,change,actual(n),rows(n*n+1),cols(n*n+1)
    real(dp) :: base(n,n),dense(n,n),values(n*n+1),known(n,2),rhs(n,2)
    complex(dp) :: zdense(n,n),zvalues(n*n+1),zknown(n,2),zrhs(n,2)
    integer(i8) :: bytes,zbytes
    logical :: matched

    base=0
    base(1,1)=2
    base(2,:2)=[100.0_dp,-3.0_dp]
    base(3,2:3)=[50.0_dp,5.0_dp]
    base(4,:)=[0.0_dp,2.0_dp,1.0_dp,7.0_dp,1.0_dp,2.0_dp]
    base(5,:)=[3.0_dp,0.0_dp,-1.0_dp,2.0_dp,8.0_dp,1.0_dp]
    base(6,:)=[2.0_dp,4.0_dp,0.0_dp,1.0_dp,3.0_dp,9.0_dp]
    do j=1,2
        do i=1,n
            known(i,j)=i+2*j
            zknown(i,j)=cmplx(i+2*j,2*i-j,dp)
        end do
    end do
    do prune=0,1
        do change=0,1
            ! Numerical activation of an explicitly stored zero changes the
            ! singleton structure while retaining exactly the same CSC graph.
            if(change==1) base(1,4)=1
            do j=1,n
                do i=1,n
                    dense(pr(i),pc(j))=base(i,j)
                    zdense(pr(i),pc(j))=base(i,j)*cmplx(1.0_dp,0.1_dp*i,dp)* &
                        cmplx(1.0_dp,-0.2_dp*j,dp)
                end do
            end do
            k=0
            do j=1,n
                do i=1,n
                    k=k+1
                    rows(k)=i;cols(k)=j
                    values(k)=dense(i,j);zvalues(k)=zdense(i,j)
                end do
            end do
            ! A duplicate COO zero coalesces; explicit zero CSC entries remain.
            rows(k+1)=1;cols(k+1)=1;values(k+1)=0;zvalues(k+1)=0
            call csc_from_triplet(n,n,rows,cols,values,a,status)
            call checked(status)
            call csc_from_triplet(n,n,rows,cols,zvalues,za,status)
            call checked(status)
            call f%factor(a,status,column_permutation=requested,pruning=prune==1)
            call checked(status)
            call zf%factor(za,status,column_permutation=requested,pruning=prune==1)
            call checked(status)
            if(change==0) then
                call f%copy_column_order(n,a%col_ptr,a%row_idx,actual,matched)
                if(.not.matched.or.any(actual/=expected_order)) &
                    error stop 'singleton prefix did not retain requested core order'
            else
                if(.not.f%last_reused_pattern.or..not.zf%last_reused_pattern) &
                    error stop 'same-pattern numerical singleton change not reused'
            end if
            bytes=f%memory_bytes();zbytes=zf%memory_bytes()
            do operation=0,2
                if(operation==0) then
                    rhs=matmul(dense,known)
                    zrhs=matmul(zdense,zknown)
                else
                    rhs=matmul(transpose(dense),known)
                    if(operation==1) then
                        zrhs=matmul(transpose(zdense),zknown)
                    else
                        zrhs=matmul(conjg(transpose(zdense)),zknown)
                    end if
                end if
                call f%solve(rhs,status,operation)
                call checked(status)
                call zf%solve(zrhs,status,operation)
                call checked(status)
                if(maxval(abs(rhs-known))>1e-10_dp) error stop 'real singleton/core NTH block'
                if(maxval(abs(zrhs-zknown))>1e-10_dp) error stop 'complex singleton/core NTH block'
                if(f%memory_bytes()/=bytes.or.zf%memory_bytes()/=zbytes) &
                    error stop 'singleton solve allocated new storage'
            end do
        end do
        ! Invalid raw CSC never reaches singleton preprocessing; failure is
        ! recoverable and cannot expose the previous numerical factors.
        a%row_idx(2)=a%row_idx(1)
        call f%factor(a,status)
        if(status%code/=FORTSPARSE_INVALID_MATRIX) error stop 'duplicate raw CSC accepted'
        a%row_idx(2)=2
        a%val=0
        call f%factor(a,status)
        if(status%code/=FORTSPARSE_SINGULAR) error stop 'zero rows accepted'
        rhs=0
        call f%solve(rhs,status)
        if(status%code/=FORTSPARSE_NOT_FACTORED) error stop 'failed factors remained usable'
        call f%free();call zf%free()
        base(1,4)=0
    end do
contains
    subroutine checked(status)
        type(fortsparse_status_t), intent(in) :: status
        if(.not.status_ok(status)) then
            print *,status%code,status%msg
            error stop 'singleton factor/solve failure'
        end if
    end subroutine
end program
