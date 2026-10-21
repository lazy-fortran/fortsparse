program test_fortsparse_mixed_lu
    use fortsparse
    use fortsparse_mixed_lu
    use fortsparse_native_lu, only: lu_factor_t
    use fortnum_precision, only: FP32,FP64
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
    implicit none
    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda,*), b(ldb,*)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface
    integer, parameter :: n=19
    type(csc_t) :: a,ill,matrices(2)
    type(csc_z_t) :: az
    type(mixed_lu_t) :: mixed,wide,systems(2)
    type(mixed_lu_z_t) :: zmixed,zsystems(2)
    type(lu_factor_t) :: raw
    type(fortsparse_status_t) :: status,statuses(2)
    type(mixed_solve_report_t) :: report,reports(2)
    real(dp) :: dense(n,n),known(n),b(n),x(n),bb(n,2),xx(n,2)
    real(dp) :: batch_b(2*n),batch_x(2*n),saved(n),inverse(n,n)
    complex(dp) :: zdense(n,n),zknown(n),zb(n),zx(n),zbb(n,2),zxx(n,2)
    complex(dp) :: zbatch_b(2*n),zbatch_x(2*n)
    integer :: rows(n*n),cols(n*n),q(n),i,j,k,op
    real(dp) :: expected,initial_error
    integer(i8) :: full_bytes
    k=0
    dense=0
    do j=1,n
        q(j)=n+1-j
        do i=1,n
            k=k+1
            rows(k)=i
            cols(k)=j
            ! Strict diagonal dominance, nonsymmetry, non-binary coefficients.
            dense(i,j)=real(mod(7*i+3*j,11)-5,dp)/71
            if (i==j) dense(i,j)=4+real(i,dp)/31
            zdense(i,j)=cmplx(dense(i,j),real(mod(i+2*j,7)-3,dp)/53,dp)
        end do
    end do
    known=[(sin(real(i,dp)),i=1,n)]
    zknown=cmplx(known,cos(2*known),dp)
    call csc_from_triplet(n,n,rows,cols,reshape(dense,[n*n]),a,status)
    call require(status%code==FORTSPARSE_OK,'real CSC')
    call csc_from_triplet(n,n,rows,cols,reshape(zdense,[n*n]),az,status)
    call require(status%code==FORTSPARSE_OK,'complex CSC')
    call wide%factor(a,status)
    call require(status%code==FORTSPARSE_OK,'default FP64 factor')
    full_bytes=wide%bytes()
    call mixed%factor(a,status,FP32,column_order=q)
    call require(status%code==FORTSPARSE_OK,'fp32 storage factor')
    call require(mixed%bytes()<full_bytes,'factor payload reduction')
    call zmixed%factor(az,status,FP32,column_order=q)
    call require(status%code==FORTSPARSE_OK,'complex fp32 storage')
    do op=0,2
        if (op==0) then
            b=matmul(dense,known)
            zb=matmul(zdense,zknown)
        else if (op==1) then
            b=matmul(transpose(dense),known)
            zb=matmul(transpose(zdense),zknown)
        else
            b=matmul(transpose(dense),known)
            zb=matmul(conjg(transpose(zdense)),zknown)
        end if
        saved=b
        call mixed%apply(saved,status,op)
        call require(status%code==FORTSPARSE_OK,'preconditioner apply')
        initial_error=maxval(abs(saved-known))
        call mixed%solve(b,x,report,status,op)
        call require(status%code==FORTSPARSE_OK.and.report%converged,'refined solve')
        call require(report%storage==FP32.and.report%steps>0,'actual mixed correction')
        call require(maxval(abs(x-known))<3e-13_dp,'real manufactured solution')
        call require(maxval(abs(x-known))<initial_error*1e-3_dp,'defect reduction')
        call require(report%backward_error<32*epsilon(1.0_dp),'backward error')
        call require(report%componentwise_error>=0,'componentwise error available')
        call require(report%rcond_estimate>1e-3_dp,'condition estimate')
        call zmixed%solve(zb,zx,report,status,op)
        call require(status%code==FORTSPARSE_OK.and.report%converged,'complex refined solve')
        call require(report%storage==FP32.and.report%steps>0,'actual complex mixed correction')
        call require(maxval(abs(zx-zknown))<5e-13_dp,'complex N/T/H manufactured solution')
        bb(:,1)=b
        bb(:,2)=2*b
        call mixed%solve_block(bb,xx,reports,statuses,op)
        call require(all(statuses%code==FORTSPARSE_OK),'real block solve')
        call require(maxval(abs(xx(:,2)-2*known))<1e-12_dp,'real block solution')
        zbb(:,1)=zb
        zbb(:,2)=2*zb
        call zmixed%solve_block(zbb,zxx,reports,statuses,op)
        call require(all(statuses%code==FORTSPARSE_OK),'complex block solve')
        call require(maxval(abs(zxx(:,2)-2*zknown))<1e-12_dp,'complex block solution')
    end do
    ! Exact inverse through independent dense LAPACK as a secondary oracle.
    inverse=0
    do i=1,n
        inverse(i,i)=1
    end do
    call dense_inverse(dense,inverse)
    expected=1/(maxval(sum(abs(dense),dim=2))*maxval(sum(abs(inverse),dim=2)))
    b=matmul(dense,known)
    call mixed%solve(b,x,report,status)
    call require(report%rcond_estimate>=expected*(1-1e-10_dp),'estimate is a lower norm estimate')
    call require(report%rcond_estimate<5*expected,'estimate oracle factor-five accuracy')
    do i=1,2
        call systems(i)%factor(a,status,merge(FP32,FP64,i==1))
        call require(status%code==FORTSPARSE_OK,'batch factor')
        call zsystems(i)%factor(az,status,FP32)
        call require(status%code==FORTSPARSE_OK,'complex batch factor')
    end do
    batch_b(:n)=b
    batch_b(n+1:)=2*b
    call mixed_solve_batch(systems,[1,n+1,2*n+1],batch_b,batch_x,reports,statuses)
    call require(all(statuses%code==FORTSPARSE_OK),'batched per-n solve')
    call require(maxval(abs(batch_x(n+1:)-2*known))<1e-12_dp,'batch oracle')
    zb=matmul(zdense,zknown)
    zbatch_b(:n)=zb
    zbatch_b(n+1:)=2*zb
    call mixed_solve_batch(zsystems,[1,n+1,2*n+1],zbatch_b,zbatch_x,reports,statuses)
    call require(all(statuses%code==FORTSPARSE_OK),'complex batch solve')
    call require(maxval(abs(zbatch_x(n+1:)-2*zknown))<1e-12_dp,'complex batch oracle')
    call block_jacobi_apply(systems,[1,n+1,2*n+1],batch_b,batch_x,statuses)
    call require(all(statuses%code==FORTSPARSE_OK),'block Jacobi apply')
    call require(maxval(abs(batch_x(:n)-known))<1e-6_dp,'fp32 block Jacobi oracle')
    call require(maxval(abs(batch_x(n+1:)-2*known))<1e-12_dp,'FP64 block Jacobi oracle')
    call block_jacobi_apply(zsystems,[1,n+1,2*n+1],zbatch_b,zbatch_x,statuses)
    call require(all(statuses%code==FORTSPARSE_OK),'complex block Jacobi')
    call require(maxval(abs(zbatch_x(:n)-zknown))<1e-6_dp,'complex Jacobi oracle')
    call mixed_solve_batch(systems,[1,0,2*n+1],batch_b,batch_x,reports,statuses)
    call require(all(statuses%code==FORTSPARSE_INVALID_MATRIX),'invalid offsets')
    call mixed%factor(a,status,FP32,memory_limit_bytes=1_i8)
    call require(status%code==FORTSPARSE_RESOURCE_LIMIT,'checker cap')
    call mixed%solve(b,x,report,status)
    call require(status%code==FORTSPARSE_NOT_FACTORED,'failed factor not usable')
    call mixed%factor(a,status,FP32,memory_limit_bytes=100000_i8)
    call require(status%code==FORTSPARSE_OK,'cap recovery')
    call mixed%factor(a,status,FP32,device=0)
    call require(status%code==FORTSPARSE_BACKEND_UNAVAILABLE,'GPU seam explicit')
    call mixed%factor(a,status,FP32)
    saved=b
    saved(1)=ieee_value(1.0_dp,ieee_quiet_nan)
    call mixed%solve(saved,x,report,status)
    call require(status%code==FORTSPARSE_INVALID_MATRIX,'nonfinite RHS')
    ! A known ill-conditioned diagonal must stay FP64 under adaptive storage.
    call csc_from_triplet(2,2,[1,2],[1,2],[1.0_dp,1e-10_dp],ill,status)
    call mixed%factor(ill,status,FP32)
    call require(status%code==FORTSPARSE_OK,'ill-conditioned factor')
    call mixed%solve([1.0_dp,1e-10_dp],x(:2),report,status)
    call require(status%code==FORTSPARSE_OK.and.report%storage==FP64,'adaptive wide storage')
    call require(abs(report%rcond_estimate-1e-10_dp)<1e-20_dp,'diagonal condition oracle')
    matrices(1)=a
    matrices(2)=ill
    call mixed_factor_batch(systems,matrices,statuses,FP32)
    call require(all(statuses%code==FORTSPARSE_OK),'different-size factor batch')
    batch_b(n+1:n+2)=[1.0_dp,1e-10_dp]
    call mixed_solve_batch(systems,[1,n+1,n+3],batch_b(:n+2),batch_x(:n+2),reports,statuses)
    call require(all(statuses%code==FORTSPARSE_OK),'different-size batch')
    call require(all(abs(batch_x(n+1:n+2)-1)<1e-13_dp),'different-size oracle')
    ! Low-level compression is atomic on cap failure and refactor-safe.
    call raw%factor(a,status)
    full_bytes=raw%memory_bytes()
    call raw%factor(a,status,memory_limit_bytes=full_bytes)
    call require(status%code==FORTSPARSE_OK,'bounded retained factors')
    call raw%compress(status)
    call require(status%code==FORTSPARSE_RESOURCE_LIMIT,'narrowing overlap cap')
    saved=b
    call raw%solve(saved,status)
    call require(status%code==FORTSPARSE_OK,'wide factors survive rejected narrowing')
    call raw%free()
    call raw%factor(a,status)
    call raw%compress(status)
    call require(status%code==FORTSPARSE_OK,'raw compress')
    call raw%factor(a,status)
    call require(status%code==FORTSPARSE_OK,'refactor after compress')
    call raw%free()
    call mixed%free()
    call zmixed%free()
    print *, 'mixed LU real/complex N/T/H, block/batch, caps and diagnostics PASS'
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not.ok) then
            print *, 'FAIL: ',message
            error stop 1
        end if
    end subroutine
    subroutine dense_inverse(a,inverse)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(inout) :: inverse(:,:)
        real(dp) :: copy(size(a,1),size(a,2))
        integer :: piv(size(a,1)),info
        copy=a
        call dgesv(size(a,1),size(a,1),copy,size(a,1),piv,inverse,size(a,1),info)
        call require(info==0,'LAPACK inverse oracle')
    end subroutine
end program
