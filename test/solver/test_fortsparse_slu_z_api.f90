program test_fortsparse_slu_z_api
    ! Public complex static LU: known solutions of A, A^T and A^H, raw and
    ! object inputs, retained block RHS, refactor and context-copy lifetimes.
    use fortsparse
    implicit none
    integer, parameter :: n=24
    type(csc_z_t) :: a, da
    type(sparse_solver_t) :: solver, copied
    type(fortsparse_status_t) :: status
    type(sparse_factor_stats_t) :: stats
    complex(dp) :: dense(n,n), known(n,2), rhs(n,2), x(n,2), op(n,n)
    complex(dp) :: dvalues(3*n-2), db(n), dx(n), bbar(n), abar(3*n-2)
    real(dp) :: lhs, adjoint
    integer :: i,j,k,operation,raw
    logical :: reused
    dense=(0.0_dp,0.0_dp)
    do j=1,n
        dense(j,j)=cmplx(5.0_dp,0.125_dp*j,dp)
        if (j<n) then
            dense(j+1,j)=cmplx(-1.0_dp,0.25_dp,dp)
            dense(j,j+1)=cmplx(0.5_dp,-0.75_dp,dp)
        end if
        known(j,1)=cmplx(real(j,dp)/n,-0.5_dp,dp)
        known(j,2)=cmplx(0.25_dp,real(n-j,dp)/n,dp)
    end do
    a%nrow=n; a%ncol=n; a%nnz=count(dense/=(0.0_dp,0.0_dp))
    allocate(a%col_ptr(n+1),a%row_idx(a%nnz),a%val(a%nnz))
    k=1; a%col_ptr(1)=1
    do j=1,n
        do i=1,n
            if (dense(i,j)==(0.0_dp,0.0_dp)) cycle
            a%row_idx(k)=i; a%val(k)=dense(i,j); k=k+1
        end do
        a%col_ptr(j+1)=k
    end do
    do raw=0,1
        solver%backend_id=FORTSPARSE_BACKEND_NATIVE
        solver%options%lu_method=FORTSPARSE_LU_SUPERNODAL
        solver%options%nthreads=2
        if (raw==0) then
            call sparse_factor(solver,a,status)
        else
            call sparse_factor(solver,n,n,a%nnz,a%col_ptr,a%row_idx,a%val,status)
        end if
        call check(status_ok(status),'factor complex public supernodal')
        call sparse_factor_stats(solver,stats,status)
        call check(stats%lu_method_used==FORTSPARSE_LU_SUPERNODAL,'uses requested method')
        call check(.not.stats%lu_fallback,'no fallback for well-conditioned matrix')
        call check(stats%active_threads>=1,'reports actual team')
        do operation=0,2
            op=dense
            if (operation==1) op=transpose(dense)
            if (operation==2) op=conjg(transpose(dense))
            rhs=matmul(op,known)
            call sparse_solve_block_op(solver,rhs,x,operation,status)
            call check(status_ok(status),'block solve status')
            call check(maxval(abs(x-known))<2.0e-13_dp,'known normal/transpose/adjoint solution')
        end do
        rhs=matmul(dense,known)
        dvalues=0.125_dp*a%val
        da=a
        da%val=dvalues
        db=rhs(:,2)+0.125_dp*rhs(:,1)
        call sparse_solve_jvp(solver,da,known(:,1),db,dx,status)
        call check(status_ok(status),'complex solve JVP status')
        call check(maxval(abs(dx-known(:,2)))<2.0e-13_dp,'manufactured solve JVP')
        call sparse_solve_vjp(solver,a,known(:,1),known(:,2),bbar,abar,status, &
            use_primal_factors=.true.)
        call check(status_ok(status),'complex solve VJP status')
        lhs=real(sum(conjg(known(:,2))*dx),dp)
        adjoint=real(sum(conjg(bbar)*db)+sum(conjg(abar)*dvalues),dp)
        call check(abs(lhs-adjoint)<2.0e-12_dp,'real-pairing sparse adjoint identity')
        if (raw==0) then
            call sparse_refactor(solver,a,status,reused)
        else
            call sparse_refactor(solver,n,n,a%nnz,a%col_ptr,a%row_idx,a%val,status,reused)
        end if
        call check(status_ok(status).and.reused,'reuse analyzed complex panels')
        copied=solver
        call sparse_destroy(solver)
        rhs=matmul(dense,known)
        call sparse_solve_block(copied,rhs,x,status)
        call check(status_ok(status),'copied factor survives source destruction')
        call check(maxval(abs(x-known))<2.0e-13_dp,'copied complex solution')
        call sparse_destroy(copied)
    end do
    solver%backend_id=FORTSPARSE_BACKEND_NATIVE
    solver%options%lu_method=FORTSPARSE_LU_SUPERNODAL
    solver%refine=.true.
    solver%options%host_memory_limit_bytes=1_i8
    call sparse_factor(solver,a,status)
    call check(status%code==FORTSPARSE_RESOURCE_LIMIT,'complex owned-payload cap')
    solver%options%host_memory_limit_bytes=0_i8
    solver%options%pivot_perturbation=0.99_dp
    call sparse_factor(solver,a,status)
    call check(status_ok(status),'complex scalar fallback factor')
    call sparse_factor_stats(solver,stats,status)
    call check(stats%lu_fallback,'perturbed complex attempt reports fallback')
    rhs=matmul(dense,known)
    call sparse_solve_block(solver,rhs,x,status)
    call check(status_ok(status),'complex fallback solve')
    call check(maxval(abs(x-known))<2.0e-13_dp,'complex fallback known solution')
    call sparse_destroy(solver)
    print *, 'complex supernodal public API passed'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not.ok) then
            print *, label,': ',status%code,trim(status%msg)
            error stop 1
        end if
    end subroutine
end program
