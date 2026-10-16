program test_fortsparse_native_general_api
    use fortsparse
    implicit none
    integer, parameter :: n=5
    type(sparse_solver_t) :: solver
    type(fortsparse_status_t) :: status
    type(sparse_factor_stats_t) :: stats
    type(csc_t) :: a, adot
    type(csc_z_t) :: z, zdot
    real(dp) :: dense(n,n), values(n*n), exact(n), rhs(n), x(n), dx(n)
    real(dp) :: seed(n), bbar(n), abar(n*n), zeros(n), block(n,9), result(n,9)
    complex(dp) :: cdense(n,n), cvalues(n*n), cexact(n), crhs(n), cx(n), cdx(n)
    complex(dp) :: cseed(n), cbbar(n), cabar(n*n), czero(n), cblock(n,9), cresult(n,9)
    integer :: rows(n*n),cols(n*n),k,i,j,op,unit,nr,nc,nz,ios
    integer, allocatable :: fr(:),fc(:)
    real(dp), allocatable :: fv(:)
    integer(i8) :: peak
    logical :: reused
    character(256) :: line

    k=0
    dense=0.0_dp
    do j=1,n
        do i=1,n
            k=k+1
            rows(k)=i
            cols(k)=j
            dense(i,j)=0.03_dp*real(i-2*j,dp)
        end do
        dense(j,j)=4.0_dp+real(j,dp)
    end do
    dense(1,1)=0.0_dp
    dense(2,1)=3.0_dp
    values=reshape(dense,[n*n])
    call csc_from_triplet(n,n,rows,cols,values,a,status)
    call check(status_ok(status),'real construction')
    solver%backend_id=FORTSPARSE_BACKEND_NATIVE
    exact=[(sin(real(i,dp)),i=1,n)]
    seed=[(cos(real(i,dp)),i=1,n)]
    zeros=0.0_dp
    do op=0,2
        call sparse_factor(solver,a,status)
        call check(status_ok(status),'public real native LU factor')
        if(op==0) then
            rhs=matmul(dense,exact)
        else
            rhs=matmul(transpose(dense),exact)
        end if
        call sparse_solve_op(solver,rhs,x,op,status)
        call check(status_ok(status),'public real operation')
        call check(maxval(abs(x-exact))<5e-13_dp,'public real dense oracle')
        do j=1,9
            block(:,j)=real(j,dp)*rhs
        end do
        call sparse_solve_block_op(solver,block,result,op,status)
        call check(status_ok(status),'public real block operation')
        call check(maxval(abs(result(:,9)-9.0_dp*exact))<5e-12_dp,'real block oracle')
    end do
    call sparse_factor_stats(solver,stats,status)
    peak=stats%host_peak_bytes
    call check(stats%strict_memory_supported,'LU budget capability')
    solver%options%host_memory_limit_bytes=peak
    call sparse_factor(solver,n,n,a%nnz,a%col_ptr,a%row_idx,a%val,status)
    call check(status_ok(status),'raw capped LU borrows CSC')
    call sparse_refactor(solver,a,status,reused)
    call check(status_ok(status),'real refactor')
    call check(reused,'real pattern retained')
    rhs=matmul(dense,exact)
    call sparse_solve(solver,rhs,x,status)
    adot=a
    adot%val=0.125_dp*a%val
    call sparse_solve_jvp(solver,adot,x,zeros,dx,status)
    call check(status_ok(status),'real JVP under cap')
    call check(maxval(abs(dx+0.125_dp*exact))<1e-12_dp,'real JVP oracle')
    call sparse_solve_vjp(solver,a,x,seed,bbar,abar,status,use_primal_factors=.true.)
    call check(status_ok(status),'real VJP primal LU factors')
    call check(abs(dot_product(seed,dx)-dot_product(abar,adot%val))<1e-12_dp,'real adjoint pairing')
    call sparse_destroy(solver)
    solver%options%host_memory_limit_bytes=0_i8

    cdense=cmplx(dense,0.0_dp,dp)
    do j=1,n
        do i=1,n
            cdense(i,j)=cdense(i,j)+cmplx(0.0_dp,0.07_dp*real(i+j,dp),dp)
        end do
    end do
    cvalues=reshape(cdense,[n*n])
    call csc_from_triplet(n,n,rows,cols,cvalues,z,status)
    cexact=cmplx(exact,seed,dp)
    cseed=cmplx(seed,-exact,dp)
    czero=cmplx(0.0_dp,0.0_dp,dp)
    solver%backend_id=FORTSPARSE_BACKEND_NATIVE
    call sparse_factor(solver,z,status)
    call check(status_ok(status),'public complex native LU')
    do op=0,2
        select case(op)
        case(0)
            crhs=matmul(cdense,cexact)
        case(1)
            crhs=matmul(transpose(cdense),cexact)
        case(2)
            crhs=matmul(conjg(transpose(cdense)),cexact)
        end select
        call sparse_solve_op(solver,crhs,cx,op,status)
        call check(status_ok(status),'public complex operation')
        call check(maxval(abs(cx-cexact))<1e-12_dp,'complex dense oracle')
        do j=1,9
            cblock(:,j)=real(j,dp)*crhs
        end do
        call sparse_solve_block_op(solver,cblock,cresult,op,status)
        call check(status_ok(status),'public complex block operation')
        call check(maxval(abs(cresult(:,9)-9.0_dp*cexact))<1e-11_dp,'complex block oracle')
    end do
    call sparse_factor_stats(solver,stats,status)
    solver%options%host_memory_limit_bytes=stats%host_peak_bytes
    call sparse_factor(solver,n,n,z%nnz,z%col_ptr,z%row_idx,z%val,status)
    call check(status_ok(status),'raw capped complex LU')
    crhs=matmul(cdense,cexact)
    call sparse_solve(solver,crhs,cx,status)
    zdot=z
    zdot%val=0.125_dp*z%val
    call sparse_solve_jvp(solver,zdot,cx,czero,cdx,status)
    call check(status_ok(status),'complex JVP at cap')
    call check(maxval(abs(cdx+0.125_dp*cexact))<1e-12_dp,'complex JVP oracle')
    call sparse_solve_vjp(solver,z,cx,cseed,cbbar,cabar,status,use_primal_factors=.true.)
    call check(status_ok(status),'complex retained adjoint')
    call check(abs(real(dot_product(cseed,cdx)-dot_product(cabar,zdot%val),dp)) &
        <1e-12_dp,'complex real-pairing oracle')
    call sparse_destroy(solver)
    solver%options%host_memory_limit_bytes=0_i8

    ! Native complex HPD dispatch has different transpose mathematics from LU.
    cdense=cmplx(0.0_dp,0.0_dp,dp)
    do j=1,n
        cdense(j,j)=cmplx(3.0_dp+real(j,dp),0.0_dp,dp)
        if(j<n) then
            cdense(j,j+1)=cmplx(-0.2_dp,0.1_dp,dp)
            cdense(j+1,j)=conjg(cdense(j,j+1))
        end if
    end do
    cvalues=reshape(cdense,[n*n])
    call csc_from_triplet(n,n,rows,cols,cvalues,z,status)
    solver%backend_id=FORTSPARSE_BACKEND_NATIVE
    solver%options%symmetry=FORTSPARSE_SYM_SPD
    call sparse_factor(solver,z,status)
    call check(status_ok(status),'public HPD factor')
    do op=0,2
        if(op==1) then
            crhs=matmul(transpose(cdense),cexact)
        else
            crhs=matmul(cdense,cexact)
        end if
        call sparse_solve_op(solver,crhs,cx,op,status)
        call check(status_ok(status),'public HPD operation')
        call check(maxval(abs(cx-cexact))<1e-13_dp,'HPD independent oracle')
    end do
    crhs=matmul(cdense,cexact)
    call sparse_solve(solver,crhs,cx,status)
    zdot=z
    zdot%val=0.125_dp*z%val
    call sparse_solve_jvp(solver,zdot,cx,czero,cdx,status)
    call check(status_ok(status),'HPD allocation-free JVP')
    call check(maxval(abs(cdx+0.125_dp*cexact))<1e-13_dp,'HPD JVP oracle')
    call sparse_refactor(solver,z,status,reused)
    call check(status_ok(status),'HPD refactor')
    call check(reused,'HPD analysis retained')
    call sparse_destroy(solver)
    solver%options%host_memory_limit_bytes=0_i8

    dense=0.0_dp
    dense(1,2)=1.0_dp
    dense(2,1)=1.0_dp
    dense(3,3)=-2.0_dp
    dense(4,4)=3.0_dp
    dense(5,5)=-4.0_dp
    values=reshape(dense,[n*n])
    call csc_from_triplet(n,n,rows,cols,values,a,status)
    solver%options%symmetry=FORTSPARSE_SYM_INDEFINITE
    call sparse_factor(solver,a,status)
    call check(status_ok(status),'public pivoted indefinite factor')
    rhs=matmul(dense,exact)
    call sparse_solve_op(solver,rhs,x,FORTSPARSE_OP_TRANSPOSE,status)
    call check(status_ok(status),'public indefinite transpose')
    call check(maxval(abs(x-exact))<1e-13_dp,'indefinite known solution')
    adot=a
    adot%val=0.125_dp*a%val
    call sparse_solve_jvp(solver,adot,x,zeros,dx,status)
    call check(status_ok(status),'LDLT allocation-free JVP')
    call check(maxval(abs(dx+0.125_dp*exact))<1e-13_dp,'LDLT JVP oracle')
    call sparse_refactor(solver,a,status,reused)
    call check(status_ok(status),'public LDLT refactor')
    call check(reused,'public LDLT retains pattern')
    call sparse_destroy(solver)
    solver%options%symmetry=FORTSPARSE_SYM_GENERAL

    ! Exact rational solution of an actual retained NEO-2 matrix, independent
    ! of the LU algorithm and its own generated arithmetic.
    open(newunit=unit,file='bench/fixtures/neo2-mini5/matrix.mtx',status='old',iostat=ios)
    call check(ios==0,'NEO-2 fixture available')
    do
        read(unit,'(a)') line
        if(line(1:1)/='%') exit
    end do
    read(line,*) nr,nc,nz
    allocate(fr(nz),fc(nz),fv(nz))
    do i=1,nz
        read(unit,*) fr(i),fc(i),fv(i)
    end do
    close(unit)
    call csc_from_triplet(nr,nc,fr,fc,fv,a,status)
    solver%backend_id=FORTSPARSE_BACKEND_NATIVE
    call sparse_factor(solver,a,status)
    call check(status_ok(status),'NEO-2 general native factor')
    rhs=1.0_dp
    call sparse_solve(solver,rhs,x,status)
    exact=[-45.0_dp/22.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp/308.0_dp]
    call check(status_ok(status),'NEO-2 solve')
    call check(maxval(abs(x-exact))<2e-12_dp,'NEO-2 exact rational oracle')
    solver%options%host_memory_limit_bytes=1_i8
    call sparse_factor(solver,a,status)
    call check(status%code==FORTSPARSE_RESOURCE_LIMIT,'LU public budget failure')
    call sparse_factor_stats(solver,stats,status)
    call check(stats%budget_required_bytes>stats%budget_available_bytes,'LU budget diagnostics')
    call sparse_destroy(solver)
    solver%options%host_memory_limit_bytes=0_i8
    print '(a)', 'native general public LU/raw/transpose/adjoint/derivative/NEO-2 tests pass'
contains
    subroutine check(ok,message)
        logical,intent(in)::ok
        character(*),intent(in)::message
        if(.not.ok) then
            print '(a)',message
            error stop 1
        end if
    end subroutine check
end program test_fortsparse_native_general_api
