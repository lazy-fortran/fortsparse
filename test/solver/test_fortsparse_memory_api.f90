program test_fortsparse_memory_api
    use fortsparse
    implicit none
    integer, parameter :: n = 73
    type(csc_t) :: a, tangent
    type(sparse_solver_t) :: solver, other
    type(sparse_memory_plan_t) :: plan
    type(sparse_factor_stats_t) :: stats
    type(fortsparse_status_t) :: status
    integer :: i, j, k, rows(3*n-2), cols(3*n-2)
    real(dp) :: vals(3*n-2), exact(n), rhs(n), x(n), dx(n), zeros(n)
    real(dp) :: block(n,9), result(n,9), seed(n), bbar(n), abar(3*n-2)
    real(dp), pointer :: vector(:)
    integer(i8) :: required
    logical :: reused
    k = 0
    do j = 1, n
        do i = max(1,j-1), min(n,j+1)
            k = k+1
            rows(k) = i
            cols(k) = j
            vals(k) = -1.0_dp
            if (i == j) vals(k) = 4.0_dp
        end do
    end do
    call csc_from_triplet(n,n,rows,cols,vals,a,status)
    call check(status_ok(status), 'construct')
    solver%backend_id = FORTSPARSE_BACKEND_NATIVE
    solver%options%symmetry = FORTSPARSE_SYM_SPD
    solver%options%memory_policy = FORTSPARSE_MEMORY_LOW
    call sparse_analyze(solver,a,plan,status,nrhs=9)
    call check(status_ok(status), 'analysis')
    call check(plan%prediction_quality == MEMORY_PLAN_EXACT_LAYOUT, 'layout quality')
    call check(plan%symbolic_nnz == 2*n-1, 'independent tridiagonal fill')
    required = 12_i8*(2*n-1)+4_i8*(3*n-2)+36_i8*n+8_i8
    call check(plan%peak_host_bytes == required, 'independent payload count')
    solver%options%host_memory_limit_bytes = required-1_i8
    call sparse_factor(solver,a,status)
    call check(status%code == FORTSPARSE_RESOURCE_LIMIT, 'one byte below cap')
    call check(.not. solver%factored, 'failure invalidates')
    call sparse_factor_stats(solver,stats,status)
    call check(status_ok(status), 'failed attempt statistics')
    call check(stats%budget_required_bytes > stats%budget_available_bytes, 'cap diagnostics')
    solver%options%host_memory_limit_bytes = required
    call sparse_factor(solver,a,status)
    call check(status_ok(status), 'exact cap succeeds')
    call sparse_factor_stats(solver,stats,status)
    call check(stats%strict_memory_supported, 'strict support')
    call check(stats%host_current_bytes == required, 'current payload')
    call check(stats%host_peak_bytes == required, 'peak payload')
    call sparse_analyze(solver,a,plan,status,nrhs=9)
    call check(status%code == FORTSPARSE_RESOURCE_LIMIT, 'planning scratch shares live cap')
    call check(solver%factored, 'planning failure preserves factors')
    exact = [(sin(real(i,dp)),i=1,n)]
    rhs = 4.0_dp*exact
    rhs(2:) = rhs(2:)-exact(:n-1)
    rhs(:n-1) = rhs(:n-1)-exact(2:)
    call sparse_solve(solver,rhs,x,status)
    call check(status_ok(status), 'bounded solve')
    call check(maxval(abs(x-exact)) < 2.0e-14_dp, 'independent solution')
    do j = 1, 9
        block(:,j) = j*rhs
    end do
    call sparse_solve_block(solver,block,result,status)
    call check(status_ok(status), 'nine RHS without extra allocation')
    do j = 1, 9
        call check(maxval(abs(result(:,j)-j*exact)) < 2.0e-13_dp, 'block oracle')
    end do
    call sparse_solve_op(solver,rhs,x,FORTSPARSE_OP_TRANSPOSE,status)
    call check(status_ok(status), 'retained transpose solve')
    call check(maxval(abs(x-exact)) < 2.0e-14_dp, 'transpose solution')
    call sparse_solve_block_op(solver,block,result,FORTSPARSE_OP_ADJOINT,status)
    call check(status_ok(status), 'retained block adjoint solve')
    call sparse_solve_op(solver,rhs,x,-1,status)
    call check(status%code == FORTSPARSE_INVALID_MATRIX, 'invalid operation')
    call sparse_solve(solver,rhs,x,status)
    seed = [(cos(real(i,dp)),i=1,n)]
    call sparse_solve_vjp(solver,a,x,seed,bbar,abar,status,use_primal_factors=.true.)
    call check(status_ok(status), 'VJP reuses primal factors')
    tangent = a
    tangent%val = 0.25_dp*a%val
    zeros = 0.0_dp
    call sparse_solve_jvp(solver,tangent,x,zeros,dx,status)
    call check(status_ok(status), 'JVP under exact cap')
    call check(maxval(abs(dx+0.25_dp*exact)) < 2.0e-14_dp, 'JVP oracle')
    call check(abs(dot_product(seed,dx)-dot_product(abar,tangent%val)) < &
        3.0e-13_dp, 'retained factor adjoint pairing')
    vector => sparse_vector(solver,n)
    call check(.not. associated(vector), 'vector growth rejected at cap')
    call sparse_factor_stats(solver,stats,status)
    call check(stats%host_current_bytes == required, 'rejected growth preserves ledger')
    a%val = 2.0_dp*a%val
    do i = 1, 4
        call sparse_refactor(solver,a,status,reused)
        call check(status_ok(status), 'bounded refactor')
        call check(reused, 'analysis reused')
        call sparse_solve(solver,2.0_dp*rhs,x,status)
        call check(maxval(abs(x-exact)) < 2.0e-14_dp, 'refactor oracle')
        call sparse_factor_stats(solver,stats,status)
        call check(stats%host_peak_bytes == required, 'stable lifetime peak')
    end do
    other%backend_id = FORTSPARSE_BACKEND_NATIVE
    other%options = solver%options
    other%options%host_memory_limit_bytes = required+8_i8*n
    call sparse_factor(other,a,status)
    call check(status_ok(status), 'independent context')
    vector => sparse_vector(other,n)
    call check(associated(vector), 'reserved vector fits')
    call sparse_factor_stats(other,stats,status)
    call check(stats%host_current_bytes == required+8_i8*n, 'vector accounted')
    solver%options%host_memory_limit_bytes = required-1_i8
    call sparse_refactor(solver,a,status,reused)
    call check(status%code == FORTSPARSE_RESOURCE_LIMIT, 'lower cap rechecked')
    call check(.not. reused, 'lower cap resets plan')
    solver%options%host_memory_limit_bytes = required
    call sparse_factor(solver,a,status)
    call check(status_ok(status), 'recovery after exhausted budget')
    solver%options%ordering = FORTSPARSE_ORDER_ND
    call sparse_factor(solver,a,status)
    call check(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, 'unbounded ND fail closed')
    solver%options%ordering = FORTSPARSE_ORDER_DEFAULT
    solver%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC
    call sparse_factor(solver,a,status)
    call check(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, 'external hard cap unsupported')
    solver%backend_id = FORTSPARSE_BACKEND_NATIVE
    solver%options%host_memory_limit_bytes = -1_i8
    call sparse_factor(solver,a,status)
    call check(status%code == FORTSPARSE_INVALID_MATRIX, 'negative cap')
    call sparse_destroy(solver)
    call sparse_destroy(other)
    call blocked_budget()
    print '(a)', 'bounded native SPD analysis, solve, JVP, reuse and recovery pass'
contains
    subroutine blocked_budget()
        integer, parameter :: nb=96
        type(sparse_solver_t) :: bounded
        type(sparse_memory_plan_t) :: predicted
        type(sparse_factor_stats_t) :: measured
        type(csc_t) :: matrix
        integer :: ii,jj,kk,rr(nb*nb),cc(nb*nb)
        real(dp) :: vv(nb*nb),dd(nb,nb),answer(nb),input(nb),output(nb)
        dd=0.01_dp
        do ii=1,nb
            dd(ii,ii)=2.0_dp
        end do
        kk=0
        do jj=1,nb
            do ii=1,nb
                kk=kk+1
                rr(kk)=ii
                cc(kk)=jj
                vv(kk)=dd(ii,jj)
            end do
        end do
        call csc_from_triplet(nb,nb,rr,cc,vv,matrix,status)
        bounded%backend_id=FORTSPARSE_BACKEND_NATIVE
        bounded%options%symmetry=FORTSPARSE_SYM_SPD
        bounded%options%nthreads=2
        bounded%options%host_memory_limit_bytes=4_i8*1024_i8*1024_i8
        call sparse_factor(bounded,matrix,status)
        call check(status_ok(status),'public bounded blocked factor')
        call sparse_factor_stats(bounded,measured,status)
        call check(trim(measured%native_factor_family)=='CHOLESKY_BLOCK','blocked schedule reported')
        call check(measured%host_peak_bytes<=bounded%options%host_memory_limit_bytes,'blocked cap')
        call check(measured%host_current_bytes>measured%factor_storage_bytes,'metadata counted')
        answer=[(sin(real(ii,dp)),ii=1,nb)]
        input=matmul(dd,answer)
        call sparse_solve(bounded,input,output,status)
        call check(status_ok(status),'bounded blocked solve')
        call check(maxval(abs(output-answer))<2e-14_dp,'blocked independent oracle')
        call sparse_refactor(bounded,matrix,status,reused)
        call check(status_ok(status),'bounded blocked refactor')
        call check(reused,'blocked analysis reuse')
        call sparse_destroy(bounded)
        bounded%options%memory_policy=FORTSPARSE_MEMORY_LOW
        call sparse_analyze(bounded,matrix,predicted,status)
        call check(status_ok(status),'fallback compact analysis')
        bounded%options%memory_policy=FORTSPARSE_MEMORY_AUTO
        bounded%options%host_memory_limit_bytes=predicted%peak_host_bytes
        call sparse_factor(bounded,matrix,status)
        call check(status_ok(status),'AUTO fits compact layout when panels cannot fit')
        call sparse_factor_stats(bounded,measured,status)
        call check(measured%compact_memory_used,'compact fallback explicitly reported')
        call check(measured%host_peak_bytes<=predicted%peak_host_bytes,'fallback cap respected')
        call sparse_solve(bounded,input,output,status)
        call check(maxval(abs(output-answer))<2e-14_dp,'fallback independent oracle')
        call sparse_destroy(bounded)
    end subroutine blocked_budget

    subroutine check(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print '(a)', message
            error stop 1
        end if
    end subroutine check
end program test_fortsparse_memory_api
