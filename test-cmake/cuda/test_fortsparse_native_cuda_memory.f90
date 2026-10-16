program test_fortsparse_native_cuda_memory
    use fortsparse
    implicit none
    integer, parameter :: n=6
    type(csc_t) :: a, changed
    type(sparse_solver_t) :: solver
    type(sparse_factor_stats_t) :: stats
    type(fortsparse_status_t) :: status
    real(dp) :: dense(n,n), exact(n), b(n), x(n)
    integer :: rows(n*n), cols(n*n), i, j, p
    real(dp) :: values(n*n)
    integer(i8) :: required
    logical :: reused

    dense=0.0_dp
    do i=1,n
        dense(i,i)=4.0_dp
        if(i<n) then
            dense(i,i+1)=-1.0_dp
            dense(i+1,i)=-1.0_dp
        end if
    end do
    p=0
    do j=1,n
        do i=1,n
            if(dense(i,j)==0) cycle
            p=p+1;rows(p)=i;cols(p)=j;values(p)=dense(i,j)
        end do
    end do
    call csc_from_triplet(n,n,rows(:p),cols(:p),values(:p),a,status)
    call require(status_ok(status),'original CSC')
    exact=[(real(i,dp)/8.0_dp,i=1,n)]
    b=matmul(dense,exact)
    solver%backend_id=FORTSPARSE_BACKEND_NATIVE
    solver%options%symmetry=FORTSPARSE_SYM_SPD
    solver%options%device=0
    call sparse_factor(solver,a,status)
    call require(status_ok(status),'uncapped GPU factor')
    call sparse_factor_stats(solver,stats,status)
    call require(stats%strict_device_memory_supported,'device payload tracked')
    call require(.not.stats%strict_memory_supported,'host cap not claimed')
    call require(stats%host_peak_bytes==-1_i8,'host peak remains unknown')
    required=stats%device_peak_bytes
    call require(required>0_i8.and.stats%device_bytes<=required,'device peak')
    solver%options%device_memory_limit_bytes=required-1_i8
    call sparse_refactor(solver,a,status,reused)
    call require(status%code==FORTSPARSE_RESOURCE_LIMIT,'cap minus one refuses')
    call require(.not.reused,'changed device cap rebuilds')
    call sparse_factor_stats(solver,stats,status)
    call require(stats%device_bytes==0_i8,'failed build frees device payload')
    call require(stats%device_peak_bytes<=required-1_i8,'successful live peak within cap')
    call require(stats%device_budget_required_bytes> &
        stats%device_budget_available_bytes,'device refusal diagnostics')
    call sparse_solve(solver,b,x,status)
    call require(status%code==FORTSPARSE_NOT_FACTORED,'refused factors unusable')
    solver%options%device_memory_limit_bytes=required
    call sparse_factor(solver,a,status)
    call require(status_ok(status),'exact device cap recovers')
    call verify()
    changed=a;changed%val=2.0_dp*a%val
    call sparse_refactor(solver,changed,status,reused)
    call require(status_ok(status).and.reused,'same policy refactor reuses plan')
    call sparse_solve(solver,2.0_dp*b,x,status)
    call require(status_ok(status),'refactor solve')
    call require(maxval(abs(x-exact))<1e-13_dp,'refactor known x')
    call sparse_factor_stats(solver,stats,status)
    call require(stats%device_peak_bytes==required,'refactor no new device allocation')
    solver%options%host_memory_limit_bytes=100000_i8
    call sparse_factor(solver,a,status)
    call require(status%code==FORTSPARSE_BACKEND_UNAVAILABLE,'GPU host cap unsupported')
    solver%options%host_memory_limit_bytes=0_i8
    solver%options%memory_policy=FORTSPARSE_MEMORY_LOW
    call sparse_factor(solver,a,status)
    call require(status%code==FORTSPARSE_BACKEND_UNAVAILABLE,'GPU LOW unsupported')
    call sparse_destroy(solver)
    print '(a)', 'Public CUDA exact device payload cap, reuse and recovery pass'
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if(.not.ok) then
            print *,trim(message),status%code,trim(status%msg)
            error stop 1
        end if
    end subroutine
    subroutine verify()
        call sparse_solve(solver,b,x,status)
        call require(status_ok(status),'known x solve')
        call require(maxval(abs(x-exact))<1e-13_dp,'independent known x')
        call require(maxval(abs(matmul(dense,x)-b))<1e-12_dp,'original A residual')
    end subroutine
end program
