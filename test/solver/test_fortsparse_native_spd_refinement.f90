program test_fortsparse_native_spd_refinement
    use fortsparse
    use fortsparse_native, only: native_backend_t
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    type(native_backend_t) :: reference,blocked
    type(csc_t) :: a,large,changed
    type(sparse_solver_t) :: solver
    type(sparse_memory_plan_t) :: plan
    type(fortsparse_status_t) :: status
    type(sparse_factor_stats_t) :: plain,refined,stats
    real(dp) :: b(2),x(2),rhs(2,2),answer(2,2),dense(2,2),delta
    real(dp), allocatable :: lb(:),lx(:)
    real(dp), pointer :: user(:)
    integer, allocatable :: rows(:),cols(:)
    real(dp), allocatable :: values(:)
    integer :: i,j,k,op
    logical :: reused
    reference%options%symmetry=FORTSPARSE_SYM_SPD
    reference%reference_only=.true.
    reference%options%memory_policy=FORTSPARSE_MEMORY_LOW
    dense=reshape([4._dp,1._dp,1._dp,3._dp],[2,2])
    b=matmul(dense,[1._dp,-2._dp])
    call csc_from_triplet(2,2,[1,2,1,2],[1,1,2,2],reshape(dense,[4]),a,status)
    call checked()
    solver%backend_id=FORTSPARSE_BACKEND_NATIVE
    solver%options%symmetry=FORTSPARSE_SYM_SPD
    solver%options%memory_policy=FORTSPARSE_MEMORY_LOW
    solver%refine=.true.
    call sparse_analyze(solver,a,plan,status)
    call checked()
    if(.not.plan%fits.or.plan%prediction_quality/=MEMORY_PLAN_EXACT_LAYOUT) error stop 'refined plan unavailable'
    call sparse_factor(solver,a,status)
    call checked()
    call sparse_analyze(solver,a,plan,status)
    call checked()
    if(.not.solver%factored) error stop 'planning changed retained factors'
    call sparse_solve(solver,b,x,status)
    call checked()
    call oracle()
    call sparse_destroy(solver)
    solver%refine=.false.
    call sparse_analyze(solver,a,plan,status)
    call checked()
    if(.not.plan%fits.or.plan%prediction_quality/=MEMORY_PLAN_EXACT_LAYOUT) error stop 'ordinary plan changed'
    call reference%factor_real(a,.false.,status)
    call checked()
    call reference%factor_stats(plain,status)
    call checked()
    call reference%factor_real(a,.true.,status)
    call checked()
    call reference%factor_stats(refined,status)
    call checked()
    if(refined%host_current_bytes-plain%host_current_bytes/=92_i8) error stop 'original CSC/two-vector payload'
    a%val=9*a%val ! Caller values must not become the later residual operator.
    do op=0,2
        call reference%solve_real_op(b,x,op,status)
        call checked()
        call oracle()
    end do
    a%val=reshape(dense,[4])
    ! Treat the fast factor as untrusted: a retained inaccurate diagonal must
    ! trigger actual correction, checked against the original dense operator.
    delta=1e-4_dp
    reference%lv(reference%lp(1))=reference%lv(reference%lp(1))*(1+delta)
    call reference%solve_real(b,x,status)
    call checked()
    call oracle()
    call reference%factor_stats(stats,status)
    call checked()
    if(stats%refinement_steps<1) error stop 'correction fixture did not correct'
    rhs(:,1)=b
    rhs(:,2)=2*b
    call reference%solve_real_block_op(rhs,answer,2,status)
    call checked()
    if(maxval(abs(matmul(dense,answer)-rhs))>1e-12_dp) error stop 'original block residual'
    ! Larger corruption cannot meet the unchanged three-correction budget.
    reference%lv(reference%lp(1))=2*reference%lv(reference%lp(1))
    call reference%solve_real(b,x,status)
    if(status%code/=FORTSPARSE_SINGULAR) error stop 'untrusted factor must fail target'
    if(reference%stats%refinement_steps/=3.or.reference%stats%backward_error<=32*epsilon(1._dp)) &
        error stop 'failed refinement diagnostics'
    call reference%refactor_real(a,.true.,status,reused)
    call checked()
    if(.not.reused) error stop 'same graph did not reuse'
    user=>reference%vector(2)
    if(.not.associated(user)) error stop 'user vector allocation'
    user=17
    changed=a
    changed%val=2*a%val
    call reference%refactor_real(changed,.true.,status,reused)
    call checked()
    if(.not.reused.or.any(user/=17)) error stop 'refactor/user vector lifetime'
    call reference%solve_real(2*b,x,status)
    call checked()
    call oracle()
    changed%val(1)=-changed%val(1)
    call reference%refactor_real(changed,.true.,status,reused)
    if(status%code/=FORTSPARSE_NOT_POSITIVE_DEFINITE) error stop 'failed pivot accepted'
    call reference%solve_real(b,x,status)
    if(status%code/=FORTSPARSE_NOT_FACTORED) error stop 'failed refactor retained ready state'
    changed%val=2*a%val
    call reference%refactor_real(changed,.true.,status,reused)
    call checked()
    if(.not.reused.or.any(user/=17)) error stop 'failed refactor recovery/user lifetime'
    x=b
    x(1)=ieee_value(1._dp,ieee_quiet_nan)
    call reference%solve_real_inplace(x,status)
    if(status%code/=FORTSPARSE_SINGULAR) error stop 'nonfinite residual accepted'
    reference%options%host_memory_limit_bytes=refined%host_current_bytes-1
    call reference%factor_real(a,.true.,status)
    if(status%code/=FORTSPARSE_RESOURCE_LIMIT) error stop 'tight refinement cap accepted'
    call reference%factor_stats(stats,status)
    call checked()
    if(stats%budget_required_bytes<=stats%budget_available_bytes) error stop 'missing cap diagnostics'
    reference%options%host_memory_limit_bytes=refined%host_peak_bytes
    call reference%factor_real(a,.true.,status)
    call checked()
    call reference%solve_real(b,x,status)
    call checked()
    call oracle()
    reference%options%ordering=FORTSPARSE_ORDER_AMD
    call reference%factor_real(a,.true.,status)
    if(status%code/=FORTSPARSE_BACKEND_UNAVAILABLE) error stop 'opaque ordering accepted'
    reference%options%ordering=FORTSPARSE_ORDER_DEFAULT
    reference%options%host_memory_limit_bytes=0
    reference%options%device=0
    call reference%factor_real(a,.true.,status)
    if(status%code/=FORTSPARSE_BACKEND_UNAVAILABLE) error stop 'GPU refinement accepted'
    reference%options%device=-1
    ! Dense 64x64 SPD exercises the accounted supernodal path, independently
    ! manufactured RHS, changed values and changed graph fallback.
    allocate(rows(64*64),cols(64*64),values(64*64),lb(64),lx(64))
    k=0
    do j=1,64
        do i=1,64
            k=k+1
            rows(k)=i
            cols(k)=j
            values(k)=0.01_dp
            if(i==j) values(k)=2.0_dp
        end do
    end do
    call csc_from_triplet(64,64,rows,cols,values,large,status)
    call checked()
    blocked%options%symmetry=FORTSPARSE_SYM_SPD
    blocked%options%nthreads=1
    call blocked%factor_real(large,.true.,status)
    call checked()
    if(.not.blocked%blocked) error stop 'blocked fixture did not use panels'
    lb=2.63_dp
    call blocked%solve_real(lb,lx,status)
    call checked()
    if(maxval(abs(lx-1._dp))>1e-12_dp) error stop 'blocked solution'
    call blocked%factor_stats(stats,status)
    call checked()
    blocked%options%host_memory_limit_bytes=stats%host_peak_bytes
    call blocked%factor_real(large,.true.,status)
    call checked()
    call blocked%refactor_real(a,.true.,status,reused)
    call checked()
    if(reused) error stop 'changed pattern incorrectly reused'
    call blocked%solve_real(b,x,status)
    call checked()
    call oracle()
    call blocked%refactor_real(a,.false.,status,reused)
    call checked()
    if(reused) error stop 'refinement toggle must rebuild'
    call blocked%solve_real(b,x,status)
    call checked()
    if(blocked%stats%backward_error/=-1._dp) error stop 'disabled refinement retained estimate'
    call blocked%refactor_real(a,.true.,status,reused)
    call checked()
    if(reused) error stop 'enabling refinement must rebuild'
contains
    subroutine checked()
        if(.not.status_ok(status)) error stop status%msg
    end subroutine
    subroutine oracle()
        real(dp) :: error
        error=maxval(abs(b-matmul(dense,x)))/(maxval(sum(abs(dense),dim=2))*maxval(abs(x))+maxval(abs(b)))
        if(error>32*epsilon(1._dp)) error stop 'independent original backward error'
    end subroutine
end program
