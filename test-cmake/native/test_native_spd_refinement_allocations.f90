program test_native_spd_refinement_allocations
    use fortsparse
    use fortsparse_native, only: native_backend_t
    use, intrinsic :: iso_c_binding, only: c_int64_t
    implicit none
    interface
        subroutine begin() bind(C,name='fsparse_refinement_alloc_begin')
        end subroutine
        subroutine finish(count,bytes) bind(C,name='fsparse_refinement_alloc_end')
            import c_int64_t
            integer(c_int64_t), intent(out) :: count,bytes
        end subroutine
    end interface
    type(native_backend_t) :: factor
    type(sparse_solver_t) :: solver
    type(csc_t) :: a
    type(fortsparse_status_t) :: status
    integer(c_int64_t) :: count,bytes,baseline_bytes
    real(dp) :: b(2),x(2),rhs(2,2),answer(2,2),zero(2),primal(2)
    logical :: reused
    call csc_from_triplet(2,2,[1,2,1,2],[1,1,2,2],[4._dp,1._dp,1._dp,3._dp],a,status)
    call checked()
    factor%options%symmetry=FORTSPARSE_SYM_SPD
    factor%reference_only=.true.
    factor%options%memory_policy=FORTSPARSE_MEMORY_LOW
    call begin()
    call factor%factor_real(a,.false.,status)
    call finish(count,baseline_bytes)
    call checked()
    call begin()
    call factor%factor_real(a,.true.,status)
    call finish(count,bytes)
    call checked()
    if(bytes-baseline_bytes/=92) error stop 'measured SPD refinement payload differs from accounting'
    factor%lv(factor%lp(1))=factor%lv(factor%lp(1))*(1+1e-4_dp)
    b=[2._dp,-5._dp]
    call begin()
    call factor%solve_real(b,x,status)
    call finish(count,bytes)
    call checked()
    if(count/=0.or.bytes/=0) error stop 'SPD correction allocated'
    if(factor%stats%refinement_steps<1) error stop 'allocation fixture did not correct'
    if(maxval(abs(x-[1._dp,-2._dp]))>1e-12_dp) error stop 'SPD correction values'
    rhs(:,1)=b
    rhs(:,2)=2*b
    call begin()
    call factor%solve_real_block_op(rhs,answer,2,status)
    call finish(count,bytes)
    call checked()
    if(count/=0.or.bytes/=0) error stop 'SPD block correction allocated'
    if(maxval(abs(answer(:,2)-2*x))>1e-12_dp) error stop 'SPD block correction values'
    call begin()
    call factor%refactor_real(a,.true.,status,reused)
    call finish(count,bytes)
    call checked()
    if(.not.reused.or.count/=0.or.bytes/=0) error stop 'SPD retained refactor allocated'
    solver%backend_id=FORTSPARSE_BACKEND_NATIVE
    solver%options%symmetry=FORTSPARSE_SYM_SPD
    solver%refine=.true.
    call sparse_factor(solver,a,status)
    call checked()
    primal=[1._dp,-2._dp]
    zero=0
    call begin()
    call sparse_solve_jvp(solver,a,primal,zero,x,status)
    call finish(count,bytes)
    call checked()
    if(count/=0.or.bytes/=0) error stop 'SPD refined public JVP allocated'
    if(maxval(abs(x+primal))>1e-12_dp) error stop 'SPD refined public JVP values'
contains
    subroutine checked()
        if(.not.status_ok(status)) error stop status%msg
    end subroutine
end program
