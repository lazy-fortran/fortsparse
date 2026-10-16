program test_native_refinement_allocations
    use fortsparse
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
    type(sparse_solver_t) :: solver
    type(csc_t) :: a
    type(fortsparse_status_t) :: status
    type(sparse_factor_stats_t) :: before,after
    integer(c_int64_t) :: count,bytes,expected
    real(dp) :: rhs(2),x(2),dx(2),zero(2),primal(2),block(2,2),result(2,2)
    call csc_from_triplet(2,2,[1,2,1,2],[1,1,2,2],[1e-12_dp,1._dp,1._dp,1._dp],a,status)
    call checked(status)
    solver%backend_id=FORTSPARSE_BACKEND_NATIVE
    solver%options%pivot_tolerance=0
    call sparse_factor(solver,a,status)
    call checked(status)
    call sparse_factor_stats(solver,before,status)
    call checked(status)
    solver%refine=.true.
    call begin()
    call sparse_factor(solver,a,status)
    call finish(count,bytes)
    call checked(status)
    expected=int(storage_size(0)/8,c_int64_t)*7+8_c_int64_t*8
    if(bytes/=expected) error stop 'measured original CSC/scratch allocations differ from planned payload'
    call sparse_factor_stats(solver,after,status)
    call checked(status)
    if(after%host_current_bytes-before%host_current_bytes/=bytes) error stop 'reported added payload misses allocations'
    rhs=[1._dp+1e-12_dp,2._dp]
    call begin()
    call sparse_solve(solver,rhs,x,status)
    call finish(count,bytes)
    call checked(status)
    if(count/=0.or.bytes/=0) error stop 'native refined solve performed a linked heap allocation'
    if(maxval(abs(x-1._dp))>1e-11_dp) error stop 'allocation-free correction values'
    primal=1.0_dp
    zero=0.0_dp
    call begin()
    call sparse_solve_jvp(solver,a,primal,zero,dx,status)
    call finish(count,bytes)
    call checked(status)
    if(count/=0.or.bytes/=0) error stop 'native refined JVP performed a linked heap allocation'
    if(maxval(abs(dx+primal))>1e-11_dp) error stop 'allocation-free JVP values'
    block(:,1)=rhs
    block(:,2)=2*rhs
    call begin()
    call sparse_solve_block_op(solver,block,result,FORTSPARSE_OP_TRANSPOSE,status)
    call finish(count,bytes)
    call checked(status)
    if(count/=0.or.bytes/=0) error stop 'native refined block performed a linked heap allocation'
    if(maxval(abs(result(:,2)-2*primal))>1e-11_dp) error stop 'allocation-free block values'
    print '(a,i8)', 'measured retained refinement allocation bytes:',expected
    call sparse_free(solver)
contains
    subroutine checked(s)
        type(fortsparse_status_t), intent(in) :: s
        if(.not.status_ok(s)) error stop s%msg
    end subroutine
end program test_native_refinement_allocations
