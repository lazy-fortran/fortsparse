program test_general_amd_bounded
    use fortsparse
    use, intrinsic :: iso_c_binding, only: c_int,c_int64_t
    implicit none
    integer, parameter :: n=9
    type(sparse_solver_t) :: solver
    type(csc_t) :: a
    type(csc_z_t) :: z
    type(fortsparse_status_t) :: status
    type(sparse_factor_stats_t) :: stats,fresh,previous
    real(dp) :: dense(n,n),exact(n),rhs(n),x(n),block(n,8),answer(n,8)
    real(dp), pointer :: user(:)
    complex(dp) :: zdense(n,n),zexact(n),zrhs(n),zx(n),zblock(n,8),zanswer(n,8)
    integer :: i,j,before,k,rr(n*n),cc(n*n)
    real(dp) :: vv(n*n)
    complex(dp) :: zz(n*n)
    integer(i8) :: quota,cap,old_live,order_phase
    integer(c_int64_t) :: arena_bytes,output_bytes
    logical :: reused
    interface
        integer(c_int) function available() bind(C,name='fsparse_native_amd_bounded_available')
            import c_int
        end function
        integer(c_int) function calls() bind(C,name='fsparse_test_bounded_amd_calls')
            import c_int
        end function
        integer(c_int) function query(n,nz,ptr,row,low,arena_bytes,output_bytes) &
                bind(C,name='fsparse_native_amd_bounded_size')
            import c_int,c_int64_t
            integer(c_int), value :: n,nz,low
            integer(c_int), intent(in) :: ptr(*),row(*)
            integer(c_int64_t), intent(out) :: arena_bytes,output_bytes
        end function
    end interface
    dense=0
    do i=1,n
        dense(i,i)=4._dp+0.1_dp*i
        dense(i,1)=0.2_dp*i
        exact(i)=sin(real(i,dp))
    end do
    dense(1,1)=0 ! Requires numerical row pivoting, independent of AMD columns.
    dense(1,2)=1
    dense(2,1)=2
    call matrix()
    solver%backend_id=FORTSPARSE_BACKEND_NATIVE
    solver%options%ordering=FORTSPARSE_ORDER_AMD
    solver%options%memory_policy=FORTSPARSE_MEMORY_LOW
    call sparse_factor(solver,a,status)
    if(available()/=1) then
        call expect(status%code==FORTSPARSE_BACKEND_UNAVAILABLE,'audited ABI unavailable status')
        call expect(calls()==0,'unavailable bounded path invoked')
        stop
    end if
    call checked()
    call operations()
    call sparse_factor_stats(solver,fresh,status)
    call checked()
    call expect(fresh%strict_memory_supported.and.fresh%host_peak_bytes>0,'bounded exact phase statistics')
    call expect(fresh%ordering_used==FORTSPARSE_ORDER_AMD,'actual bounded order')
    cap=fresh%host_peak_bytes
    call sparse_destroy(solver)
    solver%options%host_memory_limit_bytes=cap-1_i8
    call sparse_factor(solver,a,status)
    call expect(status%code==FORTSPARSE_RESOURCE_LIMIT.and..not.solver%factored,'one byte under fresh LOW peak')
    call sparse_destroy(solver)
    solver%options%host_memory_limit_bytes=cap
    call sparse_factor(solver,a,status)
    call checked()
    call sparse_factor_stats(solver,stats,status)
    call checked()
    call expect(stats%host_peak_bytes==cap,'fresh exact cap matches tracked allocation peak')
    call operations()
    ! Same-pattern numerical refactor needs only the charged permutation copy.
    call sparse_factor_stats(solver,previous,status)
    old_live=previous%host_current_bytes
    solver%options%host_memory_limit_bytes=old_live+4_i8*n-1
    before=calls()
    call sparse_refactor(solver,a,status,reused)
    call expect(status%code==FORTSPARSE_RESOURCE_LIMIT,'one byte below retained order-copy cap')
    call expect(calls()==before,'same graph cap failure invoked AMD')
    call sparse_factor_stats(solver,stats,status)
    call checked()
    call expect(stats%host_current_bytes==old_live,'copy cap refusal preserves retained payload')
    solver%options%host_memory_limit_bytes=old_live+4_i8*n
    call sparse_refactor(solver,a,status,reused)
    call checked()
    call expect(reused.and.calls()==before,'same pattern reuses bounded column order')
    call sparse_factor_stats(solver,stats,status)
    call checked()
    call expect(stats%host_peak_bytes==old_live+4_i8*n,'retained refactor exact q-copy peak')
    call operations()
    ! Existing old factors must coexist with the full arena on a strategy change.
    solver%options%host_memory_limit_bytes=0
    solver%options%ordering=FORTSPARSE_ORDER_DEFAULT
    call sparse_refactor(solver,a,status,reused)
    call checked()
    call expect(.not.reused,'changed strategy must invalidate reuse')
    call sparse_factor_stats(solver,stats,status)
    old_live=stats%host_current_bytes
    i=query(n,a%nnz,a%col_ptr,a%row_idx,1,arena_bytes,output_bytes)
    call expect(i==0,'independent workspace query')
    order_phase=arena_bytes+output_bytes
    solver%options%ordering=FORTSPARSE_ORDER_AMD
    solver%options%host_memory_limit_bytes=old_live+order_phase-1
    before=calls()
    call sparse_factor(solver,a,status)
    call expect(status%code==FORTSPARSE_RESOURCE_LIMIT,'old factors plus arena exceed cap')
    call expect(calls()==before,'arena refused before vendor call')
    call sparse_factor_stats(solver,stats,status)
    call checked()
    call expect(stats%host_current_bytes==old_live,'ordering failure preserves old allocated arrays')
    call expect(stats%budget_required_bytes==order_phase,'exact complete ordering reservation diagnostic')
    call expect(stats%budget_available_bytes==order_phase-1,'old factors charged in available budget')
    solver%options%host_memory_limit_bytes=0
    call sparse_factor(solver,a,status)
    call checked()
    call operations()
    dense(n,2)=0.03125_dp
    call matrix()
    before=calls()
    call sparse_refactor(solver,a,status,reused)
    call checked()
    call expect(.not.reused.and.calls()==before+1,'changed pattern recomputes bounded AMD')
    call operations()
    ! A retained user vector and original-matrix checker share the same cap.
    solver%refine=.true.
    call sparse_factor(solver,a,status)
    call checked()
    user=>sparse_vector(solver,n)
    call expect(associated(user),'retained user vector')
    user=17
    call sparse_factor_stats(solver,stats,status)
    old_live=stats%host_current_bytes
    quota=4_i8*(n+1+a%nnz)+8_i8*(a%nnz+2*n)
    call expect(old_live>quota+8_i8*n,'checker and user payload counted')
    solver%options%host_memory_limit_bytes=old_live+4_i8*n
    before=calls()
    call sparse_refactor(solver,a,status,reused)
    call checked()
    call expect(reused.and.calls()==before.and.all(user==17),'bounded refinement and user lifetime')
    call operations()
    ! Raw real path borrows input storage and uses the same bounded ordering.
    call sparse_destroy(solver)
    solver%options%host_memory_limit_bytes=0
    before=calls()
    call sparse_factor(solver,n,n,a%nnz,a%col_ptr,a%row_idx,a%val,status)
    call checked()
    call expect(calls()==before+1,'raw real bounded ordering')
    call operations()
    ! Ordinary bounded (non-LOW) uses the recommended larger elbow workspace.
    call sparse_destroy(solver)
    solver%refine=.false.
    solver%options%memory_policy=FORTSPARSE_MEMORY_AUTO
    solver%options%host_memory_limit_bytes=1024_i8*1024_i8
    call sparse_factor(solver,a,status)
    call checked()
    call operations()
    call sparse_factor_stats(solver,stats,status)
    call checked()
    call expect(stats%strict_memory_supported.and.stats%host_peak_bytes>0,'bounded ordinary mode')
    call sparse_destroy(solver)
    ! Complex N/T/H and eight RHS use retained factors and original CSC checks.
    zdense=cmplx(dense,0._dp,dp)
    where(dense/=0._dp) zdense=cmplx(dense,0.125_dp*dense,dp)
    k=0
    do j=1,n
        do i=1,n
            if(zdense(i,j)==(0._dp,0._dp)) cycle
            k=k+1;rr(k)=i;cc(k)=j;zz(k)=zdense(i,j)
        end do
        zexact(j)=cmplx(cos(real(j,dp)),sin(real(j,dp)),dp)
    end do
    call csc_from_triplet(n,n,rr(:k),cc(:k),zz(:k),z,status)
    call checked()
    solver%options%memory_policy=FORTSPARSE_MEMORY_LOW
    solver%options%host_memory_limit_bytes=0
    solver%refine=.true.
    call sparse_factor(solver,z,status)
    call checked()
    call zoperations()
    before=calls()
    call sparse_refactor(solver,z,status,reused)
    call checked()
    call expect(reused.and.calls()==before,'complex same pattern skips AMD')
    call sparse_factor_stats(solver,stats,status)
    call checked()
    solver%options%host_memory_limit_bytes=stats%host_current_bytes+4_i8*n
    call sparse_refactor(solver,z,status,reused)
    call checked()
    call zoperations()
    call sparse_destroy(solver)
    solver%options%host_memory_limit_bytes=0
    call sparse_factor(solver,n,n,z%nnz,z%col_ptr,z%row_idx,z%val,status)
    call checked()
    call zoperations()
    print '(a)', 'bounded real/complex AMD LU: cap/reuse/strategy/refinement/raw and N/T/H block oracles pass'
contains
    subroutine matrix()
        k=0
        do j=1,n
            do i=1,n
                if(dense(i,j)==0._dp) cycle
                k=k+1;rr(k)=i;cc(k)=j;vv(k)=dense(i,j)
            end do
        end do
        call csc_from_triplet(n,n,rr(:k),cc(:k),vv(:k),a,status)
        call checked()
    end subroutine
    subroutine operations()
        integer :: operation,column
        do operation=0,2
            if(operation==0) then
                rhs=matmul(dense,exact)
            else
                rhs=matmul(transpose(dense),exact)
            end if
            call sparse_solve_op(solver,rhs,x,operation,status)
            call checked()
            call expect(maxval(abs(x-exact))<2e-12_dp,'independent real operation solution')
            do column=1,8
                block(:,column)=column*rhs
            end do
            call sparse_solve_block_op(solver,block,answer,operation,status)
            call checked()
            call expect(maxval(abs(answer(:,8)-8*exact))<2e-11_dp,'independent real block solution')
        end do
    end subroutine
    subroutine zoperations()
        integer :: operation,column
        do operation=0,2
            select case(operation)
            case(0)
                zrhs=matmul(zdense,zexact)
            case(1)
                zrhs=matmul(transpose(zdense),zexact)
            case default
                zrhs=matmul(conjg(transpose(zdense)),zexact)
            end select
            call sparse_solve_op(solver,zrhs,zx,operation,status)
            call checked()
            call expect(maxval(abs(zx-zexact))<3e-12_dp,'independent complex operation solution')
            do column=1,8
                zblock(:,column)=column*zrhs
            end do
            call sparse_solve_block_op(solver,zblock,zanswer,operation,status)
            call checked()
            call expect(maxval(abs(zanswer(:,8)-8*zexact))<3e-11_dp,'independent complex block solution')
        end do
    end subroutine
    subroutine expect(condition,message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if(.not.condition) error stop message
    end subroutine
    subroutine checked()
        if(.not.status_ok(status)) error stop status%msg
    end subroutine
end program
