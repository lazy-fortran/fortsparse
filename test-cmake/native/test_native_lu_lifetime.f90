program test_fortsparse_native_lu_lifetime
    use fortsparse
    use fortsparse_native_lu, only: lu_factor_t
    use fortsparse_native_lu_z, only: lu_factor_z_t
    use, intrinsic :: iso_c_binding, only: c_int,c_int64_t
    implicit none
    interface
        subroutine advice_mode(denied) bind(C,name='fsparse_lu_discard_test_mode')
            import c_int
            integer(c_int), value :: denied
        end subroutine
        integer(c_int) function allocation_begin(nth) &
            bind(C,name='fsparse_lu_alloc_test_begin')
            import c_int
            integer(c_int), value :: nth
        end function
        integer(c_int) function allocation_end() &
            bind(C,name='fsparse_lu_alloc_test_end')
            import c_int
        end function
        integer(c_int64_t) function live_bytes() &
            bind(C,name='fsparse_lu_alloc_test_current')
            import c_int64_t
        end function
        integer(c_int64_t) function peak_bytes() &
            bind(C,name='fsparse_lu_alloc_test_peak')
            import c_int64_t
        end function
        integer(c_int) function advice_calls() &
            bind(C,name='fsparse_lu_discard_test_calls')
            import c_int
        end function
    end interface
    integer, parameter :: ptr(4)=[1,4,7,10],rows(9)=[1,2,3,1,2,3,1,2,3]
    real(dp), parameter :: v(9)=[4._dp,1._dp,.25_dp,1._dp,5._dp,1._dp, &
        .25_dp,1._dp,6._dp]
    integer, parameter :: n=140000
    integer(c_int) :: denied,real_calls
    character(8) :: argument
    denied=0
    call get_command_argument(1,argument)
    if(len_trim(argument)>0) read(argument,*) denied
    call advice_mode(denied)
    call real_copy()
    real_calls=advice_calls()
    if(real_calls==0) error stop 'large real LOW factor never advised'
    call complex_copy()
    if(advice_calls()/=real_calls) error stop 'complex factors unexpectedly advised'
    call allocation_failure()
    print '(a)', 'LOW factors preserve copy, self assignment and scoped lifetime'
contains
    subroutine checked(status)
        type(fortsparse_status_t), intent(in) :: status
        if(.not.status_ok(status)) then
            print *,status%code,trim(status%msg)
            error stop 'unexpected LOW lifetime status'
        end if
    end subroutine checked
    subroutine allocation_failure()
        type(lu_factor_t) :: f
        type(lu_factor_z_t) :: zf
        type(fortsparse_status_t) :: status
        real(dp), parameter :: solution(3)=[1._dp,2._dp,3._dp]
        real(dp) :: rhs(3),dense(3,3)
        complex(dp) :: zrhs(3)
        integer :: op,supported,attempts,nth,seen,kind
        dense=reshape(v,[3,3])
        do kind=1,2
            ! Discover this revision's actual allocation sequence, then fail
            ! every observed ordinal. No stale fixed-array count or byte formula.
            supported=allocation_begin(0_c_int)
            if(supported==0) return
            call allocation_attempt(kind,f,zf,status)
            attempts=allocation_end()
            print '(a,i0,a,i0)', 'Allocation oracle scalar kind ',kind, &
                ': injected ordinals ',attempts
            call checked(status)
            call payload_oracle(kind,f,zf)
            call f%free();call zf%free()
            if(live_bytes()/=0) error stop 'baseline allocation leak'
            do nth=1,attempts
                supported=allocation_begin(int(nth,c_int))
                call allocation_attempt(kind,f,zf,status)
                seen=allocation_end()
                if(status%code/=FORTSPARSE_OUT_OF_MEMORY) &
                    error stop 'injected LU allocation failure status'
                if(seen<nth) error stop 'allocation fault ordinal not reached'
                call payload_oracle(kind,f,zf)
                rhs=1._dp;zrhs=(1._dp,0._dp)
                if(kind==1) then
                    call f%solve(rhs,status)
                else
                    call zf%solve(zrhs,status)
                end if
                if(status%code/=FORTSPARSE_NOT_FACTORED) &
                    error stop 'allocation failure retained usable factor'
                supported=allocation_begin(0_c_int)
                call allocation_attempt(kind,f,zf,status)
                seen=allocation_end()
                call checked(status)
                call payload_oracle(kind,f,zf)
                do op=FORTSPARSE_OP_NORMAL,FORTSPARSE_OP_ADJOINT
                    rhs=matmul(dense,solution);zrhs=cmplx(rhs,0._dp,dp)
                    if(kind==1) then
                        call f%solve(rhs,status,op)
                        call checked(status)
                        if(maxval(abs(rhs-solution))>1e-13_dp) &
                            error stop 'real fault recovery known-x N/T/H'
                    else
                        call zf%solve(zrhs,status,op)
                        call checked(status)
                        if(maxval(abs(zrhs-cmplx(solution,0._dp,dp)))>1e-13_dp) &
                            error stop 'complex fault recovery known-x N/T/H'
                    end if
                end do
                call f%free();call zf%free()
                if(live_bytes()/=0) error stop 'allocation fault/recovery leak'
            end do
        end do
    end subroutine allocation_failure
    subroutine allocation_attempt(kind,f,zf,status)
        integer, intent(in) :: kind
        type(lu_factor_t), intent(inout) :: f
        type(lu_factor_z_t), intent(inout) :: zf
        type(fortsparse_status_t), intent(out) :: status
        if(kind==1) then
            call f%factor_raw(3,ptr,rows,v,status,low_memory=.true.,pruning=.true.)
        else
            call zf%factor_raw(3,ptr,rows,cmplx(v,0._dp,dp),status, &
                low_memory=.true.,pruning=.true.)
        end if
    end subroutine
    subroutine payload_oracle(kind,f,zf)
        integer, intent(in) :: kind
        type(lu_factor_t), intent(in) :: f
        type(lu_factor_z_t), intent(in) :: zf
        integer(i8) :: current,peak
        if(kind==1) then
            current=f%memory_bytes();peak=f%peak_memory_bytes
        else
            current=zf%memory_bytes();peak=zf%peak_memory_bytes
        end if
        if(current/=live_bytes().or.peak/=peak_bytes()) then
            print *,kind,current,live_bytes(),peak,peak_bytes()
            error stop 'external LU live/peak allocation mismatch'
        end if
    end subroutine
    subroutine real_copy()
        type(sparse_solver_t) :: a,b
        type(csc_t) :: matrix
        type(fortsparse_status_t) :: status
        real(dp) :: rhs(n),x(n)
        integer :: i,op,before
        logical :: reused
        matrix%nrow=n;matrix%ncol=n;matrix%nnz=n
        allocate(matrix%col_ptr(n+1),matrix%row_idx(n),matrix%val(n))
        matrix%col_ptr=[(i,i=1,n+1)];matrix%row_idx=[(i,i=1,n)];matrix%val=2._dp
        a%backend_id=FORTSPARSE_BACKEND_NATIVE
        a%options%memory_policy=FORTSPARSE_MEMORY_LOW
        call sparse_factor(a,matrix,status)
        call checked(status)
        ! This oracle owns only allocatable factors. Public sparse_vector
        ! pointer-copy ownership is a separate reported defect, not waived here.
        b=a
        call sparse_free(a)
        b=b
        rhs=2._dp
        do op=FORTSPARSE_OP_NORMAL,FORTSPARSE_OP_ADJOINT
            call sparse_solve_op(b,rhs,x,op,status)
            call checked(status)
            if(maxval(abs(x-1._dp))>1e-14_dp) &
                error stop 'public real copied factor alias'
        end do
        a=b
        matrix%val=4._dp
        call sparse_refactor(a,matrix,status,reused)
        call checked(status)
        call sparse_solve(a,rhs,x,status)
        call checked(status)
        if(maxval(abs(x-.5_dp))>1e-14_dp) error stop 'public real copied refactor'
        call sparse_free(a)
        call sparse_solve(b,rhs,x,status)
        call checked(status)
        if(maxval(abs(x-1._dp))>1e-14_dp) error stop 'public real copy independence'
        call sparse_free(b)
        matrix%val=2._dp
        call real_scoped_copy(a,matrix)
        call sparse_solve(a,rhs,x,status)
        call checked(status)
        if(maxval(abs(x-1._dp))>1e-14_dp) error stop 'real scoped backend finalization'
        call sparse_free(a)
        before=advice_calls()
        a%options%memory_policy=FORTSPARSE_MEMORY_AUTO
        call sparse_factor(a,matrix,status)
        call checked(status)
        call sparse_solve(a,rhs,x,status)
        call checked(status)
        if(maxval(abs(x-1._dp))>1e-14_dp) error stop 'real AUTO known x'
        call sparse_free(a)
        if(advice_calls()/=before) error stop 'real AUTO factors unexpectedly advised'
    end subroutine real_copy
    subroutine real_scoped_copy(result,matrix)
        type(sparse_solver_t), intent(out) :: result
        type(csc_t), intent(in) :: matrix
        type(sparse_solver_t) :: local
        type(fortsparse_status_t) :: status
        local%backend_id=FORTSPARSE_BACKEND_NATIVE
        local%options%memory_policy=FORTSPARSE_MEMORY_LOW
        call sparse_factor(local,matrix,status)
        call checked(status)
        result=local
        ! Existing backend finalization frees local; result must remain usable.
    end subroutine real_scoped_copy
    subroutine complex_copy()
        type(sparse_solver_t) :: a,b
        type(csc_z_t) :: matrix
        type(fortsparse_status_t) :: status
        complex(dp) :: rhs(n),x(n),diagonal
        integer :: i,op
        logical :: reused
        matrix%nrow=n;matrix%ncol=n;matrix%nnz=n
        allocate(matrix%col_ptr(n+1),matrix%row_idx(n),matrix%val(n))
        matrix%col_ptr=[(i,i=1,n+1)];matrix%row_idx=[(i,i=1,n)]
        diagonal=cmplx(2._dp,.125_dp,dp);matrix%val=diagonal
        a%backend_id=FORTSPARSE_BACKEND_NATIVE
        a%options%memory_policy=FORTSPARSE_MEMORY_LOW
        call sparse_factor(a,matrix,status)
        call checked(status)
        ! This oracle owns only allocatable factors. Public sparse_vector
        ! pointer-copy ownership is a separate reported defect, not waived here.
        b=a
        call sparse_free(a)
        b=b
        do op=FORTSPARSE_OP_NORMAL,FORTSPARSE_OP_ADJOINT
            rhs=diagonal
            if(op==FORTSPARSE_OP_ADJOINT) rhs=conjg(diagonal)
            call sparse_solve_op(b,rhs,x,op,status)
            call checked(status)
            if(maxval(abs(x-1._dp))>1e-14_dp) &
                error stop 'public complex copied factor alias'
        end do
        a=b
        matrix%val=2*diagonal
        call sparse_refactor(a,matrix,status,reused)
        call checked(status)
        rhs=diagonal
        call sparse_solve(a,rhs,x,status)
        call checked(status)
        if(maxval(abs(x-.5_dp))>1e-14_dp) error stop 'public complex copied refactor'
        call sparse_free(a)
        call sparse_solve(b,rhs,x,status)
        call checked(status)
        if(maxval(abs(x-1._dp))>1e-14_dp) error stop 'public complex copy independence'
        call sparse_free(b)
        matrix%val=diagonal
        call complex_scoped_copy(a,matrix)
        call sparse_solve(a,rhs,x,status)
        call checked(status)
        if(maxval(abs(x-1._dp))>1e-14_dp) &
                error stop 'complex scoped backend finalization'
        call sparse_free(a)
    end subroutine complex_copy
    subroutine complex_scoped_copy(result,matrix)
        type(sparse_solver_t), intent(out) :: result
        type(csc_z_t), intent(in) :: matrix
        type(sparse_solver_t) :: local
        type(fortsparse_status_t) :: status
        local%backend_id=FORTSPARSE_BACKEND_NATIVE
        local%options%memory_policy=FORTSPARSE_MEMORY_LOW
        call sparse_factor(local,matrix,status)
        call checked(status)
        result=local
    end subroutine complex_scoped_copy
end program test_fortsparse_native_lu_lifetime
