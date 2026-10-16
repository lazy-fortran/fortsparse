program test_fortsparse_refined_memory_plan
    use fortsparse
    use fortsparse_memory, only: memory_budget_t, plan_natural_spd_memory, refinement_payload_bytes
    implicit none
    type(csc_t) :: a
    type(sparse_solver_t) :: solver
    type(sparse_memory_plan_t) :: plan,ordinary
    type(sparse_factor_stats_t) :: before,after
    type(fortsparse_status_t) :: status
    type(memory_budget_t) :: budget
    integer, parameter :: maxn=17
    integer :: n,kind,backend,nz,nnz,ii,jj,kk,rr(maxn*maxn),cc(maxn*maxn),batch
    real(dp) :: dense(maxn,maxn),vv(maxn*maxn),exact(maxn),rhs(maxn),x(maxn)
    real(dp) :: block(maxn,8),answer(maxn,8)
    logical :: graph(maxn,maxn),fill(maxn,maxn)
    integer(i8) :: checker,base,peak,value
    integer, parameter :: batches(4)=[0,1,8,23]
    ! Checked shared quota helper must not wrap or overwrite output on failure.
    value=37_i8
    call refinement_payload_bytes(huge(0_i8),0_i8,8_i8,value,status)
    call expect(status%code==FORTSPARSE_INDEX_OVERFLOW.and.value==37,'quota endpoint overflow')
    call refinement_payload_bytes(0_i8,huge(0_i8)/4,16_i8,value,status)
    call expect(status%code==FORTSPARSE_INDEX_OVERFLOW.and.value==37,'quota byte overflow')
    call refinement_payload_bytes(-1_i8,0_i8,8_i8,value,status)
    call expect(status%code==FORTSPARSE_INVALID_MATRIX.and.value==37,'negative quota dimensions')
    call refinement_payload_bytes(0_i8,0_i8,3_i8,value,status)
    call expect(status%code==FORTSPARSE_INVALID_MATRIX.and.value==37,'unsupported value width')
    call refinement_payload_bytes(3_i8,9_i8,16_i8,value,status)
    call checked()
    call expect(value==292_i8,'independent complex quota')
    ! The general adapter adds one real n-vector (componentwise denominator).
    call refinement_payload_bytes(3_i8,9_i8,16_i8,value,status,componentwise=.true.)
    call checked()
    call expect(value==292_i8+8_i8*3_i8,'independent componentwise complex quota')
    value=37
    call refinement_payload_bytes(huge(0_i8)/16,0_i8,8_i8,value,status, &
        componentwise=.true.)
    call expect(status%code==FORTSPARSE_INDEX_OVERFLOW.and.value==37, &
        'componentwise quota overflow')
    do kind=1,4
        select case(kind)
        case(1)
            n=1
        case(2)
            n=2
        case(3)
            n=9
        case default
            n=17
        end select
        graph=.false.
        do ii=1,n
            graph(ii,ii)=.true.
        end do
        do ii=2,n
            if(kind==3) then
                graph(1,ii)=.true.
                graph(ii,1)=.true.
            else
                graph(ii-1,ii)=.true.
                graph(ii,ii-1)=.true.
            end if
        end do
        if(kind==4) then
            do ii=1,n-4
                graph(ii,ii+4)=.true.
                graph(ii+4,ii)=.true.
            end do
        end if
        ! Independent dense boolean elimination counts natural-order L fill.
        fill=graph
        do kk=1,n
            do jj=kk+1,n
                do ii=kk+1,n
                    if(fill(ii,kk).and.fill(jj,kk)) fill(ii,jj)=.true.
                end do
            end do
        end do
        nz=0
        do jj=1,n
            nz=nz+count(fill(jj:n,jj))
        end do
        nnz=count(graph(:n,:n))
        dense=0
        kk=0
        do jj=1,n
            do ii=1,n
                if(.not.graph(ii,jj)) cycle
                kk=kk+1
                rr(kk)=ii
                cc(kk)=jj
                vv(kk)=-0.125_dp
                if(ii==jj) vv(kk)=2._dp+real(n,dp)
                dense(ii,jj)=vv(kk)
            end do
        end do
        call csc_from_triplet(n,n,rr(:nnz),cc(:nnz),vv(:nnz),a,status)
        call checked()
        checker=4_i8*(n+1+nnz)+8_i8*(nnz+2*n)
        base=12_i8*nz+4_i8*nnz+36_i8*n+8_i8
        peak=base+checker
        call refinement_payload_bytes(int(n,i8),int(nnz,i8),8_i8,value,status)
        call checked()
        call expect(value==checker,'independent real quota')
        do backend=1,2
            solver%backend_id=FORTSPARSE_BACKEND_NATIVE
            solver%options%memory_policy=FORTSPARSE_MEMORY_LOW
            if(backend==2) then
                solver%backend_id=FORTSPARSE_BACKEND_NATIVE_REFERENCE
                solver%options%memory_policy=FORTSPARSE_MEMORY_AUTO
            end if
            solver%options%symmetry=FORTSPARSE_SYM_SPD
            solver%options%host_memory_limit_bytes=0_i8
            solver%refine=.false.
            call sparse_analyze(solver,a,ordinary,status)
            call checked()
            call expect(ordinary%peak_host_bytes==base,'unrefined exact numbers')
            solver%refine=.true.
            do batch=1,size(batches)
                call sparse_analyze(solver,a,plan,status,nrhs=batches(batch))
                call checked()
                call expect(plan%fits.and.plan%prediction_quality==MEMORY_PLAN_EXACT_LAYOUT,'refined exact quality')
                call expect(plan%factor_host_bytes==8_i8*nz,'unchanged factor payload')
                call expect(plan%persistent_host_bytes==base-8_i8*nz+checker,'original CSC and two vectors')
                call expect(plan%peak_host_bytes==peak.and.plan%symbolic_nnz==nz,'independent refined peak/fill')
                call expect(plan%max_rhs_batch==batches(batch),'bounded many-RHS plan')
            end do
            call sparse_analyze(solver,a,plan,status,nrhs=-1)
            call expect(status%code==FORTSPARSE_INVALID_MATRIX,'negative RHS count')
            solver%options%host_memory_limit_bytes=peak-1_i8
            call sparse_analyze(solver,a,plan,status)
            call expect(status%code==FORTSPARSE_RESOURCE_LIMIT.and..not.plan%fits,'one-byte-under prediction')
            call sparse_factor(solver,a,status)
            call expect(status%code==FORTSPARSE_RESOURCE_LIMIT.and..not.solver%factored,'one-byte-under factor')
            call sparse_destroy(solver)
            solver%options%host_memory_limit_bytes=peak
            call sparse_analyze(solver,a,plan,status,nrhs=8)
            call checked()
            call sparse_factor(solver,a,status)
            call checked()
            call sparse_factor_stats(solver,before,status)
            call checked()
            call expect(before%host_current_bytes==peak.and.before%host_peak_bytes==peak,'predicted exact allocation')
            call sparse_analyze(solver,a,plan,status)
            call expect(status%code==FORTSPARSE_RESOURCE_LIMIT.and.solver%factored,'retained overlap rejects safely')
            call sparse_factor_stats(solver,after,status)
            call checked()
            call expect(after%host_current_bytes==before%host_current_bytes,'planning preserves live bytes')
            call expect(after%host_peak_bytes==before%host_peak_bytes,'planning forecast is not live allocation')
            do ii=1,n
                exact(ii)=real(ii,dp)/real(n+1,dp)
            end do
            rhs(:n)=matmul(dense(:n,:n),exact(:n))
            call sparse_solve(solver,rhs(:n),x(:n),status)
            call checked()
            call expect(maxval(abs(x(:n)-exact(:n)))<1e-13_dp,'exact-cap independently manufactured solve')
            do jj=1,8
                block(:n,jj)=jj*rhs(:n)
            end do
            call sparse_solve_block(solver,block(:n,:),answer(:n,:),status)
            call checked()
            call expect(maxval(abs(answer(:n,8)-8*exact(:n)))<1e-12_dp,'many RHS exact-cap solve')
            solver%options%host_memory_limit_bytes=2*peak
            call sparse_analyze(solver,a,plan,status)
            call checked()
            call expect(plan%fits.and.plan%peak_host_bytes==peak,'retained full-layout overlap admitted')
            call sparse_destroy(solver)
        end do
        ! Forecast checker+planning scratch and final-factor admission must not
        ! charge arrays that this standalone analysis did not allocate.
        budget=memory_budget_t()
        call budget%reserve(16_i8,status)
        call checked()
        budget%limit=16_i8+checker+12_i8*n-1_i8
        call plan_natural_spd_memory(a,plan,status,budget,refine=.true.)
        call expect(status%code==FORTSPARSE_RESOURCE_LIMIT,'checker/symbolic overlap gate')
        call expect(budget%current==16_i8.and.budget%peak==16_i8,'forecast overlap not allocated')
        budget%limit=16_i8+peak-1_i8
        call plan_natural_spd_memory(a,plan,status,budget,refine=.true.)
        call expect(status%code==FORTSPARSE_RESOURCE_LIMIT,'retained final layout gate')
        call expect(budget%current==16_i8.and.budget%peak==16_i8+12_i8*n,'only actual graph scratch charged')
        budget%limit=16_i8+peak
        call plan_natural_spd_memory(a,plan,status,budget,refine=.true.)
        call checked()
        call expect(plan%fits.and.budget%current==16_i8,'successful forecast leaves ledger live bytes')
    end do
    ! A structural layout prediction is not a positive-definiteness proof.
    a%val(1)=-abs(a%val(1))
    call plan_natural_spd_memory(a,plan,status,refine=.true.)
    call checked()
    call expect(plan%fits,'layout should not claim an SPD proof')
    print '(a)', 'exact refined compact SPD plans: 4 patterns, native LOW/reference, RHS0/1/8/23, caps and live overlap'
contains
    subroutine checked()
        if(.not.status_ok(status)) error stop status%msg
    end subroutine
    subroutine expect(condition,message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if(.not.condition) error stop message
    end subroutine
end program
