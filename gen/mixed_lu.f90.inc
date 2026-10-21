    subroutine factor@S@(self,a,status,storage,column_order,device,memory_limit_bytes)
        class(mixed_lu@S@_t), intent(inout) :: self
        type(@CSC@), intent(in) :: a
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: storage,column_order(:),device
        integer(i8), intent(in), optional :: memory_limit_bytes
        integer :: selected,n,stat
        real(dp) :: condition
        integer(i8) :: cap,reserved,factor_cap
        call self%free()
        cap=0
        if (present(memory_limit_bytes)) cap=memory_limit_bytes
        if (cap<0) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'negative mixed LU cap')
            return
        end if
        if (.not.csc_is_valid(a)) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'invalid mixed LU CSC')
            return
        end if
        n=a%ncol
        if (n<1.or.a%nrow/=n) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'mixed LU needs nonempty square CSC')
            return
        end if
        ! Matrix copy and all retained checker vectors are reserved beside
        ! low-level factor peak, including narrowing replacement overlap.
        reserved=int(a%nnz,i8)*(@VB@_i8+storage_size(0)/8)+ &
            (int(n,i8)+1_i8)*(storage_size(0)/8)+ &
            int(n,i8)*(@VB@_i8*2+8_i8)
        factor_cap=0
        if (cap>0) then
            if (reserved>=cap) then
                call status_set(status,FORTSPARSE_RESOURCE_LIMIT,'mixed LU checker reservation')
                return
            end if
            factor_cap=cap-reserved
        end if
        selected=FP64
        if (present(storage)) selected=storage
        if (present(device)) then
            if (device/=-1) then
                call status_set(status,FORTSPARSE_BACKEND_UNAVAILABLE, &
                    'mixed LU currently supports CPU device -1 only')
                return
            end if
        end if
        if (selected/=FP32.and.selected/=FP64) then
            call status_set(status,FORTSPARSE_BACKEND_UNAVAILABLE, &
                'mixed LU supports FP32 or FP64 factor storage')
            return
        end if
        call self%factors%factor(a,status,column_permutation=column_order, &
            memory_limit_bytes=factor_cap)
        if (status%code/=FORTSPARSE_OK) return
        n=a%ncol
        allocate(self%residual(n),self%delta(n),self%denominator(n),stat=stat)
        if (stat/=0) then
            call self%free()
            call status_set(status,FORTSPARSE_OUT_OF_MEMORY,'mixed LU scratch')
            return
        end if
        allocate(self%original%col_ptr(size(a%col_ptr)), &
            self%original%row_idx(size(a%row_idx)),self%original%val(size(a%val)),stat=stat)
        if (stat/=0) then
            call self%free()
            call status_set(status,FORTSPARSE_OUT_OF_MEMORY,'mixed LU original CSC')
            return
        end if
        self%original%nrow=a%nrow
        self%original%ncol=a%ncol
        self%original%nnz=a%nnz
        self%original%col_ptr=a%col_ptr
        self%original%row_idx=a%row_idx
        self%original%val=a%val
        call native_original_norms(n,a%col_ptr,a%row_idx,a%val,self%delta, &
            self%norm(1),self%norm(2),status)
        if (status%code/=FORTSPARSE_OK) then
            call self%free()
            return
        end if
        ! Estimate both infinity-norm inverse operators before narrowing.
        call estimate@S@(self,0,condition,status)
        self%rcond(1)=condition
        if (status%code==FORTSPARSE_OK) then
            call estimate@S@(self,1,condition,status)
            self%rcond(2)=condition
        end if
        if (status%code/=FORTSPARSE_OK) then
            call self%free()
            return
        end if
        self%storage=FP64
        ! Adaptive precision is conservative: ill-conditioned factors stay dp.
        if (selected==FP32.and.minval(self%rcond)>100*epsilon(1.0_sp)) then
            call self%factors%compress(status)
            if (status%code==FORTSPARSE_OK) self%storage=FP32
            if (status%code==FORTSPARSE_BACKEND_UNAVAILABLE) &
                call status_set(status,FORTSPARSE_OK,'range requires FP64 factors')
            if (status%code/=FORTSPARSE_OK) then
                call self%free()
                return
            end if
        end if
    end subroutine factor@S@

    subroutine estimate@S@(self,operation,rcond,status)
        class(mixed_lu@S@_t), intent(inout) :: self
        integer, intent(in) :: operation
        real(dp), intent(out) :: rcond
        type(fortsparse_status_t), intent(out) :: status
        integer :: k,j,old,n,forward,reverse
        real(dp) :: inverse_norm,previous,magnitude
        rcond=0
        n=self%original%ncol
        forward=0
        reverse=2
        if (operation/=0) then
            forward=2
            reverse=0
        end if
        self%delta=1.0_dp/real(n,dp)
        inverse_norm=0
        old=0
        do k=1,5
            call self%factors%solve(self%delta,status,reverse)
            if (status%code/=FORTSPARSE_OK) return
            previous=inverse_norm
            inverse_norm=sum(abs(self%delta))
            if (k>1.and.inverse_norm<=previous) then
                inverse_norm=previous
                exit
            end if
            do j=1,n
                magnitude=abs(self%delta(j))
                if (magnitude>0) then
                    self%delta(j)=self%delta(j)/magnitude
                else
                    self%delta(j)=1
                end if
            end do
            call self%factors%solve(self%delta,status,forward)
            if (status%code/=FORTSPARSE_OK) return
            j=maxloc(abs(self%delta),dim=1)
            if (j==old) exit
            old=j
            self%delta=0
            self%delta(j)=1
        end do
        rcond=0
        if (inverse_norm>0.and.ieee_is_finite(inverse_norm)) &
            rcond=min(1.0_dp,(1.0_dp/self%norm(operation+1))/inverse_norm)
    end subroutine estimate@S@

    subroutine solve@S@(self,b,x,report,status,operation)
        class(mixed_lu@S@_t), intent(inout) :: self
        @SCALAR@, intent(in) :: b(:)
        @SCALAR@, intent(out) :: x(:)
        type(mixed_solve_report_t), intent(out) :: report
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: operation
        integer :: op,index,k
        real(dp) :: previous
        report=mixed_solve_report_t()
        op=0
        if (present(operation)) op=operation
        call status_set(status,FORTSPARSE_NOT_FACTORED,'mixed LU is not ready')
        if (.not.allocated(self%residual)) return
        call status_set(status,FORTSPARSE_INVALID_MATRIX,'mixed LU dimensions or op')
        if (size(b)/=self%original%ncol.or.size(x)/=size(b)) return
        if (op<0.or.op>2) return
        if (.not.all(@FINITE@(b))) return
        index=1
        if (op/=0) index=2
        report%storage=self%storage
        report%rcond_estimate=self%rcond(index)
        x=b
        call self%factors%solve(x,status,op)
        if (status%code/=FORTSPARSE_OK) return
        previous=huge(1.0_dp)
        do k=0,10
            call native_residual_error(self%original%col_ptr,self%original%row_idx, &
                self%original%val,b,x,op,self%norm(index),self%residual, &
                report%backward_error,self%denominator,report%componentwise_error)
            report%steps=k
            if (report%backward_error<0.or.report%componentwise_error<0) exit
            if (report%backward_error<=32*epsilon(1.0_dp)) then
                report%converged=.true.
                return
            end if
            if (k==10.or.report%backward_error>=previous) exit
            previous=report%backward_error
            self%delta=self%residual
            call self%factors%solve(self%delta,status,op)
            if (status%code/=FORTSPARSE_OK) return
            x=@UPDATE@(x,@MINUS_ONE@,self%delta)
        end do
        call status_set(status,FORTSPARSE_SINGULAR, &
            'mixed LU refinement failed; refactor with FP64 storage')
    end subroutine solve@S@

    subroutine block@S@(self,b,x,reports,statuses,operation)
        class(mixed_lu@S@_t), intent(inout) :: self
        @SCALAR@, intent(in) :: b(:,:)
        @SCALAR@, intent(out) :: x(:,:)
        type(mixed_solve_report_t), intent(out) :: reports(:)
        type(fortsparse_status_t), intent(out) :: statuses(:)
        integer, intent(in), optional :: operation
        integer :: j
        reports=mixed_solve_report_t()
        statuses%code=FORTSPARSE_INVALID_MATRIX
        statuses%msg='mixed LU block dimensions'
        if (size(b,1)/=size(x,1).or.size(b,2)/=size(x,2)) return
        if (size(reports)/=size(b,2).or.size(statuses)/=size(b,2)) return
        do j=1,size(b,2)
            call self%solve(b(:,j),x(:,j),reports(j),statuses(j),operation)
        end do
    end subroutine block@S@

    subroutine apply@S@(self,rhs,status,operation)
        class(mixed_lu@S@_t), intent(inout) :: self
        @SCALAR@, intent(inout) :: rhs(:)
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: operation
        ! Allocation-free approximate inverse seam for Krylov preconditioners.
        call self%factors%solve(rhs,status,operation)
    end subroutine apply@S@

    subroutine factor_batch@S@(systems,matrices,statuses,storage)
        type(mixed_lu@S@_t), intent(inout) :: systems(:)
        type(@CSC@), intent(in) :: matrices(:)
        type(fortsparse_status_t), intent(out) :: statuses(:)
        integer, intent(in), optional :: storage
        integer :: j,n
        n=size(systems)
        statuses%code=FORTSPARSE_INVALID_MATRIX
        statuses%msg='mixed LU factor batch dimensions'
        if (size(matrices)/=n.or.size(statuses)/=n) return
        !$omp parallel do
        do j=1,n
            call systems(j)%factor(matrices(j),statuses(j),storage)
        end do
    end subroutine factor_batch@S@

    subroutine batch@S@(systems,offsets,b,x,reports,statuses)
        type(mixed_lu@S@_t), intent(inout) :: systems(:)
        integer, intent(in) :: offsets(:)
        @SCALAR@, intent(in) :: b(:)
        @SCALAR@, intent(out) :: x(:)
        type(mixed_solve_report_t), intent(out) :: reports(:)
        type(fortsparse_status_t), intent(out) :: statuses(:)
        integer :: j,n
        n=size(systems)
        reports=mixed_solve_report_t()
        statuses%code=FORTSPARSE_INVALID_MATRIX
        statuses%msg='mixed LU batch offsets or dimensions'
        if (size(offsets)/=n+1.or.size(reports)/=n.or.size(statuses)/=n) return
        if (size(b)/=size(x)) return
        if (offsets(1)/=1.or.offsets(n+1)/=size(b)+1) return
        if (any(offsets<1).or.any(offsets>size(b)+1)) return
        if (any(offsets(2:)<offsets(:n))) return
        !$omp parallel do
        do j=1,n
            call systems(j)%solve(b(offsets(j):offsets(j+1)-1), &
                x(offsets(j):offsets(j+1)-1),reports(j),statuses(j))
        end do
    end subroutine batch@S@

    subroutine jacobi@S@(systems,offsets,b,x,statuses)
        ! The contexts are independently factored diagonal blocks. Full-system
        ! off-block entries are deliberately absent from this approximate inverse.
        type(mixed_lu@S@_t), intent(inout) :: systems(:)
        integer, intent(in) :: offsets(:)
        @SCALAR@, intent(in) :: b(:)
        @SCALAR@, intent(out) :: x(:)
        type(fortsparse_status_t), intent(out) :: statuses(:)
        integer :: j,n
        n=size(systems)
        statuses%code=FORTSPARSE_INVALID_MATRIX
        statuses%msg='block Jacobi offsets or dimensions'
        if (size(offsets)/=n+1.or.size(statuses)/=n) return
        if (size(b)/=size(x)) return
        if (offsets(1)/=1.or.offsets(n+1)/=size(b)+1) return
        if (any(offsets<1).or.any(offsets>size(b)+1)) return
        if (any(offsets(2:)<offsets(:n))) return
        x=b
        !$omp parallel do
        do j=1,n
            call systems(j)%apply(x(offsets(j):offsets(j+1)-1),statuses(j))
        end do
    end subroutine jacobi@S@

    integer(i8) function bytes@S@(self) result(payload)
        class(mixed_lu@S@_t), intent(in) :: self
        payload=self%factors%memory_bytes()
        if (allocated(self%residual)) payload=payload+ &
            int(size(self%residual),i8)*(@VB@_i8*2+8_i8)
        if (allocated(self%original%val)) payload=payload+ &
            int(size(self%original%val),i8)*@VB@_i8+ &
            (size(self%original%row_idx,kind=i8)+size(self%original%col_ptr,kind=i8))* &
            int(storage_size(0)/8,i8)
    end function bytes@S@

    subroutine release@S@(self)
        class(mixed_lu@S@_t), intent(inout) :: self
        type(@CSC@) :: empty
        call self%factors%free()
        self%original=empty
        if (allocated(self%residual)) deallocate(self%residual)
        if (allocated(self%delta)) deallocate(self%delta)
        if (allocated(self%denominator)) deallocate(self%denominator)
        self%norm=0
        self%rcond=0
        self%storage=FP64
    end subroutine release@S@
