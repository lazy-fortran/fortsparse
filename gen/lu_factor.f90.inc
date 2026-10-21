! Generated indexed sparse schedule from gen/lu_factor.f90.inc; DO NOT EDIT.
module @MODULE@
    use fortsparse_kinds, only: dp, sp, i8
    use fortsparse_csc, only: @CSC@, csc_is_valid
    use fortsparse_status, only: fortsparse_status_t, status_set, FORTSPARSE_OK, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_SINGULAR, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_OUT_OF_MEMORY, FORTSPARSE_RESOURCE_LIMIT, FORTSPARSE_INDEX_OVERFLOW, &
        FORTSPARSE_BACKEND_UNAVAILABLE
    use @LEAVES@, only: divide_value => @DIVIDE@, &
        scatter_update => @SCATTER@, gather_update => @GATHER@
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
@REAL_ONLY_BEGIN@
    use, intrinsic :: iso_c_binding, only: c_size_t
@REAL_ONLY_END@
    @EXTRA_USE@
    use fortsparse_memory, only: memory_checked_add, memory_checked_mul
@REAL_PANEL_BEGIN@
    use fortsparse_spd_leaves, only: spd_panel_scatter
@REAL_PANEL_END@
    implicit none
    private
    public :: lu_prune_column
@REAL_PANEL_BEGIN@
    public :: lu_panel_prefix
@REAL_PANEL_END@
    integer, parameter, public :: LU_NORMAL=0, LU_TRANSPOSE=1, LU_ADJOINT=2
    integer(i8), parameter :: ib=int(storage_size(0)/8,i8), vb=@VALUE_BYTES@_i8
@REAL_ONLY_BEGIN@
    interface
        subroutine discard_payload(payload,bytes) &
            bind(C,name='fsparse_native_lu_discard_payload')
            import :: c_size_t
            type(*), intent(inout) :: payload(*)
            integer(c_size_t), value :: bytes
        end subroutine discard_payload
    end interface
@REAL_ONLY_END@
    type, public :: @TYPE@
        integer :: n=0, panel_width=1
        integer(i8) :: factor_nnz=0_i8, peak_memory_bytes=0_i8
        integer(i8) :: failed_required_bytes=0_i8, failed_available_bytes=0_i8
        logical :: last_reused_pattern=.false.
        logical, private :: ready=.false., low=.false., pruned=.false.
        integer(i8), private :: limit=0_i8
        integer, private :: ln=0, un=0
        integer, allocatable, private :: lp(:), up(:), lr(:), ur(:), p(:), pinv(:), q(:)
        integer, allocatable, private :: ap(:), ai(:), mark(:), stack(:), edge(:), reach(:)
        integer, allocatable, private :: prune_end(:)
        @SCALAR@, allocatable, private :: lv(:), uv(:), work(:), correction(:)
        @LOW_SCALAR@, allocatable, private :: lv_low(:), uv_low(:)
        @SCALAR@, allocatable, private :: column_work(:)
        real(dp), allocatable, private :: row_norm(:)
@REAL_PANEL_BEGIN@
        real(dp), allocatable, private :: panel(:,:), panel_coefficient(:)
        integer, allocatable, private :: panel_mask(:)
@REAL_PANEL_END@
    contains
        procedure :: factor, factor_raw
        procedure, private :: solve_vector, solve_block
        generic :: solve => solve_vector, solve_block
        procedure :: free
        procedure :: memory_bytes, copy_column_order, compress
    end type @TYPE@
contains
    subroutine copy_column_order(self,n,col_ptr,row_idx,order,matched)
        ! Metadata only: reuse a validated order from successful same-pattern
        ! factors without duplicating the original pattern in the adapter.
        class(@TYPE@), intent(in) :: self
        integer, intent(in) :: n,col_ptr(:),row_idx(:)
        integer, intent(out) :: order(:)
        logical, intent(out) :: matched
        integer :: i
        matched=.false.
        if (.not.self%ready .or. n/=self%n .or. size(order)/=n) return
        if (.not.allocated(self%ap)) return
        if (.not.allocated(self%ai)) return
        if (.not.allocated(self%q)) return
        if (size(col_ptr)/=size(self%ap) .or. size(row_idx)/=size(self%ai)) return
        do i=1,size(col_ptr)
            if (col_ptr(i)/=self%ap(i)) return
        end do
        do i=1,size(row_idx)
            if (row_idx(i)/=self%ai(i)) return
        end do
        order=self%q
        matched=.true.
    end subroutine copy_column_order

    subroutine factor(self,a,status,column_permutation,pivot_threshold, &
        memory_limit_bytes,low_memory,pruning,row_preference,panel_width)
        class(@TYPE@), intent(inout) :: self
        type(@CSC@), intent(in) :: a
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: column_permutation(:)
        real(dp), intent(in), optional :: pivot_threshold
        integer(i8), intent(in), optional :: memory_limit_bytes
        logical, intent(in), optional :: low_memory,pruning
        integer, intent(in), optional :: row_preference(:), panel_width
        if (allocated(self%lv_low)) call self%free()
        self%ready=.false.
        self%factor_nnz=0_i8
        self%failed_required_bytes=0_i8
        self%failed_available_bytes=0_i8
        self%last_reused_pattern=.false.
        self%peak_memory_bytes=self%memory_bytes()
        if (.not.csc_is_valid(a) .or. a%nrow/=a%ncol) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU requires valid square CSC')
            return
        end if
        call self%factor_raw(a%ncol,a%col_ptr,a%row_idx,a%val,status, &
            column_permutation,pivot_threshold,memory_limit_bytes,low_memory,pruning,row_preference,panel_width)
    end subroutine factor

    ! Threshold partial row pivoting of D*A(:,q). Sparse triangular reach is
    ! Gilbert--Peierls style; no numerical or symbolic dense n-by-n storage.
    ! Caller input is borrowed, never copied into another CSC container.
    subroutine factor_raw(self,n,col_ptr,row_idx,values,status, &
        column_permutation,pivot_threshold,memory_limit_bytes,low_memory,pruning,row_preference,panel_width)
        class(@TYPE@), intent(inout) :: self
        integer, intent(in) :: n,col_ptr(:),row_idx(:)
        @SCALAR@, intent(in) :: values(:)
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: column_permutation(:)
        real(dp), intent(in), optional :: pivot_threshold
        integer(i8), intent(in), optional :: memory_limit_bytes
        logical, intent(in), optional :: low_memory,pruning
        integer, intent(in), optional :: row_preference(:), panel_width
        integer :: j,i,k,t,e,row,col,pivot,count,lc,uc,ierr,width
@REAL_PANEL_BEGIN@
        integer :: panel_first, panel_last, lanes, lane, active
@REAL_PANEL_END@
        integer(i8) :: cap,fixed
        real(dp) :: threshold,best,magnitude
        @SCALAR@ :: value,diagonal
        logical :: same,use_pruning,exact_growth,lower_hint,upper_hint,attempt_slack

        if (allocated(self%lv_low)) call self%free()
        self%ready=.false.
        self%factor_nnz=0_i8
        self%failed_required_bytes=0_i8
        self%failed_available_bytes=0_i8
        self%last_reused_pattern=.false.
        self%peak_memory_bytes=self%memory_bytes()
        width=1
        if (present(panel_width)) width=panel_width
        if (width/=1 .and. width/=8) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU panel width must be 1 or 8')
            return
        end if
        if (width==8) then
            if (@COMPLEX@) then
                call status_set(status,FORTSPARSE_BACKEND_UNAVAILABLE, &
                    'LU target panels currently support real arithmetic only')
                return
            end if
            if (present(low_memory)) then
                if (low_memory) then
                    call status_set(status,FORTSPARSE_BACKEND_UNAVAILABLE, &
                        'LU target panels do not support LOW_MEMORY')
                    return
                end if
            end if
        end if
        threshold=1.0_dp
        if (present(pivot_threshold)) threshold=pivot_threshold
        cap=0_i8
        if (present(memory_limit_bytes)) cap=memory_limit_bytes
        if (.not.ieee_is_finite(threshold)) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU pivot threshold must be finite')
            return
        end if
        if (threshold<0 .or. threshold>1 .or. cap<0) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU invalid pivot threshold or memory cap')
            return
        end if
        if (n<0 .or. n==huge(n) .or. size(values)==huge(0)) then
            call status_set(status,FORTSPARSE_INDEX_OVERFLOW,'LU invalid or overflowing raw dimensions')
            return
        end if
        if (size(col_ptr)/=n+1 .or. size(row_idx)/=size(values)) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU raw CSC shape mismatch')
            return
        end if
        if (col_ptr(1)/=1 .or. col_ptr(n+1)/=size(values)+1) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU raw CSC pointer endpoints')
            return
        end if
        do j=1,n
            if (col_ptr(j)<1 .or. col_ptr(j+1)<col_ptr(j) .or. col_ptr(j+1)>size(values)+1) then
                call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU raw CSC pointer bounds')
                return
            end if
        end do
        do j=1,n
            row=0
            do e=col_ptr(j),col_ptr(j+1)-1
                if (row_idx(e)<=row .or. row_idx(e)>n .or. .not.@FINITE@(values(e))) then
                    call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU needs sorted unique finite CSC entries')
                    return
                end if
                row=row_idx(e)
            end do
        end do
        if (present(column_permutation)) then
            if (size(column_permutation)/=n) then
                call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU column permutation size mismatch')
                return
            end if
        end if
        if (present(row_preference)) then
            if (size(row_preference)/=n) then
                call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU row preference size mismatch')
                return
            end if
        end if
        use_pruning=.false.
        if(present(pruning)) use_pruning=pruning
        same=.false.
        if (allocated(self%ap)) then
            if (self%n==n .and. size(self%ai)==size(values)) then
                same=.true.
                do j=1,n+1
                    if (self%ap(j)/=col_ptr(j)) same=.false.
                end do
                do e=1,size(values)
                    if (self%ai(e)/=row_idx(e)) same=.false.
                end do
            else
                call self%free()
            end if
        end if
        if (cap>0 .and. self%memory_bytes()>cap) then
            call self%free()
            same=.false.
        end if
        self%n=n
        self%limit=cap
        self%low=.false.
        if (present(low_memory)) self%low=low_memory
        self%peak_memory_bytes=self%memory_bytes()
        self%panel_width=width
@REAL_PANEL_BEGIN@
        if (width==1) then
            if (allocated(self%panel)) deallocate(self%panel)
            if (allocated(self%panel_coefficient)) deallocate(self%panel_coefficient)
            if (allocated(self%panel_mask)) deallocate(self%panel_mask)
        end if
@REAL_PANEL_END@
        self%pruned=.false.
        if(.not.use_pruning) then
            if(allocated(self%prune_end)) deallocate(self%prune_end)
        end if
        if (.not.allocated(self%ap)) then
            fixed=ib*(10_i8*n+3_i8+int(size(values),i8))+(2_i8*vb+8_i8)*n
            if(use_pruning) fixed=fixed+ib*n
            if (cap>0 .and. fixed>cap) then
                self%failed_required_bytes=fixed
                self%failed_available_bytes=max(0_i8,cap-self%memory_bytes())
                call status_set(status,FORTSPARSE_RESOURCE_LIMIT,'LU fixed workspace exceeds memory cap')
                return
            end if
            allocate(self%lp(n+1),self%up(n+1),self%p(n),self%pinv(n),self%q(n), &
                self%ap(n+1),self%ai(size(values)),self%mark(n),self%stack(n),self%edge(n), &
                self%reach(n),self%work(n),self%correction(n),self%row_norm(n),stat=ierr)
            if (ierr/=0) then
                fixed=self%memory_bytes()
                call self%free()
                self%peak_memory_bytes=fixed
                call status_set(status,FORTSPARSE_OUT_OF_MEMORY,'LU fixed workspace allocation failed')
                return
            end if
        end if
        self%peak_memory_bytes=max(self%peak_memory_bytes,self%memory_bytes())
        if(use_pruning) then
            if(.not.allocated(self%prune_end)) then
                fixed=ib*n
                if(cap>0_i8) then
                    if(self%memory_bytes()+fixed>cap) then
                        self%failed_required_bytes=fixed
                        self%failed_available_bytes=max(0_i8,cap-self%memory_bytes())
                        call status_set(status,FORTSPARSE_RESOURCE_LIMIT, &
                            'LU pruning workspace exceeds memory cap including reuse')
                        return
                    end if
                end if
                allocate(self%prune_end(n),stat=ierr)
                if(ierr/=0) then
                    call status_set(status,FORTSPARSE_OUT_OF_MEMORY, &
                        'LU pruning workspace allocation failed')
                    return
                end if
            end if
            self%prune_end=0
            self%pruned=.true.
        end if
        self%peak_memory_bytes=max(self%peak_memory_bytes,self%memory_bytes())
        self%ap=col_ptr
        self%ai=row_idx
@REAL_PANEL_BEGIN@
        if (width==8) then
            call panel_storage(self,status)
            if (status%code/=FORTSPARSE_OK) return
        end if
@REAL_PANEL_END@
        self%last_reused_pattern=same
        self%pinv=0
        self%mark=0
        do j=1,n
            col=j
            if (present(column_permutation)) col=column_permutation(j)
            if (col<1 .or. col>n) then
                call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU column permutation index out of range')
                return
            end if
            if (self%mark(col)/=0) then
                call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU column permutation is not bijective')
                return
            end if
            self%mark(col)=1
            self%q(j)=col
        end do
        self%mark=0
        ! Borrow an original-column -> original-row hint during factor only.
        ! Use existing bounded marks; no retained permutation or extra payload.
        if (present(row_preference)) then
            do j=1,n
                row=row_preference(j)
                if (row<1 .or. row>n) then
                    self%last_reused_pattern=.false.
                    call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU row preference index out of range')
                    return
                end if
                if (self%mark(row)/=0) then
                    self%last_reused_pattern=.false.
                    call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU row preference is not bijective')
                    return
                end if
                self%mark(row)=1
            end do
        end if
        self%mark=0
        self%row_norm=0
        do e=1,size(values)
            row=row_idx(e)
            self%row_norm(row)=max(self%row_norm(row),@NORM@(values(e)))
        end do
        if (any(self%row_norm==0.0_dp)) then
            call status_set(status,FORTSPARSE_SINGULAR,'LU structurally or numerically zero row')
            return
        end if
        exact_growth=.false.
        factor_attempt: do
        ! A cap failure caused by this attempt's spare capacity gets exactly
        ! one retry with exact growth. The old numerical factors are already
        ! invalid; only their storage is released, never live entry factors to
        ! evade the initial cap gate. Pattern/order/scaling workspace is kept.
        lower_hint=.false.
        upper_hint=.false.
        call status_set(status,FORTSPARSE_OK,'')
        self%ln=0
        self%un=0
        self%lp(1)=1
        self%up(1)=1
        call order_row_singletons(self,col_ptr,row_idx,values)
        self%up(1)=1
        self%pinv=0
        self%mark=0
        factor_columns: do j=1,n
@REAL_PANEL_BEGIN@
            if (width==8) then
                if (mod(j-1,8)==0) then
                    panel_first=j
                    panel_last=j+min(8,n-j+1)-1
                    lanes=panel_last-panel_first+1
                    count=0
                    do lane=1,lanes
                        col=self%q(panel_first+lane-1)
                        do e=col_ptr(col),col_ptr(col+1)-1
                            call visit(self,row_idx(e),j,count)
                        end do
                    end do
                    do t=1,count
                        row=self%reach(t)
                        self%panel(:,row)=0.0_dp
                        self%panel_mask(row)=0
                    end do
                    do lane=1,lanes
                        col=self%q(panel_first+lane-1)
                        do e=col_ptr(col),col_ptr(col+1)-1
                            row=row_idx(e)
                            value=divide_value(values(e),self%row_norm(row))
                            if (.not.ieee_is_finite(value) .or. &
                                (value==0.0_dp .and. values(e)/=0.0_dp)) then
                                call status_set(status,FORTSPARSE_SINGULAR, &
                                    'LU row scaling exceeds numerical range')
                                return
                            end if
                            self%panel(lane,row)=value
                            self%panel_mask(row)=ibset(self%panel_mask(row),lane-1)
                        end do
                    end do
                    ! Empty L need not be allocated (diagonal / singleton input).
                    if (self%ln>0) then
                        call lu_panel_prefix(self%panel,self%panel_mask, &
                            self%panel_coefficient,lanes,self%reach,count, &
                            self%pinv,self%lp,self%lr,self%lv)
                    end if
                    ! Union used generation j, also the first ordinary visit.
                    do t=1,count
                        self%mark(self%reach(t))=0
                    end do
                end if
            end if
@REAL_PANEL_END@
            col=self%q(j)
            count=0
            do e=col_ptr(col),col_ptr(col+1)-1
                call visit(self,row_idx(e),j,count)
            end do
@REAL_PANEL_BEGIN@
            if (width==8) then
                lane=j-panel_first+1
                do t=1,count
                    row=self%reach(t)
                    self%work(row)=self%panel(lane,row)
                end do
            else
@REAL_PANEL_END@
            do e=col_ptr(col),col_ptr(col+1)-1
                row=row_idx(e)
                value=divide_value(values(e),scalar_from_real(self%row_norm(row)))
                if (.not.@FINITE@(value) .or. (value==@ZERO@ .and. values(e)/=@ZERO@)) then
                    call status_set(status,FORTSPARSE_SINGULAR,'LU row scaling exceeds numerical range')
                    return
                end if
                self%work(row)=value
            end do
            do t=count,1,-1
                row=self%reach(t)
                k=self%pinv(row)
                if (k==0) cycle
                value=self%work(row)
                ! L columns contain distinct rows; value is copied before scatter.
                e=self%lp(k)
                if(self%lp(k+1)>e) then
                    call scatter_update(self%work,self%lr(e:),self%lv(e:), &
                        self%lp(k+1)-e,value)
                end if
            end do
@REAL_PANEL_BEGIN@
            end if
@REAL_PANEL_END@
            best=0
            pivot=0
            lc=0
            uc=1
            do t=1,count
                row=self%reach(t)
                if (.not.@FINITE@(self%work(row))) then
                    call status_set(status,FORTSPARSE_SINGULAR,'LU nonfinite arithmetic during elimination')
                    return
                end if
                if (self%pinv(row)==0) then
                    if (self%work(row)/=@ZERO@) lc=lc+1
                    magnitude=@NORM@(self%work(row))
                    if (magnitude>best) then
                        best=magnitude
                        pivot=row
                    end if
                else if (self%work(row)/=@ZERO@) then
                    uc=uc+1
                end if
            end do
            if (pivot==0) then
                call status_set(status,FORTSPARSE_SINGULAR,'LU zero pivot after sparse elimination')
                return
            end if
            row=self%p(j)
            if(row>0) then
                ! A structural singleton row has no entries in later columns,
                ! so this pivot creates no Schur update or numerical fill.
                if(self%work(row)==@ZERO@) then
                    call status_set(status,FORTSPARSE_SINGULAR,'LU zero singleton pivot')
                    return
                end if
                pivot=row
            else
                row=self%q(j)
                if (present(row_preference)) row=row_preference(self%q(j))
                if (self%pinv(row)==0 .and. self%mark(row)==j) then
                    magnitude=@NORM@(self%work(row))
                    if (magnitude>0 .and. magnitude>=threshold*best) pivot=row
                end if
            end if
            diagonal=self%work(pivot)
            ! A valid pivot is one of the nonzero unpivoted entries counted above.
            lc=lc-1
            if (int(self%ln,i8)+lc>=huge(0) .or. int(self%un,i8)+uc>=huge(0)) then
                call status_set(status,FORTSPARSE_INDEX_OVERFLOW,'LU factor index overflow')
                return
            end if
            call reserve(self,.true.,self%ln+lc,status,exact_growth,lower_hint)
            if (status%code/=FORTSPARSE_OK) exit factor_columns
            call reserve(self,.false.,self%un+uc,status,exact_growth,upper_hint)
            if (status%code/=FORTSPARSE_OK) exit factor_columns
            do t=count,1,-1
                row=self%reach(t)
                k=self%pinv(row)
                if (k/=0) then
                    if (self%work(row)==@ZERO@) cycle
                    self%un=self%un+1
                    self%ur(self%un)=k
                    self%uv(self%un)=self%work(row)
                else if (row/=pivot .and. self%work(row)/=@ZERO@) then
                    value=divide_value(self%work(row),diagonal)
                    if (.not.@FINITE@(value) .or. value==@ZERO@) then
                        call status_set(status,FORTSPARSE_SINGULAR,'LU multiplier exceeds numerical range')
                        return
                    end if
                    self%ln=self%ln+1
                    self%lr(self%ln)=row
                    self%lv(self%ln)=value
                    if(self%pruned) self%mark(row)=-j
                end if
            end do
            self%un=self%un+1
            self%ur(self%un)=j
            self%uv(self%un)=diagonal
            self%p(j)=pivot
            self%pinv(pivot)=j
            self%lp(j+1)=self%ln+1
            self%up(j+1)=self%un+1
            if(self%pruned) then
                self%prune_end(j)=self%lp(j+1)
                ! Negative marks certify actual retained children of new pivot j.
                ! Numerical factors keep every entry; only DFS adjacency shrinks.
                if(self%lp(j+1)>self%lp(j)) then
                    do e=self%up(j),self%up(j+1)-2
                        k=self%ur(e)
                        if(self%prune_end(k)<0) cycle
                        if(self%lp(k)==self%lp(k+1)) cycle
                        call lu_prune_column(self%lr,self%lv,self%lp(k),self%lp(k+1), &
                            pivot,self%mark,j,self%prune_end(k))
                    end do
                end if
            end if
@REAL_PANEL_BEGIN@
            if (width==8) then
                if (j<panel_last) then
                    active=0
                    do lane=1,lanes
                        if (lane>j-panel_first+1) active=ibset(active,lane-1)
                        self%panel_coefficient(lane)=self%panel(lane,pivot)
                    end do
                    e=self%lp(j)
                    if (self%lp(j+1)>e) then
                        call spd_panel_scatter(self%panel,self%panel_mask, &
                            self%lr(e:),self%lv(e:),self%lp(j+1)-e, &
                            self%panel_coefficient,lanes,active)
                    end if
                end if
            end if
@REAL_PANEL_END@
        end do factor_columns
        if(status%code==FORTSPARSE_OK) exit factor_attempt
        attempt_slack=.false.
        if(allocated(self%lr)) attempt_slack=lower_hint.and.size(self%lr)>self%ln
        if(allocated(self%ur)) attempt_slack=attempt_slack.or.(upper_hint.and.size(self%ur)>self%un)
        if(status%code/=FORTSPARSE_RESOURCE_LIMIT.or..not.self%low.or.exact_growth.or..not.attempt_slack) return
        ! No simultaneous old/new retry factors: preserve the measured first
        ! attempt peak, then free all invalid numerical payload before retry.
@REAL_ONLY_BEGIN@
        call discard_factors(self)
@REAL_ONLY_END@
        if(allocated(self%lr)) deallocate(self%lr,self%lv)
        if(allocated(self%ur)) deallocate(self%ur,self%uv)
        self%pinv=0
        self%mark=0
        if(self%pruned) self%prune_end=0
        exact_growth=.true.
        self%last_reused_pattern=.false.
        self%failed_required_bytes=0_i8
        self%failed_available_bytes=0_i8
        end do factor_attempt
        do e=1,self%ln
            self%lr(e)=self%pinv(self%lr(e))
        end do
        self%factor_nnz=int(self%ln,i8)+self%un+n
        self%ready=.true.
        call status_set(status,FORTSPARSE_OK,'')
    end subroutine factor_raw

@REAL_PANEL_BEGIN@
    subroutine lu_panel_prefix(panel,masks,coefficient,lanes,reach,count, &
            pinv,lp,rows,values)
        real(dp), intent(inout) :: panel(8,*), coefficient(8)
        integer, intent(inout) :: masks(*)
        integer, intent(in) :: lanes, reach(*), count, pinv(*), lp(*), rows(*)
        real(dp), intent(in) :: values(*)
        integer :: t,row,k,first,active,lane
        do t=count,1,-1
            row=reach(t)
            k=pinv(row)
            if (k==0) cycle
            active=masks(row)
            do lane=1,lanes
                coefficient(lane)=panel(lane,row)
            end do
            first=lp(k)
            if (lp(k+1)>first) then
                call spd_panel_scatter(panel,masks,rows(first:lp(k+1)-1),values(first:lp(k+1)-1), &
                    lp(k+1)-first,coefficient,lanes,active)
            end if
        end do
    end subroutine lu_panel_prefix

    subroutine panel_storage(self,status)
        class(@TYPE@), intent(inout) :: self
        type(fortsparse_status_t), intent(out) :: status
        integer(i8) :: extra,part,total
        integer :: ierr
        extra=0_i8
        part=0_i8
        total=0_i8
        if (.not.allocated(self%panel)) then
            call memory_checked_mul(64_i8,int(self%n,i8),part,status)
            if (status%code/=FORTSPARSE_OK) goto 90
            call memory_checked_add(extra,part,total,status)
            if (status%code/=FORTSPARSE_OK) goto 90
            extra=total
        end if
        if (.not.allocated(self%panel_mask)) then
            call memory_checked_mul(ib,int(self%n,i8),part,status)
            if (status%code/=FORTSPARSE_OK) goto 90
            call memory_checked_add(extra,part,total,status)
            if (status%code/=FORTSPARSE_OK) goto 90
            extra=total
        end if
        if (.not.allocated(self%panel_coefficient)) then
            call memory_checked_add(extra,64_i8,total,status)
            if (status%code/=FORTSPARSE_OK) goto 90
            extra=total
        end if
        call memory_checked_add(self%memory_bytes(),extra,total,status)
        if (status%code/=FORTSPARSE_OK) goto 90
        if (self%limit>0_i8) then
            if (total>self%limit) then
                self%failed_required_bytes=extra
                self%failed_available_bytes=max(0_i8,self%limit-self%memory_bytes())
                call status_set(status,FORTSPARSE_RESOURCE_LIMIT, &
                    'LU target-panel workspace exceeds cap including reuse')
                return
            end if
        end if
        ! Allocations are setup-only and checked against the whole simultaneous
        ! payload first. Partial failures remain owned and exactly observable.
        ierr=0
        if (.not.allocated(self%panel)) allocate(self%panel(8,self%n),stat=ierr)
        if (ierr/=0) goto 80
        if (.not.allocated(self%panel_mask)) allocate(self%panel_mask(self%n),stat=ierr)
        if (ierr/=0) goto 80
        if (.not.allocated(self%panel_coefficient)) allocate(self%panel_coefficient(8),stat=ierr)
        if (ierr/=0) goto 80
        self%peak_memory_bytes=max(self%peak_memory_bytes,self%memory_bytes())
        call status_set(status,FORTSPARSE_OK,'')
        return
80      self%peak_memory_bytes=max(self%peak_memory_bytes,self%memory_bytes())
        call status_set(status,FORTSPARSE_OUT_OF_MEMORY,'LU target-panel allocation failed')
        return
90      call status_set(status,FORTSPARSE_INDEX_OVERFLOW,'LU target-panel byte overflow')
    end subroutine panel_storage
@REAL_PANEL_END@

    subroutine order_row_singletons(self,col_ptr,row_idx,values)
        ! Standard singleton preprocessing. A row with one remaining original
        ! nonzero fixes its column without mixing independent equations. Such a
        ! row has no later-column entries, so its elimination creates no fill.
        ! Remove that column, update row degrees/XORs, and append the untouched
        ! core in the requested relative order. Every entry is visited O(1) times.
        ! Existing factor/DFS work arrays suffice; no additional allocation.
        class(@TYPE@), intent(inout) :: self
        integer, intent(in) :: col_ptr(:),row_idx(:)
        @SCALAR@, intent(in) :: values(:)
        integer :: row,col,e,head,tail,done,j,child
        self%reach=self%q
        self%mark=0
        self%edge=0
        self%up=0
        self%pinv=0
        self%p=0
        do col=1,self%n
            do e=col_ptr(col),col_ptr(col+1)-1
                if(values(e)==@ZERO@) cycle
                row=row_idx(e)
                self%mark(row)=self%mark(row)+1
                self%edge(row)=ieor(self%edge(row),col)
            end do
        end do
        tail=0
        do row=1,self%n
            if(self%mark(row)/=1) cycle
            tail=tail+1
            self%stack(tail)=row
        end do
        head=1
        done=0
        do while(head<=tail)
            row=self%stack(head)
            head=head+1
            if(self%pinv(row)/=0) cycle
            if(self%mark(row)/=1) cycle
            col=self%edge(row)
            done=done+1
            self%q(done)=col
            self%p(done)=row
            self%up(col)=1
            self%pinv(row)=-1
            do e=col_ptr(col),col_ptr(col+1)-1
                if(values(e)==@ZERO@) cycle
                child=row_idx(e)
                if(self%pinv(child)/=0) cycle
                self%mark(child)=self%mark(child)-1
                self%edge(child)=ieor(self%edge(child),col)
                if(self%mark(child)/=1) cycle
                tail=tail+1
                self%stack(tail)=child
            end do
        end do
        do j=1,self%n
            col=self%reach(j)
            if(self%up(col)/=0) cycle
            done=done+1
            self%q(done)=col
        end do
    end subroutine order_row_singletons

    ! Iterative DFS gives postorder of the current L graph in original row
    ! identities. Each visited numerical workspace entry is cleared once.
    subroutine visit(self,start,generation,count)
        class(@TYPE@), intent(inout) :: self
        integer, intent(in) :: start,generation
        integer, intent(inout) :: count
        integer :: depth,row,k,child,e,last
        if (self%mark(start)==generation) return
        depth=1
        self%stack(1)=start
        self%edge(1)=0
        self%mark(start)=generation
        self%work(start)=@ZERO@
        do while (depth>0)
            row=self%stack(depth)
            k=self%pinv(row)
            if (k/=0) then
                e=self%edge(depth)
                if (e==0) e=self%lp(k)
                last=self%lp(k+1)
                if(self%pruned) last=abs(self%prune_end(k))
                do while (e<last)
                    child=self%lr(e)
                    if (self%mark(child)/=generation) exit
                    e=e+1
                end do
                self%edge(depth)=e+1
                if (e<last) then
                    depth=depth+1
                    self%stack(depth)=child
                    self%edge(depth)=0
                    self%mark(child)=generation
                    self%work(child)=@ZERO@
                    cycle
                end if
            end if
            count=count+1
            self%reach(count)=row
            depth=depth-1
        end do
    end subroutine visit

    ! Graph-only transitive reduction, with a direct retained length-two witness.
    ! Source/pivot are original row identities. Rows are unique and in range.
    ! marks(r)==-generation certifies a retained edge pivot -> r in the new column.
    ! Keep source -> pivot and remove ONLY source -> r with that certificate.
    ! Paired swaps preserve numerical L, whose full span remains in every update.
    subroutine lu_prune_column(rows,values,first,last,pivot,marks,generation,endpoint)
        integer, intent(inout) :: rows(*)
        @SCALAR@, intent(inout) :: values(*)
        integer, intent(in) :: first,last,pivot,marks(*),generation
        integer, intent(inout) :: endpoint
        integer :: e,left,right,row_tmp
        @SCALAR@ :: value_tmp
        logical :: witness
        if(endpoint<0.or.first==last) return
        witness=.false.
        do e=first,last-1
            if(rows(e)==pivot) then
                witness=.true.
                exit
            end if
        end do
        if(.not.witness) return
        left=first
        right=last-1
        do while(left<=right)
            if(rows(left)==pivot) then
                left=left+1
            else if(marks(rows(left))/=-generation) then
                left=left+1
            else
                row_tmp=rows(left)
                rows(left)=rows(right)
                rows(right)=row_tmp
                value_tmp=values(left)
                values(left)=values(right)
                values(right)=value_tmp
                right=right-1
            end if
        end do
        ! All pointer endpoints are one-based; negative means processed once.
        endpoint=-left
    end subroutine lu_prune_column

    subroutine reserve(self,lower,required,status,exact_growth,hinted)
        class(@TYPE@), intent(inout) :: self
        logical, intent(in) :: lower
        integer, intent(in) :: required
        logical, intent(in) :: exact_growth
        logical, intent(inout) :: hinted
        type(fortsparse_status_t), intent(out) :: status
        integer, allocatable :: new_indices(:)
        @SCALAR@, allocatable :: new_values(:)
        integer :: old,capacity,ierr,used
        integer(i8) :: wanted,peak,current
        call status_set(status,FORTSPARSE_OK,'')
        old=0
        if (lower) then
            if (allocated(self%lr)) old=size(self%lr)
            used=self%ln
        else
            if (allocated(self%ur)) old=size(self%ur)
            used=self%un
        end if
        if (required<=old) return
        wanted=required
        if (self%low) then
            ! LOW amortizes prefix copies with at most 1/16 spare capacity at
            ! each growth. Tiny buffers and a cap-recovery retry grow exactly.
            if (.not.exact_growth) wanted=max(wanted,min(int(huge(0),i8), &
                int(old,i8)+max(1_i8,int(old,i8)/16_i8)))
        else
            wanted=max(wanted,min(int(huge(0),i8),max(16_i8,2_i8*old)))
        end if
        current=self%memory_bytes()
        peak=current+wanted*(ib+vb)
        if (self%limit>0 .and. peak>self%limit) then
            wanted=required
            peak=current+wanted*(ib+vb)
        end if
        if (self%limit>0 .and. peak>self%limit) then
            self%failed_required_bytes=wanted*(ib+vb)
            self%failed_available_bytes=max(0_i8,self%limit-current)
            call status_set(status,FORTSPARSE_RESOURCE_LIMIT,'LU factor growth exceeds memory cap including copy')
            return
        end if
        capacity=int(wanted)
        allocate(new_indices(capacity),new_values(capacity),stat=ierr)
        if (ierr/=0) then
            peak=current
            if (allocated(new_indices)) peak=peak+ib*size(new_indices,kind=i8)
            if (allocated(new_values)) peak=peak+vb*size(new_values,kind=i8)
            self%peak_memory_bytes=max(self%peak_memory_bytes,peak)
            call status_set(status,FORTSPARSE_OUT_OF_MEMORY,'LU factor growth allocation failed')
            return
        end if
        self%peak_memory_bytes=max(self%peak_memory_bytes,peak)
        hinted=wanted>required
        if (lower) then
            if (used>0) then
                new_indices(:used)=self%lr(:used)
                new_values(:used)=self%lv(:used)
            end if
@REAL_ONLY_BEGIN@
            if(self%low.and.old>0) then
                ! BOTH prefix copies are finished; the old payload is dead.
                call discard_payload(self%lr,int(int(old,i8)*ib,c_size_t))
                call discard_payload(self%lv,int(int(old,i8)*vb,c_size_t))
            end if
@REAL_ONLY_END@
            call move_alloc(new_indices,self%lr)
            call move_alloc(new_values,self%lv)
        else
            if (used>0) then
                new_indices(:used)=self%ur(:used)
                new_values(:used)=self%uv(:used)
            end if
@REAL_ONLY_BEGIN@
            if(self%low.and.old>0) then
                call discard_payload(self%ur,int(int(old,i8)*ib,c_size_t))
                call discard_payload(self%uv,int(int(old,i8)*vb,c_size_t))
            end if
@REAL_ONLY_END@
            call move_alloc(new_indices,self%ur)
            call move_alloc(new_values,self%uv)
        end if
    end subroutine reserve

    subroutine compress(self,status)
        ! Setup-only narrowing. Factors are built in dp; triangular arithmetic
        ! stays dp. No allocation in retained solves, including N/T/H.
        class(@TYPE@), intent(inout) :: self
        type(fortsparse_status_t), intent(out) :: status
        @LOW_SCALAR@, allocatable :: lower(:),upper(:)
        @SCALAR@, allocatable :: column(:)
        integer :: e,j,ierr
        integer(i8) :: extra,required
        call status_set(status,FORTSPARSE_NOT_FACTORED,'LU factor is not ready')
        if (.not.self%ready) return
        call status_set(status,FORTSPARSE_OK,'')
        if (allocated(self%lv_low)) return
        ! Validate before converting: reject overflow and vanished nonzeros.
        do e=1,self%ln
            if (@NORM@(self%lv(e))>real(huge(0.0_sp),dp)) goto 810
            if (@NORM@(self%lv(e))>0.0_dp) then
                if (@NORM@(self%lv(e))<real(tiny(0.0_sp),dp)) goto 810
            end if
        end do
        do e=1,self%un
            if (@NORM@(self%uv(e))>real(huge(0.0_sp),dp)) goto 810
            if (@NORM@(self%uv(e))>0.0_dp) then
                if (@NORM@(self%uv(e))<real(tiny(0.0_sp),dp)) goto 810
            end if
        end do
        extra=(int(self%ln,i8)+int(self%un,i8))*(vb/2)+int(self%n,i8)*vb
        call memory_checked_add(self%memory_bytes(),extra,required,status)
        if (status%code/=FORTSPARSE_OK) return
        if (self%limit>0_i8) then
            if (required>self%limit) goto 820
        end if
        allocate(lower(self%ln),upper(self%un),column(self%n),stat=ierr)
        if (ierr/=0) then
            call status_set(status,FORTSPARSE_OUT_OF_MEMORY,'LU narrow allocation')
            return
        end if
        lower=self%lv(:self%ln)
        upper=self%uv(:self%un)
        ! Componentwise complex underflow may be harmless, but a vanished
        ! entire U diagonal never defines a valid triangular solve.
        do j=1,self%n
            if (@NORM@(scalar_from_low(upper(self%up(j+1)-1)))==0.0_dp) goto 810
        end do
        self%peak_memory_bytes=max(self%peak_memory_bytes,required)
        call move_alloc(lower,self%lv_low)
        call move_alloc(upper,self%uv_low)
        call move_alloc(column,self%column_work)
        deallocate(self%lv,self%uv)
        return
810     call status_set(status,FORTSPARSE_BACKEND_UNAVAILABLE, &
            'LU factors outside safe fp32 storage range; dp factors retained')
        return
820     self%failed_required_bytes=required
        self%failed_available_bytes=self%limit
        call status_set(status,FORTSPARSE_RESOURCE_LIMIT,'LU narrowing overlap exceeds cap')
    end subroutine compress

    pure function scalar_from_low(value) result(wide)
        @LOW_SCALAR@, intent(in) :: value
        @SCALAR@ :: wide
        wide=value
    end function scalar_from_low

    subroutine expand_column(self,lower,first,last)
        class(@TYPE@), intent(inout) :: self
        logical, intent(in) :: lower
        integer, intent(in) :: first,last
        integer :: k
        if (lower) then
            !$omp simd
            do k=first,last
                self%column_work(k-first+1)=self%lv_low(k)
            end do
        else
            !$omp simd
            do k=first,last
                self%column_work(k-first+1)=self%uv_low(k)
            end do
        end if
    end subroutine expand_column

    subroutine solve_vector(self,rhs,status,operation)
        class(@TYPE@), intent(inout) :: self
        @SCALAR@, intent(inout) :: rhs(:)
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: operation
        integer :: op,j,e,row
        @SCALAR@ :: value,coefficient
        op=LU_NORMAL
        if (present(operation)) op=operation
        if (.not.self%ready) then
            call status_set(status,FORTSPARSE_NOT_FACTORED,'LU factor is not ready')
            return
        end if
        if (size(rhs)/=self%n .or. op<LU_NORMAL .or. op>LU_ADJOINT) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU solve dimension or operation mismatch')
            return
        end if
        do j=1,self%n
            if (.not.@FINITE@(rhs(j))) then
                call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU solve needs finite RHS')
                return
            end if
        end do
        self%correction=@ZERO@
        if (op==LU_NORMAL) then
            do j=1,self%n
                row=self%p(j)
                self%work(j)=divide_value(rhs(row),scalar_from_real(self%row_norm(row)))
            end do
            do j=1,self%n
                e=self%lp(j)
                self%work(j)=corrected_value(self%work(j),self%correction(j))
                self%correction(j)=@ZERO@
                value=self%work(j)
                if(self%lp(j+1)>e) then
                    if (allocated(self%lv_low)) then
                        call expand_column(self,.true.,e,self%lp(j+1)-1)
                        call solve_scatter(self%work,self%correction,self%lr(e:), &
                            self%column_work,self%lp(j+1)-e,value)
                    else
                        call solve_scatter(self%work,self%correction,self%lr(e:), &
                            self%lv(e:),self%lp(j+1)-e,value)
                    end if
                end if
            end do
            do j=self%n,1,-1
                if (allocated(self%uv_low)) then
                    coefficient=self%uv_low(self%up(j+1)-1)
                else
                    coefficient=self%uv(self%up(j+1)-1)
                end if
                self%work(j)=divide_value(corrected_value(self%work(j), &
                    self%correction(j)),coefficient)
                self%correction(j)=@ZERO@
                e=self%up(j)
                value=self%work(j)
                if(self%up(j+1)-1>e) then
                    if (allocated(self%uv_low)) then
                        call expand_column(self,.false.,e,self%up(j+1)-2)
                        call solve_scatter(self%work,self%correction,self%ur(e:), &
                            self%column_work,self%up(j+1)-e-1,value)
                    else
                        call solve_scatter(self%work,self%correction,self%ur(e:), &
                            self%uv(e:),self%up(j+1)-e-1,value)
                    end if
                end if
            end do
            do j=1,self%n
                rhs(self%q(j))=self%work(j)
            end do
        else
            do j=1,self%n
                self%work(j)=rhs(self%q(j))
            end do
            do j=1,self%n
                value=self%work(j)
                e=self%up(j)
                if(self%up(j+1)-1>e) then
                    if (allocated(self%uv_low)) then
                        call expand_column(self,.false.,e,self%up(j+1)-2)
                        call solve_gather(value,self%ur(e:),self%column_work,self%work, &
                            self%up(j+1)-e-1@GATHER_ADJOINT@)
                    else
                        call solve_gather(value,self%ur(e:),self%uv(e:),self%work, &
                            self%up(j+1)-e-1@GATHER_ADJOINT@)
                    end if
                end if
                if (allocated(self%uv_low)) then
                    coefficient=self%uv_low(self%up(j+1)-1)
                else
                    coefficient=self%uv(self%up(j+1)-1)
                end if
                coefficient=adjoint_value(coefficient,op)
                self%work(j)=divide_value(value,coefficient)
            end do
            do j=self%n,1,-1
                value=self%work(j)
                e=self%lp(j)
                if(self%lp(j+1)>e) then
                    if (allocated(self%lv_low)) then
                        call expand_column(self,.true.,e,self%lp(j+1)-1)
                        call solve_gather(value,self%lr(e:),self%column_work,self%work, &
                            self%lp(j+1)-e@GATHER_ADJOINT@)
                    else
                        call solve_gather(value,self%lr(e:),self%lv(e:),self%work, &
                            self%lp(j+1)-e@GATHER_ADJOINT@)
                    end if
                end if
                self%work(j)=value
            end do
            do j=1,self%n
                row=self%p(j)
                rhs(row)=divide_value(self%work(j),scalar_from_real(self%row_norm(row)))
            end do
        end if
        do j=1,self%n
            if (.not.@FINITE@(rhs(j))) then
                call status_set(status,FORTSPARSE_SINGULAR,'LU nonfinite solution')
                return
            end if
        end do
        call status_set(status,FORTSPARSE_OK,'')
    end subroutine solve_vector

    subroutine solve_block(self,rhs,status,operation)
        class(@TYPE@), intent(inout) :: self
        @SCALAR@, intent(inout) :: rhs(:,:)
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: operation
        integer :: j
        if (.not.self%ready) then
            call status_set(status,FORTSPARSE_NOT_FACTORED,'LU factor is not ready')
            return
        end if
        if (size(rhs,1)/=self%n) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU block solve dimension mismatch')
            return
        end if
        if (present(operation)) then
            if (operation<LU_NORMAL .or. operation>LU_ADJOINT) then
                call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU block solve operation mismatch')
                return
            end if
        end if
        call status_set(status,FORTSPARSE_OK,'')
        do j=1,size(rhs,2)
            call self%solve_vector(rhs(:,j),status,operation)
            if (status%code/=FORTSPARSE_OK) return
        end do
    end subroutine solve_block

    pure elemental function scalar_from_real(a) result(value)
        real(dp), intent(in) :: a
        @SCALAR@ :: value
        value=a
    end function scalar_from_real

    pure function adjoint_value(a,op) result(value)
        @SCALAR@, intent(in) :: a
        integer, intent(in) :: op
        @SCALAR@ :: value
        value=a
        if (op==LU_ADJOINT) value=@ADJOINT@
    end function adjoint_value

    integer(i8) function memory_bytes(self) result(bytes)
        class(@TYPE@), intent(in) :: self
        ! Exact allocated-array payload; object descriptors and allocator
        ! overhead remain process-RSS quantities rather than payload claims.
        bytes=0_i8
@MEMORY@
        if(allocated(self%lv_low)) bytes=bytes+(vb/2)*size(self%lv_low,kind=i8)
        if(allocated(self%uv_low)) bytes=bytes+(vb/2)*size(self%uv_low,kind=i8)
        if(allocated(self%column_work)) bytes=bytes+vb*size(self%column_work,kind=i8)
@REAL_PANEL_BEGIN@
        if (allocated(self%panel)) bytes=bytes+8_i8*size(self%panel,kind=i8)
        if (allocated(self%panel_coefficient)) &
            bytes=bytes+8_i8*size(self%panel_coefficient,kind=i8)
        if (allocated(self%panel_mask)) bytes=bytes+ib*size(self%panel_mask,kind=i8)
@REAL_PANEL_END@
    end function memory_bytes

@REAL_ONLY_BEGIN@
    subroutine discard_factors(self)
        class(@TYPE@), intent(inout) :: self
        ! Keep allocatable ownership. Only dead real LOW factors are advised.
        if(.not.self%low) return
        if(allocated(self%lr)) &
            call discard_payload(self%lr,int(size(self%lr,kind=i8)*ib,c_size_t))
        if(allocated(self%lv)) &
            call discard_payload(self%lv,int(size(self%lv,kind=i8)*vb,c_size_t))
        if(allocated(self%ur)) &
            call discard_payload(self%ur,int(size(self%ur,kind=i8)*ib,c_size_t))
        if(allocated(self%uv)) &
            call discard_payload(self%uv,int(size(self%uv,kind=i8)*vb,c_size_t))
    end subroutine discard_factors

@REAL_ONLY_END@
    subroutine free(self)
        class(@TYPE@), intent(inout) :: self
@REAL_ONLY_BEGIN@
        call discard_factors(self)
@REAL_ONLY_END@
@FREE@
        if(allocated(self%lv_low)) deallocate(self%lv_low)
        if(allocated(self%uv_low)) deallocate(self%uv_low)
        if(allocated(self%column_work)) deallocate(self%column_work)
@REAL_PANEL_BEGIN@
        if (allocated(self%panel)) deallocate(self%panel)
        if (allocated(self%panel_coefficient)) deallocate(self%panel_coefficient)
        if (allocated(self%panel_mask)) deallocate(self%panel_mask)
@REAL_PANEL_END@
        self%panel_width=1
        self%n=0
        self%ln=0
        self%un=0
        self%factor_nnz=0_i8
        self%peak_memory_bytes=0_i8
        self%failed_required_bytes=0_i8
        self%failed_available_bytes=0_i8
        self%last_reused_pattern=.false.
        self%ready=.false.
        self%pruned=.false.
        self%limit=0_i8
    end subroutine free
end module @MODULE@
