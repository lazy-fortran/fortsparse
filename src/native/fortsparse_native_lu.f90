! Generated indexed sparse schedule from gen/lu_factor.f90.inc; DO NOT EDIT.
module fortsparse_native_lu
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csc, only: csc_t, csc_is_valid
    use fortsparse_status, only: fortsparse_status_t, status_set, FORTSPARSE_OK, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_SINGULAR, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_OUT_OF_MEMORY, FORTSPARSE_RESOURCE_LIMIT, FORTSPARSE_INDEX_OVERFLOW
    use fortsparse_spd_leaves, only: update_value => spd_update, divide_value => spd_divide
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none
    private
    integer, parameter, public :: LU_NORMAL=0, LU_TRANSPOSE=1, LU_ADJOINT=2
    integer(i8), parameter :: ib=int(storage_size(0)/8,i8), vb=8_i8
    type, public :: lu_factor_t
        integer :: n=0
        integer(i8) :: factor_nnz=0_i8, peak_memory_bytes=0_i8
        integer(i8) :: failed_required_bytes=0_i8, failed_available_bytes=0_i8
        logical :: last_reused_pattern=.false.
        logical, private :: ready=.false., low=.false.
        integer(i8), private :: limit=0_i8
        integer, private :: ln=0, un=0
        integer, allocatable, private :: lp(:), up(:), lr(:), ur(:), p(:), pinv(:), q(:)
        integer, allocatable, private :: ap(:), ai(:), mark(:), stack(:), edge(:), reach(:)
        real(dp), allocatable, private :: lv(:), uv(:), work(:)
        real(dp), allocatable, private :: row_norm(:)
    contains
        procedure :: factor, factor_raw
        procedure, private :: solve_vector, solve_block
        generic :: solve => solve_vector, solve_block
        procedure :: free
        procedure :: memory_bytes
    end type lu_factor_t
contains
    subroutine factor(self,a,status,column_permutation,pivot_threshold,memory_limit_bytes,low_memory)
        class(lu_factor_t), intent(inout) :: self
        type(csc_t), intent(in) :: a
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: column_permutation(:)
        real(dp), intent(in), optional :: pivot_threshold
        integer(i8), intent(in), optional :: memory_limit_bytes
        logical, intent(in), optional :: low_memory
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
            column_permutation,pivot_threshold,memory_limit_bytes,low_memory)
    end subroutine factor

    ! Threshold partial row pivoting of D*A(:,q). Sparse triangular reach is
    ! Gilbert--Peierls style; no numerical or symbolic dense n-by-n storage.
    ! Caller input is borrowed, never copied into another CSC container.
    subroutine factor_raw(self,n,col_ptr,row_idx,values,status, &
        column_permutation,pivot_threshold,memory_limit_bytes,low_memory)
        class(lu_factor_t), intent(inout) :: self
        integer, intent(in) :: n,col_ptr(:),row_idx(:)
        real(dp), intent(in) :: values(:)
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: column_permutation(:)
        real(dp), intent(in), optional :: pivot_threshold
        integer(i8), intent(in), optional :: memory_limit_bytes
        logical, intent(in), optional :: low_memory
        integer :: j,i,k,t,e,row,col,pivot,count,lc,uc,ierr
        integer(i8) :: cap,fixed
        real(dp) :: threshold,best,magnitude
        real(dp) :: value,diagonal
        logical :: same

        self%ready=.false.
        self%factor_nnz=0_i8
        self%failed_required_bytes=0_i8
        self%failed_available_bytes=0_i8
        self%last_reused_pattern=.false.
        self%peak_memory_bytes=self%memory_bytes()
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
                if (row_idx(e)<=row .or. row_idx(e)>n .or. .not.ieee_is_finite(values(e))) then
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
        if (.not.allocated(self%ap)) then
            fixed=ib*(10_i8*n+3_i8+int(size(values),i8))+(vb+8_i8)*n
            if (cap>0 .and. fixed>cap) then
                self%failed_required_bytes=fixed
                self%failed_available_bytes=max(0_i8,cap-self%memory_bytes())
                call status_set(status,FORTSPARSE_RESOURCE_LIMIT,'LU fixed workspace exceeds memory cap')
                return
            end if
            allocate(self%lp(n+1),self%up(n+1),self%p(n),self%pinv(n),self%q(n), &
                self%ap(n+1),self%ai(size(values)),self%mark(n),self%stack(n),self%edge(n), &
                self%reach(n),self%work(n),self%row_norm(n),stat=ierr)
            if (ierr/=0) then
                fixed=self%memory_bytes()
                call self%free()
                self%peak_memory_bytes=fixed
                call status_set(status,FORTSPARSE_OUT_OF_MEMORY,'LU fixed workspace allocation failed')
                return
            end if
        end if
        self%peak_memory_bytes=max(self%peak_memory_bytes,self%memory_bytes())
        self%ap=col_ptr
        self%ai=row_idx
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
        self%row_norm=0
        do e=1,size(values)
            row=row_idx(e)
            self%row_norm(row)=max(self%row_norm(row),abs(values(e)))
        end do
        if (any(self%row_norm==0.0_dp)) then
            call status_set(status,FORTSPARSE_SINGULAR,'LU structurally or numerically zero row')
            return
        end if
        self%ln=0
        self%un=0
        self%lp(1)=1
        self%up(1)=1
        do j=1,n
            col=self%q(j)
            count=0
            do e=col_ptr(col),col_ptr(col+1)-1
                call visit(self,row_idx(e),j,count)
            end do
            do e=col_ptr(col),col_ptr(col+1)-1
                row=row_idx(e)
                value=divide_value(values(e),scalar_from_real(self%row_norm(row)))
                if (.not.ieee_is_finite(value) .or. (value==0.0_dp .and. values(e)/=0.0_dp)) then
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
                do e=self%lp(k),self%lp(k+1)-1
                    i=self%lr(e)
                    self%work(i)=update_value(self%work(i),self%lv(e),value)
                end do
            end do
            best=0
            pivot=0
            lc=0
            uc=1
            do t=1,count
                row=self%reach(t)
                if (.not.ieee_is_finite(self%work(row))) then
                    call status_set(status,FORTSPARSE_SINGULAR,'LU nonfinite arithmetic during elimination')
                    return
                end if
                if (self%pinv(row)==0) then
                    magnitude=abs(self%work(row))
                    if (magnitude>best) then
                        best=magnitude
                        pivot=row
                    end if
                else if (self%work(row)/=0.0_dp) then
                    uc=uc+1
                end if
            end do
            if (pivot==0) then
                call status_set(status,FORTSPARSE_SINGULAR,'LU zero pivot after sparse elimination')
                return
            end if
            row=self%q(j)
            if (self%pinv(row)==0 .and. self%mark(row)==j) then
                magnitude=abs(self%work(row))
                if (magnitude>0 .and. magnitude>=threshold*best) pivot=row
            end if
            diagonal=self%work(pivot)
            do t=1,count
                row=self%reach(t)
                if (self%pinv(row)==0 .and. row/=pivot .and. self%work(row)/=0.0_dp) lc=lc+1
            end do
            if (int(self%ln,i8)+lc>=huge(0) .or. int(self%un,i8)+uc>=huge(0)) then
                call status_set(status,FORTSPARSE_INDEX_OVERFLOW,'LU factor index overflow')
                return
            end if
            call reserve(self,.true.,self%ln+lc,status)
            if (status%code/=FORTSPARSE_OK) return
            call reserve(self,.false.,self%un+uc,status)
            if (status%code/=FORTSPARSE_OK) return
            do t=count,1,-1
                row=self%reach(t)
                k=self%pinv(row)
                if (k/=0) then
                    if (self%work(row)==0.0_dp) cycle
                    self%un=self%un+1
                    self%ur(self%un)=k
                    self%uv(self%un)=self%work(row)
                else if (row/=pivot .and. self%work(row)/=0.0_dp) then
                    value=divide_value(self%work(row),diagonal)
                    if (.not.ieee_is_finite(value) .or. value==0.0_dp) then
                        call status_set(status,FORTSPARSE_SINGULAR,'LU multiplier exceeds numerical range')
                        return
                    end if
                    self%ln=self%ln+1
                    self%lr(self%ln)=row
                    self%lv(self%ln)=value
                end if
            end do
            self%un=self%un+1
            self%ur(self%un)=j
            self%uv(self%un)=diagonal
            self%p(j)=pivot
            self%pinv(pivot)=j
            self%lp(j+1)=self%ln+1
            self%up(j+1)=self%un+1
        end do
        do e=1,self%ln
            self%lr(e)=self%pinv(self%lr(e))
        end do
        self%factor_nnz=int(self%ln,i8)+self%un+n
        self%ready=.true.
        call status_set(status,FORTSPARSE_OK,'')
    end subroutine factor_raw

    ! Iterative DFS gives postorder of the current L graph in original row
    ! identities. Each visited numerical workspace entry is cleared once.
    subroutine visit(self,start,generation,count)
        class(lu_factor_t), intent(inout) :: self
        integer, intent(in) :: start,generation
        integer, intent(inout) :: count
        integer :: depth,row,k,child,e
        if (self%mark(start)==generation) return
        depth=1
        self%stack(1)=start
        self%edge(1)=0
        self%mark(start)=generation
        self%work(start)=0.0_dp
        do while (depth>0)
            row=self%stack(depth)
            k=self%pinv(row)
            if (k/=0) then
                e=self%edge(depth)
                if (e==0) e=self%lp(k)
                do while (e<self%lp(k+1))
                    child=self%lr(e)
                    if (self%mark(child)/=generation) exit
                    e=e+1
                end do
                self%edge(depth)=e+1
                if (e<self%lp(k+1)) then
                    depth=depth+1
                    self%stack(depth)=child
                    self%edge(depth)=0
                    self%mark(child)=generation
                    self%work(child)=0.0_dp
                    cycle
                end if
            end if
            count=count+1
            self%reach(count)=row
            depth=depth-1
        end do
    end subroutine visit

    subroutine reserve(self,lower,required,status)
        class(lu_factor_t), intent(inout) :: self
        logical, intent(in) :: lower
        integer, intent(in) :: required
        type(fortsparse_status_t), intent(out) :: status
        integer, allocatable :: new_indices(:)
        real(dp), allocatable :: new_values(:)
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
        if (.not.self%low) wanted=max(wanted,min(int(huge(0),i8),max(16_i8,2_i8*old)))
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
        if (lower) then
            if (used>0) then
                new_indices(:used)=self%lr(:used)
                new_values(:used)=self%lv(:used)
            end if
            call move_alloc(new_indices,self%lr)
            call move_alloc(new_values,self%lv)
        else
            if (used>0) then
                new_indices(:used)=self%ur(:used)
                new_values(:used)=self%uv(:used)
            end if
            call move_alloc(new_indices,self%ur)
            call move_alloc(new_values,self%uv)
        end if
    end subroutine reserve

    subroutine solve_vector(self,rhs,status,operation)
        class(lu_factor_t), intent(inout) :: self
        real(dp), intent(inout) :: rhs(:)
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: operation
        integer :: op,j,e,row
        real(dp) :: value,coefficient
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
            if (.not.ieee_is_finite(rhs(j))) then
                call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU solve needs finite RHS')
                return
            end if
        end do
        if (op==LU_NORMAL) then
            do j=1,self%n
                row=self%p(j)
                self%work(j)=divide_value(rhs(row),scalar_from_real(self%row_norm(row)))
            end do
            do j=1,self%n
                do e=self%lp(j),self%lp(j+1)-1
                    row=self%lr(e)
                    self%work(row)=update_value(self%work(row),self%lv(e),self%work(j))
                end do
            end do
            do j=self%n,1,-1
                self%work(j)=divide_value(self%work(j),self%uv(self%up(j+1)-1))
                do e=self%up(j),self%up(j+1)-2
                    row=self%ur(e)
                    self%work(row)=update_value(self%work(row),self%uv(e),self%work(j))
                end do
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
                do e=self%up(j),self%up(j+1)-2
                    coefficient=adjoint_value(self%uv(e),op)
                    value=update_value(value,coefficient,self%work(self%ur(e)))
                end do
                coefficient=adjoint_value(self%uv(self%up(j+1)-1),op)
                self%work(j)=divide_value(value,coefficient)
            end do
            do j=self%n,1,-1
                value=self%work(j)
                do e=self%lp(j),self%lp(j+1)-1
                    coefficient=adjoint_value(self%lv(e),op)
                    value=update_value(value,coefficient,self%work(self%lr(e)))
                end do
                self%work(j)=value
            end do
            do j=1,self%n
                row=self%p(j)
                rhs(row)=divide_value(self%work(j),scalar_from_real(self%row_norm(row)))
            end do
        end if
        do j=1,self%n
            if (.not.ieee_is_finite(rhs(j))) then
                call status_set(status,FORTSPARSE_SINGULAR,'LU nonfinite solution')
                return
            end if
        end do
        call status_set(status,FORTSPARSE_OK,'')
    end subroutine solve_vector

    subroutine solve_block(self,rhs,status,operation)
        class(lu_factor_t), intent(inout) :: self
        real(dp), intent(inout) :: rhs(:,:)
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
        real(dp) :: value
        value=a
    end function scalar_from_real

    pure function adjoint_value(a,op) result(value)
        real(dp), intent(in) :: a
        integer, intent(in) :: op
        real(dp) :: value
        value=a
        if (op==LU_ADJOINT) value=a
    end function adjoint_value

    integer(i8) function memory_bytes(self) result(bytes)
        class(lu_factor_t), intent(in) :: self
        ! Exact allocated-array payload; object descriptors and allocator
        ! overhead remain process-RSS quantities rather than payload claims.
        bytes=0_i8
        if(allocated(self%lp)) bytes=bytes+ib*size(self%lp,kind=i8)
        if(allocated(self%up)) bytes=bytes+ib*size(self%up,kind=i8)
        if(allocated(self%lr)) bytes=bytes+ib*size(self%lr,kind=i8)
        if(allocated(self%ur)) bytes=bytes+ib*size(self%ur,kind=i8)
        if(allocated(self%p)) bytes=bytes+ib*size(self%p,kind=i8)
        if(allocated(self%pinv)) bytes=bytes+ib*size(self%pinv,kind=i8)
        if(allocated(self%q)) bytes=bytes+ib*size(self%q,kind=i8)
        if(allocated(self%ap)) bytes=bytes+ib*size(self%ap,kind=i8)
        if(allocated(self%ai)) bytes=bytes+ib*size(self%ai,kind=i8)
        if(allocated(self%mark)) bytes=bytes+ib*size(self%mark,kind=i8)
        if(allocated(self%stack)) bytes=bytes+ib*size(self%stack,kind=i8)
        if(allocated(self%edge)) bytes=bytes+ib*size(self%edge,kind=i8)
        if(allocated(self%reach)) bytes=bytes+ib*size(self%reach,kind=i8)
        if(allocated(self%lv)) bytes=bytes+vb*size(self%lv,kind=i8)
        if(allocated(self%uv)) bytes=bytes+vb*size(self%uv,kind=i8)
        if(allocated(self%work)) bytes=bytes+vb*size(self%work,kind=i8)
        if(allocated(self%row_norm)) bytes=bytes+8_i8*size(self%row_norm,kind=i8)
    end function memory_bytes

    subroutine free(self)
        class(lu_factor_t), intent(inout) :: self
        if(allocated(self%lp)) deallocate(self%lp)
        if(allocated(self%up)) deallocate(self%up)
        if(allocated(self%lr)) deallocate(self%lr)
        if(allocated(self%ur)) deallocate(self%ur)
        if(allocated(self%p)) deallocate(self%p)
        if(allocated(self%pinv)) deallocate(self%pinv)
        if(allocated(self%q)) deallocate(self%q)
        if(allocated(self%ap)) deallocate(self%ap)
        if(allocated(self%ai)) deallocate(self%ai)
        if(allocated(self%mark)) deallocate(self%mark)
        if(allocated(self%stack)) deallocate(self%stack)
        if(allocated(self%edge)) deallocate(self%edge)
        if(allocated(self%reach)) deallocate(self%reach)
        if(allocated(self%lv)) deallocate(self%lv)
        if(allocated(self%uv)) deallocate(self%uv)
        if(allocated(self%work)) deallocate(self%work)
        if(allocated(self%row_norm)) deallocate(self%row_norm)
        self%n=0
        self%ln=0
        self%un=0
        self%factor_nnz=0_i8
        self%peak_memory_bytes=0_i8
        self%failed_required_bytes=0_i8
        self%failed_available_bytes=0_i8
        self%last_reused_pattern=.false.
        self%ready=.false.
        self%limit=0_i8
    end subroutine free
end module fortsparse_native_lu
! Generated indexed sparse schedule from gen/lu_factor.f90.inc; DO NOT EDIT.
module fortsparse_native_lu_z
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csc, only: csc_z_t, csc_is_valid
    use fortsparse_status, only: fortsparse_status_t, status_set, FORTSPARSE_OK, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_SINGULAR, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_OUT_OF_MEMORY, FORTSPARSE_RESOURCE_LIMIT, FORTSPARSE_INDEX_OVERFLOW
    use fortsparse_lu_leaves, only: update_value => lu_z_update, divide_value => lu_z_divide
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortsparse_lu_leaves, only: lu_z_finite, lu_z_norm
    implicit none
    private
    integer, parameter, public :: LU_NORMAL=0, LU_TRANSPOSE=1, LU_ADJOINT=2
    integer(i8), parameter :: ib=int(storage_size(0)/8,i8), vb=16_i8
    type, public :: lu_factor_z_t
        integer :: n=0
        integer(i8) :: factor_nnz=0_i8, peak_memory_bytes=0_i8
        integer(i8) :: failed_required_bytes=0_i8, failed_available_bytes=0_i8
        logical :: last_reused_pattern=.false.
        logical, private :: ready=.false., low=.false.
        integer(i8), private :: limit=0_i8
        integer, private :: ln=0, un=0
        integer, allocatable, private :: lp(:), up(:), lr(:), ur(:), p(:), pinv(:), q(:)
        integer, allocatable, private :: ap(:), ai(:), mark(:), stack(:), edge(:), reach(:)
        complex(dp), allocatable, private :: lv(:), uv(:), work(:)
        real(dp), allocatable, private :: row_norm(:)
    contains
        procedure :: factor, factor_raw
        procedure, private :: solve_vector, solve_block
        generic :: solve => solve_vector, solve_block
        procedure :: free
        procedure :: memory_bytes
    end type lu_factor_z_t
contains
    subroutine factor(self,a,status,column_permutation,pivot_threshold,memory_limit_bytes,low_memory)
        class(lu_factor_z_t), intent(inout) :: self
        type(csc_z_t), intent(in) :: a
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: column_permutation(:)
        real(dp), intent(in), optional :: pivot_threshold
        integer(i8), intent(in), optional :: memory_limit_bytes
        logical, intent(in), optional :: low_memory
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
            column_permutation,pivot_threshold,memory_limit_bytes,low_memory)
    end subroutine factor

    ! Threshold partial row pivoting of D*A(:,q). Sparse triangular reach is
    ! Gilbert--Peierls style; no numerical or symbolic dense n-by-n storage.
    ! Caller input is borrowed, never copied into another CSC container.
    subroutine factor_raw(self,n,col_ptr,row_idx,values,status, &
        column_permutation,pivot_threshold,memory_limit_bytes,low_memory)
        class(lu_factor_z_t), intent(inout) :: self
        integer, intent(in) :: n,col_ptr(:),row_idx(:)
        complex(dp), intent(in) :: values(:)
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: column_permutation(:)
        real(dp), intent(in), optional :: pivot_threshold
        integer(i8), intent(in), optional :: memory_limit_bytes
        logical, intent(in), optional :: low_memory
        integer :: j,i,k,t,e,row,col,pivot,count,lc,uc,ierr
        integer(i8) :: cap,fixed
        real(dp) :: threshold,best,magnitude
        complex(dp) :: value,diagonal
        logical :: same

        self%ready=.false.
        self%factor_nnz=0_i8
        self%failed_required_bytes=0_i8
        self%failed_available_bytes=0_i8
        self%last_reused_pattern=.false.
        self%peak_memory_bytes=self%memory_bytes()
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
                if (row_idx(e)<=row .or. row_idx(e)>n .or. .not.lu_z_finite(values(e))) then
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
        if (.not.allocated(self%ap)) then
            fixed=ib*(10_i8*n+3_i8+int(size(values),i8))+(vb+8_i8)*n
            if (cap>0 .and. fixed>cap) then
                self%failed_required_bytes=fixed
                self%failed_available_bytes=max(0_i8,cap-self%memory_bytes())
                call status_set(status,FORTSPARSE_RESOURCE_LIMIT,'LU fixed workspace exceeds memory cap')
                return
            end if
            allocate(self%lp(n+1),self%up(n+1),self%p(n),self%pinv(n),self%q(n), &
                self%ap(n+1),self%ai(size(values)),self%mark(n),self%stack(n),self%edge(n), &
                self%reach(n),self%work(n),self%row_norm(n),stat=ierr)
            if (ierr/=0) then
                fixed=self%memory_bytes()
                call self%free()
                self%peak_memory_bytes=fixed
                call status_set(status,FORTSPARSE_OUT_OF_MEMORY,'LU fixed workspace allocation failed')
                return
            end if
        end if
        self%peak_memory_bytes=max(self%peak_memory_bytes,self%memory_bytes())
        self%ap=col_ptr
        self%ai=row_idx
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
        self%row_norm=0
        do e=1,size(values)
            row=row_idx(e)
            self%row_norm(row)=max(self%row_norm(row),lu_z_norm(values(e)))
        end do
        if (any(self%row_norm==0.0_dp)) then
            call status_set(status,FORTSPARSE_SINGULAR,'LU structurally or numerically zero row')
            return
        end if
        self%ln=0
        self%un=0
        self%lp(1)=1
        self%up(1)=1
        do j=1,n
            col=self%q(j)
            count=0
            do e=col_ptr(col),col_ptr(col+1)-1
                call visit(self,row_idx(e),j,count)
            end do
            do e=col_ptr(col),col_ptr(col+1)-1
                row=row_idx(e)
                value=divide_value(values(e),scalar_from_real(self%row_norm(row)))
                if (.not.lu_z_finite(value) .or. (value==(0.0_dp,0.0_dp) .and. values(e)/=(0.0_dp,0.0_dp))) then
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
                do e=self%lp(k),self%lp(k+1)-1
                    i=self%lr(e)
                    self%work(i)=update_value(self%work(i),self%lv(e),value)
                end do
            end do
            best=0
            pivot=0
            lc=0
            uc=1
            do t=1,count
                row=self%reach(t)
                if (.not.lu_z_finite(self%work(row))) then
                    call status_set(status,FORTSPARSE_SINGULAR,'LU nonfinite arithmetic during elimination')
                    return
                end if
                if (self%pinv(row)==0) then
                    magnitude=lu_z_norm(self%work(row))
                    if (magnitude>best) then
                        best=magnitude
                        pivot=row
                    end if
                else if (self%work(row)/=(0.0_dp,0.0_dp)) then
                    uc=uc+1
                end if
            end do
            if (pivot==0) then
                call status_set(status,FORTSPARSE_SINGULAR,'LU zero pivot after sparse elimination')
                return
            end if
            row=self%q(j)
            if (self%pinv(row)==0 .and. self%mark(row)==j) then
                magnitude=lu_z_norm(self%work(row))
                if (magnitude>0 .and. magnitude>=threshold*best) pivot=row
            end if
            diagonal=self%work(pivot)
            do t=1,count
                row=self%reach(t)
                if (self%pinv(row)==0 .and. row/=pivot .and. self%work(row)/=(0.0_dp,0.0_dp)) lc=lc+1
            end do
            if (int(self%ln,i8)+lc>=huge(0) .or. int(self%un,i8)+uc>=huge(0)) then
                call status_set(status,FORTSPARSE_INDEX_OVERFLOW,'LU factor index overflow')
                return
            end if
            call reserve(self,.true.,self%ln+lc,status)
            if (status%code/=FORTSPARSE_OK) return
            call reserve(self,.false.,self%un+uc,status)
            if (status%code/=FORTSPARSE_OK) return
            do t=count,1,-1
                row=self%reach(t)
                k=self%pinv(row)
                if (k/=0) then
                    if (self%work(row)==(0.0_dp,0.0_dp)) cycle
                    self%un=self%un+1
                    self%ur(self%un)=k
                    self%uv(self%un)=self%work(row)
                else if (row/=pivot .and. self%work(row)/=(0.0_dp,0.0_dp)) then
                    value=divide_value(self%work(row),diagonal)
                    if (.not.lu_z_finite(value) .or. value==(0.0_dp,0.0_dp)) then
                        call status_set(status,FORTSPARSE_SINGULAR,'LU multiplier exceeds numerical range')
                        return
                    end if
                    self%ln=self%ln+1
                    self%lr(self%ln)=row
                    self%lv(self%ln)=value
                end if
            end do
            self%un=self%un+1
            self%ur(self%un)=j
            self%uv(self%un)=diagonal
            self%p(j)=pivot
            self%pinv(pivot)=j
            self%lp(j+1)=self%ln+1
            self%up(j+1)=self%un+1
        end do
        do e=1,self%ln
            self%lr(e)=self%pinv(self%lr(e))
        end do
        self%factor_nnz=int(self%ln,i8)+self%un+n
        self%ready=.true.
        call status_set(status,FORTSPARSE_OK,'')
    end subroutine factor_raw

    ! Iterative DFS gives postorder of the current L graph in original row
    ! identities. Each visited numerical workspace entry is cleared once.
    subroutine visit(self,start,generation,count)
        class(lu_factor_z_t), intent(inout) :: self
        integer, intent(in) :: start,generation
        integer, intent(inout) :: count
        integer :: depth,row,k,child,e
        if (self%mark(start)==generation) return
        depth=1
        self%stack(1)=start
        self%edge(1)=0
        self%mark(start)=generation
        self%work(start)=(0.0_dp,0.0_dp)
        do while (depth>0)
            row=self%stack(depth)
            k=self%pinv(row)
            if (k/=0) then
                e=self%edge(depth)
                if (e==0) e=self%lp(k)
                do while (e<self%lp(k+1))
                    child=self%lr(e)
                    if (self%mark(child)/=generation) exit
                    e=e+1
                end do
                self%edge(depth)=e+1
                if (e<self%lp(k+1)) then
                    depth=depth+1
                    self%stack(depth)=child
                    self%edge(depth)=0
                    self%mark(child)=generation
                    self%work(child)=(0.0_dp,0.0_dp)
                    cycle
                end if
            end if
            count=count+1
            self%reach(count)=row
            depth=depth-1
        end do
    end subroutine visit

    subroutine reserve(self,lower,required,status)
        class(lu_factor_z_t), intent(inout) :: self
        logical, intent(in) :: lower
        integer, intent(in) :: required
        type(fortsparse_status_t), intent(out) :: status
        integer, allocatable :: new_indices(:)
        complex(dp), allocatable :: new_values(:)
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
        if (.not.self%low) wanted=max(wanted,min(int(huge(0),i8),max(16_i8,2_i8*old)))
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
        if (lower) then
            if (used>0) then
                new_indices(:used)=self%lr(:used)
                new_values(:used)=self%lv(:used)
            end if
            call move_alloc(new_indices,self%lr)
            call move_alloc(new_values,self%lv)
        else
            if (used>0) then
                new_indices(:used)=self%ur(:used)
                new_values(:used)=self%uv(:used)
            end if
            call move_alloc(new_indices,self%ur)
            call move_alloc(new_values,self%uv)
        end if
    end subroutine reserve

    subroutine solve_vector(self,rhs,status,operation)
        class(lu_factor_z_t), intent(inout) :: self
        complex(dp), intent(inout) :: rhs(:)
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: operation
        integer :: op,j,e,row
        complex(dp) :: value,coefficient
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
            if (.not.lu_z_finite(rhs(j))) then
                call status_set(status,FORTSPARSE_INVALID_MATRIX,'LU solve needs finite RHS')
                return
            end if
        end do
        if (op==LU_NORMAL) then
            do j=1,self%n
                row=self%p(j)
                self%work(j)=divide_value(rhs(row),scalar_from_real(self%row_norm(row)))
            end do
            do j=1,self%n
                do e=self%lp(j),self%lp(j+1)-1
                    row=self%lr(e)
                    self%work(row)=update_value(self%work(row),self%lv(e),self%work(j))
                end do
            end do
            do j=self%n,1,-1
                self%work(j)=divide_value(self%work(j),self%uv(self%up(j+1)-1))
                do e=self%up(j),self%up(j+1)-2
                    row=self%ur(e)
                    self%work(row)=update_value(self%work(row),self%uv(e),self%work(j))
                end do
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
                do e=self%up(j),self%up(j+1)-2
                    coefficient=adjoint_value(self%uv(e),op)
                    value=update_value(value,coefficient,self%work(self%ur(e)))
                end do
                coefficient=adjoint_value(self%uv(self%up(j+1)-1),op)
                self%work(j)=divide_value(value,coefficient)
            end do
            do j=self%n,1,-1
                value=self%work(j)
                do e=self%lp(j),self%lp(j+1)-1
                    coefficient=adjoint_value(self%lv(e),op)
                    value=update_value(value,coefficient,self%work(self%lr(e)))
                end do
                self%work(j)=value
            end do
            do j=1,self%n
                row=self%p(j)
                rhs(row)=divide_value(self%work(j),scalar_from_real(self%row_norm(row)))
            end do
        end if
        do j=1,self%n
            if (.not.lu_z_finite(rhs(j))) then
                call status_set(status,FORTSPARSE_SINGULAR,'LU nonfinite solution')
                return
            end if
        end do
        call status_set(status,FORTSPARSE_OK,'')
    end subroutine solve_vector

    subroutine solve_block(self,rhs,status,operation)
        class(lu_factor_z_t), intent(inout) :: self
        complex(dp), intent(inout) :: rhs(:,:)
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
        complex(dp) :: value
        value=a
    end function scalar_from_real

    pure function adjoint_value(a,op) result(value)
        complex(dp), intent(in) :: a
        integer, intent(in) :: op
        complex(dp) :: value
        value=a
        if (op==LU_ADJOINT) value=conjg(a)
    end function adjoint_value

    integer(i8) function memory_bytes(self) result(bytes)
        class(lu_factor_z_t), intent(in) :: self
        ! Exact allocated-array payload; object descriptors and allocator
        ! overhead remain process-RSS quantities rather than payload claims.
        bytes=0_i8
        if(allocated(self%lp)) bytes=bytes+ib*size(self%lp,kind=i8)
        if(allocated(self%up)) bytes=bytes+ib*size(self%up,kind=i8)
        if(allocated(self%lr)) bytes=bytes+ib*size(self%lr,kind=i8)
        if(allocated(self%ur)) bytes=bytes+ib*size(self%ur,kind=i8)
        if(allocated(self%p)) bytes=bytes+ib*size(self%p,kind=i8)
        if(allocated(self%pinv)) bytes=bytes+ib*size(self%pinv,kind=i8)
        if(allocated(self%q)) bytes=bytes+ib*size(self%q,kind=i8)
        if(allocated(self%ap)) bytes=bytes+ib*size(self%ap,kind=i8)
        if(allocated(self%ai)) bytes=bytes+ib*size(self%ai,kind=i8)
        if(allocated(self%mark)) bytes=bytes+ib*size(self%mark,kind=i8)
        if(allocated(self%stack)) bytes=bytes+ib*size(self%stack,kind=i8)
        if(allocated(self%edge)) bytes=bytes+ib*size(self%edge,kind=i8)
        if(allocated(self%reach)) bytes=bytes+ib*size(self%reach,kind=i8)
        if(allocated(self%lv)) bytes=bytes+vb*size(self%lv,kind=i8)
        if(allocated(self%uv)) bytes=bytes+vb*size(self%uv,kind=i8)
        if(allocated(self%work)) bytes=bytes+vb*size(self%work,kind=i8)
        if(allocated(self%row_norm)) bytes=bytes+8_i8*size(self%row_norm,kind=i8)
    end function memory_bytes

    subroutine free(self)
        class(lu_factor_z_t), intent(inout) :: self
        if(allocated(self%lp)) deallocate(self%lp)
        if(allocated(self%up)) deallocate(self%up)
        if(allocated(self%lr)) deallocate(self%lr)
        if(allocated(self%ur)) deallocate(self%ur)
        if(allocated(self%p)) deallocate(self%p)
        if(allocated(self%pinv)) deallocate(self%pinv)
        if(allocated(self%q)) deallocate(self%q)
        if(allocated(self%ap)) deallocate(self%ap)
        if(allocated(self%ai)) deallocate(self%ai)
        if(allocated(self%mark)) deallocate(self%mark)
        if(allocated(self%stack)) deallocate(self%stack)
        if(allocated(self%edge)) deallocate(self%edge)
        if(allocated(self%reach)) deallocate(self%reach)
        if(allocated(self%lv)) deallocate(self%lv)
        if(allocated(self%uv)) deallocate(self%uv)
        if(allocated(self%work)) deallocate(self%work)
        if(allocated(self%row_norm)) deallocate(self%row_norm)
        self%n=0
        self%ln=0
        self%un=0
        self%factor_nnz=0_i8
        self%peak_memory_bytes=0_i8
        self%failed_required_bytes=0_i8
        self%failed_available_bytes=0_i8
        self%last_reused_pattern=.false.
        self%ready=.false.
        self%limit=0_i8
    end subroutine free
end module fortsparse_native_lu_z
