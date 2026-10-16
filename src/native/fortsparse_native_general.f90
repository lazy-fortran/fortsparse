module fortsparse_native_general
    ! Native general sparse factors; orchestration only, arithmetic is generated.
    use fortsparse_kinds, only: dp, i8
    use fortsparse_backend, only: sparse_backend_t, sparse_factor_stats_t, &
        FORTSPARSE_SYM_GENERAL, FORTSPARSE_SYM_SPD, FORTSPARSE_SYM_INDEFINITE, &
        FORTSPARSE_SYM_HERMITIAN, FORTSPARSE_ORDER_DEFAULT, FORTSPARSE_MEMORY_LOW
    use fortsparse_csc, only: csc_t, csc_z_t
    use fortsparse_native_hpd, only: native_hpd_t
    use fortsparse_native_ldlt, only: native_ldlt_t
    use fortsparse_native_lu, only: lu_factor_t
    use fortsparse_native_lu_z, only: lu_factor_z_t
    use fortsparse_status, only: fortsparse_status_t, status_set, status_ok, &
        FORTSPARSE_OK, FORTSPARSE_INVALID_MATRIX, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_RESOURCE_LIMIT
    implicit none
    private
    public :: native_general_backend_t
    type, extends(sparse_backend_t) :: native_general_backend_t
        private
        type(native_hpd_t) :: hpd_factor
        type(native_ldlt_t) :: ldlt_factor
        integer :: n = 0
        integer :: mode = 0 ! 0 real LU, 1 complex LU, 2 HPD, 3 real LDLT
        logical :: reused_special = .false.
        type(lu_factor_t) :: real_factor
        type(lu_factor_z_t) :: complex_factor
        logical :: ready = .false., is_complex = .false.
        real(dp), pointer :: user_vector(:) => null()
        type(sparse_factor_stats_t) :: stats
    contains
        procedure :: factor_real => general_factor_real
        procedure :: refactor_real => general_refactor_real
        procedure :: solve_real => general_solve_real
        procedure :: factor_real_raw => general_factor_real_raw
        procedure :: refactor_real_raw => general_refactor_real_raw
        procedure :: solve_real_inplace => general_solve_real_inplace
        procedure :: solve_real_block => general_solve_real_block
        procedure :: solve_real_op => general_solve_real_op
        procedure :: solve_real_block_op => general_solve_real_block_op
        procedure :: factor_complex => general_factor_complex
        procedure :: refactor_complex => general_refactor_complex
        procedure :: solve_complex => general_solve_complex
        procedure :: factor_complex_raw => general_factor_complex_raw
        procedure :: refactor_complex_raw => general_refactor_complex_raw
        procedure :: solve_complex_inplace => general_solve_complex_inplace
        procedure :: solve_complex_block => general_solve_complex_block
        procedure :: solve_complex_op => general_solve_complex_op
        procedure :: solve_complex_block_op => general_solve_complex_block_op
        procedure :: factor_stats => general_stats
        procedure :: vector => general_vector
        procedure :: free => general_free
        final :: general_finalize
    end type native_general_backend_t
contains
    subroutine prepare(self, refine, cap, status)
        class(native_general_backend_t), intent(inout) :: self
        logical, intent(in) :: refine
        integer(i8), intent(out) :: cap
        type(fortsparse_status_t), intent(out) :: status
        self%ready = .false.
        self%stats = sparse_factor_stats_t()
        cap = self%options%host_memory_limit_bytes
        call status_set(status,FORTSPARSE_BACKEND_UNAVAILABLE, &
            'native LU: CPU general CSC, natural ordering, no refinement yet')
        if (self%options%device /= -1 .or. refine) return
        if (self%options%symmetry < FORTSPARSE_SYM_GENERAL .or. &
            self%options%symmetry > FORTSPARSE_SYM_HERMITIAN) return
        if (self%options%ordering /= FORTSPARSE_ORDER_DEFAULT) return
        if (cap > 0_i8 .and. associated(self%user_vector)) then
            cap = cap-8_i8*size(self%user_vector,kind=i8)
            if (cap <= 0_i8) then
                call status_set(status,FORTSPARSE_RESOURCE_LIMIT, &
                    'native LU: retained user vector consumes the host budget')
                return
            end if
        end if
        call status_set(status,FORTSPARSE_OK,'')
    end subroutine prepare

    subroutine general_factor_real(self,A,refine,status)
        class(native_general_backend_t), intent(inout) :: self
        type(csc_t), intent(in) :: A
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        integer(i8) :: cap
        real(dp) :: threshold, started
        call prepare(self,refine,cap,status)
        if (.not. status_ok(status)) return
        call self%complex_factor%free()
        self%reused_special = .false.
        self%is_complex = .false.
        threshold = self%options%pivot_tolerance
        if (threshold < 0.0_dp) threshold = 1.0_dp
        self%n = A%ncol
        started = seconds()
        if (self%options%symmetry == FORTSPARSE_SYM_INDEFINITE .or. &
            self%options%symmetry == FORTSPARSE_SYM_HERMITIAN) then
            call self%real_factor%free()
            call self%hpd_factor%free()
            self%mode = 3
            call self%ldlt_factor%factor(A,status,host_limit_bytes=cap, &
                low_memory=self%options%memory_policy == FORTSPARSE_MEMORY_LOW)
            self%stats%factor_s = seconds()-started
            self%ready = status_ok(status)
            call update_stats(self)
            return
        end if
        self%mode = 0
        call self%hpd_factor%free()
        call self%ldlt_factor%free()
        call self%real_factor%factor(A,status, &
            pivot_threshold=threshold,memory_limit_bytes=cap, &
            low_memory=self%options%memory_policy == FORTSPARSE_MEMORY_LOW)
        self%stats%factor_s = seconds()-started
        self%ready = status_ok(status)
        call update_stats(self)
    end subroutine general_factor_real

    subroutine general_refactor_real(self,A,refine,status,reused)
        class(native_general_backend_t), intent(inout) :: self
        type(csc_t), intent(in) :: A
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(out) :: reused
        call self%factor_real(A,refine,status)
        reused = .false.
        if (status_ok(status)) reused = self%real_factor%last_reused_pattern
    end subroutine general_refactor_real

    subroutine general_factor_real_raw(self,nrow,ncol,nz,col_ptr,row_idx,val,refine,status)
        class(native_general_backend_t), intent(inout) :: self
        integer, intent(in) :: nrow,ncol,nz,col_ptr(:),row_idx(:)
        real(dp), intent(in) :: val(:)
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        integer(i8) :: cap
        real(dp) :: threshold, started
        self%ready = .false.
        if (self%options%symmetry /= FORTSPARSE_SYM_GENERAL) then
            call status_set(status,FORTSPARSE_BACKEND_UNAVAILABLE, &
                'native structured raw CSC factorization not yet supported')
            return
        end if
        call prepare(self,refine,cap,status)
        if (.not. status_ok(status)) return
        if (nrow /= ncol .or. nz /= size(val) .or. nz /= size(row_idx)) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'native LU: raw shape mismatch')
            return
        end if
        call self%complex_factor%free()
        call self%hpd_factor%free()
        call self%ldlt_factor%free()
        self%mode = 0
        self%is_complex = .false.
        threshold = self%options%pivot_tolerance
        if (threshold < 0.0_dp) threshold = 1.0_dp
        started = seconds()
        call self%real_factor%factor_raw(ncol,col_ptr,row_idx,val,status, &
            pivot_threshold=threshold,memory_limit_bytes=cap, &
            low_memory=self%options%memory_policy == FORTSPARSE_MEMORY_LOW)
        self%stats%factor_s = seconds()-started
        self%ready = status_ok(status)
        call update_stats(self)
    end subroutine general_factor_real_raw

    subroutine general_refactor_real_raw(self,nrow,ncol,nz,col_ptr,row_idx,val,refine,status,reused)
        class(native_general_backend_t), intent(inout) :: self
        integer, intent(in) :: nrow,ncol,nz,col_ptr(:),row_idx(:)
        real(dp), intent(in) :: val(:)
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(out) :: reused
        call self%factor_real_raw(nrow,ncol,nz,col_ptr,row_idx,val,refine,status)
        reused = .false.
        if (status_ok(status)) reused = self%real_factor%last_reused_pattern
    end subroutine general_refactor_real_raw

    subroutine general_solve_real(self,b,x,status)
        class(native_general_backend_t), intent(inout) :: self
        real(dp), target, contiguous, intent(in) :: b(:)
        real(dp), target, contiguous, intent(out) :: x(:)
        type(fortsparse_status_t), intent(out) :: status
        call self%solve_real_op(b,x,0,status)
    end subroutine general_solve_real

    subroutine general_solve_real_op(self,b,x,operation,status)
        class(native_general_backend_t), intent(inout) :: self
        real(dp), target, contiguous, intent(in) :: b(:)
        real(dp), target, contiguous, intent(out) :: x(:)
        integer, intent(in) :: operation
        type(fortsparse_status_t), intent(out) :: status
        call solve_shape(self,size(b),size(x),.false.,status)
        if (.not. status_ok(status)) return
        if (self%mode == 3) then
            call self%ldlt_factor%solve(b,x,status,transpose=operation /= 0)
            return
        end if
        x = b
        call self%real_factor%solve(x,status,operation)
    end subroutine general_solve_real_op

    subroutine general_solve_real_inplace(self,b,status)
        class(native_general_backend_t), intent(inout) :: self
        real(dp), intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out) :: status
        call solve_shape(self,size(b),size(b),.false.,status)
        if (.not. status_ok(status)) return
        if (self%mode == 3) then
            call self%ldlt_factor%solve_inplace(b,status)
            return
        end if
        call self%real_factor%solve(b,status)
    end subroutine general_solve_real_inplace

    subroutine general_solve_real_block(self,b,x,status)
        class(native_general_backend_t), intent(inout) :: self
        real(dp), contiguous, intent(in) :: b(:,:)
        real(dp), contiguous, intent(out) :: x(:,:)
        type(fortsparse_status_t), intent(out) :: status
        call self%solve_real_block_op(b,x,0,status)
    end subroutine general_solve_real_block

    subroutine general_solve_real_block_op(self,b,x,operation,status)
        class(native_general_backend_t), intent(inout) :: self
        real(dp), contiguous, intent(in) :: b(:,:)
        real(dp), contiguous, intent(out) :: x(:,:)
        integer, intent(in) :: operation
        type(fortsparse_status_t), intent(out) :: status
        call solve_shape(self,size(b,1),size(x,1),.false.,status)
        if (.not. status_ok(status)) return
        if (size(b,2) /= size(x,2)) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'native LU: RHS count mismatch')
            return
        end if
        if (self%mode == 3) then
            call self%ldlt_factor%solve_block(b,x,status,transpose=operation /= 0)
            return
        end if
        x = b
        call self%real_factor%solve(x,status,operation)
    end subroutine general_solve_real_block_op

    subroutine general_factor_complex(self,A,refine,status)
        class(native_general_backend_t), intent(inout) :: self
        type(csc_z_t), intent(in) :: A
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        integer(i8) :: cap
        real(dp) :: threshold, started
        call prepare(self,refine,cap,status)
        if (.not. status_ok(status)) return
        call self%real_factor%free()
        self%reused_special = .false.
        self%is_complex = .true.
        threshold = self%options%pivot_tolerance
        if (threshold < 0.0_dp) threshold = 1.0_dp
        self%n = A%ncol
        started = seconds()
        if (self%options%symmetry == FORTSPARSE_SYM_SPD) then
            call self%complex_factor%free()
            call self%ldlt_factor%free()
            self%mode = 2
            call self%hpd_factor%refactor(A,status,self%reused_special, &
                memory_limit_bytes=cap, &
                low_memory=self%options%memory_policy == FORTSPARSE_MEMORY_LOW)
            self%stats%factor_s = seconds()-started
            self%ready = status_ok(status)
            call update_stats(self)
            return
        end if
        self%mode = 1
        call self%hpd_factor%free()
        call self%ldlt_factor%free()
        call self%complex_factor%factor(A,status, &
            pivot_threshold=threshold,memory_limit_bytes=cap, &
            low_memory=self%options%memory_policy == FORTSPARSE_MEMORY_LOW)
        self%stats%factor_s = seconds()-started
        self%ready = status_ok(status)
        call update_stats(self)
    end subroutine general_factor_complex

    subroutine general_refactor_complex(self,A,refine,status,reused)
        class(native_general_backend_t), intent(inout) :: self
        type(csc_z_t), intent(in) :: A
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(out) :: reused
        call self%factor_complex(A,refine,status)
        reused = .false.
        if (status_ok(status)) reused = self%complex_factor%last_reused_pattern
        if (self%mode == 2) reused = self%reused_special
    end subroutine general_refactor_complex

    subroutine general_factor_complex_raw(self,nrow,ncol,nz,col_ptr,row_idx,val,refine,status)
        class(native_general_backend_t), intent(inout) :: self
        integer, intent(in) :: nrow,ncol,nz,col_ptr(:),row_idx(:)
        complex(dp), intent(in) :: val(:)
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        integer(i8) :: cap
        real(dp) :: threshold, started
        self%ready = .false.
        if (self%options%symmetry /= FORTSPARSE_SYM_GENERAL) then
            call status_set(status,FORTSPARSE_BACKEND_UNAVAILABLE, &
                'native structured raw CSC factorization not yet supported')
            return
        end if
        call prepare(self,refine,cap,status)
        if (.not. status_ok(status)) return
        if (nrow /= ncol .or. nz /= size(val) .or. nz /= size(row_idx)) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'native LU: raw shape mismatch')
            return
        end if
        call self%real_factor%free()
        call self%hpd_factor%free()
        call self%ldlt_factor%free()
        self%mode = 1
        self%is_complex = .true.
        threshold = self%options%pivot_tolerance
        if (threshold < 0.0_dp) threshold = 1.0_dp
        started = seconds()
        call self%complex_factor%factor_raw(ncol,col_ptr,row_idx,val,status, &
            pivot_threshold=threshold,memory_limit_bytes=cap, &
            low_memory=self%options%memory_policy == FORTSPARSE_MEMORY_LOW)
        self%stats%factor_s = seconds()-started
        self%ready = status_ok(status)
        call update_stats(self)
    end subroutine general_factor_complex_raw

    subroutine general_refactor_complex_raw(self,nrow,ncol,nz,col_ptr,row_idx,val,refine,status,reused)
        class(native_general_backend_t), intent(inout) :: self
        integer, intent(in) :: nrow,ncol,nz,col_ptr(:),row_idx(:)
        complex(dp), intent(in) :: val(:)
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(out) :: reused
        call self%factor_complex_raw(nrow,ncol,nz,col_ptr,row_idx,val,refine,status)
        reused = .false.
        if (status_ok(status)) reused = self%complex_factor%last_reused_pattern
    end subroutine general_refactor_complex_raw

    subroutine general_solve_complex(self,b,x,status)
        class(native_general_backend_t), intent(inout) :: self
        complex(dp), intent(in) :: b(:)
        complex(dp), intent(out) :: x(:)
        type(fortsparse_status_t), intent(out) :: status
        call self%solve_complex_op(b,x,0,status)
    end subroutine general_solve_complex

    subroutine general_solve_complex_op(self,b,x,operation,status)
        class(native_general_backend_t), intent(inout) :: self
        complex(dp), intent(in) :: b(:)
        complex(dp), intent(out) :: x(:)
        integer, intent(in) :: operation
        type(fortsparse_status_t), intent(out) :: status
        call solve_shape(self,size(b),size(x),.true.,status)
        if (.not. status_ok(status)) return
        if (self%mode == 2) then
            call self%hpd_factor%solve(b,x,status,operation)
            return
        end if
        x = b
        call self%complex_factor%solve(x,status,operation)
    end subroutine general_solve_complex_op

    subroutine general_solve_complex_inplace(self,b,status)
        class(native_general_backend_t), intent(inout) :: self
        complex(dp), intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out) :: status
        call solve_shape(self,size(b),size(b),.true.,status)
        if (.not. status_ok(status)) return
        if (self%mode == 2) then
            call self%hpd_factor%solve_inplace(b,status)
            return
        end if
        call self%complex_factor%solve(b,status)
    end subroutine general_solve_complex_inplace

    subroutine general_solve_complex_block(self,b,x,status)
        class(native_general_backend_t), intent(inout) :: self
        complex(dp), contiguous, intent(in) :: b(:,:)
        complex(dp), contiguous, intent(out) :: x(:,:)
        type(fortsparse_status_t), intent(out) :: status
        call self%solve_complex_block_op(b,x,0,status)
    end subroutine general_solve_complex_block

    subroutine general_solve_complex_block_op(self,b,x,operation,status)
        class(native_general_backend_t), intent(inout) :: self
        complex(dp), contiguous, intent(in) :: b(:,:)
        complex(dp), contiguous, intent(out) :: x(:,:)
        integer, intent(in) :: operation
        type(fortsparse_status_t), intent(out) :: status
        call solve_shape(self,size(b,1),size(x,1),.true.,status)
        if (.not. status_ok(status)) return
        if (size(b,2) /= size(x,2)) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'native LU: RHS count mismatch')
            return
        end if
        if (self%mode == 2) then
            call self%hpd_factor%solve(b,x,status,operation)
            return
        end if
        x = b
        call self%complex_factor%solve(x,status,operation)
    end subroutine general_solve_complex_block_op

    subroutine solve_shape(self,nb,nx,is_complex,status)
        class(native_general_backend_t), intent(in) :: self
        integer, intent(in) :: nb,nx
        logical, intent(in) :: is_complex
        type(fortsparse_status_t), intent(out) :: status
        integer :: n
        call status_set(status,FORTSPARSE_NOT_FACTORED,'native LU: no factors for this value type')
        if (.not. self%ready) return
        if (self%is_complex .neqv. is_complex) return
        n = self%real_factor%n
        if (self%is_complex) n = self%complex_factor%n
        if (self%mode == 2) n = self%hpd_factor%dimension()
        if (self%mode == 3) n = self%n
        call status_set(status,FORTSPARSE_INVALID_MATRIX,'native LU: RHS dimension mismatch')
        if (nb /= n .or. nx /= n) return
        call status_set(status,FORTSPARSE_OK,'')
    end subroutine solve_shape

    subroutine update_stats(self)
        class(native_general_backend_t), intent(inout) :: self
        integer(i8) :: vector_bytes
        vector_bytes = 0_i8
        if (associated(self%user_vector)) vector_bytes = 8_i8*size(self%user_vector,kind=i8)
        self%stats%host_current_bytes = self%real_factor%memory_bytes() &
            +self%complex_factor%memory_bytes()+self%hpd_factor%memory_current_bytes() &
            +self%ldlt_factor%retained_bytes()+vector_bytes
        self%stats%host_peak_bytes = max(self%stats%host_peak_bytes, &
            self%real_factor%peak_memory_bytes+self%complex_factor%peak_memory_bytes &
            +self%hpd_factor%memory_peak_bytes()+self%ldlt_factor%peak_bytes()+vector_bytes)
        self%stats%factor_nnz = self%real_factor%factor_nnz+self%complex_factor%factor_nnz &
            +self%hpd_factor%factor_nnz()+self%ldlt_factor%factor_nnz()
        self%stats%budget_required_bytes = self%real_factor%failed_required_bytes &
            +self%complex_factor%failed_required_bytes
        self%stats%budget_available_bytes = self%real_factor%failed_available_bytes &
            +self%complex_factor%failed_available_bytes
        if (self%mode == 2) call self%hpd_factor%memory_failure_bytes( &
            self%stats%budget_required_bytes,self%stats%budget_available_bytes)
        if (self%mode == 3) then
            self%stats%budget_required_bytes = self%ldlt_factor%required_bytes()
            self%stats%budget_available_bytes = self%ldlt_factor%available_bytes()
        end if
        self%stats%factor_mem_mb = real(self%stats%host_peak_bytes,dp)/1.0e6_dp
        self%stats%strict_memory_supported = .true.
        self%stats%active_device = -1
        self%stats%device_bytes = 0_i8
        self%stats%active_threads = 1
        self%stats%nthreads = self%options%nthreads
        self%stats%ordering_used = FORTSPARSE_ORDER_DEFAULT
    end subroutine update_stats

    subroutine general_stats(self,stats,status)
        class(native_general_backend_t), intent(inout) :: self
        type(sparse_factor_stats_t), intent(out) :: stats
        type(fortsparse_status_t), intent(out) :: status
        call update_stats(self)
        stats = self%stats
        call status_set(status,FORTSPARSE_OK,'')
    end subroutine general_stats

    function general_vector(self,n) result(p)
        class(native_general_backend_t), intent(inout) :: self
        integer, intent(in) :: n
        real(dp), pointer :: p(:)
        integer(i8) :: required, cap
        integer :: ierr
        nullify(p)
        if (n < 0) return
        if (associated(self%user_vector)) deallocate(self%user_vector)
        required = self%real_factor%memory_bytes()+self%complex_factor%memory_bytes() &
            +self%hpd_factor%memory_current_bytes()+self%ldlt_factor%retained_bytes() &
            +8_i8*int(n,i8)
        cap = self%options%host_memory_limit_bytes
        if (cap > 0_i8 .and. required > cap) then
            self%stats%budget_required_bytes = 8_i8*int(n,i8)
            self%stats%budget_available_bytes = max(0_i8, cap-required+8_i8*int(n,i8))
            call update_stats(self)
            return
        end if
        allocate(self%user_vector(n),stat=ierr)
        if (ierr /= 0) return
        p => self%user_vector
        call update_stats(self)
    end function general_vector

    subroutine general_free(self)
        class(native_general_backend_t), intent(inout) :: self
        call self%real_factor%free()
        call self%complex_factor%free()
        call self%hpd_factor%free()
        call self%ldlt_factor%free()
        self%n = 0
        if (associated(self%user_vector)) deallocate(self%user_vector)
        nullify(self%user_vector)
        self%ready = .false.
        self%stats = sparse_factor_stats_t()
    end subroutine general_free

    subroutine general_finalize(self)
        type(native_general_backend_t), intent(inout) :: self
        call self%free()
    end subroutine general_finalize

    real(dp) function seconds()
        integer(i8) :: clock,rate
        call system_clock(clock,rate)
        seconds = real(clock,dp)/real(rate,dp)
    end function seconds
end module fortsparse_native_general
