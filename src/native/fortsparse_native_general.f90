module fortsparse_native_general
    ! Native general sparse factors; orchestration only, arithmetic is generated.
    use fortsparse_kinds, only: dp, i8
    use fortsparse_memory, only: refinement_payload_bytes, memory_checked_add, memory_checked_mul
    use, intrinsic :: iso_c_binding, only: c_int
    use fortsparse_backend, only: sparse_backend_t, sparse_factor_stats_t, &
        backend_factor_real_raw, backend_refactor_real_raw, &
        FORTSPARSE_SYM_GENERAL, FORTSPARSE_SYM_SPD, FORTSPARSE_SYM_INDEFINITE, &
        FORTSPARSE_SYM_HERMITIAN, FORTSPARSE_ORDER_DEFAULT, FORTSPARSE_ORDER_AMD, &
        FORTSPARSE_ORDER_ND, &
        FORTSPARSE_MEMORY_LOW
    use fortsparse_csc, only: csc_t, csc_z_t
    use fortsparse_spd_leaves, only: spd_update
    use fortsparse_lu_leaves, only: lu_z_update
    use fortsparse_native_hpd, only: native_hpd_t
    use fortsparse_native_ldlt, only: native_ldlt_t, ldlt_check_matrix
    use fortsparse_native_ldlt_z, only: native_ldlt_z_t
    use fortsparse_native_lu, only: lu_factor_t
    use fortsparse_native_lu_z, only: lu_factor_z_t
    use fortsparse_status, only: fortsparse_status_t, status_set, status_ok, &
        FORTSPARSE_OK, FORTSPARSE_INVALID_MATRIX, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_RESOURCE_LIMIT, &
        FORTSPARSE_OUT_OF_MEMORY, FORTSPARSE_SINGULAR, FORTSPARSE_INDEX_OVERFLOW
    use fortsparse_native_amd, only: bounded_amd_available, bounded_amd_size, bounded_amd_order
    use fortsparse_native_matching, only: structural_row_matching, matching_graph_rows
    use fortsparse_native_residual, only: native_original_norms, native_residual_error
    implicit none
    private
    public :: native_general_backend_t
    integer(i8), parameter :: matching_integer_bytes=int(storage_size(0)/8,i8)
    interface
        integer(c_int) function order_amd_general(n,nnz,ptr,row,perm,inverse) &
            bind(C,name='fsparse_native_order_amd_general')
            import :: c_int
            integer(c_int), value :: n,nnz
            integer(c_int), intent(in) :: ptr(*),row(*)
            integer(c_int), intent(out) :: perm(*),inverse(*)
        end function order_amd_general
        integer(c_int) function order_nd_general(n,nnz,ptr,row,perm,inverse) &
            bind(C,name='fsparse_native_order_general')
            import :: c_int
            integer(c_int), value :: n,nnz
            integer(c_int), intent(in) :: ptr(*),row(*)
            integer(c_int), intent(out) :: perm(*),inverse(*)
        end function order_nd_general
    end interface
    type, extends(sparse_backend_t) :: native_general_backend_t
        private
        logical, public :: reference_only = .false.
        integer :: selected_order = FORTSPARSE_ORDER_DEFAULT
        integer :: retained_order = -1
        logical :: same_ordering_strategy = .false.
        logical :: bounded_amd = .false., factor_invoked = .false.
        integer(i8) :: refinement_quota = 0_i8
        logical :: previous_matching = .false.
        integer, allocatable :: matching_preference(:), matching_inverse(:)
        integer, allocatable :: matching_work(:), matching_rows(:)
        type(native_hpd_t) :: hpd_factor
        type(native_ldlt_t) :: ldlt_factor
        type(native_ldlt_z_t) :: complex_ldlt_factor
        integer :: n = 0
        integer :: mode = 0 ! 0 real LU, 1 complex LU, 2 HPD, 3 real LDLT, 4 LDLH, 5 complex LDLT
        logical :: reused_special = .false.
        type(lu_factor_t) :: real_factor
        type(lu_factor_z_t) :: complex_factor
        logical :: ready = .false., is_complex = .false., converted_raw = .false.
        real(dp), pointer :: user_vector(:) => null()
        type(sparse_factor_stats_t) :: stats
        logical :: refining = .false.
        integer, allocatable :: original_ptr(:), original_rows(:)
        real(dp), allocatable :: original_values(:), saved_rhs(:), residual(:)
        complex(dp), allocatable :: original_zvalues(:), saved_zrhs(:), zresidual(:)
        real(dp) :: matrix_norm = 0.0_dp, transpose_norm = 0.0_dp
        integer(i8) :: factor_external_live = 0_i8
        integer(i8) :: failed_required = 0_i8, failed_available = 0_i8
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
    subroutine column_order(self,n,ptr,rows,is_complex,order,cap,status)
        class(native_general_backend_t), intent(inout) :: self
        integer, intent(in) :: n,ptr(:),rows(:)
        logical, intent(in) :: is_complex
        integer, allocatable, intent(out) :: order(:)
        integer(i8), intent(inout) :: cap
        type(fortsparse_status_t), intent(out) :: status
        integer, allocatable :: inverse(:)
        integer(c_int), allocatable :: arena(:)
        integer :: ierr,code,low
        integer(i8) :: order_bytes,arena_bytes,output_bytes,total
        logical :: matched
        real(dp) :: started
        call status_set(status,FORTSPARSE_OK,'')
        if (self%options%ordering == FORTSPARSE_ORDER_DEFAULT) return
        if (n <= 0 .or. n == huge(n)) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'native LU: ordering needs positive dimension')
            return
        end if
        if (size(ptr) /= n+1 .or. size(rows) == huge(n)) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'native LU: ordering CSC extent mismatch')
            return
        end if
        if (self%bounded_amd .and. (n>=huge(0_c_int)/4 .or. size(rows)>=huge(0_c_int)/4)) then
            call status_set(status,FORTSPARSE_INDEX_OVERFLOW,'native LU: bounded AMD exceeds audited index range')
            return
        end if
        started = seconds()
        order_bytes=0_i8
        call memory_checked_mul(4_i8,int(n,i8),order_bytes,status)
        if (.not.status_ok(status)) return
        if (self%bounded_amd) then
            call admit_order(self,order_bytes,status)
            if (.not.status_ok(status)) goto 900
        end if
        allocate(order(n),stat=ierr)
        if (ierr /= 0) then
            call status_set(status,FORTSPARSE_OUT_OF_MEMORY,'native LU: ordering output allocation failed')
            goto 900
        end if
        if (self%bounded_amd) call record_order_peak(self,order_bytes)
        matched = .false.
        if (self%same_ordering_strategy) then
            if (is_complex) then
                call self%complex_factor%copy_column_order(n,ptr,rows,order,matched)
            else if (self%mode == 3) then
                call self%ldlt_factor%copy_column_order(n,ptr,rows,order,matched)
            else
                call self%real_factor%copy_column_order(n,ptr,rows,order,matched)
            end if
        end if
        if (.not.matched) then
            if (self%bounded_amd) then
                low=0
                if (self%options%memory_policy==FORTSPARSE_MEMORY_LOW) low=1
                code=bounded_amd_size(n,size(rows),ptr,rows,low,arena_bytes,output_bytes)
                call status_set(status,code,'native LU: bounded AMD graph/ABI is unavailable or invalid')
                if (.not.status_ok(status)) goto 900
                total=0_i8
                call memory_checked_add(arena_bytes,output_bytes,total,status)
                if (.not.status_ok(status)) goto 900
                ! Includes both outputs and the complete construction arena.
                call admit_order(self,total,status)
                if (.not.status_ok(status)) goto 900
            end if
            allocate(inverse(n),stat=ierr)
            if (ierr /= 0) then
                call status_set(status,FORTSPARSE_OUT_OF_MEMORY,'native LU: ordering inverse allocation failed')
                goto 900
            end if
            if (self%bounded_amd) then
                call record_order_peak(self,output_bytes)
                ! The same admitted total is checked before this allocation too.
                call admit_order(self,total,status)
                if (.not.status_ok(status)) goto 900
                allocate(arena(arena_bytes/4_i8),stat=ierr)
                if (ierr/=0) then
                    call status_set(status,FORTSPARSE_OUT_OF_MEMORY,'native LU: bounded AMD arena allocation failed')
                    goto 900
                end if
                call record_order_peak(self,total)
                code=bounded_amd_order(n,size(rows),ptr,rows,low,arena,arena_bytes,order,inverse)
            else if (self%options%ordering == FORTSPARSE_ORDER_ND) then
                code=order_nd_general(n,size(rows),ptr,rows,order,inverse)
            else
                code=order_amd_general(n,size(rows),ptr,rows,order,inverse)
            end if
            call status_set(status,code,'native LU: column ordering failed or unavailable')
            if (.not.status_ok(status)) goto 900
        end if
        if (self%bounded_amd) then
            ! Keep the permutation input charged during LU growth. Exhaustion
            ! must never become the factor's zero-means-unlimited cap.
            if (cap>0_i8) then
                if (cap<=order_bytes) then
                    self%failed_required=order_bytes+1_i8
                    self%failed_available=cap
                    call status_set(status,FORTSPARSE_RESOURCE_LIMIT,'native LU: AMD output leaves no factor budget')
                    goto 900
                end if
                cap=cap-order_bytes
            end if
            self%factor_external_live=self%factor_external_live+order_bytes
        end if
        self%selected_order=self%options%ordering
        call status_set(status,FORTSPARSE_OK,'')
900     continue
        if (allocated(arena)) deallocate(arena)
        if (allocated(inverse)) deallocate(inverse)
        if (.not.status_ok(status)) then
            if (allocated(order)) deallocate(order)
            call update_stats(self)
        end if
        self%stats%analysis_s=seconds()-started
    end subroutine column_order

    subroutine general_order(self,n,ptr,rows,is_complex,order,cap,status,values,zvalues)
        class(native_general_backend_t), intent(inout) :: self
        integer, intent(in) :: n,ptr(:),rows(:)
        logical, intent(in) :: is_complex
        integer, allocatable, intent(out) :: order(:)
        integer(i8), intent(inout) :: cap
        type(fortsparse_status_t), intent(out) :: status
        real(dp), optional, intent(in) :: values(:)
        complex(dp), optional, intent(in) :: zvalues(:)
        integer(i8) :: entries,bytes,preference_bytes
        integer :: ierr
        real(dp) :: started
        if(.not.self%options%structural_matching) then
            call column_order(self,n,ptr,rows,is_complex,order,cap,status)
            return
        end if
        started=seconds()
        call status_set(status,FORTSPARSE_INVALID_MATRIX,'matching: dimension exceeds workspace range')
        if(n<0.or.n>huge(0)/7) return
        entries=7_i8*int(n,i8)
        if(self%options%ordering/=FORTSPARSE_ORDER_DEFAULT) entries=entries+size(rows,kind=i8)
        bytes=0_i8
        call memory_checked_mul(matching_integer_bytes,entries,bytes,status)
        if(.not.status_ok(status)) return
        call admit_order(self,bytes,status)
        if(.not.status_ok(status)) goto 900
        allocate(self%matching_preference(n),self%matching_inverse(n), &
            self%matching_work(5*n),stat=ierr)
        if(ierr/=0) then
            call status_set(status,FORTSPARSE_OUT_OF_MEMORY,'matching workspace allocation failed')
            goto 900
        end if
        call record_order_peak(self,0_i8)
        if(is_complex) then
            call structural_row_matching(n,ptr,rows,zvalues,self%matching_preference, &
                self%matching_inverse,self%matching_work,status)
        else
            call structural_row_matching(n,ptr,rows,values,self%matching_preference, &
                self%matching_inverse,self%matching_work,status)
        end if
        if(.not.status_ok(status)) goto 900
        if(self%options%ordering/=FORTSPARSE_ORDER_DEFAULT) then
            allocate(self%matching_rows(size(rows)),stat=ierr)
            if(ierr/=0) then
                call status_set(status,FORTSPARSE_OUT_OF_MEMORY,'matching graph allocation failed')
                goto 900
            end if
            call record_order_peak(self,0_i8)
            call matching_graph_rows(ptr,rows,self%matching_inverse,self%matching_rows)
            call column_order(self,n,ptr,self%matching_rows,is_complex,order,cap,status)
        else
            call column_order(self,n,ptr,rows,is_complex,order,cap,status)
        end if
900     continue
        call record_order_peak(self,0_i8)
        if(allocated(self%matching_inverse)) deallocate(self%matching_inverse)
        if(allocated(self%matching_work)) deallocate(self%matching_work)
        if(allocated(self%matching_rows)) deallocate(self%matching_rows)
        if(status_ok(status)) then
            preference_bytes=matching_integer_bytes*int(n,i8)
            if(cap>0_i8) then
                if(cap<=preference_bytes) then
                    self%failed_required=preference_bytes+1_i8
                    self%failed_available=cap
                    call status_set(status,FORTSPARSE_RESOURCE_LIMIT,'matching preference leaves no factor budget')
                else
                    cap=cap-preference_bytes
                end if
            end if
            if(status_ok(status)) self%factor_external_live=self%factor_external_live+preference_bytes
        end if
        if(.not.status_ok(status)) then
            call release_matching(self)
            call update_stats(self)
        end if
        self%stats%analysis_s=seconds()-started
    end subroutine general_order

    integer(i8) function matching_bytes(self) result(bytes)
        class(native_general_backend_t), intent(in) :: self
        bytes=0_i8
        if(allocated(self%matching_preference)) bytes=bytes+matching_integer_bytes*size(self%matching_preference,kind=i8)
        if(allocated(self%matching_inverse)) bytes=bytes+matching_integer_bytes*size(self%matching_inverse,kind=i8)
        if(allocated(self%matching_work)) bytes=bytes+matching_integer_bytes*size(self%matching_work,kind=i8)
        if(allocated(self%matching_rows)) bytes=bytes+matching_integer_bytes*size(self%matching_rows,kind=i8)
    end function matching_bytes

    subroutine release_matching(self)
        class(native_general_backend_t), intent(inout) :: self
        if(allocated(self%matching_preference)) deallocate(self%matching_preference)
        if(allocated(self%matching_inverse)) deallocate(self%matching_inverse)
        if(allocated(self%matching_work)) deallocate(self%matching_work)
        if(allocated(self%matching_rows)) deallocate(self%matching_rows)
    end subroutine release_matching

    subroutine admit_order(self,bytes,status)
        class(native_general_backend_t), intent(inout) :: self
        integer(i8), intent(in) :: bytes
        type(fortsparse_status_t), intent(out) :: status
        integer(i8) :: base,total,reserved
        reserved=max(0_i8,self%refinement_quota-original_bytes(self))
        base=0_i8
        total=0_i8
        call memory_checked_add(current_bytes(self),reserved,base,status)
        if (.not.status_ok(status)) return
        call memory_checked_add(base,bytes,total,status)
        if (.not.status_ok(status)) return
        if (self%options%host_memory_limit_bytes>0_i8) then
            if (total>self%options%host_memory_limit_bytes) then
                self%failed_required=bytes
                self%failed_available=max(0_i8,self%options%host_memory_limit_bytes-base)
                call status_set(status,FORTSPARSE_RESOURCE_LIMIT,'native LU: ordering/matching workspace exceeds host cap')
            end if
        end if
    end subroutine admit_order

    subroutine record_order_peak(self,bytes)
        class(native_general_backend_t), intent(inout) :: self
        integer(i8), intent(in) :: bytes
        self%stats%host_peak_bytes=max(self%stats%host_peak_bytes,current_bytes(self)+bytes)
    end subroutine record_order_peak

    subroutine prepare(self, refine, cap, status, real_structured)
        class(native_general_backend_t), intent(inout) :: self
        logical, intent(in) :: refine
        logical, intent(in), optional :: real_structured
        logical :: supported_structure
        integer(i8), intent(out) :: cap
        type(fortsparse_status_t), intent(out) :: status
        self%ready = .false.
        call release_matching(self)
        self%converted_raw = .false.
        self%same_ordering_strategy = self%selected_order == self%options%ordering
        self%selected_order = FORTSPARSE_ORDER_DEFAULT
        self%stats = sparse_factor_stats_t()
        self%bounded_amd=.false.
        self%factor_invoked=.false.
        self%refinement_quota=0_i8
        self%stats%host_peak_bytes = current_bytes(self)
        self%failed_required = 0_i8
        self%failed_available = 0_i8
        self%refining = .false.
        cap = self%options%host_memory_limit_bytes
        if (cap < 0_i8) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'native CPU: negative host budget')
            return
        end if
        call status_set(status,FORTSPARSE_BACKEND_UNAVAILABLE, &
            'native CPU: requested factorization or ordering is unavailable')
        if (self%options%device /= -1) return
        if (self%options%symmetry < FORTSPARSE_SYM_GENERAL .or. &
            self%options%symmetry > FORTSPARSE_SYM_HERMITIAN) return
        if (self%options%ordering /= FORTSPARSE_ORDER_DEFAULT) then
            if (self%options%ordering /= FORTSPARSE_ORDER_AMD .and. &
                self%options%ordering /= FORTSPARSE_ORDER_ND) return
            if (self%options%ordering == FORTSPARSE_ORDER_ND .and. &
                self%options%symmetry /= FORTSPARSE_SYM_GENERAL) return
            supported_structure = self%options%symmetry == FORTSPARSE_SYM_GENERAL
            if (present(real_structured)) then
                if (real_structured) supported_structure = supported_structure .or. &
                    self%options%symmetry == FORTSPARSE_SYM_INDEFINITE .or. &
                    self%options%symmetry == FORTSPARSE_SYM_HERMITIAN
            end if
            if (.not. supported_structure) return
            if (self%reference_only) return
            if (cap > 0_i8 .or. self%options%memory_policy == FORTSPARSE_MEMORY_LOW) then
                if (self%options%ordering == FORTSPARSE_ORDER_ND) then
                    call status_set(status,FORTSPARSE_BACKEND_UNAVAILABLE, &
                        'native CPU: opaque METIS workspace cannot enforce a strict cap or LOW')
                    return
                end if
                if (bounded_amd_available()/=1) then
                    call status_set(status,FORTSPARSE_BACKEND_UNAVAILABLE, &
                        'native LU: audited AMD caller-workspace ABI unavailable')
                    return
                end if
                self%bounded_amd=.true.
                self%same_ordering_strategy=self%retained_order==self%options%ordering
            end if
        end if
        if(self%options%structural_matching.or.self%previous_matching) self%same_ordering_strategy=.false.
        self%previous_matching=self%options%structural_matching
        if(self%options%structural_matching.and.self%options%symmetry/=FORTSPARSE_SYM_GENERAL) return
        if (cap > 0_i8 .and. associated(self%user_vector)) then
            cap = cap-8_i8*size(self%user_vector,kind=i8)
            if (cap <= 0_i8) then
                self%failed_required=1_i8
                self%failed_available=0_i8
                call update_stats(self)
                call status_set(status,FORTSPARSE_RESOURCE_LIMIT, &
                    'native LU: retained user vector consumes the host budget')
                return
            end if
        end if
        if (.not.refine) call release_original(self)
        call status_set(status,FORTSPARSE_OK,'')
    end subroutine prepare

    subroutine general_factor_real(self,A,refine,status)
        class(native_general_backend_t), intent(inout) :: self
        type(csc_t), intent(in) :: A
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        integer(i8) :: cap
        integer, allocatable :: order(:)
        real(dp) :: threshold, started
        call prepare(self,refine,cap,status,real_structured=.true.)
        if (.not. status_ok(status)) return
        if (.not.allocated(A%col_ptr) .or. .not.allocated(A%row_idx) .or. &
            .not.allocated(A%val)) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX, &
                'native CPU: CSC storage is unallocated')
            return
        end if
        call self%complex_factor%free()
        call self%complex_ldlt_factor%free()
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
            ! Validate numerical symmetry before constructing an ordering graph.
            call ldlt_check_matrix(A,status)
            if (.not.status_ok(status)) return
            call reserve_original(self,A%ncol,A%nnz,refine,.false.,cap,status)
            if (.not.status_ok(status)) return
            call column_order(self,A%ncol,A%col_ptr,A%row_idx,.false.,order,cap,status)
            if (.not.status_ok(status)) return
            self%factor_invoked=.true.
            self%retained_order=-1
            call self%ldlt_factor%factor(A,status,host_limit_bytes=cap, &
                low_memory=self%options%memory_policy == FORTSPARSE_MEMORY_LOW, &
                column_permutation=order)
            if (allocated(order)) deallocate(order)
            self%stats%factor_s = seconds()-started
            self%ready = status_ok(status)
            if (self%ready) self%retained_order=self%selected_order
            if (self%ready) call retain_real(self,A%ncol,A%col_ptr,A%row_idx,A%val,refine,status)
            self%ready = status_ok(status)
            call update_stats(self)
            return
        end if
        self%mode = 0
        call self%hpd_factor%free()
        call self%ldlt_factor%free()
        call reserve_original(self,A%ncol,A%nnz,refine,.false.,cap,status)
        if (.not.status_ok(status)) return
        call general_order(self,A%ncol,A%col_ptr,A%row_idx,.false., &
            order,cap,status, values=A%val)
        if (.not.status_ok(status)) return
        started = seconds()
        self%factor_invoked=.true.
        self%retained_order=-1
        call self%real_factor%factor(A,status, &
            column_permutation=order,pivot_threshold=threshold,memory_limit_bytes=cap, &
            pruning=.not.self%reference_only, row_preference=self%matching_preference, &
            low_memory=self%options%memory_policy == FORTSPARSE_MEMORY_LOW)
        self%stats%factor_s = seconds()-started
        if (allocated(order)) deallocate(order)
        call release_matching(self)
        self%ready = status_ok(status)
        if (self%ready) self%retained_order=self%selected_order
        if (self%ready) call retain_real(self,A%ncol,A%col_ptr,A%row_idx,A%val,refine,status)
        self%ready = status_ok(status)
        call update_stats(self)
    end subroutine general_factor_real

    subroutine general_refactor_real(self,A,refine,status,reused)
        class(native_general_backend_t), intent(inout) :: self
        type(csc_t), intent(in) :: A
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(out) :: reused
        integer(i8) :: cap
        integer, allocatable :: order(:)
        real(dp) :: started
        integer(i8) :: attempt_peak
        attempt_peak=0_i8
        if (self%mode == 3 .and. &
            (self%options%symmetry == FORTSPARSE_SYM_INDEFINITE .or. &
             self%options%symmetry == FORTSPARSE_SYM_HERMITIAN)) then
            reused = .false.
            call prepare(self,refine,cap,status,real_structured=.true.)
            if (.not. status_ok(status)) return
            if (.not.self%same_ordering_strategy) then
                call self%factor_real(A,refine,status)
                return
            end if
            call ldlt_check_matrix(A,status)
            if (.not.status_ok(status)) return
            started = seconds()
            call reserve_original(self,A%ncol,A%nnz,refine,.false.,cap,status)
            if (.not.status_ok(status)) return
            call column_order(self,A%ncol,A%col_ptr,A%row_idx,.false.,order,cap,status)
            if (.not.status_ok(status)) return
            self%factor_invoked=.true.
            self%retained_order=-1
            call self%ldlt_factor%refactor(A,status,host_limit_bytes=cap, &
                low_memory=self%options%memory_policy == FORTSPARSE_MEMORY_LOW, &
                column_permutation=order)
            if (allocated(order)) deallocate(order)
            if (status%code /= FORTSPARSE_INVALID_MATRIX) then
                self%stats%factor_s = seconds()-started
                self%ready = status_ok(status)
                if (self%ready) self%retained_order=self%selected_order
                if (self%ready) call retain_real(self,A%ncol,A%col_ptr,A%row_idx,A%val,refine,status)
                self%ready = status_ok(status)
                reused = self%ready
                call update_stats(self)
                return
            end if
            attempt_peak=self%stats%host_peak_bytes
        end if
        call self%factor_real(A,refine,status)
        if (self%selected_order /= FORTSPARSE_ORDER_AMD .or. self%bounded_amd) then
            self%stats%host_peak_bytes=max(self%stats%host_peak_bytes,attempt_peak)
        end if
        reused = .false.
        if (status_ok(status)) reused = self%real_factor%last_reused_pattern &
            .and. self%same_ordering_strategy
    end subroutine general_refactor_real

    subroutine structured_raw_options(self, status)
        class(native_general_backend_t), intent(in) :: self
        type(fortsparse_status_t), intent(out) :: status
        call status_set(status,FORTSPARSE_INVALID_MATRIX,'invalid structured raw memory options')
        if (self%options%host_memory_limit_bytes < 0_i8 .or. &
            self%options%device_memory_limit_bytes < 0_i8) return
        call status_set(status,FORTSPARSE_BACKEND_UNAVAILABLE, &
            'strict memory: structured raw CSC conversion is not budgeted')
        if (self%options%host_memory_limit_bytes > 0_i8 .or. &
            self%options%device_memory_limit_bytes > 0_i8 .or. &
            self%options%memory_policy == FORTSPARSE_MEMORY_LOW .or. &
            self%options%device /= -1) return
        call status_set(status,FORTSPARSE_OK,'')
    end subroutine structured_raw_options

    subroutine general_factor_real_raw(self,nrow,ncol,nz,col_ptr,row_idx,val,refine,status)
        class(native_general_backend_t), intent(inout) :: self
        integer, intent(in) :: nrow,ncol,nz,col_ptr(:),row_idx(:)
        real(dp), intent(in) :: val(:)
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        integer(i8) :: cap
        integer, allocatable :: order(:)
        real(dp) :: threshold, started
        self%ready = .false.
        if (self%options%symmetry == FORTSPARSE_SYM_INDEFINITE .or. &
            self%options%symmetry == FORTSPARSE_SYM_HERMITIAN) then
            call structured_raw_options(self, status)
            if (.not.status_ok(status)) return
            call backend_factor_real_raw(self,nrow,ncol,nz,col_ptr,row_idx,val,refine,status)
            self%converted_raw = .true.
            call update_stats(self)
            return
        end if
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
        call self%complex_ldlt_factor%free()
        call self%hpd_factor%free()
        call self%ldlt_factor%free()
        self%mode = 0
        self%is_complex = .false.
        threshold = self%options%pivot_tolerance
        if (threshold < 0.0_dp) threshold = 1.0_dp
        self%n = ncol
        started = seconds()
        call reserve_original(self,ncol,nz,refine,.false.,cap,status)
        if (.not.status_ok(status)) return
        call general_order(self,ncol,col_ptr,row_idx,.false.,order,cap,status, values=val)
        if (.not.status_ok(status)) return
        started = seconds()
        self%factor_invoked=.true.
        self%retained_order=-1
        call self%real_factor%factor_raw(ncol,col_ptr,row_idx,val,status, &
            column_permutation=order,pivot_threshold=threshold,memory_limit_bytes=cap, &
            pruning=.not.self%reference_only, row_preference=self%matching_preference, &
            low_memory=self%options%memory_policy == FORTSPARSE_MEMORY_LOW)
        self%stats%factor_s = seconds()-started
        if (allocated(order)) deallocate(order)
        call release_matching(self)
        self%ready = status_ok(status)
        if (self%ready) self%retained_order=self%selected_order
        if (self%ready) call retain_real(self,ncol,col_ptr,row_idx,val,refine,status)
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
        reused = .false.
        self%ready = .false.
        if (self%options%symmetry == FORTSPARSE_SYM_INDEFINITE .or. &
            self%options%symmetry == FORTSPARSE_SYM_HERMITIAN) then
            call structured_raw_options(self, status)
            if (.not.status_ok(status)) return
            call backend_refactor_real_raw(self,nrow,ncol,nz,col_ptr,row_idx,val, &
                refine,status,reused)
            self%converted_raw = .true.
            call update_stats(self)
            return
        end if
        call self%factor_real_raw(nrow,ncol,nz,col_ptr,row_idx,val,refine,status)
        reused = .false.
        if (status_ok(status)) reused = self%real_factor%last_reused_pattern &
            .and. self%same_ordering_strategy
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
        self%stats%refinement_steps = 0
        self%stats%backward_error = -1.0_dp
        call solve_shape(self,size(b),size(x),.false.,status,operation)
        if (.not.status_ok(status)) return
        if (self%refining) self%saved_rhs = b
        x = b
        call solve_real_values(self%real_factor,self%ldlt_factor, &
            self%mode,x,operation,status)
        if (status_ok(status) .and. self%refining) call refine_real(self,x,operation,status)
    end subroutine general_solve_real_op

    subroutine general_solve_real_inplace(self,b,status)
        class(native_general_backend_t), intent(inout) :: self
        real(dp), intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out) :: status
        self%stats%refinement_steps = 0
        self%stats%backward_error = -1.0_dp
        call solve_shape(self,size(b),size(b),.false.,status)
        if (.not.status_ok(status)) return
        if (self%refining) self%saved_rhs = b
        call solve_real_values(self%real_factor,self%ldlt_factor, &
            self%mode,b,0,status)
        if (status_ok(status) .and. self%refining) call refine_real(self,b,0,status)
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
        integer :: j,steps
        real(dp) :: error
        self%stats%refinement_steps = 0
        self%stats%backward_error = -1.0_dp
        call solve_shape(self,size(b,1),size(x,1),.false.,status,operation)
        if (.not.status_ok(status)) return
        if (size(b,2) /= size(x,2)) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'native CPU: RHS count mismatch')
            return
        end if
        steps = 0
        error = -1.0_dp
        do j = 1,size(b,2)
            call self%solve_real_op(b(:,j),x(:,j),operation,status)
            steps = max(steps,self%stats%refinement_steps)
            error = max(error,self%stats%backward_error)
            if (.not.status_ok(status)) return
        end do
        self%stats%refinement_steps = steps
        self%stats%backward_error = error
    end subroutine general_solve_real_block_op

    subroutine general_factor_complex(self,A,refine,status)
        class(native_general_backend_t), intent(inout) :: self
        type(csc_z_t), intent(in) :: A
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        integer(i8) :: cap
        integer, allocatable :: order(:)
        real(dp) :: threshold, started
        call prepare(self,refine,cap,status)
        if (.not. status_ok(status)) return
        if (.not.allocated(A%col_ptr) .or. .not.allocated(A%row_idx) .or. &
            .not.allocated(A%val)) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX, &
                'native CPU: CSC storage is unallocated')
            return
        end if
        call self%real_factor%free()
        self%reused_special = .false.
        self%is_complex = .true.
        threshold = self%options%pivot_tolerance
        if (threshold < 0.0_dp) threshold = 1.0_dp
        self%n = A%ncol
        started = seconds()
        if (self%options%symmetry == FORTSPARSE_SYM_HERMITIAN .or. &
            self%options%symmetry == FORTSPARSE_SYM_INDEFINITE) then
            call self%complex_factor%free()
            call self%hpd_factor%free()
            call self%ldlt_factor%free()
            self%mode = 4
            if (self%options%symmetry == FORTSPARSE_SYM_INDEFINITE) self%mode = 5
            call reserve_original(self,A%ncol,A%nnz,refine,.true.,cap,status)
            if (.not.status_ok(status)) return
            call self%complex_ldlt_factor%refactor(A,status,hermitian=self%mode==4, &
                host_limit_bytes=cap,low_memory=self%options%memory_policy==FORTSPARSE_MEMORY_LOW)
            self%reused_special=status_ok(status)
            if (status%code==FORTSPARSE_INVALID_MATRIX) then
                call self%complex_ldlt_factor%factor(A,status,hermitian=self%mode==4, &
                    host_limit_bytes=cap,low_memory=self%options%memory_policy==FORTSPARSE_MEMORY_LOW)
            end if
            self%stats%factor_s=seconds()-started
            self%ready=status_ok(status)
            if (self%ready) call retain_complex(self,A%ncol,A%col_ptr,A%row_idx,A%val,refine,status)
            self%ready=status_ok(status)
            call update_stats(self)
            return
        end if
        call self%complex_ldlt_factor%free()
        if (self%options%symmetry == FORTSPARSE_SYM_SPD) then
            call self%complex_factor%free()
            call self%ldlt_factor%free()
            self%mode = 2
            call reserve_original(self,A%ncol,A%nnz,refine,.true.,cap,status)
            if (.not.status_ok(status)) return
            call self%hpd_factor%refactor(A,status,self%reused_special, &
                memory_limit_bytes=cap, &
                low_memory=self%options%memory_policy == FORTSPARSE_MEMORY_LOW)
            self%stats%factor_s = seconds()-started
            self%ready = status_ok(status)
            if (self%ready) call retain_complex(self,A%ncol,A%col_ptr,A%row_idx,A%val,refine,status)
            self%ready = status_ok(status)
            call update_stats(self)
            return
        end if
        self%mode = 1
        call self%hpd_factor%free()
        call self%ldlt_factor%free()
        call reserve_original(self,A%ncol,A%nnz,refine,.true.,cap,status)
        if (.not.status_ok(status)) return
        call general_order(self,A%ncol,A%col_ptr,A%row_idx,.true., &
            order,cap,status, zvalues=A%val)
        if (.not.status_ok(status)) return
        started = seconds()
        self%factor_invoked=.true.
        self%retained_order=-1
        call self%complex_factor%factor(A,status, &
            column_permutation=order,pivot_threshold=threshold,memory_limit_bytes=cap, &
            pruning=.not.self%reference_only, row_preference=self%matching_preference, &
            low_memory=self%options%memory_policy == FORTSPARSE_MEMORY_LOW)
        self%stats%factor_s = seconds()-started
        if (allocated(order)) deallocate(order)
        call release_matching(self)
        self%ready = status_ok(status)
        if (self%ready) self%retained_order=self%selected_order
        if (self%ready) call retain_complex(self,A%ncol,A%col_ptr,A%row_idx,A%val,refine,status)
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
        if (.not.status_ok(status)) return
        reused = self%complex_factor%last_reused_pattern .and. self%same_ordering_strategy
        if (self%mode == 2 .or. self%mode == 4 .or. self%mode == 5) reused = self%reused_special
    end subroutine general_refactor_complex

    subroutine general_factor_complex_raw(self,nrow,ncol,nz,col_ptr,row_idx,val,refine,status)
        class(native_general_backend_t), intent(inout) :: self
        integer, intent(in) :: nrow,ncol,nz,col_ptr(:),row_idx(:)
        complex(dp), intent(in) :: val(:)
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        integer(i8) :: cap
        integer, allocatable :: order(:)
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
        call self%complex_ldlt_factor%free()
        call self%hpd_factor%free()
        call self%ldlt_factor%free()
        self%mode = 1
        self%is_complex = .true.
        threshold = self%options%pivot_tolerance
        if (threshold < 0.0_dp) threshold = 1.0_dp
        self%n = ncol
        started = seconds()
        call reserve_original(self,ncol,nz,refine,.true.,cap,status)
        if (.not.status_ok(status)) return
        call general_order(self,ncol,col_ptr,row_idx,.true.,order,cap,status, zvalues=val)
        if (.not.status_ok(status)) return
        started = seconds()
        self%factor_invoked=.true.
        self%retained_order=-1
        call self%complex_factor%factor_raw(ncol,col_ptr,row_idx,val,status, &
            column_permutation=order,pivot_threshold=threshold,memory_limit_bytes=cap, &
            pruning=.not.self%reference_only, row_preference=self%matching_preference, &
            low_memory=self%options%memory_policy == FORTSPARSE_MEMORY_LOW)
        self%stats%factor_s = seconds()-started
        if (allocated(order)) deallocate(order)
        call release_matching(self)
        self%ready = status_ok(status)
        if (self%ready) self%retained_order=self%selected_order
        if (self%ready) call retain_complex(self,ncol,col_ptr,row_idx,val,refine,status)
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
        if (status_ok(status)) reused = self%complex_factor%last_reused_pattern &
            .and. self%same_ordering_strategy
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
        self%stats%refinement_steps = 0
        self%stats%backward_error = -1.0_dp
        call solve_shape(self,size(b),size(x),.true.,status,operation)
        if (.not.status_ok(status)) return
        if (self%refining) self%saved_zrhs = b
        x = b
        call solve_complex_values(self%complex_factor,self%hpd_factor,self%complex_ldlt_factor, &
            self%mode,x,operation,status)
        if (status_ok(status) .and. self%refining) call refine_complex(self,x,operation,status)
    end subroutine general_solve_complex_op

    subroutine general_solve_complex_inplace(self,b,status)
        class(native_general_backend_t), intent(inout) :: self
        complex(dp), intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out) :: status
        self%stats%refinement_steps = 0
        self%stats%backward_error = -1.0_dp
        call solve_shape(self,size(b),size(b),.true.,status)
        if (.not.status_ok(status)) return
        if (self%refining) self%saved_zrhs = b
        call solve_complex_values(self%complex_factor,self%hpd_factor,self%complex_ldlt_factor, &
            self%mode,b,0,status)
        if (status_ok(status) .and. self%refining) call refine_complex(self,b,0,status)
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
        integer :: j,steps
        real(dp) :: error
        self%stats%refinement_steps = 0
        self%stats%backward_error = -1.0_dp
        call solve_shape(self,size(b,1),size(x,1),.true.,status,operation)
        if (.not.status_ok(status)) return
        if (size(b,2) /= size(x,2)) then
            call status_set(status,FORTSPARSE_INVALID_MATRIX,'native CPU: RHS count mismatch')
            return
        end if
        steps = 0
        error = -1.0_dp
        do j = 1,size(b,2)
            call self%solve_complex_op(b(:,j),x(:,j),operation,status)
            steps = max(steps,self%stats%refinement_steps)
            error = max(error,self%stats%backward_error)
            if (.not.status_ok(status)) return
        end do
        self%stats%refinement_steps = steps
        self%stats%backward_error = error
    end subroutine general_solve_complex_block_op

    subroutine solve_shape(self,nb,nx,is_complex,status,operation)
        class(native_general_backend_t), intent(in) :: self
        integer, intent(in) :: nb,nx
        logical, intent(in) :: is_complex
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: operation
        integer :: n
        if (present(operation)) then
            if (operation < 0 .or. operation > 2) then
                call status_set(status,FORTSPARSE_INVALID_MATRIX,'native CPU: invalid operation')
                return
            end if
        end if
        call status_set(status,FORTSPARSE_NOT_FACTORED,'native LU: no factors for this value type')
        if (.not. self%ready) return
        if (self%is_complex .neqv. is_complex) return
        n = self%real_factor%n
        if (self%is_complex) n = self%complex_factor%n
        if (self%mode == 2) n = self%hpd_factor%dimension()
        if (self%mode >= 3) n = self%n
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
            +self%ldlt_factor%retained_bytes() &
            +self%complex_ldlt_factor%retained_bytes()+vector_bytes &
            +original_bytes(self)+matching_bytes(self)
        if (.not.self%bounded_amd .or. self%factor_invoked) then
            self%stats%host_peak_bytes = max(self%stats%host_peak_bytes, &
                self%real_factor%peak_memory_bytes+self%complex_factor%peak_memory_bytes &
                +self%hpd_factor%memory_peak_bytes()+self%ldlt_factor%peak_bytes() &
                +self%complex_ldlt_factor%peak_bytes()+self%factor_external_live, &
                self%stats%host_current_bytes)
        else
            self%stats%host_peak_bytes=max(self%stats%host_peak_bytes,self%stats%host_current_bytes)
        end if
        self%stats%factor_nnz = self%real_factor%factor_nnz+self%complex_factor%factor_nnz &
            +self%hpd_factor%factor_nnz()+self%ldlt_factor%factor_nnz()+self%complex_ldlt_factor%factor_nnz()
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
        if (self%mode == 4 .or. self%mode == 5) then
            self%stats%budget_required_bytes = self%complex_ldlt_factor%required_bytes()
            self%stats%budget_available_bytes = self%complex_ldlt_factor%available_bytes()
        end if
        if (self%failed_required > 0_i8) then
            self%stats%budget_required_bytes = self%failed_required
            self%stats%budget_available_bytes = self%failed_available
        end if
        self%stats%factor_mem_mb = real(self%stats%host_peak_bytes,dp)/1.0e6_dp
        self%stats%strict_memory_supported = .true.
        self%stats%active_device = -1
        self%stats%device_bytes = 0_i8
        select case(self%mode)
        case(0)
            self%stats%native_factor_family = 'LU_REAL'
        case(1)
            self%stats%native_factor_family = 'LU_COMPLEX'
        case(2)
            self%stats%native_factor_family = 'CHOLESKY_HPD'
        case(3)
            self%stats%native_factor_family = 'LDLT_REAL'
        case(4)
            self%stats%native_factor_family = 'LDLH_COMPLEX'
        case(5)
            self%stats%native_factor_family = 'LDLT_COMPLEX'
        end select
        self%stats%active_threads = 1
        self%stats%nthreads = self%options%nthreads
        self%stats%ordering_used = self%selected_order
        if (self%selected_order == FORTSPARSE_ORDER_ND .or. &
            (self%selected_order == FORTSPARSE_ORDER_AMD .and. .not.self%bounded_amd) .or. &
            self%converted_raw) then
            ! Opaque graph ordering or uncapped raw conversion prevents a complete peak or
            ! strict-cap claim. Retained owned payload remains exactly counted.
            self%stats%host_peak_bytes = -1_i8
            self%stats%factor_mem_mb = -1.0_dp
            self%stats%strict_memory_supported = .false.
        end if
    end subroutine update_stats

    integer(i8) function original_bytes(self) result(bytes)
        class(native_general_backend_t), intent(in) :: self
        integer(i8), parameter :: ib=int(storage_size(0)/8,i8)
        bytes=0_i8
        if (allocated(self%original_ptr)) bytes=bytes+ib*size(self%original_ptr,kind=i8)
        if (allocated(self%original_rows)) bytes=bytes+ib*size(self%original_rows,kind=i8)
        if (allocated(self%original_values)) bytes=bytes+8_i8*size(self%original_values,kind=i8)
        if (allocated(self%saved_rhs)) bytes=bytes+8_i8*size(self%saved_rhs,kind=i8)
        if (allocated(self%residual)) bytes=bytes+8_i8*size(self%residual,kind=i8)
        if (allocated(self%original_zvalues)) bytes=bytes+16_i8*size(self%original_zvalues,kind=i8)
        if (allocated(self%saved_zrhs)) bytes=bytes+16_i8*size(self%saved_zrhs,kind=i8)
        if (allocated(self%zresidual)) bytes=bytes+16_i8*size(self%zresidual,kind=i8)
    end function original_bytes

    integer(i8) function current_bytes(self) result(bytes)
        class(native_general_backend_t), intent(in) :: self
        bytes=self%real_factor%memory_bytes()+self%complex_factor%memory_bytes() &
            +self%hpd_factor%memory_current_bytes()+self%ldlt_factor%retained_bytes() &
            +self%complex_ldlt_factor%retained_bytes()+original_bytes(self)+matching_bytes(self)
        if (associated(self%user_vector)) bytes=bytes+8_i8*size(self%user_vector,kind=i8)
    end function current_bytes

    subroutine release_original(self)
        class(native_general_backend_t), intent(inout) :: self
        if (allocated(self%original_ptr)) deallocate(self%original_ptr)
        if (allocated(self%original_rows)) deallocate(self%original_rows)
        if (allocated(self%original_values)) deallocate(self%original_values)
        if (allocated(self%saved_rhs)) deallocate(self%saved_rhs)
        if (allocated(self%residual)) deallocate(self%residual)
        if (allocated(self%original_zvalues)) deallocate(self%original_zvalues)
        if (allocated(self%saved_zrhs)) deallocate(self%saved_zrhs)
        if (allocated(self%zresidual)) deallocate(self%zresidual)
        self%refining=.false.
        self%matrix_norm=0.0_dp
        self%transpose_norm=0.0_dp
    end subroutine release_original

    subroutine reserve_original(self,n,nz,refine,is_complex,cap,status)
        class(native_general_backend_t), intent(inout) :: self
        integer, intent(in) :: n,nz
        logical, intent(in) :: refine,is_complex
        integer(i8), intent(inout) :: cap
        type(fortsparse_status_t), intent(out) :: status
        integer(i8) :: planned,value_bytes
        logical :: same
        self%refining=.false.
        call status_set(status,FORTSPARSE_INVALID_MATRIX,'native refinement: invalid CSC dimensions')
        if (n<0 .or. n==huge(0) .or. nz<0) return
        same=.false.
        if (allocated(self%original_ptr)) then
            if (size(self%original_ptr)==n+1 .and. size(self%original_rows)==nz) then
                if (is_complex) then
                    same=allocated(self%original_zvalues)
                else
                    same=allocated(self%original_values)
                end if
            end if
        end if
        if (.not.refine .or. .not.same) call release_original(self)
        planned=0_i8
        if (refine) then
            value_bytes=8_i8
            if (is_complex) value_bytes=16_i8
            call refinement_payload_bytes(int(n,i8),int(nz,i8),value_bytes,planned,status)
            if (.not.status_ok(status)) return
        end if
        self%refinement_quota=planned
        ! Reserving a positive quota must never turn an exhausted capped budget
        ! into the factors' zero-means-unlimited setting.
        if (cap>0_i8) then
            if (cap<=planned) then
                self%failed_required=planned+1_i8
                self%failed_available=cap
                call status_set(status,FORTSPARSE_RESOURCE_LIMIT, &
                    'native refinement: original CSC/workspace leaves no factor budget')
                call update_stats(self)
                return
            end if
            cap=cap-planned
        end if
        self%factor_external_live=original_bytes(self)
        if (associated(self%user_vector)) &
            self%factor_external_live=self%factor_external_live+8_i8*size(self%user_vector,kind=i8)
        call status_set(status,FORTSPARSE_OK,'')
    end subroutine reserve_original

    subroutine retain_real(self,n,ptr,rows,values,refine,status)
        class(native_general_backend_t), intent(inout) :: self
        integer, intent(in) :: n,ptr(:),rows(:)
        real(dp), intent(in) :: values(:)
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        integer :: ierr
        call status_set(status,FORTSPARSE_OK,'')
        if (.not.refine) return
        if (.not.allocated(self%original_ptr)) then
            allocate(self%original_ptr(n+1),self%original_rows(size(rows)), &
                self%original_values(size(values)),self%saved_rhs(n),self%residual(n),stat=ierr)
            if (ierr/=0) then
                self%stats%host_peak_bytes=max(self%stats%host_peak_bytes,current_bytes(self))
                call release_original(self)
                call status_set(status,FORTSPARSE_OUT_OF_MEMORY,'native refinement: original CSC allocation failed')
                return
            end if
        end if
        self%original_ptr=ptr
        self%original_rows=rows
        self%original_values=values
        call native_original_norms(n,ptr,rows,values,self%residual, &
            self%matrix_norm,self%transpose_norm,status)
        if (.not.status_ok(status)) return
        self%refining=.true.
    end subroutine retain_real

    subroutine retain_complex(self,n,ptr,rows,values,refine,status)
        class(native_general_backend_t), intent(inout) :: self
        integer, intent(in) :: n,ptr(:),rows(:)
        complex(dp), intent(in) :: values(:)
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        integer :: ierr
        call status_set(status,FORTSPARSE_OK,'')
        if (.not.refine) return
        if (.not.allocated(self%original_ptr)) then
            allocate(self%original_ptr(n+1),self%original_rows(size(rows)), &
                self%original_zvalues(size(values)),self%saved_zrhs(n),self%zresidual(n),stat=ierr)
            if (ierr/=0) then
                self%stats%host_peak_bytes=max(self%stats%host_peak_bytes,current_bytes(self))
                call release_original(self)
                call status_set(status,FORTSPARSE_OUT_OF_MEMORY,'native refinement: original CSC allocation failed')
                return
            end if
        end if
        self%original_ptr=ptr
        self%original_rows=rows
        self%original_zvalues=values
        call native_original_norms(n,ptr,rows,values,self%zresidual, &
            self%matrix_norm,self%transpose_norm,status)
        if (.not.status_ok(status)) return
        self%refining=.true.
    end subroutine retain_complex

    subroutine solve_real_values(lu,ldlt,mode,x,operation,status)
        type(lu_factor_t), intent(inout) :: lu
        type(native_ldlt_t), intent(inout) :: ldlt
        integer, intent(in) :: mode,operation
        real(dp), intent(inout) :: x(:)
        type(fortsparse_status_t), intent(out) :: status
        if (mode==3) then
            call ldlt%solve_inplace(x,status,transpose=operation/=0)
        else
            call lu%solve(x,status,operation)
        end if
    end subroutine solve_real_values

    subroutine solve_complex_values(lu,hpd,ldlt,mode,x,operation,status)
        type(lu_factor_z_t), intent(inout) :: lu
        type(native_hpd_t), intent(inout) :: hpd
        type(native_ldlt_z_t), intent(inout) :: ldlt
        integer, intent(in) :: mode,operation
        complex(dp), intent(inout) :: x(:)
        type(fortsparse_status_t), intent(out) :: status
        if (mode==2) then
            call hpd%solve_inplace(x,status,operation)
        else if (mode==4 .or. mode==5) then
            call ldlt%solve_inplace(x,status,operation)
        else
            call lu%solve(x,status,operation)
        end if
    end subroutine solve_complex_values

    subroutine refine_real(self,x,operation,status)
        class(native_general_backend_t), intent(inout) :: self
        real(dp), intent(inout) :: x(:)
        integer, intent(in) :: operation
        type(fortsparse_status_t), intent(out) :: status
        integer :: step,i
        real(dp) :: an,error
        an=self%matrix_norm
        if (operation/=0) an=self%transpose_norm
        do step=0,3
            self%stats%refinement_steps=step
            call native_residual_error(self%original_ptr,self%original_rows,self%original_values, &
                self%saved_rhs,x,operation,an,self%residual,error)
            self%stats%backward_error=error
            if (error<0.0_dp) then
                call status_set(status,FORTSPARSE_SINGULAR,'native refinement: backward error exceeds numerical range')
                return
            end if
            if (error<=32.0_dp*epsilon(1.0_dp)) then
                call status_set(status,FORTSPARSE_OK,'')
                return
            end if
            if (step==3) exit
            call solve_real_values(self%real_factor,self%ldlt_factor,self%mode,self%residual,operation,status)
            if (.not.status_ok(status)) return
            do i=1,size(x)
                x(i)=spd_update(x(i),-1.0_dp,self%residual(i))
            end do
        end do
        call status_set(status,FORTSPARSE_SINGULAR,'native refinement: backward-error target not met in three corrections')
    end subroutine refine_real

    subroutine refine_complex(self,x,operation,status)
        class(native_general_backend_t), intent(inout) :: self
        complex(dp), intent(inout) :: x(:)
        integer, intent(in) :: operation
        type(fortsparse_status_t), intent(out) :: status
        integer :: step,i
        real(dp) :: an,error
        an=self%matrix_norm
        if (operation/=0) an=self%transpose_norm
        do step=0,3
            self%stats%refinement_steps=step
            call native_residual_error(self%original_ptr,self%original_rows,self%original_zvalues, &
                self%saved_zrhs,x,operation,an,self%zresidual,error)
            self%stats%backward_error=error
            if (error<0.0_dp) then
                call status_set(status,FORTSPARSE_SINGULAR,'native refinement: backward error exceeds numerical range')
                return
            end if
            if (error<=32.0_dp*epsilon(1.0_dp)) then
                call status_set(status,FORTSPARSE_OK,'')
                return
            end if
            if (step==3) exit
            call solve_complex_values(self%complex_factor,self%hpd_factor,self%complex_ldlt_factor, &
                self%mode,self%zresidual,operation,status)
            if (.not.status_ok(status)) return
            do i=1,size(x)
                x(i)=lu_z_update(x(i),cmplx(-1.0_dp,0.0_dp,dp),self%zresidual(i))
            end do
        end do
        call status_set(status,FORTSPARSE_SINGULAR,'native refinement: backward-error target not met in three corrections')
    end subroutine refine_complex

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
            +self%complex_ldlt_factor%retained_bytes()+original_bytes(self)+8_i8*int(n,i8)
        self%failed_required = 0_i8
        self%failed_available = 0_i8
        cap = self%options%host_memory_limit_bytes
        if (cap > 0_i8 .and. required > cap) then
            self%failed_required = 8_i8*int(n,i8)
            self%failed_available = max(0_i8, cap-required+8_i8*int(n,i8))
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
        call self%complex_ldlt_factor%free()
        self%n = 0
        call release_matching(self)
        self%previous_matching=.false.
        self%selected_order = FORTSPARSE_ORDER_DEFAULT
        self%same_ordering_strategy = .false.
        call release_original(self)
        self%factor_external_live = 0_i8
        self%converted_raw=.false.
        self%retained_order=-1
        self%bounded_amd=.false.
        self%factor_invoked=.false.
        self%refinement_quota=0_i8
        self%failed_required = 0_i8
        self%failed_available = 0_i8
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
