! MIT. Thin ownership and validation boundary; arithmetic stays in FortSparse.
module fortsparse_native_c_api
    use, intrinsic :: iso_c_binding
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortsparse
    use fortsparse_memory, only: memory_checked_add, memory_checked_mul
    implicit none
    private

    type, bind(C) :: memory_info_t
        integer(c_int64_t) :: host_limit_bytes = 0, device_limit_bytes = 0
        integer(c_int64_t) :: host_current_bytes = 0, host_peak_bytes = 0
        integer(c_int64_t) :: device_current_bytes = 0, device_peak_bytes = 0
        integer(c_int64_t) :: budget_required_bytes = 0, budget_available_bytes = 0
        integer(c_int32_t) :: memory_policy = 0, strict_memory_supported = 1, peak_is_bound = 0
    end type

    type :: context_t
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        integer :: scalar_kind = 0
        integer :: n = 0
        type(memory_info_t) :: memory
    end type

contains

    function fs_native_create(scalar_kind, handle_arg) bind(C) result(code)
        integer(c_int32_t), value :: scalar_kind
        type(c_ptr), value :: handle_arg
        type(c_ptr), pointer :: handle
        integer(c_int32_t) :: code
        type(context_t), pointer :: ctx
        integer :: alloc_status

        code = FORTSPARSE_INVALID_MATRIX
        if (.not. c_associated(handle_arg)) return
        call c_f_pointer(handle_arg, handle)
        if (c_associated(handle)) return
        if (scalar_kind /= 0 .and. scalar_kind /= 1) return
        allocate(ctx, stat=alloc_status)
        code = FORTSPARSE_OUT_OF_MEMORY
        if (alloc_status /= 0) return
        ctx%scalar_kind = scalar_kind
        ctx%solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        handle = c_loc(ctx)
        code = FORTSPARSE_OK
    end function

    function fs_native_free(handle_arg) bind(C) result(code)
        type(c_ptr), value :: handle_arg
        type(c_ptr), pointer :: handle
        integer(c_int32_t) :: code
        type(context_t), pointer :: ctx

        code = FORTSPARSE_OK
        if (.not. c_associated(handle_arg)) return
        call c_f_pointer(handle_arg, handle)
        if (c_associated(handle)) then
            call c_f_pointer(handle, ctx)
            call sparse_free(ctx%solver)
            deallocate(ctx)
            handle = c_null_ptr
        end if
        code = FORTSPARSE_OK
    end function

    function fs_native_set_refine(handle, refine) bind(C) result(code)
        type(c_ptr), value :: handle
        integer(c_int32_t), value :: refine
        integer(c_int32_t) :: code
        type(context_t), pointer :: ctx

        code = FORTSPARSE_INVALID_MATRIX
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, ctx)
        call status_set(ctx%status, FORTSPARSE_INVALID_MATRIX, 'refine must be 0 or 1')
        if (refine == 0 .or. refine == 1) then
            if (ctx%solver%refine .neqv. (refine == 1)) then
                call sparse_free(ctx%solver)
                ctx%n = 0
                call begin_memory(ctx)
            end if
            ctx%solver%refine = refine == 1
            ctx%memory%budget_required_bytes=0
            ctx%memory%budget_available_bytes=0
            call status_set(ctx%status, FORTSPARSE_OK, '')
        end if
        code = ctx%status%code
    end function

    function fs_native_ordering_available(ordering) bind(C) result(code)
        integer(c_int32_t), value :: ordering
        integer(c_int32_t) :: code

        code = FORTSPARSE_INVALID_MATRIX
        if (ordering < FORTSPARSE_ORDER_DEFAULT .or. &
            ordering > FORTSPARSE_ORDER_NATURAL) return
        code = FORTSPARSE_BACKEND_UNAVAILABLE
        if (native_ordering_available(int(ordering))) code = FORTSPARSE_OK
    end function

    function fs_native_set_ordering(handle, ordering) bind(C) result(code)
        type(c_ptr), value :: handle
        integer(c_int32_t), value :: ordering
        integer(c_int32_t) :: code
        type(context_t), pointer :: ctx

        code = FORTSPARSE_INVALID_MATRIX
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, ctx)
        code = fs_native_ordering_available(ordering)
        if (code == FORTSPARSE_INVALID_MATRIX) then
            call status_set(ctx%status, code, 'invalid native ordering enum')
            return
        end if
        if (code == FORTSPARSE_BACKEND_UNAVAILABLE) then
            call status_set(ctx%status, code, 'native general LU ordering unavailable')
            return
        end if
        if (ctx%solver%options%ordering /= ordering) then
            call sparse_free(ctx%solver)
            ctx%n = 0
            call begin_memory(ctx)
        end if
        ctx%solver%options%ordering = ordering
        ctx%memory%budget_required_bytes=0
        ctx%memory%budget_available_bytes=0
        call status_set(ctx%status, FORTSPARSE_OK, '')
    end function

    function fs_native_set_memory(handle, host_limit, device_limit, policy) bind(C) result(code)
        type(c_ptr), value :: handle
        integer(c_int64_t), value :: host_limit, device_limit
        integer(c_int32_t), value :: policy
        integer(c_int32_t) :: code
        type(context_t), pointer :: ctx
        code = FORTSPARSE_INVALID_MATRIX
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, ctx)
        call status_set(ctx%status, code, 'invalid native memory limits or policy')
        if (host_limit < 0 .or. device_limit < 0) return
        if (policy /= FORTSPARSE_MEMORY_AUTO .and. policy /= FORTSPARSE_MEMORY_LOW) return
        if (host_limit /= ctx%memory%host_limit_bytes .or. &
            device_limit /= ctx%memory%device_limit_bytes .or. policy /= ctx%memory%memory_policy) then
            call sparse_free(ctx%solver)
            ctx%n = 0
            ctx%memory = memory_info_t()
        end if
        ctx%memory%host_limit_bytes = host_limit
        ctx%memory%device_limit_bytes = device_limit
        ctx%memory%memory_policy = policy
        ctx%memory%budget_required_bytes = 0
        ctx%memory%budget_available_bytes = 0
        ctx%solver%options%host_memory_limit_bytes = host_limit
        ctx%solver%options%device_memory_limit_bytes = device_limit
        ctx%solver%options%memory_policy = policy
        call status_set(ctx%status, FORTSPARSE_OK, '')
        code = FORTSPARSE_OK
    end function

    function fs_native_get_memory_stats(handle, output_arg) bind(C) result(code)
        type(c_ptr), value :: handle, output_arg
        integer(c_int32_t) :: code
        type(context_t), pointer :: ctx
        type(memory_info_t), pointer :: output
        type(sparse_factor_stats_t) :: stats
        code = FORTSPARSE_INVALID_MATRIX
        if (.not. c_associated(handle) .or. .not. c_associated(output_arg)) return
        call c_f_pointer(handle, ctx)
        call c_f_pointer(output_arg, output)
        call backend_memory(ctx, stats)
        ctx%memory%host_current_bytes = stats%host_current_bytes
        ! ND remains opaque even if a strict request failed before graph setup.
        if (ctx%solver%options%ordering == FORTSPARSE_ORDER_ND) &
            ctx%memory%strict_memory_supported = 0
        output = ctx%memory
        code = FORTSPARSE_OK
    end function

    subroutine backend_memory(ctx, stats)
        type(context_t), intent(inout) :: ctx
        type(sparse_factor_stats_t), intent(out) :: stats
        type(fortsparse_status_t) :: query_status
        stats = sparse_factor_stats_t()
        stats%host_current_bytes = 0
        stats%host_peak_bytes = 0
        stats%strict_memory_supported = .true.
        if (allocated(ctx%solver%backend)) call sparse_factor_stats(ctx%solver, stats, query_status)
    end subroutine

    subroutine begin_memory(ctx)
        type(context_t), intent(inout) :: ctx
        type(sparse_factor_stats_t) :: stats
        call backend_memory(ctx, stats)
        ctx%memory%host_current_bytes = stats%host_current_bytes
        ctx%memory%host_peak_bytes = stats%host_current_bytes
        ctx%memory%peak_is_bound = 0
        ctx%memory%budget_required_bytes = 0
        ctx%memory%budget_available_bytes = 0
        ctx%memory%strict_memory_supported = 0
        if (stats%strict_memory_supported) ctx%memory%strict_memory_supported = 1
    end subroutine

    subroutine staging_quota(ctx, bytes, backend_cap)
        type(context_t), intent(inout) :: ctx
        integer(i8), intent(in) :: bytes
        integer(i8), intent(out) :: backend_cap
        integer(i8) :: total, available
        backend_cap = 0
        call memory_checked_add(ctx%memory%host_current_bytes, bytes, total, ctx%status)
        if (.not. status_ok(ctx%status)) return
        if (ctx%memory%host_limit_bytes == 0) return
        available = max(0_i8, ctx%memory%host_limit_bytes - ctx%memory%host_current_bytes)
        if (bytes > available) then
            ctx%memory%budget_required_bytes = bytes
            ctx%memory%budget_available_bytes = available
            call status_set(ctx%status, FORTSPARSE_RESOURCE_LIMIT, 'native C staging exceeds host cap')
            return
        end if
        backend_cap = ctx%memory%host_limit_bytes - bytes
        if (backend_cap == 0) then
            ctx%memory%budget_required_bytes = 1
            ctx%memory%budget_available_bytes = 0
            call status_set(ctx%status, FORTSPARSE_RESOURCE_LIMIT, 'native C staging leaves no backend budget')
        end if
    end subroutine

    subroutine note_staging(ctx, bytes)
        type(context_t), intent(inout) :: ctx
        integer(i8), intent(in) :: bytes
        ctx%memory%host_peak_bytes = max(ctx%memory%host_peak_bytes, ctx%memory%host_current_bytes + bytes)
    end subroutine

    subroutine finish_factor_memory(ctx, bytes)
        type(context_t), intent(inout) :: ctx
        integer(i8), intent(in) :: bytes
        type(sparse_factor_stats_t) :: stats
        call backend_memory(ctx, stats)
        ctx%memory%host_current_bytes = stats%host_current_bytes
        ctx%memory%peak_is_bound = 1
        ctx%memory%strict_memory_supported = 0
        if (stats%strict_memory_supported) ctx%memory%strict_memory_supported = 1
        if (stats%host_peak_bytes < 0) then
            ctx%memory%host_peak_bytes = -1
        else
            ctx%memory%host_peak_bytes = max(ctx%memory%host_peak_bytes, stats%host_peak_bytes + bytes)
            ! Historical backend admission peaks can be conservative. An enforced
            ! complete cap is also a valid upper bound, not an observed RSS peak.
            if (stats%strict_memory_supported .and. ctx%memory%host_limit_bytes > 0) &
                ctx%memory%host_peak_bytes = min(ctx%memory%host_peak_bytes, ctx%memory%host_limit_bytes)
        end if
        ctx%memory%budget_required_bytes = stats%budget_required_bytes
        ctx%memory%budget_available_bytes = stats%budget_available_bytes
    end subroutine

    function fs_native_last_error(handle, message, capacity) bind(C) result(code)
        type(c_ptr), value :: handle, message
        integer(c_int32_t), value :: capacity
        integer(c_int32_t) :: code
        type(context_t), pointer :: ctx
        character(c_char), pointer :: output(:)
        integer :: i, count

        code = FORTSPARSE_INVALID_MATRIX
        if (.not. c_associated(handle)) return
        if (capacity < 1 .or. .not. c_associated(message)) return
        call c_f_pointer(handle, ctx)
        call c_f_pointer(message, output, [capacity])
        count = min(len_trim(ctx%status%msg), capacity-1)
        do i = 1, count
            output(i) = ctx%status%msg(i:i)
        end do
        output(count+1) = c_null_char
        code = ctx%status%code
    end function

    logical function check_csc(ctx, n, nnz, colptr, rowind, values, cp, ri)
        type(context_t), intent(inout) :: ctx
        integer(c_int32_t), intent(in) :: n, nnz
        type(c_ptr), intent(in) :: colptr, rowind, values
        integer(c_int32_t), pointer, contiguous, intent(out) :: cp(:), ri(:)
        integer :: j, k

        check_csc = .false.
        ! Every failure invalidates the old factors, including validation errors.
        ctx%solver%factored = .false.
        ctx%n = 0
        call status_set(ctx%status, FORTSPARSE_INVALID_MATRIX, 'invalid CSC dimensions/pointers')
        if (n < 1 .or. nnz < 0) return
        if (n == huge(n) .or. nnz == huge(nnz)) then
            call status_set(ctx%status, FORTSPARSE_INDEX_OVERFLOW, 'CSC index +1 overflows I4')
            return
        end if
        if (.not. c_associated(colptr)) return
        ! Zero-entry matrices still require non-NULL data pointers (see header).
        if (.not. c_associated(rowind) .or. .not. c_associated(values)) return
        call c_f_pointer(colptr, cp, [n+1])
        call c_f_pointer(rowind, ri, [nnz])
        call status_set(ctx%status, FORTSPARSE_INVALID_MATRIX, 'invalid noncanonical zero-based CSC')
        if (cp(1) /= 0 .or. cp(n+1) /= nnz) return
        do j = 1, n
            if (cp(j) < 0 .or. cp(j+1) > nnz) return
            if (cp(j) > cp(j+1)) return
            do k = cp(j)+1, cp(j+1)
                if (ri(k) < 0 .or. ri(k) >= n) return
                if (k > cp(j)+1) then
                    if (ri(k) <= ri(k-1)) return
                end if
            end do
        end do
        check_csc = .true.
    end function

    function fs_native_factor_d(handle, n, nnz, colptr, rowind, values) &
            bind(C) result(code)
        type(c_ptr), value :: handle, colptr, rowind, values
        integer(c_int32_t), value :: n, nnz
        integer(c_int32_t) :: code
        code = factor_d(handle, n, nnz, colptr, rowind, values, .false.)
    end function

    function fs_native_refactor_d(handle, n, nnz, colptr, rowind, values) &
            bind(C) result(code)
        type(c_ptr), value :: handle, colptr, rowind, values
        integer(c_int32_t), value :: n, nnz
        integer(c_int32_t) :: code
        code = factor_d(handle, n, nnz, colptr, rowind, values, .true.)
    end function

    function factor_d(handle, n, nnz, colptr, rowind, values, reuse) result(code)
        type(c_ptr), intent(in) :: handle, colptr, rowind, values
        integer(c_int32_t), intent(in) :: n, nnz
        logical, intent(in) :: reuse
        integer(c_int32_t) :: code
        type(context_t), pointer :: ctx
        integer(c_int32_t), pointer, contiguous :: cp(:), ri(:)
        real(c_double), pointer, contiguous :: v(:)
        type(csc_t) :: matrix
        integer :: alloc_status
        integer(i8) :: staging_bytes, backend_cap, part

        code = FORTSPARSE_INVALID_MATRIX
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, ctx)
        call begin_memory(ctx)
        ctx%solver%factored = .false.
        ctx%n = 0
        call status_set(ctx%status, FORTSPARSE_INVALID_MATRIX, 'wrong scalar kind')
        if (ctx%scalar_kind /= 0) return
        if (.not. check_csc(ctx, n, nnz, colptr, rowind, values, cp, ri)) then
            code = ctx%status%code
            return
        end if
        call c_f_pointer(values, v, [nnz])
        if (any(.not. ieee_is_finite(v))) then
            call status_set(ctx%status, FORTSPARSE_INVALID_MATRIX, 'nonfinite CSC value')
            return
        end if
        call memory_checked_mul(int(nnz, i8), 12_i8, part, ctx%status)
        if (.not. status_ok(ctx%status)) then
            code = ctx%status%code
            return
        end if
        call memory_checked_add(4_i8*(int(n, i8) + 1_i8), part, staging_bytes, ctx%status)
        if (.not. status_ok(ctx%status)) then
            code = ctx%status%code
            return
        end if
        call staging_quota(ctx, staging_bytes, backend_cap)
        if (.not. status_ok(ctx%status)) then
            code = ctx%status%code
            return
        end if
        allocate (matrix%col_ptr(n + 1), matrix%row_idx(nnz), matrix%val(nnz), stat=alloc_status)
        code = FORTSPARSE_OUT_OF_MEMORY
        if (alloc_status /= 0) then
            part = 0
            if (allocated(matrix%col_ptr)) part = part + 4_i8*size(matrix%col_ptr, kind=i8)
            if (allocated(matrix%row_idx)) part = part + 4_i8*size(matrix%row_idx, kind=i8)
            if (allocated(matrix%val)) part = part + 8_i8*size(matrix%val, kind=i8)
            call note_staging(ctx, part)
            call status_set(ctx%status, code, 'CSC conversion allocation failed')
            return
        end if
        call note_staging(ctx, staging_bytes)
        ctx%solver%options%host_memory_limit_bytes = backend_cap
        matrix%nrow = n
        matrix%ncol = n
        matrix%nnz = nnz
        matrix%col_ptr = cp + 1
        matrix%row_idx = ri + 1
        matrix%val = v
        if (reuse) then
            call sparse_refactor(ctx%solver, matrix, ctx%status)
        else
            call sparse_factor(ctx%solver, matrix, ctx%status)
        end if
        ctx%solver%options%host_memory_limit_bytes = ctx%memory%host_limit_bytes
        call finish_factor_memory(ctx, staging_bytes)
        code = ctx%status%code
        if (code == FORTSPARSE_OK) ctx%n = n
    end function

    function fs_native_solve_d(handle, operation, n, nrhs, b, x) bind(C) result(code)
        type(c_ptr), value :: handle, b, x
        integer(c_int32_t), value :: operation, n, nrhs
        integer(c_int32_t) :: code
        type(context_t), pointer :: ctx
        real(c_double), pointer, contiguous :: rhs(:, :), output(:, :)
        real(c_double), allocatable :: solution(:, :)
        integer :: alloc_status
        integer(i8) :: staging_bytes, unused_cap

        code = FORTSPARSE_INVALID_MATRIX
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, ctx)
        call begin_memory(ctx)
        call status_set(ctx%status, FORTSPARSE_INVALID_MATRIX, 'invalid solve shape, operation, or scalar kind')
        if (ctx%scalar_kind /= 0) return
        if (.not. ctx%solver%factored) then
            call status_set(ctx%status, FORTSPARSE_NOT_FACTORED, 'native handle is not factored')
            code = ctx%status%code
            return
        end if
        if (n /= ctx%n .or. nrhs < 1) return
        if (operation < 0 .or. operation > 2) return
        if (int(n, c_int64_t)*int(nrhs, c_int64_t) > huge(0_c_int32_t)) then
            call status_set(ctx%status, FORTSPARSE_INDEX_OVERFLOW, 'RHS shape exceeds I4 array extent')
            code = ctx%status%code
            return
        end if
        if (.not. c_associated(b) .or. .not. c_associated(x)) return
        call c_f_pointer(b, rhs, [n, nrhs])
        call c_f_pointer(x, output, [n, nrhs])
        if (any(.not. ieee_is_finite(rhs))) then
            call status_set(ctx%status, FORTSPARSE_INVALID_MATRIX, 'nonfinite RHS')
            return
        end if
        call memory_checked_mul(8_i8*int(n, i8), int(nrhs, i8), staging_bytes, ctx%status)
        if (status_ok(ctx%status)) call staging_quota(ctx, staging_bytes, unused_cap)
        if (.not. status_ok(ctx%status)) then
            code = ctx%status%code
            return
        end if
        allocate (solution(n, nrhs), stat=alloc_status)
        code = FORTSPARSE_OUT_OF_MEMORY
        if (alloc_status /= 0) then
            call status_set(ctx%status, code, 'solution staging allocation failed')
            return
        end if
        call note_staging(ctx, staging_bytes)
        call sparse_solve_block_op(ctx%solver, rhs, solution, operation, ctx%status)
        code = ctx%status%code
        if (code == FORTSPARSE_OK) then
            if (any(.not. ieee_is_finite(solution))) then
                call status_set(ctx%status, FORTSPARSE_INTERNAL_ERROR, 'nonfinite solution')
                code = ctx%status%code
                return
            end if
            output = solution
        end if
    end function

    function fs_native_factor_z(handle, n, nnz, colptr, rowind, values) &
            bind(C) result(code)
        type(c_ptr), value :: handle, colptr, rowind, values
        integer(c_int32_t), value :: n, nnz
        integer(c_int32_t) :: code
        code = factor_z(handle, n, nnz, colptr, rowind, values, .false.)
    end function

    function fs_native_refactor_z(handle, n, nnz, colptr, rowind, values) &
            bind(C) result(code)
        type(c_ptr), value :: handle, colptr, rowind, values
        integer(c_int32_t), value :: n, nnz
        integer(c_int32_t) :: code
        code = factor_z(handle, n, nnz, colptr, rowind, values, .true.)
    end function

    function factor_z(handle, n, nnz, colptr, rowind, values, reuse) result(code)
        type(c_ptr), intent(in) :: handle, colptr, rowind, values
        integer(c_int32_t), intent(in) :: n, nnz
        logical, intent(in) :: reuse
        integer(c_int32_t) :: code
        type(context_t), pointer :: ctx
        integer(c_int32_t), pointer, contiguous :: cp(:), ri(:)
        complex(c_double_complex), pointer, contiguous :: v(:)
        type(csc_z_t) :: matrix
        integer :: alloc_status
        integer(i8) :: staging_bytes, backend_cap, part

        code = FORTSPARSE_INVALID_MATRIX
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, ctx)
        call begin_memory(ctx)
        ctx%solver%factored = .false.
        ctx%n = 0
        call status_set(ctx%status, FORTSPARSE_INVALID_MATRIX, 'wrong scalar kind')
        if (ctx%scalar_kind /= 1) return
        if (.not. check_csc(ctx, n, nnz, colptr, rowind, values, cp, ri)) then
            code = ctx%status%code
            return
        end if
        call c_f_pointer(values, v, [nnz])
        if (any(.not. ieee_is_finite(real(v))) .or. &
            any(.not. ieee_is_finite(aimag(v)))) then
            call status_set(ctx%status, FORTSPARSE_INVALID_MATRIX, 'nonfinite CSC value')
            return
        end if
        call memory_checked_mul(int(nnz, i8), 20_i8, part, ctx%status)
        if (.not. status_ok(ctx%status)) then
            code = ctx%status%code
            return
        end if
        call memory_checked_add(4_i8*(int(n, i8) + 1_i8), part, staging_bytes, ctx%status)
        if (.not. status_ok(ctx%status)) then
            code = ctx%status%code
            return
        end if
        call staging_quota(ctx, staging_bytes, backend_cap)
        if (.not. status_ok(ctx%status)) then
            code = ctx%status%code
            return
        end if
        allocate (matrix%col_ptr(n + 1), matrix%row_idx(nnz), matrix%val(nnz), stat=alloc_status)
        code = FORTSPARSE_OUT_OF_MEMORY
        if (alloc_status /= 0) then
            part = 0
            if (allocated(matrix%col_ptr)) part = part + 4_i8*size(matrix%col_ptr, kind=i8)
            if (allocated(matrix%row_idx)) part = part + 4_i8*size(matrix%row_idx, kind=i8)
            if (allocated(matrix%val)) part = part + 16_i8*size(matrix%val, kind=i8)
            call note_staging(ctx, part)
            call status_set(ctx%status, code, 'CSC conversion allocation failed')
            return
        end if
        call note_staging(ctx, staging_bytes)
        ctx%solver%options%host_memory_limit_bytes = backend_cap
        matrix%nrow = n
        matrix%ncol = n
        matrix%nnz = nnz
        matrix%col_ptr = cp + 1
        matrix%row_idx = ri + 1
        matrix%val = v
        if (reuse) then
            call sparse_refactor(ctx%solver, matrix, ctx%status)
        else
            call sparse_factor(ctx%solver, matrix, ctx%status)
        end if
        ctx%solver%options%host_memory_limit_bytes = ctx%memory%host_limit_bytes
        call finish_factor_memory(ctx, staging_bytes)
        code = ctx%status%code
        if (code == FORTSPARSE_OK) ctx%n = n
    end function

    function fs_native_solve_z(handle, operation, n, nrhs, b, x) bind(C) result(code)
        type(c_ptr), value :: handle, b, x
        integer(c_int32_t), value :: operation, n, nrhs
        integer(c_int32_t) :: code
        type(context_t), pointer :: ctx
        complex(c_double_complex), pointer, contiguous :: rhs(:, :), output(:, :)
        complex(c_double_complex), allocatable :: solution(:, :)
        integer :: alloc_status
        integer(i8) :: staging_bytes, unused_cap

        code = FORTSPARSE_INVALID_MATRIX
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, ctx)
        call begin_memory(ctx)
        call status_set(ctx%status, FORTSPARSE_INVALID_MATRIX, 'invalid solve shape, operation, or scalar kind')
        if (ctx%scalar_kind /= 1) return
        if (.not. ctx%solver%factored) then
            call status_set(ctx%status, FORTSPARSE_NOT_FACTORED, 'native handle is not factored')
            code = ctx%status%code
            return
        end if
        if (n /= ctx%n .or. nrhs < 1) return
        if (operation < 0 .or. operation > 2) return
        if (int(n, c_int64_t)*int(nrhs, c_int64_t) > huge(0_c_int32_t)) then
            call status_set(ctx%status, FORTSPARSE_INDEX_OVERFLOW, 'RHS shape exceeds I4 array extent')
            code = ctx%status%code
            return
        end if
        if (.not. c_associated(b) .or. .not. c_associated(x)) return
        call c_f_pointer(b, rhs, [n, nrhs])
        call c_f_pointer(x, output, [n, nrhs])
        if (any(.not. ieee_is_finite(real(rhs))) .or. &
            any(.not. ieee_is_finite(aimag(rhs)))) then
            call status_set(ctx%status, FORTSPARSE_INVALID_MATRIX, 'nonfinite RHS')
            return
        end if
        call memory_checked_mul(16_i8*int(n, i8), int(nrhs, i8), staging_bytes, ctx%status)
        if (status_ok(ctx%status)) call staging_quota(ctx, staging_bytes, unused_cap)
        if (.not. status_ok(ctx%status)) then
            code = ctx%status%code
            return
        end if
        allocate (solution(n, nrhs), stat=alloc_status)
        code = FORTSPARSE_OUT_OF_MEMORY
        if (alloc_status /= 0) then
            call status_set(ctx%status, code, 'solution staging allocation failed')
            return
        end if
        call note_staging(ctx, staging_bytes)
        call sparse_solve_block_op(ctx%solver, rhs, solution, operation, ctx%status)
        code = ctx%status%code
        if (code == FORTSPARSE_OK) then
            if (any(.not. ieee_is_finite(real(solution))) .or. &
                any(.not. ieee_is_finite(aimag(solution)))) then
                call status_set(ctx%status, FORTSPARSE_INTERNAL_ERROR, 'nonfinite solution')
                code = ctx%status%code
                return
            end if
            output = solution
        end if
    end function

end module fortsparse_native_c_api
