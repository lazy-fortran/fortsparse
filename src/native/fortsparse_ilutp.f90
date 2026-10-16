module fortsparse_ilutp
    ! Independently written real CSR ILUT with dual dropping and column pivots.
    ! Mathematical/provenance and payload contracts: docs/ilutp.md.
    use fortsparse_kinds, only: dp, i8
    use fortsparse_status, only: fortsparse_status_t, status_set, &
                                 FORTSPARSE_OK, &
                                 FORTSPARSE_INVALID_MATRIX, FORTSPARSE_SINGULAR, &
                                 FORTSPARSE_NOT_FACTORED, &
                                 FORTSPARSE_OUT_OF_MEMORY, &
                                 FORTSPARSE_RESOURCE_LIMIT, &
                                 FORTSPARSE_INDEX_OVERFLOW, &
                                 FORTSPARSE_INTERNAL_ERROR
    use fortsparse_spd_leaves, only: divide_value => spd_divide, &
                                     update_value => spd_update
    use fortsparse_ilutp_selection, only: lower_push, lower_pop, &
                                          select_largest
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: sparse_ilutp_t, sparse_ilutp_stats_t
    integer(i8), parameter :: ib = storage_size(0)/8, &
                              vb = storage_size(0.0_dp)/8

    type :: sparse_ilutp_stats_t
        logical :: ready = .false.
        integer :: n = 0, lower_nnz = 0, upper_nnz = 0
        integer :: column_swaps = 0, regularized_pivots = 0
        integer(i8) :: duplicate_entries = 0, stored_zero_entries = 0
        integer(i8) :: lower_threshold_drops = 0, upper_threshold_drops = 0
        integer(i8) :: lower_fill_drops = 0, upper_fill_drops = 0
        integer(i8) :: current_bytes = 0, peak_owned_bytes = 0
        integer(i8) :: setup_required_bytes = 0, failed_required_bytes = 0
        integer(i8) :: available_bytes = 0
    end type

    type :: sparse_ilutp_t
        private
        type(sparse_ilutp_stats_t) :: info
        integer, allocatable :: lp(:), up(:), lc(:), uc(:), &
                                perm(:), inverse(:)
        integer, allocatable :: marker(:), touched(:), heap(:)
        real(dp), allocatable :: lv(:), uv(:), diagonal(:), work(:)
    contains
        procedure :: setup, apply, stats
        procedure :: free => release_context
        procedure, private :: memory_bytes
    end type
contains
    subroutine setup(self, n, row_ptr, columns, values, &
                     lower_fill, upper_fill, &
                     drop_tolerance, pivot_tolerance, status, &
                     host_memory_limit_bytes)
        class(sparse_ilutp_t), intent(inout) :: self
        integer, intent(in) :: n, row_ptr(:), columns(:), &
                               lower_fill, upper_fill
        real(dp), intent(in) :: values(:), drop_tolerance, pivot_tolerance
        type(fortsparse_status_t), intent(out) :: status
        integer(i8), intent(in), optional :: host_memory_limit_bytes
        integer(i8) :: lower_capacity, upper_capacity, slots, &
                       limit, required
        integer :: i, j, k, p, column, pivot, old_column, io, nz
        integer :: count, pending, kept, eligible, nl, nu
        real(dp) :: scale, mean, multiplier, threshold, value, &
                    diagonal_value

        ! Rebuild is destructive: release every old array before validation or
        ! new allocation. There is deliberately no hidden old/new factor overlap.
        call self%free()
        limit = 0_i8
        if (present(host_memory_limit_bytes)) limit = &
            host_memory_limit_bytes
        self%info%available_bytes = limit
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        'invalid ILUTP CSR/options')
        if (limit < 0_i8) return
        if (n < 1 .or. lower_fill < 0 .or. upper_fill < 0) return
        if (n == huge(n)) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                            'ILUTP row pointer overflow')
            return
        end if
        if (.not. ieee_is_finite(drop_tolerance)) return
        if (.not. ieee_is_finite(pivot_tolerance)) return
        if (drop_tolerance < 0.0_dp .or. drop_tolerance > 1.0_dp) return
        if (pivot_tolerance < 0.0_dp .or. pivot_tolerance > 1.0_dp) return
        if (size(row_ptr, kind=i8) /= int(n, i8) + 1_i8) return
        if (size(columns, kind=i8) /= size(values, kind=i8)) return
        if (size(values, kind=i8) >= 2147483647_i8) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                            'ILUTP input pointer overflow')
            return
        end if
        nz = int(size(values, kind=i8))
        if (row_ptr(1) /= 1 .or. row_ptr(n + 1) /= nz + 1) return
        do i = 1, n
            if (row_ptr(i) < 1 .or. row_ptr(i + 1) > nz + 1) return
            if (row_ptr(i) > row_ptr(i + 1)) return
        end do
        do p = 1, nz
            if (columns(p) < 1 .or. columns(p) > n) return
            if (.not. ieee_is_finite(values(p))) return
        end do
        ! Default-integer pointer capacity is checked before any allocation.
        ! n,k are bounded by I4 below, making the i8 products nonoverflowing.
        if (int(n, i8) > 2147483646_i8) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                            'ILUTP requires I4 dimensions')
            return
        end if
        k = min(lower_fill, n - 1)
        lower_capacity = int(k, i8)*int(n, i8) - int(k, &
                                                     i8)*(int(k, i8) + 1)/2
        k = min(upper_fill, n - 1)
        upper_capacity = int(k, i8)*int(n, i8) - int(k, &
                                                     i8)*(int(k, i8) + 1)/2
        if (max(lower_capacity, upper_capacity) >= 2147483647_i8) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                            'ILUTP factor pointer overflow')
            return
        end if
        slots = lower_capacity + upper_capacity
        required = ib*(slots + 7_i8*n + 2_i8) + vb*(slots + 2_i8*n)
        self%info%setup_required_bytes = required
        if (limit > 0_i8) then
            if (required > limit) then
                self%info%failed_required_bytes = required
                call status_set(status, FORTSPARSE_RESOURCE_LIMIT, &
                                'ILUTP setup array payload exceeds host limit')
                return
            end if
        end if
        allocate (self%lp(n + 1), self%up(n + 1), self%lc(lower_capacity), &
                  self%uc(upper_capacity), self%lv(lower_capacity), &
                  self%uv(upper_capacity), self%diagonal(n), self%work(n), &
                  self%perm(n), self%inverse(n), self%marker(n), &
                  self%touched(n), self%heap(n), stat=io)
        if (io /= 0) then
            call fail(FORTSPARSE_OUT_OF_MEMORY, 'ILUTP allocation failed')
            return
        end if
        self%info%peak_owned_bytes = self%memory_bytes()
        self%info%n = n
        self%marker = 0
        self%lp(1) = 1
        self%up(1) = 1
        nl = 0
        nu = 0
        do i = 1, n
            self%perm(i) = i
            self%inverse(i) = i
        end do
        do i = 1, n
            count = 0
            pending = 0
            scale = 0.0_dp
            do p = row_ptr(i), row_ptr(i + 1) - 1
                scale = max(scale, abs(values(p)))
            end do
            if (scale <= 0.0_dp) then
                call fail(FORTSPARSE_SINGULAR, &
                          'ILUTP encountered an empty/all-zero input row')
                return
            end if
            mean = 0.0_dp
            do p = row_ptr(i), row_ptr(i + 1) - 1
                ! Scale the row reduction to avoid overflow of sum(abs(a)).
                mean = update_value(mean, &
                                    -divide_value(abs(values(p)), scale), 1.0_dp)
                column = self%inverse(columns(p))
                if (self%marker(column) == i) then
                    self%info%duplicate_entries = &
                        self%info%duplicate_entries + 1_i8
                else
                    call touch(column)
                end if
                if (abs(values(p)) <= 0.0_dp) &
                    self%info%stored_zero_entries = &
                    self%info%stored_zero_entries + 1_i8
                self%work(column) = &
                    update_value(self%work(column), -values(p), 1.0_dp)
                if (.not. ieee_is_finite(self%work(column))) then
                    call fail(FORTSPARSE_INVALID_MATRIX, &
                              'ILUTP duplicate accumulation overflow')
                    return
                end if
            end do
            mean = update_value(0.0_dp, -scale, &
                                divide_value(mean, &
                                             real(row_ptr(i + 1) - row_ptr(i), dp)))
            call touch(i)
            do while (pending > 0)
                j = lower_pop(self%heap, pending)
                if (abs(self%work(j)) <= 0.0_dp) cycle
                multiplier = divide_value(self%work(j), self%diagonal(j))
                if (.not. ieee_is_finite(multiplier)) then
                    call fail(FORTSPARSE_SINGULAR, &
                              'nonfinite ILUTP elimination multiplier')
                    return
                end if
                if (abs(multiplier) <= drop_tolerance) then
                    self%work(j) = 0.0_dp
                    self%info%lower_threshold_drops = &
                        self%info%lower_threshold_drops + 1_i8
                    cycle
                end if
                self%work(j) = multiplier
                do p = self%up(j), self%up(j + 1) - 1
                    column = self%inverse(self%uc(p))
                    call touch(column)
                    self%work(column) = update_value(self%work(column), &
                                                     multiplier, self%uv(p))
                    if (.not. ieee_is_finite(self%work(column))) then
                        call fail(FORTSPARSE_SINGULAR, &
                                  'nonfinite ILUTP Schur update')
                        return
                    end if
                end do
            end do
            call select_largest(self%work, self%touched, count, i, .true., &
                                lower_fill, self%heap, kept, eligible)
            self%info%lower_fill_drops = &
                self%info%lower_fill_drops + int(eligible - kept, i8)
            do k = 1, kept
                nl = nl + 1
                self%lc(nl) = self%heap(k)
                self%lv(nl) = self%work(self%heap(k))
            end do
            self%lp(i + 1) = nl + 1
            threshold = update_value(0.0_dp, -drop_tolerance, mean)
            do k = 1, count
                column = self%touched(k)
                if (column <= i) cycle
                if (abs(self%work(column)) <= 0.0_dp) cycle
                if (abs(self%work(column)) <= threshold) then
                    self%work(column) = 0.0_dp
                    self%info%upper_threshold_drops = &
                        self%info%upper_threshold_drops + 1_i8
                end if
            end do
            call select_largest(self%work, self%touched, count, i, &
                                .false., &
                                upper_fill, self%heap, kept, eligible)
            self%info%upper_fill_drops = &
                self%info%upper_fill_drops + int(eligible - kept, i8)
            pivot = i
            value = abs(self%work(i))
            do k = 1, kept
                column = self%heap(k)
                if (update_value(0.0_dp, -abs(self%work(column)), &
                                 pivot_tolerance) &
                    <= abs(self%work(i))) cycle
                if (abs(self%work(column)) > value) then
                    pivot = column
                    value = abs(self%work(column))
                else if (.not. (abs(self%work(column)) < value)) then
                    pivot = min(pivot, column)
                end if
            end do
            if (pivot /= i) then
                diagonal_value = self%work(i)
                self%work(i) = self%work(pivot)
                self%work(pivot) = diagonal_value
                old_column = self%perm(i)
                self%perm(i) = self%perm(pivot)
                self%perm(pivot) = old_column
                self%inverse(self%perm(i)) = i
                self%inverse(self%perm(pivot)) = pivot
                self%info%column_swaps = self%info%column_swaps + 1
            end if
            do k = 1, kept
                nu = nu + 1
                self%uc(nu) = self%perm(self%heap(k))
                self%uv(nu) = self%work(self%heap(k))
            end do
            self%up(i + 1) = nu + 1
            diagonal_value = self%work(i)
            if (abs(diagonal_value) <= 0.0_dp) then
                diagonal_value = update_value(0.0_dp, -mean, &
                                              update_value(1.e-4_dp, &
                                                           -drop_tolerance, 1.0_dp))
                self%info%regularized_pivots = &
                    self%info%regularized_pivots + 1
            end if
            if (.not. ieee_is_finite(diagonal_value)) then
                call fail(FORTSPARSE_SINGULAR, &
                          'nonfinite ILUTP pivot replacement')
                return
            end if
            if (abs(diagonal_value) <= 0.0_dp) then
                call fail(FORTSPARSE_SINGULAR, &
                          'ILUTP pivot replacement underflow')
                return
            end if
            self%diagonal(i) = diagonal_value
        end do
        do p = 1, nu
            self%uc(p) = self%inverse(self%uc(p))
        end do
        deallocate (self%marker, self%touched, self%heap)
        self%info%lower_nnz = nl
        self%info%upper_nnz = nu
        self%info%ready = .true.
        call status_set(status, FORTSPARSE_OK, '')
    contains
        subroutine touch(index)
            integer, intent(in) :: index
            if (self%marker(index) == i) return
            self%marker(index) = i
            self%work(index) = 0.0_dp
            count = count + 1
            self%touched(count) = index
            if (index < i) call lower_push(self%heap, pending, index)
        end subroutine

        subroutine fail(code, message)
            integer, intent(in) :: code
            character(*), intent(in) :: message
            type(sparse_ilutp_stats_t) :: failed
            failed = self%info
            failed%peak_owned_bytes = max(failed%peak_owned_bytes, &
                                          self%memory_bytes())
            failed%failed_required_bytes = required
            call self%free()
            self%info = failed
            self%info%ready = .false.
            call status_set(status, code, message)
        end subroutine
    end subroutine setup

    subroutine apply(self, b, x, status)
        class(sparse_ilutp_t), intent(inout) :: self
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: x(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: i, p, n
        real(dp) :: value
        call status_set(status, FORTSPARSE_NOT_FACTORED, &
                        'ILUTP context is not ready')
        if (.not. self%info%ready) return
        n = self%info%n
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        'invalid ILUTP apply RHS/shape')
        if (size(b) /= n .or. size(x) /= n) return
        do i = 1, n
            if (.not. ieee_is_finite(b(i))) return
            self%work(i) = b(i)
        end do
        call status_set(status, FORTSPARSE_INTERNAL_ERROR, &
                        'nonfinite ILUTP apply arithmetic')
        do i = 1, n
            value = self%work(i)
            do p = self%lp(i), self%lp(i + 1) - 1
                value = update_value(value, self%lv(p), &
                                     self%work(self%lc(p)))
            end do
            if (.not. ieee_is_finite(value)) return
            self%work(i) = value
        end do
        do i = n, 1, -1
            value = self%work(i)
            do p = self%up(i), self%up(i + 1) - 1
                value = update_value(value, self%uv(p), &
                                     self%work(self%uc(p)))
            end do
            value = divide_value(value, self%diagonal(i))
            if (.not. ieee_is_finite(value)) return
            self%work(i) = value
        end do
        do i = 1, n
            x(self%perm(i)) = self%work(i)
        end do
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine

    function stats(self) result(info)
        class(sparse_ilutp_t), intent(in) :: self
        type(sparse_ilutp_stats_t) :: info
        info = self%info
        info%current_bytes = self%memory_bytes()
    end function

    integer(i8) function memory_bytes(self) result(bytes)
        class(sparse_ilutp_t), intent(in) :: self
        bytes = 0_i8
        if (allocated(self%lp)) bytes = bytes + ib*size(self%lp, kind=i8)
        if (allocated(self%up)) bytes = bytes + ib*size(self%up, kind=i8)
        if (allocated(self%lc)) bytes = bytes + ib*size(self%lc, kind=i8)
        if (allocated(self%uc)) bytes = bytes + ib*size(self%uc, kind=i8)
        if (allocated(self%perm)) bytes = &
            bytes + ib*size(self%perm, kind=i8)
        if (allocated(self%inverse)) bytes = &
            bytes + ib*size(self%inverse, kind=i8)
        if (allocated(self%marker)) bytes = &
            bytes + ib*size(self%marker, kind=i8)
        if (allocated(self%touched)) bytes = &
            bytes + ib*size(self%touched, kind=i8)
        if (allocated(self%heap)) bytes = &
            bytes + ib*size(self%heap, kind=i8)
        if (allocated(self%lv)) bytes = bytes + vb*size(self%lv, kind=i8)
        if (allocated(self%uv)) bytes = bytes + vb*size(self%uv, kind=i8)
        if (allocated(self%diagonal)) bytes = &
            bytes + vb*size(self%diagonal, kind=i8)
        if (allocated(self%work)) bytes = &
            bytes + vb*size(self%work, kind=i8)
    end function

    subroutine release_context(self)
        class(sparse_ilutp_t), intent(inout) :: self
        if (allocated(self%lp)) deallocate (self%lp)
        if (allocated(self%up)) deallocate (self%up)
        if (allocated(self%lc)) deallocate (self%lc)
        if (allocated(self%uc)) deallocate (self%uc)
        if (allocated(self%perm)) deallocate (self%perm)
        if (allocated(self%inverse)) deallocate (self%inverse)
        if (allocated(self%marker)) deallocate (self%marker)
        if (allocated(self%touched)) deallocate (self%touched)
        if (allocated(self%heap)) deallocate (self%heap)
        if (allocated(self%lv)) deallocate (self%lv)
        if (allocated(self%uv)) deallocate (self%uv)
        if (allocated(self%diagonal)) deallocate (self%diagonal)
        if (allocated(self%work)) deallocate (self%work)
        self%info = sparse_ilutp_stats_t()
    end subroutine
end module fortsparse_ilutp
