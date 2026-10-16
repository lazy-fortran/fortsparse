module fortsparse_memory
    ! Payload accounting only: allocator/runtime overhead and RSS are separate.
    use fortsparse_kinds, only: i8
    use fortsparse_csc, only: csc_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortsparse_status, only: fortsparse_status_t, status_set, &
        FORTSPARSE_OK, FORTSPARSE_INVALID_MATRIX, FORTSPARSE_INDEX_OVERFLOW, &
        FORTSPARSE_RESOURCE_LIMIT, FORTSPARSE_OUT_OF_MEMORY
    implicit none
    private
    public :: memory_budget_t, memory_checked_add, memory_checked_mul
    public :: sparse_memory_plan_t, plan_natural_spd_memory
    public :: MEMORY_PLAN_EXACT_LAYOUT, MEMORY_PLAN_BOUND, MEMORY_PLAN_ESTIMATE
    integer, parameter :: MEMORY_PLAN_EXACT_LAYOUT = 1
    integer, parameter :: MEMORY_PLAN_BOUND = 2, MEMORY_PLAN_ESTIMATE = 3

    type :: sparse_memory_plan_t
        ! Disjoint payload classes; peak includes factor + persistent + workspace.
        ! Caller A/B/X, object descriptors, allocator and runtime/RSS are excluded.
        integer(i8) :: factor_host_bytes = 0_i8, factor_device_bytes = 0_i8
        integer(i8) :: persistent_host_bytes = 0_i8, persistent_device_bytes = 0_i8
        integer(i8) :: workspace_host_bytes = 0_i8, workspace_device_bytes = 0_i8
        integer(i8) :: peak_host_bytes = 0_i8, peak_device_bytes = 0_i8
        integer(i8) :: symbolic_nnz = 0_i8
        integer :: max_rhs_batch = 0, prediction_quality = MEMORY_PLAN_ESTIMATE
        logical :: fits = .false.
    end type sparse_memory_plan_t

    type :: memory_budget_t
        integer(i8) :: limit = 0_i8 ! Zero means unlimited.
        integer(i8) :: current = 0_i8, peak = 0_i8
        integer(i8) :: failed_required = 0_i8, failed_available = 0_i8
    contains
        procedure :: reserve => budget_reserve
        procedure :: release => budget_release
    end type memory_budget_t

contains

    subroutine plan_natural_spd_memory(A, plan, status, budget, nrhs)
        ! Exact natural-order simplicial reference layout, NOT blocked/padded
        ! panels or an SPD proof. Structural fill analysis allocates only 3*n
        ! integers, O(n); traversals may visit every filled graph entry.
        type(csc_t), intent(in) :: A
        type(sparse_memory_plan_t), intent(out) :: plan
        type(fortsparse_status_t), intent(out) :: status
        type(memory_budget_t), intent(inout), optional :: budget
        integer, intent(in), optional :: nrhs
        type(memory_budget_t) :: tracker
        type(fortsparse_status_t) :: ignored
        integer, allocatable :: parent(:), count(:), mark(:)
        integer :: n, rhs, j, p, i, ierr
        integer(i8) :: scratch, nz, part, total

        plan = sparse_memory_plan_t()
        rhs = 1
        if (present(nrhs)) rhs = nrhs
        if (rhs < 0) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'memory plan: negative RHS count')
            return
        end if
        call validate_spd_pattern(A, status)
        if (status%code /= FORTSPARSE_OK) return
        if (present(budget)) tracker = budget
        call tracker%reserve(0_i8, status)
        if (status%code /= FORTSPARSE_OK) return
        n = A%ncol
        scratch = 0_i8
        call memory_checked_mul(int(n,i8), 12_i8, scratch, status)
        if (status%code /= FORTSPARSE_OK) return
        call tracker%reserve(scratch, status, 'SPD planning graph scratch')
        if (status%code /= FORTSPARSE_OK) then
            if (present(budget)) budget = tracker
            return
        end if
        allocate (parent(n), count(n), mark(n), stat=ierr)
        if (ierr /= 0) then
            if (allocated(parent)) deallocate (parent)
            if (allocated(count)) deallocate (count)
            if (allocated(mark)) deallocate (mark)
            call tracker%release(scratch, ignored)
            if (present(budget)) budget = tracker
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, &
                'memory plan: graph scratch allocation failed')
            return
        end if
        parent = 0
        count = 0
        mark = 0
        do j = 1, n
            do p = A%col_ptr(j), A%col_ptr(j+1)-1
                i = A%row_idx(p)
                if (i >= j) exit
                do while (i /= 0)
                    if (mark(i) == j) exit
                    mark(i) = j
                    count(i) = count(i)+1
                    if (parent(i) == 0) parent(i) = j
                    i = parent(i)
                    if (i == j) exit
                end do
            end do
        end do
        nz = int(n,i8)
        do j = 1, n
            nz = nz+int(count(j),i8)
        end do
        deallocate (parent, count, mark)
        call tracker%release(scratch, ignored)
        if (present(budget)) budget = tracker
        if (nz >= int(huge(n),i8)) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                'memory plan: factor exceeds supported index range')
            return
        end if
        plan%symbolic_nnz = nz
        plan%prediction_quality = MEMORY_PLAN_EXACT_LAYOUT
        ! lv(L), li(L), ap(n+1), ai(nnz), lp(n+1), five integer n
        ! arrays, and real work(n). Reference block solves add no workspace.
        call memory_checked_mul(nz, 8_i8, plan%factor_host_bytes, status)
        if (status%code /= FORTSPARSE_OK) return
        total = 0_i8
        call memory_checked_mul(nz, 4_i8, total, status)
        if (status%code /= FORTSPARSE_OK) return
        part = 0_i8
        call memory_checked_mul(int(A%nnz,i8), 4_i8, part, status)
        if (status%code /= FORTSPARSE_OK) return
        call memory_checked_add(total, part, plan%persistent_host_bytes, status)
        if (status%code /= FORTSPARSE_OK) return
        call memory_checked_mul(int(n,i8), 36_i8, part, status)
        if (status%code /= FORTSPARSE_OK) return
        call memory_checked_add(plan%persistent_host_bytes, part, total, status)
        if (status%code /= FORTSPARSE_OK) return
        call memory_checked_add(total, 8_i8, plan%persistent_host_bytes, status)
        if (status%code /= FORTSPARSE_OK) return
        call memory_checked_add(plan%factor_host_bytes, &
            plan%persistent_host_bytes, plan%peak_host_bytes, status)
        if (status%code /= FORTSPARSE_OK) return
        ! Reserve only on a copy: predicted factor bytes are not allocated here
        ! and must not inflate the caller's measured live/peak ledger.
        call tracker%reserve(plan%peak_host_bytes, status, 'SPD reference plan')
        if (status%code /= FORTSPARSE_OK) then
            if (present(budget)) then
                budget%failed_required = tracker%failed_required
                budget%failed_available = tracker%failed_available
            end if
            return
        end if
        plan%fits = .true.
        plan%max_rhs_batch = rhs
    end subroutine plan_natural_spd_memory

    subroutine validate_spd_pattern(A, status)
        type(csc_t), intent(in) :: A
        type(fortsparse_status_t), intent(out) :: status
        integer :: n, j, p, row, previous, lo, hi, mid
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            'memory plan: require finite full symmetric canonical square CSC')
        if (A%ncol < 0 .or. A%nrow /= A%ncol .or. A%nnz < 0) return
        n = A%ncol
        if (n == huge(n) .or. A%nnz == huge(n)) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                'memory plan: CSC endpoints exceed supported index range')
            return
        end if
        if (.not. allocated(A%col_ptr)) return
        if (.not. allocated(A%row_idx)) return
        if (.not. allocated(A%val)) return
        if (size(A%col_ptr) /= n+1 .or. size(A%row_idx) /= A%nnz) return
        if (size(A%val) /= A%nnz) return
        if (A%col_ptr(1) /= 1 .or. A%col_ptr(n+1) /= A%nnz+1) return
        do p = 1, A%nnz
            if (.not. ieee_is_finite(A%val(p))) return
        end do
        do j = 1, n
            if (A%col_ptr(j+1) < A%col_ptr(j)) return
            if (A%col_ptr(j+1) > A%nnz+1) return
        end do
        do j = 1, n
            previous = 0
            do p = A%col_ptr(j), A%col_ptr(j+1)-1
                row = A%row_idx(p)
                if (row <= previous .or. row > n) return
                previous = row
            end do
        end do
        do j = 1, n
            do p = A%col_ptr(j), A%col_ptr(j+1)-1
                row = A%row_idx(p)
                lo = A%col_ptr(row)
                hi = A%col_ptr(row+1)-1
                do while (lo <= hi)
                    mid = lo+(hi-lo)/2
                    if (A%row_idx(mid) < j) then
                        lo = mid+1
                    else if (A%row_idx(mid) > j) then
                        hi = mid-1
                    else
                        exit
                    end if
                end do
                if (lo > hi) return
                ! Exact symmetry is required; both values are already finite.
                if (A%val(p) < A%val(mid) .or. A%val(p) > A%val(mid)) return
            end do
        end do
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine validate_spd_pattern

    pure subroutine memory_checked_add(a, b, value, status)
        integer(i8), intent(in) :: a, b
        integer(i8), intent(inout) :: value
        type(fortsparse_status_t), intent(out) :: status
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            'memory arithmetic requires nonnegative byte counts')
        if (a < 0_i8 .or. b < 0_i8) return
        if (a > huge(a)-b) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                'memory byte addition overflow')
            return
        end if
        value = a+b
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine memory_checked_add

    pure subroutine memory_checked_mul(a, b, value, status)
        integer(i8), intent(in) :: a, b
        integer(i8), intent(inout) :: value
        type(fortsparse_status_t), intent(out) :: status
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            'memory arithmetic requires nonnegative byte counts')
        if (a < 0_i8 .or. b < 0_i8) return
        if (a > 0_i8) then
            if (b > huge(a)/a) then
                call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                    'memory byte multiplication overflow')
                return
            end if
        end if
        value = a*b
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine memory_checked_mul

    pure logical function budget_valid(self)
        class(memory_budget_t), intent(in) :: self
        budget_valid = .false.
        if (self%limit < 0_i8 .or. self%current < 0_i8) return
        if (self%peak < self%current) return
        if (self%limit > 0_i8) then
            if (self%current > self%limit) return
        end if
        budget_valid = .true.
    end function budget_valid

    pure subroutine budget_reserve(self, bytes, status, phase)
        class(memory_budget_t), intent(inout) :: self
        integer(i8), intent(in) :: bytes
        type(fortsparse_status_t), intent(out) :: status
        character(*), intent(in), optional :: phase
        integer(i8) :: required
        character(120) :: message

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            'invalid memory budget or reservation')
        if (.not. budget_valid(self)) return
        required = 0_i8
        call memory_checked_add(self%current, bytes, required, status)
        if (status%code /= FORTSPARSE_OK) return
        if (self%limit > 0_i8) then
            if (required > self%limit) then
                ! Live/peak accounting is transactional; these two diagnostics
                ! describe the last refused increment, not a live allocation.
                self%failed_required = bytes
                self%failed_available = self%limit-self%current
                message = 'memory cap exceeded'
                if (present(phase)) then
                    message = 'memory cap exceeded: '
                    message(22:) = phase
                end if
                call status_set(status, FORTSPARSE_RESOURCE_LIMIT, message)
                return
            end if
        end if
        self%current = required
        self%peak = max(self%peak, required)
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine budget_reserve

    pure subroutine budget_release(self, bytes, status)
        class(memory_budget_t), intent(inout) :: self
        integer(i8), intent(in) :: bytes
        type(fortsparse_status_t), intent(out) :: status
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            'invalid memory budget or release')
        if (.not. budget_valid(self)) return
        if (bytes < 0_i8 .or. bytes > self%current) return
        self%current = self%current-bytes
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine budget_release
end module fortsparse_memory
