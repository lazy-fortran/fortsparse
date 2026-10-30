module fortsparse_native_ldlt_z
    ! Serial sparse complex Bunch--Kaufman reference: P^* A P = L D L^*.
    ! * is adjoint for Hermitian input, transpose for transpose-symmetric input.
    ! Full canonical symmetric CSC is required. Schur rows and L columns use
    ! original vertex labels: symmetric pivoting changes maps, never n^2 data.
    ! Numeric factorization grows sparse rows under an explicit payload cap.
    ! This reference is not a production performance/parallelism admission.
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csc, only: csc_z_t, csc_is_valid
    use fortsparse_status, only: fortsparse_status_t, status_set, status_ok, &
                        FORTSPARSE_OK, FORTSPARSE_SINGULAR, FORTSPARSE_INVALID_MATRIX, &
                                 FORTSPARSE_NOT_FACTORED, FORTSPARSE_OUT_OF_MEMORY, &
                                 FORTSPARSE_INDEX_OVERFLOW, FORTSPARSE_RESOURCE_LIMIT, &
                                 FORTSPARSE_INTERNAL_ERROR
    use fortsparse_lu_leaves, only: lu_z_divide, lu_z_update, lu_z_norm, lu_z_finite
    use fortsparse_ldlt_z_leaves, only: ldlt_z_det, ldlt_z_numerators, ldlt_z_update
    use fortsparse_inertia, only: inertia_count_t, inertia_add_pivot, &
                                  inertia_status_of, inertia_tolerance, &
                                  FORTSPARSE_INERTIA_UNAVAILABLE
    implicit none
    private
    public :: native_ldlt_z_t
    integer(i8), parameter :: ib = storage_size(0)/8
    integer(i8), parameter :: rb = storage_size(cmplx(0.0_dp, 0.0_dp, dp))/8
    integer(i8), parameter :: rsb = storage_size(0.0_dp)/8
    real(dp), parameter :: alpha = (1.0_dp + sqrt(17.0_dp))/8.0_dp
    ! Threshold u of the ordered pivot tests: every multiplier is bounded by 1/u.
    real(dp), parameter :: delay_threshold = 0.1_dp

    type :: sparse_row_t
        integer :: used = 0
        integer, allocatable :: index(:)
        complex(dp), allocatable :: value(:)
    end type sparse_row_t

    type :: native_ldlt_z_t
        private
        integer :: n = 0
        logical :: ready = .false., low_memory = .false.
        logical :: is_hermitian = .true.
        integer(i8) :: limit = 0_i8, bytes = 0_i8, peak = 0_i8
        integer(i8) :: required = 0_i8
        ! Hermitian inertia of the last completed elimination (real pivots).
        logical :: inertia_ready = .false.
        type(inertia_count_t) :: pivots
        ! Per-vertex max |entry| of A and every Schur entry stored in its row.
        real(dp), allocatable :: row_scale(:)
        integer, allocatable :: ap(:), ai(:), perm(:), pos(:), block(:), nbr(:)
        integer, allocatable :: initial_order(:)
        ! Ordered factors only: elimination step at which a vertex was last delayed.
        integer, allocatable :: delayed(:)
        complex(dp), allocatable :: diagonal(:), offdiagonal(:), u(:), v(:), work(:)
        type(sparse_row_t), allocatable :: schur(:), lower(:)
    contains
        procedure :: factor => ldlt_factor
        procedure :: refactor => ldlt_refactor
        procedure :: solve => ldlt_solve
        procedure :: solve_inplace => ldlt_solve_inplace
        procedure :: solve_block => ldlt_solve_block
        procedure :: free => ldlt_free
        procedure :: retained_bytes => ldlt_bytes
        procedure :: peak_bytes => ldlt_peak
        procedure :: required_bytes => ldlt_required
        procedure :: available_bytes => ldlt_available
        procedure :: factor_nnz => ldlt_nnz
        procedure :: permutation => ldlt_permutation
        procedure :: copy_column_order => ldlt_copy_column_order
        procedure :: inertia => ldlt_inertia
        final :: ldlt_finalize
    end type native_ldlt_z_t
contains
    subroutine ldlt_factor(self, a, status, hermitian, host_limit_bytes, low_memory, &
            column_permutation)
        class(native_ldlt_z_t), intent(inout) :: self
        type(csc_z_t), intent(in) :: a
        type(fortsparse_status_t), intent(out) :: status
        integer(i8), intent(in), optional :: host_limit_bytes
        logical, intent(in), optional :: low_memory, hermitian
        integer, intent(in), optional :: column_permutation(:)
        call self%free()
        if (present(hermitian)) self%is_hermitian = hermitian
        if (present(host_limit_bytes)) self%limit = host_limit_bytes
        if (present(low_memory)) self%low_memory = low_memory
        call validate(a, self%is_hermitian, status)
        if (.not. status_ok(status)) return
        if (self%limit < 0_i8) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, 'LDLT: negative cap')
            return
        end if
        call initialize(self, a, status, column_permutation)
        if (.not. status_ok(status)) return
        call numeric(self, status)
    end subroutine ldlt_factor

    subroutine ldlt_refactor(self, a, status, hermitian, host_limit_bytes, low_memory, &
            column_permutation)
        class(native_ldlt_z_t), intent(inout) :: self
        type(csc_z_t), intent(in) :: a
        type(fortsparse_status_t), intent(out) :: status
        integer(i8), intent(in), optional :: host_limit_bytes
        logical, intent(in), optional :: low_memory, hermitian
        integer, intent(in), optional :: column_permutation(:)
        integer :: i
        if (present(hermitian)) then
            if (hermitian .neqv. self%is_hermitian) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                                'complex LDL: refactor cannot change symmetry class')
                return
            end if
        end if
        call validate(a, self%is_hermitian, status)
        if (.not. status_ok(status)) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'LDLT: pattern changed')
        if (.not. allocated(self%ap)) return
        if (.not. allocated(self%ai)) return
        if (a%ncol /= self%n) return
        if (size(a%row_idx) /= size(self%ai)) return
        do i = 1, self%n + 1
            if (a%col_ptr(i) /= self%ap(i)) return
        end do
        do i = 1, a%nnz
            if (a%row_idx(i) /= self%ai(i)) return
        end do
        ! The final pivot permutation is value-dependent; reuse only the same
        ! immutable initial order, then choose numerical pivots afresh.
        if (present(column_permutation)) then
            if (.not. allocated(self%initial_order)) return
            if (size(column_permutation) /= self%n) return
            if (any(column_permutation /= self%initial_order)) return
        else
            if (allocated(self%initial_order)) return
        end if
        if (present(host_limit_bytes)) then
            if (host_limit_bytes < 0_i8) return
            self%limit = host_limit_bytes
        end if
        if (present(low_memory)) self%low_memory = low_memory
        ! Destructive numeric refactor: drop old numeric factors before building
        ! new ones, retaining only the immutable pattern. No old/new overlap.
        call clear_numeric(self)
        self%peak = self%bytes
        self%required = self%bytes
        call initialize(self, a, status)
        if (.not. status_ok(status)) return
        call numeric(self, status)
    end subroutine ldlt_refactor

    subroutine validate(a, hermitian, status)
        type(csc_z_t), intent(in) :: a
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(in) :: hermitian
        complex(dp) :: opposite
        integer :: j, p, i, previous
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
               'complex LDL: full canonical square CSC with declared symmetry required')
        if (a%ncol <= 0 .or. a%ncol == huge(0)) return
        if (a%nrow /= a%ncol) return
        if (a%nnz < 0) return
        if (a%nnz == huge(0)) then
          call status_set(status, FORTSPARSE_INDEX_OVERFLOW, 'LDLT: CSC index overflow')
            return
        end if
        if (.not. allocated(a%col_ptr)) return
        if (.not. allocated(a%row_idx)) return
        if (.not. allocated(a%val)) return
        if (size(a%val) /= a%nnz) return
        if (.not. csc_is_valid(a)) return
        do j = 1, a%ncol
            previous = 0
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                i = a%row_idx(p)
                if (i <= previous) return
                previous = i
                if (.not. lu_z_finite(a%val(p))) return
                if (.not. norm_in_range(a%val(p))) then
                    call arithmetic_error(status)
                    return
                end if
                opposite = mate(csc_get(a, j, i), hermitian)
                if (real(a%val(p), dp) < real(opposite, dp)) return
                if (real(a%val(p), dp) > real(opposite, dp)) return
                if (aimag(a%val(p)) < aimag(opposite)) return
                if (aimag(a%val(p)) > aimag(opposite)) return
                if (hermitian .and. i == j) then
                    if (abs(aimag(a%val(p))) > 0.0_dp) return
                end if
            end do
        end do
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine validate

    pure complex(dp) function csc_get(a, row, column) result(value)
        type(csc_z_t), intent(in) :: a
        integer, intent(in) :: row, column
        integer :: lo, hi, mid
        lo = a%col_ptr(column)
        hi = a%col_ptr(column + 1) - 1
        value = 0.0_dp
        do while (lo <= hi)
            mid = lo + (hi - lo)/2
            if (a%row_idx(mid) == row) then
                value = a%val(mid)
                return
            else if (a%row_idx(mid) < row) then
                lo = mid + 1
            else
                hi = mid - 1
            end if
        end do
    end function csc_get

    subroutine initialize(self, a, status, column_permutation)
        class(native_ldlt_z_t), intent(inout) :: self
        type(csc_z_t), intent(in) :: a
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: column_permutation(:)
        integer :: n, ierr, j, count, p
        integer(i8) :: request
        n = a%ncol
        self%n = n
        if (.not. allocated(self%ap)) then
            request = ib*(int(n, i8) + 1_i8 + int(a%nnz, i8))
            call admit(self, request, status)
            if (.not. status_ok(status)) return
            allocate (self%ap(n + 1), self%ai(a%nnz), stat=ierr)
            if (ierr /= 0) then
                call allocation_error(status)
                call self%free()
                return
            end if
            self%ap = a%col_ptr
            self%ai = a%row_idx
        end if
        ! Includes descriptor-array payload, all scratch, permutation and D.
        request = int(n, i8)*(2_i8*int(storage_size(sparse_row_t())/8, i8) &
                              + 4_i8*ib + 5_i8*rb + rsb)
        call admit(self, request, status)
        if (.not. status_ok(status)) return
        allocate (self%schur(n), self%lower(n), self%perm(n), self%pos(n), &
                  self%block(n), self%nbr(n), self%diagonal(n), self%offdiagonal(n), &
                  self%u(n), self%v(n), self%work(n), self%row_scale(n), &
                  stat=ierr)
        if (ierr /= 0) then
            call allocation_error(status)
            call self%free()
            return
        end if
        if (present(column_permutation)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'LDLT: invalid initial order')
            if (size(column_permutation) /= n) return
            self%pos = 0
            do j = 1, n
                p = column_permutation(j)
                if (p < 1 .or. p > n) return
                if (self%pos(p) /= 0) return
                self%pos(p) = j
            end do
            if (.not. allocated(self%initial_order)) then
                call admit(self, ib*int(n, i8), status)
                if (.not. status_ok(status)) return
                allocate (self%initial_order(n), stat=ierr)
                if (ierr /= 0) then
                    call allocation_error(status)
                    call self%free()
                    return
                end if
                self%initial_order = column_permutation
            end if
        end if
        if (allocated(self%initial_order)) then
            call admit(self, ib*int(n, i8), status)
            if (.not. status_ok(status)) return
            allocate (self%delayed(n), stat=ierr)
            if (ierr /= 0) then
                call allocation_error(status)
                call self%free()
                return
            end if
            self%delayed = 0
        end if
        do j = 1, n
            self%perm(j) = j
            if (allocated(self%initial_order)) self%perm(j) = self%initial_order(j)
            self%pos(self%perm(j)) = j
        end do
        do j = 1, n
            count = a%col_ptr(j + 1) - a%col_ptr(j)
            call reserve_row(self%schur(j), count, .true., self%limit, &
                             self%bytes, self%peak, self%required, status)
            if (.not. status_ok(status)) return
            self%schur(j)%used = count
            do p = 1, count
                self%schur(j)%index(p) = a%row_idx(a%col_ptr(j) + p - 1)
                self%schur(j)%value(p) = mate(a%val(a%col_ptr(j) + p - 1), &
                                              self%is_hermitian)
            end do
        end do
        self%block = 0
        self%diagonal = 0.0_dp
        self%offdiagonal = 0.0_dp
        ! Componentwise max avoids complex-modulus overflow near huge().
        do j = 1, n
            self%row_scale(j) = 0.0_dp
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                self%row_scale(j) = max(self%row_scale(j), &
                    abs(real(a%val(p), dp)), abs(aimag(a%val(p))))
            end do
        end do
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine initialize

    subroutine numeric(self, status)
        class(native_ldlt_z_t), intent(inout) :: self
        type(fortsparse_status_t), intent(out) :: status
        integer :: k, step, pivot, first, second, m, i, j, p, q
        complex(dp) :: a, b, c, value, l1, l2
        self%pivots = inertia_count_t()
        k = 1
        do while (k <= self%n)
            if (allocated(self%delayed)) then
                call ordered_pivot(self, k, pivot, step, status)
                if (.not. status_ok(status)) return
                call rotate(self, pivot, k + step - 1)
            else
                call choose_pivot(self, k, pivot, step, status)
                if (.not. status_ok(status)) return
                call interchange(self, k + step - 1, pivot)
            end if
            first = self%perm(k)
            second = first
            a = row_get(self%schur(first), first)
            b = 0.0_dp
            c = 0.0_dp
            if (step == 2) then
                second = self%perm(k + 1)
                b = row_get(self%schur(first), second)
                c = row_get(self%schur(second), second)
            end if
            self%diagonal(k) = a
            self%offdiagonal(k) = b
            self%block(k) = step
            if (step == 2) self%diagonal(k + 1) = c
            call add_pivot_inertia(self, step, a, b, c, first, second)
            if (step == 2) then
                call pair_solve(a, b, c, (0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), &
                                self%is_hermitian, l1, l2, status)
                if (.not. status_ok(status)) return
            end if
            ! Gather the sparse union of both pivot rows into bounded scratch.
            m = 0
            p = 1
            q = 1
            if (step == 1) q = self%schur(second)%used + 1
            do
                i = self%n + 1
                j = self%n + 1
                if (p <= self%schur(first)%used) i = self%schur(first)%index(p)
                if (q <= self%schur(second)%used) j = self%schur(second)%index(q)
                i = min(i, j)
                if (i > self%n) exit
                value = 0.0_dp
                l2 = 0.0_dp
                if (p <= self%schur(first)%used) then
                    if (self%schur(first)%index(p) == i) then
                        value = mate(self%schur(first)%value(p), &
                                     self%is_hermitian)
                        p = p + 1
                    end if
                end if
                if (q <= self%schur(second)%used) then
                    if (self%schur(second)%index(q) == i) then
                        l2 = mate(self%schur(second)%value(q), &
                                  self%is_hermitian)
                        q = q + 1
                    end if
                end if
                if (self%pos(i) < k + step) cycle
                if (abs(value) <= 0.0_dp .and. abs(l2) <= 0.0_dp) cycle
                m = m + 1
                self%nbr(m) = i
                self%u(m) = value
                self%v(m) = l2
            end do
            do j = k, k + step - 1
                call reserve_row(self%lower(j), m, .true., self%limit, &
                                 self%bytes, self%peak, self%required, status)
                if (.not. status_ok(status)) return
                self%lower(j)%used = m
            end do
            do p = 1, m
                l2 = 0.0_dp
                if (step == 1) then
                    l1 = lu_z_divide(self%u(p), a)
                else
                    call pair_solve(a, b, c, self%u(p), self%v(p), &
                                    self%is_hermitian, l1, l2, status)
                    if (.not. status_ok(status)) return
                end if
                if (.not. lu_z_finite(l1) .or. .not. lu_z_finite(l2)) then
                    call arithmetic_error(status)
                    return
                end if
                self%lower(k)%index(p) = self%nbr(p)
                self%lower(k)%value(p) = l1
                if (step == 2) then
                    self%lower(k + 1)%index(p) = self%nbr(p)
                    self%lower(k + 1)%value(p) = l2
                end if
            end do
            do p = 1, m
                i = self%nbr(p)
                l1 = self%lower(k)%value(p)
                l2 = 0.0_dp
                if (step == 2) l2 = self%lower(k + 1)%value(p)
                do q = p, m
                    j = self%nbr(q)
                    value = ldlt_z_update(row_get(self%schur(i), j), &
                                          l1, mate(self%u(q), self%is_hermitian), &
                                          l2, mate(self%v(q), self%is_hermitian))
                    if (self%is_hermitian .and. i == j) &
                        value = cmplx(real(value, dp), 0.0_dp, dp)
                    if (.not. lu_z_finite(value)) then
                        call arithmetic_error(status)
                        return
                    end if
                    if (.not. norm_in_range(value)) then
                        call arithmetic_error(status)
                        return
                    end if
                    self%row_scale(i) = max(self%row_scale(i), &
                        abs(real(value, dp)), abs(aimag(value)))
                    self%row_scale(j) = max(self%row_scale(j), &
                        abs(real(value, dp)), abs(aimag(value)))
                    call row_put(self%schur(i), j, value, self%low_memory, &
                               self%limit, self%bytes, self%peak, self%required, status)
                    if (.not. status_ok(status)) return
                    if (i /= j) then
                        call row_put(self%schur(j), i, mate(value, self%is_hermitian), &
                                     self%low_memory, &
                               self%limit, self%bytes, self%peak, self%required, status)
                        if (.not. status_ok(status)) return
                    end if
                end do
                call row_erase(self%schur(i), first)
                if (step == 2) call row_erase(self%schur(i), second)
            end do
            call release_row(self%schur(first), self%bytes)
            if (step == 2) call release_row(self%schur(second), self%bytes)
            k = k + step
        end do
        call release_scratch(self)
        self%inertia_ready = self%is_hermitian
        if (self%pivots%zero > 0_i8) then
            ! Zero pivot columns were eliminated as exact zero 1x1 pivots so
            ! that Hermitian inertia is complete; solves are refused.
            call status_set(status, FORTSPARSE_SINGULAR, &
                            'LDLT: zero numerical pivot column; inertia only')
            return
        end if
        self%ready = .true.
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine numeric

    subroutine add_pivot_inertia(self, step, a, b, c, first, second)
        ! Hermitian pivots are real with a local rounding threshold from the
        ! pivot rows. Complex symmetric pivots carry no inertia; only exact
        ! zero pivots are counted there so that SINGULAR stays explicit.
        class(native_ldlt_z_t), intent(inout) :: self
        integer, intent(in) :: step, first, second
        complex(dp), intent(in) :: a, b, c
        if (.not. self%is_hermitian) then
            if (step == 1 .and. a == (0.0_dp, 0.0_dp)) &
                self%pivots%zero = self%pivots%zero + 1_i8
            return
        end if
        call inertia_add_pivot(self%pivots, step, real(a, dp), &
            hypot(real(b, dp), aimag(b)), real(c, dp), &
            inertia_tolerance(self%n, max(self%row_scale(first), &
                                          self%row_scale(second))))
    end subroutine add_pivot_inertia

    ! Hermitian inertia of the last completed factorization; UNAVAILABLE for
    ! complex symmetric input or when no elimination completed.
    subroutine ldlt_inertia(self, count, state, tolerance)
        class(native_ldlt_z_t), intent(in) :: self
        type(inertia_count_t), intent(out) :: count
        integer, intent(out) :: state
        real(dp), intent(out), optional :: tolerance
        state = FORTSPARSE_INERTIA_UNAVAILABLE
        if (present(tolerance)) tolerance = -1.0_dp
        if (.not. self%inertia_ready) return
        count = self%pivots
        state = inertia_status_of(count)
        if (present(tolerance)) tolerance = count%tolerance
    end subroutine ldlt_inertia

    subroutine choose_pivot(self, k, pivot, step, status)
        class(native_ldlt_z_t), intent(in) :: self
        integer, intent(in) :: k
        integer, intent(out) :: pivot, step
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: absakk, colmax, rowmax, magnitude
        integer :: first, largest, i, p
        first = self%perm(k)
        absakk = pivot_norm(row_get(self%schur(first), first))
        largest = first
        colmax = 0.0_dp
        do p = 1, self%schur(first)%used
            i = self%schur(first)%index(p)
            if (self%pos(i) <= k) cycle
            magnitude = pivot_norm(self%schur(first)%value(p))
            if (magnitude > colmax) then
                colmax = magnitude
                largest = i
            end if
        end do
        pivot = k
        step = 1
        call status_set(status, FORTSPARSE_OK, '')
        ! A zero column is an exact zero 1x1 pivot; numeric completes the
        ! elimination for inertia and then reports SINGULAR.
        if (max(absakk, colmax) <= 0.0_dp) return
        if (absakk >= alpha*colmax) return
        rowmax = 0.0_dp
        do p = 1, self%schur(largest)%used
            i = self%schur(largest)%index(p)
            if (self%pos(i) < k .or. i == largest) cycle
            rowmax = max(rowmax, pivot_norm(self%schur(largest)%value(p)))
        end do
        ! rowmax >= colmax > 0: quotient is safe and never exceeds one.
        if (absakk > 0.0_dp) then
            if ((absakk/colmax) >= &
                alpha*(colmax/rowmax)) return
        end if
        pivot = self%pos(largest)
        if (pivot_norm(row_get(self%schur(largest), largest)) >= alpha*rowmax) return
        step = 2
    end subroutine choose_pivot

    ! Pivoting restricted to the fill-reducing order. A 1x1 pivot in order or
    ! a 2x2 pivot with the first remaining neighbour (the elimination-tree
    ! parent) must pass the threshold test; otherwise the vertex is delayed to
    ! just after that parent, as a multifrontal solver passes a failed pivot
    ! to the parent front. A vertex delayed again before any elimination takes
    ! the unrestricted Bunch--Kaufman pivot, which always exists.
    subroutine ordered_pivot(self, k, pivot, step, status)
        class(native_ldlt_z_t), intent(inout) :: self
        integer, intent(in) :: k
        integer, intent(out) :: pivot, step
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: absakk, colmax, a, b, c, det, first_max, parent_max, block_norm
        complex(dp) :: an, bn, cn
        integer :: first, parent, i, p, block_exponent
        do
            first = self%perm(k)
            absakk = pivot_norm(row_get(self%schur(first), first))
            colmax = 0.0_dp
            parent = 0
            do p = 1, self%schur(first)%used
                i = self%schur(first)%index(p)
                if (self%pos(i) <= k) cycle
                colmax = max(colmax, pivot_norm(self%schur(first)%value(p)))
                if (parent == 0) then
                    parent = i
                else if (self%pos(i) < self%pos(parent)) then
                    parent = i
                end if
            end do
            pivot = k
            step = 1
            call status_set(status, FORTSPARSE_OK, '')
            if (max(absakk, colmax) <= 0.0_dp) return
            if (absakk >= delay_threshold*colmax) return
            first_max = 0.0_dp
            do p = 1, self%schur(first)%used
                i = self%schur(first)%index(p)
                if (self%pos(i) <= k .or. i == parent) cycle
                first_max = max(first_max, pivot_norm(self%schur(first)%value(p)))
            end do
            parent_max = 0.0_dp
            do p = 1, self%schur(parent)%used
                i = self%schur(parent)%index(p)
                if (self%pos(i) <= k .or. i == parent) cycle
                parent_max = max(parent_max, pivot_norm(self%schur(parent)%value(p)))
            end do
            an = row_get(self%schur(first), first)
            bn = row_get(self%schur(first), parent)
            cn = row_get(self%schur(parent), parent)
            det = 0.0_dp
            block_exponent = 0
            block_norm = max(lu_z_norm(an), lu_z_norm(bn), lu_z_norm(cn))
            if (block_norm > 0.0_dp) then
                block_exponent = exponent(block_norm)
                an = scale_complex(an, -block_exponent)
                bn = scale_complex(bn, -block_exponent)
                cn = scale_complex(cn, -block_exponent)
                det = abs(ldlt_z_det(an, bn, mate(bn, self%is_hermitian), cn))
            end if
            a = pivot_norm(an)
            b = pivot_norm(bn)
            c = pivot_norm(cn)
            ! |P^{-1}| (first_max, parent_max)^T <= 1/u with P scaled by 2^-e.
            if (det > 0.0_dp) then
                if (delay_threshold*scale(c*first_max + b*parent_max, -block_exponent) &
                    <= det .and. &
                    delay_threshold*scale(b*first_max + a*parent_max, -block_exponent) &
                    <= det) then
                    pivot = self%pos(parent)
                    step = 2
                    return
                end if
            end if
            if (self%delayed(first) == k) then
                call choose_pivot(self, k, pivot, step, status)
                return
            end if
            self%delayed(first) = k
            call rotate(self, k, self%pos(parent))
        end do
    end subroutine ordered_pivot

    ! Move the vertex at position from to position to; the vertices between
    ! keep their relative order, so the fill-reducing order is preserved.
    subroutine rotate(self, from, to)
        class(native_ldlt_z_t), intent(inout) :: self
        ! By value: callers pass entries of self%pos, which this routine updates.
        integer, value :: from, to
        integer :: vertex, i
        vertex = self%perm(from)
        if (from > to) then
            do i = from, to + 1, -1
                self%perm(i) = self%perm(i - 1)
                self%pos(self%perm(i)) = i
            end do
        else
            do i = from, to - 1
                self%perm(i) = self%perm(i + 1)
                self%pos(self%perm(i)) = i
            end do
        end if
        self%perm(to) = vertex
        self%pos(vertex) = to
    end subroutine rotate

    subroutine interchange(self, first, second)
        class(native_ldlt_z_t), intent(inout) :: self
        integer, intent(in) :: first, second
        integer :: vertex
        vertex = self%perm(first)
        self%perm(first) = self%perm(second)
        self%perm(second) = vertex
        self%pos(self%perm(first)) = first
        self%pos(self%perm(second)) = second
    end subroutine interchange

    subroutine ldlt_solve(self, b, x, status, operation)
        class(native_ldlt_z_t), intent(inout) :: self
        complex(dp), intent(in) :: b(:)
        complex(dp), intent(out) :: x(:)
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: operation
        call status_set(status, FORTSPARSE_NOT_FACTORED, 'LDLT: no usable factors')
        if (.not. self%ready) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'LDLT: solution shape')
        if (size(x) /= self%n) return
        call prepare_rhs(self, b, status, operation)
        if (.not. status_ok(status)) return
        call solve_workspace(self, status)
        if (status_ok(status)) then
            if (conjugate_operation(self, operation)) call conjugate_work(self)
            x = self%work
        end if
    end subroutine ldlt_solve

    subroutine ldlt_solve_inplace(self, b, status, operation)
        class(native_ldlt_z_t), intent(inout) :: self
        complex(dp), intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: operation
        call prepare_rhs(self, b, status, operation)
        if (.not. status_ok(status)) return
        call solve_workspace(self, status)
        if (status_ok(status)) then
            if (conjugate_operation(self, operation)) call conjugate_work(self)
            b = self%work
        end if
    end subroutine ldlt_solve_inplace

    subroutine prepare_rhs(self, b, status, operation)
        class(native_ldlt_z_t), intent(inout) :: self
        complex(dp), intent(in) :: b(:)
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: operation
        integer :: p
        call status_set(status, FORTSPARSE_NOT_FACTORED, 'LDLT: no usable factors')
        if (.not. self%ready) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'LDLT: RHS shape/value')
        if (present(operation)) then
            if (operation < 0 .or. operation > 2) return
        end if
        if (size(b) /= self%n) return
        do p = 1, self%n
            if (.not. lu_z_finite(b(p))) return
        end do
        ! N/T/H share factors through the declared class involution. Copy
        ! before output is touched, without passing aliasing dummies.
        self%work = b
        if (conjugate_operation(self, operation)) call conjugate_work(self)
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine prepare_rhs

    subroutine solve_workspace(self, status)
        class(native_ldlt_z_t), intent(inout) :: self
        type(fortsparse_status_t), intent(out) :: status
        integer :: k, p, first, second, row
        complex(dp) :: a, c, off, value, other
        do k = 1, self%n
            first = self%perm(k)
            do p = 1, self%lower(k)%used
                row = self%lower(k)%index(p)
                self%work(row) = lu_z_update(self%work(row), &
                                             self%lower(k)%value(p), self%work(first))
                if (.not. lu_z_finite(self%work(row))) then
                    call arithmetic_error(status)
                    return
                end if
            end do
        end do
        k = 1
        do while (k <= self%n)
            first = self%perm(k)
            a = self%diagonal(k)
            if (self%block(k) == 1) then
                self%work(first) = lu_z_divide(self%work(first), a)
                if (.not. lu_z_finite(self%work(first))) then
                    call arithmetic_error(status)
                    return
                end if
                k = k + 1
            else
                second = self%perm(k + 1)
                c = self%diagonal(k + 1)
                off = self%offdiagonal(k)
                call pair_solve(a, mate(off, self%is_hermitian), c, self%work(first), &
                             self%work(second), self%is_hermitian, value, other, status)
                if (.not. status_ok(status)) return
                self%work(first) = value
                self%work(second) = other
                k = k + 2
            end if
        end do
        do k = self%n, 1, -1
            first = self%perm(k)
            do p = 1, self%lower(k)%used
                row = self%lower(k)%index(p)
                self%work(first) = lu_z_update(self%work(first), &
                                      mate(self%lower(k)%value(p), self%is_hermitian), &
                                               self%work(row))
                if (.not. lu_z_finite(self%work(first))) then
                    call arithmetic_error(status)
                    return
                end if
            end do
        end do
        do p = 1, self%n
            if (.not. lu_z_finite(self%work(p))) then
                call arithmetic_error(status)
                return
            end if
        end do
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine solve_workspace

    subroutine ldlt_solve_block(self, b, x, status, operation)
        class(native_ldlt_z_t), intent(inout) :: self
        complex(dp), intent(in) :: b(:, :)
        complex(dp), intent(out) :: x(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: operation
        integer :: j
        call status_set(status, FORTSPARSE_NOT_FACTORED, 'LDLT: no usable factors')
        if (.not. self%ready) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'LDLT: block RHS shape')
        if (size(b, 1) /= self%n .or. size(x, 1) /= self%n) return
        if (size(b, 2) /= size(x, 2)) return
        if (present(operation)) then
            if (operation < 0 .or. operation > 2) return
        end if
        call status_set(status, FORTSPARSE_OK, '')
        do j = 1, size(b, 2)
            call self%solve(b(:, j), x(:, j), status, operation)
            if (.not. status_ok(status)) return
        end do
    end subroutine ldlt_solve_block

    pure integer function row_position(row, column) result(p)
        type(sparse_row_t), intent(in) :: row
        integer, intent(in) :: column
        integer :: lo, hi, mid
        lo = 1
        hi = row%used
        do while (lo <= hi)
            mid = lo + (hi - lo)/2
            if (row%index(mid) < column) then
                lo = mid + 1
            else
                hi = mid - 1
            end if
        end do
        p = lo
    end function row_position

    pure complex(dp) function row_get(row, column) result(value)
        type(sparse_row_t), intent(in) :: row
        integer, intent(in) :: column
        integer :: p
        value = 0.0_dp
        p = row_position(row, column)
        if (p > row%used) return
        if (row%index(p) == column) value = row%value(p)
    end function row_get

    subroutine row_put(row, column, value, low, limit, bytes, peak, required, status)
        type(sparse_row_t), intent(inout) :: row
        integer, intent(in) :: column
        complex(dp), intent(in) :: value
        logical, intent(in) :: low
        integer(i8), intent(in) :: limit
        integer(i8), intent(inout) :: bytes, peak, required
        type(fortsparse_status_t), intent(out) :: status
        integer :: p, q
        call status_set(status, FORTSPARSE_OK, '')
        p = row_position(row, column)
        if (p <= row%used) then
            if (row%index(p) == column) then
                row%value(p) = value
                return
            end if
        end if
        if (abs(value) <= 0.0_dp) return
        if (row%used == huge(0)) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, 'LDLT: row too long')
            return
        end if
        call reserve_row(row, row%used + 1, low, limit, bytes, peak, required, status)
        if (.not. status_ok(status)) return
        do q = row%used, p, -1
            row%index(q + 1) = row%index(q)
            row%value(q + 1) = row%value(q)
        end do
        row%index(p) = column
        row%value(p) = value
        row%used = row%used + 1
    end subroutine row_put

    subroutine row_erase(row, column)
        type(sparse_row_t), intent(inout) :: row
        integer, intent(in) :: column
        integer :: p, q
        p = row_position(row, column)
        if (p > row%used) return
        if (row%index(p) /= column) return
        do q = p, row%used - 1
            row%index(q) = row%index(q + 1)
            row%value(q) = row%value(q + 1)
        end do
        row%used = row%used - 1
    end subroutine row_erase

    subroutine reserve_row(row, need, low, limit, bytes, peak, required, status)
        type(sparse_row_t), intent(inout) :: row
        integer, intent(in) :: need
        logical, intent(in) :: low
        integer(i8), intent(in) :: limit
        integer(i8), intent(inout) :: bytes, peak, required
        type(fortsparse_status_t), intent(out) :: status
        integer :: capacity, wanted, ierr
        integer(i8) :: new_capacity, request
        integer, allocatable :: indices(:)
        complex(dp), allocatable :: values(:)
        capacity = 0
        if (allocated(row%index)) capacity = size(row%index)
        call status_set(status, FORTSPARSE_OK, '')
        if (need <= capacity) return
        new_capacity = int(need, i8)
        if (.not. low) new_capacity = max(new_capacity, &
                    min(int(huge(0), i8), int(capacity, i8) + max(4_i8, capacity/2_i8)))
        wanted = int(new_capacity)
        request = (ib + rb)*new_capacity
        ! Account both old and new payloads while copying a growing sparse row.
        call admit_bytes(request, limit, bytes, peak, required, status)
        if (.not. status_ok(status)) return
        allocate (indices(wanted), values(wanted), stat=ierr)
        if (ierr /= 0) then
            bytes = bytes - request
            call allocation_error(status)
            return
        end if
        if (row%used > 0) then
            indices(:row%used) = row%index(:row%used)
            values(:row%used) = row%value(:row%used)
        end if
        call release_row_payload(row, bytes)
        call move_alloc(indices, row%index)
        call move_alloc(values, row%value)
    end subroutine reserve_row

    subroutine admit(self, request, status)
        class(native_ldlt_z_t), intent(inout) :: self
        integer(i8), intent(in) :: request
        type(fortsparse_status_t), intent(out) :: status
        call admit_bytes(request, self%limit, self%bytes, self%peak, &
                         self%required, status)
    end subroutine admit

    subroutine admit_bytes(request, limit, bytes, peak, required, status)
        integer(i8), intent(in) :: request, limit
        integer(i8), intent(inout) :: bytes, peak, required
        type(fortsparse_status_t), intent(out) :: status
        call status_set(status, FORTSPARSE_OK, '')
        if (request < 0_i8) then
            required = huge(required)
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                            'LDLT: allocation byte count overflow')
            return
        end if
        if (bytes > huge(bytes) - request) then
            required = huge(required)
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                            'LDLT: allocation byte count overflow')
            return
        end if
        required = max(required, bytes + request)
        if (limit > 0_i8) then
            if (request > limit) then
                call status_set(status, FORTSPARSE_RESOURCE_LIMIT, &
                                'LDLT: solver-owned host memory cap exceeded')
                return
            end if
            if (bytes > limit - request) then
                call status_set(status, FORTSPARSE_RESOURCE_LIMIT, &
                                'LDLT: solver-owned host memory cap exceeded')
                return
            end if
        end if
        bytes = bytes + request
        peak = max(peak, bytes)
    end subroutine admit_bytes

    subroutine release_row_payload(row, bytes)
        type(sparse_row_t), intent(inout) :: row
        integer(i8), intent(inout) :: bytes
        if (allocated(row%index)) then
            bytes = bytes - ib*size(row%index, kind=i8)
            deallocate (row%index)
        end if
        if (allocated(row%value)) then
            bytes = bytes - rb*size(row%value, kind=i8)
            deallocate (row%value)
        end if
    end subroutine release_row_payload

    subroutine release_row(row, bytes)
        type(sparse_row_t), intent(inout) :: row
        integer(i8), intent(inout) :: bytes
        call release_row_payload(row, bytes)
        row%used = 0
    end subroutine release_row

    subroutine release_scratch(self)
        class(native_ldlt_z_t), intent(inout) :: self
        integer :: j
        if (allocated(self%schur)) then
            do j = 1, size(self%schur)
                call release_row(self%schur(j), self%bytes)
            end do
            self%bytes = self%bytes &
                     - int(storage_size(sparse_row_t())/8, i8)*size(self%schur, kind=i8)
            deallocate (self%schur)
        end if
        if (allocated(self%nbr)) then
            self%bytes = self%bytes - ib*size(self%nbr, kind=i8)
            deallocate (self%nbr)
        end if
        if (allocated(self%u)) then
            self%bytes = self%bytes - rb*size(self%u, kind=i8)
            deallocate (self%u)
        end if
        if (allocated(self%v)) then
            self%bytes = self%bytes - rb*size(self%v, kind=i8)
            deallocate (self%v)
        end if
        if (allocated(self%row_scale)) then
            self%bytes = self%bytes - rsb*size(self%row_scale, kind=i8)
            deallocate (self%row_scale)
        end if
        if (allocated(self%delayed)) then
            self%bytes = self%bytes - ib*size(self%delayed, kind=i8)
            deallocate (self%delayed)
        end if
    end subroutine release_scratch

    subroutine clear_numeric(self)
        class(native_ldlt_z_t), intent(inout) :: self
        self%ready = .false.
        self%inertia_ready = .false.
        if (allocated(self%schur)) deallocate (self%schur)
        if (allocated(self%lower)) deallocate (self%lower)
        if (allocated(self%perm)) deallocate (self%perm)
        if (allocated(self%pos)) deallocate (self%pos)
        if (allocated(self%block)) deallocate (self%block)
        if (allocated(self%nbr)) deallocate (self%nbr)
        if (allocated(self%diagonal)) deallocate (self%diagonal)
        if (allocated(self%offdiagonal)) deallocate (self%offdiagonal)
        if (allocated(self%u)) deallocate (self%u)
        if (allocated(self%v)) deallocate (self%v)
        if (allocated(self%work)) deallocate (self%work)
        if (allocated(self%row_scale)) deallocate (self%row_scale)
        if (allocated(self%delayed)) deallocate (self%delayed)
        self%bytes = 0_i8
        if (allocated(self%ap)) self%bytes = ib*size(self%ap, kind=i8)
        if (allocated(self%ai)) self%bytes = self%bytes + ib*size(self%ai, kind=i8)
        if (allocated(self%initial_order)) &
            self%bytes = self%bytes + ib*size(self%initial_order, kind=i8)
    end subroutine clear_numeric

    subroutine ldlt_free(self)
        class(native_ldlt_z_t), intent(inout) :: self
        call clear_numeric(self)
        if (allocated(self%ap)) deallocate (self%ap)
        if (allocated(self%ai)) deallocate (self%ai)
        if (allocated(self%initial_order)) deallocate (self%initial_order)
        self%n = 0
        self%bytes = 0_i8
        self%peak = 0_i8
        self%required = 0_i8
        self%limit = 0_i8
        self%low_memory = .false.
        self%is_hermitian = .true.
    end subroutine ldlt_free

    subroutine ldlt_finalize(self)
        type(native_ldlt_z_t), intent(inout) :: self
        call self%free()
    end subroutine ldlt_finalize

    subroutine ldlt_copy_column_order(self, n, ptr, rows, order, matched)
        class(native_ldlt_z_t), intent(in) :: self
        integer, intent(in) :: n, ptr(:), rows(:)
        integer, intent(out) :: order(:)
        logical, intent(out) :: matched
        matched = .false.
        if (.not. self%ready .or. n /= self%n) return
        if (.not. allocated(self%initial_order)) return
        if (.not. allocated(self%ap) .or. .not. allocated(self%ai)) return
        if (size(ptr) /= size(self%ap) .or. size(rows) /= size(self%ai)) return
        if (size(order) /= n) return
        if (any(ptr /= self%ap) .or. any(rows /= self%ai)) return
        order = self%initial_order
        matched = .true.
    end subroutine ldlt_copy_column_order

    pure integer(i8) function ldlt_bytes(self) result(bytes)
        class(native_ldlt_z_t), intent(in) :: self
        bytes = self%bytes
    end function ldlt_bytes

    pure integer(i8) function ldlt_peak(self) result(bytes)
        class(native_ldlt_z_t), intent(in) :: self
        bytes = self%peak
    end function ldlt_peak

    pure integer(i8) function ldlt_required(self) result(bytes)
        class(native_ldlt_z_t), intent(in) :: self
        ! Highest complete requested payload, including failed transient growth.
        bytes = self%required
    end function ldlt_required

    pure integer(i8) function ldlt_available(self) result(bytes)
        class(native_ldlt_z_t), intent(in) :: self
        ! Configured attributed-payload budget. Zero means unlimited.
        bytes = self%limit
    end function ldlt_available

    pure integer(i8) function ldlt_nnz(self) result(nnz)
        class(native_ldlt_z_t), intent(in) :: self
        integer :: j
        nnz = 0_i8
        if (.not. self%ready) return
        nnz = int(self%n, i8)
        do j = 1, self%n
            nnz = nnz + int(self%lower(j)%used, i8)
            if (self%block(j) == 2) nnz = nnz + 1_i8
        end do
    end function ldlt_nnz

    subroutine ldlt_permutation(self, permutation, status)
        class(native_ldlt_z_t), intent(in) :: self
        integer, intent(out) :: permutation(:)
        type(fortsparse_status_t), intent(out) :: status
        call status_set(status, FORTSPARSE_NOT_FACTORED, 'LDLT: no usable factors')
        if (.not. self%ready) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'LDLT: permutation size')
        if (size(permutation) /= self%n) return
        permutation = self%perm
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine ldlt_permutation

    pure elemental function mate(value, hermitian) result(result)
        complex(dp), intent(in) :: value
        logical, intent(in) :: hermitian
        complex(dp) :: result
        result = value
        if (hermitian) result = conjg(value)
    end function mate

    pure elemental real(dp) function pivot_norm(value) result(result)
        complex(dp), intent(in) :: value
        ! LAPACK complex Bunch--Kaufman uses component one-norm comparisons.
        ! norm_in_range guards this sum on all input and Schur entries.
        result = abs(real(value, dp)) + abs(aimag(value))
    end function pivot_norm

    pure elemental logical function norm_in_range(value) result(ok)
        complex(dp), intent(in) :: value
        ok = abs(real(value, dp)) <= huge(1.0_dp) - abs(aimag(value))
    end function norm_in_range

    subroutine conjugate_work(self)
        class(native_ldlt_z_t), intent(inout) :: self
        integer :: i
        do i = 1, self%n
            self%work(i) = conjg(self%work(i))
        end do
    end subroutine conjugate_work

    pure logical function conjugate_operation(self, operation) result(needed)
        class(native_ldlt_z_t), intent(in) :: self
        integer, intent(in), optional :: operation
        integer :: op
        op = 0
        if (present(operation)) op = operation
        needed = .false.
        if (self%is_hermitian) then
            needed = op == 1
        else
            needed = op == 2
        end if
    end function conjugate_operation

    subroutine pair_solve(a, b, c, u, v, hermitian, first, second, status)
        ! Row solve [first,second]*[[a,b],[mate(b),c]]=[u,v].
        ! Normalize block and RHS separately by exact powers of two; generated
        ! Cartesian polynomial leaves see bounded components. Do not form a
        ! determinant at the original physical scale.
        complex(dp), intent(in) :: a, b, c, u, v
        logical, intent(in) :: hermitian
        complex(dp), intent(out) :: first, second
        type(fortsparse_status_t), intent(out) :: status
        complex(dp) :: an, bn, cn, un, vn, det, n1, n2
        integer :: block_exponent, rhs_exponent, difference
        real(dp) :: block_norm, rhs_norm
        if (.not. lu_z_finite(u) .or. .not. lu_z_finite(v)) then
            call arithmetic_error(status)
            return
        end if
        block_norm = max(lu_z_norm(a), lu_z_norm(b), lu_z_norm(c))
        if (block_norm <= 0.0_dp) then
            call status_set(status, FORTSPARSE_SINGULAR, 'complex LDL: zero block')
            return
        end if
        block_exponent = exponent(block_norm)
        rhs_norm = max(lu_z_norm(u), lu_z_norm(v))
        rhs_exponent = exponent(rhs_norm)
        an = scale_complex(a, -block_exponent)
        bn = scale_complex(b, -block_exponent)
        cn = scale_complex(c, -block_exponent)
        un = scale_complex(u, -rhs_exponent)
        vn = scale_complex(v, -rhs_exponent)
        det = ldlt_z_det(an, bn, mate(bn, hermitian), cn)
        if (lu_z_norm(det) <= 0.0_dp) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                            'complex LDL: singular numerical 2x2 block')
            return
        end if
        call ldlt_z_numerators(an, bn, mate(bn, hermitian), cn, un, vn, n1, n2)
        difference = rhs_exponent - block_exponent
        first = scale_complex(lu_z_divide(n1, det), difference)
        second = scale_complex(lu_z_divide(n2, det), difference)
        if (.not. lu_z_finite(first) .or. .not. lu_z_finite(second)) then
            call arithmetic_error(status)
            return
        end if
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine pair_solve

    pure elemental function scale_complex(value, power) result(result)
        complex(dp), intent(in) :: value
        integer, intent(in) :: power
        complex(dp) :: result
        result = cmplx(scale(real(value, dp), power), scale(aimag(value), power), dp)
    end function scale_complex

    subroutine allocation_error(status)
        type(fortsparse_status_t), intent(out) :: status
        call status_set(status, FORTSPARSE_OUT_OF_MEMORY, &
                        'LDLT: host allocation failed')
    end subroutine allocation_error

    subroutine arithmetic_error(status)
        type(fortsparse_status_t), intent(out) :: status
        call status_set(status, FORTSPARSE_INTERNAL_ERROR, &
                        'LDLT: nonfinite arithmetic; factors or solution unavailable')
    end subroutine arithmetic_error
end module fortsparse_native_ldlt_z
