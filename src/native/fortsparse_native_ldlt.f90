module fortsparse_native_ldlt
    ! Serial sparse Bunch--Kaufman reference: P^T A P = L D L^T.
    ! Full canonical symmetric CSC is required. Schur rows and L columns use
    ! original vertex labels: symmetric pivoting changes maps, never n^2 data.
    ! Numeric factorization grows sparse rows under an explicit payload cap.
    ! This reference is not a production performance/parallelism admission.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csc, only: csc_t, csc_is_valid
    use fortsparse_status, only: fortsparse_status_t, status_set, status_ok, &
                        FORTSPARSE_OK, FORTSPARSE_SINGULAR, FORTSPARSE_INVALID_MATRIX, &
                                 FORTSPARSE_NOT_FACTORED, FORTSPARSE_OUT_OF_MEMORY, &
                                 FORTSPARSE_INDEX_OVERFLOW, FORTSPARSE_RESOURCE_LIMIT, &
                                 FORTSPARSE_INTERNAL_ERROR
    use fortsparse_spd_leaves, only: spd_divide, spd_update
    use fortsparse_ldlt_leaves, only: ldlt_det_value, ldlt_pair_value, &
                                      ldlt_update_value
    use fortsparse_inertia, only: inertia_count_t, inertia_add_pivot, &
                                  inertia_status_of, inertia_tolerance, &
                                  FORTSPARSE_INERTIA_UNAVAILABLE
    !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
    implicit none
    private
    public :: native_ldlt_t, ldlt_check_matrix
    integer(i8), parameter :: ib = storage_size(0)/8, rb = storage_size(0.0_dp)/8
    integer(i8), parameter :: rsb = storage_size(0.0_dp)/8
    real(dp), parameter :: alpha = (1.0_dp + sqrt(17.0_dp))/8.0_dp

    integer, parameter :: front_limit = 512, panel_width = 32

    type :: sparse_row_t
        integer :: used = 0
        integer, allocatable :: index(:)
        real(dp), allocatable :: value(:)
    end type sparse_row_t

    type :: native_ldlt_t
        private
        logical :: use_panels = .false.
        integer :: panels = 0, teams = 1, active_threads = 1, requested_threads = 0
        real(dp), allocatable :: front(:, :), panel_work(:, :), merge_value(:, :)
        integer, allocatable :: front_map(:), vertices(:), ipiv(:), merge_index(:, :)
        integer :: n = 0
        logical :: ready = .false., low_memory = .false.
        integer(i8) :: limit = 0_i8, bytes = 0_i8, peak = 0_i8
        integer(i8) :: required = 0_i8
        ! Inertia of the last completed elimination (also with zero pivots).
        logical :: inertia_ready = .false.
        type(inertia_count_t) :: pivots
        ! Per-vertex max |entry| of A and every Schur entry stored in its row.
        real(dp), allocatable :: row_scale(:)
        integer, allocatable :: ap(:), ai(:), perm(:), pos(:), block(:), nbr(:), initial_order(:)
        real(dp), allocatable :: diagonal(:), offdiagonal(:), u(:), v(:), work(:)
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
        procedure :: blocked_panels => ldlt_blocked_panels
        procedure :: numeric_threads => ldlt_numeric_threads
        final :: ldlt_finalize
    end type native_ldlt_t
contains
    subroutine ldlt_factor(self, a, status, host_limit_bytes, low_memory, column_permutation, blocked, nthreads)
        class(native_ldlt_t), intent(inout) :: self
        type(csc_t), intent(in) :: a
        type(fortsparse_status_t), intent(out) :: status
        integer(i8), intent(in), optional :: host_limit_bytes
        logical, intent(in), optional :: low_memory
        integer, intent(in), optional :: column_permutation(:)
        logical, intent(in), optional :: blocked
        integer, intent(in), optional :: nthreads
        call self%free()
        if (present(host_limit_bytes)) self%limit = host_limit_bytes
        if (present(low_memory)) self%low_memory = low_memory
        call ldlt_check_matrix(a, status)
        if (.not. status_ok(status)) return
        if (self%limit < 0_i8) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, 'LDLT: negative cap')
            return
        end if
        if (present(blocked)) self%use_panels = blocked
        if (present(nthreads)) self%requested_threads = nthreads
        call initialize(self, a, status, column_permutation)
        if (.not. status_ok(status)) return
        call panel_workspace(self, status)
        if (.not. status_ok(status)) return
        call numeric(self, status)
    end subroutine ldlt_factor

    subroutine ldlt_refactor(self, a, status, host_limit_bytes, low_memory, column_permutation, blocked, nthreads)
        class(native_ldlt_t), intent(inout) :: self
        type(csc_t), intent(in) :: a
        type(fortsparse_status_t), intent(out) :: status
        integer(i8), intent(in), optional :: host_limit_bytes
        logical, intent(in), optional :: low_memory
        integer, intent(in), optional :: column_permutation(:)
        logical, intent(in), optional :: blocked
        integer, intent(in), optional :: nthreads
        integer :: i
        call ldlt_check_matrix(a, status)
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
        if (present(blocked)) self%use_panels = blocked
        if (present(nthreads)) self%requested_threads = nthreads
        self%peak = self%bytes
        self%required = self%bytes
        call initialize(self, a, status)
        if (.not. status_ok(status)) return
        call panel_workspace(self, status)
        if (.not. status_ok(status)) return
        call numeric(self, status)
    end subroutine ldlt_refactor

    subroutine ldlt_check_matrix(a, status)
        type(csc_t), intent(in) :: a
        type(fortsparse_status_t), intent(out) :: status
        integer :: j, p, i, previous
        real(dp) :: counterpart
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        'LDLT: finite, full canonical symmetric square CSC required')
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
                if (.not. ieee_is_finite(a%val(p))) return
                ! A stored zero matches a missing transpose entry. Any
                ! asymmetric nonzero is still checked from its own column.
                if (abs(a%val(p)) <= 0.0_dp) cycle
                counterpart = csc_get(a, j, i)
                if (a%val(p) < counterpart .or. a%val(p) > counterpart) return
            end do
        end do
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine ldlt_check_matrix

    pure real(dp) function csc_get(a, row, column) result(value)
        type(csc_t), intent(in) :: a
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
        class(native_ldlt_t), intent(inout) :: self
        type(csc_t), intent(in) :: a
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
            call status_set(status, FORTSPARSE_INVALID_MATRIX, 'LDLT: invalid initial order')
            if (size(column_permutation) /= n) return
            self%pos = 0
            do j = 1, n
                p = column_permutation(j)
                if (p < 1 .or. p > n) return
                if (self%pos(p) /= 0) return
                self%pos(p) = j
            end do
            if (.not. allocated(self%initial_order)) then
                call admit(self, ib*int(n,i8), status)
                if (.not. status_ok(status)) return
                allocate(self%initial_order(n), stat=ierr)
                if (ierr /= 0) then
                    call allocation_error(status)
                    call self%free()
                    return
                end if
                self%initial_order = column_permutation
            end if
        end if
        do j = 1, n
            self%perm(j) = j
            if (allocated(self%initial_order)) self%perm(j) = self%initial_order(j)
            self%pos(self%perm(j)) = j
        end do
        do j = 1, n
            ! Stored zeros belong to the immutable CSC pattern/reuse key, not
            ! to the active numeric graph. Keeping them here enlarges relaxed
            ! fronts and makes every pivot scan/erase traverse irrelevant rows.
            count = 0
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                if (abs(a%val(p)) > 0.0_dp) count = count + 1
            end do
            call reserve_row(self%schur(j), count, .true., self%limit, &
                             self%bytes, self%peak, self%required, status)
            if (.not. status_ok(status)) return
            self%schur(j)%used = count
            count = 0
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                if (abs(a%val(p)) <= 0.0_dp) cycle
                count = count + 1
                self%schur(j)%index(count) = a%row_idx(p)
                self%schur(j)%value(count) = a%val(p)
            end do
        end do
        self%block = 0
        self%diagonal = 0.0_dp
        self%offdiagonal = 0.0_dp
        do j = 1, n
            self%row_scale(j) = 0.0_dp
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                self%row_scale(j) = max(self%row_scale(j), abs(a%val(p)))
            end do
        end do
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine initialize

    ! A relaxed supernode contains consecutive candidate pivots and the complete
    ! union of their active rows. DLASYF performs genuine Bunch--Kaufman tests
    ! against that complete support and uses BLAS-3 for its Schur update. An
    ! interchange rejects this private trial before any authoritative state is
    ! changed; the unrestricted sparse reference then eliminates one pivot.
    subroutine panel_workspace(self, status)
        class(native_ldlt_t), intent(inout) :: self
        type(fortsparse_status_t), intent(out) :: status
        integer :: ierr
        integer(i8) :: request
        call status_set(status, FORTSPARSE_OK, '')
        if (.not. self%use_panels) return
        ! Opaque LAPACK/BLAS workspace is not a hard allocation-cap contract.
        ! Keep the existing audited scalar algorithm for strict/LOW requests.
        if (self%limit > 0_i8 .or. self%low_memory .or. self%n < 8) then
            self%use_panels = .false.
            return
        end if
        self%teams = 1
        !$ self%teams = omp_get_max_threads()
        if (self%requested_threads > 0) self%teams = self%requested_threads
        self%teams = max(1, min(self%teams, self%n))
        request = rb*int(front_limit, i8)*(front_limit + panel_width) &
            + ib*(int(self%n, i8) + 2_i8*front_limit) &
            + (ib + rb)*int(self%n, i8)*self%teams
        call admit(self, request, status)
        if (.not. status_ok(status)) return
        allocate(self%front(front_limit, front_limit), &
            self%panel_work(front_limit, panel_width), &
            self%vertices(front_limit), self%ipiv(front_limit), &
            self%front_map(self%n), self%merge_index(self%n, self%teams), &
            self%merge_value(self%n, self%teams), stat=ierr)
        if (ierr /= 0) then
            call allocation_error(status)
            call clear_numeric(self)
        end if
    end subroutine panel_workspace

    subroutine free_panel_workspace(self)
        class(native_ldlt_t), intent(inout) :: self
        if (allocated(self%front)) then
            self%bytes = self%bytes - rb*size(self%front, kind=i8)
            deallocate(self%front)
        end if
        if (allocated(self%panel_work)) then
            self%bytes = self%bytes - rb*size(self%panel_work, kind=i8)
            deallocate(self%panel_work)
        end if
        if (allocated(self%front_map)) then
            self%bytes = self%bytes - ib*size(self%front_map, kind=i8)
            deallocate(self%front_map)
        end if
        if (allocated(self%vertices)) then
            self%bytes = self%bytes - ib*size(self%vertices, kind=i8)
            deallocate(self%vertices)
        end if
        if (allocated(self%ipiv)) then
            self%bytes = self%bytes - ib*size(self%ipiv, kind=i8)
            deallocate(self%ipiv)
        end if
        if (allocated(self%merge_index)) then
            self%bytes = self%bytes - ib*size(self%merge_index, kind=i8)
            deallocate(self%merge_index)
        end if
        if (allocated(self%merge_value)) then
            self%bytes = self%bytes - rb*size(self%merge_value, kind=i8)
            deallocate(self%merge_value)
        end if
    end subroutine free_panel_workspace

    subroutine panel_factor(self, k, eliminated, status)
        class(native_ldlt_t), intent(inout) :: self
        integer, intent(in) :: k
        integer, intent(out) :: eliminated
        type(fortsparse_status_t), intent(out) :: status
        integer :: width, m, p, q, i, j, r, vertex, nb, kb, info, step, count
        real(dp) :: a, b, c, lm1, lm2, bound1, bound2, l1, l2
        interface
            subroutine dlasyf(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info)
                import dp
                character, intent(in) :: uplo
                integer, intent(in) :: n, nb, lda, ldw
                integer, intent(out) :: kb, ipiv(*), info
                real(dp), intent(inout) :: a(lda, *), w(ldw, *)
            end subroutine dlasyf
        end interface
        call status_set(status, FORTSPARSE_OK, '')
        eliminated = 0
        width = min(panel_width, self%n - k + 1)
        if (width < 8) return
        self%front_map = 0
        m = width
        do p = 1, width
            vertex = self%perm(k + p - 1)
            self%vertices(p) = vertex
            self%front_map(vertex) = p
        end do
        do p = 1, width
            i = self%vertices(p)
            do q = 1, self%schur(i)%used
                j = self%schur(i)%index(q)
                if (self%pos(j) < k .or. self%front_map(j) /= 0) cycle
                if (m == front_limit) return
                m = m + 1
                self%vertices(m) = j
                self%front_map(j) = m
            end do
        end do
        ! Each assembly column is private. No numerical-loop allocation.
        !$omp parallel num_threads(self%teams) private(p,i,q,j,r)
        !$omp single
        !$ self%active_threads = omp_get_num_threads()
        !$omp end single
        !$omp do schedule(static)
        do p = 1, m
            self%front(:m, p) = 0.0_dp
            i = self%vertices(p)
            do q = 1, self%schur(i)%used
                j = self%schur(i)%index(q)
                r = self%front_map(j)
                if (r > 0) self%front(r, p) = self%schur(i)%value(q)
            end do
        end do
        !$omp end do
        !$omp end parallel
        nb = width
        call dlasyf('L', m, nb, kb, self%front, front_limit, &
            self%ipiv, self%panel_work, front_limit, info)
        if (info /= 0 .or. kb < 1) return
        ! A trial can inspect a noncandidate pivot row with incomplete support,
        ! but acceptance requires that NO interchange occurred, at any step.
        ! Thus every accepted pivot search used a fully assembled candidate row.
        p = 1
        do while (p <= kb)
            if (self%ipiv(p) > 0) then
                if (self%ipiv(p) /= p) return
                p = p + 1
            else
                if (p == kb) return
                if (self%ipiv(p) /= -p - 1) return
                if (self%ipiv(p + 1) /= self%ipiv(p)) return
                p = p + 2
            end if
        end do
        do p = 1, m
            do q = 1, p
                if (.not. ieee_is_finite(self%front(p, q))) return
            end do
        end do
        ! Reserve all persistent growth before committing the trial or entering
        ! parallel scatter. Allocation failures leave factors unavailable.
        p = 1
        do while (p <= kb)
            step = 1
            if (self%ipiv(p) < 0) step = 2
            do j = p, p + step - 1
                call reserve_row(self%lower(k + j - 1), m - p - step + 1, &
                    .true., self%limit, self%bytes, self%peak, self%required, status)
                if (.not. status_ok(status)) return
            end do
            p = p + step
        end do
        do p = kb + 1, m
            i = self%vertices(p)
            call reserve_row(self%schur(i), &
                int(min(int(self%n, i8), int(self%schur(i)%used, i8) + m - kb)), &
                .false., self%limit, self%bytes, self%peak, self%required, status)
            if (.not. status_ok(status)) return
        end do
        p = 1
        do while (p <= kb)
            step = 1
            if (self%ipiv(p) < 0) step = 2
            a = self%front(p, p)
            b = 0.0_dp
            c = 0.0_dp
            if (step == 2) then
                b = self%front(p + 1, p)
                c = self%front(p + 1, p + 1)
            end if
            i = self%vertices(p)
            j = self%vertices(p + step - 1)
            call inertia_add_pivot(self%pivots, step, a, abs(b), c, &
                inertia_tolerance(self%n, max(self%row_scale(i), self%row_scale(j))))
            self%diagonal(k + p - 1) = a
            self%offdiagonal(k + p - 1) = b
            self%block(k + p - 1) = step
            if (step == 2) self%diagonal(k + p) = c
            ! Conservative history bound: |S_ij - (L D L^T)_ij| <=
            ! |S_ij| + |L_i1| max_j |(D L^T)_1j| + |L_i2| max_j |(D L^T)_2j|.
            ! This includes intermediate updates hidden by LAPACK's panel and
            ! can only make the existing reliability test more conservative.
            lm1 = 0.0_dp
            lm2 = 0.0_dp
            if (p + step <= m) then
                lm1 = maxval(abs(self%front(p + step:m, p)))
                if (step == 2) lm2 = maxval(abs(self%front(p + step:m, p + 1)))
            end if
            bound1 = abs(a)*lm1 + abs(b)*lm2
            bound2 = abs(b)*lm1 + abs(c)*lm2
            count = 0
            do r = p + step, m
                l1 = self%front(r, p)
                l2 = 0.0_dp
                if (step == 2) l2 = self%front(r, p + 1)
                vertex = self%vertices(r)
                self%row_scale(vertex) = self%row_scale(vertex) &
                    + abs(l1)*bound1 + abs(l2)*bound2
                if (abs(l1) <= 0.0_dp .and. abs(l2) <= 0.0_dp) cycle
                count = count + 1
                self%lower(k + p - 1)%index(count) = vertex
                self%lower(k + p - 1)%value(count) = l1
                if (step == 2) then
                    self%lower(k + p)%index(count) = vertex
                    self%lower(k + p)%value(count) = l2
                end if
            end do
            self%lower(k + p - 1)%used = count
            if (step == 2) self%lower(k + p)%used = count
            p = p + step
        end do
        ! Sorted trailing vertices support a linear union with each Schur row.
        count = 0
        do i = 1, self%n
            if (self%front_map(i) <= kb) cycle
            count = count + 1
            self%nbr(count) = i
        end do
        !$omp parallel do num_threads(self%teams) schedule(static)
        do p = kb + 1, m
            call panel_scatter_row(self, p, kb, count)
        end do
        !$omp end parallel do
        do p = 1, kb
            call release_row(self%schur(self%vertices(p)), self%bytes)
        end do
        self%panels = self%panels + 1
        eliminated = kb
    end subroutine panel_factor

    subroutine panel_scatter_row(self, p, kb, count)
        class(native_ldlt_t), intent(inout) :: self
        integer, intent(in) :: p, kb, count
        integer :: tid, i, j, q, r, old, fresh, used, column, local
        real(dp) :: value
        logical :: existed
        tid = 1
        !$ tid = omp_get_thread_num() + 1
        i = self%vertices(p)
        old = 1
        fresh = 1
        used = 0
        do while (old <= self%schur(i)%used .or. fresh <= count)
            j = self%n + 1
            q = self%n + 1
            if (old <= self%schur(i)%used) j = self%schur(i)%index(old)
            if (fresh <= count) q = self%nbr(fresh)
            column = min(j, q)
            existed = j == column
            value = 0.0_dp
            if (existed) then
                value = self%schur(i)%value(old)
                old = old + 1
            end if
            if (q == column) fresh = fresh + 1
            local = self%front_map(column)
            if (local > 0 .and. local <= kb) cycle
            if (local > kb) value = self%front(max(p, local), min(p, local))
            if (.not. existed .and. abs(value) <= 0.0_dp) cycle
            used = used + 1
            self%merge_index(used, tid) = column
            self%merge_value(used, tid) = value
        end do
        do r = 1, used
            self%schur(i)%index(r) = self%merge_index(r, tid)
            self%schur(i)%value(r) = self%merge_value(r, tid)
        end do
        self%schur(i)%used = used
    end subroutine panel_scatter_row

    integer function ldlt_blocked_panels(self) result(panels)
        class(native_ldlt_t), intent(in) :: self
        panels = self%panels
    end function ldlt_blocked_panels

    integer function ldlt_numeric_threads(self) result(threads)
        class(native_ldlt_t), intent(in) :: self
        threads = self%active_threads
    end function ldlt_numeric_threads

    subroutine numeric(self, status)
        class(native_ldlt_t), intent(inout) :: self
        type(fortsparse_status_t), intent(out) :: status
        integer :: k, step, pivot, first, second, m, i, j, p, q, eliminated
        real(dp) :: a, b, c, scale, determinant, value, l1, l2
        self%pivots = inertia_count_t()
        k = 1
        do while (k <= self%n)
            if (self%use_panels) then
                call panel_factor(self, k, eliminated, status)
                if (.not. status_ok(status)) return
                if (eliminated > 0) then
                    k = k + eliminated
                    cycle
                end if
            end if
            call choose_pivot(self, k, pivot, step, status)
            if (.not. status_ok(status)) return
            call interchange(self, k + step - 1, pivot)
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
            ! Local rounding threshold from the pivot rows' entry history.
            call inertia_add_pivot(self%pivots, step, a, abs(b), c, &
                inertia_tolerance(self%n, max(self%row_scale(first), &
                                              self%row_scale(second))))
            scale = max(abs(a), abs(b), abs(c))
            if (step == 2) then
                call ldlt_det_value(spd_divide(a, scale), &
                                spd_divide(b, scale), spd_divide(c, scale), determinant)
                if (abs(determinant) <= 0.0_dp) then
                    call status_set(status, FORTSPARSE_SINGULAR, &
                                    'LDLT: singular numerical 2x2 pivot')
                    return
                end if
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
                        value = self%schur(first)%value(p)
                        p = p + 1
                    end if
                end if
                if (q <= self%schur(second)%used) then
                    if (self%schur(second)%index(q) == i) then
                        l2 = self%schur(second)%value(q)
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
                    l1 = spd_divide(self%u(p), a)
                else
                    call ldlt_pair_value(spd_divide(a, scale), &
                                         spd_divide(b, scale), spd_divide(c, scale), &
                           spd_divide(self%u(p), scale), spd_divide(self%v(p), scale), &
                                         l1, l2)
                end if
                if (.not. ieee_is_finite(l1) .or. .not. ieee_is_finite(l2)) then
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
                    call ldlt_update_value(row_get(self%schur(i), j), &
                                           l1, self%u(q), l2, self%v(q), value)
                    if (.not. ieee_is_finite(value)) then
                        call arithmetic_error(status)
                        return
                    end if
                    self%row_scale(i) = max(self%row_scale(i), abs(value))
                    self%row_scale(j) = max(self%row_scale(j), abs(value))
                    call row_put(self%schur(i), j, value, self%low_memory, &
                               self%limit, self%bytes, self%peak, self%required, status)
                    if (.not. status_ok(status)) return
                    if (i /= j) then
                        call row_put(self%schur(j), i, value, self%low_memory, &
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
        self%inertia_ready = .true.
        if (self%pivots%zero > 0_i8) then
            ! Zero pivot columns were eliminated as exact zero 1x1 pivots so
            ! that inertia is complete; D is singular and solves are refused.
            call status_set(status, FORTSPARSE_SINGULAR, &
                            'LDLT: zero numerical pivot column; inertia only')
            return
        end if
        self%ready = .true.
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine numeric

    ! Inertia of the last completed factorization. state is a
    ! FORTSPARSE_INERTIA_* code; UNAVAILABLE when no elimination completed.
    subroutine ldlt_inertia(self, count, state, tolerance)
        class(native_ldlt_t), intent(in) :: self
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
        class(native_ldlt_t), intent(in) :: self
        integer, intent(in) :: k
        integer, intent(out) :: pivot, step
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: absakk, colmax, rowmax, magnitude
        integer :: first, largest, i, p
        first = self%perm(k)
        absakk = abs(row_get(self%schur(first), first))
        largest = first
        colmax = 0.0_dp
        do p = 1, self%schur(first)%used
            i = self%schur(first)%index(p)
            if (self%pos(i) <= k) cycle
            magnitude = abs(self%schur(first)%value(p))
            if (magnitude > colmax) then
                colmax = magnitude
                largest = i
            end if
        end do
        pivot = k
        step = 1
        call status_set(status, FORTSPARSE_OK, '')
        ! A zero column is an exact zero 1x1 pivot (counted in the inertia);
        ! numeric reports SINGULAR after completing the elimination.
        if (max(absakk, colmax) <= 0.0_dp) return
        if (absakk >= alpha*colmax) return
        rowmax = 0.0_dp
        do p = 1, self%schur(largest)%used
            i = self%schur(largest)%index(p)
            if (self%pos(i) < k .or. i == largest) cycle
            rowmax = max(rowmax, abs(self%schur(largest)%value(p)))
        end do
        ! rowmax >= colmax > 0: quotient is safe and never exceeds one.
        if (absakk > 0.0_dp) then
            if (spd_divide(absakk, colmax) >= &
                alpha*spd_divide(colmax, rowmax)) return
        end if
        pivot = self%pos(largest)
        if (abs(row_get(self%schur(largest), largest)) >= alpha*rowmax) return
        step = 2
    end subroutine choose_pivot

    subroutine interchange(self, first, second)
        class(native_ldlt_t), intent(inout) :: self
        integer, intent(in) :: first, second
        integer :: vertex
        vertex = self%perm(first)
        self%perm(first) = self%perm(second)
        self%perm(second) = vertex
        self%pos(self%perm(first)) = first
        self%pos(self%perm(second)) = second
    end subroutine interchange

    subroutine ldlt_solve(self, b, x, status, transpose)
        class(native_ldlt_t), intent(inout) :: self
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: x(:)
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(in), optional :: transpose
        call status_set(status, FORTSPARSE_NOT_FACTORED, 'LDLT: no usable factors')
        if (.not. self%ready) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'LDLT: solution shape')
        if (size(x) /= self%n) return
        call prepare_rhs(self, b, status, transpose)
        if (.not. status_ok(status)) return
        call solve_workspace(self, status)
        if (status_ok(status)) x = self%work
    end subroutine ldlt_solve

    subroutine ldlt_solve_inplace(self, b, status, transpose)
        class(native_ldlt_t), intent(inout) :: self
        real(dp), intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(in), optional :: transpose
        call prepare_rhs(self, b, status, transpose)
        if (.not. status_ok(status)) return
        call solve_workspace(self, status)
        if (status_ok(status)) b = self%work
    end subroutine ldlt_solve_inplace

    subroutine prepare_rhs(self, b, status, transpose)
        class(native_ldlt_t), intent(inout) :: self
        real(dp), intent(in) :: b(:)
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(in), optional :: transpose
        integer :: p
        call status_set(status, FORTSPARSE_NOT_FACTORED, 'LDLT: no usable factors')
        if (.not. self%ready) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'LDLT: RHS shape/value')
        if (present(transpose)) then
            if (transpose) call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                                           'LDLT: transpose RHS shape/value')
        end if
        if (size(b) /= self%n) return
        do p = 1, self%n
            if (.not. ieee_is_finite(b(p))) return
        end do
        ! Real symmetry makes normal/transpose solves identical. The input is
        ! copied before output is touched, without passing aliasing dummies.
        self%work = b
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine prepare_rhs

    subroutine solve_workspace(self, status)
        class(native_ldlt_t), intent(inout) :: self
        type(fortsparse_status_t), intent(out) :: status
        integer :: k, p, first, second, row
        real(dp) :: a, c, off, scale, value, other
        do k = 1, self%n
            first = self%perm(k)
            do p = 1, self%lower(k)%used
                row = self%lower(k)%index(p)
                self%work(row) = spd_update(self%work(row), &
                                            self%lower(k)%value(p), self%work(first))
            end do
        end do
        k = 1
        do while (k <= self%n)
            first = self%perm(k)
            a = self%diagonal(k)
            if (self%block(k) == 1) then
                self%work(first) = spd_divide(self%work(first), a)
                k = k + 1
            else
                second = self%perm(k + 1)
                c = self%diagonal(k + 1)
                off = self%offdiagonal(k)
                scale = max(abs(a), abs(c), abs(off))
                call ldlt_pair_value(spd_divide(a, scale), &
                                     spd_divide(off, scale), spd_divide(c, scale), &
                                     spd_divide(self%work(first), scale), &
                                     spd_divide(self%work(second), scale), value, other)
                self%work(first) = value
                self%work(second) = other
                k = k + 2
            end if
        end do
        do k = self%n, 1, -1
            first = self%perm(k)
            do p = 1, self%lower(k)%used
                row = self%lower(k)%index(p)
                self%work(first) = spd_update(self%work(first), &
                                              self%lower(k)%value(p), self%work(row))
            end do
        end do
        do p = 1, self%n
            if (.not. ieee_is_finite(self%work(p))) then
                call arithmetic_error(status)
                return
            end if
        end do
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine solve_workspace

    subroutine ldlt_solve_block(self, b, x, status, transpose)
        class(native_ldlt_t), intent(inout) :: self
        real(dp), intent(in) :: b(:, :)
        real(dp), intent(out) :: x(:, :)
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(in), optional :: transpose
        integer :: j
        call status_set(status, FORTSPARSE_NOT_FACTORED, 'LDLT: no usable factors')
        if (.not. self%ready) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'LDLT: block RHS shape')
        if (size(b, 1) /= self%n .or. size(x, 1) /= self%n) return
        if (size(b, 2) /= size(x, 2)) return
        call status_set(status, FORTSPARSE_OK, '')
        do j = 1, size(b, 2)
            call self%solve(b(:, j), x(:, j), status, transpose)
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

    pure real(dp) function row_get(row, column) result(value)
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
        real(dp), intent(in) :: value
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
        real(dp), allocatable :: values(:)
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
        class(native_ldlt_t), intent(inout) :: self
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
        class(native_ldlt_t), intent(inout) :: self
        integer :: j
        call free_panel_workspace(self)
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
    end subroutine release_scratch

    subroutine clear_numeric(self)
        class(native_ldlt_t), intent(inout) :: self
        self%ready = .false.
        self%inertia_ready = .false.
        call free_panel_workspace(self)
        self%panels = 0
        self%active_threads = 1
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
        self%bytes = 0_i8
        if (allocated(self%ap)) self%bytes = ib*size(self%ap, kind=i8)
        if (allocated(self%ai)) self%bytes = self%bytes + ib*size(self%ai, kind=i8)
        if (allocated(self%initial_order)) &
            self%bytes = self%bytes + ib*size(self%initial_order, kind=i8)
    end subroutine clear_numeric

    subroutine ldlt_free(self)
        class(native_ldlt_t), intent(inout) :: self
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
        self%use_panels = .false.
        self%requested_threads = 0
    end subroutine ldlt_free

    subroutine ldlt_finalize(self)
        type(native_ldlt_t), intent(inout) :: self
        call self%free()
    end subroutine ldlt_finalize

    subroutine ldlt_copy_column_order(self, n, ptr, rows, order, matched)
        class(native_ldlt_t), intent(in) :: self
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
        class(native_ldlt_t), intent(in) :: self
        bytes = self%bytes
    end function ldlt_bytes

    pure integer(i8) function ldlt_peak(self) result(bytes)
        class(native_ldlt_t), intent(in) :: self
        bytes = self%peak
    end function ldlt_peak

    pure integer(i8) function ldlt_required(self) result(bytes)
        class(native_ldlt_t), intent(in) :: self
        ! Highest complete requested payload, including failed transient growth.
        bytes = self%required
    end function ldlt_required

    pure integer(i8) function ldlt_available(self) result(bytes)
        class(native_ldlt_t), intent(in) :: self
        ! Configured attributed-payload budget. Zero means unlimited.
        bytes = self%limit
    end function ldlt_available

    pure integer(i8) function ldlt_nnz(self) result(nnz)
        class(native_ldlt_t), intent(in) :: self
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
        class(native_ldlt_t), intent(in) :: self
        integer, intent(out) :: permutation(:)
        type(fortsparse_status_t), intent(out) :: status
        call status_set(status, FORTSPARSE_NOT_FACTORED, 'LDLT: no usable factors')
        if (.not. self%ready) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'LDLT: permutation size')
        if (size(permutation) /= self%n) return
        permutation = self%perm
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine ldlt_permutation

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
end module fortsparse_native_ldlt
