module fortsparse_native
    ! Project-owned simplicial Cholesky. Graph reach determines fill;
    ! generated leaves own continuous arithmetic. Full canonical symmetric
    ! CSC is required; the caller's matrix is never modified.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_c_binding, only: c_int, c_int32_t, c_size_t, &
        c_ptr, c_null_ptr, c_associated, c_f_pointer, c_sizeof
    use fortsparse_kinds, only: dp, i8
    use fortsparse_backend, only: sparse_backend_t, sparse_factor_stats_t, &
        FORTSPARSE_SYM_SPD, FORTSPARSE_ORDER_DEFAULT, FORTSPARSE_ORDER_ND, &
        FORTSPARSE_ORDER_AMD
    use fortsparse_csc, only: csc_t, csc_z_t, csc_is_valid
    use fortsparse_status, only: fortsparse_status_t, status_set, status_ok, &
        FORTSPARSE_OK, FORTSPARSE_INVALID_MATRIX, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_OUT_OF_MEMORY, &
        FORTSPARSE_INDEX_OVERFLOW, FORTSPARSE_NOT_POSITIVE_DEFINITE
    use fortsparse_native_supernodes, only: supernodes_t
    use fortsparse_native_symbolic, only: native_symbolic_counts
    use fortsparse_native_cuda, only: native_cuda_available
    use fortsparse_spd_leaves, only: spd_root, spd_divide, spd_update, &
        spd_scatter_update, spd_gather_update
    implicit none
    private
    public :: native_backend_t
    interface
        function symbolic_alloc(bytes) result(owner) &
                bind(C, name='fsparse_native_symbolic_alloc')
            import :: c_ptr, c_size_t
            integer(c_size_t), value :: bytes
            type(c_ptr) :: owner
        end function symbolic_alloc
        subroutine symbolic_free(owner, bytes) bind(C, name='fsparse_native_symbolic_free')
            import :: c_ptr, c_size_t
            type(c_ptr), value :: owner
            integer(c_size_t), value :: bytes
        end subroutine symbolic_free
        integer(c_int) function native_order(n, nnz, ptr, row, perm, inverse) &
                bind(C, name='fsparse_native_order')
            import :: c_int
            integer(c_int), value :: n, nnz
            integer(c_int), intent(in) :: ptr(*), row(*)
            integer(c_int), intent(out) :: perm(*), inverse(*)
        end function native_order
        integer(c_int) function native_order_amd(n, nnz, ptr, row, perm, inverse) &
                bind(C, name='fsparse_native_order_amd')
            import :: c_int
            integer(c_int), value :: n, nnz
            integer(c_int), intent(in) :: ptr(*), row(*)
            integer(c_int), intent(out) :: perm(*), inverse(*)
        end function native_order_amd
    end interface

    type, extends(sparse_backend_t) :: native_backend_t
        integer :: n = 0
        logical :: reference_only = .false.
        integer :: selected_order = FORTSPARSE_ORDER_DEFAULT
        integer :: analyzed_order = 0, analyzed_threads = 0, analyzed_device = -1
        integer, allocatable :: perm(:), inverse(:), source(:), op(:), oi(:)
        type(csc_t) :: ordered
        logical :: ready = .false., blocked = .false.
        type(supernodes_t) :: panels
        integer, allocatable :: ap(:), ai(:), parent(:), lp(:)
        integer(c_int32_t), contiguous, pointer :: li(:) => null()
        type(c_ptr), private :: li_owner = c_null_ptr
        integer(c_size_t), private :: li_bytes = 0_c_size_t
        integer, allocatable :: count(:), mark(:), reach(:), path(:)
        real(dp), allocatable :: lv(:), work(:)
        real(dp), pointer :: user_vector(:) => null()
        type(sparse_factor_stats_t) :: stats
    contains
        procedure :: factor_real => native_factor
        procedure :: refactor_real => native_refactor
        procedure :: factor_complex => native_factor_complex
        procedure :: solve_real => native_solve
        procedure :: solve_real_inplace => native_solve_inplace
        procedure :: solve_real_block => native_solve_block
        procedure :: solve_complex => native_solve_complex
        procedure :: solve_complex_inplace => native_solve_complex_inplace
        procedure :: factor_stats => native_stats
        procedure :: vector => native_vector
        procedure :: free => native_free
        final :: native_finalize
    end type native_backend_t

contains

    subroutine native_factor(self, A, refine, status)
        class(native_backend_t), intent(inout) :: self
        type(csc_t), intent(in) :: A
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: start

        call self%free()
        call validate(self, A, status)
        if (.not. status_ok(status)) return
        start = wall_time()
        self%analyzed_order = self%options%ordering
        self%analyzed_threads = self%options%nthreads
        self%analyzed_device = self%options%device
        call prepare_ordering(self, A, status)
        if (.not. status_ok(status)) return
        if (allocated(self%perm)) then
            call analyze(self, self%ordered, status)
        else
            call analyze(self, A, status)
        end if
        if (.not. status_ok(status)) return
        if (self%blocked) then
            if (allocated(self%perm)) then
                ! Panels assemble directly from the caller's canonical CSC.
                ! Retain ordering structure and reuse keys, not a second
                ! numeric matrix or a per-entry numeric permutation map.
                deallocate (self%ordered%val, self%source)
            end if
            call numeric_factor(self, A, status)
        else if (allocated(self%perm)) then
            call numeric_factor(self, self%ordered, status)
        else
            call numeric_factor(self, A, status)
        end if
        if (status_ok(status)) then
            if (allocated(self%perm)) self%stats%factor_mem_mb = &
                self%stats%factor_mem_mb + real( &
                32_i8*A%nnz + 12_i8*self%n + 8_i8, dp)/1.0e6_dp
            self%stats%analysis_s = wall_time() - start - self%stats%factor_s
        end if
    end subroutine native_factor

    subroutine prepare_ordering(self, A, status)
        class(native_backend_t), intent(inout) :: self
        type(csc_t), intent(in) :: A
        type(fortsparse_status_t), intent(out) :: status
        integer :: bandwidth, j, i, p, q, ierr, code, chosen, old
        integer, allocatable :: cursor(:)

        call status_set(status, FORTSPARSE_OK, '')
        bandwidth = 0
        do j = 1, A%ncol
            do p = A%col_ptr(j), A%col_ptr(j + 1) - 1
                bandwidth = max(bandwidth, abs(A%row_idx(p) - j))
            end do
        end do
        chosen = self%options%ordering
        if (chosen == FORTSPARSE_ORDER_DEFAULT) then
            if (.not. self%reference_only .and. A%ncol >= 4096 .and. &
                bandwidth > 12) then
                chosen = FORTSPARSE_ORDER_AMD
            end if
        end if
        self%selected_order = FORTSPARSE_ORDER_DEFAULT
        if (chosen == FORTSPARSE_ORDER_DEFAULT) return
        allocate (self%perm(A%ncol), self%inverse(A%ncol), stat=ierr)
        if (ierr /= 0) then
            call memory_error(status)
            return
        end if
        call compute_order(chosen, code)
        if (code == FORTSPARSE_BACKEND_UNAVAILABLE .and. &
            self%options%ordering == FORTSPARSE_ORDER_DEFAULT) then
            if (chosen == FORTSPARSE_ORDER_AMD) then
                chosen = FORTSPARSE_ORDER_ND
            else
                chosen = FORTSPARSE_ORDER_AMD
            end if
            call compute_order(chosen, code)
        end if
        if (code /= FORTSPARSE_OK) then
            deallocate (self%perm, self%inverse)
            if (code == FORTSPARSE_BACKEND_UNAVAILABLE .and. &
                self%options%ordering == FORTSPARSE_ORDER_DEFAULT) return
            call status_set(status, code, 'native SPD: ordering failed')
            return
        end if
        self%selected_order = chosen
        allocate (self%source(A%nnz), self%op(A%ncol + 1), self%oi(A%nnz), &
            self%ordered%col_ptr(A%ncol + 1), self%ordered%row_idx(A%nnz), &
            self%ordered%val(A%nnz), cursor(A%ncol), stat=ierr)
        if (ierr /= 0) then
            call memory_error(status)
            return
        end if
        self%op = A%col_ptr
        self%oi = A%row_idx
        self%ordered%nrow = A%nrow
        self%ordered%ncol = A%ncol
        self%ordered%nnz = A%nnz
        self%ordered%col_ptr(1) = 1
        do j = 1, A%ncol
            old = self%perm(j)
            self%ordered%col_ptr(j + 1) = self%ordered%col_ptr(j) + &
                A%col_ptr(old + 1) - A%col_ptr(old)
        end do
        cursor = self%ordered%col_ptr(:A%ncol)
        ! A is full and exactly symmetric. Walk original columns in new ROW
        ! order; scattering their transposes creates sorted permuted CSC in
        ! O(nnz+n), without triplets, sorting or a duplicate numeric layout.
        do i = 1, A%ncol
            old = self%perm(i)
            do p = A%col_ptr(old), A%col_ptr(old + 1) - 1
                j = self%inverse(A%row_idx(p))
                q = cursor(j)
                self%ordered%row_idx(q) = i
                self%ordered%val(q) = A%val(p)
                self%source(q) = p
                cursor(j) = q + 1
            end do
        end do
    contains
        subroutine compute_order(ordering, result_code)
            integer, intent(in) :: ordering
            integer, intent(out) :: result_code
            if (ordering == FORTSPARSE_ORDER_AMD) then
                result_code = native_order_amd(A%ncol, A%nnz, A%col_ptr, &
                    A%row_idx, self%perm, self%inverse)
            else
                result_code = native_order(A%ncol, A%nnz, A%col_ptr, &
                    A%row_idx, self%perm, self%inverse)
            end if
        end subroutine compute_order
    end subroutine prepare_ordering

    subroutine analyze(self, A, status)
        class(native_backend_t), intent(inout) :: self
        type(csc_t), intent(in) :: A
        type(fortsparse_status_t), intent(out) :: status
        integer :: j, i, p, k, q, top, ierr
        integer(i8) :: nz
        real(dp) :: start

        start = wall_time()
        self%n = A%ncol
        allocate (self%parent(self%n), self%lp(self%n + 1), self%count(self%n), &
            self%mark(self%n), self%reach(self%n), self%path(self%n), &
            self%work(self%n), stat=ierr)
        if (ierr /= 0) then
            call self%free()
            call memory_error(status)
            return
        end if
        if (allocated(self%perm)) then
            ! The retained original pattern and permuted CSC already own
            ! both reuse keys; a third pattern copy has no consumer.
            allocate (self%ap(0), self%ai(0), stat=ierr)
        else
            allocate (self%ap(self%n + 1), self%ai(A%nnz), stat=ierr)
        end if
        if (ierr /= 0) then
            call memory_error(status)
            return
        end if
        if (.not. allocated(self%perm)) then
            self%ap = A%col_ptr
            self%ai = A%row_idx
        end if
        self%parent = 0
        self%mark = 0
        self%count = 0
        if (self%n >= 4096 .and. allocated(self%perm)) then
            ! Dense filled paths need not be enumerated just to count them.
            ! Keep the simple low-fill/small setup schedule below.
            call native_symbolic_counts(self%n, A%col_ptr, A%row_idx, &
                self%parent, self%count, status)
            if (.not. status_ok(status)) return
        else
            ! Every previously unseen ancestor receives row j in its L column.
            do j = 1, self%n
                do p = A%col_ptr(j), A%col_ptr(j + 1) - 1
                    i = A%row_idx(p)
                    if (i >= j) exit
                    do while (i /= 0)
                        if (self%mark(i) == j) exit
                        self%mark(i) = j
                        self%count(i) = self%count(i) + 1
                        if (self%parent(i) == 0) self%parent(i) = j
                        i = self%parent(i)
                        if (i == j) exit
                    end do
                end do
            end do
        end if
        nz = int(self%n, i8) + sum(int(self%count, i8))
        if (nz >= int(huge(1), i8)) then
            call self%free()
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                'native SPD: factor exceeds supported index range')
            return
        end if
        self%lp(1) = 1
        do j = 1, self%n
            self%lp(j + 1) = self%lp(j) + self%count(j) + 1
        end do
        if (nz > int(huge(0_c_size_t)/c_sizeof(0_c_int32_t), i8)) then
            call self%free()
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                'native SPD: symbolic storage exceeds supported byte range')
            return
        end if
        self%li_bytes = int(nz, c_size_t)*c_sizeof(0_c_int32_t)
        self%li_owner = symbolic_alloc(self%li_bytes)
        if (.not. c_associated(self%li_owner)) then
            call self%free()
            call memory_error(status)
            return
        end if
        call c_f_pointer(self%li_owner, self%li, [int(nz)])
        self%mark = 0
        self%count = 1
        do j = 1, self%n
            self%li(self%lp(j)) = j
            call row_reach(self, A, j, top)
            do q = top, self%n
                i = self%reach(q)
                k = self%lp(i) + self%count(i)
                self%li(k) = j
                self%count(i) = self%count(i) + 1
            end do
        end do
        self%blocked = .not. self%reference_only .and. &
            ((self%n >= 64 .and. maxval(self%count) > 12) .or. &
            self%options%device >= 0)
        if (self%blocked) then
            call self%panels%build(self%n, self%lp, self%li, &
                self%options%nthreads, status, self%options%device)
            if (.not. status_ok(status)) return
        else
            allocate (self%lv(int(nz)), stat=ierr)
            if (ierr /= 0) then
                call memory_error(status)
                return
            end if
        end if
        self%stats%ordering_used = self%selected_order
        self%stats%symbolic_nnz = nz
        self%stats%factor_nnz = nz
        self%stats%factor_mem_mb = real( &
            4_i8*(size(self%ap, kind=i8) + size(self%ai, kind=i8) + &
            size(self%lp, kind=i8) + size(self%li, kind=i8) + &
            5_i8*self%n) + 8_i8*(nz + self%n), dp)/1.0e6_dp
        if (self%blocked) then
            self%stats%factor_nnz = self%panels%factor_nnz
            self%stats%factor_mem_mb = self%stats%factor_mem_mb + &
                real(self%panels%memory_bytes() - 8_i8*nz, dp)/1.0e6_dp
        end if
        if (self%blocked) then
            self%stats%factor_storage_bytes = self%panels%factor_storage_bytes
            self%stats%workspace_bytes = self%panels%workspace_bytes + &
                4_i8*(size(self%ap, kind=i8) + size(self%ai, kind=i8)) + &
                8_i8*self%n
            ! Panel row/index structures now own the analysis result. Release
            ! the temporary simplicial structure before touching numeric pages.
            call release_symbolic_rows(self)
            deallocate (self%lp, self%parent, self%count, self%mark, &
                self%reach, self%path)
        else
            self%stats%factor_storage_bytes = 8_i8*nz
            self%stats%workspace_bytes = 4_i8*nz + 4_i8*(self%n + 1) + 4_i8*( &
                size(self%ap, kind=i8) + size(self%ai, kind=i8) + &
                5_i8*self%n) + 8_i8*self%n
        end if
        if (allocated(self%perm)) then
            if (self%blocked) then
                self%stats%workspace_bytes = self%stats%workspace_bytes &
                    + 8_i8*A%nnz + 16_i8*self%n + 8_i8
            else
                self%stats%workspace_bytes = self%stats%workspace_bytes &
                    + 20_i8*A%nnz + 16_i8*self%n + 8_i8
            end if
        end if
        self%stats%analysis_s = wall_time() - start
        self%stats%nthreads = 1
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine analyze

    subroutine validate(self, A, status)
        class(native_backend_t), intent(in) :: self
        type(csc_t), intent(in) :: A
        type(fortsparse_status_t), intent(out) :: status
        integer :: j, p, i, lo, hi, mid, previous

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            'native SPD: require finite, full symmetric, canonical square CSC')
        if (self%options%device >= 0 .and. &
            (.not. native_cuda_available .or. self%reference_only)) then
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                'native SPD: CUDA sparse factorization not enabled')
            return
        end if
        if (self%options%symmetry /= FORTSPARSE_SYM_SPD) then
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                'native backend supports real SPD matrices only')
            return
        end if
        if (self%options%ordering /= FORTSPARSE_ORDER_DEFAULT .and. &
            self%options%ordering /= FORTSPARSE_ORDER_ND .and. &
            self%options%ordering /= FORTSPARSE_ORDER_AMD) then
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                'native SPD: unsupported ordering option')
            return
        end if
        if (.not. allocated(A%col_ptr)) return
        if (.not. allocated(A%row_idx)) return
        if (.not. allocated(A%val)) return
        if (.not. csc_is_valid(A)) return
        if (A%nrow /= A%ncol .or. size(A%val) /= A%nnz) return
        if (.not. all(ieee_is_finite(A%val))) return
        do j = 1, A%ncol
            previous = 0
            do p = A%col_ptr(j), A%col_ptr(j + 1) - 1
                i = A%row_idx(p)
                if (i <= previous) return
                previous = i
                lo = A%col_ptr(i)
                hi = A%col_ptr(i + 1) - 1
                do while (lo <= hi)
                    mid = lo + (hi - lo)/2
                    if (A%row_idx(mid) < j) then
                        lo = mid + 1
                    else if (A%row_idx(mid) > j) then
                        hi = mid - 1
                    else
                        exit
                    end if
                end do
                if (lo > hi) return
                if (A%val(p) /= A%val(mid)) return
            end do
        end do
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine validate

    subroutine row_reach(self, A, j, top)
        class(native_backend_t), intent(inout) :: self
        type(csc_t), intent(in) :: A
        integer, intent(in) :: j
        integer, intent(out) :: top
        integer :: p, i, length

        top = self%n + 1
        self%mark(j) = j
        do p = A%col_ptr(j), A%col_ptr(j + 1) - 1
            i = A%row_idx(p)
            if (i >= j) exit
            length = 0
            do while (self%mark(i) /= j)
                length = length + 1
                self%path(length) = i
                self%mark(i) = j
                i = self%parent(i)
            end do
            do while (length > 0)
                top = top - 1
                self%reach(top) = self%path(length)
                length = length - 1
            end do
        end do
    end subroutine row_reach

    subroutine numeric_factor(self, A, status)
        class(native_backend_t), intent(inout) :: self
        type(csc_t), intent(in) :: A
        type(fortsparse_status_t), intent(out) :: status
        integer :: j, p, q, i, top
        real(dp) :: d, v, start, flops

        start = wall_time()
        self%ready = .false.
        self%stats%device_bytes = 0_i8
        self%stats%active_device = -1
        if (self%blocked) then
            if (allocated(self%perm)) then
                call self%panels%factor(A, status, self%perm, self%inverse)
            else
                call self%panels%factor(A, status)
            end if
            self%stats%factor_s = self%panels%numeric_s
            self%stats%flops = self%panels%flops
            self%stats%nthreads = self%options%nthreads
            self%stats%active_threads = self%panels%active_threads
            self%stats%device_bytes = self%panels%device_bytes
            self%stats%active_device = self%options%device
            self%ready = status_ok(status)
            return
        end if
        self%stats%active_threads = 1
        self%work = 0.0_dp
        self%mark = 0
        self%count = 1
        flops = 0.0_dp
        do j = 1, self%n
            call row_reach(self, A, j, top)
            do p = A%col_ptr(j), A%col_ptr(j + 1) - 1
                i = A%row_idx(p)
                if (i > j) exit
                self%work(i) = A%val(p)
            end do
            d = self%work(j)
            self%work(j) = 0.0_dp
            do q = top, self%n
                i = self%reach(q)
                v = spd_divide(self%work(i), self%lv(self%lp(i)))
                self%work(i) = 0.0_dp
                p = self%lp(i) + 1
                call spd_scatter_update(self%work, self%li(p:), self%lv(p:), &
                    self%count(i) - 1, v)
                d = spd_update(d, v, v)
                self%lv(self%lp(i) + self%count(i)) = v
                flops = flops + 2.0_dp*self%count(i)
                self%count(i) = self%count(i) + 1
            end do
            if (.not. ieee_is_finite(d) .or. d <= 0.0_dp) then
                call status_set(status, FORTSPARSE_NOT_POSITIVE_DEFINITE, &
                    'native SPD: nonpositive or nonfinite Cholesky pivot')
                return
            end if
            self%lv(self%lp(j)) = spd_root(d)
        end do
        self%stats%factor_s = wall_time() - start
        self%stats%flops = flops + self%n
        self%ready = .true.
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine numeric_factor

    subroutine native_refactor(self, A, refine, status, reused)
        class(native_backend_t), intent(inout) :: self
        type(csc_t), intent(in) :: A
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(out) :: reused
        integer :: p

        reused = .false.
        call validate(self, A, status)
        self%ready = .false.
        if (.not. status_ok(status)) return
        if (self%analyzed_order /= self%options%ordering .or. &
            self%analyzed_threads /= self%options%nthreads .or. &
            self%analyzed_device /= self%options%device) then
            call self%factor_real(A, refine, status)
            return
        end if
        if (allocated(self%op)) then
            if (self%n == A%ncol .and. size(self%oi) == A%nnz) then
                if (all(self%op == A%col_ptr) .and. &
                    all(self%oi == A%row_idx)) then
                    reused = .true.
                    self%stats%analysis_s = 0.0_dp
                    if (self%blocked) then
                        call numeric_factor(self, A, status)
                    else
                        do p = 1, A%nnz
                            self%ordered%val(p) = A%val(self%source(p))
                        end do
                        call numeric_factor(self, self%ordered, status)
                    end if
                    return
                end if
            end if
            call self%factor_real(A, refine, status)
            return
        end if
        if (allocated(self%ap)) then
            if (self%n == A%ncol .and. size(self%ai) == A%nnz) then
                if (all(self%ap == A%col_ptr) .and. &
                    all(self%ai == A%row_idx)) then
                    reused = .true.
                    self%stats%analysis_s = 0.0_dp
                    call numeric_factor(self, A, status)
                    return
                end if
            end if
        end if
        call self%factor_real(A, refine, status)
    end subroutine native_refactor

    subroutine native_solve(self, b, x, status)
        class(native_backend_t), intent(inout) :: self
        real(dp), target, contiguous, intent(in) :: b(:)
        real(dp), target, contiguous, intent(out) :: x(:)
        type(fortsparse_status_t), intent(out) :: status

        call solve_shape(self, size(b), size(x), status)
        if (.not. status_ok(status)) return
        x = b
        call self%solve_real_inplace(x, status)
    end subroutine native_solve

    subroutine native_solve_inplace(self, b, status)
        class(native_backend_t), intent(inout) :: self
        real(dp), intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out) :: status

        call solve_shape(self, size(b), size(b), status)
        if (.not. status_ok(status)) return
        if (allocated(self%perm)) then
            self%work = b(self%perm)
            if (self%blocked) then
                call self%panels%solve(self%work, status)
            else
                call triangular_solve(self, self%work)
            end if
            b(self%perm) = self%work
        else if (self%blocked) then
            self%work = b
            call self%panels%solve(self%work, status)
            b = self%work
        else
            call triangular_solve(self, b)
        end if
    end subroutine native_solve_inplace

    subroutine triangular_solve(self, b)
        class(native_backend_t), intent(in) :: self
        real(dp), intent(inout) :: b(:)
        integer :: j, p
        real(dp) :: v

        do j = 1, self%n
            v = spd_divide(b(j), self%lv(self%lp(j)))
            b(j) = v
            p = self%lp(j) + 1
            call spd_scatter_update(b, self%li(p:), self%lv(p:), &
                self%lp(j + 1) - p, v)
        end do
        do j = self%n, 1, -1
            v = b(j)
            p = self%lp(j) + 1
            call spd_gather_update(v, self%li(p:), self%lv(p:), b, &
                self%lp(j + 1) - p)
            b(j) = spd_divide(v, self%lv(self%lp(j)))
        end do
    end subroutine triangular_solve

    subroutine native_solve_block(self, B, X, status)
        class(native_backend_t), intent(inout) :: self
        real(dp), contiguous, intent(in) :: B(:, :)
        real(dp), contiguous, intent(out) :: X(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: j

        call solve_shape(self, size(B, 1), size(X, 1), status)
        if (.not. status_ok(status)) return
        if (size(B, 2) /= size(X, 2)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'native SPD: incompatible RHS count')
            return
        end if
        X = B
        do j = 1, size(X, 2)
            call self%solve_real_inplace(X(:, j), status)
        end do
    end subroutine native_solve_block

    subroutine solve_shape(self, nb, nx, status)
        class(native_backend_t), intent(in) :: self
        integer, intent(in) :: nb, nx
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_OK, '')
        if (.not. self%ready) then
            call status_set(status, FORTSPARSE_NOT_FACTORED, &
                'native SPD: no successful factorization')
        else if (nb /= self%n .or. nx /= self%n) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'native SPD: incompatible RHS dimensions')
        end if
    end subroutine solve_shape

    subroutine native_stats(self, stats, status)
        class(native_backend_t), intent(inout) :: self
        type(sparse_factor_stats_t), intent(out) :: stats
        type(fortsparse_status_t), intent(out) :: status

        stats = self%stats
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine native_stats

    function native_vector(self, n) result(p)
        class(native_backend_t), intent(inout) :: self
        integer, intent(in) :: n
        real(dp), pointer :: p(:)

        if (associated(self%user_vector)) deallocate (self%user_vector)
        allocate (self%user_vector(n))
        p => self%user_vector
    end function native_vector

    subroutine release_symbolic_rows(self)
        class(native_backend_t), intent(inout) :: self

        if (c_associated(self%li_owner)) &
            call symbolic_free(self%li_owner, self%li_bytes)
        nullify (self%li)
        self%li_owner = c_null_ptr
        self%li_bytes = 0_c_size_t
    end subroutine release_symbolic_rows

    subroutine native_free(self)
        class(native_backend_t), intent(inout) :: self

        call self%panels%free()
        self%blocked = .false.
        if (allocated(self%perm)) deallocate (self%perm)
        if (allocated(self%inverse)) deallocate (self%inverse)
        if (allocated(self%source)) deallocate (self%source)
        if (allocated(self%op)) deallocate (self%op)
        if (allocated(self%oi)) deallocate (self%oi)
        self%ordered = csc_t()
        if (allocated(self%ap)) deallocate (self%ap)
        if (allocated(self%ai)) deallocate (self%ai)
        if (allocated(self%parent)) deallocate (self%parent)
        if (allocated(self%lp)) deallocate (self%lp)
        call release_symbolic_rows(self)
        if (allocated(self%lv)) deallocate (self%lv)
        if (allocated(self%count)) deallocate (self%count)
        if (allocated(self%mark)) deallocate (self%mark)
        if (allocated(self%reach)) deallocate (self%reach)
        if (allocated(self%path)) deallocate (self%path)
        if (allocated(self%work)) deallocate (self%work)
        if (associated(self%user_vector)) deallocate (self%user_vector)
        nullify (self%user_vector)
        self%n = 0
        self%ready = .false.
        self%stats = sparse_factor_stats_t()
    end subroutine native_free

    subroutine native_finalize(self)
        type(native_backend_t), intent(inout) :: self
        call self%free()
    end subroutine native_finalize

    subroutine native_factor_complex(self, A, refine, status)
        class(native_backend_t), intent(inout) :: self
        type(csc_z_t), intent(in) :: A
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        call self%free()
        call unsupported(status)
    end subroutine native_factor_complex

    subroutine native_solve_complex(self, b, x, status)
        class(native_backend_t), intent(inout) :: self
        complex(dp), intent(in) :: b(:)
        complex(dp), intent(out) :: x(:)
        type(fortsparse_status_t), intent(out) :: status
        call unsupported(status)
    end subroutine native_solve_complex

    subroutine native_solve_complex_inplace(self, b, status)
        class(native_backend_t), intent(inout) :: self
        complex(dp), intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out) :: status
        call unsupported(status)
    end subroutine native_solve_complex_inplace

    subroutine unsupported(status)
        type(fortsparse_status_t), intent(out) :: status
        call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
            'native backend supports real SPD matrices only')
    end subroutine unsupported

    subroutine memory_error(status)
        type(fortsparse_status_t), intent(out) :: status
        call status_set(status, FORTSPARSE_OUT_OF_MEMORY, &
            'native SPD: workspace allocation failed')
    end subroutine memory_error

    real(dp) function wall_time() result(t)
        integer(i8) :: ticks, rate
        call system_clock(ticks, rate)
        t = real(ticks, dp)/real(rate, dp)
    end function wall_time
end module fortsparse_native
