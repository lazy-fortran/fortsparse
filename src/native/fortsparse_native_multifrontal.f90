module fortsparse_native_multifrontal
    ! Real symmetric indefinite multifrontal LDL^T, P^T A P = L D L^T, with
    ! threshold 1x1/2x2 pivots among the fully summed variables of each front
    ! and delayed pivots: a variable that fails the test moves to the parent
    ! front (Duff and Reid, ACM TOMS 9 (1983) 302; MA57, MUMPS). Variables
    ! whose diagonal fails the threshold are ordered next to a matched
    ! neighbour, so most of them meet their 2x2 partner in the same front
    ! (Duff and Pralet, SIAM J. Matrix Anal. Appl. 27 (2005) 313).
    use, intrinsic :: iso_c_binding, only: c_int
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csc, only: csc_t, csc_z_t, csc_is_valid
    use fortsparse_status, only: fortsparse_status_t, status_set, status_ok, &
        FORTSPARSE_OK, FORTSPARSE_SINGULAR, FORTSPARSE_INVALID_MATRIX, &
        FORTSPARSE_NOT_FACTORED, FORTSPARSE_BACKEND_UNAVAILABLE, &
        FORTSPARSE_INTERNAL_ERROR
    use fortsparse_backend, only: sparse_backend_t, sparse_factor_stats_t, &
        FORTSPARSE_ORDER_DEFAULT, FORTSPARSE_ORDER_AMD, FORTSPARSE_ORDER_ND, &
        FORTSPARSE_ORDER_NATURAL, FORTSPARSE_OP_NORMAL, FORTSPARSE_OP_TRANSPOSE, &
        FORTSPARSE_OP_ADJOINT
    use fortsparse_inertia, only: inertia_count_t, inertia_add_pivot, &
        inertia_add_scaled_block, inertia_status_of, inertia_tolerance
    use fortsparse_native_residual, only: native_original_norms, native_residual_error
    use fortnum_blas, only: dgemm => fortnum_dgemm
    use fortnum_blas12, only: dgemv => fortnum_dgemv, dtrsv => fortnum_dtrsv
    implicit none
    private
    public :: native_multifrontal_backend_t

    real(dp), parameter :: bk_alpha = (1.0_dp + sqrt(17.0_dp))/8.0_dp
    ! Pivot threshold u: every multiplier is bounded by 1/u (MUMPS default).
    real(dp), parameter :: default_threshold = 0.01_dp
    ! Columns searched for pivots before a BLAS-3 update of the front.
    integer, parameter :: window = 32

    type :: multifrontal_ldlt_t
        integer :: n = 0, nsn = 0
        ! FORTSPARSE_ORDER_* actually used (DEFAULT: no fill-reducing order)
        integer :: order_used = FORTSPARSE_ORDER_DEFAULT
        real(dp) :: u = default_threshold
        ! perm(new) = old, inverse(old) = new: postordered fill-reducing order
        integer, allocatable :: perm(:), inverse(:)
        ! Supernode s owns the new columns first(s):first(s+1)-1; its rows
        ! below are srow(sptr(s):sptr(s+1)-1), ascending; sparent(s) = 0 at roots.
        integer, allocatable :: first(:), sparent(:), sptr(:), srow(:), nchild(:)
        ! Lower original entries of new column j: rows arow, value index asrc.
        integer, allocatable :: aptr(:), arow(:), asrc(:)
        ! Index of the transposed entry of every stored entry (0: not stored).
        integer, allocatable :: tpos(:)
        integer(i8) :: estimate = 0, lnnz = 0
        ! Front s: m rows fidx(fidx_ptr(s):), npiv pivots, L(m, npiv) at lptr(s).
        integer, allocatable :: fm(:), fnpiv(:), fidx_ptr(:), fidx(:)
        integer(i8), allocatable :: lptr(:)
        real(dp), allocatable :: l(:)
        ! Pivot kind by new index: 1 (1x1), 2 (first of a 2x2), 3 (second).
        ! d holds the diagonal, e the off-diagonal of a 2x2 block.
        integer, allocatable :: pk(:)
        real(dp), allocatable :: d(:), e(:)
        ! Inertia thresholds per variable: rowmax = max |a_ij| of its row,
        ! hist = sum of l_ik^2 |d_k| over the pivots of earlier fronts (a
        ! 2x2 pivot contributes |l|^T (|D| + |e| I) |l|), which bounds every
        ! Schur entry that the row received.
        real(dp), allocatable :: rowmax(:), hist(:)
        type(inertia_count_t) :: inertia
        integer :: delayed = 0, two_by_two = 0, max_front = 0
        integer(i8) :: flops = 0
        ! workspace
        integer, allocatable :: map(:), pos(:), istack(:)
        real(dp), allocatable :: stack(:), work(:), y(:)
    end type multifrontal_ldlt_t

    type, extends(sparse_backend_t) :: native_multifrontal_backend_t
        type(multifrontal_ldlt_t) :: f
        logical :: ready = .false., refining = .false., counted = .false.
        integer :: analysed_order = -1
        integer, allocatable :: ap(:), ai(:)
        real(dp), allocatable :: original(:), saved(:), residual(:), scratch(:)
        real(dp) :: matrix_norm = 0.0_dp
        type(sparse_factor_stats_t) :: stats
        real(dp), pointer :: user_vector(:) => null()
    contains
        procedure :: factor_real => mf_factor_real
        procedure :: refactor_real => mf_refactor_real
        procedure :: factor_complex => mf_factor_complex
        procedure :: solve_real => mf_solve_real
        procedure :: solve_real_op => mf_solve_real_op
        procedure :: solve_real_inplace => mf_solve_real_inplace
        procedure :: solve_real_block_op => mf_solve_real_block_op
        procedure :: solve_complex => mf_solve_complex
        procedure :: solve_complex_inplace => mf_solve_complex_inplace
        procedure :: factor_stats => mf_factor_stats
        procedure :: copy_from => mf_copy_from
        procedure :: vector => mf_vector
        procedure :: free => mf_free
        final :: mf_finalize
    end type native_multifrontal_backend_t

    interface
        integer(c_int) function order_amd(n, nnz, ptr, row, perm, inverse) &
                bind(C, name='fsparse_native_order_amd_general')
            import :: c_int
            integer(c_int), value :: n, nnz
            integer(c_int), intent(in) :: ptr(*), row(*)
            integer(c_int), intent(out) :: perm(*), inverse(*)
        end function order_amd
        integer(c_int) function order_nd(n, nnz, ptr, row, perm, inverse) &
                bind(C, name='fsparse_native_order_general')
            import :: c_int
            integer(c_int), value :: n, nnz
            integer(c_int), intent(in) :: ptr(*), row(*)
            integer(c_int), intent(out) :: perm(*), inverse(*)
        end function order_nd
    end interface

contains

    ! ---------------------------------------------------------------- backend

    subroutine mf_factor_real(self, A, refine, status)
        class(native_multifrontal_backend_t), intent(inout) :: self
        type(csc_t), intent(in) :: A
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        call factor_values(self, A, refine, .false., status)
    end subroutine mf_factor_real

    subroutine mf_refactor_real(self, A, refine, status, reused)
        class(native_multifrontal_backend_t), intent(inout) :: self
        type(csc_t), intent(in) :: A
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(out) :: reused
        reused = same_pattern(self, A)
        call factor_values(self, A, refine, reused, status)
    end subroutine mf_refactor_real

    logical function same_pattern(self, A) result(same)
        class(native_multifrontal_backend_t), intent(in) :: self
        type(csc_t), intent(in) :: A
        same = .false.
        if (.not. allocated(self%ap) .or. .not. allocated(A%col_ptr)) return
        if (self%analysed_order /= self%options%ordering) return
        if (A%ncol /= self%f%n .or. size(A%col_ptr) /= size(self%ap)) return
        if (size(A%row_idx) < A%nnz .or. size(self%ai) /= A%nnz) return
        if (any(A%col_ptr /= self%ap)) return
        same = all(A%row_idx(1:A%nnz) == self%ai)
    end function same_pattern

    subroutine factor_values(self, A, refine, reuse, status)
        class(native_multifrontal_backend_t), intent(inout) :: self
        type(csc_t), intent(in) :: A
        logical, intent(in) :: refine, reuse
        type(fortsparse_status_t), intent(out) :: status
        integer(i8) :: t0, t1, t2, rate
        integer :: info
        real(dp) :: tnorm

        self%ready = .false.
        self%counted = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            'multifrontal LDL^T: finite, full canonical symmetric square CSC required')
        if (.not. canonical(A)) return
        if (reuse) then
            if (.not. same_values_symmetric(self%f%tpos, A%val)) return
        end if
        if (self%options%device >= 0) then
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                'multifrontal LDL^T: CPU only')
            return
        end if
        call system_clock(t0, rate)
        if (.not. reuse) then
            call analyse(self%f, A%ncol, A%col_ptr, A%row_idx, A%val, &
                self%options%ordering, info)
            if (info < 0) return
            if (info /= 0) then
                call status_set(status, FORTSPARSE_INTERNAL_ERROR, &
                    'multifrontal LDL^T: ordering failed')
                return
            end if
            self%ap = A%col_ptr
            self%ai = A%row_idx(1:A%nnz)
            self%analysed_order = self%options%ordering
        end if
        call system_clock(t1)
        call factor(self%f, A%val)
        call system_clock(t2)
        self%refining = refine
        if (refine) then
            self%original = A%val(1:A%nnz)
            if (allocated(self%saved)) deallocate (self%saved, self%residual, self%scratch)
            allocate (self%saved(A%ncol), self%residual(A%ncol), self%scratch(A%ncol))
            call native_original_norms(A%ncol, self%ap, self%ai, self%original, &
                self%scratch, self%matrix_norm, tnorm, status)
            if (.not. status_ok(status)) return
        else if (allocated(self%original)) then
            deallocate (self%original)
        end if
        call fill_stats(self, real(t1 - t0, dp)/real(rate, dp), &
            real(t2 - t1, dp)/real(rate, dp), reuse)
        self%counted = .true.
        if (self%f%inertia%zero > 0_i8) then
            ! the inertia stays available; solves are refused
            call status_set(status, FORTSPARSE_SINGULAR, &
                'multifrontal LDL^T: exactly zero pivot, matrix is singular')
            return
        end if
        self%ready = .true.
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine factor_values

    logical function canonical(A) result(ok)
        ! square CSC with strictly ascending rows and finite values
        type(csc_t), intent(in) :: A
        integer :: j, p
        ok = .false.
        if (A%nrow /= A%ncol .or. A%ncol < 1) return
        if (.not. allocated(A%col_ptr) .or. .not. allocated(A%row_idx) .or. &
            .not. allocated(A%val)) return
        if (.not. csc_is_valid(A)) return
        if (size(A%val) < A%nnz) return
        do j = 1, A%ncol
            do p = A%col_ptr(j) + 1, A%col_ptr(j + 1) - 1
                if (A%row_idx(p) <= A%row_idx(p - 1)) return
            end do
        end do
        ok = all(ieee_is_finite(A%val(1:A%nnz)))
    end function canonical

    pure logical function same_values_symmetric(tpos, val) result(ok)
        integer, intent(in) :: tpos(:)
        real(dp), intent(in) :: val(:)
        integer :: p
        ok = .false.
        do p = 1, size(tpos)
            if (tpos(p) > 0) then
                if (val(p) /= val(tpos(p))) return
            else if (val(p) /= 0.0_dp) then
                return
            end if
        end do
        ok = .true.
    end function same_values_symmetric

    subroutine fill_stats(self, analysis_s, factor_s, reuse)
        class(native_multifrontal_backend_t), intent(inout) :: self
        real(dp), intent(in) :: analysis_s, factor_s
        logical, intent(in) :: reuse
        integer(i8) :: nnz
        integer :: s
        type(sparse_factor_stats_t) :: fresh
        if (reuse) fresh%analysis_s = self%stats%analysis_s
        if (.not. reuse) fresh%analysis_s = analysis_s
        self%stats = fresh
        nnz = 0
        do s = 1, self%f%nsn
            nnz = nnz + int(self%f%fm(s), i8)*self%f%fnpiv(s) - &
                int(self%f%fnpiv(s), i8)*(self%f%fnpiv(s) - 1)/2
        end do
        self%stats%factor_nnz = nnz
        self%stats%symbolic_nnz = self%f%lnnz
        self%stats%factor_s = factor_s
        self%stats%flops = real(self%f%flops, dp)
        self%stats%factor_storage_bytes = 8_i8*self%f%lptr(self%f%nsn + 1)
        self%stats%factor_mem_mb = real(8_i8*(size(self%f%l, kind=i8) + &
            size(self%f%stack, kind=i8)), dp)/1.0e6_dp
        self%stats%active_threads = 1
        self%stats%nthreads = 1
        self%stats%active_device = -1
        self%stats%ordering_used = self%f%order_used
        self%stats%native_factor_family = 'LDLT_MULTIFRONT'
        self%stats%perturbed_pivots = 0
        self%stats%block_pivot_swaps = self%f%two_by_two
        self%stats%inertia_negative = self%f%inertia%negative
        self%stats%inertia_zero = self%f%inertia%zero
        self%stats%inertia_positive = self%f%inertia%positive
        self%stats%inertia_uncertain = self%f%inertia%uncertain
        self%stats%inertia_status = inertia_status_of(self%f%inertia)
        self%stats%inertia_min_pivot = self%f%inertia%min_magnitude
        self%stats%inertia_tolerance = self%f%inertia%tolerance
        self%stats%min_pivot = self%f%inertia%min_magnitude
        self%stats%numeric_started = .true.
        ! retained payload; the AMD and BLAS workspaces are opaque, so the
        ! peak is unknown
        self%stats%host_current_bytes = payload_bytes(self)
        self%stats%host_peak_bytes = -1_i8
        self%stats%strict_memory_supported = .false.
    end subroutine fill_stats

    integer(i8) function payload_bytes(self) result(bytes)
        class(native_multifrontal_backend_t), intent(in) :: self
        bytes = 8_i8*(sz_r(self%f%l) + sz_r(self%f%stack) + sz_r(self%f%work) + &
            sz_r(self%f%y) + sz_r(self%f%d) + sz_r(self%f%e) + sz_r(self%f%rowmax) + &
            sz_r(self%f%hist) + sz_r(self%original) + sz_r(self%saved) + &
            sz_r(self%residual) + sz_r(self%scratch)) + &
            4_i8*(sz_i(self%f%perm) + sz_i(self%f%inverse) + sz_i(self%f%first) + &
            sz_i(self%f%sparent) + sz_i(self%f%sptr) + sz_i(self%f%srow) + &
            sz_i(self%f%nchild) + sz_i(self%f%aptr) + sz_i(self%f%arow) + &
            sz_i(self%f%asrc) + sz_i(self%f%tpos) + sz_i(self%f%fm) + &
            sz_i(self%f%fnpiv) + sz_i(self%f%fidx_ptr) + sz_i(self%f%fidx) + &
            sz_i(self%f%pk) + sz_i(self%f%map) + sz_i(self%f%pos) + &
            sz_i(self%f%istack) + sz_i(self%ap) + sz_i(self%ai)) + &
            8_i8*sz_i8(self%f%lptr)
    contains
        integer(i8) function sz_r(a)
            real(dp), allocatable, intent(in) :: a(:)
            sz_r = 0
            if (allocated(a)) sz_r = size(a, kind=i8)
        end function sz_r
        integer(i8) function sz_i(a)
            integer, allocatable, intent(in) :: a(:)
            sz_i = 0
            if (allocated(a)) sz_i = size(a, kind=i8)
        end function sz_i
        integer(i8) function sz_i8(a)
            integer(i8), allocatable, intent(in) :: a(:)
            sz_i8 = 0
            if (allocated(a)) sz_i8 = size(a, kind=i8)
        end function sz_i8
    end function payload_bytes

    subroutine mf_solve_real_inplace(self, b, status)
        class(native_multifrontal_backend_t), intent(inout) :: self
        real(dp), intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: step, i
        real(dp) :: error
        if (.not. solvable(self, size(b), status)) return
        if (self%refining) self%saved = b
        call solve(self%f, b)
        call status_set(status, FORTSPARSE_OK, '')
        if (.not. self%refining) return
        do step = 0, 3
            self%stats%refinement_steps = step
            call native_residual_error(self%ap, self%ai, self%original, self%saved, b, 0, &
                self%matrix_norm, self%residual, error)
            self%stats%backward_error = error
            if (error < 0.0_dp) then
                call status_set(status, FORTSPARSE_SINGULAR, &
                    'multifrontal LDL^T refinement: nonfinite iterate or residual')
                return
            end if
            if (error <= 32.0_dp*epsilon(1.0_dp)) return
            if (step == 3) exit
            call solve(self%f, self%residual)
            do i = 1, size(b)
                b(i) = b(i) + self%residual(i)
            end do
        end do
        call status_set(status, FORTSPARSE_SINGULAR, &
            'multifrontal LDL^T refinement: target not met in three corrections')
    end subroutine mf_solve_real_inplace

    logical function solvable(self, n, status)
        class(native_multifrontal_backend_t), intent(in) :: self
        integer, intent(in) :: n
        type(fortsparse_status_t), intent(out) :: status
        solvable = .false.
        if (.not. self%ready) then
            call status_set(status, FORTSPARSE_NOT_FACTORED, 'multifrontal LDL^T: not factored')
        else if (n /= self%f%n) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'multifrontal LDL^T: right-hand side length mismatch')
        else if (self%f%inertia%zero > 0_i8) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                'multifrontal LDL^T: zero pivot, matrix is singular')
        else
            solvable = .true.
        end if
    end function solvable

    subroutine mf_solve_real(self, b, x, status)
        class(native_multifrontal_backend_t), intent(inout) :: self
        real(dp), target, contiguous, intent(in) :: b(:)
        real(dp), target, contiguous, intent(out) :: x(:)
        type(fortsparse_status_t), intent(out) :: status
        if (size(x) /= size(b)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'multifrontal LDL^T: solution length mismatch')
            return
        end if
        x = b
        call self%solve_real_inplace(x, status)
    end subroutine mf_solve_real

    subroutine mf_solve_real_op(self, b, x, operation, status)
        ! A is symmetric: transpose and adjoint solves are the same solve.
        class(native_multifrontal_backend_t), intent(inout) :: self
        real(dp), target, contiguous, intent(in) :: b(:)
        real(dp), target, contiguous, intent(out) :: x(:)
        integer, intent(in) :: operation
        type(fortsparse_status_t), intent(out) :: status
        if (operation /= FORTSPARSE_OP_NORMAL .and. operation /= FORTSPARSE_OP_TRANSPOSE &
            .and. operation /= FORTSPARSE_OP_ADJOINT) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'multifrontal LDL^T: unknown solve operation')
            return
        end if
        call self%solve_real(b, x, status)
    end subroutine mf_solve_real_op

    subroutine mf_solve_real_block_op(self, b, x, operation, status)
        class(native_multifrontal_backend_t), intent(inout) :: self
        real(dp), contiguous, intent(in) :: b(:, :)
        real(dp), contiguous, intent(out) :: x(:, :)
        integer, intent(in) :: operation
        type(fortsparse_status_t), intent(out) :: status
        integer :: j
        call status_set(status, FORTSPARSE_OK, '')
        if (any(shape(b) /= shape(x))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'multifrontal LDL^T: block shape mismatch')
            return
        end if
        do j = 1, size(b, 2)
            call self%solve_real_op(b(:, j), x(:, j), operation, status)
            if (.not. status_ok(status)) return
        end do
    end subroutine mf_solve_real_block_op

    subroutine mf_factor_complex(self, A, refine, status)
        class(native_multifrontal_backend_t), intent(inout) :: self
        type(csc_z_t), intent(in) :: A
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
            'multifrontal LDL^T: real matrices only')
    end subroutine mf_factor_complex

    subroutine mf_solve_complex(self, b, x, status)
        class(native_multifrontal_backend_t), intent(inout) :: self
        complex(dp), intent(in) :: b(:)
        complex(dp), intent(out) :: x(:)
        type(fortsparse_status_t), intent(out) :: status
        x = 0
        call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
            'multifrontal LDL^T: real right-hand sides only')
    end subroutine mf_solve_complex

    subroutine mf_solve_complex_inplace(self, b, status)
        class(native_multifrontal_backend_t), intent(inout) :: self
        complex(dp), intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out) :: status
        call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
            'multifrontal LDL^T: real right-hand sides only')
    end subroutine mf_solve_complex_inplace

    subroutine mf_factor_stats(self, stats, status)
        class(native_multifrontal_backend_t), intent(inout) :: self
        type(sparse_factor_stats_t), intent(out) :: stats
        type(fortsparse_status_t), intent(out) :: status
        stats = self%stats
        call status_set(status, FORTSPARSE_OK, '')
        if (.not. self%counted) call status_set(status, FORTSPARSE_NOT_FACTORED, &
            'multifrontal LDL^T: not factored')
    end subroutine mf_factor_stats

    subroutine mf_copy_from(self, src, status)
        class(native_multifrontal_backend_t), intent(inout) :: self
        class(sparse_backend_t), intent(in) :: src
        type(fortsparse_status_t), intent(out) :: status
        select type (src)
        type is (native_multifrontal_backend_t)
            self%f = src%f
            self%ready = src%ready
            self%counted = src%counted
            self%refining = src%refining
            self%analysed_order = src%analysed_order
            if (allocated(src%ap)) self%ap = src%ap
            if (allocated(src%ai)) self%ai = src%ai
            if (allocated(src%original)) self%original = src%original
            if (allocated(src%saved)) then
                self%saved = src%saved
                self%residual = src%residual
                self%scratch = src%scratch
            end if
            self%matrix_norm = src%matrix_norm
            self%stats = src%stats
            call status_set(status, FORTSPARSE_OK, '')
            if (.not. self%ready) call status_set(status, FORTSPARSE_NOT_FACTORED, &
                'multifrontal LDL^T: source not factored')
        class default
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'backend copy: source and destination types differ')
        end select
    end subroutine mf_copy_from

    function mf_vector(self, n) result(p)
        class(native_multifrontal_backend_t), intent(inout) :: self
        integer, intent(in) :: n
        real(dp), pointer :: p(:)
        nullify (p)
        if (n < 0) return
        if (associated(self%user_vector)) deallocate (self%user_vector)
        allocate (self%user_vector(n))
        p => self%user_vector
    end function mf_vector

    subroutine mf_free(self)
        class(native_multifrontal_backend_t), intent(inout) :: self
        type(multifrontal_ldlt_t) :: empty
        type(sparse_factor_stats_t) :: none
        self%f = empty
        self%ready = .false.
        self%counted = .false.
        self%refining = .false.
        self%analysed_order = -1
        if (allocated(self%ap)) deallocate (self%ap)
        if (allocated(self%ai)) deallocate (self%ai)
        if (allocated(self%original)) deallocate (self%original)
        if (allocated(self%saved)) deallocate (self%saved, self%residual, self%scratch)
        self%stats = none
        call self%reset_raw_tracking()
    end subroutine mf_free

    subroutine mf_finalize(self)
        type(native_multifrontal_backend_t), intent(inout) :: self
        if (associated(self%user_vector)) deallocate (self%user_vector)
    end subroutine mf_finalize

    ! --------------------------------------------------------------- analysis

    subroutine analyse(self, n, ptr, row, val, ordering, info)
        type(multifrontal_ldlt_t), intent(inout) :: self
        integer, intent(in) :: n, ptr(:), row(:), ordering
        real(dp), intent(in) :: val(:)
        integer, intent(out) :: info
        integer, allocatable :: sp(:), si(:), parent(:), anc(:), post(:), head(:), &
            next(:), stk(:), cc(:), mark(:), p2(:), snode(:), ncol(:), nrow(:), &
            fsfirst(:), fsparent(:), mate(:)
        integer(i8), allocatable :: zeros(:)
        integer :: j, k, p, i, c, top, s, ns, t, nsr, q
        integer(i8) :: totsize, newz, nc0, nc1, nr0, nr1, tot
        real(dp) :: frac
        logical :: merge, symmetric

        info = 0
        call clear(self)
        self%n = n
        ! pattern of A + A^T with ascending rows (the stored pattern may
        ! carry explicit zeros on one side only)
        call symmetric_pattern(n, ptr, row, val, sp, si, self%tpos, symmetric)
        if (.not. symmetric) then
            info = -1
            return
        end if
        allocate (self%perm(n), self%inverse(n), mate(n))
        mate = 0
        if (ordering == FORTSPARSE_ORDER_NATURAL) then
            self%perm = [(j, j=1, n)]
            self%inverse = self%perm
            self%order_used = FORTSPARSE_ORDER_DEFAULT
        else
            call paired_order(n, ptr, row, val, sp, si, self%u, &
                ordering == FORTSPARSE_ORDER_ND, self%perm, self%inverse, mate, &
                self%order_used, info)
            if (info /= 0) return
        end if
        ! elimination tree of the permuted pattern (Liu)
        allocate (parent(n), anc(n))
        parent = 0; anc = 0
        do j = 1, n
            c = self%perm(j)
            do p = sp(c), sp(c + 1) - 1
                i = self%inverse(si(p))
                do while (i /= 0 .and. i < j)
                    k = anc(i)
                    anc(i) = j
                    if (k == 0) parent(i) = j
                    i = k
                end do
            end do
        end do
        ! postorder: a child that is the last one visited stays next to its parent
        allocate (head(n), next(n), post(n), stk(n))
        head = 0
        do j = n, 1, -1
            if (parent(j) /= 0) then
                next(j) = head(parent(j)); head(parent(j)) = j
            end if
        end do
        k = 0
        do j = 1, n
            if (parent(j) /= 0) cycle
            top = 1; stk(1) = j
            do while (top > 0)
                i = stk(top)
                c = head(i)
                if (c == 0) then
                    top = top - 1
                    k = k + 1; post(k) = i
                else
                    head(i) = next(c)
                    top = top + 1; stk(top) = c
                end if
            end do
        end do
        allocate (p2(n))
        do k = 1, n
            p2(k) = self%perm(post(k))
        end do
        self%perm = p2
        do k = 1, n
            self%inverse(self%perm(k)) = k
            anc(post(k)) = k
        end do
        do k = 1, n
            p2(k) = 0
            if (parent(post(k)) /= 0) p2(k) = anc(parent(post(k)))
        end do
        parent = p2
        ! column counts of L from the row subtrees
        allocate (cc(n), mark(n))
        cc = 1; mark = 0
        do i = 1, n
            mark(i) = i
            c = self%perm(i)
            do p = sp(c), sp(c + 1) - 1
                j = self%inverse(si(p))
                do while (j < i .and. mark(j) /= i)
                    cc(j) = cc(j) + 1
                    mark(j) = i
                    j = parent(j)
                end do
            end do
        end do
        self%lnnz = sum(int(cc, i8))
        ! fundamental supernodes; a matched pair is kept in one supernode
        head = 0
        do j = 1, n
            if (parent(j) /= 0) head(parent(j)) = head(parent(j)) + 1
        end do
        allocate (snode(n), fsfirst(n + 1))
        ns = 0
        do j = 1, n
            if (j > 1) then
                if (parent(j - 1) == j .and. ((cc(j - 1) == cc(j) + 1 .and. &
                    head(j) == 1) .or. mate(self%perm(j)) == self%perm(j - 1))) then
                    snode(j) = ns
                    cycle
                end if
            end if
            ns = ns + 1
            fsfirst(ns) = j
            snode(j) = ns
        end do
        fsfirst(ns + 1) = n + 1
        allocate (fsparent(ns), ncol(ns), nrow(ns), zeros(ns))
        do s = 1, ns
            j = fsfirst(s + 1) - 1
            fsparent(s) = 0
            if (parent(j) /= 0) fsparent(s) = snode(parent(j))
            ncol(s) = fsfirst(s + 1) - fsfirst(s)
            nrow(s) = cc(fsfirst(s))
            zeros(s) = 0
        end do
        ! relaxed amalgamation with the CHOLMOD rules: a supernode merges
        ! into its parent when its columns immediately precede the parent's
        do s = ns - 1, 1, -1
            t = fsparent(s)
            if (t /= s + 1) cycle
            nc0 = ncol(s); nr0 = nrow(s); nc1 = ncol(t); nr1 = nrow(t)
            tot = nc0 + nc1
            newz = nc0*(nc0 + nr1 - nr0)
            totsize = tot*(nc0 + nr1) - tot*(tot - 1)/2
            frac = real(zeros(s) + zeros(t) + newz, dp)/real(max(1_i8, totsize), dp)
            merge = newz == 0 .or. tot <= 4 .or. (tot <= 16 .and. frac < 0.8_dp) .or. &
                (tot <= 48 .and. frac < 0.1_dp) .or. frac < 0.05_dp
            if (.not. merge) cycle
            ncol(s) = int(tot)
            nrow(s) = int(nc0 + nr1)
            zeros(s) = zeros(s) + zeros(t) + newz
            fsparent(s) = fsparent(t)
            ncol(t) = 0
        end do
        self%nsn = count(ncol > 0)
        allocate (self%first(self%nsn + 1), self%sparent(self%nsn))
        s = 0
        do t = 1, ns
            if (ncol(t) == 0) cycle
            s = s + 1
            self%first(s) = fsfirst(t)
        end do
        self%first(self%nsn + 1) = n + 1
        do s = 1, self%nsn
            snode(self%first(s):self%first(s + 1) - 1) = s
        end do
        do s = 1, self%nsn
            j = self%first(s + 1) - 1
            self%sparent(s) = 0
            if (parent(j) /= 0) self%sparent(s) = snode(parent(j))
        end do
        ! supernodal symbolic factorization: rows below each supernode
        allocate (self%sptr(self%nsn + 1), self%nchild(self%nsn))
        self%nchild = 0
        deallocate (head, next)
        allocate (head(self%nsn), next(self%nsn))
        head = 0
        do s = self%nsn, 1, -1
            t = self%sparent(s)
            if (t == 0) cycle
            self%nchild(t) = self%nchild(t) + 1
            next(s) = head(t); head(t) = s
        end do
        mark = 0
        allocate (self%srow(max(16, 2*size(si))))
        self%sptr(1) = 1
        self%estimate = 0
        do s = 1, self%nsn
            nsr = 0
            k = self%first(s + 1) - 1
            do j = self%first(s), k
                c = self%perm(j)
                do p = sp(c), sp(c + 1) - 1
                    i = self%inverse(si(p))
                    if (i > k .and. mark(i) /= s) then
                        mark(i) = s
                        call push(i)
                    end if
                end do
            end do
            q = head(s)
            do while (q /= 0)
                do p = self%sptr(q), self%sptr(q + 1) - 1
                    i = self%srow(p)
                    if (i > k .and. mark(i) /= s) then
                        mark(i) = s
                        call push(i)
                    end if
                end do
                q = next(q)
            end do
            ! ascending rows keep child rows in the parent's relative order
            call sort_int(self%srow(self%sptr(s):self%sptr(s) + nsr - 1))
            self%sptr(s + 1) = self%sptr(s) + nsr
            j = self%first(s + 1) - self%first(s)
            self%estimate = self%estimate + int(j, i8)*(j + nsr)
        end do
        ! lower original entries per new column
        allocate (self%aptr(n + 1))
        self%aptr = 0
        do j = 1, n
            c = self%perm(j)
            do p = ptr(c), ptr(c + 1) - 1
                if (self%inverse(row(p)) >= j) self%aptr(j + 1) = self%aptr(j + 1) + 1
            end do
        end do
        self%aptr(1) = 1
        do j = 1, n
            self%aptr(j + 1) = self%aptr(j + 1) + self%aptr(j)
        end do
        allocate (self%arow(self%aptr(n + 1) - 1), self%asrc(self%aptr(n + 1) - 1))
        do j = 1, n
            c = self%perm(j)
            q = self%aptr(j)
            do p = ptr(c), ptr(c + 1) - 1
                i = self%inverse(row(p))
                if (i >= j) then
                    self%arow(q) = i; self%asrc(q) = p; q = q + 1
                end if
            end do
        end do
        allocate (self%map(n))
    contains
        subroutine push(v)
            integer, intent(in) :: v
            integer, allocatable :: tmp(:)
            if (self%sptr(s) + nsr > size(self%srow)) then
                allocate (tmp(2*size(self%srow)))
                tmp(1:size(self%srow)) = self%srow
                call move_alloc(tmp, self%srow)
            end if
            self%srow(self%sptr(s) + nsr) = v
            nsr = nsr + 1
        end subroutine push
    end subroutine analyse

    subroutine clear(self)
        type(multifrontal_ldlt_t), intent(inout) :: self
        real(dp) :: u
        u = self%u
        self = multifrontal_ldlt_t()
        self%u = u
    end subroutine clear

    subroutine symmetric_pattern(n, ptr, row, val, sp, si, tpos, symmetric)
        ! Ascending rows of the pattern of A + A^T without the diagonal, from
        ! a merge of each column of A with the same column of A^T. On the way
        ! every entry is checked against its transpose: values must agree,
        ! and an entry without a stored counterpart must be zero.
        integer, intent(in) :: n, ptr(:), row(:)
        real(dp), intent(in) :: val(:)
        integer, allocatable, intent(out) :: sp(:), si(:), tpos(:)
        logical, intent(out) :: symmetric
        integer, allocatable :: tp(:), ti(:), ts(:), cnt(:)
        real(dp), allocatable :: tv(:)
        integer :: j, p, i, q, a, b, ea, eb, ra, rb
        allocate (tp(n + 1), cnt(n + 1))
        cnt = 0
        do p = 1, ptr(n + 1) - 1
            cnt(row(p) + 1) = cnt(row(p) + 1) + 1
        end do
        tp(1) = 1
        do j = 1, n
            tp(j + 1) = tp(j) + cnt(j + 1)
        end do
        allocate (ti(max(1, ptr(n + 1) - 1)), tv(max(1, ptr(n + 1) - 1)), &
            ts(max(1, ptr(n + 1) - 1)), tpos(ptr(n + 1) - 1))
        tpos = 0
        cnt(1:n) = tp(1:n)
        do j = 1, n
            do p = ptr(j), ptr(j + 1) - 1
                i = row(p)
                ti(cnt(i)) = j
                tv(cnt(i)) = val(p)
                ts(cnt(i)) = p
                cnt(i) = cnt(i) + 1
            end do
        end do
        allocate (sp(n + 1), si(max(1, 2*(ptr(n + 1) - 1))))
        symmetric = .false.
        sp(1) = 1
        q = 1
        do j = 1, n
            a = ptr(j); ea = ptr(j + 1)
            b = tp(j); eb = tp(j + 1)
            do while (a < ea .or. b < eb)
                ra = huge(0); rb = huge(0)
                if (a < ea) ra = row(a)
                if (b < eb) rb = ti(b)
                if (ra == rb) then
                    if (val(a) /= tv(b)) return
                    tpos(a) = ts(b)
                    i = ra; a = a + 1; b = b + 1
                else if (ra < rb) then
                    if (val(a) /= 0.0_dp) return
                    i = ra; a = a + 1
                else
                    if (tv(b) /= 0.0_dp) return
                    i = rb; b = b + 1
                end if
                if (i /= j) then
                    si(q) = i; q = q + 1
                end if
            end do
            sp(j + 1) = q
        end do
        symmetric = .true.
    end subroutine symmetric_pattern

    subroutine paired_order(n, ptr, row, val, sp, si, u, nested, perm, inverse, mate, &
            used, info)
        ! Fill-reducing order of the graph in which a variable whose diagonal
        ! fails the pivot threshold u is merged with its largest unmatched
        ! neighbour; the pair is ordered consecutively, the member with the
        ! larger diagonal first.
        integer, intent(in) :: n, ptr(:), row(:), sp(:), si(:)
        real(dp), intent(in) :: val(:), u
        logical, intent(in) :: nested
        integer, intent(out) :: perm(:), inverse(:), mate(:), used
        integer, intent(out) :: info
        integer, allocatable :: sv(:), lead(:), gp(:), gi(:), mark(:), &
            sperm(:), sinv(:)
        real(dp), allocatable :: dg(:), cmax(:)
        real(dp) :: best
        integer :: c, p, i, v, nsv, s, k, q, m

        allocate (dg(n), cmax(n))
        mate = 0; dg = 0; cmax = 0
        do c = 1, n
            do p = ptr(c), ptr(c + 1) - 1
                if (row(p) == c) then
                    dg(c) = abs(val(p))
                else
                    cmax(c) = max(cmax(c), abs(val(p)))
                end if
            end do
        end do
        do c = 1, n
            if (mate(c) /= 0 .or. dg(c) >= u*cmax(c)) cycle
            best = 0; v = 0
            do p = ptr(c), ptr(c + 1) - 1
                i = row(p)
                if (i == c .or. mate(i) /= 0) cycle
                if (abs(val(p)) > best) then
                    best = abs(val(p)); v = i
                end if
            end do
            if (v /= 0) then
                mate(c) = v; mate(v) = c
            end if
        end do
        allocate (sv(n), lead(n))
        nsv = 0
        do c = 1, n
            if (mate(c) /= 0 .and. mate(c) < c) then
                sv(c) = sv(mate(c))
            else
                nsv = nsv + 1
                sv(c) = nsv
                lead(nsv) = c
            end if
        end do
        allocate (gp(nsv + 1), gi(max(1, sp(n + 1) - 1)), mark(nsv))
        mark = 0
        gp(1) = 1
        q = 1
        do s = 1, nsv
            mark(s) = s
            c = lead(s)
            do k = 1, 2
                if (k == 2) c = mate(lead(s))
                if (c == 0) exit
                do p = sp(c), sp(c + 1) - 1
                    i = sv(si(p))
                    if (mark(i) /= s) then
                        mark(i) = s
                        gi(q) = i; q = q + 1
                    end if
                end do
            end do
            gp(s + 1) = q
        end do
        call transpose_pattern(nsv, gp, gi)
        allocate (sperm(nsv), sinv(nsv))
        info = 1
        used = FORTSPARSE_ORDER_ND
        if (nested) info = order_nd(int(nsv, c_int), int(q - 1, c_int), gp, gi, sperm, sinv)
        if (info /= 0) then
            used = FORTSPARSE_ORDER_AMD
            info = order_amd(int(nsv, c_int), int(q - 1, c_int), gp, gi, sperm, sinv)
        end if
        if (info /= 0) return
        k = 0
        do s = 1, nsv
            c = lead(sperm(s))
            v = mate(c)
            if (v /= 0) then
                if (dg(v) > dg(c)) then
                    m = c; c = v; v = m
                end if
            end if
            k = k + 1; perm(k) = c
            if (v /= 0) then
                k = k + 1; perm(k) = v
            end if
        end do
        do k = 1, n
            inverse(perm(k)) = k
        end do
    end subroutine paired_order

    ! ---------------------------------------------------------------- numeric

    subroutine factor(self, val)
        type(multifrontal_ldlt_t), intent(inout) :: self
        real(dp), intent(in) :: val(:)
        integer :: s, ns, nd, m, nfs, npiv, k, j, p, q, mc, nstk, ncb, base, fidx_top, itop
        integer(i8) :: sp
        integer, allocatable :: cbm(:), cbn(:), csp_i(:)
        integer(i8), allocatable :: csp(:)

        ns = self%nsn
        self%inertia = inertia_count_t()
        self%delayed = 0; self%two_by_two = 0; self%max_front = 0; self%flops = 0
        if (allocated(self%fm)) deallocate (self%fm, self%fnpiv, self%fidx_ptr, self%lptr)
        allocate (self%fm(ns), self%fnpiv(ns), self%fidx_ptr(ns + 1), self%lptr(ns + 1))
        if (.not. allocated(self%fidx)) allocate (self%fidx(2*self%n + 16))
        if (.not. allocated(self%l)) allocate (self%l(max(1024_i8, 3*self%estimate/2)))
        if (allocated(self%pk)) then
            if (size(self%pk) /= self%n) deallocate (self%pk)
        end if
        if (.not. allocated(self%pk)) allocate (self%pk(self%n))
        if (.not. allocated(self%d)) allocate (self%d(self%n), self%e(self%n), &
            self%rowmax(self%n), self%hist(self%n))
        if (.not. allocated(self%stack)) allocate (self%stack(1024))
        if (.not. allocated(self%istack)) allocate (self%istack(1024))
        if (.not. allocated(self%work)) allocate (self%work(1024))
        if (.not. allocated(self%pos)) allocate (self%pos(1024))
        self%rowmax = 0
        self%hist = 0
        do j = 1, self%n
            do p = self%aptr(j), self%aptr(j + 1) - 1
                k = self%arow(p)
                self%rowmax(j) = max(self%rowmax(j), abs(val(self%asrc(p))))
                self%rowmax(k) = max(self%rowmax(k), abs(val(self%asrc(p))))
            end do
        end do
        allocate (cbm(ns), cbn(ns), csp(ns), csp_i(ns))
        nstk = 0; sp = 0; itop = 0
        self%fidx_ptr(1) = 1
        self%lptr(1) = 1
        fidx_top = 1
        do s = 1, ns
            ! fully summed: variables delayed by the children, then own columns
            nd = sum(cbn(nstk - self%nchild(s) + 1:nstk))
            ncb = self%sptr(s + 1) - self%sptr(s)
            nfs = nd + self%first(s + 1) - self%first(s)
            m = nfs + ncb
            call grow_int(self%fidx, fidx_top + m)
            base = fidx_top - 1
            q = 0
            do k = nstk - self%nchild(s) + 1, nstk
                self%fidx(base + q + 1:base + q + cbn(k)) = &
                    self%istack(csp_i(k) + 1:csp_i(k) + cbn(k))
                q = q + cbn(k)
            end do
            do j = self%first(s), self%first(s + 1) - 1
                q = q + 1
                self%fidx(base + q) = j
            end do
            self%fidx(base + q + 1:base + m) = self%srow(self%sptr(s):self%sptr(s + 1) - 1)
            do q = 1, m
                self%map(self%fidx(base + q)) = q
            end do
            ! the front is factored in place at the end of the factor storage
            call grow_real(self%l, self%lptr(s) + int(m, i8)*m)
            call reserve_int(self%pos, m)
            call assemble(self, self%l(self%lptr(s):), m, s, val, nstk, cbm, csp, csp_i)
            if (self%nchild(s) > 0) then
                k = nstk - self%nchild(s) + 1
                sp = csp(k)
                itop = csp_i(k)
                nstk = k - 1
            end if
            self%max_front = max(self%max_front, m)
            call factor_front(self, self%l(self%lptr(s):), m, nfs, self%sparent(s) == 0, &
                self%fidx(base + 1:base + m), npiv)
            self%fm(s) = m
            self%fnpiv(s) = npiv
            self%delayed = self%delayed + (nfs - npiv)
            self%flops = self%flops + int(npiv, i8)*m*m - int(npiv, i8)*npiv*m + &
                int(npiv, i8)*npiv*npiv/3
            fidx_top = fidx_top + m
            self%fidx_ptr(s + 1) = fidx_top
            ! contribution block: positions npiv+1:m (delayed variables first)
            if (self%sparent(s) /= 0) then
                mc = m - npiv
                nstk = nstk + 1
                cbm(nstk) = mc
                cbn(nstk) = nfs - npiv
                csp(nstk) = sp
                csp_i(nstk) = itop
                call grow_real(self%stack, sp + int(mc, i8)*mc)
                call grow_int(self%istack, itop + mc)
                do j = 1, mc
                    self%stack(sp + (j - 1)*int(mc, i8) + j:sp + j*int(mc, i8)) = &
                        self%l(self%lptr(s) + npiv + (npiv + j - 1)*int(m, i8) + (j - 1): &
                               self%lptr(s) + m - 1 + (npiv + j - 1)*int(m, i8))
                end do
                self%istack(itop + 1:itop + mc) = self%fidx(base + npiv + 1:base + m)
                sp = sp + int(mc, i8)*mc
                itop = itop + mc
            end if
            self%lptr(s + 1) = self%lptr(s) + int(m, i8)*npiv
        end do
        if (allocated(self%y)) deallocate (self%y)
        allocate (self%y(max(1, self%max_front)))
    end subroutine factor

    subroutine assemble(self, f, m, s, val, nstk, cbm, csp, csp_i)
        ! Zero the lower front, add the original entries of the supernode's
        ! columns and extend-add the children's contribution blocks. Child
        ! rows keep their relative order in the parent (delayed variables
        ! first, then ascending rows), so every entry stays in the lower part.
        type(multifrontal_ldlt_t), intent(inout) :: self
        integer, intent(in) :: m, s, nstk, cbm(:), csp_i(:)
        real(dp), intent(inout) :: f(m, *)
        real(dp), intent(in) :: val(:)
        integer(i8), intent(in) :: csp(:)
        integer :: j, p, k, a, b, pb, mc
        integer(i8) :: c0
        do j = 1, m
            f(j:m, j) = 0.0_dp
        end do
        do j = self%first(s), self%first(s + 1) - 1
            pb = self%map(j)
            do p = self%aptr(j), self%aptr(j + 1) - 1
                f(self%map(self%arow(p)), pb) = f(self%map(self%arow(p)), pb) + &
                    val(self%asrc(p))
            end do
        end do
        do k = nstk - self%nchild(s) + 1, nstk
            mc = cbm(k)
            do a = 1, mc
                self%pos(a) = self%map(self%istack(csp_i(k) + a))
            end do
            c0 = csp(k)
            do b = 1, mc
                pb = self%pos(b)
                do a = b, mc
                    f(self%pos(a), pb) = f(self%pos(a), pb) + &
                        self%stack(c0 + a + (b - 1)*int(mc, i8))
                end do
            end do
        end do
    end subroutine assemble

    subroutine factor_front(self, f, m, nfs, root, idx, npiv)
        ! Eliminate up to nfs fully summed variables of the m x m front f
        ! (lower triangle). Pivots end at positions 1:npiv, rejected
        ! variables at npiv+1:nfs, updated; the contribution block nfs+1:m
        ! receives the Schur update. Pivots are sought in windows of columns
        ! updated by rank-1/2 steps; the rest of the front by BLAS-3.
        type(multifrontal_ldlt_t), intent(inout) :: self
        integer, intent(in) :: m, nfs
        real(dp), intent(inout) :: f(m, *)
        logical, intent(in) :: root
        integer, intent(inout) :: idx(:)
        integer, intent(out) :: npiv
        integer :: p, q, pend, pend0, lim, k, r, kind, j, i
        real(dp) :: akk, arr, akr, amaxk, amaxr, det, u, d11, d22, d21, t1, t2, &
            gk, gr, sc
        logical :: wide, force

        u = self%u
        p = 1
        wide = .false.
        force = .false.
        do while (p <= nfs)
            pend0 = nfs
            if (.not. wide) pend0 = min(nfs, p + window - 1)
            pend = pend0
            ! untested candidates q:lim, candidates that failed in this
            ! window lim+1:pend (still updated, retested in the next window)
            lim = pend
            q = p
            do while (q <= lim)
                k = q
                if (force) then
                    ! a root front that delays nothing: Bunch-Kaufman on column q
                    force = .false.
                    call column_max(f, m, q, k, 0, amaxk, r, pend)
                    akk = f(k, k)
                    kind = 1
                    if (abs(akk) < bk_alpha*amaxk) then
                        call column_max(f, m, q, r, 0, amaxr, i, pend)
                        if (abs(akk)*amaxr < bk_alpha*amaxk*amaxk) then
                            if (abs(f(r, r)) >= bk_alpha*amaxr) then
                                k = r
                                amaxk = amaxr
                            else
                                kind = 2
                            end if
                        end if
                    end if
                else
                    call column_max(f, m, q, k, 0, amaxk, r, pend)
                    akk = f(k, k)
                    kind = 0
                    if (abs(akk) >= u*amaxk) then
                        kind = 1
                    else if (r /= 0) then
                        ! |D^-1| (amax_k, amax_r)^T <= 1/u with D scaled by its
                        ! largest entry sc, so no product over- or underflows
                        akr = f(r, k)
                        arr = f(r, r)
                        sc = max(abs(akk), abs(akr), abs(arr))
                        det = (akk/sc)*(arr/sc) - (akr/sc)**2
                        if (det /= 0.0_dp) then
                            call column_max(f, m, q, k, r, amaxk, i, pend)
                            call column_max(f, m, q, r, k, amaxr, i, pend)
                            if ((abs(arr/sc)*amaxk + abs(akr/sc)*amaxr)*u <= sc*abs(det) &
                                .and. (abs(akr/sc)*amaxk + abs(akk/sc)*amaxr)*u <= &
                                sc*abs(det)) kind = 2
                        end if
                    end if
                    if (kind == 0) then
                        call swap(f, m, idx, q, lim)
                        lim = lim - 1
                        cycle
                    end if
                end if
                if (kind == 1) then
                    call swap(f, m, idx, q, k)
                    d11 = f(q, q)
                    gk = max(self%rowmax(idx(q)) + self%hist(idx(q)) + &
                        own_history(self, f, m, idx, q, q), amaxk, abs(d11))
                    call inertia_add_pivot(self%inertia, 1, d11, 0.0_dp, 0.0_dp, &
                        inertia_tolerance(self%n, gk))
                    if (d11 /= 0.0_dp) then
                        do j = q + 1, pend
                            t1 = f(j, q)/d11
                            if (t1 /= 0.0_dp) f(j:m, j) = f(j:m, j) - t1*f(j:m, q)
                        end do
                        f(q + 1:m, q) = f(q + 1:m, q)/d11
                    end if
                    self%pk(idx(q)) = 1
                    self%d(idx(q)) = d11
                    self%e(idx(q)) = 0
                    q = q + 1
                else
                    ! k = q, r > q in the window; a failed partner joins q+1
                    ! and the untested column it displaces stays untested
                    if (r > lim) then
                        call swap(f, m, idx, lim + 1, r)
                        call swap(f, m, idx, q + 1, lim + 1)
                        lim = lim + 1
                    else
                        call swap(f, m, idx, q + 1, r)
                    end if
                    d11 = f(q, q); d21 = f(q + 1, q); d22 = f(q + 1, q + 1)
                    sc = max(abs(d11), abs(d21), abs(d22))
                    det = (d11/sc)*(d22/sc) - (d21/sc)**2
                    call column_max(f, m, q, q, q + 1, amaxk, i, q)
                    call column_max(f, m, q, q + 1, q, amaxr, i, q)
                    gk = max(self%rowmax(idx(q)) + self%hist(idx(q)) + &
                        own_history(self, f, m, idx, q, q), amaxk, abs(d11), abs(d21))
                    gr = max(self%rowmax(idx(q + 1)) + self%hist(idx(q + 1)) + &
                        own_history(self, f, m, idx, q, q + 1), amaxr, abs(d22), abs(d21))
                    call inertia_add_scaled_block(self%inertia, self%n, d11, abs(d21), &
                        d22, gk, gr)
                    self%two_by_two = self%two_by_two + 1
                    ! multipliers [t1 t2] = [a_jq a_jq+1] D^-1
                    do j = q + 2, pend
                        t1 = ((f(j, q)*(d22/sc) - f(j, q + 1)*(d21/sc))/det)/sc
                        t2 = ((f(j, q + 1)*(d11/sc) - f(j, q)*(d21/sc))/det)/sc
                        f(j:m, j) = f(j:m, j) - t1*f(j:m, q) - t2*f(j:m, q + 1)
                    end do
                    do i = q + 2, m
                        t1 = ((f(i, q)*(d22/sc) - f(i, q + 1)*(d21/sc))/det)/sc
                        t2 = ((f(i, q + 1)*(d11/sc) - f(i, q)*(d21/sc))/det)/sc
                        f(i, q) = t1; f(i, q + 1) = t2
                    end do
                    f(q + 1, q) = 0.0_dp
                    self%pk(idx(q)) = 2; self%d(idx(q)) = d11; self%e(idx(q)) = d21
                    self%pk(idx(q + 1)) = 3; self%d(idx(q + 1)) = d22
                    self%e(idx(q + 1)) = d21
                    q = q + 2
                end if
            end do
            if (q == p) then
                ! no pivot in this window: widen it, force a pivot at a root,
                ! or delay the rest to the parent
                if (pend0 < nfs) then
                    wide = .true.
                else if (root) then
                    wide = .true.
                    force = .true.
                else
                    exit
                end if
                cycle
            end if
            wide = .false.
            if (pend < nfs) call update_block(self, f, m, idx, p, q - 1, pend + 1, nfs)
            p = q
        end do
        npiv = p - 1
        if (npiv > 0 .and. m > nfs) call update_block(self, f, m, idx, 1, npiv, nfs + 1, m)
        do i = npiv + 1, m
            self%hist(idx(i)) = self%hist(idx(i)) + own_history(self, f, m, idx, npiv + 1, i)
        end do
    end subroutine factor_front

    real(dp) function own_history(self, f, m, idx, last, row) result(h)
        ! sum over the pivots 1:last-1 of this front of l^2 |d| for one row
        type(multifrontal_ldlt_t), intent(in) :: self
        integer, intent(in) :: m, last, row
        real(dp), intent(in) :: f(m, *)
        integer, intent(in) :: idx(:)
        integer :: k
        real(dp) :: l1, l2, e
        h = 0
        k = 1
        do while (k < min(last, row))
            if (self%pk(idx(k)) == 2) then
                l1 = abs(f(row, k)); l2 = abs(f(row, k + 1)); e = abs(self%e(idx(k)))
                h = h + abs(self%d(idx(k)))*l1*l1 + abs(self%d(idx(k + 1)))*l2*l2 + &
                    2*e*l1*l2 + e*(l1*l1 + l2*l2)
                k = k + 2
            else
                h = h + abs(self%d(idx(k)))*f(row, k)**2
                k = k + 1
            end if
        end do
    end function own_history

    subroutine update_block(self, f, m, idx, p0, p1, c0, c1)
        ! Lower part of columns c0:c1 (rows j:m) -= L D L^T over the pivots
        ! p0:p1, with W = L D formed once.
        type(multifrontal_ldlt_t), intent(inout) :: self
        integer, intent(in) :: m, p0, p1, c0, c1
        real(dp), intent(inout) :: f(m, *)
        integer, intent(in) :: idx(:)
        integer, parameter :: cb = 64
        integer :: j, kp, nr, j0, jb
        integer(i8) :: o
        real(dp) :: d11, d21, d22

        kp = p1 - p0 + 1
        nr = m - c0 + 1
        call reserve_real(self%work, int(nr, i8)*kp)
        j = p0
        do while (j <= p1)
            o = int(j - p0, i8)*nr
            if (self%pk(idx(j)) == 2) then
                d11 = self%d(idx(j)); d21 = self%e(idx(j)); d22 = self%d(idx(j + 1))
                self%work(o + 1:o + nr) = f(c0:m, j)*d11 + f(c0:m, j + 1)*d21
                self%work(o + nr + 1:o + 2*nr) = f(c0:m, j)*d21 + f(c0:m, j + 1)*d22
                j = j + 2
            else
                self%work(o + 1:o + nr) = f(c0:m, j)*self%d(idx(j))
                j = j + 1
            end if
        end do
        do j0 = c0, c1, cb
            jb = min(cb, c1 - j0 + 1)
            call dgemm('N', 'T', m - j0 + 1, jb, kp, -1.0_dp, self%work(1 + j0 - c0), nr, &
                f(j0, p0), m, 1.0_dp, f(j0, j0), m)
        end do
    end subroutine update_block

    subroutine column_max(f, m, p, k, skip, amax, r, nfs)
        ! Largest |a_ik| over the remaining rows i >= p, i /= k, i /= skip;
        ! r is the row of the largest fully summed (i <= nfs) one. The skipped
        ! entry is zeroed meanwhile, so the long loop over the contribution
        ! rows has no branch.
        integer, intent(in) :: m, p, k, skip, nfs
        real(dp), intent(inout) :: f(m, *)
        real(dp), intent(out) :: amax
        integer, intent(out) :: r
        integer :: i, si, sj
        real(dp) :: v, saved
        if (skip /= 0) then
            si = max(k, skip); sj = min(k, skip)
            saved = f(si, sj)
            f(si, sj) = 0.0_dp
        end if
        amax = 0; r = 0
        do i = p, k - 1
            v = abs(f(k, i))
            if (v > amax) then
                amax = v; r = i
            end if
        end do
        do i = k + 1, nfs
            v = abs(f(i, k))
            if (v > amax) then
                amax = v; r = i
            end if
        end do
        v = 0
        do i = max(k, nfs) + 1, m
            v = max(v, abs(f(i, k)))
        end do
        amax = max(amax, v)
        if (skip /= 0) f(si, sj) = saved
    end subroutine column_max

    subroutine swap(f, m, idx, i1, i2)
        ! symmetric interchange of rows/columns i1 and i2 in lower storage,
        ! including the stored L columns left of them
        integer, intent(in) :: m, i1, i2
        real(dp), intent(inout) :: f(m, *)
        integer, intent(inout) :: idx(:)
        integer :: a, b, j, t
        real(dp) :: tmp
        if (i1 == i2) return
        a = min(i1, i2); b = max(i1, i2)
        do j = 1, a - 1
            tmp = f(a, j); f(a, j) = f(b, j); f(b, j) = tmp
        end do
        tmp = f(a, a); f(a, a) = f(b, b); f(b, b) = tmp
        do j = a + 1, b - 1
            tmp = f(j, a); f(j, a) = f(b, j); f(b, j) = tmp
        end do
        do j = b + 1, m
            tmp = f(j, a); f(j, a) = f(j, b); f(j, b) = tmp
        end do
        t = idx(a); idx(a) = idx(b); idx(b) = t
    end subroutine swap

    subroutine solve(self, b)
        ! b: right-hand side in the original numbering, overwritten by x
        type(multifrontal_ldlt_t), intent(inout) :: self
        real(dp), intent(inout) :: b(:)
        real(dp), allocatable :: x(:)
        integer :: s, m, np, base, j, g
        integer(i8) :: lp
        real(dp) :: d11, d21, d22, det, t1, t2, sc
        allocate (x(self%n))
        x = b(self%perm)
        associate (y => self%y)
            do s = 1, self%nsn
                m = self%fm(s); np = self%fnpiv(s); base = self%fidx_ptr(s) - 1
                if (np == 0) cycle
                lp = self%lptr(s)
                y(1:m) = x(self%fidx(base + 1:base + m))
                call dtrsv('L', 'N', 'U', np, self%l(lp), m, y, 1)
                if (m > np) then
                    call dgemv('N', m - np, np, -1.0_dp, self%l(lp + np), m, y, 1, &
                        1.0_dp, y(np + 1), 1)
                    x(self%fidx(base + np + 1:base + m)) = y(np + 1:m)
                end if
                j = 1
                do while (j <= np)
                    g = self%fidx(base + j)
                    if (self%pk(g) == 2) then
                        d11 = self%d(g); d21 = self%e(g); d22 = self%d(self%fidx(base + j + 1))
                        sc = max(abs(d11), abs(d21), abs(d22))
                        det = (d11/sc)*(d22/sc) - (d21/sc)**2
                        t1 = (((d22/sc)*y(j) - (d21/sc)*y(j + 1))/det)/sc
                        t2 = (((d11/sc)*y(j + 1) - (d21/sc)*y(j))/det)/sc
                        y(j) = t1; y(j + 1) = t2
                        j = j + 2
                    else
                        y(j) = y(j)/self%d(g)
                        j = j + 1
                    end if
                end do
                x(self%fidx(base + 1:base + np)) = y(1:np)
            end do
            do s = self%nsn, 1, -1
                m = self%fm(s); np = self%fnpiv(s); base = self%fidx_ptr(s) - 1
                if (np == 0) cycle
                lp = self%lptr(s)
                y(1:m) = x(self%fidx(base + 1:base + m))
                if (m > np) call dgemv('T', m - np, np, -1.0_dp, self%l(lp + np), m, &
                    y(np + 1), 1, 1.0_dp, y, 1)
                call dtrsv('L', 'T', 'U', np, self%l(lp), m, y, 1)
                x(self%fidx(base + 1:base + np)) = y(1:np)
            end do
        end associate
        b(self%perm) = x
    end subroutine solve

    ! ---------------------------------------------------------------- helpers

    subroutine sort_int(a)
        integer, intent(inout) :: a(:)
        integer :: i, j, v, gap
        gap = size(a)/2
        do while (gap > 0)
            do i = gap + 1, size(a)
                v = a(i)
                j = i
                do while (j > gap)
                    if (a(j - gap) <= v) exit
                    a(j) = a(j - gap)
                    j = j - gap
                end do
                a(j) = v
            end do
            gap = gap/2
        end do
    end subroutine sort_int

    subroutine transpose_pattern(n, ptr, idx)
        ! a symmetric pattern is its own transpose: rewrite with ascending rows
        integer, intent(in) :: n
        integer, intent(inout) :: ptr(:), idx(:)
        integer, allocatable :: tp(:), ti(:)
        integer :: j, p, i
        allocate (tp(n + 1), ti(max(1, ptr(n + 1) - 1)))
        tp = 0
        do p = 1, ptr(n + 1) - 1
            tp(idx(p) + 1) = tp(idx(p) + 1) + 1
        end do
        tp(1) = 1
        do j = 1, n
            tp(j + 1) = tp(j + 1) + tp(j)
        end do
        do j = 1, n
            do p = ptr(j), ptr(j + 1) - 1
                i = idx(p)
                ti(tp(i)) = j
                tp(i) = tp(i) + 1
            end do
        end do
        idx(1:ptr(n + 1) - 1) = ti(1:ptr(n + 1) - 1)
    end subroutine transpose_pattern

    subroutine reserve_real(a, need)
        ! workspace growth without copying the old contents
        real(dp), allocatable, intent(inout) :: a(:)
        integer(i8), intent(in) :: need
        integer(i8) :: old
        old = size(a, kind=i8)
        if (old >= need) return
        deallocate (a)
        allocate (a(max(need, 2*old)))
    end subroutine reserve_real

    subroutine reserve_int(a, need)
        integer, allocatable, intent(inout) :: a(:)
        integer, intent(in) :: need
        integer :: old
        old = size(a)
        if (old >= need) return
        deallocate (a)
        allocate (a(max(need, 2*old)))
    end subroutine reserve_int

    subroutine grow_real(a, need)
        real(dp), allocatable, intent(inout) :: a(:)
        integer(i8), intent(in) :: need
        real(dp), allocatable :: t(:)
        if (size(a, kind=i8) >= need) return
        allocate (t(max(need, 2*size(a, kind=i8))))
        t(1:size(a, kind=i8)) = a
        call move_alloc(t, a)
    end subroutine grow_real

    subroutine grow_int(a, need)
        integer, allocatable, intent(inout) :: a(:)
        integer, intent(in) :: need
        integer, allocatable :: t(:)
        if (size(a) >= need) return
        allocate (t(max(need, 2*size(a))))
        t(1:size(a)) = a
        call move_alloc(t, a)
    end subroutine grow_int

end module fortsparse_native_multifrontal
