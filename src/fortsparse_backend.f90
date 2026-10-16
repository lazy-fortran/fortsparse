module fortsparse_backend
    ! Abstract sparse-solver backend. A concrete backend retains a
    ! factorization across solves and maps its native error codes onto the
    ! fortsparse status codes. The solver holds one of these polymorphically and
    ! dispatches factor/solve/free through the deferred bindings, so adding a
    ! backend never touches client code or the solver dispatch logic.
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csc, only: csc_t, csc_z_t
    use fortsparse_memory, only: memory_budget_t, memory_checked_add, memory_checked_mul
    use fortsparse_status, only: fortsparse_status_t, status_set, &
        FORTSPARSE_OK, FORTSPARSE_BACKEND_UNAVAILABLE, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OUT_OF_MEMORY, &
        FORTSPARSE_INDEX_OVERFLOW, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_RESOURCE_LIMIT
    implicit none
    private

    public :: backend_factor_real_raw, backend_refactor_real_raw
    public :: backend_factor_real_raw_budgeted, backend_refactor_real_raw_budgeted
    public :: backend_factor_complex_raw_budgeted, backend_refactor_complex_raw_budgeted
    public :: sparse_backend_t
    public :: sparse_options_t
    public :: sparse_factor_stats_t

    ! Matrix structure hint for a factorization (sparse_options_t%symmetry).
    ! SPD and INDEFINITE factor only the lower triangle of the supplied full
    ! symmetric CSC matrix (MUMPS); GENERAL factors the matrix as given.
    ! SuperLU factors the full matrix with every hint; a symmetric hint
    ! selects its symmetric-mode ordering and threshold pivoting.
    integer, parameter, public :: FORTSPARSE_SYM_GENERAL = 0
    integer, parameter, public :: FORTSPARSE_SYM_SPD = 1
    integer, parameter, public :: FORTSPARSE_SYM_INDEFINITE = 2
    ! Complex only: A = A^H, MUMPS reads the lower triangle (real treats it as
    ! INDEFINITE). For complex matrices SPD means Hermitian positive definite,
    ! and INDEFINITE means complex symmetric, A = A^T (not A^H).
    integer, parameter, public :: FORTSPARSE_SYM_HERMITIAN = 3

    ! Fill-reducing ordering (sparse_options_t%ordering).
    integer, parameter, public :: FORTSPARSE_ORDER_DEFAULT = 0
    integer, parameter, public :: FORTSPARSE_ORDER_AMD = 1
    integer, parameter, public :: FORTSPARSE_ORDER_ND = 2

    integer, parameter, public :: FORTSPARSE_OP_NORMAL = 0
    integer, parameter, public :: FORTSPARSE_OP_TRANSPOSE = 1
    integer, parameter, public :: FORTSPARSE_OP_ADJOINT = 2

    integer, parameter, public :: FORTSPARSE_MEMORY_AUTO = 0
    integer, parameter, public :: FORTSPARSE_MEMORY_LOW = 1

    ! Native CPU real GENERAL LU method (sparse_options_t%lu_method). AUTO
    ! behaves as SCALAR until its thresholds are calibrated; SUPERNODAL is
    ! the static-structure supernodal LU (docs/lu-supernodal.md).
    integer, parameter, public :: FORTSPARSE_LU_AUTO = 0
    integer, parameter, public :: FORTSPARSE_LU_SCALAR = 1
    integer, parameter, public :: FORTSPARSE_LU_SUPERNODAL = 2

    ! sparse_factor_stats_t%matching_used.
    integer, parameter, public :: FORTSPARSE_MATCHING_NONE = 0
    integer, parameter, public :: FORTSPARSE_MATCHING_STRUCTURAL = 1
    integer, parameter, public :: FORTSPARSE_MATCHING_WEIGHTED = 2

    ! Backend tuning, copied from the solver to the backend at every factor.
    ! Backends that have no use for a field ignore it.
    type :: sparse_options_t
        integer :: symmetry = FORTSPARSE_SYM_GENERAL
        integer :: ordering = FORTSPARSE_ORDER_DEFAULT
        logical :: structural_matching = .false. ! native CPU GENERAL, nonzero row matching
        integer :: nthreads = 0 ! 0: backend default (OMP_NUM_THREADS)
        integer :: device = -1 ! native: -1 CPU, >=0 explicit CUDA device
        logical :: verbose = .false.
        ! Strict native array-payload limits; allocator/runtime/RSS excluded.
        integer(i8) :: host_memory_limit_bytes = 0_i8
        integer(i8) :: device_memory_limit_bytes = 0_i8
        integer :: memory_policy = FORTSPARSE_MEMORY_AUTO
        ! UMFPACK IPC: negative keeps its native pivot defaults; [0,1] overrides.
        real(dp) :: pivot_tolerance = -1.0_dp
        ! UMFPACK estimated peak only, excluding caller/IPC storage; 0: no gate.
        integer(i8) :: memory_budget_bytes = 0_i8
        ! Native CPU real GENERAL LU method; see FORTSPARSE_LU_*.
        integer :: lu_method = FORTSPARSE_LU_SCALAR
        ! Supernodal LU: static pivots |p| < pivot_perturbation*||A~||_inf
        ! are replaced by sign(p)*pivot_perturbation*||A~||_inf; 0 < value < 1.
        real(dp) :: pivot_perturbation = 1.0e-13_dp
        ! Supernodal LU: validate every factor with one refined solve against
        ! the original matrix (always done when pivots were perturbed).
        logical :: verify_factor = .true.
    end type sparse_options_t

    ! Factorization statistics; -1 / -1.0 mark a quantity the backend does not
    ! report.
    type :: sparse_factor_stats_t
        integer(i8) :: factor_nnz = -1 ! entries stored in the factors
        real(dp)    :: factor_mem_mb = -1.0_dp ! peak factorization memory
        real(dp)    :: flops = -1.0_dp ! factorization operations
        real(dp)    :: analysis_s = -1.0_dp ! wall time of ordering + symbolic
        real(dp)    :: factor_s = -1.0_dp   ! wall time of the numerical factor
        integer     :: active_threads = -1 ! actual native factor team size
        integer     :: active_device = -2 ! -2 unknown; -1 CPU; >=0 CUDA device
        integer     :: ordering_used = -1 ! actual native ordering option
        integer(i8) :: device_bytes = -1 ! current explicit native CUDA payload
        integer(i8) :: device_peak_bytes = -1_i8
        integer(i8) :: device_budget_required_bytes = 0_i8
        integer(i8) :: device_budget_available_bytes = 0_i8
        logical :: strict_device_memory_supported = .false. ! Explicit payload only
        integer(i8) :: symbolic_nnz = -1 ! unpadded native Cholesky structure
        integer(i8) :: factor_storage_bytes = -1 ! native numeric factor payload
        integer(i8) :: workspace_bytes = -1 ! retained native auxiliary payload
        integer     :: nthreads = -1        ! threads requested
        integer(i8) :: symbolic_peak_bytes = -1_i8
        integer(i8) :: estimated_peak_bytes = -1_i8
        integer(i8) :: numeric_peak_bytes = -1_i8
        integer(i8) :: ipc_bytes = -1_i8
        real(dp) :: effective_pivot_tolerance = -1.0_dp
        logical :: numeric_started = .false.
        integer(i8) :: host_current_bytes = -1_i8, host_peak_bytes = -1_i8
        integer(i8) :: budget_required_bytes = 0_i8, budget_available_bytes = 0_i8
        logical :: strict_memory_supported = .false.
        logical :: compact_memory_used = .false., budget_reduced_threads = .false.
        character(16) :: native_factor_family = ''
        integer :: refinement_steps = 0
        real(dp) :: backward_error = -1.0_dp
        ! Native LU pivoting report (-1 when not reported). The supernodal
        ! counters describe the supernodal attempt, also after a fallback.
        integer :: perturbed_pivots = -1, block_pivot_swaps = -1
        real(dp) :: pivot_growth = -1.0_dp, min_pivot = -1.0_dp
        integer :: lu_method_used = -1 ! FORTSPARSE_LU_SCALAR or _SUPERNODAL
        ! A requested supernodal LU fell back to the scalar LU; the reason is
        ! NONFINITE, PERTURBATION, GROWTH, REFINEMENT, SINGULAR or SOLVE.
        logical :: lu_fallback = .false.
        character(16) :: lu_fallback_reason = ''
        integer :: matching_used = -1 ! FORTSPARSE_MATCHING_*
        ! max(log10 range of the row scale, of the column scale).
        real(dp) :: scaling_log10_range = -1.0_dp
    end type sparse_factor_stats_t

    type, abstract :: sparse_backend_t
        type(sparse_options_t) :: options
        ! Only temporary solver-owned RAW conversion arrays, never caller CSC.
        logical, private :: raw_active = .false., raw_tracked = .false.
        logical, private :: raw_strict = .false.
        integer(i8), private :: raw_payload = 0_i8, raw_peak = 0_i8
        integer(i8), private :: raw_old_live = 0_i8
        integer(i8), private :: raw_required = 0_i8, raw_available = 0_i8
    contains
        procedure :: reset_raw_tracking
        procedure :: apply_raw_tracking
        procedure(factor_real_i), deferred :: factor_real
        procedure(factor_complex_i), deferred :: factor_complex
        procedure :: factor_real_raw => backend_factor_real_raw
        procedure :: refactor_real => backend_refactor_real
        procedure :: refactor_complex => backend_refactor_complex
        procedure :: refactor_real_raw => backend_refactor_real_raw
        procedure :: refactor_complex_raw => backend_refactor_complex_raw
        procedure :: factor_complex_raw => backend_factor_complex_raw
        procedure(solve_real_i), deferred :: solve_real
        procedure(solve_complex_i), deferred :: solve_complex
        procedure(solve_real_inplace_i), deferred :: solve_real_inplace
        procedure(solve_complex_inplace_i), deferred :: solve_complex_inplace
        procedure :: solve_real_block => backend_solve_real_block
        procedure :: solve_complex_block => backend_solve_complex_block
        procedure :: solve_real_op => backend_solve_real_op
        procedure :: solve_complex_op => backend_solve_complex_op
        procedure :: solve_real_block_op => backend_solve_real_block_op
        procedure :: solve_complex_block_op => backend_solve_complex_block_op
        procedure :: factor_stats => backend_factor_stats
        procedure :: copy_from => backend_copy_from
        procedure(vector_i), deferred :: vector
        procedure(free_i), deferred :: free
    end type sparse_backend_t

    abstract interface

        subroutine factor_real_i(self, A, refine, status)
            import :: sparse_backend_t, csc_t, fortsparse_status_t
            class(sparse_backend_t),   intent(inout) :: self
            type(csc_t),               intent(in)    :: A
            logical,                   intent(in)    :: refine
            type(fortsparse_status_t), intent(out)   :: status
        end subroutine factor_real_i

        subroutine factor_complex_i(self, A, refine, status)
            import :: sparse_backend_t, csc_z_t, fortsparse_status_t
            class(sparse_backend_t),   intent(inout) :: self
            type(csc_z_t),             intent(in)    :: A
            logical,                   intent(in)    :: refine
            type(fortsparse_status_t), intent(out)   :: status
        end subroutine factor_complex_i

        subroutine solve_real_i(self, b, x, status)
            import :: sparse_backend_t, dp, fortsparse_status_t
            class(sparse_backend_t),       intent(inout) :: self
            real(dp), target, contiguous,  intent(in)    :: b(:)
            real(dp), target, contiguous,  intent(out)   :: x(:)
            type(fortsparse_status_t),     intent(out)   :: status
        end subroutine solve_real_i

        subroutine solve_complex_i(self, b, x, status)
            import :: sparse_backend_t, dp, fortsparse_status_t
            class(sparse_backend_t),   intent(inout) :: self
            complex(dp),               intent(in)    :: b(:)
            complex(dp),               intent(out)   :: x(:)
            type(fortsparse_status_t), intent(out)   :: status
        end subroutine solve_complex_i

        ! In-place real solve: b holds the right-hand side on entry and the
        ! solution on return. Lets a caller reuse one vector with no temporary
        ! and no extra copy at the call site.
        subroutine solve_real_inplace_i(self, b, status)
            import :: sparse_backend_t, dp, fortsparse_status_t
            class(sparse_backend_t),   intent(inout) :: self
            real(dp),                  intent(inout) :: b(:)
            type(fortsparse_status_t), intent(out)   :: status
        end subroutine solve_real_inplace_i

        ! In-place complex solve; b is the RHS on entry, the solution on return.
        subroutine solve_complex_inplace_i(self, b, status)
            import :: sparse_backend_t, dp, fortsparse_status_t
            class(sparse_backend_t),   intent(inout) :: self
            complex(dp),               intent(inout) :: b(:)
            type(fortsparse_status_t), intent(out)   :: status
        end subroutine solve_complex_inplace_i

        ! Allocate a length-n real solve vector owned by the backend. For the
        ! out-of-process backend it is a slot in the shared mapping, so a solve
        ! whose RHS and solution are such vectors copies nothing across the
        ! boundary; for an in-process backend it is a plain array. The vector
        ! lives until the backend is finalized; the caller does not free it.
        function vector_i(self, n) result(p)
            import :: sparse_backend_t, dp
            class(sparse_backend_t), intent(inout) :: self
            integer,                 intent(in)    :: n
            real(dp), pointer                      :: p(:)
        end function vector_i

        subroutine free_i(self)
            import :: sparse_backend_t
            class(sparse_backend_t), intent(inout) :: self
        end subroutine free_i

    end interface

contains

    subroutine backend_solve_real_op(self, b, x, operation, status)
        class(sparse_backend_t), intent(inout) :: self
        real(dp), target, contiguous, intent(in) :: b(:)
        real(dp), target, contiguous, intent(out) :: x(:)
        integer, intent(in) :: operation
        type(fortsparse_status_t), intent(out) :: status
        if (operation == FORTSPARSE_OP_NORMAL) then
            call self%solve_real(b, x, status)
        else
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                'selected backend does not expose retained transpose solves')
        end if
    end subroutine backend_solve_real_op

    subroutine backend_solve_complex_op(self, b, x, operation, status)
        class(sparse_backend_t), intent(inout) :: self
        complex(dp), intent(in) :: b(:)
        complex(dp), intent(out) :: x(:)
        integer, intent(in) :: operation
        type(fortsparse_status_t), intent(out) :: status
        if (operation == FORTSPARSE_OP_NORMAL) then
            call self%solve_complex(b, x, status)
        else
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                'selected backend does not expose retained adjoint solves')
        end if
    end subroutine backend_solve_complex_op

    subroutine backend_solve_real_block_op(self, b, x, operation, status)
        class(sparse_backend_t), intent(inout) :: self
        real(dp), contiguous, intent(in) :: b(:,:)
        real(dp), contiguous, intent(out) :: x(:,:)
        integer, intent(in) :: operation
        type(fortsparse_status_t), intent(out) :: status
        integer :: j
        if (operation == FORTSPARSE_OP_NORMAL) then
            call self%solve_real_block(b, x, status)
            return
        end if
        call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
            'selected backend does not expose block transpose solves')
        if (any(shape(b) /= shape(x))) return
        ! An empty valid block requires no transpose capability or solve call.
        if (size(b,2) == 0 .and. (operation == FORTSPARSE_OP_TRANSPOSE .or. &
                operation == FORTSPARSE_OP_ADJOINT)) then
            call status_set(status, FORTSPARSE_OK, '')
            return
        end if
        do j = 1, size(b,2)
            call self%solve_real_op(b(:,j), x(:,j), operation, status)
            if (status%code /= 0) return
        end do
    end subroutine backend_solve_real_block_op

    subroutine backend_solve_complex_block_op(self, b, x, operation, status)
        class(sparse_backend_t), intent(inout) :: self
        complex(dp), contiguous, intent(in) :: b(:,:)
        complex(dp), contiguous, intent(out) :: x(:,:)
        integer, intent(in) :: operation
        type(fortsparse_status_t), intent(out) :: status
        integer :: j
        if (operation == FORTSPARSE_OP_NORMAL) then
            call self%solve_complex_block(b, x, status)
            return
        end if
        call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
            'selected backend does not expose block adjoint solves')
        if (any(shape(b) /= shape(x))) return
        ! An empty valid block requires no transpose capability or solve call.
        if (size(b,2) == 0 .and. (operation == FORTSPARSE_OP_TRANSPOSE .or. &
                operation == FORTSPARSE_OP_ADJOINT)) then
            call status_set(status, FORTSPARSE_OK, '')
            return
        end if
        do j = 1, size(b,2)
            call self%solve_complex_op(b(:,j), x(:,j), operation, status)
            if (status%code /= 0) return
        end do
    end subroutine backend_solve_complex_block_op

    ! Solve for the columns of B one by one. A backend with a native multi-RHS
    ! solve overrides this.
    subroutine backend_solve_real_block(self, B, X, status)
        class(sparse_backend_t),   intent(inout) :: self
        real(dp), contiguous,      intent(in)    :: B(:, :)
        real(dp), contiguous,      intent(out)   :: X(:, :)
        type(fortsparse_status_t), intent(out)   :: status

        real(dp), allocatable :: bj(:), xj(:)
        integer :: j

        allocate (bj(size(B, 1)), xj(size(B, 1)))
        do j = 1, size(B, 2)
            bj = B(:, j)
            call self%solve_real(bj, xj, status)
            if (status%code /= 0) return
            X(:, j) = xj
        end do
    end subroutine backend_solve_real_block

    subroutine backend_solve_complex_block(self, B, X, status)
        class(sparse_backend_t),   intent(inout) :: self
        complex(dp), contiguous,   intent(in)    :: B(:, :)
        complex(dp), contiguous,   intent(out)   :: X(:, :)
        type(fortsparse_status_t), intent(out)   :: status

        complex(dp), allocatable :: bj(:), xj(:)
        integer :: j

        allocate (bj(size(B, 1)), xj(size(B, 1)))
        do j = 1, size(B, 2)
            bj = B(:, j)
            call self%solve_complex(bj, xj, status)
            if (status%code /= 0) return
            X(:, j) = xj
        end do
    end subroutine backend_solve_complex_block

    ! Context-copy contract. self is a freshly default-initialized backend of
    ! the dynamic type of src whose options the caller has already set; its
    ! host_memory_limit_bytes is the destination cap. A backend that can
    ! deep-copy its retained state overrides this. The default copies no
    ! factors: self stays unfactored with null handles, so both finalize safely.
    subroutine backend_copy_from(self, src, status)
        class(sparse_backend_t), intent(inout) :: self
        class(sparse_backend_t), intent(in) :: src
        type(fortsparse_status_t), intent(out) :: status

        if (.not. same_type_as(self, src)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'backend copy: source and destination types differ')
            return
        end if
        call status_set(status, FORTSPARSE_NOT_FACTORED, &
            'backend copy: external factors are not copied; destination is unfactored')
    end subroutine backend_copy_from

    ! Default: the backend reports no statistics.
    subroutine backend_factor_stats(self, stats, status)
        class(sparse_backend_t),    intent(inout) :: self
        type(sparse_factor_stats_t), intent(out)  :: stats
        type(fortsparse_status_t),  intent(out)   :: status

        call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
            "factor statistics are not available for this backend")
        stats%factor_nnz = -1
    end subroutine backend_factor_stats

    subroutine reset_raw_tracking(self)
        class(sparse_backend_t), intent(inout) :: self
        ! CSC factor/free may run inside the RAW transaction; preserve its owner.
        if (self%raw_active) return
        self%raw_tracked = .false.
        self%raw_payload = 0_i8
        self%raw_peak = 0_i8
        self%raw_required = 0_i8
        self%raw_available = 0_i8
        self%raw_old_live = 0_i8
        self%raw_strict = .false.
    end subroutine reset_raw_tracking

    subroutine apply_raw_tracking(self, stats)
        class(sparse_backend_t), intent(in) :: self
        type(sparse_factor_stats_t), intent(inout) :: stats
        if (.not. self%raw_tracked) return
        if (self%raw_active) then
            if (stats%host_current_bytes >= 0_i8) &
                stats%host_current_bytes = stats%host_current_bytes + self%raw_payload
            if (stats%host_peak_bytes >= 0_i8) &
                stats%host_peak_bytes = max(self%raw_peak, &
                                            stats%host_peak_bytes + self%raw_payload)
        else
            if (self%raw_peak >= 0_i8) then
                stats%host_peak_bytes = max(stats%host_peak_bytes, self%raw_peak)
            else
                stats%host_peak_bytes = -1_i8
            end if
            stats%strict_memory_supported = self%raw_strict
        end if
        if (self%raw_required > 0_i8) then
            stats%budget_required_bytes = self%raw_required
            stats%budget_available_bytes = self%raw_available
        end if
        stats%factor_mem_mb = -1.0_dp
        if (stats%host_peak_bytes >= 0_i8) &
            stats%factor_mem_mb = real(stats%host_peak_bytes, dp)/1.0e6_dp
    end subroutine apply_raw_tracking

    subroutine begin_raw_conversion(self, nrow, ncol, nz, np, nr, nv, vbytes, &
                                    bytes, status)
        class(sparse_backend_t), intent(inout) :: self
        integer, intent(in) :: nrow, ncol, nz, np, nr, nv
        integer(i8), intent(in) :: vbytes
        integer(i8), intent(out) :: bytes
        type(fortsparse_status_t), intent(out) :: status
        type(sparse_factor_stats_t) :: old
        type(fortsparse_status_t) :: ignored
        type(memory_budget_t) :: budget
        integer(i8) :: indices, values
        call self%reset_raw_tracking()
        bytes = 0_i8
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'raw CSC extent mismatch')
        if (nrow < 0 .or. ncol < 0 .or. nz < 0) return
        if (ncol == huge(ncol) .or. nz == huge(nz)) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                            'raw CSC endpoint overflow')
            return
        end if
        if (np < ncol + 1 .or. nr < nz .or. nv < nz) return
        if (self%options%host_memory_limit_bytes < 0_i8 .or. &
            self%options%device_memory_limit_bytes < 0_i8 .or. &
            self%options%memory_policy < FORTSPARSE_MEMORY_AUTO .or. &
            self%options%memory_policy > FORTSPARSE_MEMORY_LOW) return
        if (self%options%device /= -1) then
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                            'budgeted raw CSC conversion requires native CPU execution')
            return
        end if
        indices = 0_i8
        values = 0_i8
        call memory_checked_mul(int(ncol, i8) + 1_i8 + int(nz, i8), &
                                int(storage_size(0)/8, i8), indices, status)
        if (status%code /= FORTSPARSE_OK) return
        call memory_checked_mul(int(nz, i8), vbytes, values, status)
        if (status%code /= FORTSPARSE_OK) return
        call memory_checked_add(indices, values, bytes, status)
        if (status%code /= FORTSPARSE_OK) return
        call self%factor_stats(old, ignored)
        self%raw_old_live = max(0_i8, old%host_current_bytes)
        self%raw_peak = self%raw_old_live
        self%raw_active = .true.
        self%raw_tracked = .true.
        self%raw_strict = .true.
        budget%limit = self%options%host_memory_limit_bytes
        budget%current = self%raw_old_live
        budget%peak = budget%current
        if (budget%limit > 0_i8 .and. budget%current > budget%limit) then
            self%raw_required = budget%current
            self%raw_available = budget%limit
            call status_set(status, FORTSPARSE_RESOURCE_LIMIT, &
                            'native RAW conversion: retained payload exceeds new cap')
            return
        end if
        call budget%reserve(bytes, status, 'native structured RAW conversion')
        if (status%code /= FORTSPARSE_OK) then
            self%raw_required = budget%failed_required
            self%raw_available = budget%failed_available
            return
        end if
        if (budget%limit > 0_i8) then
            ! Never let an exhausted positive cap turn into zero/unlimited.
            if (budget%limit <= bytes) then
                self%raw_required = 1_i8
                self%raw_available = 0_i8
                call status_set(status, FORTSPARSE_RESOURCE_LIMIT, &
                                'native RAW conversion leaves no factor budget')
            end if
        end if
    end subroutine begin_raw_conversion

    subroutine finish_raw_conversion(self, called, status)
        class(sparse_backend_t), intent(inout) :: self
        logical, intent(in) :: called
        type(fortsparse_status_t), intent(in) :: status
        type(sparse_factor_stats_t) :: stats
        type(fortsparse_status_t) :: ignored
        if (.not. self%raw_active) then
            if (status%code /= FORTSPARSE_OK) call self%free()
            return
        end if
        if (called) then
            call self%factor_stats(stats, ignored)
            self%raw_peak = stats%host_peak_bytes
            self%raw_strict = stats%strict_memory_supported
            self%raw_required = stats%budget_required_bytes
            self%raw_available = stats%budget_available_bytes
        else
            self%raw_peak = max(self%raw_peak, self%raw_old_live + self%raw_payload)
        end if
        if (status%code /= FORTSPARSE_OK) call self%free()
        self%raw_active = .false.
        self%raw_payload = 0_i8
    end subroutine finish_raw_conversion

    subroutine backend_factor_real_raw_budgeted(self, nrow, ncol, nz, ptr, &
                                                rows, values, refine, status)
        class(sparse_backend_t), intent(inout) :: self
        integer, intent(in) :: nrow, ncol, nz, ptr(:), rows(:)
        real(dp), intent(in) :: values(:)
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        logical :: reused
        call converted_real(self, nrow, ncol, nz, ptr, rows, values, refine, status, &
                            .false., reused)
    end subroutine backend_factor_real_raw_budgeted

    subroutine backend_refactor_real_raw_budgeted(self, nrow, ncol, nz, ptr, &
                                                  rows, values, refine, status, reused)
        class(sparse_backend_t), intent(inout) :: self
        integer, intent(in) :: nrow, ncol, nz, ptr(:), rows(:)
        real(dp), intent(in) :: values(:)
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(out) :: reused
        call converted_real(self, nrow, ncol, nz, ptr, rows, values, refine, status, &
                            .true., reused)
    end subroutine backend_refactor_real_raw_budgeted

    subroutine converted_real(self, nrow, ncol, nz, ptr, rows, values, refine, &
                              status, refactor, reused)
        class(sparse_backend_t), intent(inout) :: self
        integer, intent(in) :: nrow, ncol, nz, ptr(:), rows(:)
        real(dp), intent(in) :: values(:)
        logical, intent(in) :: refine, refactor
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(out) :: reused
        type(csc_t) :: a
        integer(i8) :: bytes, saved_cap
        integer :: ierr
        logical :: called
        reused = .false.
        called = .false.
        saved_cap = self%options%host_memory_limit_bytes
        call begin_raw_conversion(self, nrow, ncol, nz, size(ptr), size(rows), &
                                  size(values), 8_i8, bytes, status)
        if (status%code /= FORTSPARSE_OK) goto 900
        allocate (a%col_ptr(ncol + 1), stat=ierr)
        if (ierr /= 0) goto 800
        self%raw_payload = int(storage_size(0)/8, i8)*size(a%col_ptr, kind=i8)
        allocate (a%row_idx(nz), stat=ierr)
        if (ierr /= 0) goto 800
        self%raw_payload = self%raw_payload + int(storage_size(0)/8, i8)*int(nz, i8)
        allocate (a%val(nz), stat=ierr)
        if (ierr /= 0) goto 800
        self%raw_payload = bytes
        self%raw_peak = max(self%raw_peak, self%raw_old_live + bytes)
        a%nrow = nrow
        a%ncol = ncol
        a%nnz = nz
        a%col_ptr = ptr(:ncol + 1)
        a%row_idx = rows(:nz)
        a%val = values(:nz)
        if (saved_cap > 0_i8) self%options%host_memory_limit_bytes = saved_cap - bytes
        called = .true.
        if (refactor) then
            call self%refactor_real(a, refine, status, reused)
        else
            call self%factor_real(a, refine, status)
        end if
        goto 900
800     call status_set(status, FORTSPARSE_OUT_OF_MEMORY, &
                        'native RAW conversion allocation failed')
900     self%options%host_memory_limit_bytes = saved_cap
        ! Snapshot the simultaneous conversion/factor peak before local arrays
        ! leave scope. Failed attempts also invalidate direct backend factors.
        call finish_raw_conversion(self, called, status)
        if (allocated(a%val)) deallocate (a%val)
        if (allocated(a%row_idx)) deallocate (a%row_idx)
        if (allocated(a%col_ptr)) deallocate (a%col_ptr)
    end subroutine converted_real

    subroutine backend_factor_complex_raw_budgeted(self, nrow, ncol, nz, ptr, &
                                                   rows, values, refine, status)
        class(sparse_backend_t), intent(inout) :: self
        integer, intent(in) :: nrow, ncol, nz, ptr(:), rows(:)
        complex(dp), intent(in) :: values(:)
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        logical :: reused
        call converted_complex(self, nrow, ncol, nz, ptr, rows, values, &
                               refine, status, &
                               .false., reused)
    end subroutine backend_factor_complex_raw_budgeted

    subroutine backend_refactor_complex_raw_budgeted(self, nrow, ncol, nz, &
                                                     ptr, rows, values, refine, &
                                                         status, reused)
        class(sparse_backend_t), intent(inout) :: self
        integer, intent(in) :: nrow, ncol, nz, ptr(:), rows(:)
        complex(dp), intent(in) :: values(:)
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(out) :: reused
        call converted_complex(self, nrow, ncol, nz, ptr, rows, values, &
                               refine, status, &
                               .true., reused)
    end subroutine backend_refactor_complex_raw_budgeted

    subroutine converted_complex(self, nrow, ncol, nz, ptr, rows, values, &
                                 refine, status, refactor, reused)
        class(sparse_backend_t), intent(inout) :: self
        integer, intent(in) :: nrow, ncol, nz, ptr(:), rows(:)
        complex(dp), intent(in) :: values(:)
        logical, intent(in) :: refine, refactor
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(out) :: reused
        type(csc_z_t) :: a
        integer(i8) :: bytes, saved_cap
        integer :: ierr
        logical :: called
        reused = .false.
        called = .false.
        saved_cap = self%options%host_memory_limit_bytes
        call begin_raw_conversion(self, nrow, ncol, nz, size(ptr), size(rows), &
                                  size(values), 16_i8, bytes, status)
        if (status%code /= FORTSPARSE_OK) goto 900
        allocate (a%col_ptr(ncol + 1), stat=ierr)
        if (ierr /= 0) goto 800
        self%raw_payload = int(storage_size(0)/8, i8)*size(a%col_ptr, kind=i8)
        allocate (a%row_idx(nz), stat=ierr)
        if (ierr /= 0) goto 800
        self%raw_payload = self%raw_payload + int(storage_size(0)/8, i8)*int(nz, i8)
        allocate (a%val(nz), stat=ierr)
        if (ierr /= 0) goto 800
        self%raw_payload = bytes
        self%raw_peak = max(self%raw_peak, self%raw_old_live + bytes)
        a%nrow = nrow
        a%ncol = ncol
        a%nnz = nz
        a%col_ptr = ptr(:ncol + 1)
        a%row_idx = rows(:nz)
        a%val = values(:nz)
        if (saved_cap > 0_i8) self%options%host_memory_limit_bytes = saved_cap - bytes
        called = .true.
        if (refactor) then
            call self%refactor_complex(a, refine, status, reused)
        else
            call self%factor_complex(a, refine, status)
        end if
        goto 900
800     call status_set(status, FORTSPARSE_OUT_OF_MEMORY, &
                        'native RAW conversion allocation failed')
900     self%options%host_memory_limit_bytes = saved_cap
        ! Snapshot the simultaneous conversion/factor peak before local arrays
        ! leave scope. Failed attempts also invalidate direct backend factors.
        call finish_raw_conversion(self, called, status)
        if (allocated(a%val)) deallocate (a%val)
        if (allocated(a%row_idx)) deallocate (a%row_idx)
        if (allocated(a%col_ptr)) deallocate (a%col_ptr)
    end subroutine converted_complex

    ! Factor directly from caller-owned 1-based CSC arrays. The default wraps
    ! the arrays into a csc_t (one copy) and dispatches to factor_real; a
    ! backend that can stream the arrays into its own storage overrides this to
    ! skip that copy. The real conversion checks requested extents before
    ! slicing; the concrete backend validates CSC contents.
    subroutine backend_factor_real_raw(self, nrow, ncol, nz, col_ptr, &
            row_idx, val, refine, status)
        class(sparse_backend_t),   intent(inout) :: self
        integer,                   intent(in)    :: nrow, ncol, nz
        integer,                   intent(in)    :: col_ptr(:), row_idx(:)
        real(dp),                  intent(in)    :: val(:)
        logical,                   intent(in)    :: refine
        type(fortsparse_status_t), intent(out)   :: status

        type(csc_t) :: A

        call copy_raw_real(nrow, ncol, nz, col_ptr, row_idx, val, A, status)
        if (status%code /= FORTSPARSE_OK) return
        call self%factor_real(A, refine, status)
    end subroutine backend_factor_real_raw

    subroutine copy_raw_real(nrow, ncol, nz, col_ptr, row_idx, val, A, status)
        ! Owned uncapped conversion shared by default/structured raw dispatch.
        integer, intent(in) :: nrow, ncol, nz, col_ptr(:), row_idx(:)
        real(dp), intent(in) :: val(:)
        type(csc_t), intent(out) :: A
        type(fortsparse_status_t), intent(out) :: status
        integer :: ierr
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'raw CSC extent mismatch')
        if (nrow < 0 .or. ncol < 0 .or. nz < 0) return
        if (ncol == huge(ncol)) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, 'raw CSC pointer extent overflow')
            return
        end if
        if (size(col_ptr) < ncol + 1) return
        if (size(row_idx) < nz .or. size(val) < nz) return
        allocate(A%col_ptr(ncol+1), A%row_idx(nz), A%val(nz), stat=ierr)
        if (ierr /= 0) then
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, 'raw CSC conversion allocation failed')
            return
        end if
        A%nrow = nrow
        A%ncol = ncol
        A%nnz = nz
        A%col_ptr = col_ptr(1:ncol+1)
        A%row_idx = row_idx(1:nz)
        A%val = val(1:nz)
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine copy_raw_real

    ! Complex twin of backend_factor_real_raw.
    subroutine backend_factor_complex_raw(self, nrow, ncol, nz, col_ptr, &
            row_idx, val, refine, status)
        class(sparse_backend_t),   intent(inout) :: self
        integer,                   intent(in)    :: nrow, ncol, nz
        integer,                   intent(in)    :: col_ptr(:), row_idx(:)
        complex(dp),               intent(in)    :: val(:)
        logical,                   intent(in)    :: refine
        type(fortsparse_status_t), intent(out)   :: status

        type(csc_z_t) :: A

        A%nrow = nrow
        A%ncol = ncol
        A%nnz = nz
        A%col_ptr = col_ptr(1:ncol + 1)
        A%row_idx = row_idx(1:nz)
        A%val = val(1:nz)
        call self%factor_complex(A, refine, status)
    end subroutine backend_factor_complex_raw

    ! Refactor: new values, ideally the same sparsity pattern as the last
    ! factorization. A backend that can keep its symbolic analysis overrides
    ! these and sets reused = .true. when it did; the default is a full factor.
    subroutine backend_refactor_real(self, A, refine, status, reused)
        class(sparse_backend_t),   intent(inout) :: self
        type(csc_t),               intent(in)    :: A
        logical,                   intent(in)    :: refine
        type(fortsparse_status_t), intent(out)   :: status
        logical,                   intent(out)   :: reused

        reused = .false.
        call self%factor_real(A, refine, status)
    end subroutine backend_refactor_real

    subroutine backend_refactor_complex(self, A, refine, status, reused)
        class(sparse_backend_t),   intent(inout) :: self
        type(csc_z_t),             intent(in)    :: A
        logical,                   intent(in)    :: refine
        type(fortsparse_status_t), intent(out)   :: status
        logical,                   intent(out)   :: reused

        reused = .false.
        call self%factor_complex(A, refine, status)
    end subroutine backend_refactor_complex

    subroutine backend_refactor_real_raw(self, nrow, ncol, nz, col_ptr, &
            row_idx, val, refine, status, reused)
        class(sparse_backend_t),   intent(inout) :: self
        integer,                   intent(in)    :: nrow, ncol, nz
        integer,                   intent(in)    :: col_ptr(:), row_idx(:)
        real(dp),                  intent(in)    :: val(:)
        logical,                   intent(in)    :: refine
        type(fortsparse_status_t), intent(out)   :: status
        logical,                   intent(out)   :: reused

        type(csc_t) :: A

        reused = .false.
        call copy_raw_real(nrow, ncol, nz, col_ptr, row_idx, val, A, status)
        if (status%code /= FORTSPARSE_OK) return
        call self%refactor_real(A, refine, status, reused)
    end subroutine backend_refactor_real_raw

    subroutine backend_refactor_complex_raw(self, nrow, ncol, nz, col_ptr, &
            row_idx, val, refine, status, reused)
        class(sparse_backend_t),   intent(inout) :: self
        integer,                   intent(in)    :: nrow, ncol, nz
        integer,                   intent(in)    :: col_ptr(:), row_idx(:)
        complex(dp),               intent(in)    :: val(:)
        logical,                   intent(in)    :: refine
        type(fortsparse_status_t), intent(out)   :: status
        logical,                   intent(out)   :: reused

        type(csc_z_t) :: A

        A%nrow = nrow
        A%ncol = ncol
        A%nnz = nz
        A%col_ptr = col_ptr(1:ncol + 1)
        A%row_idx = row_idx(1:nz)
        A%val = val(1:nz)
        call self%refactor_complex(A, refine, status, reused)
    end subroutine backend_refactor_complex_raw

end module fortsparse_backend
