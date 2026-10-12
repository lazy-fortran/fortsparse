module fortsparse_backend
    ! Abstract sparse-solver backend. A concrete backend retains a
    ! factorization across solves and maps its native error codes onto the
    ! fortsparse status codes. The solver holds one of these polymorphically and
    ! dispatches factor/solve/free through the deferred bindings, so adding a
    ! backend never touches client code or the solver dispatch logic.
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csc, only: csc_t, csc_z_t
    use fortsparse_status, only: fortsparse_status_t, status_set, &
        FORTSPARSE_BACKEND_UNAVAILABLE
    implicit none
    private

    public :: sparse_backend_t
    public :: sparse_options_t
    public :: sparse_factor_stats_t

    ! Matrix structure hint for a factorization (sparse_options_t%symmetry).
    ! SPD and INDEFINITE factor only the lower triangle of the supplied full
    ! symmetric CSC matrix; GENERAL factors the matrix as given.
    integer, parameter, public :: FORTSPARSE_SYM_GENERAL = 0
    integer, parameter, public :: FORTSPARSE_SYM_SPD = 1
    integer, parameter, public :: FORTSPARSE_SYM_INDEFINITE = 2
    ! Complex only: A = A^H, lower triangle read (a real matrix treats it as
    ! INDEFINITE). For complex matrices SPD means Hermitian positive definite,
    ! and INDEFINITE means complex symmetric, A = A^T (not A^H).
    integer, parameter, public :: FORTSPARSE_SYM_HERMITIAN = 3

    ! Fill-reducing ordering (sparse_options_t%ordering).
    integer, parameter, public :: FORTSPARSE_ORDER_DEFAULT = 0
    integer, parameter, public :: FORTSPARSE_ORDER_AMD = 1
    integer, parameter, public :: FORTSPARSE_ORDER_ND = 2

    ! Backend tuning, copied from the solver to the backend at every factor.
    ! Backends that have no use for a field ignore it.
    type :: sparse_options_t
        integer :: symmetry = FORTSPARSE_SYM_GENERAL
        integer :: ordering = FORTSPARSE_ORDER_DEFAULT
        integer :: nthreads = 0 ! 0: backend default (OMP_NUM_THREADS)
        logical :: verbose = .false.
    end type sparse_options_t

    ! Factorization statistics; -1 / -1.0 mark a quantity the backend does not
    ! report.
    type :: sparse_factor_stats_t
        integer(i8) :: factor_nnz = -1      ! entries stored in the factors
        real(dp)    :: factor_mem_mb = -1.0_dp ! peak factorization memory
        real(dp)    :: flops = -1.0_dp      ! factorization operations
        real(dp)    :: analysis_s = -1.0_dp ! wall time of ordering + symbolic
        real(dp)    :: factor_s = -1.0_dp   ! wall time of the numerical factor
        integer     :: nthreads = -1        ! threads requested
    end type sparse_factor_stats_t

    type, abstract :: sparse_backend_t
        type(sparse_options_t) :: options
    contains
        procedure(factor_real_i), deferred :: factor_real
        procedure(factor_complex_i), deferred :: factor_complex
        procedure :: factor_real_raw => backend_factor_real_raw
        procedure :: factor_complex_raw => backend_factor_complex_raw
        procedure(solve_real_i), deferred :: solve_real
        procedure(solve_complex_i), deferred :: solve_complex
        procedure(solve_real_inplace_i), deferred :: solve_real_inplace
        procedure(solve_complex_inplace_i), deferred :: solve_complex_inplace
        procedure :: solve_real_block => backend_solve_real_block
        procedure :: solve_complex_block => backend_solve_complex_block
        procedure :: factor_stats => backend_factor_stats
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

    ! Default: the backend reports no statistics.
    subroutine backend_factor_stats(self, stats, status)
        class(sparse_backend_t),    intent(inout) :: self
        type(sparse_factor_stats_t), intent(out)  :: stats
        type(fortsparse_status_t),  intent(out)   :: status

        call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
            "factor statistics are not available for this backend")
        stats%factor_nnz = -1
    end subroutine backend_factor_stats

    ! Factor directly from caller-owned 1-based CSC arrays. The default wraps
    ! the arrays into a csc_t (one copy) and dispatches to factor_real; a
    ! backend that can stream the arrays into its own storage overrides this to
    ! skip that copy. The arrays are trusted as-is: col_ptr(ncol+1),
    ! row_idx(nz), val(nz).
    subroutine backend_factor_real_raw(self, nrow, ncol, nz, col_ptr, &
            row_idx, val, refine, status)
        class(sparse_backend_t),   intent(inout) :: self
        integer,                   intent(in)    :: nrow, ncol, nz
        integer,                   intent(in)    :: col_ptr(:), row_idx(:)
        real(dp),                  intent(in)    :: val(:)
        logical,                   intent(in)    :: refine
        type(fortsparse_status_t), intent(out)   :: status

        type(csc_t) :: A

        A%nrow = nrow
        A%ncol = ncol
        A%nnz = nz
        A%col_ptr = col_ptr(1:ncol + 1)
        A%row_idx = row_idx(1:nz)
        A%val = val(1:nz)
        call self%factor_real(A, refine, status)
    end subroutine backend_factor_real_raw

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

end module fortsparse_backend
