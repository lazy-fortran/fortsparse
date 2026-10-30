module fortsparse_solver
    ! Public, backend-pluggable sparse direct solver API.
    !
    ! A client factors a matrix once, solves any number of right-hand sides
    ! reusing that factorization, then frees it. The backend is selected by an
    ! integer tag on the solver. The default, FORTSPARSE_BACKEND_AUTO, picks
    ! SuperLU (BSD) or, for large matrices with a ready MUMPS runtime, MUMPS
    ! per factorization (sparse_auto_backend); UMFPACK runs out-of-process
    ! through a separate GPL helper and is reached by selecting
    ! FORTSPARSE_BACKEND_UMFPACK_IPC. Selecting an unbuilt or unknown backend
    ! yields a documented
    ! FORTSPARSE_BACKEND_UNAVAILABLE status rather than a crash.
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated
    use fortsparse_kinds, only: dp, i8
    use fortsparse_status, only: fortsparse_status_t, status_set, status_ok, &
        FORTSPARSE_OK, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_INVALID_MATRIX
    use fortsparse_csc, only: csc_t, csc_z_t, csc_is_valid, csc_matvec
    use fortsparse_backend, only: sparse_backend_t, sparse_options_t, &
        sparse_factor_stats_t, &
        FORTSPARSE_SYM_GENERAL, FORTSPARSE_SYM_SPD, FORTSPARSE_SYM_INDEFINITE, &
        FORTSPARSE_SYM_HERMITIAN, &
        FORTSPARSE_ORDER_DEFAULT, FORTSPARSE_ORDER_AMD, FORTSPARSE_ORDER_ND, &
        FORTSPARSE_ORDER_NATURAL, &
        FORTSPARSE_MEMORY_AUTO, FORTSPARSE_MEMORY_LOW, &
        FORTSPARSE_OP_NORMAL, FORTSPARSE_OP_TRANSPOSE, FORTSPARSE_OP_ADJOINT, &
        FORTSPARSE_LU_AUTO, FORTSPARSE_LU_SUPERNODAL
    use fortsparse_inertia, only: FORTSPARSE_INERTIA_UNAVAILABLE
#ifdef FORTSPARSE_HAVE_MUMPS
    use fortsparse_mumps, only: mumps_backend_t, mumps_available
#endif
#ifdef FORTSPARSE_HAVE_SUPERLU
    use fortsparse_superlu, only: superlu_backend_t
#endif
    use fortsparse_umfpack_ipc, only: umfpack_ipc_backend_t
    use fortsparse_memory, only: sparse_memory_plan_t, memory_budget_t, &
        plan_natural_spd_memory
    use fortsparse_spd_leaves, only: spd_update
    use fortsparse_native_general, only: native_general_backend_t
    use fortsparse_lu_leaves, only: lu_z_update
    use fortsparse_native, only: native_backend_t, native_plan_ordered_spd, &
        bounded_default_prefers_amd
    implicit none
    private

    ! Backend identifiers. FORTSPARSE_BACKEND_DEFAULT below names the default
    ! in-process backend; the UMFPACK_IPC id selects the out-of-process GPL
    ! helper. KLU and PARDISO are reserved for later releases.
    integer, parameter, public :: FORTSPARSE_BACKEND_SUPERLU = 1
    integer, parameter, public :: FORTSPARSE_BACKEND_UMFPACK_IPC = 2
    ! In-process threaded MUMPS: Cholesky/LDL^T/LU, nested dissection, block
    ! solves (present only when built with FORTSPARSE_ENABLE_MUMPS).
    integer, parameter, public :: FORTSPARSE_BACKEND_MUMPS = 3
    integer, parameter, public :: FORTSPARSE_BACKEND_NATIVE = 4
    integer, parameter, public :: FORTSPARSE_BACKEND_NATIVE_REFERENCE = 5
    ! Explicit approximate inverse; excluded from AUTO and exact derivatives.
    integer, parameter, public :: FORTSPARSE_BACKEND_SUPERLU_ILU = 6
    ! Per-matrix choice between the in-process backends at factor time
    ! (sparse_auto_backend); the default of sparse_solver_t%backend_id and
    ! hence of sparse_solve_once.
    integer, parameter, public :: FORTSPARSE_BACKEND_AUTO = 0
    integer, parameter, public :: FORTSPARSE_BACKEND_DEFAULT = &
        FORTSPARSE_BACKEND_AUTO

    ! Fewest stored entries for which AUTO factors with MUMPS (measured
    ! crossover of factor + solve wall time; README, "Backend selection").
    integer, parameter :: AUTO_MUMPS_MIN_NNZ = 40000
    ! Backends AUTO may choose; MPI readiness is checked at factor time.
#ifdef FORTSPARSE_HAVE_MUMPS
    logical, parameter :: AUTO_HAS_MUMPS = .true.
#else
    logical, parameter :: AUTO_HAS_MUMPS = .false.
#endif
#ifdef FORTSPARSE_HAVE_ZMUMPS
    logical, parameter :: AUTO_HAS_ZMUMPS = .true.
#else
    logical, parameter :: AUTO_HAS_ZMUMPS = .false.
#endif
#ifdef FORTSPARSE_HAVE_SUPERLU
    logical, parameter :: AUTO_HAS_SUPERLU = .true.
#else
    logical, parameter :: AUTO_HAS_SUPERLU = .false.
#endif

    ! Re-exported option and statistics vocabulary.
    public :: sparse_options_t, sparse_factor_stats_t
    public :: FORTSPARSE_MEMORY_AUTO, FORTSPARSE_MEMORY_LOW
    public :: FORTSPARSE_OP_NORMAL, FORTSPARSE_OP_TRANSPOSE, FORTSPARSE_OP_ADJOINT
    public :: sparse_solve_op, sparse_solve_block_op
    public :: FORTSPARSE_SYM_GENERAL, FORTSPARSE_SYM_SPD, &
        FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_SYM_HERMITIAN
    public :: FORTSPARSE_ORDER_DEFAULT, FORTSPARSE_ORDER_AMD, FORTSPARSE_ORDER_ND, &
        FORTSPARSE_ORDER_NATURAL

    ! Solver handle. Holds backend selection, refinement toggle, the factored
    ! flag, and the polymorphic backend that owns the retained factorization.
    ! Assignment dst = src is a deep copy (sparse_copy); since it cannot
    ! return a status, a failed copy leaves dst a valid unfactored context.
    type, public :: sparse_solver_t
        integer :: backend_id = FORTSPARSE_BACKEND_DEFAULT
        logical :: refine = .false.
        logical :: factored = .false.
        ! Structure hint, ordering and thread count, applied at the next factor.
        type(sparse_options_t) :: options
        class(sparse_backend_t), allocatable :: backend
    contains
        procedure, private :: assign_solver
        generic :: assignment(=) => assign_solver
    end type sparse_solver_t

    public :: sparse_factor, sparse_analyze, sparse_memory_plan_t
    public :: sparse_refactor
    public :: sparse_solve
    public :: sparse_free
    public :: sparse_destroy
    public :: sparse_copy
    public :: sparse_solve_once
    public :: sparse_solve_jvp
    public :: sparse_solve_vjp
    public :: sparse_vector
    public :: sparse_solve_block
    public :: sparse_factor_stats
    public :: sparse_inertia
    public :: sparse_auto_backend

    ! Refactor after the matrix values changed. When the sparsity pattern (and
    ! the symmetry and ordering options) equal those of the last factorization
    ! and the backend supports it (MUMPS), the symbolic analysis is kept and
    ! only the numerical factorization is redone; otherwise this is a full
    ! sparse_factor. The optional reused is .true. when the analysis was kept.
    ! A numerical failure (singular) leaves the solver unfactored.
    interface sparse_refactor
        module procedure sparse_refactor_real
        module procedure sparse_refactor_complex
        module procedure sparse_refactor_real_raw
        module procedure sparse_refactor_complex_raw
    end interface sparse_refactor

    ! Factor a real or complex matrix into the solver handle. The csc_t forms
    ! take an assembled matrix; the raw forms take caller-owned 1-based CSC
    ! arrays directly, which lets a backend stream them into its own storage
    ! without an intermediate csc_t copy: for a large system that copy is the
    ! caller's peak-memory overhead.
    interface sparse_factor
        module procedure sparse_factor_real
        module procedure sparse_factor_complex
        module procedure sparse_factor_real_raw
        module procedure sparse_factor_complex_raw
    end interface sparse_factor

    ! Solve A x = b reusing the factorization; real or complex vectors. The
    ! two-vector forms write the solution to x; the in-place forms return it in
    ! b, sparing the caller a temporary and a copy in a tight solve loop.
    interface sparse_solve
        module procedure sparse_solve_real
        module procedure sparse_solve_complex
        module procedure sparse_solve_real_inplace
        module procedure sparse_solve_complex_inplace
    end interface sparse_solve

    ! Convenience: factor, solve one RHS, and free in a single call.
    interface sparse_solve_block
        module procedure sparse_solve_block_real
        module procedure sparse_solve_block_complex
    end interface sparse_solve_block

    interface sparse_solve_op
        module procedure sparse_solve_real_op
        module procedure sparse_solve_complex_op
    end interface sparse_solve_op
    interface sparse_solve_block_op
        module procedure sparse_solve_real_block_op
        module procedure sparse_solve_complex_block_op
    end interface sparse_solve_block_op

    interface sparse_solve_once
        module procedure sparse_solve_once_real
        module procedure sparse_solve_once_complex
    end interface sparse_solve_once

    interface sparse_solve_jvp
        module procedure sparse_solve_jvp_real
        module procedure sparse_solve_jvp_complex
    end interface sparse_solve_jvp

    interface sparse_solve_vjp
        module procedure sparse_solve_vjp_real
        module procedure sparse_solve_vjp_complex
    end interface sparse_solve_vjp

contains

    subroutine sparse_solve_real_op(solver, b, x, operation, status)
        type(sparse_solver_t), intent(inout) :: solver
        real(dp), target, contiguous, intent(in) :: b(:)
        real(dp), target, contiguous, intent(out) :: x(:)
        integer, intent(in) :: operation
        type(fortsparse_status_t), intent(out) :: status
        if (.not. solver%factored) then
            call not_factored(status)
            return
        end if
        if (operation < FORTSPARSE_OP_NORMAL .or. operation > FORTSPARSE_OP_ADJOINT) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, 'invalid solve operation')
            return
        end if
        if (any(shape(b) /= shape(x))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, 'incompatible solve dimensions')
            return
        end if
        call solver%backend%solve_real_op(b, x, operation, status)
    end subroutine sparse_solve_real_op

    subroutine sparse_solve_real_block_op(solver, b, x, operation, status)
        type(sparse_solver_t), intent(inout) :: solver
        real(dp), contiguous, intent(in) :: b(:,:)
        real(dp), contiguous, intent(out) :: x(:,:)
        integer, intent(in) :: operation
        type(fortsparse_status_t), intent(out) :: status
        if (.not. solver%factored) then
            call not_factored(status)
            return
        end if
        if (operation < FORTSPARSE_OP_NORMAL .or. operation > FORTSPARSE_OP_ADJOINT) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, 'invalid solve operation')
            return
        end if
        if (any(shape(b) /= shape(x))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, 'incompatible solve dimensions')
            return
        end if
        call solver%backend%solve_real_block_op(b, x, operation, status)
    end subroutine sparse_solve_real_block_op

    subroutine sparse_solve_complex_op(solver, b, x, operation, status)
        type(sparse_solver_t), intent(inout) :: solver
        complex(dp), intent(in) :: b(:)
        complex(dp), intent(out) :: x(:)
        integer, intent(in) :: operation
        type(fortsparse_status_t), intent(out) :: status
        if (.not. solver%factored) then
            call not_factored(status)
            return
        end if
        if (operation < FORTSPARSE_OP_NORMAL .or. operation > FORTSPARSE_OP_ADJOINT) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, 'invalid solve operation')
            return
        end if
        if (any(shape(b) /= shape(x))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, 'incompatible solve dimensions')
            return
        end if
        call solver%backend%solve_complex_op(b, x, operation, status)
    end subroutine sparse_solve_complex_op

    subroutine sparse_solve_complex_block_op(solver, b, x, operation, status)
        type(sparse_solver_t), intent(inout) :: solver
        complex(dp), contiguous, intent(in) :: b(:,:)
        complex(dp), contiguous, intent(out) :: x(:,:)
        integer, intent(in) :: operation
        type(fortsparse_status_t), intent(out) :: status
        if (.not. solver%factored) then
            call not_factored(status)
            return
        end if
        if (operation < FORTSPARSE_OP_NORMAL .or. operation > FORTSPARSE_OP_ADJOINT) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, 'invalid solve operation')
            return
        end if
        if (any(shape(b) /= shape(x))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, 'incompatible solve dimensions')
            return
        end if
        call solver%backend%solve_complex_block_op(b, x, operation, status)
    end subroutine sparse_solve_complex_block_op

    subroutine sparse_analyze(solver, A, plan, status, nrhs)
        ! Plans the bounded compact natural/AMD CPU SPD layout without
        ! factorizing or changing retained factors; DEFAULT plans the order the
        ! factor selects. Other schedules are explicit.
        type(sparse_solver_t), intent(in) :: solver
        type(csc_t), intent(in) :: A
        type(sparse_memory_plan_t), intent(out) :: plan
        type(fortsparse_status_t), intent(out) :: status
        integer, optional, intent(in) :: nrhs
        type(memory_budget_t) :: budget
        type(sparse_options_t) :: ordered
        plan = sparse_memory_plan_t()
        if (solver%backend_id /= FORTSPARSE_BACKEND_NATIVE .and. &
            solver%backend_id /= FORTSPARSE_BACKEND_NATIVE_REFERENCE) then
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                'memory analysis: native CPU SPD reference layout only')
            return
        end if
        if (solver%options%symmetry /= FORTSPARSE_SYM_SPD .or. &
            (solver%options%ordering /= FORTSPARSE_ORDER_DEFAULT .and. &
             solver%options%ordering /= FORTSPARSE_ORDER_NATURAL .and. &
             solver%options%ordering /= FORTSPARSE_ORDER_AMD) .or. &
            solver%options%structural_matching .or. solver%options%device /= -1) then
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                'memory analysis: require compact natural/AMD CPU SPD layout')
            return
        end if
        if (solver%backend_id == FORTSPARSE_BACKEND_NATIVE .and. &
            solver%options%memory_policy /= FORTSPARSE_MEMORY_LOW) then
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                'memory analysis: select LOW_MEMORY for the exact reference layout')
            return
        end if
        if (solver%options%device_memory_limit_bytes < 0_i8 .or. &
            solver%options%memory_policy < FORTSPARSE_MEMORY_AUTO .or. &
            solver%options%memory_policy > FORTSPARSE_MEMORY_LOW) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'memory analysis: invalid memory limit or policy')
            return
        end if
        if (solver%options%ordering == FORTSPARSE_ORDER_AMD .and. &
            solver%options%memory_policy /= FORTSPARSE_MEMORY_LOW) then
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                'memory analysis: AMD planning requires LOW_MEMORY audited ordering')
            return
        end if
        budget%limit = solver%options%host_memory_limit_bytes
        if (allocated(solver%backend)) then
            select type (backend => solver%backend)
            type is (native_backend_t)
                if (.not. backend%bounded_memory) then
                    call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                        'memory analysis: release unbudgeted factors before planning')
                    return
                end if
                budget%current = backend%host_budget%current
                budget%peak = budget%current
            class default
                call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                    'memory analysis: release other factors before strict planning')
                return
            end select
        end if
        if (solver%options%ordering == FORTSPARSE_ORDER_AMD) then
            call native_plan_ordered_spd(A, solver%options, solver%refine, &
                                         budget, plan, status, nrhs)
        else if (bounded_default_prefers_amd(solver%options, &
                solver%backend_id == FORTSPARSE_BACKEND_NATIVE_REFERENCE, A)) then
            ! Plan the AMD order that the production DEFAULT factor selects.
            ordered = solver%options
            ordered%ordering = FORTSPARSE_ORDER_AMD
            call native_plan_ordered_spd(A, ordered, solver%refine, &
                                         budget, plan, status, nrhs)
        else
            call plan_natural_spd_memory(A, plan, status, budget, nrhs, &
                refine=solver%refine)
        end if
    end subroutine sparse_analyze

    ! Factor a real matrix. Ensures the concrete backend, then dispatches.
    subroutine sparse_factor_real(solver, A, status)
        type(sparse_solver_t),     intent(inout) :: solver
        type(csc_t),               intent(in)    :: A
        type(fortsparse_status_t), intent(out)   :: status

        call ensure_backend(solver, A%nnz, .false., status)
        if (.not. status_ok(status)) return
        call solver%backend%factor_real(A, solver%refine, status)
        solver%factored = status_ok(status)
    end subroutine sparse_factor_real

    ! Factor a complex matrix. Ensures the concrete backend, then dispatches.
    subroutine sparse_factor_complex(solver, A, status)
        type(sparse_solver_t),     intent(inout) :: solver
        type(csc_z_t),             intent(in)    :: A
        type(fortsparse_status_t), intent(out)   :: status

        call ensure_backend(solver, A%nnz, .true., status)
        if (.not. status_ok(status)) return
        call solver%backend%factor_complex(A, solver%refine, status)
        solver%factored = status_ok(status)
    end subroutine sparse_factor_complex

    ! Factor a real matrix given as raw 1-based CSC arrays: col_ptr(ncol+1),
    ! row_idx(nz), val(nz).
    subroutine sparse_factor_real_raw(solver, nrow, ncol, nz, col_ptr, &
            row_idx, val, status)
        type(sparse_solver_t),     intent(inout) :: solver
        integer,                   intent(in)    :: nrow, ncol, nz
        integer,                   intent(in)    :: col_ptr(:), row_idx(:)
        real(dp),                  intent(in)    :: val(:)
        type(fortsparse_status_t), intent(out)   :: status

        solver%factored = .false.
        call ensure_backend(solver, nz, .false., status)
        if (.not. status_ok(status)) return
        call solver%backend%factor_real_raw(nrow, ncol, nz, col_ptr, row_idx, &
            val, solver%refine, status)
        solver%factored = status_ok(status)
    end subroutine sparse_factor_real_raw

    ! Factor a complex matrix given as raw 1-based CSC arrays.
    subroutine sparse_factor_complex_raw(solver, nrow, ncol, nz, col_ptr, &
            row_idx, val, status)
        type(sparse_solver_t),     intent(inout) :: solver
        integer,                   intent(in)    :: nrow, ncol, nz
        integer,                   intent(in)    :: col_ptr(:), row_idx(:)
        complex(dp),               intent(in)    :: val(:)
        type(fortsparse_status_t), intent(out)   :: status

        solver%factored = .false.
        call ensure_backend(solver, nz, .true., status)
        if (.not. status_ok(status)) return
        call solver%backend%factor_complex_raw(nrow, ncol, nz, col_ptr, &
            row_idx, val, solver%refine, status)
        solver%factored = status_ok(status)
    end subroutine sparse_factor_complex_raw

    subroutine sparse_refactor_real(solver, A, status, reused)
        type(sparse_solver_t),     intent(inout) :: solver
        type(csc_t),               intent(in)    :: A
        type(fortsparse_status_t), intent(out)   :: status
        logical, optional,         intent(out)   :: reused

        logical :: r

        r = .false.
        call ensure_backend(solver, A%nnz, .false., status)
        if (status_ok(status)) then
            call solver%backend%refactor_real(A, solver%refine, status, r)
            solver%factored = status_ok(status)
        end if
        if (present(reused)) reused = r
    end subroutine sparse_refactor_real

    subroutine sparse_refactor_complex(solver, A, status, reused)
        type(sparse_solver_t),     intent(inout) :: solver
        type(csc_z_t),             intent(in)    :: A
        type(fortsparse_status_t), intent(out)   :: status
        logical, optional,         intent(out)   :: reused

        logical :: r

        r = .false.
        call ensure_backend(solver, A%nnz, .true., status)
        if (status_ok(status)) then
            call solver%backend%refactor_complex(A, solver%refine, status, r)
            solver%factored = status_ok(status)
        end if
        if (present(reused)) reused = r
    end subroutine sparse_refactor_complex

    subroutine sparse_refactor_real_raw(solver, nrow, ncol, nz, col_ptr, &
            row_idx, val, status, reused)
        type(sparse_solver_t),     intent(inout) :: solver
        integer,                   intent(in)    :: nrow, ncol, nz
        integer,                   intent(in)    :: col_ptr(:), row_idx(:)
        real(dp),                  intent(in)    :: val(:)
        type(fortsparse_status_t), intent(out)   :: status
        logical, optional,         intent(out)   :: reused

        logical :: r

        r = .false.
        solver%factored = .false.
        call ensure_backend(solver, nz, .false., status)
        if (status_ok(status)) then
            call solver%backend%refactor_real_raw(nrow, ncol, nz, col_ptr, &
                row_idx, val, solver%refine, status, r)
            solver%factored = status_ok(status)
        end if
        if (present(reused)) reused = r
    end subroutine sparse_refactor_real_raw

    subroutine sparse_refactor_complex_raw(solver, nrow, ncol, nz, col_ptr, &
            row_idx, val, status, reused)
        type(sparse_solver_t),     intent(inout) :: solver
        integer,                   intent(in)    :: nrow, ncol, nz
        integer,                   intent(in)    :: col_ptr(:), row_idx(:)
        complex(dp),               intent(in)    :: val(:)
        type(fortsparse_status_t), intent(out)   :: status
        logical, optional,         intent(out)   :: reused

        logical :: r

        r = .false.
        solver%factored = .false.
        call ensure_backend(solver, nz, .true., status)
        if (status_ok(status)) then
            call solver%backend%refactor_complex_raw(nrow, ncol, nz, col_ptr, &
                row_idx, val, solver%refine, status, r)
            solver%factored = status_ok(status)
        end if
        if (present(reused)) reused = r
    end subroutine sparse_refactor_complex_raw

    ! Solve A x = b for a real RHS, reusing the stored factorization.
    subroutine sparse_solve_real(solver, b, x, status)
        type(sparse_solver_t),        intent(inout) :: solver
        real(dp), target, contiguous, intent(in)    :: b(:)
        real(dp), target, contiguous, intent(out)   :: x(:)
        type(fortsparse_status_t),    intent(out)   :: status

        if (.not. solver%factored) then
            call not_factored(status)
            return
        end if
        call solver%backend%solve_real(b, x, status)
    end subroutine sparse_solve_real

    ! Solve A x = b for a complex RHS, reusing the stored factorization.
    subroutine sparse_solve_complex(solver, b, x, status)
        type(sparse_solver_t),     intent(inout) :: solver
        complex(dp),               intent(in)    :: b(:)
        complex(dp),               intent(out)   :: x(:)
        type(fortsparse_status_t), intent(out)   :: status

        if (.not. solver%factored) then
            call not_factored(status)
            return
        end if
        call solver%backend%solve_complex(b, x, status)
    end subroutine sparse_solve_complex

    subroutine sparse_solve_jvp_real(solver, A_dot, x, b_dot, x_dot, status)
        ! Exact implicit tangent for A x = b:
        ! A x_dot = b_dot - A_dot x.
        !
        ! The sparsity pattern may be fixed or changing; A_dot is an ordinary
        ! CSC matrix. The retained primal factorization is reused.
        type(sparse_solver_t), intent(inout) :: solver
        type(csc_t), intent(in) :: A_dot
        real(dp), target, contiguous, intent(in) :: x(:), b_dot(:)
        real(dp), target, contiguous, intent(out) :: x_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: column, entry, row

        if (.not. solver%factored) then
            call not_factored(status)
            return
        end if
        if (has_approximate_factors(solver)) then
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                "implicit derivatives require exact factors; SuperLU ILU is approximate")
            return
        end if
        if (.not. csc_is_valid(A_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "sparse_solve_jvp: invalid matrix tangent")
            return
        end if
        if (A_dot%nrow /= size(b_dot) .or. A_dot%ncol /= size(x)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "sparse_solve_jvp: incompatible tangent dimensions")
            return
        end if
        if (size(x_dot) /= size(b_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "sparse_solve_jvp: incompatible output dimension")
            return
        end if
        x_dot = b_dot
        do column = 1, A_dot%ncol
            do entry = A_dot%col_ptr(column), A_dot%col_ptr(column+1)-1
                row = A_dot%row_idx(entry)
                x_dot(row) = spd_update(x_dot(row), A_dot%val(entry), x(column))
            end do
        end do
        call sparse_solve_real_inplace(solver, x_dot, status)
    end subroutine sparse_solve_jvp_real

    subroutine sparse_solve_vjp_real( &
            transpose_solver, A, x, x_bar, b_bar, A_values_bar, status, use_primal_factors)
        ! Exact implicit adjoint for A x = b:
        ! A^T b_bar = x_bar, (A_ij)_bar = -b_bar_i*x_j.
        !
        ! transpose_solver must retain a factorization of A^T. Returning only
        ! active CSC values preserves the caller's sparsity contract and avoids
        ! materializing a dense matrix cotangent.
        type(sparse_solver_t), intent(inout) :: transpose_solver
        type(csc_t), intent(in) :: A
        real(dp), intent(in) :: x(:)
        real(dp), target, contiguous, intent(in) :: x_bar(:)
        real(dp), target, contiguous, intent(out) :: b_bar(:)
        real(dp), intent(out) :: A_values_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        logical, optional, intent(in) :: use_primal_factors
        logical :: primal
        integer :: column, entry

        if (.not. transpose_solver%factored) then
            call not_factored(status)
            return
        end if
        if (has_approximate_factors(transpose_solver)) then
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                "implicit derivatives require exact factors; SuperLU ILU is approximate")
            return
        end if
        if (.not. csc_is_valid(A)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "sparse_solve_vjp: invalid primal matrix")
            return
        end if
        if (A%ncol /= size(x) .or. A%ncol /= size(x_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "sparse_solve_vjp: incompatible primal dimensions")
            return
        end if
        if (A%nrow /= size(b_bar) .or. A%nnz /= size(A_values_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "sparse_solve_vjp: incompatible cotangent dimensions")
            return
        end if
        primal = .false.
        if (present(use_primal_factors)) primal = use_primal_factors
        if (primal) then
            call sparse_solve_real_op(transpose_solver, x_bar, b_bar, FORTSPARSE_OP_TRANSPOSE, status)
        else
        call sparse_solve_real(transpose_solver, x_bar, b_bar, status)
        end if
        if (.not. status_ok(status)) return
        do column = 1, A%ncol
            do entry = A%col_ptr(column), A%col_ptr(column + 1) - 1
                A_values_bar(entry) = -b_bar(A%row_idx(entry))*x(column)
            end do
        end do
    end subroutine sparse_solve_vjp_real

    subroutine sparse_solve_jvp_complex( &
            solver, A_dot, x, b_dot, x_dot, status)
        ! Exact complex implicit tangent for A x = b:
        ! A x_dot = b_dot - A_dot x.
        type(sparse_solver_t), intent(inout) :: solver
        type(csc_z_t), intent(in) :: A_dot
        complex(dp), intent(in) :: x(:), b_dot(:)
        complex(dp), intent(out) :: x_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: column, entry, row

        if (.not. solver%factored) then
            call not_factored(status)
            return
        end if
        if (has_approximate_factors(solver)) then
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                "implicit derivatives require exact factors; SuperLU ILU is approximate")
            return
        end if
        if (.not. csc_is_valid(A_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "sparse_solve_jvp: invalid complex matrix tangent")
            return
        end if
        if (A_dot%nrow /= size(b_dot) .or. A_dot%ncol /= size(x)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "sparse_solve_jvp: incompatible complex tangent dimensions")
            return
        end if
        if (size(x_dot) /= size(b_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "sparse_solve_jvp: incompatible complex output dimension")
            return
        end if
        x_dot = b_dot
        do column = 1, A_dot%ncol
            do entry = A_dot%col_ptr(column), A_dot%col_ptr(column+1)-1
                row = A_dot%row_idx(entry)
                x_dot(row) = lu_z_update(x_dot(row), A_dot%val(entry), x(column))
            end do
        end do
        call sparse_solve_complex_inplace(solver,x_dot,status)
    end subroutine sparse_solve_jvp_complex

    subroutine sparse_solve_vjp_complex( &
            adjoint_solver, A, x, x_bar, b_bar, A_values_bar, status, use_primal_factors)
        ! Exact complex implicit adjoint under the real inner product:
        ! A^H b_bar = x_bar, (A_ij)_bar = -b_bar_i*conjg(x_j).
        ! adjoint_solver must retain a factorization of A^H.
        type(sparse_solver_t), intent(inout) :: adjoint_solver
        type(csc_z_t), intent(in) :: A
        complex(dp), intent(in) :: x(:), x_bar(:)
        complex(dp), intent(out) :: b_bar(:), A_values_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        logical, optional, intent(in) :: use_primal_factors
        logical :: primal
        integer :: column, entry

        if (.not. adjoint_solver%factored) then
            call not_factored(status)
            return
        end if
        if (has_approximate_factors(adjoint_solver)) then
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                "implicit derivatives require exact factors; SuperLU ILU is approximate")
            return
        end if
        if (.not. csc_is_valid(A)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "sparse_solve_vjp: invalid complex primal matrix")
            return
        end if
        if (A%ncol /= size(x) .or. A%ncol /= size(x_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "sparse_solve_vjp: incompatible complex primal dimensions")
            return
        end if
        if (A%nrow /= size(b_bar) .or. A%nnz /= size(A_values_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "sparse_solve_vjp: incompatible complex cotangent dimensions")
            return
        end if
        primal = .false.
        if (present(use_primal_factors)) primal = use_primal_factors
        if (primal) then
            call sparse_solve_complex_op(adjoint_solver, x_bar, b_bar, FORTSPARSE_OP_ADJOINT, status)
        else
        call sparse_solve_complex(adjoint_solver, x_bar, b_bar, status)
        end if
        if (.not. status_ok(status)) return
        do column = 1, A%ncol
            do entry = A%col_ptr(column), A%col_ptr(column + 1) - 1
                A_values_bar(entry) = &
                    -b_bar(A%row_idx(entry))*conjg(x(column))
            end do
        end do
    end subroutine sparse_solve_vjp_complex

    ! Solve A X = B for all columns of B at once, reusing the factorization.
    ! A backend with a native multi-RHS solve (MUMPS) does it in one call;
    ! others solve column by column. Results equal repeated single solves to
    ! rounding.
    subroutine sparse_solve_block_real(solver, B, X, status)
        type(sparse_solver_t),     intent(inout) :: solver
        real(dp), contiguous,      intent(in)    :: B(:, :)
        real(dp), contiguous,      intent(out)   :: X(:, :)
        type(fortsparse_status_t), intent(out)   :: status

        if (.not. solver%factored) then
            call not_factored(status)
            return
        end if
        if (size(B, 1) /= size(X, 1) .or. size(B, 2) /= size(X, 2)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "sparse_solve_block: B and X shapes differ")
            return
        end if
        call solver%backend%solve_real_block(B, X, status)
    end subroutine sparse_solve_block_real

    ! Complex block solve with the same contract as the real one.
    subroutine sparse_solve_block_complex(solver, B, X, status)
        type(sparse_solver_t),     intent(inout) :: solver
        complex(dp), contiguous,   intent(in)    :: B(:, :)
        complex(dp), contiguous,   intent(out)   :: X(:, :)
        type(fortsparse_status_t), intent(out)   :: status

        if (.not. solver%factored) then
            call not_factored(status)
            return
        end if
        if (size(B, 1) /= size(X, 1) .or. size(B, 2) /= size(X, 2)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "sparse_solve_block: B and X shapes differ")
            return
        end if
        call solver%backend%solve_complex_block(B, X, status)
    end subroutine sparse_solve_block_complex

    ! Last factorization attempt, including a rejected UMFPACK memory estimate.
    ! Unknown quantities remain -1; availability depends on the backend.
    subroutine sparse_factor_stats(solver, stats, status)
        type(sparse_solver_t),       intent(inout) :: solver
        type(sparse_factor_stats_t), intent(out)   :: stats
        type(fortsparse_status_t),   intent(out)   :: status

        if (.not. allocated(solver%backend)) then
            call not_factored(status)
            return
        end if
        call solver%backend%factor_stats(stats, status)
    end subroutine sparse_factor_stats

    ! Sylvester inertia (negative, zero, positive pivot eigenvalue counts) of
    ! the last native symmetric/Hermitian factorization, including one that
    ! returned FORTSPARSE_SINGULAR because of exact zero pivots. reliability is
    ! a FORTSPARSE_INERTIA_* code (docs/native-cpu.md); status is OK whenever
    ! counts are returned and BACKEND_UNAVAILABLE for LU, perturbed, complex
    ! symmetric and external factors, which do not provide inertia.
    subroutine sparse_inertia(solver, n_negative, n_zero, n_positive, &
                              reliability, status)
        type(sparse_solver_t),     intent(inout) :: solver
        integer(i8),               intent(out)   :: n_negative, n_zero, n_positive
        integer,                   intent(out)   :: reliability
        type(fortsparse_status_t), intent(out)   :: status
        type(sparse_factor_stats_t) :: stats

        n_negative = -1_i8
        n_zero = -1_i8
        n_positive = -1_i8
        reliability = FORTSPARSE_INERTIA_UNAVAILABLE
        call sparse_factor_stats(solver, stats, status)
        if (.not. status_ok(status)) return
        if (stats%inertia_status == FORTSPARSE_INERTIA_UNAVAILABLE) then
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                'inertia unavailable: no native LDL^T/LDL^H/Cholesky factor')
            return
        end if
        n_negative = stats%inertia_negative
        n_zero = stats%inertia_zero
        n_positive = stats%inertia_positive
        reliability = stats%inertia_status
    end subroutine sparse_inertia

    ! In-place real solve: b is the RHS on entry, the solution on return.
    subroutine sparse_solve_real_inplace(solver, b, status)
        type(sparse_solver_t),     intent(inout) :: solver
        real(dp),                  intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out)   :: status

        if (.not. solver%factored) then
            call not_factored(status)
            return
        end if
        call solver%backend%solve_real_inplace(b, status)
    end subroutine sparse_solve_real_inplace

    ! In-place complex solve: b is the RHS on entry, the solution on return.
    subroutine sparse_solve_complex_inplace(solver, b, status)
        type(sparse_solver_t),     intent(inout) :: solver
        complex(dp),               intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out)   :: status

        if (.not. solver%factored) then
            call not_factored(status)
            return
        end if
        call solver%backend%solve_complex_inplace(b, status)
    end subroutine sparse_solve_complex_inplace

    ! Release the current factorization, keeping the backend ready to factor
    ! again. For the out-of-process backend this also releases the helper
    ! process and its matrix-sized shared mapping, so a freed solver holds no
    ! memory; the next factor spawns a fresh helper. Vectors from sparse_vector
    ! do not survive a free.
    subroutine sparse_free(solver)
        type(sparse_solver_t), intent(inout) :: solver

        if (allocated(solver%backend)) call solver%backend%free()
        solver%factored = .false.
    end subroutine sparse_free

    ! Allocate a length-n real solve vector owned by the solver. Used as the RHS
    ! and/or solution of sparse_solve, it avoids copying that vector across the
    ! backend boundary: for the out-of-process backend it is a shared-memory
    ! slot the helper reads and writes directly. Valid after factorization (it
    ! is sized to the factorization) and until sparse_free or teardown, so the
    ! caller never frees it. Returns null if no backend is active or the pool
    ! is exhausted.
    function sparse_vector(solver, n) result(p)
        type(sparse_solver_t), intent(inout) :: solver
        integer,               intent(in)    :: n
        real(dp), pointer                    :: p(:)

        if (allocated(solver%backend)) then
            p => solver%backend%vector(n)
        else
            p => null()
        end if
    end function sparse_vector

    ! Tear the backend down, releasing its factorization and any helper process.
    ! A solver releases everything automatically when it goes out of scope (its
    ! backend component is finalized), so most code never needs this; it is for
    ! freeing a long-lived solver's resources before it goes out of scope.
    subroutine sparse_destroy(solver)
        type(sparse_solver_t), intent(inout) :: solver

        if (allocated(solver%backend)) deallocate (solver%backend)
        solver%factored = .false.
    end subroutine sparse_destroy

    ! Deep-copy the context src into dst. Native CPU factors, analysis and
    ! refinement state are duplicated into storage owned by dst, so either
    ! context may be solved, refactored or destroyed independently. Vectors
    ! from sparse_vector are not copied. dst keeps its own
    ! options%host_memory_limit_bytes, and the copied payload is charged
    ! against it (RESOURCE_LIMIT when it does not fit). External backend
    ! factors and native CUDA device factors are not copied: dst then has
    ! src's selection and options but is unfactored, and status is
    ! FORTSPARSE_NOT_FACTORED. On every failure dst is a valid unfactored
    ! context. dst's previous factors are released after the copy is built;
    ! call sparse_destroy(dst) first to lower the transient peak.
    subroutine sparse_copy(dst, src, status)
        type(sparse_solver_t), target, intent(inout) :: dst
        type(sparse_solver_t), target, intent(in) :: src
        type(fortsparse_status_t), intent(out) :: status

        call copy_context(dst, src, status)
    end subroutine sparse_copy

    ! Defined assignment; status is dropped. Use sparse_copy when it matters.
    impure elemental subroutine assign_solver(dst, src)
        class(sparse_solver_t), target, intent(inout) :: dst
        type(sparse_solver_t), target, intent(in) :: src
        type(fortsparse_status_t) :: status

        call copy_context(dst, src, status)
    end subroutine assign_solver

    subroutine copy_context(dst, src, status)
        class(sparse_solver_t), target, intent(inout) :: dst
        type(sparse_solver_t), target, intent(in) :: src
        type(fortsparse_status_t), intent(out) :: status
        class(sparse_backend_t), allocatable :: fresh
        integer(i8) :: cap

        call status_set(status, FORTSPARSE_OK, '')
        if (c_associated(c_loc(dst%backend_id), c_loc(src%backend_id))) return
        ! gfortran passes a shallow temporary for self-assignment;
        ! its backend is dst's backend, so this is also self-assignment.
        if (allocated(dst%backend) .and. allocated(src%backend)) then
            if (c_associated(c_loc(dst%backend%options), &
                             c_loc(src%backend%options))) return
        end if
        cap = dst%options%host_memory_limit_bytes
        ! Build the copy before releasing dst, so a source that shares
        ! storage with dst is read before that storage is released.
        if (allocated(src%backend)) then
            allocate (fresh, mold=src%backend)
            fresh%options = src%backend%options
            fresh%options%host_memory_limit_bytes = cap
            call fresh%copy_from(src%backend, status)
            ! An unfactored source loses nothing when no factors are copied.
            if (status%code == FORTSPARSE_NOT_FACTORED .and. .not. src%factored) &
                call status_set(status, FORTSPARSE_OK, '')
        end if
        if (allocated(dst%backend)) deallocate (dst%backend)
        if (allocated(fresh)) call move_alloc(fresh, dst%backend)
        dst%backend_id = src%backend_id
        dst%refine = src%refine
        dst%options = src%options
        dst%options%host_memory_limit_bytes = cap
        dst%factored = src%factored .and. status_ok(status)
    end subroutine copy_context

    ! Convenience real driver: factor, solve, free.
    subroutine sparse_solve_once_real(A, b, x, status)
        type(csc_t),               intent(in)  :: A
        real(dp),                  intent(in)  :: b(:)
        real(dp),                  intent(out) :: x(:)
        type(fortsparse_status_t), intent(out) :: status

        type(sparse_solver_t) :: solver

        call sparse_factor_real(solver, A, status)
        if (.not. status_ok(status)) then
            call sparse_destroy(solver)
            return
        end if
        call sparse_solve_real(solver, b, x, status)
        call sparse_destroy(solver)
    end subroutine sparse_solve_once_real

    ! Convenience complex driver: factor, solve, free.
    subroutine sparse_solve_once_complex(A, b, x, status)
        type(csc_z_t),             intent(in)  :: A
        complex(dp),               intent(in)  :: b(:)
        complex(dp),               intent(out) :: x(:)
        type(fortsparse_status_t), intent(out) :: status

        type(sparse_solver_t) :: solver

        call sparse_factor_complex(solver, A, status)
        if (.not. status_ok(status)) then
            call sparse_destroy(solver)
            return
        end if
        call sparse_solve_complex(solver, b, x, status)
        call sparse_destroy(solver)
    end subroutine sparse_solve_once_complex

    ! Backend chosen by FORTSPARSE_BACKEND_AUTO for a matrix with nnz stored
    ! entries: MUMPS from AUTO_MUMPS_MIN_NNZ on, SuperLU below (README,
    ! "Backend selection", has the measured crossover). MUMPS needs a
    ! live runtime (an MPI MUMPS needs the application to own MPI_Init)
    ! and, for a complex matrix, ZMUMPS; otherwise SuperLU. Without SuperLU,
    ! AUTO uses a ready MUMPS throughout, or the native backend if not ready.
    integer function sparse_auto_backend(nnz, is_complex) result(id)
        integer, intent(in) :: nnz
        logical, intent(in) :: is_complex
        logical :: ready

        ready = AUTO_HAS_MUMPS
#ifdef FORTSPARSE_HAVE_MUMPS
        ready = mumps_available()
#endif
        id = FORTSPARSE_BACKEND_SUPERLU
        if (.not. AUTO_HAS_SUPERLU .and. .not. ready) &
            id = FORTSPARSE_BACKEND_NATIVE
        if (ready .and. nnz >= AUTO_MUMPS_MIN_NNZ .and. &
            (AUTO_HAS_ZMUMPS .or. .not. is_complex)) id = FORTSPARSE_BACKEND_MUMPS
        if (ready .and. .not. AUTO_HAS_SUPERLU) &
            id = FORTSPARSE_BACKEND_MUMPS
    end function sparse_auto_backend

    ! Ensure solver%backend is allocated to the concrete type for backend_id,
    ! resolving FORTSPARSE_BACKEND_AUTO from the matrix about to be factored
    ! (sparse_auto_backend). An existing backend of the selected kind is kept,
    ! so its retained session (the resident helper, for the out-of-process
    ! backend; the symbolic analysis, for a MUMPS refactorization) survives
    ! across factorizations; the backend's own factor replaces any prior
    ! factorization. The backend is reallocated only when absent or when the
    ! selected kind changed, and deallocation runs the old backend's
    ! finalizer, releasing its resources. Unknown ids set the
    ! backend-unavailable status.
    subroutine ensure_backend(solver, nnz, is_complex, status)
        type(sparse_solver_t),     intent(inout) :: solver
        integer,                   intent(in)    :: nnz
        logical,                   intent(in)    :: is_complex
        type(fortsparse_status_t), intent(out)   :: status

        integer :: id

        if (solver%options%host_memory_limit_bytes < 0_i8 .or. &
            solver%options%device_memory_limit_bytes < 0_i8 .or. &
            solver%options%memory_policy < FORTSPARSE_MEMORY_AUTO .or. &
            solver%options%memory_policy > FORTSPARSE_MEMORY_LOW) then
            solver%factored = .false.
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'sparse solver: invalid memory limit or policy')
            return
        end if
        id = solver%backend_id
        if (id == FORTSPARSE_BACKEND_AUTO) id = sparse_auto_backend(nnz, &
            is_complex)
        if (solver%options%structural_matching) then
            if ((id /= FORTSPARSE_BACKEND_NATIVE .and. &
                 id /= FORTSPARSE_BACKEND_NATIVE_REFERENCE) .or. &
                solver%options%symmetry /= FORTSPARSE_SYM_GENERAL .or. &
                solver%options%device /= -1) then
                solver%factored = .false.
                call status_set(status,FORTSPARSE_BACKEND_UNAVAILABLE, &
                    'structural matching requires native CPU GENERAL factorization')
                return
            end if
        end if
        if (solver%options%lu_method < FORTSPARSE_LU_AUTO .or. &
            solver%options%lu_method > FORTSPARSE_LU_SUPERNODAL) then
            solver%factored = .false.
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'sparse solver: invalid LU method')
            return
        end if
        if (solver%options%lu_method == FORTSPARSE_LU_SUPERNODAL) then
            if (id /= FORTSPARSE_BACKEND_NATIVE .or. &
                solver%options%symmetry /= FORTSPARSE_SYM_GENERAL .or. &
                solver%options%device /= -1) then
                solver%factored = .false.
                call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                    'supernodal LU requires native CPU GENERAL factorization')
                return
            end if
        end if
        if (id /= FORTSPARSE_BACKEND_NATIVE .and. &
            id /= FORTSPARSE_BACKEND_NATIVE_REFERENCE) then
            if (solver%options%host_memory_limit_bytes > 0_i8 .or. &
                solver%options%device_memory_limit_bytes > 0_i8 .or. &
                solver%options%memory_policy == FORTSPARSE_MEMORY_LOW) then
                solver%factored = .false.
                call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                    'sparse solver: selected external backend has no strict memory contract')
                return
            end if
        end if
        if (allocated(solver%backend)) then
            if (backend_is(solver%backend, id, is_complex, solver%options%symmetry, &
                solver%options%ldlt_supernodal)) then
                solver%backend%options = solver%options
                call status_set(status, FORTSPARSE_OK, "")
                return
            end if
            deallocate (solver%backend)
        end if
        select case (id)
        case (FORTSPARSE_BACKEND_NATIVE, FORTSPARSE_BACKEND_NATIVE_REFERENCE)
            if (.not. native_cholesky_family(is_complex, solver%options%symmetry, &
                solver%options%ldlt_supernodal)) then
                allocate(native_general_backend_t :: solver%backend)
            else
                allocate(native_backend_t :: solver%backend)
            end if
            select type (backend => solver%backend)
            type is (native_backend_t)
                backend%reference_only = id == FORTSPARSE_BACKEND_NATIVE_REFERENCE
            type is (native_general_backend_t)
                backend%reference_only = id == FORTSPARSE_BACKEND_NATIVE_REFERENCE
            end select
#ifdef FORTSPARSE_HAVE_SUPERLU
        case (FORTSPARSE_BACKEND_SUPERLU, FORTSPARSE_BACKEND_SUPERLU_ILU)
            allocate (superlu_backend_t :: solver%backend)
            select type (backend => solver%backend)
            type is (superlu_backend_t)
                backend%incomplete = id == FORTSPARSE_BACKEND_SUPERLU_ILU
            end select
#endif
        case (FORTSPARSE_BACKEND_UMFPACK_IPC)
            allocate (umfpack_ipc_backend_t :: solver%backend)
#ifdef FORTSPARSE_HAVE_MUMPS
        case (FORTSPARSE_BACKEND_MUMPS)
            allocate (mumps_backend_t :: solver%backend)
#endif
        case default
            call unavailable_backend(id, status)
            return
        end select
        solver%backend%options = solver%options
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine ensure_backend

    ! Factor identity is retained in the backend: backend_id may already name
    ! a different backend for the next factorization.
    logical function has_approximate_factors(solver) result(approximate)
        type(sparse_solver_t), intent(in) :: solver
        approximate = .false.
#ifdef FORTSPARSE_HAVE_SUPERLU
        if (.not. allocated(solver%backend)) return
        select type (backend => solver%backend)
        type is (superlu_backend_t)
            approximate = backend%incomplete
        end select
#endif
    end function has_approximate_factors

    ! True when an allocated backend is the concrete kind for backend_id.
    ! The native Cholesky backend owns real SPD and, on request, the real
    ! symmetric indefinite supernodal LDL^T.
    pure logical function native_cholesky_family(is_complex, symmetry, &
            ldlt_supernodal) result(owned)
        logical, intent(in) :: is_complex, ldlt_supernodal
        integer, intent(in) :: symmetry
        owned = .not. is_complex .and. (symmetry == FORTSPARSE_SYM_SPD .or. &
            (symmetry == FORTSPARSE_SYM_INDEFINITE .and. ldlt_supernodal))
    end function native_cholesky_family

    logical function backend_is(backend, backend_id, is_complex, symmetry, &
            ldlt_supernodal)
        class(sparse_backend_t), intent(in) :: backend
        integer,                 intent(in) :: backend_id, symmetry
        logical, intent(in) :: is_complex, ldlt_supernodal

        backend_is = .false.
        select type (backend)
            type is (native_general_backend_t)
            if (backend%reference_only .neqv. &
                (backend_id == FORTSPARSE_BACKEND_NATIVE_REFERENCE)) return
            backend_is = (backend_id == FORTSPARSE_BACKEND_NATIVE .or. &
                backend_id == FORTSPARSE_BACKEND_NATIVE_REFERENCE) .and. &
                .not. native_cholesky_family(is_complex, symmetry, ldlt_supernodal)
            type is (native_backend_t)
            if (.not. native_cholesky_family(is_complex, symmetry, ldlt_supernodal)) return
            if (backend%reference_only) then
                backend_is = backend_id == FORTSPARSE_BACKEND_NATIVE_REFERENCE
            else
                backend_is = backend_id == FORTSPARSE_BACKEND_NATIVE
            end if
#ifdef FORTSPARSE_HAVE_SUPERLU
            type is (superlu_backend_t)
            if (backend%incomplete) then
                backend_is = backend_id == FORTSPARSE_BACKEND_SUPERLU_ILU
            else
                backend_is = backend_id == FORTSPARSE_BACKEND_SUPERLU
            end if
#endif
            type is (umfpack_ipc_backend_t)
            backend_is = (backend_id == FORTSPARSE_BACKEND_UMFPACK_IPC)
#ifdef FORTSPARSE_HAVE_MUMPS
            type is (mumps_backend_t)
            backend_is = (backend_id == FORTSPARSE_BACKEND_MUMPS)
#endif
        end select
    end function backend_is

    ! Set status for a solve attempted before factorization.
    subroutine not_factored(status)
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_NOT_FACTORED, &
            "sparse_solve: matrix has not been factored")
    end subroutine not_factored

    ! Set status for a backend tag that is not available in this build.
    subroutine unavailable_backend(backend_id, status)
        integer,                   intent(in)  :: backend_id
        type(fortsparse_status_t), intent(out) :: status
        character(32)                          :: tag

        write (tag, '(i0)') backend_id
        call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
            "sparse solver: backend "//trim(tag)//" is not available")
    end subroutine unavailable_backend

end module fortsparse_solver
