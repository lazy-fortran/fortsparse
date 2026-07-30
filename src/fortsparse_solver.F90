module fortsparse_solver
    ! Public, backend-pluggable sparse direct solver API.
    !
    ! A client factors a matrix once, solves any number of right-hand sides
    ! reusing that factorization, then frees it. The backend is selected by an
    ! integer tag on the solver. The default in-process backend is SuperLU
    ! (BSD); UMFPACK runs out-of-process through a separate GPL helper and is
    ! reached by selecting FORTSPARSE_BACKEND_UMFPACK_IPC. Selecting an
    ! unbuilt or unknown backend yields a documented
    ! FORTSPARSE_BACKEND_UNAVAILABLE status rather than a crash.
    use fortsparse_kinds, only: dp
    use fortsparse_status, only: fortsparse_status_t, status_set, status_ok, &
        FORTSPARSE_OK, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_INVALID_MATRIX
    use fortsparse_csc, only: csc_t, csc_z_t, csc_is_valid, csc_matvec
    use fortsparse_backend, only: sparse_backend_t
#ifdef FORTSPARSE_HAVE_SUPERLU
    use fortsparse_superlu, only: superlu_backend_t
#endif
    use fortsparse_umfpack_ipc, only: umfpack_ipc_backend_t
    implicit none
    private

    ! Backend identifiers. SuperLU is the default in-process backend; the
    ! UMFPACK_IPC id selects the out-of-process GPL helper. KLU and PARDISO are
    ! reserved for later releases.
    integer, parameter, public :: FORTSPARSE_BACKEND_SUPERLU = 1
    integer, parameter, public :: FORTSPARSE_BACKEND_UMFPACK_IPC = 2

    ! Solver handle. Holds backend selection, refinement toggle, the factored
    ! flag, and the polymorphic backend that owns the retained factorization.
    type, public :: sparse_solver_t
        integer :: backend_id = FORTSPARSE_BACKEND_SUPERLU
        logical :: refine = .false.
        logical :: factored = .false.
        class(sparse_backend_t), allocatable :: backend
    end type sparse_solver_t

    public :: sparse_factor
    public :: sparse_solve
    public :: sparse_free
    public :: sparse_destroy
    public :: sparse_solve_once
    public :: sparse_solve_jvp
    public :: sparse_solve_vjp
    public :: sparse_vector

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

    ! Factor a real matrix. Ensures the concrete backend, then dispatches.
    subroutine sparse_factor_real(solver, A, status)
        type(sparse_solver_t),     intent(inout) :: solver
        type(csc_t),               intent(in)    :: A
        type(fortsparse_status_t), intent(out)   :: status

        call ensure_backend(solver, status)
        if (.not. status_ok(status)) return
        call solver%backend%factor_real(A, solver%refine, status)
        solver%factored = status_ok(status)
    end subroutine sparse_factor_real

    ! Factor a complex matrix. Ensures the concrete backend, then dispatches.
    subroutine sparse_factor_complex(solver, A, status)
        type(sparse_solver_t),     intent(inout) :: solver
        type(csc_z_t),             intent(in)    :: A
        type(fortsparse_status_t), intent(out)   :: status

        call ensure_backend(solver, status)
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

        call ensure_backend(solver, status)
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

        call ensure_backend(solver, status)
        if (.not. status_ok(status)) return
        call solver%backend%factor_complex_raw(nrow, ncol, nz, col_ptr, &
            row_idx, val, solver%refine, status)
        solver%factored = status_ok(status)
    end subroutine sparse_factor_complex_raw

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

        real(dp), allocatable, target :: A_dot_x(:), tangent_rhs(:)

        if (.not. solver%factored) then
            call not_factored(status)
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
        A_dot_x = csc_matvec(A_dot, x)
        tangent_rhs = b_dot - A_dot_x
        call sparse_solve_real(solver, tangent_rhs, x_dot, status)
    end subroutine sparse_solve_jvp_real

    subroutine sparse_solve_vjp_real( &
            transpose_solver, A, x, x_bar, b_bar, A_values_bar, status)
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

        integer :: column, entry

        if (.not. transpose_solver%factored) then
            call not_factored(status)
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
        call sparse_solve_real(transpose_solver, x_bar, b_bar, status)
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

        complex(dp), allocatable :: A_dot_x(:), tangent_rhs(:)

        if (.not. solver%factored) then
            call not_factored(status)
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
        A_dot_x = csc_matvec(A_dot, x)
        tangent_rhs = b_dot - A_dot_x
        call sparse_solve_complex(solver, tangent_rhs, x_dot, status)
    end subroutine sparse_solve_jvp_complex

    subroutine sparse_solve_vjp_complex( &
            adjoint_solver, A, x, x_bar, b_bar, A_values_bar, status)
        ! Exact complex implicit adjoint under the real inner product:
        ! A^H b_bar = x_bar, (A_ij)_bar = -b_bar_i*conjg(x_j).
        ! adjoint_solver must retain a factorization of A^H.
        type(sparse_solver_t), intent(inout) :: adjoint_solver
        type(csc_z_t), intent(in) :: A
        complex(dp), intent(in) :: x(:), x_bar(:)
        complex(dp), intent(out) :: b_bar(:), A_values_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: column, entry

        if (.not. adjoint_solver%factored) then
            call not_factored(status)
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
        call sparse_solve_complex(adjoint_solver, x_bar, b_bar, status)
        if (.not. status_ok(status)) return
        do column = 1, A%ncol
            do entry = A%col_ptr(column), A%col_ptr(column + 1) - 1
                A_values_bar(entry) = &
                    -b_bar(A%row_idx(entry))*conjg(x(column))
            end do
        end do
    end subroutine sparse_solve_vjp_complex

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

    ! Ensure solver%backend is allocated to the concrete type for backend_id.
    ! An existing backend of the selected kind is kept, so its retained session
    ! (the resident helper, for the out-of-process backend) survives across
    ! factorizations; the backend's own factor replaces any prior factorization.
    ! The backend is reallocated only when absent or when backend_id changed, and
    ! deallocation runs the old backend's finalizer, releasing its resources.
    ! Unknown ids set the backend-unavailable status.
    subroutine ensure_backend(solver, status)
        type(sparse_solver_t),     intent(inout) :: solver
        type(fortsparse_status_t), intent(out)   :: status

        if (allocated(solver%backend)) then
            if (backend_is(solver%backend, solver%backend_id)) then
                call status_set(status, FORTSPARSE_OK, "")
                return
            end if
            deallocate (solver%backend)
        end if
        select case (solver%backend_id)
#ifdef FORTSPARSE_HAVE_SUPERLU
        case (FORTSPARSE_BACKEND_SUPERLU)
            allocate (superlu_backend_t :: solver%backend)
#endif
        case (FORTSPARSE_BACKEND_UMFPACK_IPC)
            allocate (umfpack_ipc_backend_t :: solver%backend)
        case default
            call unavailable_backend(solver%backend_id, status)
            return
        end select
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine ensure_backend

    ! True when an allocated backend is the concrete kind for backend_id.
    logical function backend_is(backend, backend_id)
        class(sparse_backend_t), intent(in) :: backend
        integer,                 intent(in) :: backend_id

        backend_is = .false.
        select type (backend)
#ifdef FORTSPARSE_HAVE_SUPERLU
            type is (superlu_backend_t)
            backend_is = (backend_id == FORTSPARSE_BACKEND_SUPERLU)
#endif
            type is (umfpack_ipc_backend_t)
            backend_is = (backend_id == FORTSPARSE_BACKEND_UMFPACK_IPC)
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
