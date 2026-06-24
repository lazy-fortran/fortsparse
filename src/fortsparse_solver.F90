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
        FORTSPARSE_OK, FORTSPARSE_NOT_FACTORED, FORTSPARSE_BACKEND_UNAVAILABLE
    use fortsparse_csc, only: csc_t, csc_z_t
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

    ! Factor a real or complex matrix into the solver handle.
    interface sparse_factor
        module procedure sparse_factor_real
        module procedure sparse_factor_complex
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

    ! Solve A x = b for a real RHS, reusing the stored factorization.
    subroutine sparse_solve_real(solver, b, x, status)
        type(sparse_solver_t),     intent(inout) :: solver
        real(dp),                  intent(in)    :: b(:)
        real(dp),                  intent(out)   :: x(:)
        type(fortsparse_status_t), intent(out)   :: status

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
    ! again. For the out-of-process backend the helper stays resident, so the
    ! next factor reuses it. The backend and any resource it owns (the helper
    ! process) are released when the solver is finalized, or at once by
    ! sparse_destroy.
    subroutine sparse_free(solver)
        type(sparse_solver_t), intent(inout) :: solver

        if (allocated(solver%backend)) call solver%backend%free()
        solver%factored = .false.
    end subroutine sparse_free

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
