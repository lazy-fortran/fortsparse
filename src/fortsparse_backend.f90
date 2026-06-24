module fortsparse_backend
    ! Abstract sparse-solver backend. A concrete backend retains a
    ! factorization across solves and maps its native error codes onto the
    ! fortsparse status codes. The solver holds one of these polymorphically and
    ! dispatches factor/solve/free through the deferred bindings, so adding a
    ! backend never touches client code or the solver dispatch logic.
    use fortsparse_kinds, only: dp
    use fortsparse_csc, only: csc_t, csc_z_t
    use fortsparse_status, only: fortsparse_status_t
    implicit none
    private

    public :: sparse_backend_t

    type, abstract :: sparse_backend_t
    contains
        procedure(factor_real_i), deferred :: factor_real
        procedure(factor_complex_i), deferred :: factor_complex
        procedure(solve_real_i), deferred :: solve_real
        procedure(solve_complex_i), deferred :: solve_complex
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
            class(sparse_backend_t),   intent(inout) :: self
            real(dp),                  intent(in)    :: b(:)
            real(dp),                  intent(out)   :: x(:)
            type(fortsparse_status_t), intent(out)   :: status
        end subroutine solve_real_i

        subroutine solve_complex_i(self, b, x, status)
            import :: sparse_backend_t, dp, fortsparse_status_t
            class(sparse_backend_t),   intent(inout) :: self
            complex(dp),               intent(in)    :: b(:)
            complex(dp),               intent(out)   :: x(:)
            type(fortsparse_status_t), intent(out)   :: status
        end subroutine solve_complex_i

        subroutine free_i(self)
            import :: sparse_backend_t
            class(sparse_backend_t), intent(inout) :: self
        end subroutine free_i

    end interface

end module fortsparse_backend
