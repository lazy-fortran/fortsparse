module fortsparse_compiler
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated
    implicit none
    private
    public :: require_safe_self_assignment

    type, abstract :: probe_backend_t
        integer :: identity = 0
    end type probe_backend_t

    type, extends(probe_backend_t) :: probe_owner_t
        integer, pointer :: finalizations => null()
    contains
        final :: finalize_probe
    end type probe_owner_t

    type :: probe_solver_t
        integer :: identity = 0
        class(probe_backend_t), allocatable :: backend
    contains
        procedure :: assign_probe
        generic :: assignment(=) => assign_probe
    end type probe_solver_t

contains

    subroutine require_safe_self_assignment()
        type(probe_solver_t) :: probe
        integer, target :: finalizations

        finalizations = 0
        allocate (probe_owner_t :: probe%backend)
        select type (owner => probe%backend)
        type is (probe_owner_t)
            owner%finalizations => finalizations
        end select
        probe = probe
        ! An affected compiler has already released the aliased backend.
        if (finalizations /= 0) error stop &
            'compiler finalizes a live backend during defined self-assignment'
    end subroutine require_safe_self_assignment

    impure elemental subroutine assign_probe(dst, src)
        class(probe_solver_t), target, intent(inout) :: dst
        type(probe_solver_t), target, intent(in) :: src

        if (c_associated(c_loc(dst%identity), c_loc(src%identity))) return
        if (c_associated(c_loc(dst%backend%identity), &
                         c_loc(src%backend%identity))) return
        error stop 'self-assignment compiler probe received different objects'
    end subroutine assign_probe

    subroutine finalize_probe(self)
        type(probe_owner_t), intent(inout) :: self

        if (associated(self%finalizations)) &
            self%finalizations = self%finalizations + 1
    end subroutine finalize_probe

end module fortsparse_compiler
