module fortsparse_status
    implicit none
    private

    ! Named error codes.  Zero always means success; positive values are
    ! domain-specific failures; keep codes stable across releases because
    ! callers may branch on them.
    integer, parameter, public :: FORTSPARSE_OK                  = 0
    integer, parameter, public :: FORTSPARSE_SINGULAR            = 1
    integer, parameter, public :: FORTSPARSE_INVALID_MATRIX      = 2
    integer, parameter, public :: FORTSPARSE_BACKEND_UNAVAILABLE = 3
    integer, parameter, public :: FORTSPARSE_NOT_FACTORED        = 4
    integer, parameter, public :: FORTSPARSE_INTERNAL_ERROR      = 5

    ! Maximum length for the human-readable status message.
    integer, parameter :: MSG_LEN = 120

    ! Opaque status carrier.  Pass by value in pure contexts; pass by reference
    ! when the callee must update it (intent(out) or intent(inout)).
    type, public :: fortsparse_status_t
        integer            :: code = FORTSPARSE_OK
        character(MSG_LEN) :: msg  = ""
    end type fortsparse_status_t

    public :: status_ok
    public :: status_set

contains

    ! Returns .true. iff the status indicates no error.
    pure logical function status_ok(s)
        type(fortsparse_status_t), intent(in) :: s
        status_ok = (s%code == FORTSPARSE_OK)
    end function status_ok

    ! Sets the code and message on an existing status object.
    pure subroutine status_set(s, code, msg)
        type(fortsparse_status_t), intent(out) :: s
        integer,                   intent(in)  :: code
        character(*),              intent(in)  :: msg
        s%code = code
        s%msg  = msg
    end subroutine status_set

end module fortsparse_status
