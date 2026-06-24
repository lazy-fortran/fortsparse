program test_fortsparse_smoke
    ! Smoke test: the umbrella module exposes a non-empty version string.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: fortsparse_version_string
    implicit none

    integer :: nfail
    nfail = 0

    call check_true("version_nonempty", len_trim(fortsparse_version_string) > 0, &
        nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    subroutine check_true(label, cond, nfail)
        character(*), intent(in)    :: label
        logical,      intent(in)    :: cond
        integer,      intent(inout) :: nfail
        if (.not. cond) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a)") "FAIL [", label, "] expected .true."
        end if
    end subroutine check_true

end program test_fortsparse_smoke
