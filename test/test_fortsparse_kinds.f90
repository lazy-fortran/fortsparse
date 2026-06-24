program test_fortsparse_kinds
    ! Verify the exported kind parameters resolve to the expected
    ! iso_fortran_env kinds and that arithmetic uses the declared kind.
    use, intrinsic :: iso_fortran_env, only: real64, real32, int32, int64, &
        error_unit
    use fortsparse, only: dp, sp, i4, i8
    implicit none

    integer :: nfail
    nfail = 0

    call check_int("dp_real64", dp, real64, nfail)
    call check_int("sp_real32", sp, real32, nfail)
    call check_int("i4_int32", i4, int32, nfail)
    call check_int("i8_int64", i8, int64, nfail)

    ! Double precision must carry at least 15 significant decimal digits.
    call check_true("dp_precision", precision(1.0_dp) >= 15, nfail)
    ! 64-bit integers must hold values past the 32-bit range.
    call check_true("i8_range", range(1_i8) > range(1_i4), nfail)

    block
        real(dp)    :: x
        real(sp)    :: y
        integer(i4) :: n4
        integer(i8) :: n8
        x  = 1.0_dp
        y  = 1.0_sp
        n4 = 1_i4
        n8 = 1_i8
        call check_true("kind_dp", kind(x) == dp, nfail)
        call check_true("kind_sp", kind(y) == sp, nfail)
        call check_true("kind_i4", kind(n4) == i4, nfail)
        call check_true("kind_i8", kind(n8) == i8, nfail)
    end block

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    subroutine check_int(label, got, expected, nfail)
        character(*), intent(in)    :: label
        integer,      intent(in)    :: got, expected
        integer,      intent(inout) :: nfail
        if (got /= expected) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,i0,a,i0)") "FAIL [", label, "] got=", &
                got, " expected=", expected
        end if
    end subroutine check_int

    subroutine check_true(label, cond, nfail)
        character(*), intent(in)    :: label
        logical,      intent(in)    :: cond
        integer,      intent(inout) :: nfail
        if (.not. cond) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a)") "FAIL [", label, "] expected .true."
        end if
    end subroutine check_true

end program test_fortsparse_kinds
