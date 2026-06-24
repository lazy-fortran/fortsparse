program test_fortsparse_status
    ! Behavioral tests for the status carrier: default OK, status_set updates
    ! code and message, non-OK codes flip status_ok to false, and the named
    ! codes are distinct and non-negative.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: fortsparse_status_t, status_ok, status_set, &
        FORTSPARSE_OK, FORTSPARSE_SINGULAR, FORTSPARSE_INVALID_MATRIX, &
        FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_INTERNAL_ERROR
    implicit none

    type(fortsparse_status_t) :: s
    integer                   :: nfail
    nfail = 0

    ! Default-constructed status is OK with an empty message.
    call check_true("default_ok", status_ok(s), nfail)
    call check_int("default_code", s%code, FORTSPARSE_OK, nfail)
    call check_true("default_msg_empty", len_trim(s%msg) == 0, nfail)

    ! status_set writes both code and message; a non-OK code is not ok.
    call status_set(s, FORTSPARSE_SINGULAR, "matrix is singular")
    call check_int("singular_code", s%code, FORTSPARSE_SINGULAR, nfail)
    call check_true("singular_msg", trim(s%msg) == "matrix is singular", nfail)
    call check_true("singular_not_ok",.not. status_ok(s), nfail)

    call status_set(s, FORTSPARSE_INVALID_MATRIX, "bad")
    call check_true("invalid_not_ok",.not. status_ok(s), nfail)
    call status_set(s, FORTSPARSE_BACKEND_UNAVAILABLE, "bad")
    call check_true("backend_not_ok",.not. status_ok(s), nfail)
    call status_set(s, FORTSPARSE_NOT_FACTORED, "bad")
    call check_true("not_factored_not_ok",.not. status_ok(s), nfail)
    call status_set(s, FORTSPARSE_INTERNAL_ERROR, "bad")
    call check_true("internal_not_ok",.not. status_ok(s), nfail)

    ! Resetting to OK clears the error state.
    call status_set(s, FORTSPARSE_OK, "")
    call check_true("reset_ok", status_ok(s), nfail)

    ! Named codes: OK is non-negative and every failure code is distinct.
    call check_true("ok_nonneg", FORTSPARSE_OK >= 0, nfail)
    call check_true("ok_is_zero", FORTSPARSE_OK == 0, nfail)
    call check_true("codes_distinct", all_distinct( &
        [FORTSPARSE_OK, FORTSPARSE_SINGULAR, FORTSPARSE_INVALID_MATRIX, &
        FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_INTERNAL_ERROR]), nfail)

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

    pure logical function all_distinct(v) result(ok)
        integer, intent(in) :: v(:)
        integer             :: i, j
        ok = .true.
        do i = 1, size(v)
            do j = i + 1, size(v)
                if (v(i) == v(j)) ok = .false.
            end do
        end do
    end function all_distinct

end program test_fortsparse_status
