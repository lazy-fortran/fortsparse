program test_fortsparse_csc
    ! CSC construction from triplets: duplicate (row, col) entries are summed,
    ! structural validity holds for a well-formed matrix and fails for a
    ! deliberately corrupted one.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: dp, csc_t, csc_from_triplet, csc_is_valid, &
        fortsparse_status_t, status_ok
    implicit none

    integer :: nfail
    nfail = 0

    call test_duplicate_sum(nfail)
    call test_valid_and_invalid(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    ! Two triplets at (1,1) must collapse to one summed entry.
    subroutine test_duplicate_sum(nfail)
        integer, intent(inout) :: nfail

        type(csc_t)               :: A
        type(fortsparse_status_t) :: status
        integer                   :: rows(4), cols(4)
        real(dp)                  :: vals(4)

        rows = [1, 1, 2, 2]
        cols = [1, 1, 2, 1]
        vals = [1.0_dp, 2.5_dp, 4.0_dp, 7.0_dp]
        call csc_from_triplet(2, 2, rows, cols, vals, A, status)
        call check_true("dup_status_ok", status_ok(status), nfail)
        call check_true("dup_nnz", A%nnz == 3, nfail)
        call check_true("dup_valid", csc_is_valid(A), nfail)
        ! Column 1 holds (1,1)=3.5 and (2,1)=7.0; row order ascending.
        call check_true("dup_colptr1", A%col_ptr(1) == 1, nfail)
        call check_close("dup_val11", A%val(1), 3.5_dp, nfail)
        call check_close("dup_val21", A%val(2), 7.0_dp, nfail)
    end subroutine test_duplicate_sum

    ! A built matrix is valid; corrupting col_ptr makes it invalid.
    subroutine test_valid_and_invalid(nfail)
        integer, intent(inout) :: nfail

        type(csc_t)               :: A
        type(fortsparse_status_t) :: status
        integer                   :: rows(3), cols(3)
        real(dp)                  :: vals(3)

        rows = [1, 2, 3]
        cols = [1, 2, 3]
        vals = [1.0_dp, 2.0_dp, 3.0_dp]
        call csc_from_triplet(3, 3, rows, cols, vals, A, status)
        call check_true("diag_valid", csc_is_valid(A), nfail)
        A%col_ptr(2) = 99
        call check_true("corrupt_invalid", .not. csc_is_valid(A), nfail)
    end subroutine test_valid_and_invalid

    subroutine check_true(label, cond, nfail)
        character(*), intent(in)    :: label
        logical,      intent(in)    :: cond
        integer,      intent(inout) :: nfail
        if (.not. cond) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a)") "FAIL [", label, "] expected .true."
        end if
    end subroutine check_true

    subroutine check_close(label, got, want, nfail)
        character(*), intent(in)    :: label
        real(dp),     intent(in)    :: got, want
        integer,      intent(inout) :: nfail
        if (abs(got - want) > 1.0e-12_dp) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,es13.6,a,es13.6)") "FAIL [", label, &
                "] got ", got, " want ", want
        end if
    end subroutine check_close

end program test_fortsparse_csc
