program test_fortsparse_csc
    ! CSC construction from triplets: duplicate (row, col) entries are summed,
    ! structural validity holds for a well-formed matrix and fails for a
    ! deliberately corrupted one.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: dp, csc_t, csc_from_triplet, csc_is_valid, &
        csc_matmul, csc_matvec, csc_transpose, fortsparse_status_t, status_ok
    implicit none

    integer :: nfail
    nfail = 0

    call test_duplicate_sum(nfail)
    call test_valid_and_invalid(nfail)
    call test_transpose_pairing(nfail)
    call test_sparse_product(nfail)

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

    subroutine test_transpose_pairing(nfail)
        integer, intent(inout) :: nfail

        type(csc_t) :: A, transpose_A
        type(fortsparse_status_t) :: status
        integer :: rows(5), cols(5)
        real(dp) :: left, right, vals(5), x(3), y(2)

        rows = [1, 2, 1, 2, 1]
        cols = [1, 1, 2, 2, 3]
        vals = [2.0_dp, -1.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
        x = [0.25_dp, -2.0_dp, 1.5_dp]
        y = [-0.75_dp, 2.5_dp]
        call csc_from_triplet(2, 3, rows, cols, vals, A, status)
        call csc_transpose(A, transpose_A, status)
        left = dot_product(csc_matvec(A, x), y)
        right = dot_product(x, csc_matvec(transpose_A, y))
        call check_true("transpose_status", status_ok(status), nfail)
        call check_true("transpose_shape", &
            transpose_A%nrow == 3 .and. transpose_A%ncol == 2, nfail)
        call check_close("transpose_pairing", left, right, nfail)
    end subroutine test_transpose_pairing

    subroutine test_sparse_product(nfail)
        integer, intent(inout) :: nfail

        type(csc_t) :: A, B, C
        type(fortsparse_status_t) :: status
        integer :: a_rows(5), a_cols(5), b_rows(4), b_cols(4)
        real(dp) :: a_vals(5), b_vals(4), composed(2), product(2), x(2)

        a_rows = [1, 2, 1, 2, 1]
        a_cols = [1, 1, 2, 2, 3]
        a_vals = [2.0_dp, -1.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
        b_rows = [1, 2, 3, 2]
        b_cols = [1, 1, 1, 2]
        b_vals = [1.0_dp, -2.0_dp, 0.5_dp, 3.0_dp]
        x = [0.6_dp, -1.25_dp]
        call csc_from_triplet(2, 3, a_rows, a_cols, a_vals, A, status)
        call csc_from_triplet(3, 2, b_rows, b_cols, b_vals, B, status)
        call csc_matmul(A, B, C, status)
        composed = csc_matvec(A, csc_matvec(B, x))
        product = csc_matvec(C, x)
        call check_true("product_status", status_ok(status), nfail)
        call check_true("product_shape", C%nrow == 2 .and. C%ncol == 2, nfail)
        call check_true("product_valid", csc_is_valid(C), nfail)
        call check_close("product_row_1", product(1), composed(1), nfail)
        call check_close("product_row_2", product(2), composed(2), nfail)
    end subroutine test_sparse_product

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
