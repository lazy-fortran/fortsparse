program test_fortsparse_csc
    ! CSC construction from triplets: duplicate (row, col) entries are summed,
    ! structural validity holds for a well-formed matrix and fails for a
    ! deliberately corrupted one.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: dp, csc_t, csc_from_triplet, csc_is_valid, &
        csc_matmul, csc_matvec, csc_transpose, fortsparse_status_t, status_ok, &
        csc_triplet_plan_t, csc_triplet_plan_build, csc_triplet_plan_apply
    implicit none

    integer :: nfail
    nfail = 0

    call test_duplicate_sum(nfail)
    call test_triplet_plan(nfail)
    call test_triplet_plan_parallel(nfail)
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

    subroutine test_triplet_plan(nfail)
        integer, intent(inout) :: nfail
        integer, parameter :: rows(7) = [3, 2, 1, 3, 3, 2, 2]
        integer, parameter :: cols(7) = [4, 2, 2, 4, 4, 2, 4]
        real(dp), parameter :: vals(7) = [1.0e16_dp, 3.0_dp, -2.0_dp, &
            -1.0e16_dp, 1.0_dp, -3.0_dp, 5.0_dp]
        type(csc_triplet_plan_t) :: plan, other
        type(csc_t) :: a
        type(fortsparse_status_t) :: status

        call csc_triplet_plan_apply(plan, vals, a, status)
        call check_true("unbuilt_plan", .not. status_ok(status), nfail)
        call csc_triplet_plan_build(3, 5, rows, cols, plan, status)
        call check_true("plan_build", status_ok(status), nfail)
        call csc_triplet_plan_build(3, 5, rows, 6 - cols, other, status)
        call check_true("other_plan_build", status_ok(status), nfail)
        call csc_triplet_plan_apply(plan, vals, a, status)
        if (.not. csc_is_valid(a)) then
            call check_true("plan_valid_matrix", .false., nfail)
            return
        end if
        call check_true("plan_apply", status_ok(status), nfail)
        call check_true("plan_size", a%nnz == 4, nfail)
        if (a%nnz /= 4) return
        call check_true("plan_columns", all(a%col_ptr == [1, 1, 3, 3, 5, 5]), nfail)
        call check_true("plan_rows", all(a%row_idx == [1, 2, 2, 3]), nfail)
        call check_true("plan_duplicate_order", &
            all(a%val == [-2.0_dp, 0.0_dp, 5.0_dp, 1.0_dp]), nfail)

        call csc_triplet_plan_apply(other, vals, a, status)
        call check_true("other_plan_apply", status_ok(status), nfail)
        call check_true("same_size_independent_plan", &
            all(a%row_idx == [2, 3, 1, 2]) .and. &
            all(a%val == [5.0_dp, 1.0_dp, -2.0_dp, 0.0_dp]), nfail)
        call csc_triplet_plan_apply(plan, 2.0_dp*vals, a, status)
        call check_true("plan_reapply", status_ok(status), nfail)
        call check_true("plan_new_values", &
            all(a%val == [-4.0_dp, 0.0_dp, 10.0_dp, 2.0_dp]), nfail)
        call csc_triplet_plan_apply(plan, vals(:6), a, status)
        call check_true("plan_wrong_value_count", .not. status_ok(status), nfail)

        call csc_triplet_plan_build(2, 5, rows, cols, plan, status)
        call check_true("plan_bad_coordinates", .not. status_ok(status), nfail)
        call csc_triplet_plan_apply(plan, vals, a, status)
        call check_true("failed_rebuild_clears_plan", .not. status_ok(status), nfail)
        call csc_triplet_plan_build(3, 5, rows(:0), cols(:0), plan, status)
        call check_true("empty_plan_build", status_ok(status), nfail)
        call csc_triplet_plan_apply(plan, vals(:0), a, status)
        call check_true("empty_plan_apply", status_ok(status), nfail)
        call check_true("empty_plan_matrix", csc_is_valid(a) .and. a%nnz == 0, nfail)
    end subroutine test_triplet_plan

    subroutine test_triplet_plan_parallel(nfail)
        integer, intent(inout) :: nfail
        integer, parameter :: n = 20021
        integer :: cols(2*n), i
        real(dp) :: vals(2*n)
        type(csc_triplet_plan_t) :: plan
        type(csc_t) :: a
        type(fortsparse_status_t) :: status

        cols(:n) = [(i, i=n, 1, -1)]
        cols(n + 1:) = [(i, i=1, n)]
        vals(:n) = real(cols(:n), dp)
        vals(n + 1:) = -0.5_dp*real(cols(n + 1:), dp)
        call csc_triplet_plan_build(1, n, spread(1, 1, 2*n), cols, plan, status)
        call check_true("parallel_plan_build", status_ok(status), nfail)
        call csc_triplet_plan_apply(plan, vals, a, status)
        call check_true("parallel_plan_apply", status_ok(status), nfail)
        call check_true("parallel_plan_size", a%nnz == n, nfail)
        if (a%nnz /= n) return
        call check_true("parallel_plan_values", all(a%val == -vals(n + 1:)), nfail)
    end subroutine test_triplet_plan_parallel

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
