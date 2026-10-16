program test_fortsparse_ilutp
    use fortsparse, only: dp, i8, sparse_ilutp_t, sparse_ilutp_stats_t, &
                          fortsparse_status_t, FORTSPARSE_RESOURCE_LIMIT, &
                          FORTSPARSE_INVALID_MATRIX, &
                          FORTSPARSE_SINGULAR, FORTSPARSE_NOT_FACTORED, &
                          FORTSPARSE_INDEX_OVERFLOW
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    type(sparse_ilutp_t) :: preconditioner
    type(sparse_ilutp_stats_t) :: info
    type(fortsparse_status_t) :: status
    integer, parameter :: ptr(5) = [1, 2, 4, 5, 7]
    integer, parameter :: col(6) = [2, 1, 2, 4, 3, 4]
    real(dp), parameter :: val(6) = [2.0_dp, 2.0_dp, 3.0_dp, &
                                     4.0_dp, 4.0_dp, 5.0_dp]
    real(dp) :: b(4), x(4), expected(4), x2(2), x3(3), bad_values(2)
    integer :: i
    integer(i8) :: required, retained
    integer, allocatable :: large_ptr(:)

    ! Two nonsingular zero-diagonal constraint blocks. Their exact inverses
    ! are [-3/4,1/2;1/2,0] and [-5/16,1/4;1/4,0].
    call preconditioner%setup(4, ptr, col, val, 3, 3, 0.0_dp, .1_dp, status)
    call success()
    info = preconditioner%stats()
    if (info%column_swaps /= 2 .or. info%regularized_pivots /= 0) &
        error stop 'constraint column pivot diagnostics'
    required = info%setup_required_bytes
    retained = info%current_bytes
    if (info%peak_owned_bytes /= required) &
        error stop 'setup payload accounting'
    if (required - retained /= 3_i8*4*(storage_size(0)/8)) &
        error stop 'scratch accounting'
    do i = 1, 5
        b = real(i, dp)*[1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
        expected = [-.75_dp*b(1) + .5_dp*b(2), .5_dp*b(1), &
                    -.3125_dp*b(3) + .25_dp*b(4), .25_dp*b(3)]
        call preconditioner%apply(b, x, status)
        call success()
        if (maxval(abs(x - expected)) > 8*epsilon(1.0_dp)) &
            error stop 'exact block inverse'
        info = preconditioner%stats()
        if (info%current_bytes /= retained) &
            error stop 'apply grew context'
    end do
    b(1) = ieee_value(1.0_dp, ieee_quiet_nan)
    call preconditioner%apply(b, x, status)
    if (status%code /= FORTSPARSE_INVALID_MATRIX) &
        error stop 'nonfinite RHS'
    call preconditioner%apply([1.0_dp], x, status)
    if (status%code /= FORTSPARSE_INVALID_MATRIX) error stop 'RHS shape'

    ! Unordered rows, additive duplicates and explicit zero are accepted;
    ! borrowed arrays remain unchanged, including the caller's column order.
    call preconditioner%setup(2, [1, 4, 6], [2, 1, 2, 2, 1], &
                              [1.0_dp, 0.0_dp, 1.0_dp, 3.0_dp, 2.0_dp], 1, 1, &
                              0.0_dp, .1_dp, status)
    call success()
    call preconditioner%apply([1.0_dp, 2.0_dp], x2, status)
    call success()
    if (maxval(abs(x2 - [.25_dp, .5_dp])) > 8*epsilon(1.0_dp)) &
        error stop 'CSR duplicates'
    info = preconditioner%stats()
    if (info%duplicate_entries /= 1 .or. info%stored_zero_entries /= 1) &
        error stop 'duplicate/zero diagnostics'

    ! Two swaps form a non-involutory 3-cycle, detecting perm/inverse mixups.
    call preconditioner%setup(3, [1, 4, 6, 9], [1, 2, 3, 1, 3, 1, 2, 3], &
                              [.125_dp, .25_dp, 2.0_dp, 1.0_dp, &
                               1.0_dp, 1.0_dp, 3.0_dp, 1.0_dp], &
                              2, 2, 0.0_dp, 1.0_dp, status)
    call success()
    call preconditioner%apply([2.5_dp, 1.75_dp, 1.0_dp], x3, status)
    call success()
    if (maxval(abs(x3 - [.5_dp, -.25_dp, 1.25_dp])) > 1.e-14_dp) &
        error stop 'non-involutory column map'

    ! A previous upper row creates a lower entry between pending columns.
    ! Original dense row formulas give this independent dyadic RHS.
    call preconditioner%setup(4, [1, 4, 7, 10, 13], &
                              [3, 1, 2, 4, 2, 1, 4, 1, 3, 3, 4, 2], &
                              [1.0_dp, .125_dp, 2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, &
                               1.0_dp, 1.0_dp, 4.0_dp, 1.0_dp, 3.0_dp, 1.0_dp], &
                              3, 3, 0.0_dp, .1_dp, status)
    call success()
    call preconditioner%apply([.3125_dp, 1.5_dp, 4.75_dp, &
                               4.25_dp], x, status)
    call success()
    if (maxval(abs(x - [.5_dp, -.25_dp, .75_dp, 1.25_dp])) > 1.e-14_dp) &
        error stop 'sparse lower reach and Schur updates'

    ! Lower multiplier equality is dropped before Schur updates.
    call preconditioner%setup(2, [1, 3, 5], [1, 2, 1, 2], &
                              [1.0_dp, .25_dp, .25_dp, 1.0_dp], 1, 1, &
                              .25_dp, 0.0_dp, status)
    call success()
    call preconditioner%apply([1.0_dp, 1.0_dp], x2, status)
    call success()
    if (maxval(abs(x2 - [.75_dp, 1.0_dp])) > 0.0_dp) &
        error stop 'lower threshold equality'
    info = preconditioner%stats()
    if (info%lower_threshold_drops /= 1) &
        error stop 'lower threshold count'

    ! Upper equality is dropped using the ORIGINAL stored row mean.
    call preconditioner%setup(2, [1, 3, 4], [1, 2, 2], &
                              [3.0_dp, 1.0_dp, 2.0_dp], 0, 1, &
                              .5_dp, 0.0_dp, status)
    call success()
    call preconditioner%apply([3.0_dp, 2.0_dp], x2, status)
    call success()
    if (maxval(abs(x2 - [1.0_dp, 1.0_dp])) > 0.0_dp) &
        error stop 'upper threshold equality'
    info = preconditioner%stats()
    if (info%upper_threshold_drops /= 1) &
        error stop 'upper threshold count'
    ! An explicit zero duplicate changes the mean's count before coalescing.
    call preconditioner%setup(2, [1, 4, 5], [1, 2, 1, 2], &
                              [3.0_dp, 1.0_dp, 0.0_dp, 2.0_dp], 0, 1, &
                              .5_dp, 0.0_dp, status)
    call success()
    call preconditioner%apply([3.0_dp, 2.0_dp], x2, status)
    call success()
    if (maxval(abs(x2 - [2.0_dp/3, 1.0_dp])) > epsilon(1.0_dp)) &
        error stop 'stored zero must count in pre-coalescing mean'

    ! Upper fill ties keep the smaller current column, independently of input order.
    call preconditioner%setup(4, [1, 5, 6, 7, 8], [4, 3, 2, 1, 2, 3, 4], &
                              [1.0_dp, 2.0_dp, 2.0_dp, 8.0_dp, &
                               2.0_dp, 4.0_dp, 8.0_dp], &
                              0, 1, 0.0_dp, 0.0_dp, status)
    call success()
    call preconditioner%apply([8.0_dp, 2.0_dp, 8.0_dp, 8.0_dp], x, status)
    call success()
    if (maxval(abs(x - [.75_dp, 1.0_dp, 2.0_dp, 1.0_dp])) > 0.0_dp) &
        error stop 'upper fill selection'
    info = preconditioner%stats()
    if (info%upper_fill_drops /= 2) error stop 'upper fill count'
    call preconditioner%setup(3, [1, 2, 3, 6], [1, 2, 1, 2, 3], &
                              [1.0_dp, 1.0_dp, .25_dp, .5_dp, 1.0_dp], 1, 0, &
                              0.0_dp, 0.0_dp, status)
    call success()
    call preconditioner%apply([1.0_dp, 2.0_dp, 3.0_dp], x3, status)
    call success()
    if (maxval(abs(x3 - [1.0_dp, 2.0_dp, 2.0_dp])) > 0.0_dp) &
        error stop 'lower fill selection'
    info = preconditioner%stats()
    if (info%lower_fill_drops /= 1) error stop 'lower fill count'

    ! Dropping may produce a zero pivot even for a nonzero input row.
    call preconditioner%setup(2, [1, 2, 3], [2, 2], [1.0_dp, 1.0_dp], &
                              0, 0, 0.0_dp, .1_dp, status)
    call success()
    info = preconditioner%stats()
    if (info%regularized_pivots /= 1) &
        error stop 'regularized pivot diagnostic'
    call preconditioner%apply([1.0_dp, 1.0_dp], x2, status)
    call success()
    if (maxval(abs(x2 - [10000.0_dp, 1.0_dp])) > 1.e-10_dp) &
        error stop 'regularized preconditioner action'

    ! A rebuild frees the old factor BEFORE a rejected allocation request.
    call preconditioner%setup(4, ptr, col, val, 3, 3, 0.0_dp, .1_dp, &
                              status, required - 1)
    if (status%code /= FORTSPARSE_RESOURCE_LIMIT) &
        error stop 'setup cap exhaustion'
    info = preconditioner%stats()
    if (info%current_bytes /= 0 .or. info%ready) &
        error stop 'old factor survived failed rebuild'
    if (info%failed_required_bytes /= required) &
        error stop 'required setup payload'
    call preconditioner%apply(b, x, status)
    if (status%code /= FORTSPARSE_NOT_FACTORED) &
        error stop 'failed setup remained ready'
    call preconditioner%setup(4, ptr, col, val, 3, 3, 0.0_dp, &
                              .1_dp, status, required)
    call success()

    ! Late numerical failure also clears partial factors; subsequent setup recovers.
    call preconditioner%setup(2, [1, 2, 3], [1, 2], [1.0_dp, 0.0_dp], &
                              1, 1, 0.0_dp, .1_dp, status)
    if (status%code /= FORTSPARSE_SINGULAR) error stop 'zero input row'
    info = preconditioner%stats()
    if (info%ready .or. info%current_bytes /= 0) &
        error stop 'partial factor retained'
    call preconditioner%setup(2, [1, 3, 3], [1, 1], [huge(1.0_dp), &
                                                     huge(1.0_dp)], &
                              1, 1, 0.0_dp, .1_dp, status)
    if (status%code /= FORTSPARSE_INVALID_MATRIX) &
        error stop 'duplicate overflow'
    bad_values = [1.0_dp, ieee_value(1.0_dp, ieee_quiet_nan)]
    call preconditioner%setup(2, [1, 2, 3], [1, 2], bad_values, &
                              1, 1, 0.0_dp, .1_dp, status)
    if (status%code /= FORTSPARSE_INVALID_MATRIX) &
        error stop 'nonfinite CSR'
    call preconditioner%setup(2, [1, 3, 2], [1], [1.0_dp], &
                              1, 1, 0.0_dp, .1_dp, status)
    if (status%code /= FORTSPARSE_INVALID_MATRIX) &
        error stop 'invalid row pointers'
    call preconditioner%setup(2, [1, 2, 3], [1, 3], [1.0_dp, 1.0_dp], &
                              1, 1, 0.0_dp, .1_dp, status)
    if (status%code /= FORTSPARSE_INVALID_MATRIX) &
        error stop 'invalid column'
    allocate (large_ptr(65538))
    large_ptr = 1
    call preconditioner%setup(65537, large_ptr, [integer ::], &
                              [real(dp) ::], &
                              65536, 65536, 0.0_dp, .1_dp, status)
    if (status%code /= FORTSPARSE_INDEX_OVERFLOW) &
        error stop 'I4 factor capacity'
    call preconditioner%setup(4, ptr, col, val, 3, 3, 0.0_dp, .1_dp, status)
    call success()
    call preconditioner%free()
    call preconditioner%free()
    info = preconditioner%stats()
    if (info%current_bytes /= 0 .or. info%ready) error stop 'free'
    print '(a)', 'PASS CSR ILUTP action, pivot/drop/fill, caps and recovery'
contains
    subroutine success()
        if (status%code /= 0) then
            print *, status%code, status%msg
            error stop 'unexpected ILUTP failure'
        end if
    end subroutine
end program
