program test_structural_matching
    use fortsparse, only: dp, fortsparse_status_t, status_ok, &
                          FORTSPARSE_SINGULAR, &
                          FORTSPARSE_INVALID_MATRIX
    use fortsparse_native_matching, only: structural_row_matching, &
                                          matching_graph_rows
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    integer :: mask, subset, i, j, e, neighbors, columns, nmask, previous
    integer :: ptr(4), rows(9), pref(3), inverse(3), work(15), graph(9)
    integer :: p4(5), r4(16), m4(4), i4(4), w4(20)
    real(dp) :: v4(16)
    integer :: p8(9), r8(15), m8(8), i8(8), w8(40), g8(15)
    real(dp) :: values(9), v8(15)
    complex(dp) :: z(9)
    logical :: perfect
    type(fortsparse_status_t) :: status
    ptr = [1, 4, 7, 10]
    rows = [1, 2, 3, 1, 2, 3, 1, 2, 3]
    ! Exhaust every three-by-three graph. Hall's subset criterion is an
    ! independent mathematical oracle, not another augmenting-path algorithm.
    do mask = 0, 511
        do e = 1, 9
            values(e) = 0
            if (btest(mask, e - 1)) values(e) = 1
        end do
        perfect = .true.
        do subset = 1, 7
            nmask = 0; columns = popcnt(subset)
            do j = 1, 3
                if (.not. btest(subset, j - 1)) cycle
                do i = 1, 3
                    if (values(3*(j - 1) + i) /= 0) nmask = ibset(nmask, i - 1)
                end do
            end do
            neighbors = popcnt(nmask)
            if (neighbors < columns) perfect = .false.
        end do
        call structural_row_matching(3, ptr, rows, values, pref, inverse, &
                                     work, status)
        call check(status_ok(status) .eqv. perfect, 'Hall oracle for real graph')
        if (perfect) then
            do j = 1, 3
                call check(values(ptr(j) + pref(j) - 1) /= 0, &
                           'matched edge nonzero')
                call check(inverse(pref(j)) == j, 'inverse bijection')
            end do
            call matching_graph_rows(ptr, rows, inverse, graph)
            do j = 1, 3
                call check(all(graph(ptr(j):ptr(j + 1) - 1) == [1, 2, 3]), &
                           'sorted graph labels')
            end do
        else
            call check(status%code == FORTSPARSE_SINGULAR, &
                       'Hall failure explicit singular')
        end if
        z = cmplx(0.0_dp, values, dp)
        call structural_row_matching(3, ptr, rows, z, pref, inverse, work, status)
        call check(status_ok(status) .eqv. perfect, &
                   'imaginary-only edges remain active')
    end do
    ! Requires an eight-column alternating augmentation after the initial
    ! shortest-path phase: columns1..7 have rows j,j+1; column8 only row1.
    e = 1; p8(1) = 1
    do j = 1, 7
        r8(e:e + 1) = [j, j + 1]; e = e + 2; p8(j + 1) = e
    end do
    r8(e) = 1; p8(9) = e + 1; v8 = 1
    call structural_row_matching(8, p8, r8, v8, m8, i8, w8, status)
    call check(status_ok(status), 'deep iterative augmenting path')
    call check(all(m8 == [2, 3, 4, 5, 6, 7, 8, 1]), &
               'deterministic noninvolutory match')
    call matching_graph_rows(p8, r8, i8, g8)
    do j = 1, 8
        previous = 0
        do e = p8(j), p8(j + 1) - 1
            call check(g8(e) > previous, 'permuted graph remains canonical')
            previous = g8(e)
        end do
    end do
    ! Hall deficiency with no empty row or column: first three columns
    ! share only rows1,2, while column4 reaches both rows3,4.
    p4 = [1, 5, 9, 13, 17]; v4 = 0
    do j = 1, 4
        r4(4*(j - 1) + 1:4*j) = [1, 2, 3, 4]
        if (j <= 3) then
            v4(4*(j - 1) + 1:4*(j - 1) + 2) = 1
        else
            v4(15:16) = 1
        end if
    end do
    call structural_row_matching(4, p4, r4, v4, m4, i4, w4, status)
    call check(status%code == FORTSPARSE_SINGULAR, &
               'nonempty Hall-deficient graph')
    values = 1; values(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call structural_row_matching(3, ptr, rows, values, pref, inverse, work, &
                                 status)
    call check(status%code == FORTSPARSE_INVALID_MATRIX, &
               'nonfinite real edge rejected')
    z = cmplx(1.0_dp, 0.0_dp, dp); z(1) = cmplx(1.0_dp, values(1), dp)
    call structural_row_matching(3, ptr, rows, z, pref, inverse, work, status)
    call check(status%code == FORTSPARSE_INVALID_MATRIX, &
               'nonfinite imaginary edge rejected')
    values = 1; rows(2) = 1
    call structural_row_matching(3, ptr, rows, values, pref, inverse, work, &
                                 status)
    call check(status%code == FORTSPARSE_INVALID_MATRIX, &
               'duplicate CSC rows rejected')
contains
    subroutine check(ok, message)
        logical, intent(in)::ok
        character(*), intent(in)::message
        if (.not. ok) then
            print *, message
            error stop 1
        end if
    end subroutine
end program
