module fortsparse_native_matching
    ! Discrete Hopcroft--Karp analysis of CSC columns versus original rows.
    ! All storage is supplied by the caller; numerical values are borrowed.
    use fortsparse_kinds, only: dp
    use fortsparse_status, only: fortsparse_status_t, status_set, &
                                 FORTSPARSE_OK, FORTSPARSE_INVALID_MATRIX, &
                                 FORTSPARSE_SINGULAR
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: structural_row_matching, matching_graph_rows
    interface structural_row_matching
        module procedure match_real, match_complex
    end interface
contains
    subroutine match_real(n, ptr, rows, values, preference, inverse, work, &
                          status)
        integer, intent(in) :: n, ptr(:), rows(:)
        real(dp), intent(in) :: values(:)
        integer, intent(out) :: preference(:), inverse(:), work(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: e
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        'matching: invalid real values')
        if (size(values) /= size(rows)) return
        do e = 1, size(values)
            if (.not. ieee_is_finite(values(e))) return
        end do
        call match_graph(n, ptr, rows, preference, inverse, work, status, &
                         real_values=values)
    end subroutine

    subroutine match_complex(n, ptr, rows, values, preference, inverse, work, &
                             status)
        integer, intent(in) :: n, ptr(:), rows(:)
        complex(dp), intent(in) :: values(:)
        integer, intent(out) :: preference(:), inverse(:), work(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: e
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        'matching: invalid complex values')
        if (size(values) /= size(rows)) return
        do e = 1, size(values)
            if (.not. ieee_is_finite(real(values(e), dp))) return
            if (.not. ieee_is_finite(aimag(values(e)))) return
        end do
        call match_graph(n, ptr, rows, preference, inverse, work, status, &
                         complex_values=values)
    end subroutine

    subroutine match_graph(n, ptr, rows, preference, inverse, work, status, &
                           real_values, complex_values)
        integer, intent(in) :: n, ptr(:), rows(:)
        integer, intent(out) :: preference(:), inverse(:)
        integer, target, intent(out) :: work(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), optional, intent(in) :: real_values(:)
        complex(dp), optional, intent(in) :: complex_values(:)
        integer, pointer :: level(:), queue(:), cursor(:), stack(:), chosen(:)
        integer :: j, e, row, previous, head, tail, shortest, mate, root, top, &
                   d, matched
        logical :: advanced, augmented
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        'matching: require canonical square CSC and workspace')
        if (n < 0 .or. n > huge(0)/5) return
        if (size(ptr) /= n + 1 .or. size(preference) /= n .or. &
            size(inverse) /= n) return
        if (size(work) /= 5*n) return
        if (ptr(n + 1) < 1) return
        if (ptr(1) /= 1 .or. ptr(n + 1) - 1 /= size(rows)) return
        do j = 1, n
            if (ptr(j) < 1) return
            if (ptr(j + 1) < ptr(j)) return
            if (ptr(j + 1) - 1 > size(rows)) return
            previous = 0
            do e = ptr(j), ptr(j + 1) - 1
                row = rows(e)
                if (row <= previous .or. row > n) return
                previous = row
            end do
        end do
        level => work(1:n); queue => work(n + 1:2*n); cursor => work(2*n + 1:3*n)
        stack => work(3*n + 1:4*n); chosen => work(4*n + 1:5*n)
        preference = 0; inverse = 0; matched = 0
        do while (matched < n)
            level = -1; head = 1; tail = 0; shortest = huge(0)
            do j = 1, n
                if (preference(j) /= 0) cycle
                tail = tail + 1; queue(tail) = j; level(j) = 0
            end do
            do while (head <= tail)
                j = queue(head); head = head + 1
                if (level(j) >= shortest) cycle
                do e = ptr(j), ptr(j + 1) - 1
                    if (.not. nonzero(e)) cycle
                    mate = inverse(rows(e))
                    if (mate == 0) then
                        shortest = level(j) + 1
                    else if (level(mate) < 0) then
                        level(mate) = level(j) + 1
                        tail = tail + 1; queue(tail) = mate
                    end if
                end do
            end do
            if (shortest == huge(0)) exit
            cursor = ptr(:n)
            do root = 1, n
                if (preference(root) /= 0 .or. level(root) /= 0) cycle
                top = 1; stack(1) = root; augmented = .false.
                do while (top > 0)
                    j = stack(top); advanced = .false.
                    do while (cursor(j) < ptr(j + 1))
                        e = cursor(j); cursor(j) = e + 1
                        if (.not. nonzero(e)) cycle
                        row = rows(e); mate = inverse(row)
                        if (mate == 0) then
                            if (level(j) + 1 /= shortest) cycle
                            chosen(top) = row
                            do d = top, 1, -1
                                preference(stack(d)) = chosen(d)
                                inverse(chosen(d)) = stack(d)
                            end do
                            matched = matched + 1; augmented = .true.
                            exit
                        else if (level(mate) == level(j) + 1) then
                            chosen(top) = row; top = top + 1; stack(top) = mate
                            advanced = .true.
                            exit
                        end if
                    end do
                    if (augmented) exit
                    if (advanced) cycle
                    level(j) = -1; top = top - 1
                end do
            end do
        end do
        if (matched /= n) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                            'matching: nonzero graph has no perfect matching')
            return
        end if
        call status_set(status, FORTSPARSE_OK, '')
    contains
        logical function nonzero(entry)
            integer, intent(in) :: entry
            if (present(real_values)) then
                nonzero = real_values(entry) /= 0.0_dp
            else
                nonzero = complex_values(entry) /= cmplx(0.0_dp, 0.0_dp, dp)
            end if
        end function
    end subroutine

    subroutine matching_graph_rows(ptr, rows, inverse, graph)
        ! AMD needs sorted row labels of P*A. Heap sort each borrowed column;
        ! no numerical CSC copy and no degree-dependent automatic workspace.
        integer, intent(in) :: ptr(:), rows(:), inverse(:)
        integer, intent(out) :: graph(:)
        integer :: e, j, first, count, start, finish, root, child, value
        do e = 1, size(rows)
            graph(e) = inverse(rows(e))
        end do
        do j = 1, size(ptr) - 1
            first = ptr(j) - 1; count = ptr(j + 1) - ptr(j)
            do start = count/2, 1, -1
                call sift(start, count)
            end do
            do finish = count, 2, -1
                value = graph(first + finish)
                graph(first + finish) = graph(first + 1)
                graph(first + 1) = value
                call sift(1, finish - 1)
            end do
        end do
    contains
        subroutine sift(begin, last)
            integer, intent(in) :: begin, last
            root = begin; value = graph(first + root)
            do while (root <= last/2)
                child = 2*root
                if (child < last) then
                    if (graph(first + child) < graph(first + child + 1)) then
                        child = child + 1
                    end if
                end if
                if (value >= graph(first + child)) exit
                graph(first + root) = graph(first + child); root = child
            end do
            graph(first + root) = value
        end subroutine
    end subroutine
end module
