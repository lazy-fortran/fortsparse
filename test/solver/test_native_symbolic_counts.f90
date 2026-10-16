program test_native_symbolic_counts
    use fortsparse_native_symbolic, only: native_symbolic_counts
    use fortsparse_status, only: fortsparse_status_t, status_ok, &
        FORTSPARSE_INVALID_MATRIX
    implicit none
    logical, allocatable :: graph(:, :), reversed(:, :)
    integer :: n, mask, bit, i, j, trial, cases
    type(fortsparse_status_t) :: status
    integer :: empty_parent(0), empty_count(0), one_parent(1), one_count(1)

    cases = 0
    ! Exhaustive independent structural oracle, including every disconnected
    ! graph and permutation through five vertices.
    do n = 1, 5
        allocate (graph(n, n))
        do mask = 0, 2**(n*(n - 1)/2) - 1
            graph = .false.
            bit = 0
            do j = 1, n
                do i = 1, j - 1
                    graph(i, j) = btest(mask, bit)
                    graph(j, i) = graph(i, j)
                    bit = bit + 1
                end do
            end do
            call check_graph(graph)
        end do
        deallocate (graph)
    end do

    ! Deep paths exercise binary lifting beyond the tiny exhaustive corpus.
    ! Dense, star, ladder and disconnected components exercise overlapping
    ! paths/LCAs, fill and multiple forest roots. Reverse every larger graph.
    do n = 17, 129, 112
        allocate (graph(n, n), reversed(n, n))
        do trial = 1, 7
            graph = .false.
            do j = 2, n
                do i = 1, j - 1
                    select case (trial)
                    case (1)
                        graph(i, j) = j == i + 1
                    case (2)
                        graph(i, j) = .true.
                    case (3)
                        graph(i, j) = i == 1
                    case (4)
                        graph(i, j) = j - i == 1 .or. j - i == 3
                    case (5)
                        graph(i, j) = (i - 1)/8 == (j - 1)/8
                    case (6)
                        graph(i, j) = mod(17*i + 31*j + i*j, 53) < 3
                    case (7)
                        graph(i, j) = mod(7*i + 13*j + i*j, 29) < 11
                    end select
                    graph(j, i) = graph(i, j)
                end do
            end do
            call check_graph(graph)
            reversed = graph(n:1:-1, n:1:-1)
            call check_graph(reversed)
        end do
        deallocate (graph, reversed)
    end do
    call native_symbolic_counts(0, [1], [integer ::], &
        empty_parent, empty_count, status)
    if (.not. status_ok(status)) error stop 'empty symbolic forest failed'
    call native_symbolic_counts(1, [1], [1], one_parent, one_count, status)
    if (status%code /= FORTSPARSE_INVALID_MATRIX) &
        error stop 'invalid pointer extent accepted'
    call native_symbolic_counts(1, [0, 2], [1], one_parent, one_count, status)
    if (status%code /= FORTSPARSE_INVALID_MATRIX) &
        error stop 'non-one-based pointer accepted'
    call native_symbolic_counts(1, [1, 3], [1], one_parent, one_count, status)
    if (status%code /= FORTSPARSE_INVALID_MATRIX) &
        error stop 'incompatible row extent accepted'
    print '(a,i0,a)', 'symbolic counts and parents match dense elimination: ', &
        cases, ' graphs'

contains

    subroutine check_graph(adjacency)
        logical, intent(in) :: adjacency(:, :)
        logical :: filled(size(adjacency, 1), size(adjacency, 1))
        integer :: ptr(size(adjacency, 1) + 1), row(size(adjacency))
        integer :: parent(size(adjacency, 1)), counts(size(adjacency, 1))
        integer :: expected_parent(size(adjacency, 1))
        integer :: expected_counts(size(adjacency, 1))
        integer :: order, column, r, c, k, nz
        type(fortsparse_status_t) :: result

        order = size(adjacency, 1)
        nz = 0
        ptr(1) = 1
        do column = 1, order
            do r = 1, order
                if (.not. adjacency(r, column) .and. r /= column) cycle
                nz = nz + 1
                row(nz) = r
            end do
            ptr(column + 1) = nz + 1
        end do
        call native_symbolic_counts(order, ptr, row(:nz), parent, counts, result)
        if (.not. status_ok(result)) then
            print *, 'symbolic status:', result%code, result%msg
            error stop 'symbolic graph failed'
        end if

        ! Independent dense graph elimination: connect all surviving neighbors
        ! into a clique. First surviving neighbor is the elimination parent.
        filled = adjacency
        expected_parent = 0
        expected_counts = 0
        do k = 1, order
            do r = k + 1, order
                if (.not. filled(r, k)) cycle
                expected_counts(k) = expected_counts(k) + 1
                if (expected_parent(k) == 0) expected_parent(k) = r
                do c = r + 1, order
                    if (.not. filled(c, k)) cycle
                    filled(c, r) = .true.
                    filled(r, c) = .true.
                end do
            end do
        end do
        if (any(parent /= expected_parent)) then
            print *, 'actual parent:', parent
            print *, 'oracle parent:', expected_parent
            error stop 'elimination forest mismatch'
        end if
        if (any(counts /= expected_counts)) then
            print *, 'actual counts:', counts
            print *, 'oracle counts:', expected_counts
            error stop 'factor column count mismatch'
        end if
        cases = cases + 1
    end subroutine check_graph
end program test_native_symbolic_counts
