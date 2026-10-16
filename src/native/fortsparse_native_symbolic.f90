module fortsparse_native_symbolic
    ! Graph topology only; no numerical values or generated arithmetic leaves.
    use fortsparse_kinds, only: i8
    use fortsparse_status, only: fortsparse_status_t, status_set, &
        FORTSPARSE_OK, FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OUT_OF_MEMORY, &
        FORTSPARSE_INTERNAL_ERROR
    implicit none
    private
    public :: native_symbolic_counts, native_symbolic_workspace_bytes

contains

    pure integer function symbolic_levels(n) result(levels)
        integer, intent(in) :: n
        integer :: span
        levels = 0
        span = n
        do while (span > 1)
            levels = levels + 1
            span = span/2
        end do
    end function symbolic_levels

    pure integer(i8) function native_symbolic_workspace_bytes(n) result(bytes)
        ! Actual temporary array payload below; no caller arrays or RSS overhead.
        integer, intent(in) :: n
        integer(i8) :: columns
        bytes = 0_i8
        if (n <= 0) return
        columns = int(n, i8)
        bytes = 4_i8*(5_i8*columns + 2_i8 + &
            (columns + 1_i8)*(symbolic_levels(n) + 1_i8)) + &
            8_i8*(columns + 1_i8)
    end function native_symbolic_workspace_bytes

    subroutine native_symbolic_counts(n, col_ptr, row_idx, parent, count, status)
        ! Caller owns canonical, sorted, symmetric, one-based CSC validation.
        ! parent is the elimination forest (zero roots); count excludes L's
        ! diagonal. All temporary storage is released on return.
        integer, intent(in) :: n, col_ptr(:), row_idx(:)
        integer, intent(out) :: parent(:), count(:)
        type(fortsparse_status_t), intent(out) :: status
        integer, allocatable :: ancestor(:), sibling(:), preorder(:), depth(:)
        integer, allocatable :: up(:, :), starts(:)
        integer(i8), allocatable :: delta(:)
        integer :: levels, ierr, j, p, i, previous, level, top, clock
        integer :: child, number, common

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            'native symbolic: incompatible array extents')
        if (n < 0 .or. n == huge(n)) return
        if (size(col_ptr) /= n + 1) return
        if (size(parent) /= n .or. size(count) /= n) return
        if (col_ptr(1) /= 1) return
        if (col_ptr(n + 1) < 1) return
        if (size(row_idx) /= col_ptr(n + 1) - 1) return
        parent = 0
        count = 0
        if (n == 0) then
            call status_set(status, FORTSPARSE_OK, '')
            return
        end if
        levels = symbolic_levels(n)
        allocate (ancestor(0:n), sibling(n), preorder(n), depth(0:n), &
            up(0:n, 0:levels), starts(n), delta(0:n), stat=ierr)
        if (ierr /= 0) then
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, &
                'native symbolic: setup allocation failed')
            return
        end if

        ! Standard upper-triangle elimination tree with ancestor compression.
        ancestor = 0
        do j = 1, n
            do p = col_ptr(j), col_ptr(j + 1) - 1
                i = row_idx(p)
                if (i >= j) exit
                do while (i /= 0)
                    if (i >= j) exit
                    previous = ancestor(i)
                    ancestor(i) = j
                    if (previous == 0) parent(i) = j
                    i = previous
                end do
            end do
        end do

        ! Reuse ancestor as child-list heads. Iterative DFS handles deep trees.
        ancestor = 0
        do i = 1, n
            sibling(i) = ancestor(parent(i))
            ancestor(parent(i)) = i
        end do
        top = 0
        child = ancestor(0)
        do while (child /= 0)
            top = top + 1
            starts(top) = child
            child = sibling(child)
        end do
        clock = 0
        do while (top > 0)
            i = starts(top)
            top = top - 1
            clock = clock + 1
            preorder(i) = clock
            child = ancestor(i)
            do while (child /= 0)
                top = top + 1
                starts(top) = child
                child = sibling(child)
            end do
        end do

        ! Parent indices strictly increase. Row zero is valid sentinel padding
        ! for all lifting levels; no indexed expression needs short circuiting.
        depth = 0
        up = 0
        do i = n, 1, -1
            depth(i) = depth(parent(i)) + 1
            up(i, 0) = parent(i)
        end do
        do level = 1, levels
            do i = 1, n
                up(i, level) = up(up(i, level - 1), level - 1)
            end do
        end do

        ! Each column contributes one unit along the union of paths from its
        ! upper neighbors to j, excluding j. Adjacent DFS LCAs remove overlap.
        delta = 0_i8
        do j = 1, n
            number = 0
            do p = col_ptr(j), col_ptr(j + 1) - 1
                i = row_idx(p)
                if (i >= j) exit
                number = number + 1
                starts(number) = i
                delta(i) = delta(i) + 1_i8
            end do
            if (number == 0) cycle
            call sort_preorder(starts, number, preorder)
            do p = 2, number
                common = lowest_common_ancestor(starts(p - 1), starts(p), &
                    depth, up, levels)
                delta(common) = delta(common) - 1_i8
            end do
            delta(j) = delta(j) - 1_i8
        end do
        do i = 1, n
            if (delta(i) < 0_i8 .or. delta(i) > int(n - i, i8)) then
                call status_set(status, FORTSPARSE_INTERNAL_ERROR, &
                    'native symbolic: inconsistent path-union count')
                return
            end if
            count(i) = int(delta(i))
            delta(parent(i)) = delta(parent(i)) + delta(i)
        end do
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine native_symbolic_counts

    integer function lowest_common_ancestor(first, second, depth, up, levels) &
            result(common)
        integer, intent(in) :: first, second, depth(0:), up(0:, 0:), levels
        integer :: a, b, difference, level

        a = first
        b = second
        if (depth(a) < depth(b)) then
            a = second
            b = first
        end if
        difference = depth(a) - depth(b)
        do level = 0, levels
            if (btest(difference, level)) a = up(a, level)
        end do
        if (a == b) then
            common = a
            return
        end if
        do level = levels, 0, -1
            if (up(a, level) /= up(b, level)) then
                a = up(a, level)
                b = up(b, level)
            end if
        end do
        common = up(a, 0)
    end function lowest_common_ancestor

    subroutine sort_preorder(vertices, number, preorder)
        ! In-place heapsort; no allocation, recursion or quadratic worst case.
        integer, intent(inout) :: vertices(:)
        integer, intent(in) :: number, preorder(:)
        integer :: first, last, temporary

        do first = number/2, 1, -1
            call sift(first, number)
        end do
        do last = number, 2, -1
            temporary = vertices(1)
            vertices(1) = vertices(last)
            vertices(last) = temporary
            call sift(1, last - 1)
        end do

    contains

        subroutine sift(first, last)
            integer, intent(in) :: first, last
            integer :: root, child, value

            root = first
            value = vertices(root)
            do while (root <= last/2)
                child = 2*root
                if (child < last) then
                    if (preorder(vertices(child)) < &
                        preorder(vertices(child + 1))) child = child + 1
                end if
                if (preorder(value) >= preorder(vertices(child))) exit
                vertices(root) = vertices(child)
                root = child
            end do
            vertices(root) = value
        end subroutine sift
    end subroutine sort_preorder
end module fortsparse_native_symbolic
