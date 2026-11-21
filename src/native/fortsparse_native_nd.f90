module fortsparse_native_nd
    ! The native nested-dissection ordering behind the C entry points
    ! fsparse_native_order (symmetric pattern, every off-diagonal entry
    ! mirrored) and fsparse_native_order_general (the pattern of A + A^T):
    ! canonical square CSC with one-based indices in, perm(new) = old and
    ! inverse(old) = new, one-based, out. Status: 0 OK, 2 invalid matrix,
    ! 6 out of memory.
    use, intrinsic :: iso_c_binding, only: c_int
    use, intrinsic :: iso_fortran_env, only: int32, int64
    use fortsparse_nd, only: nd_order
    implicit none
    private

    public :: native_nd_order

    integer(c_int), parameter :: ORDER_OK = 0, ORDER_INVALID = 2, ORDER_NO_MEMORY = 6

contains

    integer(c_int) function fsparse_native_order_available() &
            bind(C, name='fsparse_native_order_available')
        fsparse_native_order_available = 1
    end function fsparse_native_order_available

    integer(c_int) function fsparse_native_order(n, nnz, ptr, row, perm, inverse) &
            bind(C, name='fsparse_native_order')
        integer(c_int), value :: n, nnz
        integer(c_int), intent(in) :: ptr(*), row(*)
        integer(c_int), intent(out) :: perm(*), inverse(*)
        fsparse_native_order = native_nd_order(n, nnz, ptr, row, perm, inverse, .false.)
    end function fsparse_native_order

    integer(c_int) function fsparse_native_order_general(n, nnz, ptr, row, perm, inverse) &
            bind(C, name='fsparse_native_order_general')
        integer(c_int), value :: n, nnz
        integer(c_int), intent(in) :: ptr(*), row(*)
        integer(c_int), intent(out) :: perm(*), inverse(*)
        fsparse_native_order_general = native_nd_order(n, nnz, ptr, row, perm, inverse, .true.)
    end function fsparse_native_order_general

    !> The ordering of the graph of A (general: of A + A^T), self edges
    !> dropped.
    integer(c_int) function native_nd_order(n, nnz, ptr, row, perm, inverse, general) &
            result(status)
        integer(c_int), intent(in) :: n, nnz, ptr(*), row(*)
        integer(c_int), intent(out) :: perm(*), inverse(*)
        logical, intent(in) :: general
        integer(int32), allocatable :: xadj(:), adj(:), next(:), p(:), ip(:)
        integer(int32) :: col, entry, r, count, previous, begin, finish
        integer :: io

        status = ORDER_INVALID
        if (.not. csc_valid(n, nnz, ptr, row)) return
        if (.not. general) then
            if (.not. mirrored(n, ptr, row)) return
        end if
        if (int(nnz, int64)*merge(2, 1, general) > int(huge(1_int32), int64)) then
            status = ORDER_NO_MEMORY
            return
        end if
        allocate (xadj(0:n), adj(0:max(1, merge(2, 1, general)*nnz) - 1), p(0:n - 1), &
            ip(0:n - 1), stat=io)
        if (io /= 0) then
            status = ORDER_NO_MEMORY
            return
        end if
        count = 0
        if (general) then
            allocate (next(0:n - 1), stat=io)
            if (io /= 0) then
                status = ORDER_NO_MEMORY
                return
            end if
            next = 0
            do col = 0, n - 1
                do entry = ptr(col + 1), ptr(col + 2) - 1
                    r = row(entry) - 1
                    if (r == col) cycle
                    next(col) = next(col) + 1
                    next(r) = next(r) + 1
                end do
            end do
            xadj(0) = 0
            do col = 0, n - 1
                xadj(col + 1) = xadj(col) + next(col)
                next(col) = xadj(col)
            end do
            do col = 0, n - 1
                do entry = ptr(col + 1), ptr(col + 2) - 1
                    r = row(entry) - 1
                    if (r == col) cycle
                    adj(next(col)) = r
                    next(col) = next(col) + 1
                    adj(next(r)) = col
                    next(r) = next(r) + 1
                end do
            end do
            ! Sorted, without repeats.
            do col = 0, n - 1
                begin = xadj(col)
                finish = xadj(col + 1)
                call sort_ints(adj(begin:finish - 1))
                xadj(col) = count
                previous = -1
                do entry = begin, finish - 1
                    r = adj(entry)
                    if (r /= previous) then
                        adj(count) = r
                        count = count + 1
                    end if
                    previous = r
                end do
            end do
        else
            do col = 0, n - 1
                xadj(col) = count
                do entry = ptr(col + 1), ptr(col + 2) - 1
                    if (row(entry) /= col + 1) then
                        adj(count) = row(entry) - 1
                        count = count + 1
                    end if
                end do
            end do
        end if
        xadj(n) = count
        call nd_order(n, xadj, adj, p, ip)
        perm(1:n) = p + 1
        inverse(1:n) = ip + 1
        status = ORDER_OK
    end function native_nd_order

    !> Ascending in place (Shell sort, Ciura's gaps extended by 2.25; no workspace).
    subroutine sort_ints(a)
        integer(int32), intent(inout) :: a(:)
        integer, parameter :: gaps(19) = [5818032, 2585792, 1149241, 510774, 227011, 100894, &
            44842, 19930, 8858, 3937, 1750, 701, 301, 132, 57, 23, 10, 4, 1]
        integer :: g, gap, i, j
        integer(int32) :: t
        do g = 1, size(gaps)
            gap = gaps(g)
            do i = gap + 1, size(a)
                t = a(i)
                j = i
                do while (j > gap)
                    if (a(j - gap) <= t) exit
                    a(j) = a(j - gap)
                    j = j - gap
                end do
                a(j) = t
            end do
        end do
    end subroutine sort_ints

    !> Canonical square CSC with one-based, strictly increasing rows.
    logical function csc_valid(n, nnz, ptr, row) result(valid)
        integer(c_int), intent(in) :: n, nnz, ptr(*), row(*)
        integer(c_int) :: col, entry, previous
        valid = .false.
        if (n <= 0 .or. nnz < 0 .or. nnz == huge(1_c_int)) return
        if (ptr(1) /= 1 .or. ptr(n + 1) /= nnz + 1) return
        do col = 1, n
            if (ptr(col) < 1 .or. ptr(col) > ptr(col + 1) .or. ptr(col + 1) > nnz + 1) return
            previous = 0
            do entry = ptr(col), ptr(col + 1) - 1
                if (row(entry) < 1 .or. row(entry) > n .or. row(entry) <= previous) return
                previous = row(entry)
            end do
        end do
        valid = .true.
    end function csc_valid

    !> Every off-diagonal entry (r, c) has its mirror (c, r).
    logical function mirrored(n, ptr, row)
        integer(c_int), intent(in) :: n, ptr(*), row(*)
        integer(c_int) :: col, entry, lo, hi, mid, wanted, other
        mirrored = .false.
        do col = 1, n
            do entry = ptr(col), ptr(col + 1) - 1
                other = row(entry)
                if (other == col) cycle
                ! Binary search for col in column other.
                wanted = col
                lo = ptr(other)
                hi = ptr(other + 1)
                do while (lo < hi)
                    mid = lo + (hi - lo)/2
                    if (row(mid) < wanted) then
                        lo = mid + 1
                    else
                        hi = mid
                    end if
                end do
                if (lo >= ptr(other + 1)) return
                if (row(lo) /= wanted) return
            end do
        end do
        mirrored = .true.
    end function mirrored
end module fortsparse_native_nd
