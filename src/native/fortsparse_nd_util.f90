module fortsparse_nd_util
    ! Support of the nested-dissection port (fortsparse_nd): the random
    ! numbers, sort and priority queue of GKlib as METIS 5.2.1 uses them, so
    ! that the port reproduces METIS_NodeND orderings exactly. GKlib is
    ! Copyright 1995-2018 Regents of the University of Minnesota, Apache
    ! License 2.0; this is a translation into Fortran.
    !
    ! Random numbers: METIS 5.2.1 builds GKlib without USE_GKRAND, so its
    ! generator is the C library's rand(). This is glibc's (TYPE_3 additive
    ! feedback, r(i) = r(i-31) + r(i-3) mod 2^32, output r(i)/2), which makes
    ! the orderings those of METIS on Linux on every platform.
    use, intrinsic :: iso_fortran_env, only: int32, int64, real32
    implicit none
    private

    public :: nd_rng_t, nd_srand, nd_rand, nd_rand_in_range, nd_rand_permute
    public :: nd_sort_kv_inc, nd_sort_inc, nd_sort_dec
    public :: rpq_t, rpq_create, rpq_reset, rpq_length, rpq_insert, rpq_delete
    public :: rpq_update, rpq_get_top, rpq_see_top_val, rpq_see_top_key

    integer, parameter :: ik = int32

    type :: nd_rng_t
        integer(int64) :: buf(0:33) = 0
        integer :: next = 0
    end type nd_rng_t

    !> Max-heap of vertices keyed by gains in single precision (GKlib rpq_t).
    type :: rpq_t
        integer(ik) :: nnodes = 0
        real(real32), allocatable :: key(:)
        integer(ik), allocatable :: val(:)
        integer(ik), allocatable :: locator(:)
    end type rpq_t

    integer(int64), parameter :: MASK32 = 4294967295_int64

contains

    !> srand(seed) of glibc: the 31-word table from the minimal standard
    !> generator, then 310 outputs discarded.
    subroutine nd_srand(rng, seed)
        type(nd_rng_t), intent(inout) :: rng
        integer(ik), intent(in) :: seed
        integer(int64) :: word, hi, lo, r(0:343)
        integer :: i
        word = seed
        if (word == 0) word = 1
        r(0) = word
        do i = 1, 30
            hi = word/127773
            lo = mod(word, 127773_int64)
            word = 16807*lo - 2836*hi
            if (word < 0) word = word + 2147483647_int64
            r(i) = word
        end do
        do i = 31, 33
            r(i) = r(i - 31)
        end do
        do i = 34, 343
            r(i) = iand(r(i - 31) + r(i - 3), MASK32)
        end do
        do i = 310, 343
            rng%buf(mod(i, 34)) = r(i)
        end do
        rng%next = 344
    end subroutine nd_srand

    integer(ik) function nd_rand(rng)
        type(nd_rng_t), intent(inout) :: rng
        integer(int64) :: x
        x = iand(rng%buf(mod(rng%next - 31, 34)) + rng%buf(mod(rng%next - 3, 34)), MASK32)
        rng%buf(mod(rng%next, 34)) = x
        rng%next = rng%next + 1
        if (rng%next >= 344 + 34*1000000) rng%next = rng%next - 34*1000000
        nd_rand = int(shiftr(x, 1), ik)
    end function nd_rand

    !> irandInRange: a number in [0, max).
    integer(ik) function nd_rand_in_range(rng, max)
        type(nd_rng_t), intent(inout) :: rng
        integer(ik), intent(in) :: max
        nd_rand_in_range = mod(nd_rand(rng), max)
    end function nd_rand_in_range

    !> irandArrayPermute(n, p, nshuffles, flag) on p(0:n-1).
    subroutine nd_rand_permute(rng, n, p, nshuffles, flag)
        type(nd_rng_t), intent(inout) :: rng
        integer(ik), intent(in) :: n, nshuffles
        integer(ik), intent(inout) :: p(0:)
        integer, intent(in) :: flag
        integer(ik) :: i, u, v, tmp
        if (flag == 1) then
            do i = 0, n - 1
                p(i) = i
            end do
        end if
        if (n < 10) then
            do i = 0, n - 1
                v = nd_rand_in_range(rng, n)
                u = nd_rand_in_range(rng, n)
                tmp = p(v)
                p(v) = p(u)
                p(u) = tmp
            end do
        else
            do i = 0, nshuffles - 1
                v = nd_rand_in_range(rng, n - 3)
                u = nd_rand_in_range(rng, n - 3)
                tmp = p(v)
                p(v) = p(u + 2)
                p(u + 2) = tmp
                tmp = p(v + 1)
                p(v + 1) = p(u + 3)
                p(u + 3) = tmp
                tmp = p(v + 2)
                p(v + 2) = p(u)
                p(u) = tmp
                tmp = p(v + 3)
                p(v + 3) = p(u + 1)
                p(u + 1) = tmp
            end do
        end if
    end subroutine nd_rand_permute

    !> ikvsorti: (key, val)(0:n-1) by increasing key, GKlib's quicksort (not
    !> stable; the order of equal keys is that of METIS).
    subroutine nd_sort_kv_inc(n, key, val)
        integer(ik), intent(in) :: n
        integer(ik), intent(inout) :: key(0:), val(0:)
        call gk_qsort(n, key, val, 1)
    end subroutine nd_sort_kv_inc

    !> isorti on a(0:n-1).
    subroutine nd_sort_inc(n, a)
        integer(ik), intent(in) :: n
        integer(ik), intent(inout) :: a(0:)
        integer(ik), allocatable :: dummy(:)
        allocate (dummy(0:max(0, n - 1)))
        dummy = 0
        call gk_qsort(n, a, dummy, 1)
    end subroutine nd_sort_inc

    !> isortd on a(0:n-1).
    subroutine nd_sort_dec(n, a)
        integer(ik), intent(in) :: n
        integer(ik), intent(inout) :: a(0:)
        integer(ik), allocatable :: dummy(:)
        allocate (dummy(0:max(0, n - 1)))
        dummy = 0
        call gk_qsort(n, a, dummy, -1)
    end subroutine nd_sort_dec

    !> GK_MKQSORT on the pairs (key, val), compared by key: increasing for
    !> dir 1, decreasing for dir -1. Indices stand for the C pointers.
    subroutine gk_qsort(n, key, val, dir)
        integer(ik), intent(in) :: n
        integer(ik), intent(inout) :: key(0:), val(0:)
        integer, intent(in) :: dir
        integer, parameter :: THRESH = 8
        integer(ik) :: lo, hi, mid, left, right, top, run, tmp, thresh_end, trav, h
        integer(ik) :: stack_lo(64), stack_hi(64), hold_k, hold_v
        if (n < 1) return
        if (n > THRESH) then
            lo = 0
            hi = n - 1
            top = 1
            do while (top > 0)
                mid = lo + shiftr(hi - lo, 1)
                if (lt(mid, lo)) call swap(mid, lo)
                if (lt(hi, mid)) then
                    call swap(mid, hi)
                    if (lt(mid, lo)) call swap(mid, lo)
                end if
                left = lo + 1
                right = hi - 1
                do
                    do while (lt(left, mid))
                        left = left + 1
                    end do
                    do while (lt(mid, right))
                        right = right - 1
                    end do
                    if (left < right) then
                        call swap(left, right)
                        if (mid == left) then
                            mid = right
                        else if (mid == right) then
                            mid = left
                        end if
                        left = left + 1
                        right = right - 1
                    else if (left == right) then
                        left = left + 1
                        right = right - 1
                        exit
                    end if
                    if (left > right) exit
                end do
                if (right - lo <= THRESH) then
                    if (hi - left <= THRESH) then
                        ! pop; the bottom entry (top = 1) ends the loop
                        top = top - 1
                        if (top > 0) then
                            lo = stack_lo(top)
                            hi = stack_hi(top)
                        end if
                    else
                        lo = left
                    end if
                else if (hi - left <= THRESH) then
                    hi = right
                else if (right - lo > hi - left) then
                    stack_lo(top) = lo
                    stack_hi(top) = right
                    top = top + 1
                    lo = left
                else
                    stack_lo(top) = left
                    stack_hi(top) = hi
                    top = top + 1
                    hi = right
                end if
            end do
        end if
        ! Insertion sort, after the smallest of the first THRESH+1 elements
        ! is moved to the front.
        thresh_end = min(int(THRESH, ik), n - 1)
        tmp = 0
        do run = 1, thresh_end
            if (lt(run, tmp)) tmp = run
        end do
        if (tmp /= 0) call swap(tmp, 0_ik)
        run = 1
        do
            run = run + 1
            if (run > n - 1) exit
            tmp = run - 1
            do while (lt(run, tmp))
                tmp = tmp - 1
            end do
            tmp = tmp + 1
            if (tmp /= run) then
                trav = run
                hold_k = key(trav)
                hold_v = val(trav)
                h = trav
                do while (h - 1 >= tmp)
                    key(h) = key(h - 1)
                    val(h) = val(h - 1)
                    h = h - 1
                end do
                key(h) = hold_k
                val(h) = hold_v
            end if
        end do
    contains
        logical function lt(a, b)
            integer(ik), intent(in) :: a, b
            if (dir > 0) then
                lt = key(a) < key(b)
            else
                lt = key(a) > key(b)
            end if
        end function lt
        subroutine swap(a, b)
            integer(ik), intent(in) :: a, b
            integer(ik) :: t
            t = key(a)
            key(a) = key(b)
            key(b) = t
            t = val(a)
            val(a) = val(b)
            val(b) = t
        end subroutine swap
    end subroutine gk_qsort

    subroutine rpq_create(q, maxnodes)
        type(rpq_t), intent(inout) :: q
        integer(ik), intent(in) :: maxnodes
        if (allocated(q%key)) deallocate (q%key, q%val, q%locator)
        allocate (q%key(0:max(0, maxnodes - 1)), q%val(0:max(0, maxnodes - 1)), &
            q%locator(0:max(0, maxnodes - 1)))
        q%locator = -1
        q%nnodes = 0
    end subroutine rpq_create

    subroutine rpq_reset(q)
        type(rpq_t), intent(inout) :: q
        integer(ik) :: i
        do i = q%nnodes - 1, 0, -1
            q%locator(q%val(i)) = -1
        end do
        q%nnodes = 0
    end subroutine rpq_reset

    integer(ik) function rpq_length(q)
        type(rpq_t), intent(in) :: q
        rpq_length = q%nnodes
    end function rpq_length

    subroutine rpq_insert(q, node, key)
        type(rpq_t), intent(inout) :: q
        integer(ik), intent(in) :: node
        real(real32), intent(in) :: key
        call heap_insert(q%nnodes, q%key, q%val, q%locator, node, key)
    end subroutine rpq_insert

    subroutine rpq_delete(q, node)
        type(rpq_t), intent(inout) :: q
        integer(ik), intent(in) :: node
        call heap_delete(q%nnodes, q%key, q%val, q%locator, node)
    end subroutine rpq_delete

    subroutine rpq_update(q, node, newkey)
        type(rpq_t), intent(inout) :: q
        integer(ik), intent(in) :: node
        real(real32), intent(in) :: newkey
        call heap_update(q%nnodes, q%key, q%val, q%locator, node, newkey)
    end subroutine rpq_update

    integer(ik) function rpq_get_top(q) result(vtx)
        type(rpq_t), intent(inout) :: q
        vtx = heap_get_top(q%nnodes, q%key, q%val, q%locator)
    end function rpq_get_top

    ! The heap operations on plain arrays, so that the compiler keeps the
    ! array addresses in registers.

    subroutine heap_insert(nnodes, key, val, locator, node, k)
        integer(ik), intent(inout) :: nnodes, val(0:*), locator(0:*)
        real(real32), intent(inout) :: key(0:*)
        integer(ik), intent(in) :: node
        real(real32), intent(in) :: k
        integer(ik) :: i, j
        i = nnodes
        nnodes = nnodes + 1
        do while (i > 0)
            j = shiftr(i - 1, 1)
            if (k > key(j)) then
                key(i) = key(j)
                val(i) = val(j)
                locator(val(i)) = i
                i = j
            else
                exit
            end if
        end do
        key(i) = k
        val(i) = node
        locator(node) = i
    end subroutine heap_insert

    subroutine heap_delete(nnodes, key, val, locator, node_in)
        integer(ik), intent(inout) :: nnodes, val(0:*), locator(0:*)
        real(real32), intent(inout) :: key(0:*)
        integer(ik), intent(in) :: node_in
        integer(ik) :: i, node
        node = node_in
        i = locator(node)
        locator(node) = -1
        nnodes = nnodes - 1
        if (nnodes > 0) then
            if (val(nnodes) /= node) then
                node = val(nnodes)
                call heap_place(nnodes, key, val, locator, i, node, key(nnodes), key(i))
            end if
        end if
    end subroutine heap_delete

    subroutine heap_update(nnodes, key, val, locator, node, newkey)
        integer(ik), intent(inout) :: nnodes, val(0:*), locator(0:*)
        real(real32), intent(inout) :: key(0:*)
        integer(ik), intent(in) :: node
        real(real32), intent(in) :: newkey
        real(real32) :: oldkey
        oldkey = key(locator(node))
        if (.not. (newkey > oldkey) .and. .not. (oldkey > newkey)) return
        call heap_place(nnodes, key, val, locator, locator(node), node, newkey, oldkey)
    end subroutine heap_update

    !> Places node with newkey at heap position i, filtering up when newkey
    !> beats oldkey, else down (GKlib's Delete and Update). newkey and
    !> oldkey are copied first, because they may name heap entries.
    subroutine heap_place(nnodes, key, val, locator, i_in, node, newkey_in, oldkey_in)
        integer(ik), intent(in) :: nnodes, i_in, node
        integer(ik), intent(inout) :: val(0:*), locator(0:*)
        real(real32), intent(inout) :: key(0:*)
        real(real32), value :: newkey_in, oldkey_in
        real(real32) :: newkey
        integer(ik) :: i, j
        i = i_in
        newkey = newkey_in
        if (newkey > oldkey_in) then
            do while (i > 0)
                j = shiftr(i - 1, 1)
                if (newkey > key(j)) then
                    key(i) = key(j)
                    val(i) = val(j)
                    locator(val(i)) = i
                    i = j
                else
                    exit
                end if
            end do
        else
            do
                j = 2*i + 1
                if (j >= nnodes) exit
                if (key(j) > newkey) then
                    if (j + 1 < nnodes) then
                        if (key(j + 1) > key(j)) j = j + 1
                    end if
                else if (j + 1 < nnodes) then
                    if (key(j + 1) > newkey) then
                        j = j + 1
                    else
                        exit
                    end if
                else
                    exit
                end if
                key(i) = key(j)
                val(i) = val(j)
                locator(val(i)) = i
                i = j
            end do
        end if
        key(i) = newkey
        val(i) = node
        locator(node) = i
    end subroutine heap_place

    integer(ik) function heap_get_top(nnodes, key, val, locator) result(vtx)
        integer(ik), intent(inout) :: nnodes, val(0:*), locator(0:*)
        real(real32), intent(inout) :: key(0:*)
        integer(ik) :: i, j, node
        real(real32) :: k
        if (nnodes == 0) then
            vtx = -1
            return
        end if
        nnodes = nnodes - 1
        vtx = val(0)
        locator(vtx) = -1
        i = nnodes
        if (i > 0) then
            k = key(i)
            node = val(i)
            i = 0
            do
                j = 2*i + 1
                if (j >= nnodes) exit
                if (key(j) > k) then
                    if (j + 1 < nnodes) then
                        if (key(j + 1) > key(j)) j = j + 1
                    end if
                else if (j + 1 < nnodes) then
                    if (key(j + 1) > k) then
                        j = j + 1
                    else
                        exit
                    end if
                else
                    exit
                end if
                key(i) = key(j)
                val(i) = val(j)
                locator(val(i)) = i
                i = j
            end do
            key(i) = k
            val(i) = node
            locator(node) = i
        end if
    end function heap_get_top

    integer(ik) function rpq_see_top_val(q)
        type(rpq_t), intent(in) :: q
        if (q%nnodes == 0) then
            rpq_see_top_val = -1
        else
            rpq_see_top_val = q%val(0)
        end if
    end function rpq_see_top_val

    real(real32) function rpq_see_top_key(q)
        type(rpq_t), intent(in) :: q
        if (q%nnodes == 0) then
            rpq_see_top_key = huge(1.0_real32)
        else
            rpq_see_top_key = q%key(0)
        end if
    end function rpq_see_top_key
end module fortsparse_nd_util
