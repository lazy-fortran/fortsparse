module fortsparse_amd_i32
    ! The AMD port (fortsparse_amd_kernel.inc) on 32-bit indices.
    use, intrinsic :: iso_c_binding, only: c_int
    use, intrinsic :: iso_fortran_env, only: int32, int64
    implicit none
    private
    public :: amd_aat, amd_graph, amd_order_core
    integer, parameter :: ik = int32
    include 'fortsparse_amd_kernel.inc'
end module fortsparse_amd_i32

module fortsparse_amd_i64
    ! The AMD port on 64-bit indices, for graphs beyond the 32-bit range.
    use, intrinsic :: iso_c_binding, only: c_int
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    private
    public :: amd_aat, amd_graph, amd_order_core
    integer, parameter :: ik = int64
    include 'fortsparse_amd_kernel.inc'
end module fortsparse_amd_i64

module fortsparse_amd
    ! Native AMD ordering of one-based CSC patterns: the ordering of A+A^T,
    ! equal to SuiteSparse amd_order's. The bind(C) names are the ordering
    ! ABI shared with the C METIS adapter and the tests.
    use, intrinsic :: iso_c_binding, only: c_int, c_int32_t, c_int64_t
    use, intrinsic :: iso_fortran_env, only: int32, int64
    use fortsparse_amd_i32, only: aat32 => amd_aat, graph32 => amd_graph, &
        core32 => amd_order_core
    use fortsparse_amd_i64, only: aat64 => amd_aat, graph64 => amd_graph, &
        core64 => amd_order_core
    implicit none
    private
    public :: csc_graph_valid

    ! Status codes of the ordering ABI.
    integer(c_int), parameter :: FS_OK = 0, FS_INVALID = 2, FS_INTERNAL = 5, &
        FS_MEMORY = 6, FS_INDEX = 7
contains

    ! Canonical square CSC with one-based, strictly increasing rows.
    logical function csc_graph_valid(n, nnz, ptr, row) result(valid)
        integer(c_int), intent(in) :: n, nnz
        integer(c_int), intent(in) :: ptr(*), row(*)
        integer :: col, entry, previous

        valid = .false.
        if (n <= 0 .or. nnz < 0 .or. nnz == huge(nnz)) return
        if (ptr(1) /= 1 .or. ptr(n + 1) /= nnz + 1) return
        do col = 1, n
            if (ptr(col) < 1 .or. ptr(col) > ptr(col + 1) .or. &
                ptr(col + 1) > nnz + 1) return
            previous = 0
            do entry = ptr(col), ptr(col + 1) - 1
                if (row(entry) < 1 .or. row(entry) > n .or. &
                    row(entry) <= previous) return
                previous = row(entry)
            end do
        end do
        valid = .true.
    end function csc_graph_valid

    logical function has_row(ptr, row, column, wanted)
        integer(c_int), intent(in) :: ptr(*), row(*), column, wanted
        integer :: lo, hi, mid

        lo = ptr(column)
        hi = ptr(column + 1)
        do while (lo < hi)
            mid = lo + (hi - lo) / 2
            if (row(mid) < wanted) then
                lo = mid + 1
            else
                hi = mid
            end if
        end do
        has_row = lo < ptr(column + 1)
        if (has_row) has_row = row(lo) == wanted
    end function has_row

    ! The SPD contract: every off-diagonal edge is mirrored.
    logical function graph_valid(n, nnz, ptr, row) result(valid)
        integer(c_int), intent(in) :: n, nnz, ptr(*), row(*)
        integer :: col, entry

        valid = csc_graph_valid(n, nnz, ptr, row)
        if (.not. valid) return
        do col = 1, n
            do entry = ptr(col), ptr(col + 1) - 1
                if (row(entry) /= col) then
                    if (.not. has_row(ptr, row, row(entry), col)) then
                        valid = .false.
                        return
                    end if
                end if
            end do
        end do
    end function graph_valid

    ! Arena words (7n + iwlen) and the A+A^T entry count m, or a status.
    integer(c_int) function plan(n, nnz, ptr, row, low, m, iwlen, words) &
            result(code)
        integer(c_int), intent(in) :: n, nnz, low
        integer(c_int), intent(in) :: ptr(*), row(*)
        integer(int64), intent(out) :: m, iwlen, words
        integer(int64) :: pairs
        integer :: j, p, i

        m = 0
        iwlen = 0
        words = 0
        code = FS_INVALID
        if (n <= 0 .or. nnz < 0 .or. (low /= 0 .and. low /= 1)) return
        code = FS_INDEX
        ! AMD's front-end size restriction; every arena offset is checked too.
        if (n >= huge(n) / 4 .or. nnz >= huge(nnz) / 4) return
        code = FS_INVALID
        if (.not. csc_graph_valid(n, nnz, ptr, row)) return
        pairs = 0
        do j = 1, n
            do p = ptr(j), ptr(j + 1) - 1
                i = row(p)
                if (i == j) cycle
                if (i < j) then
                    pairs = pairs + 1
                else if (.not. has_row(ptr, row, i, j)) then
                    pairs = pairs + 1
                end if
            end do
        end do
        m = 2 * pairs
        iwlen = m + n
        if (low == 0) iwlen = iwlen + m / 5
        words = 7 * int(n, int64) + iwlen
        code = FS_INDEX
        if (m >= huge(0_int32) .or. iwlen >= huge(0_int32) .or. &
            words > huge(0_int64) / 4) return
        code = FS_OK
    end function plan

    integer(c_int) function bounded_available() result(available) &
            bind(C, name='fsparse_native_amd_bounded_available')
        available = 1
    end function bounded_available

    ! Arena and output bytes for fsparse_native_amd_bounded; outputs are left
    ! unchanged on failure.
    integer(c_int) function bounded_size(n, nnz, ptr, row, low, &
            workspace_bytes, output_bytes) result(code) &
            bind(C, name='fsparse_native_amd_bounded_size')
        integer(c_int), value :: n, nnz, low
        integer(c_int), intent(in) :: ptr(*), row(*)
        integer(c_int64_t), intent(inout) :: workspace_bytes, output_bytes
        integer(int64) :: m, iwlen, words

        code = plan(n, nnz, ptr, row, low, m, iwlen, words)
        if (code /= FS_OK) return
        workspace_bytes = 4 * words
        output_bytes = 8 * int(n, int64)
    end function bounded_size

    ! AMD on a caller arena of int32 words; allocates nothing. low = 1 trades
    ! the elbow room for more graph compressions. perm(new) = old and
    ! inverse(old) = new, one-based.
    integer(c_int) function bounded_order(n, nnz, ptr, row, low, workspace, &
            workspace_bytes, perm, inverse) result(code) &
            bind(C, name='fsparse_native_amd_bounded')
        integer(c_int), value :: n, nnz, low
        integer(c_int), intent(in) :: ptr(*), row(*)
        integer(c_int32_t), intent(inout) :: workspace(*)
        integer(c_int64_t), value :: workspace_bytes
        integer(c_int), intent(out) :: perm(*), inverse(*)

        code = FS_INVALID
        if (n <= 0 .or. nnz < 0 .or. (low /= 0 .and. low /= 1)) return
        code = FS_INDEX
        if (n >= huge(n) / 4 .or. nnz >= huge(nnz) / 4) return
        code = FS_INVALID
        if (.not. csc_graph_valid(n, nnz, ptr, row)) return
        if (workspace_bytes < 28 * int(n, int64)) return
        ! Element arguments associate each array with its arena segment.
        code = arena_order(n, ptr, row, low, workspace_bytes / 4, &
            workspace(1), workspace(n + 1), workspace(2 * n + 1), &
            workspace(3 * n + 1), workspace(4 * n + 1), workspace(5 * n + 1), &
            workspace(6 * n + 1), workspace(7 * n + 1), perm, inverse)
    end function bounded_order

    ! The arena path on validated input: the A+A^T counts first, then the
    ! elbow room the arena must hold.
    integer(c_int) function arena_order(n, ptr, row, low, words, pe, len, nv, &
            head, elen, degree, w, iw, perm, inverse) result(code)
        integer(c_int), intent(in) :: n, ptr(*), row(*), low
        integer(int64), intent(in) :: words
        integer(int32), intent(inout) :: pe(n), len(n), nv(n), head(n), &
            elen(n), degree(n), w(n), iw(*)
        integer(c_int), intent(out) :: perm(n), inverse(n)
        integer(int64) :: m, iwlen

        m = aat32(n, ptr, row, len, w)
        iwlen = m + n
        if (low == 0) iwlen = iwlen + m / 5
        code = FS_INDEX
        if (m >= huge(0_int32) .or. iwlen >= huge(0_int32)) return
        code = FS_INVALID
        if (words < 7 * int(n, int64) + iwlen) return
        call graph32(n, ptr, row, len, pe, iw, nv, w)
        code = order_graph32(n, int(m, int32), int(iwlen, int32), pe, len, nv, &
            head, elen, degree, w, iw, perm, inverse)
    end function arena_order

    integer(c_int) function order_graph32(n, m, iwlen, pe, len, nv, head, elen, &
            degree, w, iw, perm, inverse) result(code)
        integer(c_int), intent(in) :: n
        integer(int32), intent(in) :: m, iwlen
        integer(int32), intent(inout) :: pe(n), len(n), nv(n), head(n), &
            elen(n), degree(n), w(n), iw(iwlen)
        integer(c_int), intent(out) :: perm(n), inverse(n)
        integer :: i

        code = FS_INTERNAL
        if (core32(n, pe, iw, len, iwlen, m, nv, inverse, perm, head, elen, &
            degree, w) /= 0) return
        do i = 1, n
            perm(i) = perm(i) + 1
            inverse(i) = inverse(i) + 1
        end do
        do i = 1, n
            if (perm(i) < 1 .or. perm(i) > n .or. inverse(i) < 1 .or. &
                inverse(i) > n) return
            if (perm(inverse(i)) /= i) return
        end do
        code = FS_OK
    end function order_graph32

    ! 64-bit indices, for graphs whose A+A^T pattern exceeds the 32-bit range.
    integer(c_int) function order_wide(n, ptr, row, perm, inverse) result(code)
        integer(c_int), intent(in) :: n, ptr(*), row(*)
        integer(c_int), intent(out) :: perm(n), inverse(n)
        integer(int64), allocatable :: len(:), tp(:), pe(:), iw(:), nv(:), &
            next(:), last(:), head(:), elen(:), degree(:), w(:)
        integer(int64) :: nzaat, iwlen, n64
        integer :: status

        code = FS_MEMORY
        n64 = n
        allocate (len(n), tp(n), stat=status)
        if (status /= 0) return
        nzaat = aat64(n64, ptr, row, len, tp)
        iwlen = nzaat + nzaat / 5 + n64
        allocate (pe(n), iw(iwlen), nv(n), next(n), last(n), head(n), elen(n), &
            degree(n), w(n), stat=status)
        if (status /= 0) return
        call graph64(n64, ptr, row, len, pe, iw, nv, tp)
        code = FS_INTERNAL
        if (core64(n64, pe, iw, len, iwlen, nzaat, nv, next, last, head, elen, &
            degree, w) /= 0) return
        perm = int(last + 1, c_int)
        inverse = int(next + 1, c_int)
        code = FS_OK
    end function order_wide

    ! Validated input on heap arrays: one pass for the A+A^T counts, one to
    ! list the graph. perm(new) = old; inverse(old) = new; all one-based.
    integer(c_int) function order_heap(n, ptr, row, perm, inverse) result(code)
        integer(c_int), intent(in) :: n, ptr(*), row(*)
        integer(c_int), intent(out) :: perm(*), inverse(*)
        integer(int32), allocatable :: pe(:), len(:), nv(:), head(:), elen(:), &
            degree(:), w(:), iw(:)
        integer(int64) :: m, iwlen
        integer :: status

        code = FS_INDEX
        if (n >= huge(n) / 4) return
        code = FS_MEMORY
        allocate (len(n), w(n), stat=status)
        if (status /= 0) return
        m = aat32(n, ptr, row, len, w)
        iwlen = m + m / 5 + n
        if (m >= huge(0_int32) .or. iwlen >= huge(0_int32)) then
            code = order_wide(n, ptr, row, perm, inverse)
            return
        end if
        allocate (pe(n), nv(n), head(n), elen(n), degree(n), iw(iwlen), &
            stat=status)
        if (status /= 0) return
        call graph32(n, ptr, row, len, pe, iw, nv, w)
        code = order_graph32(n, int(m, int32), int(iwlen, int32), pe, len, nv, &
            head, elen, degree, w, iw, perm, inverse)
    end function order_heap

    integer(c_int) function order_available() result(available) &
            bind(C, name='fsparse_native_order_amd_available')
        available = 1
    end function order_available

    ! Structurally symmetric input, as the SPD factorizations supply.
    integer(c_int) function order_symmetric(n, nnz, ptr, row, perm, inverse) &
            result(code) bind(C, name='fsparse_native_order_amd')
        integer(c_int), value :: n, nnz
        integer(c_int), intent(in) :: ptr(*), row(*)
        integer(c_int), intent(out) :: perm(*), inverse(*)

        code = FS_INVALID
        if (.not. graph_valid(n, nnz, ptr, row)) return
        code = order_heap(n, ptr, row, perm, inverse)
    end function order_symmetric

    ! Unsymmetric input; LU uses the ordering of A+A^T as its column order and
    ! keeps its numerical row pivots.
    integer(c_int) function order_general(n, nnz, ptr, row, perm, inverse) &
            result(code) bind(C, name='fsparse_native_order_amd_general')
        integer(c_int), value :: n, nnz
        integer(c_int), intent(in) :: ptr(*), row(*)
        integer(c_int), intent(out) :: perm(*), inverse(*)

        code = FS_INVALID
        if (.not. csc_graph_valid(n, nnz, ptr, row)) return
        code = order_heap(n, ptr, row, perm, inverse)
    end function order_general

    ! The 64-bit path on any size, for tests: it orders like the 32-bit one.
    integer(c_int) function order_general_wide(n, nnz, ptr, row, perm, inverse) &
            result(code) bind(C, name='fsparse_native_order_amd_general_wide')
        integer(c_int), value :: n, nnz
        integer(c_int), intent(in) :: ptr(*), row(*)
        integer(c_int), intent(out) :: perm(*), inverse(*)

        code = FS_INVALID
        if (.not. csc_graph_valid(n, nnz, ptr, row)) return
        code = order_wide(n, ptr, row, perm, inverse)
    end function order_general_wide
end module fortsparse_amd
