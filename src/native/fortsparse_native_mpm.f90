module fortsparse_native_mpm
    ! Maximum-product transversal with optimal LP duals and the derived
    ! row/column scaling (the MC64 job-5 problem of Duff & Koster, 2001),
    ! implemented independently as shortest augmenting paths (Dijkstra with a
    ! binary heap) on the costs
    !     c(i,j) = log(colmax(j)) - log|a(i,j)| >= 0,
    ! taken over numerically nonzero entries only; complex entries use |a|.
    !
    ! Output convention (column -> row): cmatch(j) is the original row matched
    ! to column j. The row permutation Pm of design section 2 places the
    ! matched entries on the diagonal, (Pm*A)(j,:) = A(cmatch(j),:). On return
    ! iwork(1:n) holds the inverse, iwork(i) = column matched to row i.
    !
    ! Certificate (duals u, v indexed by original row and column):
    !     c(i,j) - u(i) - v(j) >= 0 on every nonzero entry,
    !     c(cmatch(j),j) - u(cmatch(j)) - v(j) = 0,
    ! hence sum_j c(cmatch(j),j) = sum(u) + sum(v) and the matching is optimal
    ! by LP duality (up to rounding).
    !
    ! All storage is supplied by the caller: iwork(mpm_iwork_size(n)) and
    ! rwork(mpm_rwork_size(n, nnz, cache_costs)). A real workspace that also
    ! holds nnz entries caches the costs; otherwise costs are recomputed, with
    ! bitwise identical results. Nothing is allocated; the procedures are
    ! reentrant.
    use fortsparse_kinds, only: dp, i8
    use fortsparse_memory, only: memory_checked_add, memory_checked_mul
    use fortsparse_status, only: fortsparse_status_t, status_set, &
                                 FORTSPARSE_OK, FORTSPARSE_INVALID_MATRIX, &
                                 FORTSPARSE_SINGULAR
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: mpm_stats_t, mpm_iwork_size, mpm_rwork_size, mpm_workspace_bytes
    public :: max_product_matching, mpm_scaling

    type :: mpm_stats_t
        ! Columns matched by the initial tight-edge pass.
        integer :: greedy_matched = 0
        ! Shortest augmenting paths found (= n - greedy_matched on success).
        integer :: augmentations = 0
        ! Rows finalized (heap pops) and edge relaxations over all searches.
        integer(i8) :: heap_pops = 0_i8, relaxations = 0_i8
        ! log10(max/min) of dr = exp(u) and dc = exp(v)/colmax.
        real(dp) :: log_range_dr = 0.0_dp, log_range_dc = 0.0_dp
    end type

    interface max_product_matching
        module procedure mpm_real, mpm_complex
    end interface

    interface mpm_scaling
        module procedure scaling_real, scaling_complex
    end interface

    integer, parameter :: IWORK_PER_N = 6, RWORK_PER_N = 2
contains
    pure integer(i8) function mpm_iwork_size(n) result(count)
        ! Integer workspace elements; -1 when n is invalid for 32-bit indexing.
        integer, intent(in) :: n
        count = -1_i8
        if (n < 0 .or. n > huge(0)/IWORK_PER_N) return
        count = int(IWORK_PER_N, i8)*int(n, i8)
    end function

    pure integer(i8) function mpm_rwork_size(n, nnz, cache_costs) result(count)
        ! Real workspace elements (dist, log colmax and optionally the costs);
        ! -1 when n or nnz is invalid for 32-bit indexing.
        integer, intent(in) :: n, nnz
        logical, intent(in) :: cache_costs
        count = -1_i8
        if (mpm_iwork_size(n) < 0_i8 .or. nnz < 0) return
        count = int(RWORK_PER_N, i8)*int(n, i8)
        if (cache_costs) count = count + int(nnz, i8)
        if (count > int(huge(0), i8)) count = -1_i8
    end function

    pure integer(i8) function mpm_workspace_bytes(n, nnz, cache_costs) &
        result(bytes)
        ! Exact caller workspace payload in bytes, 4*iwork + 8*rwork, or -1 for
        ! invalid sizes. The outputs cmatch(n), u(n), v(n) add 4n + 16n bytes,
        ! giving the 7n*I + 4n*R (+ nnz*R) MPM phase of the design (section 7).
        integer, intent(in) :: n, nnz
        logical, intent(in) :: cache_costs
        integer(i8) :: ints, reals, ibytes, rbytes
        type(fortsparse_status_t) :: status
        bytes = -1_i8
        ints = mpm_iwork_size(n)
        reals = mpm_rwork_size(n, nnz, cache_costs)
        if (ints < 0_i8 .or. reals < 0_i8) return
        ibytes = 0_i8; rbytes = 0_i8
        call memory_checked_mul(ints, int(storage_size(0)/8, i8), ibytes, status)
        if (status%code /= FORTSPARSE_OK) return
        call memory_checked_mul(reals, int(storage_size(0.0_dp)/8, i8), rbytes, &
                                status)
        if (status%code /= FORTSPARSE_OK) return
        call memory_checked_add(ibytes, rbytes, bytes, status)
        if (status%code /= FORTSPARSE_OK) bytes = -1_i8
    end function

    subroutine mpm_real(n, ptr, rows, values, cmatch, u, v, iwork, rwork, &
                        status, stats)
        integer, intent(in) :: n, ptr(:), rows(:)
        real(dp), intent(in) :: values(:)
        integer, intent(out) :: cmatch(:)
        real(dp), intent(out) :: u(:), v(:)
        integer, contiguous, intent(out) :: iwork(:)
        real(dp), contiguous, intent(out) :: rwork(:)
        type(fortsparse_status_t), intent(out) :: status
        type(mpm_stats_t), intent(out) :: stats
        integer :: e
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        'mpm: invalid real values')
        if (size(values) /= size(rows)) return
        do e = 1, size(values)
            if (.not. ieee_is_finite(values(e))) return
        end do
        call mpm_driver(n, ptr, rows, cmatch, u, v, iwork, rwork, status, stats, &
                        real_values=values)
    end subroutine

    subroutine mpm_complex(n, ptr, rows, values, cmatch, u, v, iwork, rwork, &
                           status, stats)
        integer, intent(in) :: n, ptr(:), rows(:)
        complex(dp), intent(in) :: values(:)
        integer, intent(out) :: cmatch(:)
        real(dp), intent(out) :: u(:), v(:)
        integer, contiguous, intent(out) :: iwork(:)
        real(dp), contiguous, intent(out) :: rwork(:)
        type(fortsparse_status_t), intent(out) :: status
        type(mpm_stats_t), intent(out) :: stats
        integer :: e
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        'mpm: invalid complex values')
        if (size(values) /= size(rows)) return
        do e = 1, size(values)
            if (.not. ieee_is_finite(real(values(e), dp))) return
            if (.not. ieee_is_finite(aimag(values(e)))) return
            ! |a| must be representable for a finite cost; the halved modulus
            ! is tested so that the check itself cannot overflow.
            if (abs(0.5_dp*values(e)) > 0.5_dp*huge(1.0_dp)) return
        end do
        call mpm_driver(n, ptr, rows, cmatch, u, v, iwork, rwork, status, stats, &
                        complex_values=values)
    end subroutine

    subroutine mpm_driver(n, ptr, rows, cmatch, u, v, iwork, rwork, status, &
                          stats, real_values, complex_values)
        integer, intent(in) :: n, ptr(:), rows(:)
        integer, intent(out) :: cmatch(:)
        real(dp), intent(out) :: u(:), v(:)
        integer, contiguous, intent(out) :: iwork(:)
        real(dp), contiguous, intent(out) :: rwork(:)
        type(fortsparse_status_t), intent(out) :: status
        type(mpm_stats_t), intent(out) :: stats
        real(dp), optional, intent(in) :: real_values(:)
        complex(dp), optional, intent(in) :: complex_values(:)
        integer :: j, e, row, previous, nnz
        logical :: cache
        stats = mpm_stats_t()
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        'mpm: require canonical square CSC and workspace')
        if (mpm_iwork_size(n) < 0_i8) return
        if (size(ptr) /= n + 1 .or. size(cmatch) /= n) return
        if (size(u) /= n .or. size(v) /= n) return
        if (ptr(1) /= 1 .or. ptr(n + 1) - 1 /= size(rows)) return
        nnz = size(rows)
        if (size(iwork, kind=i8) /= mpm_iwork_size(n)) return
        if (size(rwork, kind=i8) == mpm_rwork_size(n, nnz, .true.) .and. &
            nnz > 0) then
            cache = .true.
        else if (size(rwork, kind=i8) == mpm_rwork_size(n, nnz, .false.)) then
            cache = .false.
        else
            return
        end if
        do j = 1, n
            if (ptr(j + 1) < ptr(j)) return
            previous = 0
            do e = ptr(j), ptr(j + 1) - 1
                row = rows(e)
                if (row <= previous .or. row > n) return
                previous = row
            end do
        end do
        if (n == 0) then
            call status_set(status, FORTSPARSE_OK, '')
            return
        end if
        if (cache) then
            call mpm_core(n, ptr, rows, cmatch, u, v, iwork(1:n), &
                          iwork(n + 1:2*n), iwork(2*n + 1:3*n), &
                          iwork(3*n + 1:4*n), iwork(4*n + 1:5*n), &
                          iwork(5*n + 1:6*n), rwork(1:n), rwork(n + 1:2*n), &
                          rwork(2*n + 1:2*n + nnz), .true., status, stats, &
                          real_values, complex_values)
        else
            call mpm_core(n, ptr, rows, cmatch, u, v, iwork(1:n), &
                          iwork(n + 1:2*n), iwork(2*n + 1:3*n), &
                          iwork(3*n + 1:4*n), iwork(4*n + 1:5*n), &
                          iwork(5*n + 1:6*n), rwork(1:n), rwork(n + 1:2*n), &
                          rwork(1:0), .false., status, stats, real_values, &
                          complex_values)
        end if
    end subroutine

    subroutine mpm_core(n, ptr, rows, cmatch, u, v, rmatch, heap, hpos, &
                        touched, pred, state, dist, lcmax, cost, cache, status, &
                        stats, real_values, complex_values)
        integer, intent(in) :: n, ptr(:), rows(:)
        integer, intent(out) :: cmatch(:)
        real(dp), intent(out) :: u(:), v(:)
        integer, intent(out) :: rmatch(n), heap(n), hpos(n), touched(n), &
                                pred(n), state(n)
        real(dp), intent(out) :: dist(n), lcmax(n), cost(*)
        logical, intent(in) :: cache
        type(fortsparse_status_t), intent(inout) :: status
        type(mpm_stats_t), intent(inout) :: stats
        real(dp), optional, intent(in) :: real_values(:)
        complex(dp), optional, intent(in) :: complex_values(:)
        ! state: 0 untouched, 1 queued in the heap, 2 distance final.
        integer, parameter :: UNSEEN = 0, QUEUED = 1, FINAL = 2
        real(dp), parameter :: INF = huge(1.0_dp)
        real(dp) :: colmax, c, bound, slack, lo, hi
        integer :: i, j, e, j0, hn, nt, ifree, col, row, matched

        ! Column maxima in log form; numerically zero entries are excluded.
        do j = 1, n
            colmax = 0.0_dp
            do e = ptr(j), ptr(j + 1) - 1
                colmax = max(colmax, magnitude(e))
            end do
            if (colmax == 0.0_dp) then
                call singular(status, 'mpm: column ', j, &
                              ' has no numerically nonzero entry')
                return
            end if
            lcmax(j) = log(colmax)
        end do
        if (cache) then
            do j = 1, n
                do e = ptr(j), ptr(j + 1) - 1
                    cost(e) = edge_cost(e, j)
                end do
            end do
        end if

        ! Feasible initial duals: u(i) = row minimum, v(j) = reduced column
        ! minimum. The edge attaining v(j) is tight in exact arithmetic of the
        ! same expression, so the greedy pass matches it when its row is free.
        u = INF
        do j = 1, n
            do e = ptr(j), ptr(j + 1) - 1
                c = cost_of(e, j)
                if (c < INF) u(rows(e)) = min(u(rows(e)), c)
            end do
        end do
        do i = 1, n
            if (u(i) == INF) then
                call singular(status, 'mpm: row ', i, &
                              ' has no numerically nonzero entry')
                return
            end if
        end do
        rmatch = 0; cmatch = 0; matched = 0
        do j = 1, n
            v(j) = INF
            do e = ptr(j), ptr(j + 1) - 1
                c = cost_of(e, j)
                if (c < INF) v(j) = min(v(j), c - u(rows(e)))
            end do
            do e = ptr(j), ptr(j + 1) - 1
                i = rows(e)
                if (rmatch(i) /= 0) cycle
                c = cost_of(e, j)
                if (c == INF) cycle
                if (c - u(i) - v(j) <= 0.0_dp) then
                    rmatch(i) = j; cmatch(j) = i; matched = matched + 1
                    exit
                end if
            end do
        end do
        stats%greedy_matched = matched

        ! One shortest augmenting path per free column.
        dist = INF; state = UNSEEN; hpos = 0
        do j0 = 1, n
            if (cmatch(j0) /= 0) cycle
            hn = 0; nt = 0; ifree = 0; bound = INF
            call scan_column(j0, 0.0_dp)
            do while (hn > 0)
                i = pop()
                state(i) = FINAL
                stats%heap_pops = stats%heap_pops + 1_i8
                if (rmatch(i) == 0) then
                    ifree = i
                    exit
                end if
                call scan_column(rmatch(i), dist(i))
            end do
            if (ifree == 0) then
                call singular(status, 'mpm: structurally singular, no ' // &
                              'augmenting path from column ', j0, '', matched, &
                              ' columns matched before')
                return
            end if
            ! Dual update on finalized rows (distances <= dist(ifree)) and their
            ! matched columns keeps every reduced cost nonnegative.
            bound = dist(ifree)
            do e = 1, nt
                i = touched(e)
                if (state(i) /= FINAL) cycle
                slack = bound - dist(i)
                if (slack <= 0.0_dp) cycle
                u(i) = u(i) - slack
                if (rmatch(i) /= 0) v(rmatch(i)) = v(rmatch(i)) + slack
            end do
            v(j0) = v(j0) + bound
            i = ifree
            do
                col = pred(i); row = cmatch(col)
                cmatch(col) = i; rmatch(i) = col
                if (col == j0) exit
                i = row
            end do
            matched = matched + 1
            stats%augmentations = stats%augmentations + 1
            do e = 1, nt
                i = touched(e)
                dist(i) = INF; state(i) = UNSEEN; hpos(i) = 0
            end do
        end do

        lo = minval(u); hi = maxval(u)
        stats%log_range_dr = (hi - lo)/log(10.0_dp)
        lo = INF; hi = -INF
        do j = 1, n
            lo = min(lo, v(j) - lcmax(j)); hi = max(hi, v(j) - lcmax(j))
        end do
        stats%log_range_dc = (hi - lo)/log(10.0_dp)
        call status_set(status, FORTSPARSE_OK, '')
    contains
        real(dp) function magnitude(entry)
            integer, intent(in) :: entry
            if (present(real_values)) then
                magnitude = abs(real_values(entry))
            else
                magnitude = abs(complex_values(entry))
            end if
        end function

        real(dp) function edge_cost(entry, column)
            integer, intent(in) :: entry, column
            real(dp) :: a
            a = magnitude(entry)
            if (a == 0.0_dp) then
                edge_cost = INF
            else
                edge_cost = lcmax(column) - log(a)
            end if
        end function

        real(dp) function cost_of(entry, column)
            integer, intent(in) :: entry, column
            if (cache) then
                cost_of = cost(entry)
            else
                cost_of = edge_cost(entry, column)
            end if
        end function

        subroutine scan_column(jc, base)
            ! Relax every edge of column jc reached at distance base. Rows at
            ! or beyond the best free-row distance found so far are pruned.
            integer, intent(in) :: jc
            real(dp), intent(in) :: base
            integer :: q, r
            real(dp) :: cq, nd
            do q = ptr(jc), ptr(jc + 1) - 1
                r = rows(q)
                if (state(r) == FINAL) cycle
                cq = cost_of(q, jc)
                if (cq == INF) cycle
                stats%relaxations = stats%relaxations + 1_i8
                nd = base + max(0.0_dp, cq - u(r) - v(jc))
                if (nd >= bound) cycle
                if (state(r) == UNSEEN) then
                    state(r) = QUEUED; nt = nt + 1; touched(nt) = r
                    dist(r) = nd; pred(r) = jc
                    call push(r)
                else if (nd < dist(r)) then
                    dist(r) = nd; pred(r) = jc
                    call sift_up(hpos(r))
                end if
                if (rmatch(r) == 0) bound = min(bound, nd)
            end do
        end subroutine

        subroutine push(r)
            integer, intent(in) :: r
            hn = hn + 1; heap(hn) = r; hpos(r) = hn
            call sift_up(hn)
        end subroutine

        integer function pop() result(r)
            r = heap(1); hpos(r) = 0
            heap(1) = heap(hn); hn = hn - 1
            if (hn > 0) then
                hpos(heap(1)) = 1
                call sift_down(1)
            end if
        end function

        subroutine sift_up(k0)
            integer, intent(in) :: k0
            integer :: k, parent, t
            k = k0
            do while (k > 1)
                parent = k/2
                if (dist(heap(parent)) <= dist(heap(k))) exit
                t = heap(parent); heap(parent) = heap(k); heap(k) = t
                hpos(heap(parent)) = parent; hpos(heap(k)) = k
                k = parent
            end do
        end subroutine

        subroutine sift_down(k0)
            integer, intent(in) :: k0
            integer :: k, child, t
            k = k0
            do
                child = 2*k
                if (child > hn) exit
                if (child < hn) then
                    if (dist(heap(child + 1)) < dist(heap(child))) child = child + 1
                end if
                if (dist(heap(k)) <= dist(heap(child))) exit
                t = heap(child); heap(child) = heap(k); heap(k) = t
                hpos(heap(child)) = child; hpos(heap(k)) = k
                k = child
            end do
        end subroutine
    end subroutine

    pure subroutine singular(status, head, first, middle, second, tail)
        ! SINGULAR message with integers, formatted without runtime I/O or
        ! allocating concatenation (the no-allocation contract covers errors).
        type(fortsparse_status_t), intent(inout) :: status
        character(*), intent(in) :: head, middle
        integer, intent(in) :: first
        integer, intent(in), optional :: second
        character(*), intent(in), optional :: tail
        character(len=len(status%msg)) :: text
        integer :: at
        text = ''; at = 0
        call append(text, at, head)
        call append_int(text, at, first)
        call append(text, at, middle)
        if (present(second)) then
            call append(text, at, ' (')
            call append_int(text, at, second)
            if (present(tail)) call append(text, at, tail)
            call append(text, at, ')')
        end if
        call status_set(status, FORTSPARSE_SINGULAR, text(1:at))
    end subroutine

    pure subroutine append(text, at, piece)
        character(*), intent(inout) :: text
        integer, intent(inout) :: at
        character(*), intent(in) :: piece
        integer :: k
        do k = 1, len(piece)
            if (at >= len(text)) return
            at = at + 1; text(at:at) = piece(k:k)
        end do
    end subroutine

    pure subroutine append_int(text, at, value)
        character(*), intent(inout) :: text
        integer, intent(inout) :: at
        integer, intent(in) :: value
        character(len=12) :: digits
        integer :: k
        integer(i8) :: rest
        k = len(digits); rest = abs(int(value, i8))
        do
            digits(k:k) = achar(iachar('0') + int(mod(rest, 10_i8)))
            rest = rest/10_i8
            if (rest == 0_i8) exit
            k = k - 1
        end do
        if (value < 0) then
            k = k - 1; digits(k:k) = '-'
        end if
        call append(text, at, digits(k:))
    end subroutine

    pure subroutine scaling_real(n, ptr, rows, values, u, v, dr, dc, status)
        ! dr(i) = exp(u(i)) scales ORIGINAL row i and dc(j) = exp(v(j))/colmax(j)
        ! scales column j, so |dr(i)*a(i,j)*dc(j)| = exp(-reduced cost) <= 1 with
        ! equality on matched entries. In the permuted system Dr*Pm*A*Dc of
        ! design section 2 the diagonal of Dr is dr(cmatch(:)), i.e. the
        ! permutation is applied after the row scaling: Dr*Pm = Pm*diag(dr).
        integer, intent(in) :: n, ptr(:), rows(:)
        real(dp), intent(in) :: values(:), u(:), v(:)
        real(dp), intent(out) :: dr(:), dc(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: j, e
        real(dp) :: colmax
        call check_scaling_sizes(n, ptr, size(rows), size(values), u, v, dr, &
                                 dc, status)
        if (status%code /= FORTSPARSE_OK) return
        do j = 1, n
            colmax = 0.0_dp
            do e = ptr(j), ptr(j + 1) - 1
                colmax = max(colmax, abs(values(e)))
            end do
            call scale_column(colmax, v(j), dc(j), status)
            if (status%code /= FORTSPARSE_OK) return
        end do
        call scale_rows(u, dr, status)
    end subroutine

    pure subroutine scaling_complex(n, ptr, rows, values, u, v, dr, dc, status)
        ! Complex twin of scaling_real with magnitudes |a(i,j)|.
        integer, intent(in) :: n, ptr(:), rows(:)
        complex(dp), intent(in) :: values(:)
        real(dp), intent(in) :: u(:), v(:)
        real(dp), intent(out) :: dr(:), dc(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: j, e
        real(dp) :: colmax
        call check_scaling_sizes(n, ptr, size(rows), size(values), u, v, dr, &
                                 dc, status)
        if (status%code /= FORTSPARSE_OK) return
        do j = 1, n
            colmax = 0.0_dp
            do e = ptr(j), ptr(j + 1) - 1
                colmax = max(colmax, abs(values(e)))
            end do
            call scale_column(colmax, v(j), dc(j), status)
            if (status%code /= FORTSPARSE_OK) return
        end do
        call scale_rows(u, dr, status)
    end subroutine

    pure subroutine check_scaling_sizes(n, ptr, nrows, nvalues, u, v, dr, dc, &
                                        status)
        integer, intent(in) :: n, ptr(:), nrows, nvalues
        real(dp), intent(in) :: u(:), v(:), dr(:), dc(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: j
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        'mpm scaling: inconsistent CSC or vector sizes')
        if (n < 0 .or. size(ptr) /= n + 1 .or. nvalues /= nrows) return
        if (size(u) /= n .or. size(v) /= n) return
        if (size(dr) /= n .or. size(dc) /= n) return
        if (ptr(1) /= 1 .or. ptr(n + 1) - 1 /= nrows) return
        do j = 1, n
            if (ptr(j + 1) < ptr(j)) return
        end do
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine

    pure subroutine scale_column(colmax, vj, dcj, status)
        real(dp), intent(in) :: colmax, vj
        real(dp), intent(out) :: dcj
        type(fortsparse_status_t), intent(inout) :: status
        dcj = 0.0_dp
        if (colmax == 0.0_dp) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                            'mpm scaling: numerically zero column')
            return
        end if
        dcj = exp(vj - log(colmax))
        if (.not. ieee_is_finite(dcj) .or. dcj == 0.0_dp) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                            'mpm scaling: column scale not representable')
        end if
    end subroutine

    pure subroutine scale_rows(u, dr, status)
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: dr(:)
        type(fortsparse_status_t), intent(inout) :: status
        integer :: i
        do i = 1, size(u)
            dr(i) = exp(u(i))
            if (.not. ieee_is_finite(dr(i)) .or. dr(i) == 0.0_dp) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                                'mpm scaling: row scale not representable')
                return
            end if
        end do
    end subroutine
end module
