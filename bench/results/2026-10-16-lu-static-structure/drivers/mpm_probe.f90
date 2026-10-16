! Prototype maximum-product matching with dual scaling (independent
! implementation of the shortest-augmenting-path / Dijkstra Hungarian method,
! cf. Duff & Koster 2001 MC64 job 5; no HSL code used).
! Cost c_ij = log(max_k |a_kj|) - log|a_ij| >= 0 on numerically nonzero entries.
! Duals: c_ij - u_i - v_j >= 0, equality on matched edges.
! Scaling: Dr = exp(u), Dc = exp(v)/colmax makes |Dr A Dc| <= 1, matched = 1.
! Usage: mpm_probe matrix.mtx out_prefix
!   writes out_prefix-fpref.bin (int32, 1-based matched row per column)
!          out_prefix-fdual.bin (real64 u(1:n), v(1:n))
program mpm_probe
    use, intrinsic :: iso_fortran_env, only: real64, int64, int32
    use fortsparse_csc, only: csc_t, csc_from_triplet
    use fortsparse_status, only: fortsparse_status_t
    implicit none
    integer, parameter :: dp = real64, i8 = int64
    type(csc_t) :: a
    type(fortsparse_status_t) :: status
    integer :: n, m, nz, e, uu, i, j, p, j0, row, col, nfree, scans_total
    integer, allocatable :: rowc(:), colc(:), rmatch(:), cmatch(:)
    integer, allocatable :: heap(:), hpos(:), touched(:), pred(:), state(:)
    real(dp), allocatable :: val(:), cost(:), u(:), v(:), dist(:), colmax(:)
    real(dp) :: big, d, nd, dmax, minslack, maxms, obj, red
    integer(i8) :: t0, t1, rate, relax
    integer :: hn, nt, ifree
    character(512) :: path, prefix

    call get_command_argument(1, path)
    call get_command_argument(2, prefix)
    open(newunit=uu, file=trim(path), status='old')
    block
        character(512) :: line
        do
            read(uu,'(a)') line
            if (line(1:1) /= '%') exit
        end do
        read(line,*) n, m, nz
    end block
    allocate(rowc(nz), colc(nz), val(nz))
    do e = 1, nz
        read(uu,*) rowc(e), colc(e), val(e)
    end do
    close(uu)
    call csc_from_triplet(n, n, rowc, colc, val, a, status)
    if (status%code /= 0) error stop 'invalid matrix'
    call system_clock(t0, rate)
    big = huge(1.0_dp)
    allocate(cost(size(a%row_idx)), colmax(n), u(n), v(n), dist(n), rmatch(n), cmatch(n))
    allocate(heap(n), hpos(n), touched(n), pred(n), state(n))
    ! Cost; numerically zero entries are excluded (cost = +inf).
    do j = 1, n
        colmax(j) = 0
        do p = a%col_ptr(j), a%col_ptr(j+1)-1
            colmax(j) = max(colmax(j), abs(a%val(p)))
        end do
        if (colmax(j) == 0) error stop 'structurally/numerically zero column'
        do p = a%col_ptr(j), a%col_ptr(j+1)-1
            if (a%val(p) == 0) then
                cost(p) = big
            else
                cost(p) = log(colmax(j)) - log(abs(a%val(p)))
            end if
        end do
    end do
    ! Feasible initial duals: u = row minimum, v = column reduced minimum.
    u = big
    do j = 1, n
        do p = a%col_ptr(j), a%col_ptr(j+1)-1
            u(a%row_idx(p)) = min(u(a%row_idx(p)), cost(p))
        end do
    end do
    if (any(u == big)) error stop 'zero row'
    rmatch = 0; cmatch = 0
    do j = 1, n
        v(j) = big
        do p = a%col_ptr(j), a%col_ptr(j+1)-1
            if (cost(p) < big) v(j) = min(v(j), cost(p) - u(a%row_idx(p)))
        end do
        ! Greedy: match a free row on a tight edge.
        do p = a%col_ptr(j), a%col_ptr(j+1)-1
            i = a%row_idx(p)
            if (cost(p) < big .and. rmatch(i) == 0) then
                if (cost(p) - u(i) - v(j) <= 0) then
                    rmatch(i) = j; cmatch(j) = i; exit
                end if
            end if
        end do
    end do
    nfree = count(cmatch == 0)
    ! Shortest augmenting paths (Dijkstra on reduced costs) for free columns.
    dist = big; state = 0; hpos = 0; relax = 0; scans_total = 0
    do j0 = 1, n
        if (cmatch(j0) /= 0) cycle
        hn = 0; nt = 0; ifree = 0; dmax = big
        call scan_column(j0, 0.0_dp)
        do while (hn > 0)
            i = pop()
            state(i) = 2
            scans_total = scans_total + 1
            if (rmatch(i) == 0) then
                ifree = i; dmax = dist(i); exit
            end if
            call scan_column(rmatch(i), dist(i))
        end do
        if (ifree == 0) error stop 'structurally singular: no augmenting path'
        ! Potential update over touched vertices with dist < dmax.
        do e = 1, nt
            i = touched(e)
            if (state(i) == 2 .and. dist(i) < dmax) then
                u(i) = u(i) - (dmax - dist(i))
                if (rmatch(i) /= 0) v(rmatch(i)) = v(rmatch(i)) + (dmax - dist(i))
            end if
        end do
        v(j0) = v(j0) + dmax
        ! Augment along pred.
        i = ifree
        do
            col = pred(i)
            row = cmatch(col)
            cmatch(col) = i; rmatch(i) = col
            if (col == j0) exit
            i = row
        end do
        ! Reset touched rows.
        do e = 1, nt
            i = touched(e)
            dist(i) = big; state(i) = 0; hpos(i) = 0
        end do
    end do
    call system_clock(t1)
    ! Independent checks: permutation, dual feasibility, complementary slackness.
    if (any(cmatch == 0) .or. any(rmatch == 0)) error stop 'incomplete matching'
    minslack = big; maxms = 0; obj = 0
    do j = 1, n
        do p = a%col_ptr(j), a%col_ptr(j+1)-1
            if (cost(p) == big) cycle
            red = cost(p) - u(a%row_idx(p)) - v(j)
            minslack = min(minslack, red)
            if (a%row_idx(p) == cmatch(j)) then
                maxms = max(maxms, abs(red)); obj = obj + cost(p)
            end if
        end do
    end do
    write(*,'(a,*(g0,:,","))') 'MPM,n,greedy_free,dijkstra_row_scans,matching_s_unqualified,obj_sum_cost,dual_obj,min_reduced,max_matched_reduced: ', &
        n, nfree, scans_total, real(t1-t0,dp)/rate, obj, sum(u)+sum(v), minslack, maxms
    write(*,'(a,*(g0,:,","))') 'MPM,sum_log10_abs_matched: ', &
        sum([(log10(abs(a%val(find(j)))), j=1,n)])
    open(newunit=uu, file=trim(prefix)//'-fpref.bin', access='stream', status='replace')
    write(uu) int(cmatch, int32)
    close(uu)
    open(newunit=uu, file=trim(prefix)//'-fdual.bin', access='stream', status='replace')
    write(uu) u, v
    close(uu)
contains
    integer function find(j) result(q)
        integer, intent(in) :: j
        do q = a%col_ptr(j), a%col_ptr(j+1)-1
            if (a%row_idx(q) == cmatch(j)) return
        end do
        error stop 'matched entry missing'
    end function
    subroutine scan_column(jc, base)
        integer, intent(in) :: jc
        real(dp), intent(in) :: base
        integer :: q, r
        do q = a%col_ptr(jc), a%col_ptr(jc+1)-1
            if (cost(q) == big) cycle
            r = a%row_idx(q)
            if (state(r) == 2) cycle
            nd = base + (cost(q) - u(r) - v(jc))
            relax = relax + 1
            if (nd >= dmax) cycle
            if (state(r) == 0) then
                state(r) = 1; nt = nt + 1; touched(nt) = r
                dist(r) = nd; pred(r) = jc; call push(r)
            else if (nd < dist(r)) then
                dist(r) = nd; pred(r) = jc; call up(hpos(r))
            end if
        end do
    end subroutine
    subroutine push(r)
        integer, intent(in) :: r
        hn = hn + 1; heap(hn) = r; hpos(r) = hn
        call up(hn)
    end subroutine
    integer function pop() result(r)
        r = heap(1); hpos(r) = 0
        heap(1) = heap(hn); hn = hn - 1
        if (hn > 0) then
            hpos(heap(1)) = 1; call down(1)
        end if
    end function
    subroutine up(k0)
        integer, intent(in) :: k0
        integer :: k, pk, t
        k = k0
        do while (k > 1)
            pk = k/2
            if (dist(heap(pk)) <= dist(heap(k))) exit
            t = heap(pk); heap(pk) = heap(k); heap(k) = t
            hpos(heap(pk)) = pk; hpos(heap(k)) = k
            k = pk
        end do
    end subroutine
    subroutine down(k0)
        integer, intent(in) :: k0
        integer :: k, c, t
        k = k0
        do
            c = 2*k
            if (c > hn) exit
            if (c < hn) then
                if (dist(heap(c+1)) < dist(heap(c))) c = c + 1
            end if
            if (dist(heap(k)) <= dist(heap(c))) exit
            t = heap(c); heap(c) = heap(k); heap(k) = t
            hpos(heap(c)) = c; hpos(heap(k)) = k
            k = c
        end do
    end subroutine
end program mpm_probe
