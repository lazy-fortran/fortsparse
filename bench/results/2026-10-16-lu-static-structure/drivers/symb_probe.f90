! Static symmetrized-structure diagnostic for a PARDISO-style supernodal LU.
! Usage: symb_probe matrix.mtx out_prefix [weighted_preference.bin]
! For matching in {none, structural (native Hopcroft-Karp), weighted (file)}
! and ordering in {AMD, METIS ND} (production *_general wrappers on P*A, i.e.
! the pattern of P*A + (P*A)^T), compute the symbolic Cholesky structure of the
! symmetrized matched matrix, an independent ancestor-walk count check, etree
! postorder, the SPD supernode rule (supernodes.F90 build_storage replica on
! parent/count), and LU payload / flop / front statistics.
program symb_probe
    use, intrinsic :: iso_c_binding, only: c_int
    use, intrinsic :: iso_fortran_env, only: real64, int64, int32
    use fortsparse_csc, only: csc_t, csc_from_triplet
    use fortsparse_status, only: fortsparse_status_t
    use fortsparse_native_matching, only: structural_row_matching, matching_graph_rows
    use fortsparse_native_symbolic, only: native_symbolic_counts
    implicit none
    integer, parameter :: dp = real64, i8 = int64
    interface
        integer(c_int) function order_amd_general(n,nnz,ptr,row,perm,inverse) &
            bind(C,name='fsparse_native_order_amd_general')
            import :: c_int
            integer(c_int), value :: n,nnz
            integer(c_int), intent(in) :: ptr(*),row(*)
            integer(c_int), intent(out) :: perm(*),inverse(*)
        end function
        integer(c_int) function order_nd_general(n,nnz,ptr,row,perm,inverse) &
            bind(C,name='fsparse_native_order_general')
            import :: c_int
            integer(c_int), value :: n,nnz
            integer(c_int), intent(in) :: ptr(*),row(*)
            integer(c_int), intent(out) :: perm(*),inverse(*)
        end function
    end interface
    type(csc_t) :: a
    type(fortsparse_status_t) :: status
    integer :: n, m, nz, e, u, variant, ord, code, nvariants, i
    integer, allocatable :: row(:), col(:), pref(:), inv(:), work(:), graph(:)
    integer, allocatable :: perm(:), iperm(:), sp(:), si(:), cp(:), ci(:)
    real(dp), allocatable :: val(:)
    character(512) :: path, prefix, wpath
    character(16), parameter :: vname(3) = [character(16) :: 'none', 'structural', 'weighted']
    character(4), parameter :: oname(2) = ['AMD ', 'ND  ']
    integer, allocatable :: snode_starts(:)

    call get_command_argument(1, path)
    call get_command_argument(2, prefix)
    nvariants = 2
    if (command_argument_count() >= 3) then
        call get_command_argument(3, wpath)
        nvariants = 3
    end if
    open(newunit=u, file=trim(path), status='old')
    block
        character(512) :: line
        do
            read(u,'(a)') line
            if (line(1:1) /= '%') exit
        end do
        read(line,*) n, m, nz
    end block
    if (n /= m) error stop 'not square'
    allocate(row(nz), col(nz), val(nz))
    do e = 1, nz
        read(u,*) row(e), col(e), val(e)
    end do
    close(u)
    call csc_from_triplet(n, n, row, col, val, a, status)
    if (status%code /= 0) error stop 'invalid matrix'
    deallocate(row, col, val)
    write(*,'(a,i0,a,i0)') 'MATRIX n=', n, ' nnz=', a%nnz
    allocate(pref(n), inv(n), work(5*n), graph(size(a%row_idx)), perm(n), iperm(n))
    call structural_row_matching(n, a%col_ptr, a%row_idx, a%val, pref, inv, work, status)
    if (status%code /= 0) error stop 'structural matching failed'
    open(newunit=u, file=trim(prefix)//'-spref.bin', access='stream', status='replace')
    write(u) int(pref, int32)
    close(u)
    write(*,'(a)') 'HDR,variant,ordering,zero_diag,nnzL,nnzLU_sym,flopsLU_sym,check_nnzL,' // &
        'order,relaxed,nsuper,maxwidth,meanwidth,maxfront,padded_nnzLU,pad_overhead,' // &
        'payload_bytes,flopsLU_padded,frac_cols_w>=8,frac_flops_w>=8,frac_flops_w>=16,' // &
        'gemm_frac,gemm_frac_k>=8,diag_frac,trsm_frac,w1,w2-3,w4-7,w8-15,w16-31,w32-63,w64-128'
    do variant = 1, nvariants
        select case (variant)
        case (1)
            do i = 1, n
                inv(i) = i
            end do
        case (2)
            call structural_row_matching(n, a%col_ptr, a%row_idx, a%val, pref, inv, work, status)
        case (3)
            open(newunit=u, file=trim(wpath), access='stream', status='old')
            read(u) pref
            close(u)
            inv = 0
            do i = 1, n
                inv(pref(i)) = i
            end do
            if (any(inv == 0)) error stop 'weighted preference is not a permutation'
        end select
        call matching_graph_rows(a%col_ptr, a%row_idx, inv, graph)
        call symmetrize(n, a%col_ptr, graph, sp, si)
        do ord = 1, 2
            if (ord == 1) then
                code = order_amd_general(n, size(graph), a%col_ptr, graph, perm, iperm)
            else
                code = order_nd_general(n, size(graph), a%col_ptr, graph, perm, iperm)
            end if
            if (code /= 0) then
                write(*,'(a,i0)') 'ordering failed code=', code
                cycle
            end if
            call check_perm(n, perm, iperm)
            call permute_sym(n, sp, si, iperm, cp, ci)
            call analyse(trim(vname(variant)), trim(oname(ord)), n, cp, ci, &
                count_structural_zero_diag(n, a%col_ptr, graph), perm, trim(prefix))
        end do
    end do

contains

    integer function count_structural_zero_diag(n, ptr, rows) result(z)
        integer, intent(in) :: n, ptr(:), rows(:)
        integer :: j
        z = 0
        do j = 1, n
            if (.not. any(rows(ptr(j):ptr(j+1)-1) == j)) z = z + 1
        end do
    end function

    subroutine check_perm(n, perm, iperm)
        integer, intent(in) :: n, perm(:), iperm(:)
        integer :: k
        do k = 1, n
            if (perm(k) < 1 .or. perm(k) > n) error stop 'perm range'
            if (iperm(perm(k)) /= k) error stop 'perm/inverse mismatch'
        end do
    end subroutine

    subroutine symmetrize(n, ptr, rows, sp, si)
        ! Full sorted pattern of B + B^T + I.
        integer, intent(in) :: n, ptr(:), rows(:)
        integer, allocatable, intent(out) :: sp(:), si(:)
        integer, allocatable :: cnt(:), tp(:), ti(:), mark(:)
        integer :: j, p, i, k
        allocate(cnt(n), mark(n))
        cnt = 1
        do j = 1, n
            do p = ptr(j), ptr(j+1)-1
                i = rows(p)
                if (i /= j) then
                    cnt(i) = cnt(i) + 1
                    cnt(j) = cnt(j) + 1
                end if
            end do
        end do
        allocate(tp(n+1))
        tp(1) = 1
        do j = 1, n
            tp(j+1) = tp(j) + cnt(j)
        end do
        allocate(ti(tp(n+1)-1))
        cnt = 0
        do j = 1, n
            ti(tp(j)) = j
            cnt(j) = 1
        end do
        do j = 1, n
            do p = ptr(j), ptr(j+1)-1
                i = rows(p)
                if (i /= j) then
                    ti(tp(j)+cnt(j)) = i; cnt(j) = cnt(j) + 1
                    ti(tp(i)+cnt(i)) = j; cnt(i) = cnt(i) + 1
                end if
            end do
        end do
        ! Deduplicate and sort each column.
        allocate(sp(n+1))
        mark = 0
        sp(1) = 1
        k = 0
        do j = 1, n
            do p = tp(j), tp(j)+cnt(j)-1
                i = ti(p)
                if (mark(i) == j) cycle
                mark(i) = j
                k = k + 1
                ti(k) = i
            end do
            sp(j+1) = k + 1
            call isort(ti(sp(j):k))
        end do
        si = ti(1:k)
    end subroutine

    subroutine permute_sym(n, sp, si, iperm, cp, ci)
        ! C(iperm(i), iperm(j)) = S(i, j); rows sorted.
        integer, intent(in) :: n, sp(:), si(:), iperm(:)
        integer, allocatable, intent(out) :: cp(:), ci(:)
        integer :: j, p, newj, k
        allocate(cp(n+1), ci(size(si)))
        cp(1) = 1
        do j = 1, n
            cp(iperm(j)+1) = sp(j+1) - sp(j)
        end do
        do j = 1, n
            cp(j+1) = cp(j+1) + cp(j)
        end do
        do j = 1, n
            newj = iperm(j)
            k = cp(newj)
            do p = sp(j), sp(j+1)-1
                ci(k) = iperm(si(p)); k = k + 1
            end do
            call isort(ci(cp(newj):cp(newj+1)-1))
        end do
    end subroutine

    recursive subroutine isort(x)
        integer, intent(inout) :: x(:)
        integer :: i, j, t, pv, nn
        nn = size(x)
        if (nn < 24) then
            do i = 2, nn
                t = x(i); j = i - 1
                do while (j >= 1)
                    if (x(j) <= t) exit
                    x(j+1) = x(j); j = j - 1
                end do
                x(j+1) = t
            end do
            return
        end if
        pv = x(nn/2); i = 1; j = nn
        do
            do while (x(i) < pv); i = i + 1; end do
            do while (x(j) > pv); j = j - 1; end do
            if (i >= j) exit
            t = x(i); x(i) = x(j); x(j) = t
            i = i + 1; j = j - 1
        end do
        call isort(x(1:j))
        call isort(x(j+1:nn))
    end subroutine

    subroutine analyse(vn, on, n, cp, ci, zdiag, perm, prefix)
        character(*), intent(in) :: vn, on, prefix
        integer, intent(in) :: n, cp(:), ci(:), zdiag, perm(:)
        integer :: uu
        integer, allocatable :: parent(:), cnt(:), chk(:), mark(:), post(:), ipost(:)
        integer, allocatable :: pparent(:), pcnt(:)
        integer :: j, p, i, k
        integer(i8) :: nnzl, nnzchk
        real(dp) :: flops
        type(fortsparse_status_t) :: st
        allocate(parent(n), cnt(n), chk(n), mark(n), post(n), ipost(n), pparent(n), pcnt(n))
        call native_symbolic_counts(n, cp, ci, parent, cnt, st)
        if (st%code /= 0) error stop 'native_symbolic_counts failed'
        ! Independent oracle: row-subtree ancestor walk (counts excl. diagonal).
        chk = 0; mark = 0
        block
            integer, allocatable :: par2(:)
            allocate(par2(n)); par2 = 0
            do j = 1, n
                do p = cp(j), cp(j+1)-1
                    i = ci(p)
                    if (i >= j) exit
                    do while (i /= 0)
                        if (mark(i) == j) exit
                        mark(i) = j
                        chk(i) = chk(i) + 1
                        if (par2(i) == 0) par2(i) = j
                        i = par2(i)
                        if (i == j) exit
                    end do
                end do
            end do
            if (any(par2 /= parent)) error stop 'etree mismatch vs oracle'
        end block
        if (any(chk /= cnt)) error stop 'column counts mismatch vs oracle'
        nnzl = int(n,i8) + sum(int(cnt,i8))
        nnzchk = int(n,i8) + sum(int(chk,i8))
        flops = 0
        do j = 1, n
            flops = flops + real(cnt(j),dp) + 2.0_dp*real(cnt(j),dp)**2
        end do
        ! Etree postorder; children visited so the largest-count child is last.
        call postorder(n, parent, cnt, post)
        do k = 1, n
            ipost(post(k)) = k
        end do
        do k = 1, n
            j = post(k)
            pcnt(k) = cnt(j)
            if (parent(j) == 0) then
                pparent(k) = 0
            else
                pparent(k) = ipost(parent(j))
            end if
        end do
        call supernode_stats(vn, on, n, zdiag, nnzl, flops, nnzchk, 'given', parent, cnt, .false.)
        call supernode_stats(vn, on, n, zdiag, nnzl, flops, nnzchk, 'given', parent, cnt, .true.)
        call supernode_stats(vn, on, n, zdiag, nnzl, flops, nnzchk, 'postordered', pparent, pcnt, .false.)
        call supernode_stats(vn, on, n, zdiag, nnzl, flops, nnzchk, 'postordered', pparent, pcnt, .true.)
        ! Final order: position k holds original (matched-row-relabelled) column perm(post(k)).
        open(newunit=uu, file=prefix//'-'//vn//'-'//trim(on)//'-order.bin', access='stream', status='replace')
        write(uu) (int(perm(post(k)),int32), k=1,n)
        close(uu)
        open(newunit=uu, file=prefix//'-'//vn//'-'//trim(on)//'-snodes.bin', access='stream', status='replace')
        write(uu) int(snode_starts, int32)
        close(uu)
    end subroutine

    subroutine postorder(n, parent, cnt, post)
        integer, intent(in) :: n, parent(:), cnt(:)
        integer, intent(out) :: post(:)
        integer, allocatable :: head(:), next(:), stack(:), best(:)
        integer :: j, k, top, c, prev, b
        allocate(head(0:n), next(n), stack(n+1), best(0:n))
        ! Build child lists; place the max-count child at the list tail.
        head = 0; best = 0
        do j = 1, n
            if (best(parent(j)) == 0) then
                best(parent(j)) = j
            else if (cnt(j) > cnt(best(parent(j)))) then
                best(parent(j)) = j
            end if
        end do
        do j = n, 1, -1
            if (j == best(parent(j))) cycle
            next(j) = head(parent(j)); head(parent(j)) = j
        end do
        ! Append best child at tail of each list.
        do j = 0, n
            b = best(j)
            if (b == 0) cycle
            next(b) = 0
            if (head(j) == 0) then
                head(j) = b
            else
                c = head(j)
                do while (next(c) /= 0)
                    c = next(c)
                end do
                next(c) = b
            end if
        end do
        ! Iterative DFS postorder over forest roots (children of 0).
        k = 0
        c = head(0)
        do while (c /= 0)
            top = 1; stack(1) = c
            do while (top > 0)
                j = stack(top)
                if (head(j) /= 0) then
                    prev = head(j); head(j) = next(prev)
                    top = top + 1; stack(top) = prev
                else
                    k = k + 1; post(k) = j; top = top - 1
                end if
            end do
            c = next(c)
        end do
        if (k /= n) error stop 'postorder incomplete'
    end subroutine

    subroutine supernode_stats(vn, on, n, zdiag, nnzl, flops, nnzchk, post, parent, cnt, relaxed)
        ! Replica of supernodes.F90 build_storage merge rule using etree
        ! parent/counts: along a parent chain the union of column patterns is
        ! {first..last} plus struct(last), so m = width + count(last).
        character(*), intent(in) :: vn, on, post
        integer, intent(in) :: n, zdiag, parent(:), cnt(:)
        integer(i8), intent(in) :: nnzl, nnzchk
        real(dp), intent(in) :: flops
        logical, intent(in) :: relaxed
        integer :: first, last, w, mm, cw, cm, exact_limit, ns, maxw, maxm, k, bin
        integer(i8) :: orig, cand, padded, padlu, rowsum, cols8
        integer(i8) :: hist(7)
        real(dp) :: fpad, f8, f16, fs, gemm, gemm8, diagf, trsm, off
        exact_limit = 64
        if (nnzl > 128_i8*n) exact_limit = 128
        ns = 0; maxw = 0; maxm = 0; padlu = 0; rowsum = 0; hist = 0; cols8 = 0
        fpad = 0; f8 = 0; f16 = 0; gemm = 0; gemm8 = 0; diagf = 0; trsm = 0
        first = 1
        do while (first <= n)
            last = first
            mm = cnt(first) + 1
            orig = mm
            do while (last < n)
                if (cnt(last) < 1) exit
                if (parent(last) /= last + 1) exit
                cw = last - first + 2
                if (cw > exact_limit) exit
                cm = cw + cnt(last+1)
                cand = orig + cnt(last+1) + 1
                padded = int(cm,i8)*cw - int(cw,i8)*(cw-1)/2
                if (padded /= cand) then
                    if (.not. relaxed) exit
                    if (cw > 32) exit
                    if (4_i8*padded > 5_i8*cand) exit
                end if
                last = last + 1
                mm = cm
                orig = cand
            end do
            w = last - first + 1
            ns = ns + 1
            if (relaxed .and. post == 'postordered') then
                if (ns == 1) then
                    if (allocated(snode_starts)) deallocate(snode_starts)
                    allocate(snode_starts(0))
                end if
                snode_starts = [snode_starts, first]
            end if
            maxw = max(maxw, w); maxm = max(maxm, mm)
            padded = int(mm,i8)*w - int(w,i8)*(w-1)/2
            padlu = padlu + 2_i8*padded - w
            rowsum = rowsum + mm
            fs = 0
            do k = 0, w-1
                off = real(mm - k - 1, dp)
                fs = fs + off + 2.0_dp*off*off
            end do
            fpad = fpad + fs
            ! Split: dense diagonal-block LU, two TRSMs, external Schur GEMM.
            diagf = diagf + 2.0_dp/3.0_dp*real(w,dp)**3
            trsm = trsm + 2.0_dp*real(mm-w,dp)*real(w,dp)**2
            gemm = gemm + 2.0_dp*real(mm-w,dp)**2*real(w,dp)
            if (w >= 8) then
                f8 = f8 + fs
                gemm8 = gemm8 + 2.0_dp*real(mm-w,dp)**2*real(w,dp)
                cols8 = cols8 + w
            end if
            if (w >= 16) f16 = f16 + fs
            bin = 1
            if (w >= 2) bin = 2
            if (w >= 4) bin = 3
            if (w >= 8) bin = 4
            if (w >= 16) bin = 5
            if (w >= 32) bin = 6
            if (w >= 64) bin = 7
            hist(bin) = hist(bin) + 1
            first = last + 1
        end do
        write(*,'(*(g0,:,","))') 'CSV', vn, on, zdiag, nnzl, 2_i8*nnzl-n, flops, nnzchk, &
            post, relaxed, ns, maxw, real(n,dp)/ns, maxm, padlu, &
            real(padlu,dp)/real(2_i8*nnzl-n,dp), 8_i8*padlu + 4_i8*rowsum, fpad, &
            real(cols8,dp)/n, f8/fpad, f16/fpad, gemm/(gemm+trsm+diagf), gemm8/(gemm+trsm+diagf), &
            diagf/(gemm+trsm+diagf), trsm/(gemm+trsm+diagf), hist
    end subroutine

    function itoa(i) result(s)
        integer, intent(in) :: i
        character(24) :: s
        write(s,'(i0)') i
    end function
    function itoa8(i) result(s)
        integer(i8), intent(in) :: i
        character(24) :: s
        write(s,'(i0)') i
    end function
end program symb_probe
