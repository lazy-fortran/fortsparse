! Approximate minimum degree ordering: P. R. Amestoy, T. A. Davis and
! I. S. Duff, "An approximate minimum degree ordering algorithm", SIAM J.
! Matrix Anal. Appl. 17(4), 886-905, 1996. Ported from SuiteSparse AMD 3.3
! (amd_aat.c, amd_1.c, amd_2.c, amd_postorder.c, amd_post_tree.c) with the
! default controls (dense = 10 sqrt(n), aggressive absorption), so the
! ordering equals amd_order's; statistics and debug code are left out. Arrays
! are zero-based as in the original. Included with the index kind ik.
!
! AMD, Copyright (c) 1996-2022, Timothy A. Davis, Patrick R. Amestoy, and
! Iain S. Duff.  All Rights Reserved.
! SPDX-License-Identifier: BSD-3-Clause
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!   * Redistributions of source code must retain the above copyright notice,
!     this list of conditions and the following disclaimer.
!   * Redistributions in binary form must reproduce the above copyright
!     notice, this list of conditions and the following disclaimer in the
!     documentation and/or other materials provided with the distribution.
!   * Neither the name of the organizations to which the authors are
!     affiliated, nor the names of its contributors may be used to endorse or
!     promote products derived from this software without specific prior
!     written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
! DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
! ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    integer(ik), parameter :: EMPTY = -1_ik

contains

    elemental integer(ik) function flip(i)
        integer(ik), intent(in) :: i
        flip = -i - 2_ik
    end function flip

    ! amd_aat: len(k) = off-diagonal entries of column k of A+A^T; returns
    ! their sum. ap/ai are one-based CSC with sorted, duplicate-free columns.
    integer(ik) function amd_aat(n, ap, ai, len, tp) result(nzaat)
        integer(ik), intent(in) :: n
        integer(c_int), intent(in) :: ap(0:*), ai(0:*)
        integer(ik), intent(out) :: len(0:*), tp(0:*)
        integer(ik) :: k, p, p1, p2, i, j, pj, pj2

        len(0:n - 1) = 0
        do k = 0, n - 1
            p1 = ap(k) - 1
            p2 = ap(k + 1) - 1
            p = p1
            do while (p < p2)
                j = ai(p) - 1
                if (j < k) then
                    len(j) = len(j) + 1
                    len(k) = len(k) + 1
                    p = p + 1
                else if (j == k) then
                    p = p + 1
                    exit
                else
                    exit
                end if
                pj2 = ap(j + 1) - 1
                pj = tp(j)
                do while (pj < pj2)
                    i = ai(pj) - 1
                    if (i < k) then
                        len(i) = len(i) + 1
                        len(j) = len(j) + 1
                        pj = pj + 1
                    else if (i == k) then
                        pj = pj + 1
                        exit
                    else
                        exit
                    end if
                end do
                tp(j) = pj
            end do
            tp(k) = p
        end do
        do j = 0, n - 1
            do pj = tp(j), ap(j + 1) - 2
                i = ai(pj) - 1
                len(i) = len(i) + 1
                len(j) = len(j) + 1
            end do
        end do
        nzaat = sum(len(0:n - 1))
    end function amd_aat

    ! amd_1: the adjacency lists of A+A^T in pe/iw from the counts in len,
    ! listed in the order amd_order lists them; sp and tp are scratch.
    subroutine amd_graph(n, ap, ai, len, pe, iw, sp, tp)
        integer(ik), intent(in) :: n
        integer(c_int), intent(in) :: ap(0:*), ai(0:*)
        integer(ik), intent(in) :: len(0:*)
        integer(ik), intent(out) :: pe(0:*), sp(0:*), tp(0:*)
        integer(ik), intent(inout) :: iw(0:*)
        integer(ik) :: i, j, k, p, p1, p2, pj, pj2, pfree

        pfree = 0
        do j = 0, n - 1
            pe(j) = pfree
            sp(j) = pfree
            pfree = pfree + len(j)
        end do
        do k = 0, n - 1
            p1 = ap(k) - 1
            p2 = ap(k + 1) - 1
            p = p1
            do while (p < p2)
                j = ai(p) - 1
                if (j < k) then
                    iw(sp(j)) = k
                    sp(j) = sp(j) + 1
                    iw(sp(k)) = j
                    sp(k) = sp(k) + 1
                    p = p + 1
                else if (j == k) then
                    p = p + 1
                    exit
                else
                    exit
                end if
                pj2 = ap(j + 1) - 1
                pj = tp(j)
                do while (pj < pj2)
                    i = ai(pj) - 1
                    if (i < k) then
                        iw(sp(i)) = j
                        sp(i) = sp(i) + 1
                        iw(sp(j)) = i
                        sp(j) = sp(j) + 1
                        pj = pj + 1
                    else if (i == k) then
                        pj = pj + 1
                        exit
                    else
                        exit
                    end if
                end do
                tp(j) = pj
            end do
            tp(k) = p
        end do
        do j = 0, n - 1
            do pj = tp(j), ap(j + 1) - 2
                i = ai(pj) - 1
                iw(sp(i)) = j
                sp(i) = sp(i) + 1
                iw(sp(j)) = i
                sp(j) = sp(j) + 1
            end do
        end do
    end subroutine amd_graph

    subroutine clear_flag(wflg, wbig, w, n)
        integer(ik), intent(inout) :: wflg
        integer(ik), intent(in) :: wbig, n
        integer(ik), intent(inout) :: w(0:*)
        integer(ik) :: x

        if (wflg < 2 .or. wflg >= wbig) then
            do x = 0, n - 1
                if (w(x) /= 0) w(x) = 1
            end do
            wflg = 2
        end if
    end subroutine clear_flag

    ! amd_post_tree, iterative: postorder of the subtree at root.
    subroutine post_tree(root, k, child, sibling, order, stack)
        integer(ik), intent(in) :: root
        integer(ik), intent(inout) :: k
        integer(ik), intent(inout) :: child(0:*), order(0:*), stack(0:*)
        integer(ik), intent(in) :: sibling(0:*)
        integer(ik) :: f, head, h, i

        head = 0
        stack(0) = root
        do while (head >= 0)
            i = stack(head)
            if (child(i) /= EMPTY) then
                f = child(i)
                do while (f /= EMPTY)
                    head = head + 1
                    f = sibling(f)
                end do
                h = head
                f = child(i)
                do while (f /= EMPTY)
                    stack(h) = f
                    h = h - 1
                    f = sibling(f)
                end do
                child(i) = EMPTY
            else
                head = head - 1
                order(i) = k
                k = k + 1
            end if
        end do
    end subroutine post_tree

    ! amd_postorder: the largest front of each family is ordered last.
    subroutine postorder(nn, parent, nv, fsize, order, child, sibling, stack)
        integer(ik), intent(in) :: nn
        integer(ik), intent(in) :: parent(0:*), nv(0:*), fsize(0:*)
        integer(ik), intent(out) :: order(0:*), child(0:*), sibling(0:*)
        integer(ik), intent(inout) :: stack(0:*)
        integer(ik) :: i, j, k, par, frsize, f, fprev, maxfrsize, bigfprev, &
            bigf, fnext

        child(0:nn - 1) = EMPTY
        sibling(0:nn - 1) = EMPTY
        do j = nn - 1, 0, -1
            if (nv(j) > 0) then
                par = parent(j)
                if (par /= EMPTY) then
                    sibling(j) = child(par)
                    child(par) = j
                end if
            end if
        end do
        do i = 0, nn - 1
            if (nv(i) > 0 .and. child(i) /= EMPTY) then
                fprev = EMPTY
                maxfrsize = EMPTY
                bigfprev = EMPTY
                bigf = EMPTY
                f = child(i)
                do while (f /= EMPTY)
                    frsize = fsize(f)
                    if (frsize >= maxfrsize) then
                        maxfrsize = frsize
                        bigfprev = fprev
                        bigf = f
                    end if
                    fprev = f
                    f = sibling(f)
                end do
                fnext = sibling(bigf)
                if (fnext /= EMPTY) then
                    if (bigfprev == EMPTY) then
                        child(i) = fnext
                    else
                        sibling(bigfprev) = fnext
                    end if
                    sibling(bigf) = EMPTY
                    sibling(fprev) = bigf
                end if
            end if
        end do
        order(0:nn - 1) = EMPTY
        k = 0
        do i = 0, nn - 1
            if (parent(i) == EMPTY .and. nv(i) > 0) &
                call post_tree(i, k, child, sibling, order, stack)
        end do
    end subroutine postorder

    ! amd_2 on the graph in pe/iw/len with elbow room iw(pfree:iwlen-1),
    ! iwlen >= pfree+n; every array is destroyed. On return last(k) is the
    ! k-th pivot and next(i) the position of i. Returns 0, or nonzero for an
    ! invalid argument or a broken internal invariant.
    integer function amd_order_core(n, pe, iw, len, iwlen, pfree_in, nv, &
            next, last, head, elen, degree, w) result(code)
        integer(ik), intent(in) :: n, iwlen, pfree_in
        integer(ik), intent(inout) :: pe(0:*), iw(0:*), len(0:*)
        integer(ik), intent(out) :: nv(0:*), next(0:*), last(0:*), head(0:*), &
            elen(0:*), degree(0:*), w(0:*)
        integer(ik) :: deg, degme, dext, lemax, e, elenme, eln, i, ilast, &
            inext, j, jlast, jnext, k, knt1, knt2, knt3, lenj, ln, me, mindeg, &
            nel, nleft, nvi, nvj, nvpiv, slenme, wbig, we, wflg, wnvi, ndense, &
            dense, hash, pfree
        integer(int64) :: hsum
        integer(ik) :: p, p1, p2, p3, p4, pdst, pend, pj, pme, pme1, pme2, pn, &
            psrc
        logical :: ok

        code = 1
        if (n <= 0 .or. n >= huge(n) / 4 .or. pfree_in < 0 .or. &
            pfree_in > iwlen - n) return
        code = 2
        pfree = pfree_in
        me = EMPTY
        mindeg = 0
        nel = 0
        lemax = 0
        dense = int(10.0d0 * sqrt(real(n, kind(1.0d0))), ik)
        dense = max(16_ik, dense)
        dense = min(n, dense)
        do i = 0, n - 1
            last(i) = EMPTY
            head(i) = EMPTY
            next(i) = EMPTY
            nv(i) = 1
            w(i) = 1
            elen(i) = 0
            degree(i) = len(i)
        end do
        wbig = huge(n) - n
        wflg = 0
        call clear_flag(wflg, wbig, w, n)
        ndense = 0
        do i = 0, n - 1
            deg = degree(i)
            if (deg == 0) then
                elen(i) = flip(1_ik)
                nel = nel + 1
                pe(i) = EMPTY
                w(i) = 0
            else if (deg > dense) then
                ndense = ndense + 1
                nv(i) = 0
                elen(i) = EMPTY
                nel = nel + 1
                pe(i) = EMPTY
            else
                inext = head(deg)
                if (inext /= EMPTY) last(inext) = i
                next(i) = inext
                head(deg) = i
            end if
        end do

        do while (nel < n)
            ! Pivot of minimum approximate degree.
            do deg = mindeg, n - 1
                me = head(deg)
                if (me /= EMPTY) exit
            end do
            if (deg >= n) return
            mindeg = deg
            inext = next(me)
            if (inext /= EMPTY) last(inext) = EMPTY
            head(deg) = inext
            elenme = elen(me)
            nvpiv = nv(me)
            nel = nel + nvpiv
            nv(me) = -nvpiv
            degme = 0

            ! Construct the new element Lme.
            if (elenme == 0) then
                pme1 = pe(me)
                pme2 = pme1 - 1
                do p = pme1, pme1 + len(me) - 1
                    i = iw(p)
                    nvi = nv(i)
                    if (nvi > 0) then
                        degme = degme + nvi
                        nv(i) = -nvi
                        pme2 = pme2 + 1
                        iw(pme2) = i
                        ilast = last(i)
                        inext = next(i)
                        if (inext /= EMPTY) last(inext) = ilast
                        if (ilast /= EMPTY) then
                            next(ilast) = inext
                        else
                            head(degree(i)) = inext
                        end if
                    end if
                end do
            else
                p = pe(me)
                pme1 = pfree
                slenme = len(me) - elenme
                do knt1 = 1, elenme + 1
                    if (knt1 > elenme) then
                        e = me
                        pj = p
                        ln = slenme
                    else
                        e = iw(p)
                        p = p + 1
                        pj = pe(e)
                        ln = len(e)
                    end if
                    do knt2 = 1, ln
                        i = iw(pj)
                        pj = pj + 1
                        nvi = nv(i)
                        if (nvi > 0) then
                            if (pfree >= iwlen) then
                                ! Compress iw.
                                pe(me) = p
                                len(me) = len(me) - knt1
                                if (len(me) == 0) pe(me) = EMPTY
                                pe(e) = pj
                                len(e) = ln - knt2
                                if (len(e) == 0) pe(e) = EMPTY
                                do j = 0, n - 1
                                    pn = pe(j)
                                    if (pn >= 0) then
                                        pe(j) = iw(pn)
                                        iw(pn) = flip(j)
                                    end if
                                end do
                                psrc = 0
                                pdst = 0
                                pend = pme1 - 1
                                do while (psrc <= pend)
                                    j = flip(iw(psrc))
                                    psrc = psrc + 1
                                    if (j >= 0) then
                                        iw(pdst) = pe(j)
                                        pe(j) = pdst
                                        pdst = pdst + 1
                                        lenj = len(j)
                                        do knt3 = 0, lenj - 2
                                            iw(pdst) = iw(psrc)
                                            pdst = pdst + 1
                                            psrc = psrc + 1
                                        end do
                                    end if
                                end do
                                p1 = pdst
                                do psrc = pme1, pfree - 1
                                    iw(pdst) = iw(psrc)
                                    pdst = pdst + 1
                                end do
                                pme1 = p1
                                pfree = pdst
                                pj = pe(e)
                                p = pe(me)
                                if (pfree >= iwlen) return
                            end if
                            degme = degme + nvi
                            nv(i) = -nvi
                            iw(pfree) = i
                            pfree = pfree + 1
                            ilast = last(i)
                            inext = next(i)
                            if (inext /= EMPTY) last(inext) = ilast
                            if (ilast /= EMPTY) then
                                next(ilast) = inext
                            else
                                head(degree(i)) = inext
                            end if
                        end if
                    end do
                    if (e /= me) then
                        ! Absorb e into me.
                        pe(e) = flip(me)
                        w(e) = 0
                    end if
                end do
                pme2 = pfree - 1
            end if
            degree(me) = degme
            pe(me) = pme1
            len(me) = pme2 - pme1 + 1
            elen(me) = flip(nvpiv + degme)

            ! w(e) - wflg = |Le \ Lme| for every element e adjacent to Lme.
            call clear_flag(wflg, wbig, w, n)
            do pme = pme1, pme2
                i = iw(pme)
                eln = elen(i)
                if (eln > 0) then
                    nvi = -nv(i)
                    wnvi = wflg - nvi
                    do p = pe(i), pe(i) + eln - 1
                        e = iw(p)
                        we = w(e)
                        if (we >= wflg) then
                            we = we - nvi
                        else if (we /= 0) then
                            we = degree(e) + wnvi
                        end if
                        w(e) = we
                    end do
                end if
            end do

            ! Degree update, element absorption and supervariable hashing.
            do pme = pme1, pme2
                i = iw(pme)
                p1 = pe(i)
                p2 = p1 + elen(i) - 1
                pn = p1
                hsum = 0
                deg = 0
                do p = p1, p2
                    e = iw(p)
                    we = w(e)
                    if (we /= 0) then
                        dext = we - wflg
                        if (dext > 0) then
                            deg = deg + dext
                            iw(pn) = e
                            pn = pn + 1
                            hsum = hsum + e
                        else
                            ! Aggressive absorption.
                            pe(e) = flip(me)
                            w(e) = 0
                        end if
                    end if
                end do
                elen(i) = pn - p1 + 1
                p3 = pn
                p4 = p1 + len(i)
                do p = p2 + 1, p4 - 1
                    j = iw(p)
                    nvj = nv(j)
                    if (nvj > 0) then
                        deg = deg + nvj
                        iw(pn) = j
                        pn = pn + 1
                        hsum = hsum + j
                    end if
                end do
                if (elen(i) == 1 .and. p3 == pn) then
                    ! Mass elimination.
                    pe(i) = flip(me)
                    nvi = -nv(i)
                    degme = degme - nvi
                    nvpiv = nvpiv + nvi
                    nel = nel + nvi
                    nv(i) = 0
                    elen(i) = EMPTY
                else
                    degree(i) = min(degree(i), deg)
                    iw(pn) = iw(p3)
                    iw(p3) = iw(p1)
                    iw(p1) = me
                    len(i) = pn - p1 + 1
                    ! AMD sums the hash as an unsigned Int, which wraps.
                    if (bit_size(n) == 32) hsum = iand(hsum, int(z'FFFFFFFF', int64))
                    hash = int(modulo(hsum, int(n, int64)), ik)
                    j = head(hash)
                    if (j <= EMPTY) then
                        next(i) = flip(j)
                        head(hash) = flip(i)
                    else
                        next(i) = last(j)
                        last(j) = i
                    end if
                    last(i) = hash
                end if
            end do
            degree(me) = degme
            lemax = max(lemax, degme)
            wflg = wflg + lemax
            call clear_flag(wflg, wbig, w, n)

            ! Supervariable detection within each hash bucket.
            do pme = pme1, pme2
                i = iw(pme)
                if (nv(i) < 0) then
                    hash = last(i)
                    j = head(hash)
                    if (j == EMPTY) then
                        i = EMPTY
                    else if (j < EMPTY) then
                        i = flip(j)
                        head(hash) = EMPTY
                    else
                        i = last(j)
                        last(j) = EMPTY
                    end if
                    do while (i /= EMPTY)
                        if (next(i) == EMPTY) exit
                        ln = len(i)
                        eln = elen(i)
                        do p = pe(i) + 1, pe(i) + ln - 1
                            w(iw(p)) = wflg
                        end do
                        jlast = i
                        j = next(i)
                        do while (j /= EMPTY)
                            ok = len(j) == ln .and. elen(j) == eln
                            p = pe(j) + 1
                            do while (ok .and. p <= pe(j) + ln - 1)
                                if (w(iw(p)) /= wflg) ok = .false.
                                p = p + 1
                            end do
                            if (ok) then
                                pe(j) = flip(i)
                                nv(i) = nv(i) + nv(j)
                                nv(j) = 0
                                elen(j) = EMPTY
                                j = next(j)
                                next(jlast) = j
                            else
                                jlast = j
                                j = next(j)
                            end if
                        end do
                        wflg = wflg + 1
                        i = next(i)
                    end do
                end if
            end do

            ! Finalize degrees and restore the degree lists.
            p = pme1
            nleft = n - nel
            do pme = pme1, pme2
                i = iw(pme)
                nvi = -nv(i)
                if (nvi > 0) then
                    nv(i) = nvi
                    deg = degree(i) + degme - nvi
                    deg = min(deg, nleft - nvi)
                    inext = head(deg)
                    if (inext /= EMPTY) last(inext) = i
                    next(i) = inext
                    last(i) = EMPTY
                    head(deg) = i
                    mindeg = min(mindeg, deg)
                    degree(i) = deg
                    iw(p) = i
                    p = p + 1
                end if
            end do
            nv(me) = nvpiv
            len(me) = p - pme1
            if (len(me) == 0) then
                pe(me) = EMPTY
                w(me) = 0
            end if
            if (elenme /= 0) pfree = p
        end do

        ! Assembly tree: pe is the parent, elen the front size.
        do i = 0, n - 1
            pe(i) = flip(pe(i))
            elen(i) = flip(elen(i))
        end do
        do i = 0, n - 1
            if (nv(i) == 0) then
                j = pe(i)
                if (j == EMPTY) cycle
                do while (nv(j) == 0)
                    j = pe(j)
                end do
                e = j
                j = i
                do while (nv(j) == 0)
                    jnext = pe(j)
                    pe(j) = e
                    j = jnext
                end do
            end if
        end do
        call postorder(n, pe, nv, elen, w, head, next, last)
        head(0:n - 1) = EMPTY
        next(0:n - 1) = EMPTY
        do e = 0, n - 1
            k = w(e)
            if (k /= EMPTY) head(k) = e
        end do
        nel = 0
        do k = 0, n - 1
            e = head(k)
            if (e == EMPTY) exit
            next(e) = nel
            nel = nel + nv(e)
        end do
        do i = 0, n - 1
            if (nv(i) == 0) then
                e = pe(i)
                if (e /= EMPTY) then
                    next(i) = next(e)
                    next(e) = next(e) + 1
                else
                    next(i) = nel
                    nel = nel + 1
                end if
            end if
        end do
        do i = 0, n - 1
            k = next(i)
            if (k < 0 .or. k >= n) return
            last(k) = i
        end do
        code = 0
    end function amd_order_core
