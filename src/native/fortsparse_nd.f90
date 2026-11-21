module fortsparse_nd
    ! Nested-dissection fill-reducing ordering: a translation into Fortran of
    ! METIS_NodeND of METIS 5.2.1 with its default options and seed 0
    ! (ometis.c and mmd.c, and the modules fortsparse_nd_util, _graph and
    ! _refine), which returns the orderings of METIS. METIS is Copyright
    ! 1995-2013 Regents of the University of Minnesota, Apache License 2.0.
    ! The minimum-degree leaves (genmmd) are Liu's multiple minimum degree.
    use, intrinsic :: iso_fortran_env, only: real32, real64
    use fortsparse_nd_util, only: nd_srand
    use fortsparse_nd_graph, only: ik, nd_graph_t, nd_ctrl_t, graph_setup, graph_free, &
        graph_free_rdata, compress_graph, coarsen_graph, coarsen_graph_nlevels, &
        split_graph_order
    use fortsparse_nd_refine, only: init_separator, refine_2way_node, compute_2way_node_params
    implicit none
    private

    public :: nd_order

    integer, parameter :: dp = real64
    integer(ik), parameter :: MMDSWITCH = 120, LARGENIPARTS = 7

contains

    !> perm(new) = old and iperm(old) = new (zero-based) of the graph
    !> (xadj, adjncy), zero-based, symmetric, without self edges.
    subroutine nd_order(nvtxs, xadj, adjncy, perm, iperm)
        integer(ik), intent(in) :: nvtxs, xadj(0:), adjncy(0:)
        integer(ik), intent(out) :: perm(0:), iperm(0:)
        type(nd_ctrl_t) :: ctrl
        type(nd_graph_t), pointer :: graph
        integer(ik), allocatable :: cptr(:), cind(:)
        integer(ik) :: i, ii, j, l, nnvtxs
        real(real32) :: cfactor

        if (nvtxs <= 0) return
        ctrl%ubfactor = real(real(real(1.0_dp + 0.001_dp*200, real32), dp) + 0.0000499_dp, real32)
        call nd_srand(ctrl%rng, 0_ik)
        ctrl%compress = 1
        allocate (cptr(0:nvtxs), cind(0:nvtxs - 1))
        call compress_graph(nvtxs, xadj, adjncy, cptr, cind, graph)
        if (.not. associated(graph)) then
            ctrl%compress = 0
            allocate (graph)
            call graph_setup(graph, nvtxs, xadj, adjncy)
            nnvtxs = nvtxs
        else
            nnvtxs = graph%nvtxs
            cfactor = real(1.0_dp*real(nvtxs, dp)/real(nnvtxs, dp), real32)
            if (cfactor > 1.5_real32 .and. ctrl%nseps == 1) ctrl%nseps = 2
        end if
        call mlevel_nested_dissection(ctrl, graph, iperm, nnvtxs)
        if (ctrl%compress /= 0) then
            do i = 0, nnvtxs - 1
                perm(iperm(i)) = i
            end do
            l = 0
            do ii = 0, nnvtxs - 1
                i = perm(ii)
                do j = cptr(i), cptr(i + 1) - 1
                    iperm(cind(j)) = l
                    l = l + 1
                end do
            end do
        end if
        do i = 0, nvtxs - 1
            perm(iperm(i)) = i
        end do
    end subroutine nd_order

    !> MlevelNestedDissection: order(label) of the separator last, then the
    !> two parts recursively, small ones by minimum degree. Frees graph.
    recursive subroutine mlevel_nested_dissection(ctrl, graph, order, lastvtx_in)
        type(nd_ctrl_t), intent(inout) :: ctrl
        type(nd_graph_t), pointer, intent(inout) :: graph
        integer(ik), intent(inout) :: order(0:)
        integer(ik), intent(in) :: lastvtx_in
        type(nd_graph_t), pointer :: lgraph, rgraph
        integer(ik) :: i, lastvtx, rn

        call node_bisection_multiple(ctrl, graph)
        lastvtx = lastvtx_in
        do i = 0, graph%nbnd - 1
            lastvtx = lastvtx - 1
            order(graph%label(graph%bndind(i))) = lastvtx
        end do
        call split_graph_order(graph, lgraph, rgraph)
        call graph_free(graph)
        rn = rgraph%nvtxs
        if (lgraph%nvtxs > MMDSWITCH .and. lgraph%nedges > 0) then
            call mlevel_nested_dissection(ctrl, lgraph, order, lastvtx - rn)
        else
            call mmd_order(lgraph, order, lastvtx - rn)
            call graph_free(lgraph)
        end if
        if (rgraph%nvtxs > MMDSWITCH .and. rgraph%nedges > 0) then
            call mlevel_nested_dissection(ctrl, rgraph, order, lastvtx)
        else
            call mmd_order(rgraph, order, lastvtx)
            call graph_free(rgraph)
        end if
    end subroutine mlevel_nested_dissection

    !> MlevelNodeBisectionMultiple.
    subroutine node_bisection_multiple(ctrl, graph)
        type(nd_ctrl_t), intent(inout) :: ctrl
        type(nd_graph_t), pointer, intent(inout) :: graph
        integer(ik), allocatable :: bestwhere(:)
        integer(ik) :: i, mincut

        if (ctrl%nseps == 1 .or. graph%nvtxs < merge(1000_ik, 2000_ik, ctrl%compress /= 0)) then
            call node_bisection_l2(ctrl, graph, LARGENIPARTS)
            return
        end if
        allocate (bestwhere(0:graph%nvtxs - 1))
        mincut = graph%tvwgt
        do i = 0, ctrl%nseps - 1
            call node_bisection_l2(ctrl, graph, LARGENIPARTS)
            if (i == 0 .or. graph%mincut < mincut) then
                mincut = graph%mincut
                if (i < ctrl%nseps - 1) bestwhere = graph%where(0:graph%nvtxs - 1)
            end if
            if (mincut == 0) exit
            if (i < ctrl%nseps - 1) call graph_free_rdata(graph)
        end do
        if (mincut /= graph%mincut) then
            graph%where(0:graph%nvtxs - 1) = bestwhere
            call compute_2way_node_params(graph)
        end if
    end subroutine node_bisection_multiple

    !> MlevelNodeBisectionL2: for 5000 vertices or more, five separators of
    !> a four-level coarsening, the best refined up.
    subroutine node_bisection_l2(ctrl, graph, niparts)
        type(nd_ctrl_t), intent(inout) :: ctrl
        type(nd_graph_t), pointer, intent(inout) :: graph
        integer(ik), intent(in) :: niparts
        type(nd_graph_t), pointer :: cgraph
        integer(ik), allocatable :: bestwhere(:)
        integer(ik) :: i, mincut
        integer(ik), parameter :: nruns = 5

        if (graph%nvtxs < 5000) then
            call node_bisection_l1(ctrl, graph, niparts)
            return
        end if
        ctrl%coarsen_to = max(100_ik, graph%nvtxs/30)
        cgraph => coarsen_graph_nlevels(ctrl, graph, 4_ik)
        allocate (bestwhere(0:cgraph%nvtxs - 1))
        mincut = graph%tvwgt
        do i = 0, nruns - 1
            call node_bisection_l1(ctrl, cgraph, int(0.7_dp*real(niparts, dp), ik))
            if (i == 0 .or. cgraph%mincut < mincut) then
                mincut = cgraph%mincut
                if (i < nruns - 1) bestwhere = cgraph%where(0:cgraph%nvtxs - 1)
            end if
            if (mincut == 0) exit
            if (i < nruns - 1) call graph_free_rdata(cgraph)
        end do
        if (mincut /= cgraph%mincut) cgraph%where(0:cgraph%nvtxs - 1) = bestwhere
        call refine_2way_node(ctrl, graph, cgraph)
    end subroutine node_bisection_l2

    !> MlevelNodeBisectionL1: coarsen, initial separator, refine up.
    subroutine node_bisection_l1(ctrl, graph, niparts_in)
        type(nd_ctrl_t), intent(inout) :: ctrl
        type(nd_graph_t), pointer, intent(inout) :: graph
        integer(ik), intent(in) :: niparts_in
        type(nd_graph_t), pointer :: cgraph
        integer(ik) :: niparts

        ctrl%coarsen_to = graph%nvtxs/8
        if (ctrl%coarsen_to > 100) then
            ctrl%coarsen_to = 100
        else if (ctrl%coarsen_to < 40) then
            ctrl%coarsen_to = 40
        end if
        cgraph => coarsen_graph(ctrl, graph)
        if (cgraph%nvtxs <= ctrl%coarsen_to) then
            niparts = max(1_ik, niparts_in/2)
        else
            niparts = max(1_ik, niparts_in)
        end if
        call init_separator(ctrl, cgraph, niparts)
        call refine_2way_node(ctrl, graph, cgraph)
    end subroutine node_bisection_l1

    !> MMDOrder: minimum-degree order of graph, numbered down from lastvtx.
    subroutine mmd_order(graph, order, lastvtx)
        type(nd_graph_t), intent(in) :: graph
        integer(ik), intent(inout) :: order(0:)
        integer(ik), intent(in) :: lastvtx
        integer(ik), allocatable :: xadj(:), adjncy(:), perm(:), iperm(:), head(:), qsize(:)
        integer(ik), allocatable :: list(:), marker(:)
        integer(ik) :: i, n, nofsub, firstvtx

        n = graph%nvtxs
        if (n <= 0) return
        allocate (xadj(1:n + 1), adjncy(1:max(1, graph%nedges)), perm(1:n + 5), &
            iperm(1:n + 5), head(1:n + 5), qsize(1:n + 5), list(1:n + 5), marker(1:n + 5))
        xadj(1:n + 1) = graph%xadj(0:n) + 1
        if (graph%nedges > 0) adjncy(1:graph%nedges) = graph%adjncy(0:graph%nedges - 1) + 1
        call genmmd(n, xadj, adjncy, iperm, perm, 1_ik, head, qsize, list, marker, &
            huge(1_ik), nofsub)
        firstvtx = lastvtx - n
        do i = 0, n - 1
            order(graph%label(i)) = firstvtx + iperm(i + 1) - 1
        end do
    end subroutine mmd_order

    !> genmmd: multiple minimum degree (Liu), one-based; adjncy is
    !> overwritten. invp(node) is the new position of node.
    subroutine genmmd(neqns, xadj, adjncy, invp, perm, delta, head, qsize, list, marker, &
            maxint, ncsub)
        integer(ik), intent(in) :: neqns, delta, maxint
        integer(ik), intent(inout) :: xadj(:), adjncy(:), invp(:), perm(:), head(:), qsize(:)
        integer(ik), intent(inout) :: list(:), marker(:)
        integer(ik), intent(out) :: ncsub
        integer(ik) :: ehead, i, mdeg, mdlmt, mdeg_node, nextmd, num, tag
        logical :: done, update

        ncsub = 0
        if (neqns <= 0) return
        call mmdint(neqns, xadj, head, invp, perm, qsize, list, marker)
        num = 1
        nextmd = head(1)
        do while (nextmd > 0)
            mdeg_node = nextmd
            nextmd = invp(mdeg_node)
            marker(mdeg_node) = maxint
            invp(mdeg_node) = -num
            num = num + 1
        end do
        done = num > neqns
        tag = 1
        head(1) = 0
        mdeg = 2
        do while (.not. done)
            do while (head(mdeg) <= 0)
                mdeg = mdeg + 1
            end do
            mdlmt = min(neqns, mdeg + delta)
            ehead = 0
            update = .false.
            do
                mdeg_node = head(mdeg)
                do while (mdeg_node <= 0)
                    mdeg = mdeg + 1
                    if (mdeg > mdlmt) then
                        update = .true.
                        exit
                    end if
                    mdeg_node = head(mdeg)
                end do
                if (update) exit
                nextmd = invp(mdeg_node)
                head(mdeg) = nextmd
                if (nextmd > 0) perm(nextmd) = -mdeg
                invp(mdeg_node) = -num
                ncsub = ncsub + mdeg + qsize(mdeg_node) - 2
                if (num + qsize(mdeg_node) > neqns) then
                    done = .true.
                    exit
                end if
                tag = tag + 1
                if (tag >= maxint) then
                    tag = 1
                    do i = 1, neqns
                        if (marker(i) < maxint) marker(i) = 0
                    end do
                end if
                call mmdelm(mdeg_node, xadj, adjncy, head, invp, perm, qsize, list, marker, &
                    maxint, tag)
                num = num + qsize(mdeg_node)
                list(mdeg_node) = ehead
                ehead = mdeg_node
                if (delta < 0) then
                    update = .true.
                    exit
                end if
            end do
            if (done) exit
            if (num > neqns) exit
            call mmdupd(ehead, neqns, xadj, adjncy, delta, mdeg, head, invp, perm, qsize, &
                list, marker, maxint, tag)
        end do
        call mmdnum(neqns, perm, invp, qsize)
    end subroutine genmmd

    !> mmdint: degree lists.
    subroutine mmdint(neqns, xadj, head, forward, backward, qsize, list, marker)
        integer(ik), intent(in) :: neqns, xadj(:)
        integer(ik), intent(inout) :: head(:), forward(:), backward(:), qsize(:), list(:)
        integer(ik), intent(inout) :: marker(:)
        integer(ik) :: fnode, ndeg, node
        do node = 1, neqns
            head(node) = 0
            qsize(node) = 1
            marker(node) = 0
            list(node) = 0
        end do
        do node = 1, neqns
            ndeg = xadj(node + 1) - xadj(node) + 1
            fnode = head(ndeg)
            forward(node) = fnode
            head(ndeg) = node
            if (fnode > 0) backward(fnode) = node
            backward(node) = -ndeg
        end do
    end subroutine mmdint

    !> mmdelm: eliminates mdeg_node and transforms the quotient graph.
    subroutine mmdelm(mdeg_node, xadj, adjncy, head, forward, backward, qsize, list, marker, &
            maxint, tag)
        integer(ik), intent(in) :: mdeg_node, maxint, tag
        integer(ik), intent(inout) :: xadj(:), adjncy(:), head(:), forward(:), backward(:)
        integer(ik), intent(inout) :: qsize(:), list(:), marker(:)
        integer(ik) :: element, i, istop, istart, j, jstop, jstart, link, nabor, node, npv
        integer(ik) :: nqnbrs, nxnode, pvnode, rlmt, rloc, rnode, xqnbr
        logical :: restart

        marker(mdeg_node) = tag
        istart = xadj(mdeg_node)
        istop = xadj(mdeg_node + 1) - 1
        element = 0
        rloc = istart
        rlmt = istop
        do i = istart, istop
            nabor = adjncy(i)
            if (nabor == 0) exit
            if (marker(nabor) < tag) then
                marker(nabor) = tag
                if (forward(nabor) < 0) then
                    list(nabor) = element
                    element = nabor
                else
                    adjncy(rloc) = nabor
                    rloc = rloc + 1
                end if
            end if
        end do
        do while (element > 0)
            adjncy(rlmt) = -element
            link = element
            do
                restart = .false.
                jstart = xadj(link)
                jstop = xadj(link + 1) - 1
                do j = jstart, jstop
                    node = adjncy(j)
                    link = -node
                    if (node < 0) then
                        restart = .true.
                        exit
                    end if
                    if (node == 0) exit
                    if (marker(node) < tag .and. forward(node) >= 0) then
                        marker(node) = tag
                        do while (rloc >= rlmt)
                            link = -adjncy(rlmt)
                            rloc = xadj(link)
                            rlmt = xadj(link + 1) - 1
                        end do
                        adjncy(rloc) = node
                        rloc = rloc + 1
                    end if
                end do
                if (.not. restart) exit
            end do
            element = list(element)
        end do
        if (rloc <= rlmt) adjncy(rloc) = 0
        link = mdeg_node
        do
            restart = .false.
            istart = xadj(link)
            istop = xadj(link + 1) - 1
            do i = istart, istop
                rnode = adjncy(i)
                link = -rnode
                if (rnode < 0) then
                    restart = .true.
                    exit
                end if
                if (rnode == 0) return
                pvnode = backward(rnode)
                if (pvnode /= 0 .and. pvnode /= -maxint) then
                    nxnode = forward(rnode)
                    if (nxnode > 0) backward(nxnode) = pvnode
                    if (pvnode > 0) forward(pvnode) = nxnode
                    npv = -pvnode
                    if (pvnode < 0) head(npv) = nxnode
                end if
                jstart = xadj(rnode)
                jstop = xadj(rnode + 1) - 1
                xqnbr = jstart
                do j = jstart, jstop
                    nabor = adjncy(j)
                    if (nabor == 0) exit
                    if (marker(nabor) < tag) then
                        adjncy(xqnbr) = nabor
                        xqnbr = xqnbr + 1
                    end if
                end do
                nqnbrs = xqnbr - jstart
                if (nqnbrs <= 0) then
                    qsize(mdeg_node) = qsize(mdeg_node) + qsize(rnode)
                    qsize(rnode) = 0
                    marker(rnode) = maxint
                    forward(rnode) = -mdeg_node
                    backward(rnode) = -maxint
                else
                    forward(rnode) = nqnbrs + 1
                    backward(rnode) = 0
                    adjncy(xqnbr) = mdeg_node
                    xqnbr = xqnbr + 1
                    if (xqnbr <= jstop) adjncy(xqnbr) = 0
                end if
            end do
            if (.not. restart) exit
        end do
    end subroutine mmdelm

    !> mmdupd: degree update of the nodes adjacent to the eliminated ones.
    subroutine mmdupd(ehead, neqns, xadj, adjncy, delta, mdeg, head, forward, backward, &
            qsize, list, marker, maxint, tag)
        integer(ik), intent(in) :: ehead, neqns, delta, maxint
        integer(ik), intent(inout) :: xadj(:), adjncy(:), mdeg, head(:), forward(:)
        integer(ik), intent(inout) :: backward(:), qsize(:), list(:), marker(:), tag
        integer(ik) :: deg, deg0, element, enode, fnode, i, iq2, istop, istart, j, jstop
        integer(ik) :: jstart, link, mdeg0, mtag, nabor, node, q2head, qxhead
        logical :: restart, ended

        mdeg0 = mdeg + delta
        element = ehead
        do while (element > 0)
            mtag = tag + mdeg0
            if (mtag >= maxint) then
                tag = 1
                do i = 1, neqns
                    if (marker(i) < maxint) marker(i) = 0
                end do
                mtag = tag + mdeg0
            end if
            q2head = 0
            qxhead = 0
            deg0 = 0
            link = element
            do
                restart = .false.
                istart = xadj(link)
                istop = xadj(link + 1) - 1
                do i = istart, istop
                    enode = adjncy(i)
                    link = -enode
                    if (enode < 0) then
                        restart = .true.
                        exit
                    end if
                    if (enode == 0) exit
                    if (qsize(enode) /= 0) then
                        deg0 = deg0 + qsize(enode)
                        marker(enode) = mtag
                        if (backward(enode) == 0) then
                            if (forward(enode) /= 2) then
                                list(enode) = qxhead
                                qxhead = enode
                            else
                                list(enode) = q2head
                                q2head = enode
                            end if
                        end if
                    end if
                end do
                if (.not. restart) exit
            end do
            ! The nodes with two neighbours (iq2 = 1), then the others.
            do iq2 = 1, 0, -1
                if (iq2 == 1) then
                    enode = q2head
                else
                    enode = qxhead
                end if
                do while (enode > 0)
                    if (backward(enode) == 0) then
                        tag = tag + 1
                        deg = deg0
                        if (iq2 == 1) then
                            istart = xadj(enode)
                            nabor = adjncy(istart)
                            if (nabor == element) nabor = adjncy(istart + 1)
                            link = nabor
                            if (forward(nabor) >= 0) then
                                deg = deg + qsize(nabor)
                            else
                                ended = .false.
                                do
                                    restart = .false.
                                    istart = xadj(link)
                                    istop = xadj(link + 1) - 1
                                    do i = istart, istop
                                        node = adjncy(i)
                                        link = -node
                                        if (node /= enode) then
                                            if (node < 0) then
                                                restart = .true.
                                                exit
                                            end if
                                            if (node == 0) then
                                                ended = .true.
                                                exit
                                            end if
                                            if (qsize(node) /= 0) then
                                                if (marker(node) < tag) then
                                                    marker(node) = tag
                                                    deg = deg + qsize(node)
                                                else if (backward(node) == 0) then
                                                    if (forward(node) == 2) then
                                                        qsize(enode) = qsize(enode) + qsize(node)
                                                        qsize(node) = 0
                                                        marker(node) = maxint
                                                        forward(node) = -enode
                                                        backward(node) = -maxint
                                                    else
                                                        backward(node) = -maxint
                                                    end if
                                                end if
                                            end if
                                        end if
                                    end do
                                    if (.not. restart .or. ended) exit
                                end do
                            end if
                        else
                            istart = xadj(enode)
                            istop = xadj(enode + 1) - 1
                            do i = istart, istop
                                nabor = adjncy(i)
                                if (nabor == 0) exit
                                if (marker(nabor) < tag) then
                                    marker(nabor) = tag
                                    link = nabor
                                    if (forward(nabor) >= 0) then
                                        deg = deg + qsize(nabor)
                                    else
                                        do
                                            restart = .false.
                                            jstart = xadj(link)
                                            jstop = xadj(link + 1) - 1
                                            do j = jstart, jstop
                                                node = adjncy(j)
                                                link = -node
                                                if (node < 0) then
                                                    restart = .true.
                                                    exit
                                                end if
                                                if (node == 0) exit
                                                if (marker(node) < tag) then
                                                    marker(node) = tag
                                                    deg = deg + qsize(node)
                                                end if
                                            end do
                                            if (.not. restart) exit
                                        end do
                                    end if
                                end if
                            end do
                        end if
                        ! n2100
                        deg = deg - qsize(enode) + 1
                        fnode = head(deg)
                        forward(enode) = fnode
                        backward(enode) = -deg
                        if (fnode > 0) backward(fnode) = enode
                        head(deg) = enode
                        if (deg < mdeg) mdeg = deg
                    end if
                    enode = list(enode)
                end do
            end do
            tag = mtag
            element = list(element)
        end do
    end subroutine mmdupd

    !> mmdnum: the final numbering from the merged-node forest.
    subroutine mmdnum(neqns, perm, invp, qsize)
        integer(ik), intent(in) :: neqns
        integer(ik), intent(inout) :: perm(:), invp(:), qsize(:)
        integer(ik) :: father, nextf, node, nqsize, num, root
        do node = 1, neqns
            nqsize = qsize(node)
            if (nqsize <= 0) perm(node) = invp(node)
            if (nqsize > 0) perm(node) = -invp(node)
        end do
        do node = 1, neqns
            if (perm(node) <= 0) then
                father = node
                do while (perm(father) <= 0)
                    father = -perm(father)
                end do
                root = father
                num = perm(root) + 1
                invp(node) = -num
                perm(root) = num
                father = node
                nextf = -perm(father)
                do while (nextf > 0)
                    perm(father) = -root
                    father = nextf
                    nextf = -perm(father)
                end do
            end if
        end do
        do node = 1, neqns
            num = -invp(node)
            invp(node) = num
            perm(num) = node
        end do
    end subroutine mmdnum
end module fortsparse_nd
