module fortsparse_nd_refine
    ! Separators of the nested-dissection port (fortsparse_nd): METIS 5.2.1
    ! initpart.c (InitSeparator with edge bisection), balance.c, fm.c,
    ! refine.c, separator.c (ConstructSeparator), srefine.c and sfm.c, for
    ! one vertex weight. METIS is Copyright 1995-2013 Regents of the
    ! University of Minnesota, Apache License 2.0; this is a translation into
    ! Fortran with zero-based arrays and METIS's float and double expressions.
    use, intrinsic :: iso_fortran_env, only: real32, real64
    use fortsparse_nd_util, only: rpq_t, rpq_create, rpq_reset, rpq_insert, rpq_delete, &
        rpq_update, rpq_get_top, rpq_see_top_val, nd_rand_permute, nd_rand_in_range
    use fortsparse_nd_graph, only: ik, nd_graph_t, nd_ctrl_t, graph_free_rdata, graph_free
    implicit none
    private

    public :: init_separator, refine_2way_node, compute_2way_node_params

    integer, parameter :: dp = real64
    real(real32), parameter :: HALF = 0.5_real32

contains

    !> InitSeparator with METIS_IPTYPE_EDGE.
    subroutine init_separator(ctrl, graph, niparts)
        type(nd_ctrl_t), intent(inout) :: ctrl
        type(nd_graph_t), intent(inout) :: graph
        integer(ik), intent(in) :: niparts
        ! Setup2WayBalMultipliers with tpwgts = (0.5, 0.5)
        ctrl%pijbm = graph%invtvwgt/HALF
        if (graph%nedges == 0) then
            call random_bisection(ctrl, graph, niparts)
        else
            call grow_bisection(ctrl, graph, niparts)
        end if
        call compute_2way_params(graph)
        call construct_separator(ctrl, graph)
    end subroutine init_separator

    subroutine allocate_2way(graph)
        type(nd_graph_t), intent(inout) :: graph
        integer(ik) :: n
        n = max(1, graph%nvtxs)
        call graph_free_rdata(graph)
        allocate (graph%pwgts(0:2), graph%where(0:n - 1), graph%bndptr(0:n - 1), &
            graph%bndind(0:n - 1), graph%id(0:n - 1), graph%ed(0:n - 1))
    end subroutine allocate_2way

    !> RandomBisection.
    subroutine random_bisection(ctrl, graph, niparts)
        type(nd_ctrl_t), intent(inout) :: ctrl
        type(nd_graph_t), intent(inout) :: graph
        integer(ik), intent(in) :: niparts
        integer(ik), allocatable :: bestwhere(:), perm(:)
        integer(ik) :: i, ii, nvtxs, pwgts(0:1), zeromaxpwgt, bestcut, inbfs

        nvtxs = graph%nvtxs
        call allocate_2way(graph)
        allocate (bestwhere(0:nvtxs - 1), perm(0:nvtxs - 1))
        zeromaxpwgt = int((ctrl%ubfactor*real(graph%tvwgt, real32))*HALF, ik)
        bestcut = 0
        do inbfs = 0, niparts - 1
            graph%where(0:nvtxs - 1) = 1
            if (inbfs > 0) then
                call nd_rand_permute(ctrl%rng, nvtxs, perm, nvtxs/2, 1)
                pwgts(1) = graph%tvwgt
                pwgts(0) = 0
                do ii = 0, nvtxs - 1
                    i = perm(ii)
                    if (pwgts(0) + graph%vwgt(i) < zeromaxpwgt) then
                        graph%where(i) = 0
                        pwgts(0) = pwgts(0) + graph%vwgt(i)
                        pwgts(1) = pwgts(1) - graph%vwgt(i)
                        if (pwgts(0) > zeromaxpwgt) exit
                    end if
                end do
            end if
            call compute_2way_params(graph)
            call balance_2way(ctrl, graph)
            call fm_2way_cut_refine(ctrl, graph, 4_ik)
            if (inbfs == 0 .or. bestcut > graph%mincut) then
                bestcut = graph%mincut
                bestwhere = graph%where(0:nvtxs - 1)
                if (bestcut == 0) exit
            end if
        end do
        graph%mincut = bestcut
        graph%where(0:nvtxs - 1) = bestwhere
    end subroutine random_bisection

    !> GrowBisection: breadth-first growth from random seeds.
    subroutine grow_bisection(ctrl, graph, niparts)
        type(nd_ctrl_t), intent(inout) :: ctrl
        type(nd_graph_t), intent(inout) :: graph
        integer(ik), intent(in) :: niparts
        integer(ik), allocatable :: bestwhere(:), queue(:), touched(:)
        integer(ik) :: i, j, k, nvtxs, nleft, first, last, pwgts(0:1), oneminpwgt
        integer(ik) :: onemaxpwgt, bestcut, inbfs
        logical :: drain

        nvtxs = graph%nvtxs
        call allocate_2way(graph)
        allocate (bestwhere(0:nvtxs - 1), queue(0:nvtxs - 1), touched(0:nvtxs - 1))
        onemaxpwgt = int((ctrl%ubfactor*real(graph%tvwgt, real32))*HALF, ik)
        oneminpwgt = int((1.0_dp/real(ctrl%ubfactor, dp))*real(graph%tvwgt, dp)*0.5_dp, ik)
        bestcut = 0
        do inbfs = 0, niparts - 1
            graph%where(0:nvtxs - 1) = 1
            touched = 0
            pwgts(1) = graph%tvwgt
            pwgts(0) = 0
            queue(0) = nd_rand_in_range(ctrl%rng, nvtxs)
            touched(queue(0)) = 1
            first = 0
            last = 1
            nleft = nvtxs - 1
            drain = .false.
            do
                if (first == last) then
                    ! Disconnected: restart from a random untouched vertex.
                    if (nleft == 0 .or. drain) exit
                    k = nd_rand_in_range(ctrl%rng, nleft)
                    do i = 0, nvtxs - 1
                        if (touched(i) == 0) then
                            if (k == 0) exit
                            k = k - 1
                        end if
                    end do
                    queue(0) = i
                    touched(i) = 1
                    first = 0
                    last = 1
                    nleft = nleft - 1
                end if
                i = queue(first)
                first = first + 1
                if (pwgts(0) > 0 .and. pwgts(1) - graph%vwgt(i) < oneminpwgt) then
                    drain = .true.
                    cycle
                end if
                graph%where(i) = 0
                pwgts(0) = pwgts(0) + graph%vwgt(i)
                pwgts(1) = pwgts(1) - graph%vwgt(i)
                if (pwgts(1) <= onemaxpwgt) exit
                drain = .false.
                do j = graph%xadj(i), graph%xadj(i + 1) - 1
                    k = graph%adjncy(j)
                    if (touched(k) == 0) then
                        queue(last) = k
                        last = last + 1
                        touched(k) = 1
                        nleft = nleft - 1
                    end if
                end do
            end do
            if (pwgts(1) == 0) graph%where(nd_rand_in_range(ctrl%rng, nvtxs)) = 1
            if (pwgts(0) == 0) graph%where(nd_rand_in_range(ctrl%rng, nvtxs)) = 0
            call compute_2way_params(graph)
            call balance_2way(ctrl, graph)
            call fm_2way_cut_refine(ctrl, graph, ctrl%niter)
            if (inbfs == 0 .or. bestcut > graph%mincut) then
                bestcut = graph%mincut
                bestwhere = graph%where(0:nvtxs - 1)
                if (bestcut == 0) exit
            end if
        end do
        graph%mincut = bestcut
        graph%where(0:nvtxs - 1) = bestwhere
    end subroutine grow_bisection

    !> Compute2WayPartitionParams.
    subroutine compute_2way_params(graph)
        type(nd_graph_t), intent(inout) :: graph
        integer(ik) :: i, j, nbnd, mincut, istart, iend, tid, ted, me
        graph%pwgts(0:1) = 0
        graph%bndptr(0:graph%nvtxs - 1) = -1
        do i = 0, graph%nvtxs - 1
            graph%pwgts(graph%where(i)) = graph%pwgts(graph%where(i)) + graph%vwgt(i)
        end do
        nbnd = 0
        mincut = 0
        do i = 0, graph%nvtxs - 1
            istart = graph%xadj(i)
            iend = graph%xadj(i + 1)
            me = graph%where(i)
            tid = 0
            ted = 0
            do j = istart, iend - 1
                if (me == graph%where(graph%adjncy(j))) then
                    tid = tid + graph%adjwgt(j)
                else
                    ted = ted + graph%adjwgt(j)
                end if
            end do
            graph%id(i) = tid
            graph%ed(i) = ted
            if (ted > 0 .or. istart == iend) then
                call bnd_insert(nbnd, graph%bndind, graph%bndptr, i)
                mincut = mincut + ted
            end if
        end do
        graph%mincut = mincut/2
        graph%nbnd = nbnd
    end subroutine compute_2way_params

    subroutine bnd_insert(nbnd, bndind, bndptr, vtx)
        integer(ik), intent(inout) :: nbnd, bndind(0:*), bndptr(0:*)
        integer(ik), intent(in) :: vtx
        bndind(nbnd) = vtx
        bndptr(vtx) = nbnd
        nbnd = nbnd + 1
    end subroutine bnd_insert

    subroutine bnd_delete(nbnd, bndind, bndptr, vtx)
        integer(ik), intent(inout) :: nbnd, bndind(0:*), bndptr(0:*)
        integer(ik), intent(in) :: vtx
        nbnd = nbnd - 1
        bndind(bndptr(vtx)) = bndind(nbnd)
        bndptr(bndind(nbnd)) = bndptr(vtx)
        bndptr(vtx) = -1
    end subroutine bnd_delete

    !> Balance2Way with tpwgts = (0.5, 0.5).
    subroutine balance_2way(ctrl, graph)
        type(nd_ctrl_t), intent(inout) :: ctrl
        type(nd_graph_t), intent(inout) :: graph
        real(real32) :: imbalance, cur
        integer(ik) :: j
        ! ComputeLoadImbalanceDiff
        imbalance = -1.0_real32
        do j = 0, 1
            cur = real(graph%pwgts(j), real32)*ctrl%pijbm(j) - ctrl%ubfactor
            if (cur > imbalance) imbalance = cur
        end do
        if (imbalance <= 0) return
        if (abs(HALF*real(graph%tvwgt, real32) - real(graph%pwgts(0), real32)) &
            < real(3*graph%tvwgt/graph%nvtxs, real32)) return
        call general_balance(ctrl, graph, graph%nbnd > 0)
    end subroutine balance_2way

    !> Bnd2WayBalance (boundary = .true.) and General2WayBalance.
    subroutine general_balance(ctrl, graph, boundary)
        type(nd_ctrl_t), intent(inout) :: ctrl
        type(nd_graph_t), intent(inout) :: graph
        logical, intent(in) :: boundary
        integer(ik), allocatable :: moved(:), perm(:)
        type(rpq_t) :: queue
        integer(ik) :: i, ii, j, k, kwgt, nvtxs, nbnd, nswaps, from, to, tmp, higain
        integer(ik) :: mincut, mindiff, tpwgts(0:1)

        nvtxs = graph%nvtxs
        allocate (moved(0:nvtxs - 1), perm(0:nvtxs - 1))
        tpwgts(0) = int(real(graph%tvwgt, real32)*HALF, ik)
        tpwgts(1) = graph%tvwgt - tpwgts(0)
        mindiff = abs(tpwgts(0) - graph%pwgts(0))
        from = merge(1_ik, 0_ik, graph%pwgts(0) < tpwgts(0))
        to = mod(from + 1, 2_ik)
        call rpq_create(queue, nvtxs)
        moved = -1
        if (boundary) then
            nbnd = graph%nbnd
            call nd_rand_permute(ctrl%rng, nbnd, perm, nbnd/5, 1)
            do ii = 0, nbnd - 1
                i = graph%bndind(perm(ii))
                if (graph%where(i) == from .and. graph%vwgt(i) <= mindiff) &
                    call rpq_insert(queue, i, real(graph%ed(i) - graph%id(i), real32))
            end do
        else
            call nd_rand_permute(ctrl%rng, nvtxs, perm, nvtxs/5, 1)
            do ii = 0, nvtxs - 1
                i = perm(ii)
                if (graph%where(i) == from .and. graph%vwgt(i) <= mindiff) &
                    call rpq_insert(queue, i, real(graph%ed(i) - graph%id(i), real32))
            end do
        end if
        mincut = graph%mincut
        nbnd = graph%nbnd
        do nswaps = 0, nvtxs - 1
            higain = rpq_get_top(queue)
            if (higain == -1) exit
            if (graph%pwgts(to) + graph%vwgt(higain) > tpwgts(to)) exit
            mincut = mincut - (graph%ed(higain) - graph%id(higain))
            graph%pwgts(to) = graph%pwgts(to) + graph%vwgt(higain)
            graph%pwgts(from) = graph%pwgts(from) - graph%vwgt(higain)
            graph%where(higain) = to
            moved(higain) = nswaps
            tmp = graph%id(higain)
            graph%id(higain) = graph%ed(higain)
            graph%ed(higain) = tmp
            if (boundary) then
                if (graph%ed(higain) == 0 .and. graph%xadj(higain) < graph%xadj(higain + 1)) &
                    call bnd_delete(nbnd, graph%bndind, graph%bndptr, higain)
            else
                if (graph%ed(higain) == 0 .and. graph%bndptr(higain) /= -1 .and. &
                    graph%xadj(higain) < graph%xadj(higain + 1)) &
                    call bnd_delete(nbnd, graph%bndind, graph%bndptr, higain)
                if (graph%ed(higain) > 0 .and. graph%bndptr(higain) == -1) &
                    call bnd_insert(nbnd, graph%bndind, graph%bndptr, higain)
            end if
            do j = graph%xadj(higain), graph%xadj(higain + 1) - 1
                k = graph%adjncy(j)
                kwgt = merge(graph%adjwgt(j), -graph%adjwgt(j), to == graph%where(k))
                graph%id(k) = graph%id(k) + kwgt
                graph%ed(k) = graph%ed(k) - kwgt
                if (boundary) then
                    if (graph%bndptr(k) /= -1) then
                        if (graph%ed(k) == 0) then
                            call bnd_delete(nbnd, graph%bndind, graph%bndptr, k)
                            if (moved(k) == -1 .and. graph%where(k) == from .and. &
                                graph%vwgt(k) <= mindiff) call rpq_delete(queue, k)
                        else
                            if (moved(k) == -1 .and. graph%where(k) == from .and. &
                                graph%vwgt(k) <= mindiff) &
                                call rpq_update(queue, k, real(graph%ed(k) - graph%id(k), real32))
                        end if
                    else if (graph%ed(k) > 0) then
                        call bnd_insert(nbnd, graph%bndind, graph%bndptr, k)
                        if (moved(k) == -1 .and. graph%where(k) == from .and. &
                            graph%vwgt(k) <= mindiff) &
                            call rpq_insert(queue, k, real(graph%ed(k) - graph%id(k), real32))
                    end if
                else
                    if (moved(k) == -1 .and. graph%where(k) == from .and. &
                        graph%vwgt(k) <= mindiff) &
                        call rpq_update(queue, k, real(graph%ed(k) - graph%id(k), real32))
                    if (graph%ed(k) == 0 .and. graph%bndptr(k) /= -1) then
                        call bnd_delete(nbnd, graph%bndind, graph%bndptr, k)
                    else if (graph%ed(k) > 0 .and. graph%bndptr(k) == -1) then
                        call bnd_insert(nbnd, graph%bndind, graph%bndptr, k)
                    end if
                end if
            end do
        end do
        graph%mincut = mincut
        graph%nbnd = nbnd
    end subroutine general_balance

    !> FM_2WayCutRefine.
    subroutine fm_2way_cut_refine(ctrl, graph, niter)
        type(nd_ctrl_t), intent(inout) :: ctrl
        type(nd_graph_t), intent(inout) :: graph
        integer(ik), intent(in) :: niter
        integer(ik), allocatable :: moved(:), swaps(:), perm(:)
        type(rpq_t) :: queues(0:1)
        integer(ik) :: i, ii, j, k, kwgt, nvtxs, nbnd, nswaps, from, to, pass, limit, tmp
        integer(ik) :: higain, mincut, mindiff, origdiff, initcut, newcut, mincutorder
        integer(ik) :: avgvwgt, tpwgts(0:1)

        nvtxs = graph%nvtxs
        allocate (moved(0:nvtxs - 1), swaps(0:nvtxs - 1), perm(0:nvtxs - 1))
        tpwgts(0) = int(real(graph%tvwgt, real32)*HALF, ik)
        tpwgts(1) = graph%tvwgt - tpwgts(0)
        limit = int(min(max(0.01_dp*real(nvtxs, dp), 15.0_dp), 100.0_dp), ik)
        avgvwgt = min((graph%pwgts(0) + graph%pwgts(1))/20, &
            2*(graph%pwgts(0) + graph%pwgts(1))/nvtxs)
        call rpq_create(queues(0), nvtxs)
        call rpq_create(queues(1), nvtxs)
        origdiff = abs(tpwgts(0) - graph%pwgts(0))
        moved = -1
        do pass = 0, niter - 1
            call rpq_reset(queues(0))
            call rpq_reset(queues(1))
            mincutorder = -1
            initcut = graph%mincut
            mincut = initcut
            newcut = initcut
            mindiff = abs(tpwgts(0) - graph%pwgts(0))
            nbnd = graph%nbnd
            call nd_rand_permute(ctrl%rng, nbnd, perm, nbnd, 1)
            do ii = 0, nbnd - 1
                i = graph%bndind(perm(ii))
                call rpq_insert(queues(graph%where(i)), i, real(graph%ed(i) - graph%id(i), real32))
            end do
            nswaps = 0
            do while (nswaps < nvtxs)
                from = merge(0_ik, 1_ik, tpwgts(0) - graph%pwgts(0) < tpwgts(1) - graph%pwgts(1))
                to = mod(from + 1, 2_ik)
                higain = rpq_get_top(queues(from))
                if (higain == -1) exit
                newcut = newcut - (graph%ed(higain) - graph%id(higain))
                graph%pwgts(to) = graph%pwgts(to) + graph%vwgt(higain)
                graph%pwgts(from) = graph%pwgts(from) - graph%vwgt(higain)
                if ((newcut < mincut .and. abs(tpwgts(0) - graph%pwgts(0)) <= origdiff + avgvwgt) &
                    .or. (newcut == mincut .and. abs(tpwgts(0) - graph%pwgts(0)) < mindiff)) then
                    mincut = newcut
                    mindiff = abs(tpwgts(0) - graph%pwgts(0))
                    mincutorder = nswaps
                else if (nswaps - mincutorder > limit) then
                    newcut = newcut + (graph%ed(higain) - graph%id(higain))
                    graph%pwgts(from) = graph%pwgts(from) + graph%vwgt(higain)
                    graph%pwgts(to) = graph%pwgts(to) - graph%vwgt(higain)
                    exit
                end if
                graph%where(higain) = to
                moved(higain) = nswaps
                swaps(nswaps) = higain
                tmp = graph%id(higain)
                graph%id(higain) = graph%ed(higain)
                graph%ed(higain) = tmp
                if (graph%ed(higain) == 0 .and. graph%xadj(higain) < graph%xadj(higain + 1)) &
                    call bnd_delete(nbnd, graph%bndind, graph%bndptr, higain)
                do j = graph%xadj(higain), graph%xadj(higain + 1) - 1
                    k = graph%adjncy(j)
                    kwgt = merge(graph%adjwgt(j), -graph%adjwgt(j), to == graph%where(k))
                    graph%id(k) = graph%id(k) + kwgt
                    graph%ed(k) = graph%ed(k) - kwgt
                    if (graph%bndptr(k) /= -1) then
                        if (graph%ed(k) == 0) then
                            call bnd_delete(nbnd, graph%bndind, graph%bndptr, k)
                            if (moved(k) == -1) call rpq_delete(queues(graph%where(k)), k)
                        else if (moved(k) == -1) then
                            call rpq_update(queues(graph%where(k)), k, &
                                real(graph%ed(k) - graph%id(k), real32))
                        end if
                    else if (graph%ed(k) > 0) then
                        call bnd_insert(nbnd, graph%bndind, graph%bndptr, k)
                        if (moved(k) == -1) call rpq_insert(queues(graph%where(k)), k, &
                            real(graph%ed(k) - graph%id(k), real32))
                    end if
                end do
                nswaps = nswaps + 1
            end do
            ! Roll back to the best prefix of moves.
            do i = 0, nswaps - 1
                moved(swaps(i)) = -1
            end do
            nswaps = nswaps - 1
            do while (nswaps > mincutorder)
                higain = swaps(nswaps)
                graph%where(higain) = mod(graph%where(higain) + 1, 2_ik)
                to = graph%where(higain)
                tmp = graph%id(higain)
                graph%id(higain) = graph%ed(higain)
                graph%ed(higain) = tmp
                if (graph%ed(higain) == 0 .and. graph%bndptr(higain) /= -1 .and. &
                    graph%xadj(higain) < graph%xadj(higain + 1)) then
                    call bnd_delete(nbnd, graph%bndind, graph%bndptr, higain)
                else if (graph%ed(higain) > 0 .and. graph%bndptr(higain) == -1) then
                    call bnd_insert(nbnd, graph%bndind, graph%bndptr, higain)
                end if
                graph%pwgts(to) = graph%pwgts(to) + graph%vwgt(higain)
                graph%pwgts(mod(to + 1, 2_ik)) = graph%pwgts(mod(to + 1, 2_ik)) - graph%vwgt(higain)
                do j = graph%xadj(higain), graph%xadj(higain + 1) - 1
                    k = graph%adjncy(j)
                    kwgt = merge(graph%adjwgt(j), -graph%adjwgt(j), to == graph%where(k))
                    graph%id(k) = graph%id(k) + kwgt
                    graph%ed(k) = graph%ed(k) - kwgt
                    if (graph%bndptr(k) /= -1 .and. graph%ed(k) == 0) &
                        call bnd_delete(nbnd, graph%bndind, graph%bndptr, k)
                    if (graph%bndptr(k) == -1 .and. graph%ed(k) > 0) &
                        call bnd_insert(nbnd, graph%bndind, graph%bndptr, k)
                end do
                nswaps = nswaps - 1
            end do
            graph%mincut = mincut
            graph%nbnd = nbnd
            if (mincutorder <= 0 .or. mincut == initcut) exit
        end do
    end subroutine fm_2way_cut_refine

    !> ConstructSeparator: the boundary of the bisection becomes the
    !> separator, then node refinement.
    subroutine construct_separator(ctrl, graph)
        type(nd_ctrl_t), intent(inout) :: ctrl
        type(nd_graph_t), intent(inout) :: graph
        integer(ik), allocatable :: where(:)
        integer(ik) :: i, j
        allocate (where(0:max(1, graph%nvtxs) - 1))
        where(0:graph%nvtxs - 1) = graph%where(0:graph%nvtxs - 1)
        do i = 0, graph%nbnd - 1
            j = graph%bndind(i)
            if (graph%xadj(j + 1) - graph%xadj(j) > 0) where(j) = 2
        end do
        call allocate_node(graph)
        graph%where(0:graph%nvtxs - 1) = where(0:graph%nvtxs - 1)
        call compute_2way_node_params(graph)
        call fm_2way_node_refine_2sided(ctrl, graph, 1_ik)
        call fm_2way_node_refine_1sided(ctrl, graph, 4_ik)
    end subroutine construct_separator

    !> FreeRData then Allocate2WayNodePartitionMemory.
    subroutine allocate_node(graph)
        type(nd_graph_t), intent(inout) :: graph
        integer(ik) :: n
        n = max(1, graph%nvtxs)
        call graph_free_rdata(graph)
        allocate (graph%pwgts(0:2), graph%where(0:n - 1), graph%bndptr(0:n - 1), &
            graph%bndind(0:n - 1), graph%edeg(0:1, 0:n - 1))
    end subroutine allocate_node

    !> Compute2WayNodePartitionParams.
    subroutine compute_2way_node_params(graph)
        type(nd_graph_t), intent(inout) :: graph
        integer(ik) :: i, j, nbnd, me, other
        graph%pwgts(0:2) = 0
        graph%bndptr(0:graph%nvtxs - 1) = -1
        nbnd = 0
        do i = 0, graph%nvtxs - 1
            me = graph%where(i)
            graph%pwgts(me) = graph%pwgts(me) + graph%vwgt(i)
            if (me == 2) then
                call bnd_insert(nbnd, graph%bndind, graph%bndptr, i)
                graph%edeg(:, i) = 0
                do j = graph%xadj(i), graph%xadj(i + 1) - 1
                    other = graph%where(graph%adjncy(j))
                    if (other /= 2) graph%edeg(other, i) = graph%edeg(other, i) + &
                        graph%vwgt(graph%adjncy(j))
                end do
            end if
        end do
        graph%mincut = graph%pwgts(2)
        graph%nbnd = nbnd
    end subroutine compute_2way_node_params

    !> Refine2WayNode: projects the separator of graph (coarse) level by
    !> level up to orggraph, balancing and refining one-sided at each.
    subroutine refine_2way_node(ctrl, orggraph, graph_in)
        type(nd_ctrl_t), intent(inout) :: ctrl
        type(nd_graph_t), pointer, intent(in) :: orggraph, graph_in
        type(nd_graph_t), pointer :: graph
        graph => graph_in
        if (associated(graph, orggraph)) then
            call compute_2way_node_params(graph)
            return
        end if
        do
            graph => graph%finer
            call project_2way_node(graph)
            call fm_2way_node_balance(ctrl, graph)
            call fm_2way_node_refine_1sided(ctrl, graph, ctrl%niter)
            if (associated(graph, orggraph)) exit
        end do
    end subroutine refine_2way_node

    !> Project2WayNodePartition; frees the coarser graph.
    subroutine project_2way_node(graph)
        type(nd_graph_t), intent(inout) :: graph
        integer(ik) :: i
        call allocate_node(graph)
        do i = 0, graph%nvtxs - 1
            graph%where(i) = graph%coarser%where(graph%cmap(i))
        end do
        call graph_free(graph%coarser)
        call compute_2way_node_params(graph)
    end subroutine project_2way_node

    !> FM_2WayNodeRefine2Sided.
    subroutine fm_2way_node_refine_2sided(ctrl, graph, niter)
        type(nd_ctrl_t), intent(inout) :: ctrl
        type(nd_graph_t), intent(inout) :: graph
        integer(ik), intent(in) :: niter
        type(rpq_t) :: queues(0:1)
        integer(ik), allocatable :: mptr(:), mind(:), moved(:), swaps(:)
        integer(ik) :: i, ii, j, k, jj, kk, nvtxs, nbnd, nswaps, nmind, higain, oldgain
        integer(ik) :: mincut, initcut, mincutorder, pass, to, other, limit, badmaxpwgt
        integer(ik) :: mindiff, newdiff, u(0:1), g(0:1)
        real(real32) :: mult

        nvtxs = graph%nvtxs
        call rpq_create(queues(0), nvtxs)
        call rpq_create(queues(1), nvtxs)
        allocate (moved(0:nvtxs - 1), swaps(0:nvtxs - 1), mptr(0:nvtxs), mind(0:2*nvtxs - 1))
        mult = real(0.5_dp*real(ctrl%ubfactor, dp), real32)
        badmaxpwgt = int(mult*real(graph%pwgts(0) + graph%pwgts(1) + graph%pwgts(2), real32), ik)
        do pass = 0, niter - 1
            moved = -1
            call rpq_reset(queues(0))
            call rpq_reset(queues(1))
            mincutorder = -1
            initcut = graph%mincut
            mincut = initcut
            nbnd = graph%nbnd
            call nd_rand_permute(ctrl%rng, nbnd, swaps, nbnd, 1)
            do ii = 0, nbnd - 1
                i = graph%bndind(swaps(ii))
                call rpq_insert(queues(0), i, real(graph%vwgt(i) - graph%edeg(1, i), real32))
                call rpq_insert(queues(1), i, real(graph%vwgt(i) - graph%edeg(0, i), real32))
            end do
            if (ctrl%compress /= 0) then
                limit = min(5*nbnd, 400)
            else
                limit = min(2*nbnd, 300)
            end if
            mptr(0) = 0
            nmind = 0
            mindiff = abs(graph%pwgts(0) - graph%pwgts(1))
            to = merge(0_ik, 1_ik, graph%pwgts(0) < graph%pwgts(1))
            nswaps = 0
            do while (nswaps < nvtxs)
                u(0) = rpq_see_top_val(queues(0))
                u(1) = rpq_see_top_val(queues(1))
                if (u(0) /= -1 .and. u(1) /= -1) then
                    g(0) = graph%vwgt(u(0)) - graph%edeg(1, u(0))
                    g(1) = graph%vwgt(u(1)) - graph%edeg(0, u(1))
                    if (g(0) > g(1)) then
                        to = 0
                    else if (g(0) < g(1)) then
                        to = 1
                    else
                        to = mod(pass, 2_ik)
                    end if
                    if (graph%pwgts(to) + graph%vwgt(u(to)) > badmaxpwgt) to = mod(to + 1, 2_ik)
                else if (u(0) == -1 .and. u(1) == -1) then
                    exit
                else if (u(0) /= -1 .and. graph%pwgts(0) + graph%vwgt(max(u(0), 0_ik)) <= badmaxpwgt) then
                    to = 0
                else if (u(1) /= -1 .and. graph%pwgts(1) + graph%vwgt(max(u(1), 0_ik)) <= badmaxpwgt) then
                    to = 1
                else
                    exit
                end if
                other = mod(to + 1, 2_ik)
                higain = rpq_get_top(queues(to))
                if (moved(higain) == -1) call rpq_delete(queues(other), higain)
                if (nmind + graph%xadj(higain + 1) - graph%xadj(higain) >= 2*nvtxs - 1) exit
                graph%pwgts(2) = graph%pwgts(2) - (graph%vwgt(higain) - graph%edeg(other, higain))
                newdiff = abs(graph%pwgts(to) + graph%vwgt(higain) - &
                    (graph%pwgts(other) - graph%edeg(other, higain)))
                if (graph%pwgts(2) < mincut .or. (graph%pwgts(2) == mincut .and. newdiff < mindiff)) then
                    mincut = graph%pwgts(2)
                    mincutorder = nswaps
                    mindiff = newdiff
                else if (nswaps - mincutorder > 2*limit .or. (nswaps - mincutorder > limit .and. &
                        real(graph%pwgts(2), dp) > 1.10_dp*real(mincut, dp))) then
                    graph%pwgts(2) = graph%pwgts(2) + (graph%vwgt(higain) - graph%edeg(other, higain))
                    exit
                end if
                call bnd_delete(nbnd, graph%bndind, graph%bndptr, higain)
                graph%pwgts(to) = graph%pwgts(to) + graph%vwgt(higain)
                graph%where(higain) = to
                moved(higain) = nswaps
                swaps(nswaps) = higain
                do j = graph%xadj(higain), graph%xadj(higain + 1) - 1
                    k = graph%adjncy(j)
                    if (graph%where(k) == 2) then
                        oldgain = graph%vwgt(k) - graph%edeg(to, k)
                        graph%edeg(to, k) = graph%edeg(to, k) + graph%vwgt(higain)
                        if (moved(k) == -1 .or. moved(k) == -(2 + other)) &
                            call rpq_update(queues(other), k, real(oldgain - graph%vwgt(higain), real32))
                    else if (graph%where(k) == other) then
                        call bnd_insert(nbnd, graph%bndind, graph%bndptr, k)
                        mind(nmind) = k
                        nmind = nmind + 1
                        graph%where(k) = 2
                        graph%pwgts(other) = graph%pwgts(other) - graph%vwgt(k)
                        graph%edeg(:, k) = 0
                        do jj = graph%xadj(k), graph%xadj(k + 1) - 1
                            kk = graph%adjncy(jj)
                            if (graph%where(kk) /= 2) then
                                graph%edeg(graph%where(kk), k) = graph%edeg(graph%where(kk), k) + graph%vwgt(kk)
                            else
                                oldgain = graph%vwgt(kk) - graph%edeg(other, kk)
                                graph%edeg(other, kk) = graph%edeg(other, kk) - graph%vwgt(k)
                                if (moved(kk) == -1 .or. moved(kk) == -(2 + to)) &
                                    call rpq_update(queues(to), kk, real(oldgain + graph%vwgt(k), real32))
                            end if
                        end do
                        if (moved(k) == -1) then
                            call rpq_insert(queues(to), k, real(graph%vwgt(k) - graph%edeg(other, k), real32))
                            moved(k) = -(2 + to)
                        end if
                    end if
                end do
                mptr(nswaps + 1) = nmind
                nswaps = nswaps + 1
            end do
            nswaps = nswaps - 1
            do while (nswaps > mincutorder)
                higain = swaps(nswaps)
                to = graph%where(higain)
                other = mod(to + 1, 2_ik)
                graph%pwgts(2) = graph%pwgts(2) + graph%vwgt(higain)
                graph%pwgts(to) = graph%pwgts(to) - graph%vwgt(higain)
                graph%where(higain) = 2
                call bnd_insert(nbnd, graph%bndind, graph%bndptr, higain)
                call rollback_neighbours(graph, higain, to, other, mptr, mind, nswaps, nbnd)
                nswaps = nswaps - 1
            end do
            graph%mincut = mincut
            graph%nbnd = nbnd
            if (mincutorder == -1 .or. mincut >= initcut) exit
        end do
    end subroutine fm_2way_node_refine_2sided

    !> The rollback of one node move: the edegrees of higain and its
    !> neighbours, and the vertices that the move pulled into the separator
    !> (mind(mptr(nswaps):mptr(nswaps+1)-1)) back to side other.
    subroutine rollback_neighbours(graph, higain, to, other, mptr, mind, nswaps, nbnd)
        type(nd_graph_t), intent(inout) :: graph
        integer(ik), intent(in) :: higain, to, other, mptr(0:), mind(0:), nswaps
        integer(ik), intent(inout) :: nbnd
        integer(ik) :: j, k, jj, kk
        graph%edeg(:, higain) = 0
        do j = graph%xadj(higain), graph%xadj(higain + 1) - 1
            k = graph%adjncy(j)
            if (graph%where(k) == 2) then
                graph%edeg(to, k) = graph%edeg(to, k) - graph%vwgt(higain)
            else
                graph%edeg(graph%where(k), higain) = graph%edeg(graph%where(k), higain) + graph%vwgt(k)
            end if
        end do
        do j = mptr(nswaps), mptr(nswaps + 1) - 1
            k = mind(j)
            graph%where(k) = other
            graph%pwgts(other) = graph%pwgts(other) + graph%vwgt(k)
            graph%pwgts(2) = graph%pwgts(2) - graph%vwgt(k)
            call bnd_delete(nbnd, graph%bndind, graph%bndptr, k)
            do jj = graph%xadj(k), graph%xadj(k + 1) - 1
                kk = graph%adjncy(jj)
                if (graph%where(kk) == 2) graph%edeg(other, kk) = graph%edeg(other, kk) + graph%vwgt(k)
            end do
        end do
    end subroutine rollback_neighbours

    !> rollback_neighbours on plain arrays.
    subroutine rollback_arrays(xadj, adjncy, vwgt, where, bndind, bndptr, edeg, pwgts, &
            higain, to, other, mptr, mind, nswaps, nbnd)
        integer(ik), intent(in) :: xadj(0:*), adjncy(0:*), vwgt(0:*)
        integer(ik), intent(inout) :: where(0:*), bndind(0:*), bndptr(0:*), edeg(0:1, 0:*)
        integer(ik), intent(inout) :: pwgts(0:2), nbnd
        integer(ik), intent(in) :: higain, to, other, mptr(0:), mind(0:), nswaps
        integer(ik) :: j, k, jj, kk
        edeg(:, higain) = 0
        do j = xadj(higain), xadj(higain + 1) - 1
            k = adjncy(j)
            if (where(k) == 2) then
                edeg(to, k) = edeg(to, k) - vwgt(higain)
            else
                edeg(where(k), higain) = edeg(where(k), higain) + vwgt(k)
            end if
        end do
        do j = mptr(nswaps), mptr(nswaps + 1) - 1
            k = mind(j)
            where(k) = other
            pwgts(other) = pwgts(other) + vwgt(k)
            pwgts(2) = pwgts(2) - vwgt(k)
            call bnd_delete(nbnd, bndind, bndptr, k)
            do jj = xadj(k), xadj(k + 1) - 1
                kk = adjncy(jj)
                if (where(kk) == 2) edeg(other, kk) = edeg(other, kk) + vwgt(k)
            end do
        end do
    end subroutine rollback_arrays

    !> FM_2WayNodeRefine1Sided.
    subroutine fm_2way_node_refine_1sided(ctrl, graph, niter)
        type(nd_ctrl_t), intent(inout) :: ctrl
        type(nd_graph_t), intent(inout) :: graph
        integer(ik), intent(in) :: niter
        call node_refine_1sided(ctrl, graph%nvtxs, graph%xadj, graph%adjncy, graph%vwgt, &
            graph%where, graph%bndind, graph%bndptr, graph%edeg, graph%pwgts, graph%nbnd, &
            graph%mincut, niter)
    end subroutine fm_2way_node_refine_1sided

    subroutine node_refine_1sided(ctrl, nvtxs_g, xadj, adjncy, vwgt, where, bndind, bndptr, &
            edeg, pwgts, nbnd_g, mincut_g, niter)
        type(nd_ctrl_t), intent(inout) :: ctrl
        integer(ik), intent(in) :: nvtxs_g, niter
        integer(ik), intent(in) :: xadj(0:nvtxs_g), adjncy(0:*), vwgt(0:*)
        integer(ik), intent(inout) :: where(0:*), bndind(0:*), bndptr(0:*), edeg(0:1, 0:*)
        integer(ik), intent(inout) :: pwgts(0:2), nbnd_g, mincut_g
        type(rpq_t) :: queue
        integer(ik), allocatable :: mptr(:), mind(:), swaps(:)
        integer(ik) :: i, ii, j, k, jj, kk, nvtxs, nbnd, nswaps, nmind, higain
        integer(ik) :: mincut, initcut, mincutorder, pass, to, other, limit, badmaxpwgt
        integer(ik) :: mindiff, newdiff
        real(real32) :: mult

        nvtxs = nvtxs_g
        call rpq_create(queue, nvtxs)
        allocate (swaps(0:nvtxs - 1), mptr(0:nvtxs), mind(0:2*nvtxs - 1))
        mult = real(0.5_dp*real(ctrl%ubfactor, dp), real32)
        badmaxpwgt = int(mult*real(pwgts(0) + pwgts(1) + pwgts(2), real32), ik)
        to = merge(1_ik, 0_ik, pwgts(0) < pwgts(1))
        do pass = 0, 2*niter - 1
            other = to
            to = mod(to + 1, 2_ik)
            call rpq_reset(queue)
            mincutorder = -1
            initcut = mincut_g
            mincut = initcut
            nbnd = nbnd_g
            call nd_rand_permute(ctrl%rng, nbnd, swaps, nbnd, 1)
            do ii = 0, nbnd - 1
                i = bndind(swaps(ii))
                call rpq_insert(queue, i, real(vwgt(i) - edeg(other, i), real32))
            end do
            if (ctrl%compress /= 0) then
                limit = min(5*nbnd, 500)
            else
                limit = min(3*nbnd, 300)
            end if
            mptr(0) = 0
            nmind = 0
            mindiff = abs(pwgts(0) - pwgts(1))
            nswaps = 0
            do while (nswaps < nvtxs)
                higain = rpq_get_top(queue)
                if (higain == -1) exit
                if (nmind + xadj(higain + 1) - xadj(higain) >= 2*nvtxs - 1) exit
                if (pwgts(to) + vwgt(higain) > badmaxpwgt) exit
                pwgts(2) = pwgts(2) - (vwgt(higain) - edeg(other, higain))
                newdiff = abs(pwgts(to) + vwgt(higain) - &
                    (pwgts(other) - edeg(other, higain)))
                if (pwgts(2) < mincut .or. (pwgts(2) == mincut .and. newdiff < mindiff)) then
                    mincut = pwgts(2)
                    mincutorder = nswaps
                    mindiff = newdiff
                else if (nswaps - mincutorder > 3*limit .or. (nswaps - mincutorder > limit .and. &
                        real(pwgts(2), dp) > 1.10_dp*real(mincut, dp))) then
                    pwgts(2) = pwgts(2) + (vwgt(higain) - edeg(other, higain))
                    exit
                end if
                call bnd_delete(nbnd, bndind, bndptr, higain)
                pwgts(to) = pwgts(to) + vwgt(higain)
                where(higain) = to
                swaps(nswaps) = higain
                do j = xadj(higain), xadj(higain + 1) - 1
                    k = adjncy(j)
                    if (where(k) == 2) then
                        edeg(to, k) = edeg(to, k) + vwgt(higain)
                    else if (where(k) == other) then
                        call bnd_insert(nbnd, bndind, bndptr, k)
                        mind(nmind) = k
                        nmind = nmind + 1
                        where(k) = 2
                        pwgts(other) = pwgts(other) - vwgt(k)
                        edeg(:, k) = 0
                        do jj = xadj(k), xadj(k + 1) - 1
                            kk = adjncy(jj)
                            if (where(kk) /= 2) then
                                edeg(where(kk), k) = edeg(where(kk), k) + vwgt(kk)
                            else
                                edeg(other, kk) = edeg(other, kk) - vwgt(k)
                                call rpq_update(queue, kk, real(vwgt(kk) - edeg(other, kk), real32))
                            end if
                        end do
                        call rpq_insert(queue, k, real(vwgt(k) - edeg(other, k), real32))
                    end if
                end do
                mptr(nswaps + 1) = nmind
                nswaps = nswaps + 1
            end do
            nswaps = nswaps - 1
            do while (nswaps > mincutorder)
                higain = swaps(nswaps)
                pwgts(2) = pwgts(2) + vwgt(higain)
                pwgts(to) = pwgts(to) - vwgt(higain)
                where(higain) = 2
                call bnd_insert(nbnd, bndind, bndptr, higain)
                call rollback_arrays(xadj, adjncy, vwgt, where, bndind, bndptr, edeg, pwgts, higain, to, &
                    other, mptr, mind, nswaps, nbnd)
                nswaps = nswaps - 1
            end do
            mincut_g = mincut
            nbnd_g = nbnd
            if (mod(pass, 2_ik) == 1 .and. (mincutorder == -1 .or. mincut >= initcut)) exit
        end do
    end subroutine node_refine_1sided

    !> FM_2WayNodeBalance.
    subroutine fm_2way_node_balance(ctrl, graph)
        type(nd_ctrl_t), intent(inout) :: ctrl
        type(nd_graph_t), intent(inout) :: graph
        type(rpq_t) :: queue
        integer(ik), allocatable :: perm(:), moved(:)
        integer(ik) :: i, ii, j, k, jj, kk, nvtxs, nbnd, nswaps, gain, badmaxpwgt, higain
        integer(ik) :: oldgain, to, other
        real(real32) :: mult

        nvtxs = graph%nvtxs
        mult = real(0.5_dp*real(ctrl%ubfactor, dp), real32)
        badmaxpwgt = int(mult*real(graph%pwgts(0) + graph%pwgts(1), real32), ik)
        if (max(graph%pwgts(0), graph%pwgts(1)) < badmaxpwgt) return
        if (abs(graph%pwgts(0) - graph%pwgts(1)) < 3*graph%tvwgt/nvtxs) return
        to = merge(0_ik, 1_ik, graph%pwgts(0) < graph%pwgts(1))
        other = mod(to + 1, 2_ik)
        call rpq_create(queue, nvtxs)
        allocate (perm(0:nvtxs - 1), moved(0:nvtxs - 1))
        moved = -1
        nbnd = graph%nbnd
        call nd_rand_permute(ctrl%rng, nbnd, perm, nbnd, 1)
        do ii = 0, nbnd - 1
            i = graph%bndind(perm(ii))
            call rpq_insert(queue, i, real(graph%vwgt(i) - graph%edeg(other, i), real32))
        end do
        do nswaps = 0, nvtxs - 1
            higain = rpq_get_top(queue)
            if (higain == -1) exit
            moved(higain) = 1
            gain = graph%vwgt(higain) - graph%edeg(other, higain)
            badmaxpwgt = int(mult*real(graph%pwgts(0) + graph%pwgts(1), real32), ik)
            if (graph%pwgts(to) > graph%pwgts(other)) exit
            if (gain < 0 .and. graph%pwgts(other) < badmaxpwgt) exit
            if (graph%pwgts(to) + graph%vwgt(higain) > badmaxpwgt) cycle
            graph%pwgts(2) = graph%pwgts(2) - gain
            call bnd_delete(nbnd, graph%bndind, graph%bndptr, higain)
            graph%pwgts(to) = graph%pwgts(to) + graph%vwgt(higain)
            graph%where(higain) = to
            do j = graph%xadj(higain), graph%xadj(higain + 1) - 1
                k = graph%adjncy(j)
                if (graph%where(k) == 2) then
                    graph%edeg(to, k) = graph%edeg(to, k) + graph%vwgt(higain)
                else if (graph%where(k) == other) then
                    call bnd_insert(nbnd, graph%bndind, graph%bndptr, k)
                    graph%where(k) = 2
                    graph%pwgts(other) = graph%pwgts(other) - graph%vwgt(k)
                    graph%edeg(:, k) = 0
                    do jj = graph%xadj(k), graph%xadj(k + 1) - 1
                        kk = graph%adjncy(jj)
                        if (graph%where(kk) /= 2) then
                            graph%edeg(graph%where(kk), k) = graph%edeg(graph%where(kk), k) + graph%vwgt(kk)
                        else
                            oldgain = graph%vwgt(kk) - graph%edeg(other, kk)
                            graph%edeg(other, kk) = graph%edeg(other, kk) - graph%vwgt(k)
                            if (moved(kk) == -1) &
                                call rpq_update(queue, kk, real(oldgain + graph%vwgt(k), real32))
                        end if
                    end do
                    call rpq_insert(queue, k, real(graph%vwgt(k) - graph%edeg(other, k), real32))
                end if
            end do
        end do
        graph%mincut = graph%pwgts(2)
        graph%nbnd = nbnd
    end subroutine fm_2way_node_balance
end module fortsparse_nd_refine
