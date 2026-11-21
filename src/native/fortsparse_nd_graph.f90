module fortsparse_nd_graph
    ! Graphs of the nested-dissection port (fortsparse_nd) and their
    ! coarsening: METIS 5.2.1 graph.c, compress.c, coarsen.c, bucketsort.c and
    ! SplitGraphOrder of ometis.c, for one vertex weight. METIS is Copyright
    ! 1995-2013 Regents of the University of Minnesota, Apache License 2.0;
    ! this is a translation into Fortran. Arrays are zero-based like METIS's,
    ! and each expression keeps its C precision (real32 is METIS's real_t).
    use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
    use fortsparse_nd_util, only: nd_rng_t, nd_rand_permute, nd_sort_kv_inc
    implicit none
    private

    public :: ik, nd_graph_t, nd_ctrl_t, graph_setup, graph_free, graph_free_rdata
    public :: graph_tvwgt, compress_graph, coarsen_graph, coarsen_graph_nlevels
    public :: split_graph_order

    integer, parameter :: ik = int32
    integer, parameter :: dp = real64
    integer(ik), parameter :: UNMATCHED = -1
    integer(ik), parameter :: HTLENGTH = 8191
    real(dp), parameter :: COARSEN_FRACTION = 0.85_dp, UNMATCHEDFOR2HOP = 0.10_dp
    real(dp), parameter :: COMPRESSION_FRACTION = 0.85_dp

    type :: nd_graph_t
        integer(ik) :: nvtxs = 0, nedges = 0, tvwgt = 0, mincut = 0, nbnd = 0
        real(real32) :: invtvwgt = 0
        integer(ik), allocatable :: xadj(:), vwgt(:), adjncy(:), adjwgt(:)
        integer(ik), allocatable :: label(:), cmap(:)
        ! Partition and refinement data (FreeRData releases them).
        integer(ik), allocatable :: where(:), pwgts(:), bndptr(:), bndind(:)
        integer(ik), allocatable :: id(:), ed(:), edeg(:, :)
        type(nd_graph_t), pointer :: coarser => null(), finer => null()
    end type nd_graph_t

    type :: nd_ctrl_t
        type(nd_rng_t) :: rng
        integer(ik) :: compress = 1, nseps = 1, niter = 10, coarsen_to = 100
        integer(ik) :: maxvwgt = 0
        real(real32) :: ubfactor = 0, pijbm(0:1) = 0
    end type nd_ctrl_t

contains

    !> SetupGraph for unit vertex and edge weights, with labels.
    subroutine graph_setup(graph, nvtxs, xadj, adjncy)
        type(nd_graph_t), intent(inout) :: graph
        integer(ik), intent(in) :: nvtxs, xadj(0:), adjncy(0:)
        integer(ik) :: i
        graph%nvtxs = nvtxs
        graph%nedges = xadj(nvtxs)
        allocate (graph%xadj(0:nvtxs), graph%adjncy(0:max(0, graph%nedges - 1)), &
            graph%vwgt(0:max(0, nvtxs - 1)), graph%adjwgt(0:max(0, graph%nedges - 1)), &
            graph%label(0:max(0, nvtxs - 1)))
        graph%xadj = xadj(0:nvtxs)
        if (graph%nedges > 0) graph%adjncy = adjncy(0:graph%nedges - 1)
        graph%vwgt = 1
        graph%adjwgt = 1
        call graph_tvwgt(graph)
        do i = 0, nvtxs - 1
            graph%label(i) = i
        end do
    end subroutine graph_setup

    !> SetupGraph_tvwgt.
    subroutine graph_tvwgt(graph)
        type(nd_graph_t), intent(inout) :: graph
        integer(ik) :: t
        t = 0
        if (graph%nvtxs > 0) t = sum(graph%vwgt(0:graph%nvtxs - 1))
        graph%tvwgt = t
        graph%invtvwgt = real(1.0_dp/real(merge(t, 1_ik, t > 0), dp), real32)
    end subroutine graph_tvwgt

    !> FreeRData.
    subroutine graph_free_rdata(graph)
        type(nd_graph_t), intent(inout) :: graph
        if (allocated(graph%where)) deallocate (graph%where)
        if (allocated(graph%pwgts)) deallocate (graph%pwgts)
        if (allocated(graph%id)) deallocate (graph%id)
        if (allocated(graph%ed)) deallocate (graph%ed)
        if (allocated(graph%bndptr)) deallocate (graph%bndptr)
        if (allocated(graph%bndind)) deallocate (graph%bndind)
        if (allocated(graph%edeg)) deallocate (graph%edeg)
    end subroutine graph_free_rdata

    !> FreeGraph: releases graph and nulls the pointer.
    subroutine graph_free(graph)
        type(nd_graph_t), pointer, intent(inout) :: graph
        if (.not. associated(graph)) return
        deallocate (graph)
        graph => null()
    end subroutine graph_free

    !> CompressGraph: merges vertices with identical adjacency (and itself);
    !> returns a null graph unless the vertices shrink below 85 %.
    subroutine compress_graph(nvtxs, xadj, adjncy, cptr, cind, graph)
        integer(ik), intent(in) :: nvtxs, xadj(0:), adjncy(0:)
        integer(ik), intent(out) :: cptr(0:), cind(0:)
        type(nd_graph_t), pointer, intent(out) :: graph
        integer(ik), allocatable :: mark(:), map(:), key(:), val(:)
        integer(ik) :: i, ii, iii, j, jj, k, l, cnvtxs, cnedges
        logical :: same

        graph => null()
        allocate (mark(0:nvtxs - 1), map(0:nvtxs - 1), key(0:nvtxs - 1), val(0:nvtxs - 1))
        mark = -1
        map = -1
        do i = 0, nvtxs - 1
            k = 0
            do j = xadj(i), xadj(i + 1) - 1
                k = k + adjncy(j)
            end do
            key(i) = k + i
            val(i) = i
        end do
        call nd_sort_kv_inc(nvtxs, key, val)
        l = 0
        cptr(0) = 0
        cnvtxs = 0
        do i = 0, nvtxs - 1
            ii = val(i)
            if (map(ii) /= -1) cycle
            mark(ii) = i
            do j = xadj(ii), xadj(ii + 1) - 1
                mark(adjncy(j)) = i
            end do
            map(ii) = cnvtxs
            cind(l) = ii
            l = l + 1
            do j = i + 1, nvtxs - 1
                iii = val(j)
                if (key(i) /= key(j) .or. xadj(ii + 1) - xadj(ii) /= xadj(iii + 1) - xadj(iii)) exit
                if (map(iii) == -1) then
                    same = .true.
                    do jj = xadj(iii), xadj(iii + 1) - 1
                        if (mark(adjncy(jj)) /= i) then
                            same = .false.
                            exit
                        end if
                    end do
                    if (same) then
                        map(iii) = cnvtxs
                        cind(l) = iii
                        l = l + 1
                    end if
                end if
            end do
            cnvtxs = cnvtxs + 1
            cptr(cnvtxs) = l
        end do
        if (.not. (real(cnvtxs, dp) < COMPRESSION_FRACTION*real(nvtxs, dp))) return

        allocate (graph)
        cnedges = 0
        do i = 0, cnvtxs - 1
            ii = cind(cptr(i))
            cnedges = cnedges + xadj(ii + 1) - xadj(ii)
        end do
        allocate (graph%xadj(0:cnvtxs), graph%vwgt(0:cnvtxs - 1), &
            graph%adjncy(0:max(0, cnedges - 1)), graph%adjwgt(0:max(0, cnedges - 1)), &
            graph%label(0:cnvtxs - 1))
        graph%vwgt = 0
        graph%adjwgt = 1
        mark = -1
        l = 0
        graph%xadj(0) = 0
        do i = 0, cnvtxs - 1
            mark(i) = i
            do j = cptr(i), cptr(i + 1) - 1
                ii = cind(j)
                graph%vwgt(i) = graph%vwgt(i) + 1
                do jj = xadj(ii), xadj(ii + 1) - 1
                    k = map(adjncy(jj))
                    if (mark(k) /= i) then
                        mark(k) = i
                        graph%adjncy(l) = k
                        l = l + 1
                    end if
                end do
            end do
            graph%xadj(i + 1) = l
        end do
        graph%nvtxs = cnvtxs
        graph%nedges = l
        call graph_tvwgt(graph)
        do i = 0, cnvtxs - 1
            graph%label(i) = i
        end do
    end subroutine compress_graph

    !> CoarsenGraph: coarsens until CoarsenTo vertices or no more progress;
    !> returns the coarsest graph of the chain below graph.
    function coarsen_graph(ctrl, graph_in) result(graph)
        type(nd_ctrl_t), intent(inout) :: ctrl
        type(nd_graph_t), pointer, intent(in) :: graph_in
        type(nd_graph_t), pointer :: graph
        logical :: eqewgts

        graph => graph_in
        eqewgts = equal_edge_weights(graph)
        ctrl%maxvwgt = int(1.5_dp*real(graph%tvwgt, dp)/real(ctrl%coarsen_to, dp), ik)
        do
            if (.not. allocated(graph%cmap)) allocate (graph%cmap(0:max(0, graph%nvtxs - 1)))
            if (eqewgts .or. graph%nedges == 0) then
                call match(ctrl, graph, .false.)
            else
                call match(ctrl, graph, .true.)
            end if
            graph => graph%coarser
            eqewgts = .false.
            if (.not. (graph%nvtxs > ctrl%coarsen_to .and. &
                real(graph%nvtxs, dp) < COARSEN_FRACTION*real(graph%finer%nvtxs, dp) .and. &
                graph%nedges > graph%nvtxs/2)) exit
        end do
    end function coarsen_graph

    !> CoarsenGraphNlevels: at most nlevels coarsening steps.
    function coarsen_graph_nlevels(ctrl, graph_in, nlevels) result(graph)
        type(nd_ctrl_t), intent(inout) :: ctrl
        type(nd_graph_t), pointer, intent(in) :: graph_in
        integer(ik), intent(in) :: nlevels
        type(nd_graph_t), pointer :: graph
        logical :: eqewgts
        integer(ik) :: level

        graph => graph_in
        eqewgts = equal_edge_weights(graph)
        ctrl%maxvwgt = int(1.5_dp*real(graph%tvwgt, dp)/real(ctrl%coarsen_to, dp), ik)
        do level = 0, nlevels - 1
            if (.not. allocated(graph%cmap)) allocate (graph%cmap(0:max(0, graph%nvtxs - 1)))
            if (eqewgts .or. graph%nedges == 0) then
                call match(ctrl, graph, .false.)
            else
                call match(ctrl, graph, .true.)
            end if
            graph => graph%coarser
            eqewgts = .false.
            if (graph%nvtxs < ctrl%coarsen_to .or. &
                real(graph%nvtxs, dp) > COARSEN_FRACTION*real(graph%finer%nvtxs, dp) .or. &
                graph%nedges < graph%nvtxs/2) exit
        end do
    end function coarsen_graph_nlevels

    logical function equal_edge_weights(graph)
        type(nd_graph_t), intent(in) :: graph
        integer(ik) :: i
        equal_edge_weights = .true.
        do i = 1, graph%nedges - 1
            if (graph%adjwgt(0) /= graph%adjwgt(i)) then
                equal_edge_weights = .false.
                return
            end if
        end do
    end function equal_edge_weights

    !> Match_RM (heavy = .false.) and Match_SHEM (heavy = .true.), then
    !> CreateCoarseGraph.
    subroutine match(ctrl, graph, heavy)
        type(nd_ctrl_t), intent(inout) :: ctrl
        type(nd_graph_t), intent(inout), target :: graph
        logical, intent(in) :: heavy
        integer(ik), allocatable :: matchv(:), perm(:), tperm(:), degrees(:)
        integer(ik) :: i, nvtxs, cnvtxs, avgdegree, bnum, maxvwgt
        integer(int64) :: nunmatched

        nvtxs = graph%nvtxs
        maxvwgt = ctrl%maxvwgt
        allocate (matchv(0:nvtxs - 1), perm(0:nvtxs - 1), tperm(0:nvtxs - 1), &
            degrees(0:nvtxs - 1))
        matchv = UNMATCHED
        nunmatched = 0
        call nd_rand_permute(ctrl%rng, nvtxs, tperm, nvtxs/8, 1)
        avgdegree = int(4.0_dp*real(graph%xadj(nvtxs)/nvtxs, dp), ik)
        do i = 0, nvtxs - 1
            bnum = int(sqrt(real(1 + graph%xadj(i + 1) - graph%xadj(i), dp)), ik)
            degrees(i) = merge(avgdegree, bnum, bnum > avgdegree)
        end do
        call bucket_sort_keys_inc(nvtxs, avgdegree, degrees, tperm, perm)

        call greedy_match(nvtxs, graph%xadj, graph%adjncy, graph%adjwgt, graph%vwgt, perm, &
            matchv, graph%cmap, maxvwgt, heavy, cnvtxs, nunmatched)

        if (real(nunmatched, dp) > UNMATCHEDFOR2HOP*real(nvtxs, dp)) &
            call match_2hop(graph, perm, matchv, cnvtxs, nunmatched)

        ! Unmatched vertices match themselves; coarse vertices in fine order.
        cnvtxs = 0
        do i = 0, nvtxs - 1
            if (matchv(i) == UNMATCHED) then
                matchv(i) = i
                graph%cmap(i) = cnvtxs
                cnvtxs = cnvtxs + 1
            else if (i <= matchv(i)) then
                graph%cmap(i) = cnvtxs
                graph%cmap(matchv(i)) = cnvtxs
                cnvtxs = cnvtxs + 1
            end if
        end do
        call create_coarse_graph(graph, cnvtxs, matchv)
    end subroutine match

    !> The greedy pass of Match_RM / Match_SHEM on plain arrays.
    subroutine greedy_match(nvtxs, xadj, adjncy, adjwgt, vwgt, perm, matchv, cmap, maxvwgt, &
            heavy, cnvtxs, nunmatched)
        integer(ik), intent(in) :: nvtxs, maxvwgt
        integer(ik), intent(in) :: xadj(0:nvtxs), adjncy(0:*), adjwgt(0:*), vwgt(0:*), perm(0:*)
        integer(ik), intent(inout) :: matchv(0:*), cmap(0:*)
        logical, intent(in) :: heavy
        integer(ik), intent(out) :: cnvtxs
        integer(int64), intent(inout) :: nunmatched
        integer(ik) :: i, pi, j, k, maxidx, maxwgt, last_unmatched

        cnvtxs = 0
        last_unmatched = 0
        do pi = 0, nvtxs - 1
            i = perm(pi)
            if (matchv(i) /= UNMATCHED) cycle
            maxidx = i
            maxwgt = -1
            if (vwgt(i) < maxvwgt) then
                if (xadj(i) == xadj(i + 1)) then
                    ! An island matches the next unmatched vertex, ignoring maxvwgt.
                    last_unmatched = max(pi, last_unmatched) + 1
                    do while (last_unmatched < nvtxs)
                        j = perm(last_unmatched)
                        if (matchv(j) == UNMATCHED) then
                            maxidx = j
                            exit
                        end if
                        last_unmatched = last_unmatched + 1
                    end do
                else
                    do j = xadj(i), xadj(i + 1) - 1
                        k = adjncy(j)
                        if (heavy) then
                            if (matchv(k) == UNMATCHED .and. maxwgt < adjwgt(j) .and. &
                                vwgt(i) + vwgt(k) <= maxvwgt) then
                                maxidx = k
                                maxwgt = adjwgt(j)
                            end if
                        else
                            if (matchv(k) == UNMATCHED .and. &
                                vwgt(i) + vwgt(k) <= maxvwgt) then
                                maxidx = k
                                exit
                            end if
                        end if
                    end do
                    if (maxidx == i .and. 2*vwgt(i) < maxvwgt) then
                        nunmatched = nunmatched + 1
                        maxidx = UNMATCHED
                    end if
                end if
            end if
            if (maxidx /= UNMATCHED) then
                cmap(i) = cnvtxs
                cmap(maxidx) = cnvtxs
                cnvtxs = cnvtxs + 1
                matchv(i) = maxidx
                matchv(maxidx) = i
            end if
        end do

    end subroutine greedy_match

    !> BucketSortKeysInc.
    subroutine bucket_sort_keys_inc(n, max, keys, tperm, perm)
        integer(ik), intent(in) :: n, max, keys(0:), tperm(0:)
        integer(ik), intent(out) :: perm(0:)
        integer(ik), allocatable :: counts(:)
        integer(ik) :: i, ii
        allocate (counts(0:max + 1))
        counts = 0
        do i = 0, n - 1
            counts(keys(i)) = counts(keys(i)) + 1
        end do
        ! MAKECSR(i, max+1, counts)
        do i = 1, max
            counts(i) = counts(i) + counts(i - 1)
        end do
        do i = max + 1, 1, -1
            counts(i) = counts(i - 1)
        end do
        counts(0) = 0
        do ii = 0, n - 1
            i = tperm(ii)
            perm(counts(keys(i))) = i
            counts(keys(i)) = counts(keys(i)) + 1
        end do
    end subroutine bucket_sort_keys_inc

    !> Match_2Hop.
    subroutine match_2hop(graph, perm, matchv, cnvtxs, nunmatched)
        type(nd_graph_t), intent(inout) :: graph
        integer(ik), intent(in) :: perm(0:)
        integer(ik), intent(inout) :: matchv(0:), cnvtxs
        integer(int64), intent(inout) :: nunmatched
        call match_2hop_any(graph, perm, matchv, cnvtxs, nunmatched, 2_int64)
        call match_2hop_all(graph, perm, matchv, cnvtxs, nunmatched, 64_int64)
        if (real(nunmatched, dp) > 1.5_dp*UNMATCHEDFOR2HOP*real(graph%nvtxs, dp)) &
            call match_2hop_any(graph, perm, matchv, cnvtxs, nunmatched, 3_int64)
        if (real(nunmatched, dp) > 2.0_dp*UNMATCHEDFOR2HOP*real(graph%nvtxs, dp)) &
            call match_2hop_any(graph, perm, matchv, cnvtxs, nunmatched, &
            int(graph%nvtxs, int64))
    end subroutine match_2hop

    !> Match_2HopAny: unmatched vertices of degree below maxdegree that share
    !> a neighbour.
    subroutine match_2hop_any(graph, perm, matchv, cnvtxs, nunmatched, maxdegree)
        type(nd_graph_t), intent(inout) :: graph
        integer(ik), intent(in) :: perm(0:)
        integer(ik), intent(inout) :: matchv(0:), cnvtxs
        integer(int64), intent(inout) :: nunmatched
        integer(int64), intent(in) :: maxdegree
        integer(ik), allocatable :: colptr(:), rowind(:)
        integer(ik) :: i, pi, j, jj, nvtxs

        nvtxs = graph%nvtxs
        allocate (colptr(0:nvtxs))
        colptr = 0
        do i = 0, nvtxs - 1
            if (matchv(i) == UNMATCHED .and. &
                int(graph%xadj(i + 1) - graph%xadj(i), int64) < maxdegree) then
                do j = graph%xadj(i), graph%xadj(i + 1) - 1
                    colptr(graph%adjncy(j)) = colptr(graph%adjncy(j)) + 1
                end do
            end if
        end do
        ! MAKECSR(i, nvtxs, colptr)
        do i = 1, nvtxs - 1
            colptr(i) = colptr(i) + colptr(i - 1)
        end do
        do i = nvtxs, 1, -1
            colptr(i) = colptr(i - 1)
        end do
        colptr(0) = 0
        allocate (rowind(0:max(0, colptr(nvtxs) - 1)))
        do pi = 0, nvtxs - 1
            i = perm(pi)
            if (matchv(i) == UNMATCHED .and. &
                int(graph%xadj(i + 1) - graph%xadj(i), int64) < maxdegree) then
                do j = graph%xadj(i), graph%xadj(i + 1) - 1
                    rowind(colptr(graph%adjncy(j))) = i
                    colptr(graph%adjncy(j)) = colptr(graph%adjncy(j)) + 1
                end do
            end if
        end do
        ! SHIFTCSR(i, nvtxs, colptr)
        do i = nvtxs, 1, -1
            colptr(i) = colptr(i - 1)
        end do
        colptr(0) = 0
        do pi = 0, nvtxs - 1
            i = perm(pi)
            if (colptr(i + 1) - colptr(i) < 2) cycle
            jj = colptr(i + 1)
            j = colptr(i)
            do while (j < jj)
                if (matchv(rowind(j)) == UNMATCHED) then
                    jj = jj - 1
                    do while (jj > j)
                        if (matchv(rowind(jj)) == UNMATCHED) then
                            graph%cmap(rowind(j)) = cnvtxs
                            graph%cmap(rowind(jj)) = cnvtxs
                            cnvtxs = cnvtxs + 1
                            matchv(rowind(j)) = rowind(jj)
                            matchv(rowind(jj)) = rowind(j)
                            nunmatched = nunmatched - 2
                            exit
                        end if
                        jj = jj - 1
                    end do
                end if
                j = j + 1
            end do
        end do
    end subroutine match_2hop_any

    !> Match_2HopAll: unmatched vertices with identical adjacency lists.
    subroutine match_2hop_all(graph, perm, matchv, cnvtxs, nunmatched, maxdegree)
        type(nd_graph_t), intent(inout) :: graph
        integer(ik), intent(in) :: perm(0:)
        integer(ik), intent(inout) :: matchv(0:), cnvtxs
        integer(int64), intent(inout) :: nunmatched
        integer(int64), intent(in) :: maxdegree
        integer(ik), allocatable :: key(:), val(:), mark(:)
        integer(ik) :: i, pi, pk, j, jj, k, nvtxs, mask, idegree, ncand
        logical :: same

        nvtxs = graph%nvtxs
        mask = int(int(huge(1_ik), int64)/maxdegree, ik)
        allocate (key(0:max(0_int64, nunmatched - 1)), val(0:max(0_int64, nunmatched - 1)))
        ncand = 0
        do pi = 0, nvtxs - 1
            i = perm(pi)
            idegree = graph%xadj(i + 1) - graph%xadj(i)
            if (matchv(i) == UNMATCHED .and. idegree > 1 .and. &
                int(idegree, int64) < maxdegree) then
                k = 0
                do j = graph%xadj(i), graph%xadj(i + 1) - 1
                    k = k + mod(graph%adjncy(j), mask)
                end do
                val(ncand) = i
                key(ncand) = mod(k, mask)*int(maxdegree, ik) + idegree
                ncand = ncand + 1
            end if
        end do
        call nd_sort_kv_inc(ncand, key, val)
        allocate (mark(0:nvtxs - 1))
        mark = 0
        do pi = 0, ncand - 1
            i = val(pi)
            if (matchv(i) /= UNMATCHED) cycle
            do j = graph%xadj(i), graph%xadj(i + 1) - 1
                mark(graph%adjncy(j)) = i
            end do
            do pk = pi + 1, ncand - 1
                k = val(pk)
                if (matchv(k) /= UNMATCHED) cycle
                if (key(pi) /= key(pk)) exit
                if (graph%xadj(i + 1) - graph%xadj(i) /= graph%xadj(k + 1) - graph%xadj(k)) exit
                same = .true.
                do jj = graph%xadj(k), graph%xadj(k + 1) - 1
                    if (mark(graph%adjncy(jj)) /= i) then
                        same = .false.
                        exit
                    end if
                end do
                if (same) then
                    graph%cmap(i) = cnvtxs
                    graph%cmap(k) = cnvtxs
                    cnvtxs = cnvtxs + 1
                    matchv(i) = k
                    matchv(k) = i
                    nunmatched = nunmatched - 2
                    exit
                end if
            end do
        end do
    end subroutine match_2hop_all

    !> CreateCoarseGraph (no edge dropping, no vsize) and SetupCoarseGraph:
    !> graph%coarser becomes the contracted graph.
    subroutine create_coarse_graph(graph, cnvtxs, matchv)
        type(nd_graph_t), intent(inout), target :: graph
        integer(ik), intent(in) :: cnvtxs, matchv(0:)
        type(nd_graph_t), pointer :: cgraph

        allocate (cgraph)
        cgraph%nvtxs = cnvtxs
        cgraph%finer => graph
        graph%coarser => cgraph
        allocate (cgraph%xadj(0:cnvtxs), cgraph%adjncy(0:graph%nedges), &
            cgraph%adjwgt(0:graph%nedges), cgraph%vwgt(0:max(0, cnvtxs - 1)))
        call contract(graph%nvtxs, graph%xadj, graph%adjncy, graph%adjwgt, graph%vwgt, &
            graph%cmap, matchv, cnvtxs, cgraph%xadj, cgraph%adjncy, cgraph%adjwgt, &
            cgraph%vwgt, cgraph%nedges)
        call graph_tvwgt(cgraph)
    end subroutine create_coarse_graph

    !> The contraction of CreateCoarseGraph on plain arrays: a hash table for
    !> short merged lists, a direct table otherwise.
    subroutine contract(nvtxs, xadj, adjncy, adjwgt, vwgt, cmap, matchv, ncoarse, cxadj, &
            cadjncy, cadjwgt, cvwgt, cnedges)
        integer(ik), intent(in) :: nvtxs, ncoarse
        integer(ik), intent(in) :: xadj(0:nvtxs), adjncy(0:*), adjwgt(0:*), vwgt(0:*)
        integer(ik), intent(in) :: cmap(0:*), matchv(0:*)
        integer(ik), intent(out) :: cxadj(0:ncoarse), cadjncy(0:*), cadjwgt(0:*), cvwgt(0:*)
        integer(ik), intent(out) :: cnedges
        integer(ik), parameter :: mask = HTLENGTH
        integer(ik) :: htable(0:mask)
        integer(ik), allocatable :: dtable(:)
        integer(ik) :: v, u, w, j, jw, k, kk, m, nedges, cnvtxs, off, pass

        allocate (dtable(0:max(0, ncoarse - 1)))
        htable = -1
        dtable = -1
        cxadj(0) = 0
        cnvtxs = 0
        cnedges = 0
        do v = 0, nvtxs - 1
            u = matchv(v)
            if (u < v) cycle
            cvwgt(cnvtxs) = vwgt(v)
            if (v /= u) cvwgt(cnvtxs) = cvwgt(cnvtxs) + vwgt(u)
            off = cnedges
            if (xadj(v + 1) - xadj(v) + xadj(u + 1) - xadj(u) < shiftr(mask, 2)) then
                ! The coarse vertex itself sits first, so that it drops out.
                htable(iand(cnvtxs, mask)) = 0
                cadjncy(off) = cnvtxs
                nedges = 1
                do pass = 0, 1
                    if (pass == 0) then
                        w = v
                    else if (v /= u) then
                        w = u
                    else
                        exit
                    end if
                    do jw = xadj(w), xadj(w + 1) - 1
                        k = cmap(adjncy(jw))
                        kk = iand(k, mask)
                        do while (htable(kk) /= -1)
                            if (cadjncy(off + htable(kk)) == k) exit
                            kk = iand(kk + 1, mask)
                        end do
                        m = htable(kk)
                        if (m == -1) then
                            cadjncy(off + nedges) = k
                            cadjwgt(off + nedges) = adjwgt(jw)
                            htable(kk) = nedges
                            nedges = nedges + 1
                        else
                            cadjwgt(off + m) = cadjwgt(off + m) + adjwgt(jw)
                        end if
                    end do
                end do
                ! Reset in reverse order (LIFO), so that no removal breaks a
                ! probe chain still to be walked.
                do j = nedges - 1, 0, -1
                    k = cadjncy(off + j)
                    kk = iand(k, mask)
                    do while (cadjncy(off + htable(kk)) /= k)
                        kk = iand(kk + 1, mask)
                    end do
                    htable(kk) = -1
                end do
                nedges = nedges - 1
                cadjncy(off) = cadjncy(off + nedges)
                cadjwgt(off) = cadjwgt(off + nedges)
            else
                nedges = 0
                do pass = 0, 1
                    if (pass == 0) then
                        w = v
                    else if (v /= u) then
                        w = u
                    else
                        exit
                    end if
                    do jw = xadj(w), xadj(w + 1) - 1
                        k = cmap(adjncy(jw))
                        m = dtable(k)
                        if (m == -1) then
                            cadjncy(off + nedges) = k
                            cadjwgt(off + nedges) = adjwgt(jw)
                            dtable(k) = nedges
                            nedges = nedges + 1
                        else
                            cadjwgt(off + m) = cadjwgt(off + m) + adjwgt(jw)
                        end if
                    end do
                end do
                if (v /= u) then
                    j = dtable(cnvtxs)
                    if (j /= -1) then
                        nedges = nedges - 1
                        cadjncy(off + j) = cadjncy(off + nedges)
                        cadjwgt(off + j) = cadjwgt(off + nedges)
                        dtable(cnvtxs) = -1
                    end if
                end if
                do j = 0, nedges - 1
                    dtable(cadjncy(off + j)) = -1
                end do
            end if
            cnedges = cnedges + nedges
            cnvtxs = cnvtxs + 1
            cxadj(cnvtxs) = cnedges
        end do
    end subroutine contract

    !> SplitGraphOrder: the left and right parts of a tri-section (where 0
    !> and 1; 2 is the separator) as graphs with unit edge weights.
    subroutine split_graph_order(graph, lgraph, rgraph)
        type(nd_graph_t), intent(inout) :: graph
        type(nd_graph_t), pointer, intent(out) :: lgraph, rgraph
        integer(ik), allocatable :: rename(:)
        integer(ik) :: snvtxs(0:2), snedges(0:2), i, ii, j, k, l, mypart, istart, iend
        type(nd_graph_t), pointer :: part

        allocate (rename(0:max(0, graph%nvtxs - 1)))
        snvtxs = 0
        snedges = 0
        do i = 0, graph%nvtxs - 1
            k = graph%where(i)
            rename(i) = snvtxs(k)
            snvtxs(k) = snvtxs(k) + 1
            snedges(k) = snedges(k) + graph%xadj(i + 1) - graph%xadj(i)
        end do
        allocate (lgraph, rgraph)
        call setup_split(lgraph, snvtxs(0), snedges(0))
        call setup_split(rgraph, snvtxs(1), snedges(1))
        ! Mark the neighbours of the separator as boundary vertices.
        do ii = 0, graph%nbnd - 1
            i = graph%bndind(ii)
            do j = graph%xadj(i), graph%xadj(i + 1) - 1
                graph%bndptr(graph%adjncy(j)) = 1
            end do
        end do
        snvtxs = 0
        snedges = 0
        lgraph%xadj(0) = 0
        rgraph%xadj(0) = 0
        do i = 0, graph%nvtxs - 1
            mypart = graph%where(i)
            if (mypart == 2) cycle
            if (mypart == 0) then
                part => lgraph
            else
                part => rgraph
            end if
            istart = graph%xadj(i)
            iend = graph%xadj(i + 1)
            if (graph%bndptr(i) == -1) then
                do j = istart, iend - 1
                    part%adjncy(snedges(mypart) + j - istart) = graph%adjncy(j)
                end do
                snedges(mypart) = snedges(mypart) + iend - istart
            else
                l = snedges(mypart)
                do j = istart, iend - 1
                    k = graph%adjncy(j)
                    if (graph%where(k) == mypart) then
                        part%adjncy(l) = k
                        l = l + 1
                    end if
                end do
                snedges(mypart) = l
            end if
            part%vwgt(snvtxs(mypart)) = graph%vwgt(i)
            part%label(snvtxs(mypart)) = graph%label(i)
            snvtxs(mypart) = snvtxs(mypart) + 1
            part%xadj(snvtxs(mypart)) = snedges(mypart)
        end do
        do mypart = 0, 1
            if (mypart == 0) then
                part => lgraph
            else
                part => rgraph
            end if
            do i = 0, snedges(mypart) - 1
                part%adjwgt(i) = 1
                part%adjncy(i) = rename(part%adjncy(i))
            end do
            part%nvtxs = snvtxs(mypart)
            part%nedges = snedges(mypart)
            call graph_tvwgt(part)
        end do
    contains
        subroutine setup_split(s, n, ne)
            type(nd_graph_t), intent(inout) :: s
            integer(ik), intent(in) :: n, ne
            s%nvtxs = n
            s%nedges = ne
            allocate (s%xadj(0:n), s%vwgt(0:max(0, n - 1)), s%adjncy(0:max(0, ne - 1)), &
                s%adjwgt(0:max(0, ne - 1)), s%label(0:max(0, n - 1)))
        end subroutine setup_split
    end subroutine split_graph_order
end module fortsparse_nd_graph
