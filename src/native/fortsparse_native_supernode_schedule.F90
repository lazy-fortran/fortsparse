module fortsparse_native_supernode_schedule
    ! Structural supernode partition, row structure and the OpenMP update
    ! schedule. Graph topology and ledger charging only: no numerical values.
    ! The SPD panels use it today; the static-structure LU reuses the same
    ! partition, maps and task groups (docs/lu-supernodal.md, section 3).
    use fortsparse_kinds, only: dp, i8
    use fortsparse_memory, only: memory_budget_t, memory_checked_add, memory_checked_mul
    use fortsparse_status, only: fortsparse_status_t, status_set, FORTSPARSE_OK, &
        FORTSPARSE_OUT_OF_MEMORY, FORTSPARSE_INDEX_OVERFLOW, &
        FORTSPARSE_INTERNAL_ERROR, FORTSPARSE_RESOURCE_LIMIT
    ! Family-neutral name for the generated subtract-into-row-map scatter.
    use fortsparse_spd_leaves, only: supernode_scatter_mapped_update => &
        spd_scatter_mapped_update
    !$ use omp_lib, only: omp_get_thread_num, omp_get_max_threads
    implicit none
    private
    public :: supernode_t, supernode_schedule_t
    public :: supernode_partition, supernode_gather_rows, supernode_sort_rows
    public :: supernode_charge, supernode_scatter_mapped_update, supernode_find_row

    ! Exact (zero-padding) merges stop at COMPACT_WIDTH columns, or at
    ! EXACT_WIDTH when the fill can amortize wide panels. Relaxed merges add
    ! explicit zeros only up to RELAXED_WIDTH columns and 25% padding.
    integer, parameter, public :: SUPERNODE_COMPACT_WIDTH = 64, &
        SUPERNODE_EXACT_WIDTH = 128, SUPERNODE_RELAXED_WIDTH = 32

    ! Structural supernode: columns first..first+width-1 share sorted rows(1:m);
    ! rows(1:width) are the supernode's own columns.
    type :: supernode_t
        integer :: first = 0, width = 0, m = 0
        integer, allocatable :: rows(:)
    end type supernode_t

    ! Destination-driven update lists and subtree task groups over a partition.
    ! Components are public for panel-granularity use by the numeric owners.
    type :: supernode_schedule_t
        integer :: n = 0, node_count = 0, active_threads = 1, task_groups = 0
        logical :: budget_reduced_threads = .false.
        real(dp) :: group_target_flops = 0.0_dp
        integer, allocatable :: colnode(:), tokens(:), rowmap(:,:)
        integer, allocatable :: incoming_ptr(:), incoming_source(:)
        integer, allocatable :: incoming_first(:), incoming_q(:)
        integer, allocatable :: group_member_ptr(:), group_members(:), group_tokens(:)
        integer, allocatable :: group_ptr(:), group_sources(:)
        real(dp), allocatable :: update_work(:,:)
    contains
        procedure :: request_threads
        procedure :: setup_incoming
        procedure :: setup_groups
        procedure :: fit_thread_scratch
        procedure :: allocate_rowmap
        procedure :: allocate_update_work
        procedure :: workspace_thread
        procedure :: update_group
        procedure :: memory_bytes
        procedure :: free => schedule_free
    end type supernode_schedule_t
contains
    pure subroutine supernode_partition(n, parent, count, strict, starts, ends, m, &
            node_count)
        ! Fundamental and relaxed amalgamation along elimination-tree parent
        ! chains. count excludes the diagonal. A merge of last+1 requires
        ! parent(last) = last+1; then the union of the column structures has
        ! width+count(last) rows, so the rule needs only parent and count.
        ! strict selects the compact exact-width limit (bounded workspaces).
        integer, intent(in) :: n, parent(:), count(:)
        logical, intent(in) :: strict
        integer, intent(out) :: starts(:), ends(:), m(:), node_count
        integer :: first, last, rows, candidate_m, candidate_w, exact_limit
        integer(i8) :: original_count, candidate_count, padded

        exact_limit = SUPERNODE_COMPACT_WIDTH
        if (.not. strict) then
            if (int(n,i8)+sum(int(count(:n),i8)) > 128_i8*n) &
                exact_limit = SUPERNODE_EXACT_WIDTH
        end if
        node_count = 0
        first = 1
        do while (first <= n)
            node_count = node_count+1
            last = first
            rows = count(first)+1
            original_count = rows
            do while (last < n)
                if (count(last) < 1) exit
                if (parent(last) /= last+1) exit
                candidate_w = last-first+2
                if (candidate_w > exact_limit) exit
                candidate_m = candidate_w+count(last+1)
                candidate_count = original_count+count(last+1)+1
                padded = int(candidate_m,i8)*candidate_w &
                    -int(candidate_w,i8)*(candidate_w-1)/2
                if (padded /= candidate_count) then
                    if (candidate_w > SUPERNODE_RELAXED_WIDTH) exit
                    if (4_i8*padded > 5_i8*candidate_count) exit
                end if
                last = last+1
                rows = candidate_m
                original_count = candidate_count
            end do
            starts(node_count) = first
            ends(node_count) = last
            m(node_count) = rows
            first = last+1
        end do
    end subroutine supernode_partition

    subroutine supernode_gather_rows(first, last, col_ptr, row_idx, stamp, mark, &
            rowlist, q)
        ! Sorted union of the rows >= first in columns first..last of a
        ! one-based column structure: the factor columns (L structure), or a
        ! symmetrized pattern merged with child supernode rows by the caller.
        ! mark(row) == stamp marks rows already present in rowlist(1:q); pass
        ! q = 0 for a fresh union. Entries below first are skipped.
        integer, intent(in) :: first, last, col_ptr(:), row_idx(:), stamp
        integer, intent(inout) :: mark(:), rowlist(:), q
        integer :: j, p, row
        do j = first, last
            do p = col_ptr(j), col_ptr(j+1)-1
                row = row_idx(p)
                if (row < first) cycle
                if (mark(row) == stamp) cycle
                mark(row) = stamp
                q = q+1
                rowlist(q) = row
            end do
        end do
        call supernode_sort_rows(rowlist, 1, q)
    end subroutine supernode_gather_rows

    recursive subroutine supernode_sort_rows(rows, lo, hi)
        integer, intent(inout) :: rows(:)
        integer, intent(in) :: lo, hi
        integer :: pivot, i, j, saved
        if (lo >= hi) return
        pivot = rows(lo+(hi-lo)/2)
        i = lo
        j = hi
        do
            do while (i <= hi)
                if (rows(i) >= pivot) exit
                i = i+1
            end do
            do while (j >= lo)
                if (rows(j) <= pivot) exit
                j = j-1
            end do
            if (i > j) exit
            saved = rows(i)
            rows(i) = rows(j)
            rows(j) = saved
            i = i+1
            j = j-1
        end do
        if (lo < j) call supernode_sort_rows(rows, lo, j)
        if (i < hi) call supernode_sort_rows(rows, i, hi)
    end subroutine supernode_sort_rows

    pure integer function supernode_find_row(node, row) result(position)
        ! Position of row among the sorted node%rows(width+1:m), or 0.
        type(supernode_t), intent(in) :: node
        integer, intent(in) :: row
        integer :: lo, hi, mid
        position = 0
        lo = node%width+1
        hi = node%m
        do while (lo <= hi)
            mid = lo+(hi-lo)/2
            if (node%rows(mid) == row) then
                position = mid
                return
            else if (node%rows(mid) < row) then
                lo = mid+1
            else
                hi = mid-1
            end if
        end do
    end function supernode_find_row

    subroutine supernode_charge(ledger, count, element_bytes, status, phase)
        type(memory_budget_t), intent(inout) :: ledger
        integer(i8), intent(in) :: count, element_bytes
        type(fortsparse_status_t), intent(out) :: status
        character(*), intent(in) :: phase
        integer(i8) :: bytes
        bytes = 0_i8
        call memory_checked_mul(count, element_bytes, bytes, status)
        if (status%code /= FORTSPARSE_OK) return
        call ledger%reserve(bytes, status, phase)
    end subroutine supernode_charge

    subroutine request_threads(self, nthreads, flops, serial)
        ! nthreads: 0 = OpenMP default, >0 explicit. Small factors and
        ! serial (device) executions keep one team member.
        class(supernode_schedule_t), intent(inout) :: self
        integer, intent(in) :: nthreads
        real(dp), intent(in) :: flops
        logical, intent(in) :: serial
        self%active_threads = 1
        if (nthreads == 0) then
            !$ self%active_threads = omp_get_max_threads()
        end if
        if (nthreads > 0) then
            !$ self%active_threads = nthreads
        end if
        if (flops < 1.0e6_dp) self%active_threads = 1
        if (serial) self%active_threads = 1
    end subroutine request_threads

    subroutine fit_thread_scratch(self, maxwork, retained_bytes, ledger, status)
        ! Halve the team until per-thread row maps and update scratch, plus
        ! the caller's retained non-thread scratch, fit the ledger.
        class(supernode_schedule_t), intent(inout) :: self
        integer, intent(in) :: maxwork
        integer(i8), intent(in) :: retained_bytes
        type(memory_budget_t), intent(inout) :: ledger
        type(fortsparse_status_t), intent(out) :: status
        type(memory_budget_t) :: probe
        integer(i8) :: rows, updates, per_thread, bytes
        rows = 0_i8
        updates = 0_i8
        per_thread = 0_i8
        bytes = 0_i8
        call memory_checked_mul(int(self%n,i8), 4_i8, rows, status)
        if (status%code /= FORTSPARSE_OK) return
        call memory_checked_mul(int(maxwork,i8), 8_i8, updates, status)
        if (status%code /= FORTSPARSE_OK) return
        call memory_checked_add(rows, updates, per_thread, status)
        if (status%code /= FORTSPARSE_OK) return
        do
            call memory_checked_mul(per_thread, int(self%active_threads,i8), &
                bytes, status)
            if (status%code /= FORTSPARSE_OK) return
            call memory_checked_add(bytes, retained_bytes, rows, status)
            if (status%code /= FORTSPARSE_OK) return
            probe = ledger
            call probe%reserve(rows, status, 'CPU supernode thread workspace')
            if (status%code == FORTSPARSE_OK) return
            if (status%code /= FORTSPARSE_RESOURCE_LIMIT) return
            if (self%active_threads == 1) then
                ledger = probe
                return
            end if
            self%active_threads = max(1,self%active_threads/2)
            self%budget_reduced_threads = .true.
        end do
    end subroutine fit_thread_scratch

    subroutine allocate_rowmap(self, ledger, status)
        class(supernode_schedule_t), intent(inout) :: self
        type(memory_budget_t), intent(inout) :: ledger
        type(fortsparse_status_t), intent(out) :: status
        integer :: ierr
        call supernode_charge(ledger, int(self%n,i8)*self%active_threads, 4_i8, &
            status, 'thread row maps')
        if (status%code /= FORTSPARSE_OK) return
        allocate(self%rowmap(self%n,self%active_threads), stat=ierr)
        if (ierr /= 0) then
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, 'supernode allocation failed')
            return
        end if
        self%rowmap = 0
    end subroutine allocate_rowmap

    subroutine allocate_update_work(self, maxwork, ledger, status)
        class(supernode_schedule_t), intent(inout) :: self
        integer, intent(in) :: maxwork
        type(memory_budget_t), intent(inout) :: ledger
        type(fortsparse_status_t), intent(out) :: status
        integer :: ierr
        call supernode_charge(ledger, int(maxwork,i8)*self%active_threads, 8_i8, &
            status, 'thread Schur scratch')
        if (status%code /= FORTSPARSE_OK) return
        allocate(self%update_work(maxwork,self%active_threads), stat=ierr)
        if (ierr /= 0) &
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, 'supernode allocation failed')
    end subroutine allocate_update_work

    integer function workspace_thread(self) result(thread)
        class(supernode_schedule_t), intent(in) :: self
        thread = 1
        if (self%active_threads > 1) then
            !$ thread = 1+omp_get_thread_num()
        end if
    end function workspace_thread

    subroutine update_group(self, node, first, group, dest)
        ! Rows node%rows(first:group-1) all belong to destination supernode dest.
        class(supernode_schedule_t), intent(in) :: self
        type(supernode_t), intent(in) :: node
        integer, intent(in) :: first
        integer, intent(out) :: group, dest
        dest = self%colnode(node%rows(first))
        group = first+1
        do while (group <= node%m)
            if (self%colnode(node%rows(group)) /= dest) exit
            group = group+1
        end do
    end subroutine update_group

    subroutine setup_incoming(self, nodes, ledger, status)
        class(supernode_schedule_t), intent(inout) :: self
        ! Polymorphic, so a numeric extension array is passed without a copy.
        class(supernode_t), intent(in) :: nodes(:)
        type(memory_budget_t), intent(inout) :: ledger
        type(fortsparse_status_t), intent(out) :: status
        integer, allocatable :: next(:)
        integer :: s, first, group, dest, p, ierr
        integer(i8) :: entries
        call supernode_charge(ledger, 2_i8*self%node_count+1_i8, 4_i8, status, &
            'incoming task headers and scratch')
        if (status%code /= FORTSPARSE_OK) return
        allocate(self%incoming_ptr(self%node_count+1), next(self%node_count), stat=ierr)
        if (ierr /= 0) then
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, 'incoming task map allocation failed')
            return
        end if
        self%incoming_ptr = 0
        entries = 0_i8
        do s = 1, self%node_count
            first = nodes(s)%width+1
            do while (first <= nodes(s)%m)
                call self%update_group(nodes(s), first, group, dest)
                self%incoming_ptr(dest+1) = self%incoming_ptr(dest+1)+1
                entries = entries+1
                first = group
            end do
        end do
        if (entries >= huge(p)) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, 'incoming task map too large')
            return
        end if
        self%incoming_ptr(1) = 1
        do s = 1, self%node_count
            self%incoming_ptr(s+1) = self%incoming_ptr(s+1)+self%incoming_ptr(s)
        end do
        call supernode_charge(ledger, entries, 12_i8, status, 'incoming task entries')
        if (status%code /= FORTSPARSE_OK) return
        allocate(self%incoming_source(int(entries)), self%incoming_first(int(entries)), &
            self%incoming_q(int(entries)), stat=ierr)
        if (ierr /= 0) then
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, 'incoming task entries allocation failed')
            return
        end if
        next = self%incoming_ptr(:self%node_count)
        do s = 1, self%node_count
            first = nodes(s)%width+1
            do while (first <= nodes(s)%m)
                call self%update_group(nodes(s), first, group, dest)
                p = next(dest)
                self%incoming_source(p) = s
                self%incoming_first(p) = first
                self%incoming_q(p) = group-first
                next(dest) = p+1
                first = group
            end do
        end do
        deallocate(next)
        call ledger%release(4_i8*self%node_count, status)
    end subroutine setup_incoming

    subroutine setup_groups(self, nodes, flops, ledger, status, lu)
        ! Requires setup_incoming. flops is the factor work that sets the
        ! group target; the team is then limited to a quarter of the groups.
        ! lu selects the unsymmetric (L and U halves) subtree work estimate,
        ! consistent with LU flops; the default is the Cholesky estimate.
        class(supernode_schedule_t), intent(inout) :: self
        class(supernode_t), intent(in) :: nodes(:)
        real(dp), intent(in) :: flops
        type(memory_budget_t), intent(inout) :: ledger
        type(fortsparse_status_t), intent(out) :: status
        logical, optional, intent(in) :: lu
        integer, allocatable :: node_group(:), parent(:), root_group(:), mark(:), next(:)
        real(dp), allocatable :: subtree_work(:)
        integer :: s, p, source, g, predecessor, ierr, position, member
        integer(i8) :: entries
        logical :: unsymmetric

        unsymmetric = .false.
        if (present(lu)) unsymmetric = lu
        ! Collapse complete small elimination subtrees. Grouping arbitrary
        ! consecutive nodes can introduce false chains between independent
        ! subtrees. Large separator/critical panels remain individual tasks.
        self%group_target_flops = max(1.0e6_dp, &
            flops/(12.0_dp*self%active_threads))
        call supernode_charge(ledger, int(self%node_count,i8), 20_i8, status, &
            'subtree group scratch')
        if (status%code /= FORTSPARSE_OK) return
        allocate(node_group(self%node_count), parent(self%node_count), &
            root_group(self%node_count), subtree_work(self%node_count), stat=ierr)
        if (ierr /= 0) then
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, 'task group allocation failed')
            return
        end if
        do s = 1, self%node_count
            parent(s) = 0
            if (nodes(s)%m > nodes(s)%width) &
                parent(s) = self%colnode(nodes(s)%rows(nodes(s)%width+1))
            if (unsymmetric) then
                ! Block LU, both triangular solves and both update GEMMs.
                subtree_work(s) = 2.0_dp/3.0_dp*real(nodes(s)%width,dp)**3 &
                    +2.0_dp*real(nodes(s)%m-nodes(s)%width,dp)*nodes(s)%width**2
            else
                subtree_work(s) = real(nodes(s)%m,dp)*nodes(s)%width**2
            end if
            do p = self%incoming_ptr(s), self%incoming_ptr(s+1)-1
                source = self%incoming_source(p)
                subtree_work(s) = subtree_work(s) &
                    +2.0_dp*(nodes(source)%m-self%incoming_first(p)+1) &
                    *self%incoming_q(p)*nodes(source)%width
                if (unsymmetric) subtree_work(s) = subtree_work(s) &
                    +2.0_dp*(nodes(source)%m-self%incoming_first(p)+1 &
                    -self%incoming_q(p))*self%incoming_q(p)*nodes(source)%width
            end do
        end do
        ! Parent indices are larger: this ascending pass accumulates every
        ! completed child's work before its parent contributes farther up.
        do s = 1, self%node_count
            p = parent(s)
            if (p > 0) subtree_work(p) = subtree_work(p)+subtree_work(s)
        end do
        do s = self%node_count, 1, -1
            node_group(s) = s
            p = parent(s)
            if (p > 0) then
                if (subtree_work(p) <= self%group_target_flops) node_group(s) = node_group(p)
            end if
        end do
        root_group = 0
        do s = 1, self%node_count
            if (node_group(s) /= s) cycle
            self%task_groups = self%task_groups+1
            root_group(s) = self%task_groups
        end do
        ! Roots are each group's largest panel, so ascending root order is
        ! also a topological creation order for inter-group dependencies.
        do s = 1, self%node_count
            node_group(s) = root_group(node_group(s))
        end do
        call supernode_charge(ledger, int(self%node_count,i8)+5_i8*self%task_groups+2_i8, &
            4_i8, status, 'subtree group headers and scratch')
        if (status%code /= FORTSPARSE_OK) return
        allocate(self%group_member_ptr(self%task_groups+1), self%group_members(self%node_count), &
            self%group_tokens(self%task_groups), self%group_ptr(self%task_groups+1), &
            mark(self%task_groups), next(self%task_groups), stat=ierr)
        if (ierr /= 0) then
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, 'task group headers allocation failed')
            return
        end if
        self%group_member_ptr = 0
        do s = 1, self%node_count
            g = node_group(s)
            self%group_member_ptr(g+1) = self%group_member_ptr(g+1)+1
        end do
        self%group_member_ptr(1) = 1
        do g = 1, self%task_groups
            self%group_member_ptr(g+1) = self%group_member_ptr(g+1)+self%group_member_ptr(g)
        end do
        next = self%group_member_ptr(:self%task_groups)
        do s = 1, self%node_count
            g = node_group(s)
            self%group_members(next(g)) = s
            next(g) = next(g)+1
        end do
        self%group_tokens = 0
        self%group_ptr = 0
        mark = 0
        entries = 0_i8
        do g = 1, self%task_groups
            do member = self%group_member_ptr(g), self%group_member_ptr(g+1)-1
                s = self%group_members(member)
                do p = self%incoming_ptr(s), self%incoming_ptr(s+1)-1
                    predecessor = node_group(self%incoming_source(p))
                    if (predecessor == g) cycle
                    if (predecessor > g) then
                        call status_set(status, FORTSPARSE_INTERNAL_ERROR, &
                            'invalid elimination-subtree task dependency')
                        return
                    end if
                    if (mark(predecessor) == g) cycle
                    mark(predecessor) = g
                    self%group_ptr(g+1) = self%group_ptr(g+1)+1
                    entries = entries+1
                end do
            end do
        end do
        if (entries >= huge(position)) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, 'task group dependencies too large')
            return
        end if
        self%group_ptr(1) = 1
        do g = 1, self%task_groups
            self%group_ptr(g+1) = self%group_ptr(g+1)+self%group_ptr(g)
        end do
        call supernode_charge(ledger, entries, 4_i8, status, 'subtree group dependencies')
        if (status%code /= FORTSPARSE_OK) return
        allocate(self%group_sources(int(entries)), stat=ierr)
        if (ierr /= 0) then
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, 'task group edges allocation failed')
            return
        end if
        mark = 0
        do g = 1, self%task_groups
            position = self%group_ptr(g)
            do member = self%group_member_ptr(g), self%group_member_ptr(g+1)-1
                s = self%group_members(member)
                do p = self%incoming_ptr(s), self%incoming_ptr(s+1)-1
                    predecessor = node_group(self%incoming_source(p))
                    if (predecessor == g) cycle
                    if (mark(predecessor) == g) cycle
                    mark(predecessor) = g
                    self%group_sources(position) = predecessor
                    position = position+1
                end do
            end do
        end do
        deallocate(node_group, parent, root_group, subtree_work, mark, next)
        call ledger%release(20_i8*self%node_count+8_i8*self%task_groups, status)
        if (status%code /= FORTSPARSE_OK) return
        self%active_threads = min(self%active_threads, max(1,self%task_groups/4))
    end subroutine setup_groups

    integer(i8) function memory_bytes(self) result(bytes)
        ! Owned array payload only; descriptors and OpenMP runtime excluded.
        class(supernode_schedule_t), intent(in) :: self
        bytes = 0_i8
        if (allocated(self%colnode)) bytes = bytes+4_i8*size(self%colnode,kind=i8)
        if (allocated(self%tokens)) bytes = bytes+4_i8*size(self%tokens,kind=i8)
        if (allocated(self%rowmap)) bytes = bytes+4_i8*size(self%rowmap,kind=i8)
        if (allocated(self%incoming_ptr)) bytes = bytes+4_i8*size(self%incoming_ptr,kind=i8)
        if (allocated(self%incoming_source)) bytes = bytes+4_i8*size(self%incoming_source,kind=i8)
        if (allocated(self%incoming_first)) bytes = bytes+4_i8*size(self%incoming_first,kind=i8)
        if (allocated(self%incoming_q)) bytes = bytes+4_i8*size(self%incoming_q,kind=i8)
        if (allocated(self%group_member_ptr)) bytes = bytes+4_i8*size(self%group_member_ptr,kind=i8)
        if (allocated(self%group_members)) bytes = bytes+4_i8*size(self%group_members,kind=i8)
        if (allocated(self%group_tokens)) bytes = bytes+4_i8*size(self%group_tokens,kind=i8)
        if (allocated(self%group_ptr)) bytes = bytes+4_i8*size(self%group_ptr,kind=i8)
        if (allocated(self%group_sources)) bytes = bytes+4_i8*size(self%group_sources,kind=i8)
        if (allocated(self%update_work)) bytes = bytes+8_i8*size(self%update_work,kind=i8)
    end function memory_bytes

    subroutine schedule_free(self)
        class(supernode_schedule_t), intent(inout) :: self
        if (allocated(self%colnode)) deallocate(self%colnode)
        if (allocated(self%tokens)) deallocate(self%tokens)
        if (allocated(self%rowmap)) deallocate(self%rowmap)
        if (allocated(self%incoming_ptr)) deallocate(self%incoming_ptr)
        if (allocated(self%incoming_source)) deallocate(self%incoming_source)
        if (allocated(self%incoming_first)) deallocate(self%incoming_first)
        if (allocated(self%incoming_q)) deallocate(self%incoming_q)
        if (allocated(self%group_member_ptr)) deallocate(self%group_member_ptr)
        if (allocated(self%group_members)) deallocate(self%group_members)
        if (allocated(self%group_tokens)) deallocate(self%group_tokens)
        if (allocated(self%group_ptr)) deallocate(self%group_ptr)
        if (allocated(self%group_sources)) deallocate(self%group_sources)
        if (allocated(self%update_work)) deallocate(self%update_work)
        self%n = 0
        self%node_count = 0
        self%active_threads = 1
        self%task_groups = 0
        self%group_target_flops = 0.0_dp
        self%budget_reduced_threads = .false.
    end subroutine schedule_free
end module fortsparse_native_supernode_schedule
