module fortsparse_native_slu
    ! Static-structure supernodal LU (docs/lu-supernodal.md), real values.
    !
    ! Factor M = Q'^T Dr Pm A Dc Q' = L U with the maximum-product matching
    ! Pm and its dual scaling Dr, Dc (fortsparse_native_mpm), a fill-reducing
    ! order Q of the symmetrized matched pattern and its elimination-tree
    ! postorder (Q' = Q then postorder). The factor structure is the static
    ! symmetric structure of M + M^T, so the L and U halves share one
    ! supernode partition, row lists, row maps and update schedule.
    !
    ! Retained index vectors are in the final (new) labels k = 1..n:
    !   col_old(k)  original column of new column k (Q');
    !   row_old(k)  original row of new row k, = cmatch(col_old(k));
    !   row_new(i)  new row of original row i (inverse of row_old);
    !   drp(k), dcp(k) the row and column scale of new row/column k,
    ! so that M(k,l) = drp(k) * A(row_old(k), col_old(l)) * dcp(l).
    !
    ! analyze performs every allocation, including the numeric panels and all
    ! scratch, after charging each phase to the ledger. factor, refactor and
    ! solve allocate nothing. All components are allocatable: intrinsic
    ! assignment produces an independent deep copy.
    use fortsparse_kinds, only: dp, i8
    use fortsparse_memory, only: memory_budget_t, memory_checked_add, &
        memory_checked_mul, sparse_memory_plan_t, MEMORY_PLAN_EXACT_LAYOUT
    use fortsparse_status, only: fortsparse_status_t, status_set, status_ok, &
        FORTSPARSE_OK, FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OUT_OF_MEMORY, &
        FORTSPARSE_NOT_FACTORED, FORTSPARSE_INDEX_OVERFLOW, &
        FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_INTERNAL_ERROR, &
        FORTSPARSE_RESOURCE_LIMIT, FORTSPARSE_SINGULAR
    use fortsparse_native_mpm, only: mpm_stats_t, mpm_iwork_size, &
        mpm_rwork_size, mpm_workspace_bytes, max_product_matching, mpm_scaling
    use fortsparse_native_matching, only: matching_graph_rows
    use fortsparse_native_amd, only: bounded_amd_available, bounded_amd_size, &
        bounded_amd_order
    use fortsparse_native_symbolic, only: native_symbolic_counts, &
        native_symbolic_workspace_bytes
    use fortsparse_native_supernode_schedule, only: supernode_t, &
        supernode_schedule_t, supernode_partition, supernode_gather_rows, &
        supernode_sort_rows
    use, intrinsic :: iso_c_binding, only: c_int, c_int64_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: slu_factor_t

    integer, parameter, public :: SLU_ORDER_AMD = 1, SLU_ORDER_ND = 2
    ! Ledger phases, in analysis order (failure reports name the phase).
    integer, parameter, public :: SLU_PHASE_MPM = 1, SLU_PHASE_GRAPH = 2, &
        SLU_PHASE_ORDER = 3, SLU_PHASE_COUNTS = 4, SLU_PHASE_ROWS = 5, &
        SLU_PHASE_INCOMING = 6, SLU_PHASE_PANELS = 7, SLU_PHASE_SCRATCH = 8
    integer, parameter, public :: SLU_PHASE_COUNT = 8
    character(*), parameter :: phase_names(SLU_PHASE_COUNT) = [character(16) :: &
        'MPM', 'graph', 'AMD arena', 'counts', 'rows', 'incoming', 'panels', &
        'scratch']

    ! L_s(i,j) = values(offset_l + (j-1)*m + i), i = 1..m, j = 1..width;
    ! Ut_s(r,c) = values(offset_u + (c-1)*(m-width) + r) = U(first+c-1,
    ! rows(width+r)). offset_u = offset_l + m*width: Ut follows L.
    type, extends(supernode_t) :: slu_panel_t
        integer(i8) :: offset_l = 0_i8, offset_u = 0_i8
    end type slu_panel_t

    type :: slu_factor_t
        integer :: n = 0, node_count = 0, active_threads = 1
        ! factor_nnz: padded L+U panel entries; symbolic_nnz: exact
        ! nnz(L+U) of the symmetrized structure (L diagonal counted once).
        integer(i8) :: factor_nnz = 0_i8, symbolic_nnz = 0_i8
        real(dp) :: flops = 0.0_dp
        integer :: perturbed_pivots = 0, block_swaps = 0
        ! growth = max|U| / max|M|; min_pivot = min |u_kk| (scaled matrix).
        real(dp) :: max_growth = 0.0_dp, min_pivot = 0.0_dp
        real(dp) :: scaling_log10_range = 0.0_dp
        ! plan_bytes is the exact owned-array layout computed from the
        ! structure before the panels and scratch are allocated; it equals
        ! memory_bytes() afterwards. peak_memory_bytes is the ledger
        ! high-water of the analysis above the caller's baseline.
        integer(i8) :: plan_bytes = 0_i8, peak_memory_bytes = 0_i8
        integer(i8) :: failed_required_bytes = 0_i8, failed_available_bytes = 0_i8
        integer :: failed_phase = 0
        ! Opaque (unbounded) ordering workspace: no exact peak claim.
        logical :: opaque_ordering = .false.
        real(dp) :: analysis_s = 0.0_dp, factor_s = 0.0_dp
        logical, private :: analyzed = .false., factored = .false.
        logical, private :: strict = .false., use_blas = .true.
        integer, private :: inject = 0, maxwork = 0, maxrows = 0
        real(dp), private :: tau = 1.0e-13_dp, tau_abs = 0.0_dp, max_abs_m = 0.0_dp
        integer(i8), private :: owned_budget_bytes = 0_i8
        type(supernode_schedule_t), private :: sched
        type(slu_panel_t), allocatable, private :: panels(:)
        integer, allocatable, private :: pivots(:), failures(:)
        integer, allocatable, private :: row_old(:), row_new(:), col_old(:)
        real(dp), allocatable, private :: values(:), drp(:), dcp(:), solve_work(:)
    contains
        procedure :: analyze
        procedure :: free => slu_free
        procedure :: memory_bytes
        procedure :: plan
        procedure :: is_analyzed
        procedure :: node_structure
        procedure :: ordering
        procedure, private :: structure_bytes
    end type slu_factor_t

    interface
        integer(c_int) function order_amd_general(n, nnz, ptr, row, perm, inverse) &
            bind(C, name='fsparse_native_order_amd_general')
            import :: c_int
            integer(c_int), value :: n, nnz
            integer(c_int), intent(in) :: ptr(*), row(*)
            integer(c_int), intent(out) :: perm(*), inverse(*)
        end function order_amd_general
        integer(c_int) function order_nd_general(n, nnz, ptr, row, perm, inverse) &
            bind(C, name='fsparse_native_order_general')
            import :: c_int
            integer(c_int), value :: n, nnz
            integer(c_int), intent(in) :: ptr(*), row(*)
            integer(c_int), intent(out) :: perm(*), inverse(*)
        end function order_nd_general
    end interface
contains
    subroutine analyze(self, n, ptr, rows, values, status, ordering, low_memory, &
            pivot_perturbation, budget, inject_phase)
        ! Matching, scaling, ordering, symbolic structure, schedule and the
        ! allocation of all retained arrays. values are read only by the
        ! matching and scaling; factor assembles them later.
        ! ordering: SLU_ORDER_AMD (default) or SLU_ORDER_ND. ND (METIS) has
        ! opaque workspace and is refused under a positive budget limit or
        ! low_memory. AMD uses the bounded caller-workspace adapter when it is
        ! available; otherwise opaque AMD is refused under a limit.
        ! budget: caller ledger including its live allocations. Every phase
        ! is reserved before allocation; on failure the ledger current is
        ! restored, failed_required/available and failed_phase are set.
        ! inject_phase (tests only): refuse the first reservation of that
        ! phase as if the cap were one byte short.
        class(slu_factor_t), intent(inout) :: self
        integer, intent(in) :: n, ptr(:), rows(:)
        real(dp), intent(in) :: values(:)
        type(fortsparse_status_t), intent(out) :: status
        integer, optional, intent(in) :: ordering
        logical, optional, intent(in) :: low_memory
        real(dp), optional, intent(in) :: pivot_perturbation
        type(memory_budget_t), optional, intent(inout) :: budget
        integer, optional, intent(in) :: inject_phase
        type(memory_budget_t) :: ledger
        integer(i8) :: baseline
        integer :: order_kind
        logical :: low
        real(dp) :: started

        started = seconds()
        call self%free(budget, status)
        if (.not. status_ok(status)) return
        order_kind = SLU_ORDER_AMD
        if (present(ordering)) order_kind = ordering
        low = .false.
        if (present(low_memory)) low = low_memory
        if (present(pivot_perturbation)) self%tau = pivot_perturbation
        if (present(inject_phase)) self%inject = inject_phase
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            'supernodal LU: invalid ordering or pivot perturbation')
        if (order_kind /= SLU_ORDER_AMD .and. order_kind /= SLU_ORDER_ND) return
        if (.not. ieee_is_finite(self%tau)) return
        if (.not. self%tau > 0.0_dp .or. self%tau >= 1.0_dp) return
        if (present(budget)) ledger = budget
        call ledger%reserve(0_i8, status)
        if (.not. status_ok(status)) return
        self%strict = ledger%limit > 0_i8
        baseline = ledger%current
        call analyze_phases(self, n, ptr, rows, values, order_kind, low, ledger, &
            status)
        if (status_ok(status)) then
            self%owned_budget_bytes = self%memory_bytes()
            if (ledger%current-baseline /= self%owned_budget_bytes .or. &
                self%plan_bytes /= self%owned_budget_bytes) &
                call status_set(status, FORTSPARSE_INTERNAL_ERROR, &
                    'supernodal LU: ledger disagrees with the retained layout')
        end if
        self%peak_memory_bytes = ledger%peak-baseline
        if (.not. status_ok(status)) then
            if (status%code == FORTSPARSE_RESOURCE_LIMIT) then
                self%failed_required_bytes = ledger%failed_required
                self%failed_available_bytes = ledger%failed_available
            end if
            call release_on_failure(self)
            ledger%current = baseline
        else
            self%analyzed = .true.
        end if
        self%inject = 0
        if (present(budget)) budget = ledger
        self%analysis_s = seconds()-started
    end subroutine analyze

    subroutine release_on_failure(self)
        ! Keep the failure diagnostics; release every owned array.
        class(slu_factor_t), intent(inout) :: self
        integer(i8) :: required, available, peak
        integer :: phase
        logical :: opaque
        required = self%failed_required_bytes
        available = self%failed_available_bytes
        peak = self%peak_memory_bytes
        phase = self%failed_phase
        opaque = self%opaque_ordering
        call self%free()
        self%failed_required_bytes = required
        self%failed_available_bytes = available
        self%peak_memory_bytes = peak
        self%failed_phase = phase
        self%opaque_ordering = opaque
    end subroutine release_on_failure

    subroutine charge(self, ledger, bytes, phase, status)
        class(slu_factor_t), intent(inout) :: self
        type(memory_budget_t), intent(inout) :: ledger
        integer(i8), intent(in) :: bytes
        integer, intent(in) :: phase
        type(fortsparse_status_t), intent(out) :: status
        if (self%inject == phase .and. bytes > 0_i8) then
            ! Test injection: the cap is exactly one byte short.
            self%inject = 0
            ledger%failed_required = bytes
            ledger%failed_available = bytes-1_i8
            self%failed_phase = phase
            call status_set(status, FORTSPARSE_RESOURCE_LIMIT, &
                'memory cap exceeded: supernodal LU '//trim(phase_names(phase)))
            return
        end if
        call ledger%reserve(bytes, status, 'supernodal LU '//trim(phase_names(phase)))
        if (status%code == FORTSPARSE_RESOURCE_LIMIT) self%failed_phase = phase
    end subroutine charge

    subroutine charge_count(self, ledger, count, element_bytes, phase, status)
        class(slu_factor_t), intent(inout) :: self
        type(memory_budget_t), intent(inout) :: ledger
        integer(i8), intent(in) :: count, element_bytes
        integer, intent(in) :: phase
        type(fortsparse_status_t), intent(out) :: status
        integer(i8) :: bytes
        bytes = 0_i8
        call memory_checked_mul(count, element_bytes, bytes, status)
        if (.not. status_ok(status)) return
        call charge(self, ledger, bytes, phase, status)
    end subroutine charge_count

    subroutine release(ledger, bytes, status)
        type(memory_budget_t), intent(inout) :: ledger
        integer(i8), intent(in) :: bytes
        type(fortsparse_status_t), intent(out) :: status
        call ledger%release(bytes, status)
    end subroutine release

    subroutine analyze_phases(self, n, ptr, rows, values, order_kind, low, ledger, &
            status)
        class(slu_factor_t), intent(inout) :: self
        integer, intent(in) :: n, ptr(:), rows(:), order_kind
        real(dp), intent(in) :: values(:)
        logical, intent(in) :: low
        type(memory_budget_t), intent(inout) :: ledger
        type(fortsparse_status_t), intent(out) :: status
        integer, allocatable :: iwork(:), graph(:), order(:), finv(:)
        integer, allocatable :: sptr(:), sidx(:), parent(:), count(:), post(:)
        integer, allocatable :: head(:), next(:), stack(:)
        integer(c_int), allocatable :: arena(:)
        real(dp), allocatable :: rwork(:), u(:), v(:)
        type(mpm_stats_t) :: mstats
        integer :: nnz, ierr, code, k, len, lowflag
        integer(i8) :: bytes, mpm_bytes, arena_bytes, output_bytes, symbolic_bytes
        integer(i8) :: sym_entries
        logical :: cache

        call validate_csc(n, ptr, rows, values, status)
        if (.not. status_ok(status)) return
        nnz = ptr(n+1)-1
        self%n = n

        ! Phase MPM: retained scaling/permutation vectors (28n bytes), then
        ! the matching workspace and duals (released after scaling).
        call charge_count(self, ledger, int(n,i8), 28_i8, SLU_PHASE_MPM, status)
        if (.not. status_ok(status)) return
        allocate(self%drp(n), self%dcp(n), self%row_old(n), self%row_new(n), &
            self%col_old(n), stat=ierr)
        if (ierr /= 0) goto 910
        cache = .not. low
        mpm_bytes = mpm_workspace_bytes(n, nnz, cache)
        if (mpm_bytes < 0_i8) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                'supernodal LU: matching workspace index range')
            return
        end if
        bytes = 0_i8
        call memory_checked_add(mpm_bytes, 16_i8*int(n,i8), bytes, status)
        if (.not. status_ok(status)) return
        call charge(self, ledger, bytes, SLU_PHASE_MPM, status)
        if (.not. status_ok(status)) return
        allocate(iwork(mpm_iwork_size(n)), rwork(mpm_rwork_size(n, nnz, cache)), &
            u(n), v(n), stat=ierr)
        if (ierr /= 0) goto 910
        ! row_old temporarily holds cmatch (column -> original row) and
        ! row_new the inverse rmatch; drp/dcp hold dr/dc in original labels.
        call max_product_matching(n, ptr, rows, values, self%row_old, u, v, iwork, &
            rwork, status, mstats)
        if (.not. status_ok(status)) return
        self%row_new = iwork(1:n)
        call mpm_scaling(n, ptr, rows, values, u, v, self%drp, self%dcp, status)
        if (.not. status_ok(status)) then
            if (status%code /= FORTSPARSE_INVALID_MATRIX) &
                call status_set(status, FORTSPARSE_SINGULAR, status%msg)
            return
        end if
        self%scaling_log10_range = max(mstats%log_range_dr, mstats%log_range_dc)
        deallocate(iwork, rwork, u, v)
        call release(ledger, bytes, status)
        if (.not. status_ok(status)) return

        ! Phase graph: row labels of Pm*A (matched indices = column labels).
        call charge_count(self, ledger, int(nnz,i8), 4_i8, SLU_PHASE_GRAPH, status)
        if (.not. status_ok(status)) return
        allocate(graph(nnz), stat=ierr)
        if (ierr /= 0) goto 910
        call matching_graph_rows(ptr, rows, self%row_new, graph)

        ! Phase ordering on the pattern of Pm*A + (Pm*A)^T.
        call charge_count(self, ledger, int(n,i8), 8_i8, SLU_PHASE_ORDER, status)
        if (.not. status_ok(status)) return
        allocate(order(n), finv(n), stat=ierr)
        if (ierr /= 0) goto 910
        if (order_kind == SLU_ORDER_ND) then
            if (self%strict .or. low) then
                call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                    'supernodal LU: opaque METIS workspace cannot enforce a cap or LOW')
                return
            end if
            self%opaque_ordering = .true.
            code = order_nd_general(n, nnz, ptr, graph, order, finv)
        else if (bounded_amd_available() == 1 .and. &
                n < huge(0_c_int)/4 .and. nnz < huge(0_c_int)/4) then
            lowflag = 0
            if (low) lowflag = 1
            code = bounded_amd_size(n, nnz, ptr, graph, lowflag, arena_bytes, &
                output_bytes)
            call status_set(status, code, &
                'supernodal LU: bounded AMD graph/ABI is unavailable or invalid')
            if (.not. status_ok(status)) return
            ! The two n-vectors of output_bytes are charged above.
            call charge(self, ledger, arena_bytes, SLU_PHASE_ORDER, status)
            if (.not. status_ok(status)) return
            allocate(arena(arena_bytes/4_i8), stat=ierr)
            if (ierr /= 0) goto 910
            code = bounded_amd_order(n, nnz, ptr, graph, lowflag, arena, arena_bytes, &
                order, finv)
            deallocate(arena)
            call release(ledger, arena_bytes, status)
            if (.not. status_ok(status)) return
        else
            if (self%strict .or. low) then
                call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                    'supernodal LU: audited AMD caller-workspace ABI unavailable')
                return
            end if
            self%opaque_ordering = .true.
            code = order_amd_general(n, nnz, ptr, graph, order, finv)
        end if
        call status_set(status, code, 'supernodal LU: column ordering failed or unavailable')
        if (.not. status_ok(status)) return

        ! Phase counts: symmetrized pattern, etree/counts and postorder.
        sym_entries = 2_i8*int(nnz,i8)+int(n,i8)
        if (sym_entries >= int(huge(0),i8)) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                'supernodal LU: symmetrized pattern index range')
            return
        end if
        call charge_count(self, ledger, sym_entries+int(n,i8)+1_i8+6_i8*n, 4_i8, &
            SLU_PHASE_COUNTS, status)
        if (.not. status_ok(status)) return
        allocate(sptr(n+1), sidx(int(sym_entries)), parent(n), count(n), post(n), &
            head(n), next(n), stack(n), stat=ierr)
        if (ierr /= 0) goto 910
        call build_symmetric(n, ptr, graph, finv, sptr, sidx, stack, len)
        symbolic_bytes = native_symbolic_workspace_bytes(n)
        call charge(self, ledger, symbolic_bytes, SLU_PHASE_COUNTS, status)
        if (.not. status_ok(status)) return
        call native_symbolic_counts(n, sptr, sidx(1:len), parent, count, status)
        if (.not. status_ok(status)) return
        call release(ledger, symbolic_bytes, status)
        if (.not. status_ok(status)) return
        call etree_postorder(n, parent, count, post, head, next, stack)
        ! Final order: new k is matched index order(post(k)). Relabel the
        ! forest; a topological renumbering keeps the same filled graph.
        do k = 1, n
            self%col_old(k) = order(post(k))
            head(post(k)) = k
        end do
        do k = 1, n
            next(k) = 0
            if (parent(post(k)) > 0) next(k) = head(parent(post(k)))
            stack(k) = count(post(k))
        end do
        parent = next
        count = stack
        do k = 1, n
            finv(self%col_old(k)) = k
        end do
        call build_symmetric(n, ptr, graph, finv, sptr, sidx, stack, len)
        deallocate(graph, order)
        call release(ledger, 4_i8*int(nnz,i8)+4_i8*int(n,i8), status)
        if (.not. status_ok(status)) return
        ! The postorder is checked here; a violation is an internal error.
        do k = 1, n
            if (parent(k) /= 0 .and. parent(k) <= k) then
                call status_set(status, FORTSPARSE_INTERNAL_ERROR, &
                    'supernodal LU: elimination tree not postordered')
                return
            end if
        end do
        self%symbolic_nnz = 2_i8*sum(int(count,i8))+int(n,i8)

        call build_nodes(self, n, sptr, sidx, len, parent, count, head, next, &
            stack, finv, ledger, status)
        if (.not. status_ok(status)) return
        deallocate(sptr, sidx, parent, count, post, head, next, stack, finv)
        call release(ledger, 4_i8*(sym_entries+int(n,i8)+1_i8+7_i8*n), status)
        if (.not. status_ok(status)) return

        call build_schedule(self, ledger, status)
        if (.not. status_ok(status)) return
        call allocate_numeric(self, ledger, status)
        if (.not. status_ok(status)) return
        call finalize_labels(self)
        call status_set(status, FORTSPARSE_OK, '')
        return
910     call status_set(status, FORTSPARSE_OUT_OF_MEMORY, &
            'supernodal LU: analysis allocation failed')
    end subroutine analyze_phases

    subroutine validate_csc(n, ptr, rows, values, status)
        integer, intent(in) :: n, ptr(:), rows(:)
        real(dp), intent(in) :: values(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: j, p
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            'supernodal LU: require a square sorted CSC with finite values')
        if (n <= 0 .or. n >= huge(n)/7) return
        if (size(ptr) /= n+1) return
        if (ptr(1) /= 1 .or. ptr(n+1)-1 /= size(rows)) return
        if (size(values) /= size(rows)) return
        do j = 1, n
            if (ptr(j+1) < ptr(j)) return
            do p = ptr(j), ptr(j+1)-1
                if (rows(p) < 1 .or. rows(p) > n) return
                if (p > ptr(j)) then
                    if (rows(p) <= rows(p-1)) return
                end if
                if (.not. ieee_is_finite(values(p))) return
            end do
        end do
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine validate_csc

    subroutine build_symmetric(n, ptr, graph, inverse, sptr, sidx, fill, len)
        ! Sorted, duplicate-free full CSC of S = P (B + B^T + I) P^T, where B
        ! is the pattern (ptr, graph) and inverse maps its labels to the new
        ! ones. Capacity 2*nnz(B)+n; len returns the used entries.
        integer, intent(in) :: n, ptr(:), graph(:), inverse(:)
        integer, intent(out) :: sptr(:), sidx(:), fill(:), len
        integer :: c, e, k, l, q, start, finish, last
        sptr = 0
        do c = 1, n
            l = inverse(c)
            sptr(l+1) = sptr(l+1)+1
            do e = ptr(c), ptr(c+1)-1
                k = inverse(graph(e))
                if (k == l) cycle
                sptr(l+1) = sptr(l+1)+1
                sptr(k+1) = sptr(k+1)+1
            end do
        end do
        sptr(1) = 1
        do l = 1, n
            sptr(l+1) = sptr(l+1)+sptr(l)
        end do
        fill(1:n) = sptr(1:n)
        do c = 1, n
            l = inverse(c)
            sidx(fill(l)) = l
            fill(l) = fill(l)+1
            do e = ptr(c), ptr(c+1)-1
                k = inverse(graph(e))
                if (k == l) cycle
                sidx(fill(l)) = k
                fill(l) = fill(l)+1
                sidx(fill(k)) = l
                fill(k) = fill(k)+1
            end do
        end do
        q = 0
        start = 1
        do l = 1, n
            finish = sptr(l+1)-1
            call supernode_sort_rows(sidx, start, finish)
            sptr(l) = q+1
            last = 0
            do e = start, finish
                if (sidx(e) == last) cycle
                q = q+1
                sidx(q) = sidx(e)
                last = sidx(e)
            end do
            start = finish+1
        end do
        sptr(n+1) = q+1
        len = q
    end subroutine build_symmetric

    subroutine etree_postorder(n, parent, count, post, head, next, stack)
        ! Depth-first postorder of the elimination forest; children are
        ! visited in ascending order except that the child with the largest
        ! column count (the last such on ties) is visited last.
        integer, intent(in) :: n, parent(:), count(:)
        integer, intent(out) :: post(:), head(:), next(:), stack(:)
        integer :: j, p, child, best, previous, top, k, i
        head(1:n) = 0
        next(1:n) = 0
        do j = n, 1, -1
            p = parent(j)
            if (p == 0) cycle
            next(j) = head(p)
            head(p) = j
        end do
        do p = 1, n
            child = head(p)
            if (child == 0) cycle
            best = child
            do while (child /= 0)
                if (count(child) >= count(best)) best = child
                child = next(child)
            end do
            if (next(best) == 0) cycle
            ! Unlink best and append it at the tail.
            if (head(p) == best) then
                head(p) = next(best)
            else
                previous = head(p)
                do while (next(previous) /= best)
                    previous = next(previous)
                end do
                next(previous) = next(best)
            end if
            child = head(p)
            do while (next(child) /= 0)
                child = next(child)
            end do
            next(child) = best
            next(best) = 0
        end do
        k = 0
        do j = 1, n
            if (parent(j) /= 0) cycle
            top = 1
            stack(1) = j
            do while (top > 0)
                p = stack(top)
                i = head(p)
                if (i == 0) then
                    top = top-1
                    k = k+1
                    post(k) = p
                else
                    head(p) = next(i)
                    top = top+1
                    stack(top) = i
                end if
            end do
        end do
    end subroutine etree_postorder

    subroutine build_nodes(self, n, sptr, sidx, len, parent, count, starts, ends, &
            msize, mark, ledger, status)
        ! Partition, descriptors, column ownership and per-node row lists.
        ! starts/ends/msize/mark are caller scratch of length n.
        class(slu_factor_t), intent(inout) :: self
        integer, intent(in) :: n, sptr(:), sidx(:), len, parent(:), count(:)
        integer, intent(inout) :: starts(:), ends(:), msize(:), mark(:)
        type(memory_budget_t), intent(inout) :: ledger
        type(fortsparse_status_t), intent(out) :: status
        integer, allocatable :: rowlist(:), chead(:), cnext(:)
        integer :: s, ns, first, last, m, w, q, c, p, r, ierr
        integer(i8) :: descriptor_bytes

        call supernode_partition(n, parent, count, self%strict, starts, ends, msize, ns)
        self%node_count = ns
        descriptor_bytes = int(storage_size(slu_panel_t())/8,i8)
        ! Descriptors, task tokens and failure flags, column ownership.
        call charge_count(self, ledger, int(ns,i8), descriptor_bytes+8_i8, &
            SLU_PHASE_ROWS, status)
        if (.not. status_ok(status)) return
        call charge_count(self, ledger, int(n,i8), 4_i8, SLU_PHASE_ROWS, status)
        if (.not. status_ok(status)) return
        allocate(self%panels(ns), self%sched%tokens(ns), self%failures(ns), &
            self%sched%colnode(n), stat=ierr)
        if (ierr /= 0) goto 910
        self%sched%n = n
        self%sched%node_count = ns
        self%sched%tokens = 0
        self%failures = 0
        ! Builder scratch: row union and child lists.
        call charge_count(self, ledger, int(n,i8)+2_i8*ns, 4_i8, SLU_PHASE_ROWS, status)
        if (.not. status_ok(status)) return
        allocate(rowlist(n), chead(ns), cnext(ns), stat=ierr)
        if (ierr /= 0) goto 910
        chead = 0
        cnext = 0
        mark(1:n) = 0
        do s = 1, ns
            self%sched%colnode(starts(s):ends(s)) = s
        end do
        self%factor_nnz = 0_i8
        self%maxrows = 0
        do s = 1, ns
            first = starts(s)
            last = ends(s)
            m = msize(s)
            w = last-first+1
            q = 0
            call supernode_gather_rows(first, last, sptr, sidx(1:len), s, mark, &
                rowlist, q)
            c = chead(s)
            do while (c /= 0)
                do p = self%panels(c)%width+1, self%panels(c)%m
                    r = self%panels(c)%rows(p)
                    if (mark(r) == s) cycle
                    mark(r) = s
                    q = q+1
                    rowlist(q) = r
                end do
                c = cnext(c)
            end do
            call supernode_sort_rows(rowlist, 1, q)
            if (q /= m) then
                call status_set(status, FORTSPARSE_INTERNAL_ERROR, &
                    'supernodal LU: row union disagrees with column counts')
                return
            end if
            call charge_count(self, ledger, int(m,i8), 4_i8, SLU_PHASE_ROWS, status)
            if (.not. status_ok(status)) return
            allocate(self%panels(s)%rows(m), stat=ierr)
            if (ierr /= 0) goto 910
            self%panels(s)%rows = rowlist(1:m)
            self%panels(s)%first = first
            self%panels(s)%width = w
            self%panels(s)%m = m
            self%factor_nnz = self%factor_nnz+2_i8*int(m,i8)*w-int(w,i8)*w
            self%maxrows = max(self%maxrows, m)
            if (m > w) then
                p = self%sched%colnode(rowlist(w+1))
                if (p <= s) then
                    call status_set(status, FORTSPARSE_INTERNAL_ERROR, &
                        'supernodal LU: non-forward supernode parent')
                    return
                end if
                cnext(s) = chead(p)
                chead(p) = s
            end if
        end do
        deallocate(rowlist, chead, cnext)
        call release(ledger, 4_i8*(int(n,i8)+2_i8*ns), status)
        return
910     call status_set(status, FORTSPARSE_OUT_OF_MEMORY, &
            'supernodal LU: analysis allocation failed')
    end subroutine build_nodes

    subroutine build_schedule(self, ledger, status)
        ! Update lists, subtree groups and the exact retained plan.
        class(slu_factor_t), intent(inout) :: self
        type(memory_budget_t), intent(inout) :: ledger
        type(fortsparse_status_t), intent(out) :: status
        integer :: s, first, group, dest, m, w, q, m2
        integer(i8) :: limit, update_size, panel_count, bytes
        real(dp) :: flops

        flops = 0.0_dp
        self%maxwork = 1
        panel_count = 0_i8
        do s = 1, self%node_count
            m = self%panels(s)%m
            w = self%panels(s)%width
            flops = flops+2.0_dp/3.0_dp*real(w,dp)**3+2.0_dp*real(m-w,dp)*real(w,dp)**2
            first = w+1
            do while (first <= m)
                call self%sched%update_group(self%panels(s)%supernode_t, first, group, &
                    dest)
                q = group-first
                m2 = m-first+1
                update_size = int(m2,i8)*q
                if (update_size > int(huge(m),i8)) then
                    call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                        'supernodal LU: update too large')
                    return
                end if
                self%maxwork = max(self%maxwork, int(update_size))
                flops = flops+2.0_dp*real(update_size,dp)*w &
                    +2.0_dp*real(m2-q,dp)*q*w
                first = group
            end do
            panel_count = panel_count+int(m,i8)*w+int(m-w,i8)*w
        end do
        self%flops = flops
        self%use_blas = flops >= 1.0e7_dp .or. any(self%panels%width >= 96)
        if (self%strict) self%use_blas = .false.
        ! Stage 4 executes serially; the schedule keeps the task groups.
        call self%sched%request_threads(1, flops, .true.)
        limit = ledger%limit
        if (self%inject == SLU_PHASE_INCOMING) then
            ! Test injection: the schedule's first reservation is refused.
            ledger%limit = max(1_i8, ledger%current)
            if (limit > 0_i8) ledger%limit = min(limit, ledger%limit)
        end if
        call self%sched%setup_incoming(self%panels%supernode_t, ledger, status)
        if (status_ok(status)) call self%sched%setup_groups(self%panels%supernode_t, &
            flops, ledger, status)
        if (self%inject == SLU_PHASE_INCOMING) then
            ledger%limit = limit
            self%inject = 0
        end if
        if (.not. status_ok(status)) then
            if (status%code == FORTSPARSE_RESOURCE_LIMIT) &
                self%failed_phase = SLU_PHASE_INCOMING
            return
        end if
        self%active_threads = self%sched%active_threads
        ! Exact retained layout (design section 7), before allocation.
        bytes = 8_i8*panel_count+4_i8*self%n
        bytes = bytes+self%structure_bytes()
        bytes = bytes+self%sched%memory_bytes()
        bytes = bytes+4_i8*int(self%n,i8)*self%sched%active_threads
        bytes = bytes+8_i8*int(self%maxwork,i8)*self%sched%active_threads
        bytes = bytes+8_i8*(int(self%n,i8)+int(self%maxrows,i8))
        self%plan_bytes = bytes
    end subroutine build_schedule

    integer(i8) function structure_bytes(self) result(bytes)
        ! Descriptors, row lists, failure flags and the 28n retained vectors.
        class(slu_factor_t), intent(in) :: self
        integer :: s
        bytes = 0_i8
        if (allocated(self%panels)) then
            bytes = int(storage_size(self%panels)/8,i8)*size(self%panels,kind=i8)
            do s = 1, size(self%panels)
                if (allocated(self%panels(s)%rows)) &
                    bytes = bytes+4_i8*size(self%panels(s)%rows,kind=i8)
            end do
        end if
        if (allocated(self%failures)) bytes = bytes+4_i8*size(self%failures,kind=i8)
        if (allocated(self%drp)) bytes = bytes+8_i8*size(self%drp,kind=i8)
        if (allocated(self%dcp)) bytes = bytes+8_i8*size(self%dcp,kind=i8)
        if (allocated(self%row_old)) bytes = bytes+4_i8*size(self%row_old,kind=i8)
        if (allocated(self%row_new)) bytes = bytes+4_i8*size(self%row_new,kind=i8)
        if (allocated(self%col_old)) bytes = bytes+4_i8*size(self%col_old,kind=i8)
    end function structure_bytes

    subroutine allocate_numeric(self, ledger, status)
        ! Phase panels (L and Ut values, block pivots) and phase scratch
        ! (row maps, update and solve work). Nothing is allocated later.
        class(slu_factor_t), intent(inout) :: self
        type(memory_budget_t), intent(inout) :: ledger
        type(fortsparse_status_t), intent(out) :: status
        integer :: s, m, w, ierr
        integer(i8) :: offset, count
        offset = 0_i8
        do s = 1, self%node_count
            m = self%panels(s)%m
            w = self%panels(s)%width
            self%panels(s)%offset_l = offset
            self%panels(s)%offset_u = offset+int(m,i8)*w
            offset = self%panels(s)%offset_u+int(m-w,i8)*w
        end do
        count = offset
        call charge_count(self, ledger, count, 8_i8, SLU_PHASE_PANELS, status)
        if (.not. status_ok(status)) return
        call charge_count(self, ledger, int(self%n,i8), 4_i8, SLU_PHASE_PANELS, status)
        if (.not. status_ok(status)) return
        allocate(self%values(count), self%pivots(self%n), stat=ierr)
        if (ierr /= 0) goto 910
        self%values = 0.0_dp
        self%pivots = 0
        if (self%inject == SLU_PHASE_SCRATCH) then
            call charge(self, ledger, 1_i8, SLU_PHASE_SCRATCH, status)
            return
        end if
        call self%sched%allocate_rowmap(ledger, status)
        if (status%code == FORTSPARSE_RESOURCE_LIMIT) self%failed_phase = SLU_PHASE_SCRATCH
        if (.not. status_ok(status)) return
        call self%sched%allocate_update_work(self%maxwork, ledger, status)
        if (status%code == FORTSPARSE_RESOURCE_LIMIT) self%failed_phase = SLU_PHASE_SCRATCH
        if (.not. status_ok(status)) return
        call charge_count(self, ledger, int(self%n,i8)+self%maxrows, 8_i8, &
            SLU_PHASE_SCRATCH, status)
        if (.not. status_ok(status)) return
        allocate(self%solve_work(self%n+self%maxrows), stat=ierr)
        if (ierr /= 0) goto 910
        self%solve_work = 0.0_dp
        return
910     call status_set(status, FORTSPARSE_OUT_OF_MEMORY, &
            'supernodal LU: panel allocation failed')
    end subroutine allocate_numeric

    subroutine finalize_labels(self)
        ! row_old held cmatch (matched column -> original row), row_new held
        ! rmatch and drp/dcp the original-label scales. Compose with col_old
        ! using the (still zero) row map and the solve scratch as temporaries.
        class(slu_factor_t), intent(inout) :: self
        integer :: k
        associate (tmp => self%sched%rowmap(:,1), work => self%solve_work)
            do k = 1, self%n
                tmp(k) = self%row_old(self%col_old(k))
            end do
            self%row_old = tmp(1:self%n)
            do k = 1, self%n
                self%row_new(self%row_old(k)) = k
            end do
            do k = 1, self%n
                work(k) = self%drp(self%row_old(k))
            end do
            self%drp = work(1:self%n)
            do k = 1, self%n
                work(k) = self%dcp(self%col_old(k))
            end do
            self%dcp = work(1:self%n)
            tmp = 0
            work = 0.0_dp
        end associate
    end subroutine finalize_labels

    logical function is_analyzed(self)
        class(slu_factor_t), intent(in) :: self
        is_analyzed = self%analyzed
    end function is_analyzed

    subroutine node_structure(self, s, first, width, m, rows)
        ! Diagnostics: supernode s columns first..first+width-1, rows(1:m).
        class(slu_factor_t), intent(in) :: self
        integer, intent(in) :: s
        integer, intent(out) :: first, width, m
        integer, allocatable, intent(out) :: rows(:)
        first = 0
        width = 0
        m = 0
        if (.not. self%analyzed .or. s < 1 .or. s > self%node_count) return
        first = self%panels(s)%first
        width = self%panels(s)%width
        m = self%panels(s)%m
        rows = self%panels(s)%rows
    end subroutine node_structure

    subroutine ordering(self, row_old, col_old, row_scale, col_scale)
        ! Diagnostics: M(k,l) = row_scale(k)*A(row_old(k),col_old(l))*col_scale(l).
        class(slu_factor_t), intent(in) :: self
        integer, allocatable, intent(out) :: row_old(:), col_old(:)
        real(dp), allocatable, intent(out) :: row_scale(:), col_scale(:)
        if (.not. self%analyzed) return
        row_old = self%row_old
        col_old = self%col_old
        row_scale = self%drp
        col_scale = self%dcp
    end subroutine ordering

    integer(i8) function memory_bytes(self) result(bytes)
        ! Owned array payload only; descriptors of the top-level object,
        ! allocator rounding and the OpenMP runtime are excluded.
        class(slu_factor_t), intent(in) :: self
        bytes = self%structure_bytes()+self%sched%memory_bytes()
        if (allocated(self%values)) bytes = bytes+8_i8*size(self%values,kind=i8)
        if (allocated(self%pivots)) bytes = bytes+4_i8*size(self%pivots,kind=i8)
        if (allocated(self%solve_work)) bytes = bytes+8_i8*size(self%solve_work,kind=i8)
    end function memory_bytes

    function plan(self) result(layout)
        ! Exact layout of the analyzed factor (owned arrays only).
        class(slu_factor_t), intent(in) :: self
        type(sparse_memory_plan_t) :: layout
        layout = sparse_memory_plan_t()
        if (.not. self%analyzed) return
        layout%factor_host_bytes = 8_i8*size(self%values,kind=i8) &
            +4_i8*size(self%pivots,kind=i8)
        layout%workspace_host_bytes = 4_i8*size(self%sched%rowmap,kind=i8) &
            +8_i8*size(self%sched%update_work,kind=i8) &
            +8_i8*size(self%solve_work,kind=i8)
        layout%persistent_host_bytes = self%plan_bytes-layout%factor_host_bytes &
            -layout%workspace_host_bytes
        layout%peak_host_bytes = self%peak_memory_bytes
        layout%symbolic_nnz = self%symbolic_nnz
        layout%max_rhs_batch = 1
        layout%prediction_quality = MEMORY_PLAN_EXACT_LAYOUT
        layout%fits = .true.
    end function plan

    subroutine slu_free(self, budget, status)
        class(slu_factor_t), intent(inout) :: self
        type(memory_budget_t), optional, intent(inout) :: budget
        type(fortsparse_status_t), optional, intent(out) :: status
        type(fortsparse_status_t) :: released
        call status_set(released, FORTSPARSE_OK, '')
        if (present(budget)) then
            call budget%release(self%owned_budget_bytes, released)
            if (present(status)) status = released
            if (.not. status_ok(released)) return
        end if
        if (present(status)) status = released
        if (allocated(self%panels)) deallocate(self%panels)
        if (allocated(self%pivots)) deallocate(self%pivots)
        if (allocated(self%failures)) deallocate(self%failures)
        if (allocated(self%row_old)) deallocate(self%row_old)
        if (allocated(self%row_new)) deallocate(self%row_new)
        if (allocated(self%col_old)) deallocate(self%col_old)
        if (allocated(self%values)) deallocate(self%values)
        if (allocated(self%drp)) deallocate(self%drp)
        if (allocated(self%dcp)) deallocate(self%dcp)
        if (allocated(self%solve_work)) deallocate(self%solve_work)
        call self%sched%free()
        self%n = 0
        self%node_count = 0
        self%active_threads = 1
        self%factor_nnz = 0_i8
        self%symbolic_nnz = 0_i8
        self%flops = 0.0_dp
        self%perturbed_pivots = 0
        self%block_swaps = 0
        self%max_growth = 0.0_dp
        self%min_pivot = 0.0_dp
        self%scaling_log10_range = 0.0_dp
        self%plan_bytes = 0_i8
        self%peak_memory_bytes = 0_i8
        self%failed_required_bytes = 0_i8
        self%failed_available_bytes = 0_i8
        self%failed_phase = 0
        self%opaque_ordering = .false.
        self%analysis_s = 0.0_dp
        self%factor_s = 0.0_dp
        self%analyzed = .false.
        self%factored = .false.
        self%strict = .false.
        self%use_blas = .true.
        self%inject = 0
        self%maxwork = 0
        self%maxrows = 0
        self%tau = 1.0e-13_dp
        self%tau_abs = 0.0_dp
        self%max_abs_m = 0.0_dp
        self%owned_budget_bytes = 0_i8
    end subroutine slu_free

    real(dp) function seconds()
        integer(i8) :: clock, rate
        call system_clock(clock, rate)
        seconds = real(clock,dp)/real(rate,dp)
    end function seconds
end module fortsparse_native_slu
