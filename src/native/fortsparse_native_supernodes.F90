module fortsparse_native_supernodes
    use fortsparse_kinds, only: dp, i8
    use fortsparse_memory, only: memory_budget_t, memory_checked_add, memory_checked_mul
    use fortsparse_csc, only: csc_t, csc_is_valid
    use fortsparse_status, only: fortsparse_status_t, status_set, FORTSPARSE_OK, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OUT_OF_MEMORY, &
        FORTSPARSE_NOT_FACTORED, FORTSPARSE_INDEX_OVERFLOW, &
        FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_INTERNAL_ERROR, &
        FORTSPARSE_NOT_POSITIVE_DEFINITE, FORTSPARSE_RESOURCE_LIMIT
    use fortsparse_spd_leaves, only: spd_panel_factor, spd_panel_solve, &
        spd_scatter_mapped_update, spd_scatter_update, spd_schur_update, &
        spd_divide, spd_gather_update
    use fortsparse_native_supernode_schedule, only: supernode_t, supernode_schedule_t, &
        supernode_partition, supernode_gather_rows, supernode_charge
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, &
        c_int, c_size_t, c_loc
#ifdef FORTSPARSE_HAVE_NATIVE_CUDA
    use fortsparse_native_cuda, only: spd_cuda_create_deferred_limited, spd_cuda_destroy, &
        spd_cuda_memory, spd_cuda_memory_stats, &
        spd_cuda_prepare_blas, &
        spd_cuda_upload, spd_cuda_download, spd_cuda_synchronize, &
        spd_cuda_panel_factor, spd_cuda_factor_info, spd_cuda_panel_solve, &
        spd_cuda_schur_product, spd_cuda_set_maps, spd_cuda_begin_factor, &
        spd_cuda_scatter_update, spd_cuda_schur_task, &
        spd_cuda_host_plan_bytes, spd_cuda_level_panel, spd_cuda_level_task, &
        spd_cuda_set_levels, spd_cuda_factor_level, spd_cuda_update_level, &
        spd_cuda_set_solve_plan, spd_cuda_solve
#endif
    !$ use omp_lib, only: omp_get_num_threads
    implicit none
    private
    integer, parameter :: NOT_SPD = FORTSPARSE_NOT_POSITIVE_DEFINITE

    ! The structural supernode (first, width, m, rows) comes from the shared
    ! schedule module; the numeric views below are SPD-owned.
    type, extends(supernode_t) :: panel_t
        integer :: update_first = 0, update_last = -1
        integer(c_size_t) :: offset = 0
        real(dp), pointer, contiguous :: values(:,:) => null()
    end type panel_t

    type :: update_t
        integer :: dest, first, q
        integer(c_size_t) :: row_map, column_map
    end type update_t

    type, public :: supernodes_t
        ! active_threads, task_groups, group_target_flops and
        ! budget_reduced_threads report the composed schedule's state.
        integer :: n = 0, node_count = 0, active_threads = 1, task_groups = 0
        integer(i8) :: factor_nnz = 0_i8
        logical :: budget_reduced_threads = .false.
        logical, private :: strict_budget = .false.
        integer(i8), private :: owned_budget_bytes = 0_i8
        type(memory_budget_t), private :: build_ledger
        integer(i8) :: device_bytes = 0_i8, device_peak_bytes = 0_i8
        integer(i8) :: device_required_bytes = 0_i8, device_available_bytes = 0_i8
        logical :: device_payload_tracked = .false.
        integer(i8), private :: device_limit = 0_i8
#ifdef FORTSPARSE_HAVE_NATIVE_CUDA
        integer(i8), private :: cuda_host_bytes = 0_i8
        integer, allocatable, private :: cuda_wide_ptr(:), cuda_wide_members(:)
#endif
        integer(i8) :: symbolic_nnz = 0_i8
        integer(i8) :: factor_storage_bytes = 0_i8, workspace_bytes = 0_i8
        real(dp) :: flops = 0.0_dp, numeric_s = 0.0_dp
        real(dp) :: group_target_flops = 0.0_dp
        logical, private :: ready = .false., use_blas = .true.
        type(panel_t), allocatable, private :: panels(:)
        type(supernode_schedule_t), private :: sched
        integer, allocatable, private :: failures(:)
        real(dp), allocatable, private :: solve_work(:)
        ! One owning host assembly allocation permits one factor upload. Panel
        ! views share it without maintaining a second full host factor copy.
        real(dp), pointer, contiguous, private :: host_values(:) => null()
        integer(c_size_t), private :: numeric_count = 0
        integer, private :: device = -1
        type(c_ptr), private :: cuda = c_null_ptr
        type(update_t), allocatable, private :: updates(:)
    contains
        procedure :: build
        procedure :: factor
        procedure :: solve
        procedure :: solve_block
        procedure :: free => supernodes_free
        procedure :: memory_bytes
        procedure :: has_device_context
        procedure :: finish_copy
        final :: finalize_supernodes
    end type supernodes_t

    interface
        subroutine dgemm(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import dp
            character, intent(in) :: ta, tb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            real(dp), intent(in) :: alpha, beta, a(lda,*), b(ldb,*)
            real(dp), intent(inout) :: c(ldc,*)
        end subroutine dgemm
        subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
            import dp
            character, intent(in) :: uplo, trans
            integer, intent(in) :: n, k, lda, ldc
            real(dp), intent(in) :: alpha, beta, a(lda,*)
            real(dp), intent(inout) :: c(ldc,*)
        end subroutine dsyrk
        subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            import dp
            character, intent(in) :: trans
            integer, intent(in) :: m, n, lda, incx, incy
            real(dp), intent(in) :: alpha, beta, a(lda,*), x(*)
            real(dp), intent(inout) :: y(*)
        end subroutine dgemv
        subroutine dtrsv(uplo, trans, diag, n, a, lda, x, incx)
            import dp
            character, intent(in) :: uplo, trans, diag
            integer, intent(in) :: n, lda, incx
            real(dp), intent(in) :: a(lda,*)
            real(dp), intent(inout) :: x(*)
        end subroutine dtrsv
    end interface
contains
    subroutine charge(self, count, element_bytes, status, phase)
        class(supernodes_t), intent(inout) :: self
        integer(i8), intent(in) :: count, element_bytes
        type(fortsparse_status_t), intent(out) :: status
        character(*), intent(in) :: phase
        call supernode_charge(self%build_ledger, count, element_bytes, status, phase)
    end subroutine charge

    subroutine publish_schedule(self)
        class(supernodes_t), intent(inout) :: self
        self%active_threads = self%sched%active_threads
        self%task_groups = self%sched%task_groups
        self%group_target_flops = self%sched%group_target_flops
        self%budget_reduced_threads = self%sched%budget_reduced_threads
    end subroutine publish_schedule

    subroutine build(self, n, lp, li, nthreads, status, device, budget, device_limit)
        ! The supplied ledger already includes caller-owned live allocations.
        ! Pass the same ledger to free, or reset the whole owning context ledger.
        ! No pointer to that caller ledger survives this call.
        ! Charge variable-size panel descriptors and all owned array payload;
        ! OpenMP task/runtime storage, stacks and allocator rounding are excluded.
        class(supernodes_t), intent(inout) :: self
        integer, intent(in) :: n, lp(:), li(:), nthreads
        type(fortsparse_status_t), intent(out) :: status
        integer, optional, intent(in) :: device
        type(memory_budget_t), optional, intent(inout) :: budget
        integer(i8), optional, intent(in) :: device_limit
        type(memory_budget_t) :: ledger
        integer(i8) :: baseline, device_peak, device_required, device_available
        logical :: device_tracked

        call self%free(budget, status)
        if (status%code /= FORTSPARSE_OK) return
        if (present(device_limit)) self%device_limit = device_limit
        if (self%device_limit < 0_i8) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, 'negative device payload cap')
            return
        end if
        self%strict_budget = present(budget)
        if (present(budget)) self%build_ledger = budget
        call self%build_ledger%reserve(0_i8, status)
        if (status%code /= FORTSPARSE_OK) return
        baseline = self%build_ledger%current
        if (self%strict_budget .and. present(device)) then
            if (device >= 0) then
                call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                    'supernode CPU budget does not cover CUDA allocations')
                return
            end if
        end if
        call build_storage(self, n, lp, li, nthreads, status, device)
        call publish_schedule(self)
        if (status%code == FORTSPARSE_OK .and. self%strict_budget) then
            self%owned_budget_bytes = self%memory_bytes()
            if (self%build_ledger%current-baseline /= self%owned_budget_bytes) &
                call status_set(status, FORTSPARSE_INTERNAL_ERROR, &
                    'supernode allocation ledger disagrees with retained layout')
        end if
        ledger = self%build_ledger
        if (status%code /= FORTSPARSE_OK) then
            ! All builder-local temporary arrays are already deallocated.
            ! Keep device refusal/high-water diagnostics after releasing payload.
            device_peak = self%device_peak_bytes
            device_required = self%device_required_bytes
            device_available = self%device_available_bytes
            device_tracked = self%device_payload_tracked
            call self%free()
            self%device_peak_bytes = device_peak
            self%device_required_bytes = device_required
            self%device_available_bytes = device_available
            self%device_payload_tracked = device_tracked
            ledger%current = baseline
        end if
        if (present(budget)) budget = ledger
    end subroutine build

    subroutine build_storage(self, n, lp, li, nthreads, status, device)
        class(supernodes_t), intent(inout) :: self
        integer, intent(in) :: n, lp(:), li(:), nthreads
        type(fortsparse_status_t), intent(out) :: status
        integer, optional, intent(in) :: device
        integer, allocatable :: mark(:), starts(:), ends(:), counts(:), rowlist(:)
        integer :: s, j, p, m, first, last, group, q, dest, ierr, maxwork, maxrows
        integer(i8) :: update_size, graph_bytes, bytes, solve_bytes

        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'invalid symbolic Cholesky pattern')
        if (present(device)) self%device = device
        if (self%device < -1) return
#ifndef FORTSPARSE_HAVE_NATIVE_CUDA
        if (self%device >= 0) then
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, 'native CUDA is not built')
            return
        end if
#endif
        if (n < 0) return
        if (n == huge(n) .or. size(li,kind=i8) >= int(huge(n),i8)) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                'supernode symbolic index range overflow')
            return
        end if
        if (size(lp) /= n+1) return
        if (lp(1) /= 1 .or. lp(n+1) /= size(li)+1) return
        do j = 1, n
            if (lp(j+1) <= lp(j)) return
            if (lp(j) < 1 .or. lp(j+1) > size(li)+1) return
            if (li(lp(j)) /= j) return
            do p = lp(j)+1, lp(j+1)-1
                if (li(p) <= li(p-1) .or. li(p) > n) return
            end do
        end do
        self%n = n
        self%symbolic_nnz = size(li,kind=i8)
        graph_bytes = 0_i8
        call memory_checked_mul(int(n,i8), 20_i8, graph_bytes, status)
        if (status%code /= FORTSPARSE_OK) return
        call self%build_ledger%reserve(graph_bytes, status, 'panel graph scratch')
        if (status%code /= FORTSPARSE_OK) return
        allocate(mark(n), starts(n), ends(n), counts(n), rowlist(n), stat=ierr)
        if (ierr /= 0) then
            call allocation_error(self, status)
            return
        end if
        ! The column structure gives the elimination-tree parent (first
        ! off-diagonal row) and the column count. mark and rowlist hold them
        ! only until the partition is formed; the scratch charge is unchanged.
        do j = 1, n
            rowlist(j) = lp(j+1)-lp(j)-1
            mark(j) = 0
            if (rowlist(j) > 0) mark(j) = li(lp(j)+1)
        end do
        ! Wider dense panels pay only when the structural fill can amortize
        ! packing and diagonal storage; strict budgets keep the compact limit.
        call supernode_partition(n, mark, rowlist, self%strict_budget, starts, ends, &
            counts, self%node_count)
        call charge(self, int(self%node_count,i8), &
            int(storage_size(panel_t())/8,i8)+8_i8, status, 'panel descriptors/tokens')
        if (status%code /= FORTSPARSE_OK) return
        call charge(self, int(n,i8), 4_i8, status, 'column ownership')
        if (status%code /= FORTSPARSE_OK) return
        allocate(self%panels(self%node_count), self%sched%colnode(n), &
            self%sched%tokens(self%node_count), self%failures(self%node_count), stat=ierr)
        if (ierr /= 0) then
            call allocation_error(self, status)
            return
        end if
        self%sched%n = n
        self%sched%node_count = self%node_count
        mark = 0
        self%sched%tokens = 0
        self%failures = 0
        maxrows = 0
        do s = 1, self%node_count
            first = starts(s)
            last = ends(s)
            m = counts(s)
            self%sched%colnode(first:last) = s
            q = 0
            call supernode_gather_rows(first, last, lp, li, s, mark, rowlist, q)
            ! The parent/count rule assumes a closed elimination structure.
            if (q /= m) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    'invalid symbolic Cholesky pattern')
                return
            end if
            self%panels(s)%first = first
            self%panels(s)%width = last-first+1
            self%panels(s)%m = m
            call charge(self, int(m,i8), 4_i8, status, 'panel row indices')
            if (status%code /= FORTSPARSE_OK) return
            allocate(self%panels(s)%rows(m), stat=ierr)
            if (ierr /= 0) then
                call allocation_error(self, status)
                return
            end if
            self%panels(s)%rows = rowlist(1:m)
            self%factor_nnz = self%factor_nnz + int(m,i8)*(last-first+1) &
                -int(last-first+1,i8)*(last-first)/2
            self%panels(s)%offset = self%numeric_count
            bytes = 0_i8
            call memory_checked_mul(int(m,i8), int(last-first+1,i8), bytes, status)
            if (status%code /= FORTSPARSE_OK) return
            update_size = 0_i8
            call memory_checked_add(int(self%numeric_count,i8), bytes, &
                update_size, status)
            if (status%code /= FORTSPARSE_OK) return
            self%numeric_count = int(update_size,c_size_t)
            self%flops = self%flops+real(m,dp)*(last-first+1)**2
            maxrows = max(maxrows, m)
        end do
        call memory_checked_mul(int(self%numeric_count,i8), 8_i8, &
            self%factor_storage_bytes, status)
        if (status%code /= FORTSPARSE_OK) return
        call self%build_ledger%reserve(self%factor_storage_bytes, status, &
            'panel values')
        if (status%code /= FORTSPARSE_OK) return
        allocate(self%host_values(self%numeric_count), stat=ierr)
        if (ierr /= 0) then
            call allocation_error(self, status)
            return
        end if
        do s = 1, self%node_count
            m = self%panels(s)%m
            q = self%panels(s)%width
            update_size = self%panels(s)%offset
            self%panels(s)%values(1:m,1:q) => &
                self%host_values(update_size+1:update_size+int(m,i8)*q)
        end do
        maxwork = 1
        do s = 1, self%node_count
            first = self%panels(s)%width+1
            m = self%panels(s)%m
            do while (first <= m)
                call self%sched%update_group(self%panels(s)%supernode_t, first, group, dest)
                q = group-first
                update_size = int(m-first+1,i8)*q
                if (update_size > huge(maxwork)) then
                    call status_set(status, FORTSPARSE_INDEX_OVERFLOW, 'supernode update too large')
                    return
                end if
                maxwork = max(maxwork, int(update_size))
                self%flops = self%flops+2.0_dp*real(update_size,dp)*self%panels(s)%width
                first = group
            end do
        end do
        if (self%strict_budget) then
            deallocate(mark, starts, ends, counts, rowlist)
            call self%build_ledger%release(graph_bytes, status)
            if (status%code /= FORTSPARSE_OK) return
        end if
        ! Below this work budget, generated contractions avoid the first
        ! BLAS workspace/code-page cost. Wide panels retain their BLAS TRSM
        ! schedule, and device factors already pay runtime initialization.
        self%use_blas = self%flops >= 1.0e7_dp .or. self%device >= 0
        if (any(self%panels%width >= 96)) self%use_blas = .true.
        if (self%strict_budget) self%use_blas = .false.
        call self%sched%request_threads(nthreads, self%flops, self%device >= 0)
        if (self%device < 0) then
            call self%sched%setup_incoming(self%panels%supernode_t, self%build_ledger, status)
            if (status%code == FORTSPARSE_OK) call self%sched%setup_groups( &
                self%panels%supernode_t, self%flops, self%build_ledger, status)
            if (status%code /= FORTSPARSE_OK) then
                return
            end if
        end if
        if (self%strict_budget) then
            solve_bytes = 0_i8
            call memory_checked_mul(int(maxrows,i8), 8_i8, solve_bytes, status)
            if (status%code /= FORTSPARSE_OK) return
            call self%sched%fit_thread_scratch(maxwork, solve_bytes, self%build_ledger, status)
            if (status%code /= FORTSPARSE_OK) return
        end if
        call self%sched%allocate_rowmap(self%build_ledger, status)
        if (status%code /= FORTSPARSE_OK) return
        call charge(self, int(maxrows,i8), 8_i8, status, 'triangular solve scratch')
        if (status%code /= FORTSPARSE_OK) return
        allocate(self%solve_work(maxrows), stat=ierr)
        if (ierr /= 0) then
            call allocation_error(self, status)
            return
        end if
#ifdef FORTSPARSE_HAVE_NATIVE_CUDA
        if (self%device >= 0) then
            call setup_cuda(self, maxwork, status)
            self%workspace_bytes = self%memory_bytes()-self%factor_storage_bytes
            return
        end if
#endif
        call self%sched%allocate_update_work(maxwork, self%build_ledger, status)
        if (status%code /= FORTSPARSE_OK) return
        self%workspace_bytes = self%memory_bytes()-self%factor_storage_bytes
        if (allocated(mark)) then
            deallocate(mark, starts, ends, counts, rowlist)
            call self%build_ledger%release(graph_bytes, status)
            if (status%code /= FORTSPARSE_OK) return
        end if
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine build_storage

    subroutine factor(self, a, status, perm, inverse)
        class(supernodes_t), target, intent(inout) :: self
        type(csc_t), intent(in) :: a
        type(fortsparse_status_t), intent(out) :: status
        integer, optional, intent(in) :: perm(:), inverse(:)
        integer, pointer :: token(:)
#ifdef FORTSPARSE_HAVE_OMP_DEPEND_ITERATOR
        integer, pointer :: sources(:)
        integer :: incoming_start, incoming_end
#endif
        integer :: s, j, p, row, position, original_column, mapped_row
#ifndef FORTSPARSE_HAVE_OMP_DEPEND_ITERATOR
        integer :: group, first, q, dest
#endif
        real(dp) :: started

        self%ready = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'incompatible supernode matrix')
        if (.not. allocated(self%panels)) return
        if (.not. csc_is_valid(a)) return
        if (a%nrow /= self%n .or. a%ncol /= self%n) return
        if (.not. all(ieee_is_finite(a%val))) return
        if (present(perm) .neqv. present(inverse)) return
        if (present(perm)) then
            if (size(perm) /= self%n .or. size(inverse) /= self%n) return
            if (any(perm < 1) .or. any(perm > self%n)) return
            if (any(inverse < 1) .or. any(inverse > self%n)) return
        end if
        started = seconds()
        self%failures = 0
        self%sched%rowmap = 0
        do s = 1, self%node_count
            self%panels(s)%values = 0.0_dp
            do p = 1, self%panels(s)%m
                self%sched%rowmap(self%panels(s)%rows(p),1) = p
            end do
            do j = 1, self%panels(s)%width
                row = self%panels(s)%first+j-1
                original_column = row
                if (present(perm)) original_column = perm(row)
                do p = a%col_ptr(original_column), a%col_ptr(original_column+1)-1
                    mapped_row = a%row_idx(p)
                    if (present(inverse)) mapped_row = inverse(mapped_row)
                    if (mapped_row < row) cycle
                    position = self%sched%rowmap(mapped_row,1)
                    if (position == 0) return
                    self%panels(s)%values(position,j) = a%val(p)
                end do
            end do
            do p = 1, self%panels(s)%m
                self%sched%rowmap(self%panels(s)%rows(p),1) = 0
            end do
        end do
#ifdef FORTSPARSE_HAVE_NATIVE_CUDA
        if (self%device >= 0) then
            call factor_cuda(self, status)
            self%numeric_s = seconds()-started
            return
        end if
#endif
        if (self%sched%active_threads == 1) then
            call factor_serial(self)
        else
#ifdef FORTSPARSE_HAVE_OMP_DEPEND_ITERATOR
            token => self%sched%group_tokens
            sources => self%sched%group_sources
            !$omp parallel num_threads(self%sched%active_threads) default(shared)
            !$omp single
            !$ self%sched%active_threads = omp_get_num_threads()
            do s = 1, self%sched%task_groups
                incoming_start = self%sched%group_ptr(s)
                ! OpenMP iterator end is exclusive, including in Fortran.
                incoming_end = self%sched%group_ptr(s+1)
                !$omp task firstprivate(s) &
                !$omp depend(iterator(incoming=incoming_start:incoming_end), &
                !$omp in:token(sources(incoming))) depend(out:token(s))
                call factor_group(self, s)
                !$omp end task
            end do
            !$omp end single
            !$omp end parallel
#else
            token => self%sched%tokens
            !$omp parallel num_threads(self%sched%active_threads) default(shared)
            !$omp single
            !$ self%sched%active_threads = omp_get_num_threads()
            do s = 1, self%node_count
                !$omp task firstprivate(s) depend(inout:token(s))
                call factor_panel(self, s)
                !$omp end task
                first = self%panels(s)%width+1
                do while (first <= self%panels(s)%m)
                    call self%sched%update_group(self%panels(s)%supernode_t, &
                        first, group, dest)
                    q = group-first
                    !$omp task firstprivate(s,dest,first,q) &
                    !$omp depend(in:token(s)) depend(inout:token(dest))
                    call update_panel(self, s, dest, first, q)
                    !$omp end task
                    first = group
                end do
            end do
            !$omp end single
            !$omp end parallel
#endif
            self%active_threads = self%sched%active_threads
        end if
        self%numeric_s = seconds()-started
        if (any(self%failures /= 0)) then
            call status_set(status, NOT_SPD, 'supernode: nonpositive or nonfinite pivot')
            return
        end if
        self%ready = .true.
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine factor

    subroutine factor_serial(self)
        class(supernodes_t), intent(inout) :: self
        integer :: s
        do s = 1, self%node_count
            call factor_destination(self, s)
        end do
    end subroutine factor_serial

#ifdef FORTSPARSE_HAVE_OMP_DEPEND_ITERATOR
    subroutine factor_group(self, group)
        class(supernodes_t), intent(inout) :: self
        integer, intent(in) :: group
        integer :: member
        do member = self%sched%group_member_ptr(group), self%sched%group_member_ptr(group+1)-1
            call factor_destination(self, self%sched%group_members(member))
        end do
    end subroutine factor_group
#endif

    subroutine factor_destination(self, dest)
        class(supernodes_t), intent(inout) :: self
        integer, intent(in) :: dest
        integer :: p, thread
        thread = self%sched%workspace_thread()
        do p = 1, self%panels(dest)%m
            self%sched%rowmap(self%panels(dest)%rows(p),thread) = p
        end do
        do p = self%sched%incoming_ptr(dest), self%sched%incoming_ptr(dest+1)-1
            call accumulate_update(self, self%sched%incoming_source(p), dest, &
                self%sched%incoming_first(p), self%sched%incoming_q(p), thread)
        end do
        call factor_panel(self, dest)
        do p = 1, self%panels(dest)%m
            self%sched%rowmap(self%panels(dest)%rows(p),thread) = 0
        end do
    end subroutine factor_destination

    subroutine factor_panel(self, s)
        class(supernodes_t), intent(inout) :: self
        integer, intent(in) :: s
        integer :: info, w, m
        real(dp), pointer, contiguous :: values(:,:), tail(:)
        if (self%failures(s) /= 0) return
        w = self%panels(s)%width
        m = self%panels(s)%m
        values => self%panels(s)%values
        tail => self%host_values(self%panels(s)%offset+w+1:)
        call spd_panel_factor(values, w, m, info)
        if (info /= 0) then
            self%failures(s) = NOT_SPD
            return
        end if
        if (m > w) call spd_panel_solve(values, w, m, tail, m-w, m, self%use_blas)
        if (.not. all(ieee_is_finite(self%panels(s)%values))) self%failures(s) = NOT_SPD
    end subroutine factor_panel

#ifndef FORTSPARSE_HAVE_OMP_DEPEND_ITERATOR
    subroutine update_panel(self, source, dest, first, q)
        class(supernodes_t), intent(inout) :: self
        integer, intent(in) :: source, dest, first, q
        integer :: thread, p
        thread = self%sched%workspace_thread()
        do p = 1, self%panels(dest)%m
            self%sched%rowmap(self%panels(dest)%rows(p),thread) = p
        end do
        call accumulate_update(self, source, dest, first, q, thread)
        do p = 1, self%panels(dest)%m
            self%sched%rowmap(self%panels(dest)%rows(p),thread) = 0
        end do
    end subroutine update_panel
#endif

    subroutine accumulate_update(self, source, dest, first, q, thread)
        class(supernodes_t), intent(inout) :: self
        integer, intent(in) :: source, dest, first, q, thread
        integer :: m, w, m2, j, p, column, offset, position, dest_m
        integer(c_size_t) :: dest_offset
        logical :: contiguous_rows
        real(dp), pointer, contiguous :: source_values(:), dest_values(:)
        if (self%failures(source) /= 0) self%failures(dest) = self%failures(source)
        if (self%failures(dest) /= 0) return
        m = self%panels(source)%m
        w = self%panels(source)%width
        m2 = m-first+1
        source_values => self%host_values(self%panels(source)%offset+first:)
        ! A rank-one contraction is already the retained generated scalar
        ! Schur update. Apply it directly through the destination row map;
        ! no BLAS dispatch, product scratch or dense-intersection scan.
        if (w == 1) then
            do j = 1, q
                column = self%panels(source)%rows(first+j-1)-self%panels(dest)%first+1
                dest_values => self%panels(dest)%values(:,column)
                call spd_scatter_mapped_update(dest_values, &
                    self%panels(source)%rows(first+j-1:), self%sched%rowmap(:,thread), &
                    source_values(j:), m2-j+1, source_values(j))
            end do
            return
        end if
        if (.not. self%use_blas) then
            self%sched%update_work(:m2*q,thread) = 0.0_dp
            call spd_schur_update(self%sched%update_work(1,thread), m2, &
                self%panels(source)%values, m, self%panels(source)%values, m, &
                m2, q, w, first)
        else
            ! Exact dense intersections can update their owning panel directly.
            ! SYRK touches only the diagonal block's lower triangle; GEMM handles
            ! its rectangular tail. Irregular intersections retain mapped scatter.
            position = self%sched%rowmap(self%panels(source)%rows(first),thread)
            contiguous_rows = position > 0
            if (contiguous_rows) then
                do p = first+1, m
                    if (self%sched%rowmap(self%panels(source)%rows(p),thread) &
                        == position+p-first) cycle
                    contiguous_rows = .false.
                    exit
                end do
            end if
            if (contiguous_rows) then
                do j = 2, q
                    if (self%panels(source)%rows(first+j-1) &
                        == self%panels(source)%rows(first)+j-1) cycle
                    contiguous_rows = .false.
                    exit
                end do
            end if
            if (contiguous_rows) then
                column = self%panels(source)%rows(first)-self%panels(dest)%first+1
                dest_m = self%panels(dest)%m
                dest_offset = self%panels(dest)%offset &
                    +int(column-1,c_size_t)*dest_m+position
                dest_values => self%host_values(dest_offset:)
                call dsyrk('L', 'N', q, w, -1.0_dp, source_values, m, &
                    1.0_dp, dest_values, dest_m)
                if (m2 > q) then
                    dest_values => self%host_values(dest_offset+q:)
                    call dgemm('N', 'T', m2-q, q, w, -1.0_dp, &
                        source_values(q+1:), m, source_values, m, &
                        1.0_dp, dest_values, dest_m)
                end if
                return
            end if
            call dgemm('N', 'T', m2, q, w, -1.0_dp, &
                source_values, m, source_values, m, &
                0.0_dp, self%sched%update_work(1,thread), m2)
        end if
        do j = 1, q
            column = self%panels(source)%rows(first+j-1)-self%panels(dest)%first+1
            offset = (j-1)*m2+j
            dest_values => self%panels(dest)%values(:,column)
            call spd_scatter_mapped_update(dest_values, &
                self%panels(source)%rows(first+j-1:), self%sched%rowmap(:,thread), &
                self%sched%update_work(offset:,thread), m2-j+1, -1.0_dp)
        end do
    end subroutine accumulate_update

    subroutine solve(self, b, status)
        class(supernodes_t), intent(inout) :: self
#ifdef FORTSPARSE_HAVE_NATIVE_CUDA
        real(dp), target, contiguous, intent(inout) :: b(:)
#else
        real(dp), contiguous, intent(inout) :: b(:)
#endif
        type(fortsparse_status_t), intent(out) :: status
        integer :: s, m, w, first, p, j, column
        real(dp) :: value
        real(dp), pointer, contiguous :: values(:,:), tail(:)
        call status_set(status, FORTSPARSE_NOT_FACTORED, 'supernode factor unavailable')
        if (.not. self%ready) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'incompatible supernode RHS')
        if (size(b) /= self%n) return
#ifdef FORTSPARSE_HAVE_NATIVE_CUDA
        if (self%device >= 0) then
            call solve_cuda_tile(self, c_loc(b), 1, status)
            return
        end if
#endif
        if (.not. self%use_blas) then
            do s = 1, self%node_count
                m = self%panels(s)%m
                values => self%panels(s)%values
                do j = 1, self%panels(s)%width
                    column = self%panels(s)%first+j-1
                    b(column) = spd_divide(b(column), values(j,j))
                    call spd_scatter_update(b, self%panels(s)%rows(j+1:), &
                        values(j+1:,j), m-j, b(column))
                end do
            end do
            do s = self%node_count, 1, -1
                m = self%panels(s)%m
                values => self%panels(s)%values
                do j = self%panels(s)%width, 1, -1
                    column = self%panels(s)%first+j-1
                    value = b(column)
                    call spd_gather_update(value, self%panels(s)%rows(j+1:), &
                        values(j+1:,j), b, m-j)
                    b(column) = spd_divide(value, values(j,j))
                end do
            end do
            call status_set(status, FORTSPARSE_OK, '')
            return
        end if
        do s = 1, self%node_count
            m = self%panels(s)%m
            w = self%panels(s)%width
            first = self%panels(s)%first
            values => self%panels(s)%values
            tail => self%host_values(self%panels(s)%offset+w+1:)
            call dtrsv('L', 'N', 'N', w, values, m, b(first:), 1)
            if (m == w) cycle
            call dgemv('N', m-w, w, 1.0_dp, &
                tail, m, &
                b(first:), 1, 0.0_dp, self%solve_work, 1)
            call spd_scatter_update(b, self%panels(s)%rows(w+1:), self%solve_work, m-w, 1.0_dp)
        end do
        do s = self%node_count, 1, -1
            m = self%panels(s)%m
            w = self%panels(s)%width
            first = self%panels(s)%first
            values => self%panels(s)%values
            tail => self%host_values(self%panels(s)%offset+w+1:)
            if (m > w) then
                do p = w+1, m
                    self%solve_work(p-w) = b(self%panels(s)%rows(p))
                end do
                call dgemv('T', m-w, w, -1.0_dp, &
                    tail, m, &
                    self%solve_work, 1, 1.0_dp, b(first:), 1)
            end if
            call dtrsv('L', 'T', 'N', w, values, m, b(first:), 1)
        end do
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine solve

    subroutine solve_block(self, b, status)
        class(supernodes_t), intent(inout) :: self
        real(dp), target, contiguous, intent(inout) :: b(:,:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: j
#ifdef FORTSPARSE_HAVE_NATIVE_CUDA
        integer :: count
#endif
        call status_set(status, FORTSPARSE_NOT_FACTORED, 'supernode factor unavailable')
        if (.not. self%ready) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'incompatible supernode RHS')
        if (size(b,1) /= self%n) return
        call status_set(status, FORTSPARSE_OK, '')
#ifdef FORTSPARSE_HAVE_NATIVE_CUDA
        if (self%device >= 0) then
            do j = 1, size(b,2), 8
                count = min(8,size(b,2)-j+1)
                call solve_cuda_tile(self, c_loc(b(1,j)), count, status)
                if (status%code /= FORTSPARSE_OK) return
            end do
            return
        end if
#endif
        do j = 1, size(b,2)
            call self%solve(b(:,j), status)
            if (status%code /= FORTSPARSE_OK) return
        end do
    end subroutine solve_block

    subroutine supernodes_free(self, budget, status)
        class(supernodes_t), intent(inout) :: self
        type(memory_budget_t), optional, intent(inout) :: budget
        type(fortsparse_status_t), optional, intent(out) :: status
        type(fortsparse_status_t) :: released
#ifdef FORTSPARSE_HAVE_NATIVE_CUDA
        integer(c_int) :: code
#endif
        call status_set(released, FORTSPARSE_OK, '')
        if (present(budget)) then
            call budget%release(self%owned_budget_bytes, released)
            if (present(status)) status = released
            if (released%code /= FORTSPARSE_OK) return
        end if
        if (present(status)) status = released
#ifdef FORTSPARSE_HAVE_NATIVE_CUDA
        if (c_associated(self%cuda)) code = spd_cuda_destroy(self%cuda)
        if (allocated(self%cuda_wide_ptr)) deallocate(self%cuda_wide_ptr)
        if (allocated(self%cuda_wide_members)) deallocate(self%cuda_wide_members)
#endif
        self%cuda = c_null_ptr
        if (allocated(self%panels)) deallocate(self%panels)
        if (associated(self%host_values)) deallocate(self%host_values)
        if (allocated(self%updates)) deallocate(self%updates)
        call self%sched%free()
        if (allocated(self%failures)) deallocate(self%failures)
        if (allocated(self%solve_work)) deallocate(self%solve_work)
        self%n = 0
        self%node_count = 0
        self%active_threads = 1
        self%task_groups = 0
        self%group_target_flops = 0.0_dp
        self%factor_nnz = 0_i8
        self%symbolic_nnz = 0_i8
        self%factor_storage_bytes = 0_i8
        self%workspace_bytes = 0_i8
        self%numeric_count = 0_c_size_t
        self%device_bytes = 0_i8
        self%device_peak_bytes = 0_i8
        self%device_required_bytes = 0_i8
        self%device_available_bytes = 0_i8
        self%device_payload_tracked = .false.
        self%device_limit = 0_i8
#ifdef FORTSPARSE_HAVE_NATIVE_CUDA
        self%cuda_host_bytes = 0_i8
#endif
        self%device = -1
        self%flops = 0.0_dp
        self%use_blas = .true.
        self%numeric_s = 0.0_dp
        self%ready = .false.
        self%strict_budget = .false.
        self%budget_reduced_threads = .false.
        self%owned_budget_bytes = 0_i8
        self%build_ledger = memory_budget_t()
    end subroutine supernodes_free

    logical function has_device_context(self)
        class(supernodes_t), intent(in) :: self
        has_device_context = c_associated(self%cuda)
    end function has_device_context

    subroutine finish_copy(self, src, status)
        ! Call immediately after the intrinsic assignment self = src. It drops
        ! the aliased host/device ownership before any fallible step, then
        ! gives self its own host panel storage. Device contexts are never
        ! copied: self is then released and reported unfactored.
        class(supernodes_t), intent(inout) :: self
        type(supernodes_t), intent(in) :: src
        type(fortsparse_status_t), intent(out) :: status
        integer :: s, ierr
        integer(c_size_t) :: offset

        nullify(self%host_values)
        self%cuda = c_null_ptr
        if (allocated(self%panels)) then
            do s = 1, size(self%panels)
                nullify(self%panels(s)%values)
            end do
        end if
        call status_set(status, FORTSPARSE_OK, '')
        if (c_associated(src%cuda)) then
            call self%free()
            call status_set(status, FORTSPARSE_NOT_FACTORED, &
                'native copy: CUDA device factors are not copied; ' // &
                'destination is unfactored')
            return
        end if
        if (.not. associated(src%host_values)) return
        allocate(self%host_values(size(src%host_values, kind=i8)), stat=ierr)
        if (ierr /= 0) then
            call self%free()
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, &
                'native copy: panel value allocation failed')
            return
        end if
        self%host_values = src%host_values
        do s = 1, self%node_count
            if (.not. associated(src%panels(s)%values)) cycle
            offset = self%panels(s)%offset
            self%panels(s)%values(1:self%panels(s)%m, 1:self%panels(s)%width) => &
                self%host_values(offset+1:offset+int(self%panels(s)%m,i8)* &
                self%panels(s)%width)
        end do
    end subroutine finish_copy

    subroutine finalize_supernodes(self)
        type(supernodes_t), intent(inout) :: self
        call self%free()
    end subroutine finalize_supernodes

    integer(i8) function memory_bytes(self) result(bytes)
        class(supernodes_t), intent(in) :: self
        integer :: s
        bytes = 0_i8
#ifdef FORTSPARSE_HAVE_NATIVE_CUDA
        bytes = self%cuda_host_bytes
        if (allocated(self%cuda_wide_ptr)) bytes = bytes+4_i8*size(self%cuda_wide_ptr,kind=i8)
        if (allocated(self%cuda_wide_members)) bytes = bytes+4_i8*size(self%cuda_wide_members,kind=i8)
#endif
        if (.not. allocated(self%panels)) return
        bytes = bytes+int(storage_size(self%panels)/8,i8)*size(self%panels,kind=i8)
        do s = 1, self%node_count
            if (allocated(self%panels(s)%rows)) &
                bytes = bytes+4_i8*size(self%panels(s)%rows,kind=i8)
            if (associated(self%panels(s)%values)) &
                bytes = bytes+8_i8*size(self%panels(s)%values,kind=i8)
        end do
        bytes = bytes+self%sched%memory_bytes()
        if (allocated(self%failures)) bytes = bytes+4_i8*size(self%failures,kind=i8)
        if (allocated(self%solve_work)) bytes = bytes+8_i8*size(self%solve_work,kind=i8)
        if (allocated(self%updates)) &
            bytes = bytes+int(storage_size(self%updates)/8,i8)*size(self%updates,kind=i8)
    end function memory_bytes

#ifdef FORTSPARSE_HAVE_NATIVE_CUDA
    subroutine setup_cuda(self, maxwork, status)
        class(supernodes_t), intent(inout) :: self
        integer, intent(in) :: maxwork
        type(fortsparse_status_t), intent(out) :: status
        integer(c_int), allocatable, target :: maps(:), tail_rows(:)
        type(spd_cuda_level_task), allocatable, target :: tasks(:)
        type(spd_cuda_level_panel), allocatable, target :: panels(:)
        integer, allocatable :: levels(:), next(:)
        integer :: level_count, l, count_wide
        integer :: s, first, group, dest, q, m2, u, v, count_updates, count_thin, base, p, ierr
        integer(i8) :: count_maps, count_tail
        integer(c_int) :: code, queried
        type(spd_cuda_memory) :: report
        integer(i8) :: arena_count, rhs_count
        logical :: needs_blas

        needs_blas = .false.
        count_updates = 0
        count_thin = 0
        count_maps = 0_i8
        count_tail = 0_i8
        do s = 1, self%node_count
            count_tail = count_tail+int(self%panels(s)%m-self%panels(s)%width,i8)
            if (self%panels(s)%width > 8 .and. &
                self%panels(s)%m > self%panels(s)%width) needs_blas = .true.
            first = self%panels(s)%width+1
            do while (first <= self%panels(s)%m)
                call self%sched%update_group(self%panels(s)%supernode_t, first, group, dest)
                count_updates = count_updates+1
                if (self%panels(s)%width <= 8) count_thin = count_thin+1
                count_maps = count_maps+self%panels(s)%m-first+1+group-first
                first = group
            end do
        end do
        if (count_maps > huge(base) .or. count_tail > huge(base)) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, 'CUDA scatter maps too large')
            return
        end if
        allocate(self%updates(count_updates), maps(int(count_maps)), &
            tasks(count_thin), panels(self%node_count), levels(self%node_count), &
            tail_rows(int(count_tail)), stat=ierr)
        if (ierr /= 0) then
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, 'CUDA scatter map allocation failed')
            return
        end if
        base = 0
        u = 0
        v = 0
        levels = 0
        do s = 1, self%node_count
            self%panels(s)%update_first = u+1
            first = self%panels(s)%width+1
            do while (first <= self%panels(s)%m)
                call self%sched%update_group(self%panels(s)%supernode_t, first, group, dest)
                q = group-first
                m2 = self%panels(s)%m-first+1
                if (dest <= s .or. dest > self%node_count) then
                    call status_set(status, FORTSPARSE_INTERNAL_ERROR, 'CUDA invalid update DAG edge')
                    return
                end if
                levels(dest) = max(levels(dest),levels(s)+1)
                u = u+1
                self%updates(u)%dest = dest
                self%updates(u)%first = first
                self%updates(u)%q = q
                self%updates(u)%row_map = int(base,c_size_t)
                do p = 1, self%panels(dest)%m
                    self%sched%rowmap(self%panels(dest)%rows(p),1) = p
                end do
                do p = first, self%panels(s)%m
                    maps(base+p-first+1) = self%sched%rowmap(self%panels(s)%rows(p),1)
                end do
                base = base+m2
                self%updates(u)%column_map = int(base,c_size_t)
                maps(base+1:base+q) = self%panels(s)%rows(first:group-1) &
                    -self%panels(dest)%first+1
                base = base+q
                if (self%panels(s)%width <= 8) then
                    v = v+1
                    tasks(v)%source = int(s-1,c_size_t)
                    tasks(v)%destination = int(dest-1,c_size_t)
                    tasks(v)%task = spd_cuda_schur_task( &
                        self%panels(s)%offset+int(first-1,c_size_t), &
                        self%panels(s)%offset+int(first-1,c_size_t), &
                        self%panels(dest)%offset, self%updates(u)%row_map, &
                        self%updates(u)%column_map, int(self%panels(s)%m,c_int), &
                        int(self%panels(s)%m,c_int), int(m2,c_int), int(q,c_int), &
                        int(self%panels(s)%width,c_int), int(self%panels(dest)%m,c_int), &
                        int(self%panels(dest)%width,c_int))
                end if
                do p = 1, self%panels(dest)%m
                    self%sched%rowmap(self%panels(dest)%rows(p),1) = 0
                end do
                first = group
            end do
            self%panels(s)%update_last = u
        end do
        level_count = 0
        if (self%node_count > 0) level_count = maxval(levels)+1
        count_wide = 0
        do s = 1, self%node_count
            panels(s) = spd_cuda_level_panel(self%panels(s)%offset, &
                int(self%panels(s)%width,c_int), int(self%panels(s)%m,c_int), &
                int(self%panels(s)%m,c_int), int(levels(s),c_int))
            if (self%panels(s)%width > 8) count_wide = count_wide+1
        end do
        allocate(self%cuda_wide_ptr(level_count+1), self%cuda_wide_members(count_wide), &
            next(level_count), stat=ierr)
        if (ierr /= 0) then
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, 'CUDA level allocation failed')
            return
        end if
        self%cuda_wide_ptr = 0
        do s = 1, self%node_count
            if (self%panels(s)%width > 8) &
                self%cuda_wide_ptr(levels(s)+2) = self%cuda_wide_ptr(levels(s)+2)+1
        end do
        self%cuda_wide_ptr(1) = 1
        do l = 1, level_count
            self%cuda_wide_ptr(l+1) = self%cuda_wide_ptr(l+1)+self%cuda_wide_ptr(l)
        end do
        next = self%cuda_wide_ptr(:level_count)
        do s = 1, self%node_count
            if (self%panels(s)%width <= 8) cycle
            l = levels(s)+1
            self%cuda_wide_members(next(l)) = s
            next(l) = next(l)+1
        end do
        base = 0
        do s = 1, self%node_count
            m2 = self%panels(s)%m-self%panels(s)%width
            tail_rows(base+1:base+m2) = self%panels(s)%rows(self%panels(s)%width+1:)
            base = base+m2
        end do
        ! Check the Fortran signed-byte range and size_t conversion before C.
        call memory_checked_mul(8_i8, int(self%n,i8), rhs_count, status)
        if (status%code /= FORTSPARSE_OK) return
        call memory_checked_add(int(self%numeric_count,i8), &
            max(int(maxwork,i8),rhs_count), arena_count, status)
        if (status%code /= FORTSPARSE_OK) return
        if (arena_count > int(huge(0_c_size_t),i8) .or. &
            self%device_limit > int(huge(0_c_size_t),i8)) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, 'CUDA device byte range')
            return
        end if
        self%device_payload_tracked = .true.
        code = spd_cuda_create_deferred_limited(self%cuda, int(self%device,c_int), &
            int(arena_count,c_size_t), int(self%device_limit,c_size_t), report)
        if (code == 0 .and. needs_blas) code = spd_cuda_prepare_blas(self%cuda)
        if (code == 0) then
            if (count_maps > 0) &
                code = spd_cuda_set_maps(self%cuda, c_loc(maps), int(count_maps,c_size_t))
        end if
        if (code == 0) then
            code = spd_cuda_set_levels(self%cuda, c_loc(panels), int(self%node_count,c_size_t), &
                c_loc(tasks), int(v,c_size_t), int(level_count,c_size_t))
        end if
        if (code == 0) then
            if (count_tail > 0) then
                code = spd_cuda_set_solve_plan(self%cuda, int(self%n,c_int), c_loc(panels), &
                    int(self%node_count,c_size_t), c_loc(tail_rows), int(count_tail,c_size_t), &
                    int(level_count,c_size_t))
            else
                code = spd_cuda_set_solve_plan(self%cuda, int(self%n,c_int), c_loc(panels), &
                    int(self%node_count,c_size_t), c_null_ptr, 0_c_size_t, int(level_count,c_size_t))
            end if
        end if
        ! A create refusal has no context; its output report survives cleanup.
        ! Later setter failures keep the old/partial context until build cleanup.
        if (c_associated(self%cuda)) then
            queried = spd_cuda_memory_stats(self%cuda, report)
            if (code == 0) code = queried
            self%cuda_host_bytes = int(spd_cuda_host_plan_bytes(self%cuda),i8)
        end if
        self%device_bytes = report%current_bytes
        self%device_peak_bytes = report%peak_bytes
        self%device_required_bytes = report%required_bytes
        self%device_available_bytes = report%available_bytes
        call cuda_status(code, status)
    end subroutine setup_cuda

    subroutine factor_cuda(self, status)
        class(supernodes_t), intent(inout) :: self
        type(fortsparse_status_t), intent(out) :: status
        integer :: s, u, m, w, m2, dest, l, member
        integer(c_int) :: code, info, ignored
        integer(c_size_t) :: offset, source_offset
        real(dp), pointer, contiguous :: host(:)

        host => self%host_values
        code = spd_cuda_begin_factor(self%cuda)
        if (code == 0 .and. self%numeric_count > 0) &
            code = spd_cuda_upload(self%cuda, 0_c_size_t, c_loc(host), &
                self%numeric_count)
        ! All sources in one level are independent. Every update writes to a
        ! strictly later level. Atomic thin contributions may accumulate in a
        ! hardware-dependent FP64 order; independent residual checking applies.
        do l = 1, size(self%cuda_wide_ptr)-1
            if (code /= 0) exit
            code = spd_cuda_factor_level(self%cuda, int(l-1,c_size_t))
            do member = self%cuda_wide_ptr(l), self%cuda_wide_ptr(l+1)-1
                if (code /= 0) exit
                s = self%cuda_wide_members(member)
                m = self%panels(s)%m
                w = self%panels(s)%width
                offset = self%panels(s)%offset
                code = spd_cuda_panel_factor(self%cuda, offset, w, m)
                if (code == 0 .and. m > w) &
                    code = spd_cuda_panel_solve(self%cuda, offset, w, m, offset+w, m-w, m)
            end do
            if (code == 0) code = spd_cuda_update_level(self%cuda, int(l-1,c_size_t))
            do member = self%cuda_wide_ptr(l), self%cuda_wide_ptr(l+1)-1
                if (code /= 0) exit
                s = self%cuda_wide_members(member)
                m = self%panels(s)%m
                w = self%panels(s)%width
                offset = self%panels(s)%offset
                do u = self%panels(s)%update_first, self%panels(s)%update_last
                    if (code /= 0) exit
                    dest = self%updates(u)%dest
                    m2 = m-self%updates(u)%first+1
                    source_offset = offset+self%updates(u)%first-1
                    code = spd_cuda_schur_product(self%cuda, self%numeric_count, m2, &
                        source_offset, m, source_offset, m, m2, self%updates(u)%q, w)
                    if (code /= 0) exit
                    code = spd_cuda_scatter_update(self%cuda, self%numeric_count, &
                        m2, self%updates(u)%q, self%panels(dest)%offset, &
                        self%panels(dest)%m, self%panels(dest)%width, &
                        self%updates(u)%row_map, self%updates(u)%column_map)
                end do
            end do
        end do
        if (code == 0) code = spd_cuda_factor_info(self%cuda, info)
        if (code /= 0) then
            ignored = spd_cuda_synchronize(self%cuda)
            call cuda_status(code, status)
            return
        end if
        if (info /= 0) then
            call status_set(status, NOT_SPD, 'CUDA supernode: nonpositive or nonfinite pivot')
            return
        end if
        ! factor_info has synchronized and checked every pivot. The numeric
        ! factors stay on device; host_values remains the assembly input.
        call cuda_status(code, status)
        self%ready = code == 0
    end subroutine factor_cuda

    subroutine solve_cuda_tile(self, host, nrhs, status)
        class(supernodes_t), intent(inout) :: self
        type(c_ptr), intent(in) :: host
        integer, intent(in) :: nrhs
        type(fortsparse_status_t), intent(out) :: status
        integer(c_int) :: code, waited
        integer(c_size_t) :: count
        count = int(self%n,c_size_t)*int(nrhs,c_size_t)
        code = spd_cuda_upload(self%cuda, self%numeric_count, host, count)
        if (code == 0) code = spd_cuda_solve(self%cuda, self%numeric_count, &
            int(self%n,c_int), int(nrhs,c_int))
        if (code == 0) code = spd_cuda_download(self%cuda, self%numeric_count, host, count)
        ! Always finish pending copies before a private TARGET dummy expires.
        waited = spd_cuda_synchronize(self%cuda)
        if (code == 0) code = waited
        if (code /= 0) self%ready = .false.
        call cuda_status(code, status)
    end subroutine solve_cuda_tile

    subroutine cuda_status(code, status)
        integer(c_int), intent(in) :: code
        type(fortsparse_status_t), intent(out) :: status
        integer :: status_code
        character(120) :: message
        status_code = FORTSPARSE_INTERNAL_ERROR
        if (code == 0) then
            call status_set(status, FORTSPARSE_OK, '')
            return
        end if
        if (code == -2) status_code = FORTSPARSE_RESOURCE_LIMIT
        if (code == -3) status_code = FORTSPARSE_INDEX_OVERFLOW
        if (code == 2 .or. code == 1003) status_code = FORTSPARSE_OUT_OF_MEMORY
        if (code == 10) status_code = FORTSPARSE_BACKEND_UNAVAILABLE
        write(message, '(a,i0)') 'native CUDA operation failed: ', code
        call status_set(status, status_code, message)
    end subroutine cuda_status
#endif

    subroutine allocation_error(self, status)
        class(supernodes_t), intent(inout) :: self
        type(fortsparse_status_t), intent(out) :: status
        self%ready = .false.
        call status_set(status, FORTSPARSE_OUT_OF_MEMORY, 'supernode allocation failed')
    end subroutine allocation_error

    real(dp) function seconds()
        integer(i8) :: clock, rate
        call system_clock(clock, rate)
        seconds = real(clock,dp)/real(rate,dp)
    end function seconds
end module fortsparse_native_supernodes
