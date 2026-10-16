module fortsparse_native_supernodes
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csc, only: csc_t, csc_is_valid
    use fortsparse_status, only: fortsparse_status_t, status_set, FORTSPARSE_OK, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OUT_OF_MEMORY, &
        FORTSPARSE_NOT_FACTORED, FORTSPARSE_INDEX_OVERFLOW, &
        FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_INTERNAL_ERROR
    use fortsparse_spd_leaves, only: spd_panel_factor, spd_panel_solve, &
        spd_scatter_mapped_update, spd_scatter_update
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, &
        c_int, c_size_t, c_loc
#ifdef FORTSPARSE_HAVE_NATIVE_CUDA
    use fortsparse_native_cuda, only: spd_cuda_create, spd_cuda_destroy, &
        spd_cuda_upload, spd_cuda_download, spd_cuda_synchronize, &
        spd_cuda_panel_factor, spd_cuda_factor_info, spd_cuda_panel_solve, &
        spd_cuda_schur_product, spd_cuda_set_maps, spd_cuda_begin_factor, &
        spd_cuda_scatter_update, spd_cuda_allocated_bytes
#endif
    !$ use omp_lib, only: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
    implicit none
    private
    integer, parameter :: EXACT_WIDTH = 64, RELAXED_WIDTH = 32, NOT_SPD = 8

    type :: panel_t
        integer :: first = 0, width = 0, m = 0
        integer :: update_first = 0, update_last = -1
        integer(c_size_t) :: offset = 0
        integer, allocatable :: rows(:)
        real(dp), pointer, contiguous :: values(:,:) => null()
    end type panel_t

    type :: update_t
        integer :: dest, first, q
        integer(c_size_t) :: row_map, column_map
    end type update_t

    type, public :: supernodes_t
        integer :: n = 0, node_count = 0, active_threads = 1, task_groups = 0
        integer(i8) :: factor_nnz = 0_i8
        integer(i8) :: device_bytes = 0_i8
        integer(i8) :: symbolic_nnz = 0_i8
        integer(i8) :: factor_storage_bytes = 0_i8, workspace_bytes = 0_i8
        real(dp) :: flops = 0.0_dp, numeric_s = 0.0_dp
        real(dp) :: group_target_flops = 0.0_dp
        logical, private :: ready = .false.
        type(panel_t), allocatable, private :: panels(:)
        integer, allocatable, private :: colnode(:), tokens(:), failures(:), rowmap(:,:)
        integer, allocatable, private :: incoming_ptr(:), incoming_source(:)
        integer, allocatable, private :: incoming_first(:), incoming_q(:)
        integer, allocatable, private :: group_first(:), group_last(:), group_tokens(:)
        integer, allocatable, private :: group_ptr(:), group_sources(:)
        real(dp), allocatable, private :: update_work(:,:), solve_work(:)
        ! One owning allocation permits one factor upload/download. Panel
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
        procedure :: free => supernodes_free
        procedure :: memory_bytes
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
    subroutine build(self, n, lp, li, nthreads, status, device)
        class(supernodes_t), intent(inout) :: self
        integer, intent(in) :: n, lp(:), li(:), nthreads
        type(fortsparse_status_t), intent(out) :: status
        integer, optional, intent(in) :: device
        integer, allocatable :: mark(:), starts(:), ends(:), counts(:), rowlist(:)
        integer :: s, j, p, row, m, first, last, group, q, dest, ierr, maxwork, maxrows
        integer :: candidate_m, candidate_w
        integer(i8) :: update_size, original_count, candidate_count, padded

        call self%free()
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'invalid symbolic Cholesky pattern')
        if (present(device)) self%device = device
        if (self%device < -1) return
#ifndef FORTSPARSE_HAVE_NATIVE_CUDA
        if (self%device >= 0) then
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, 'native CUDA is not built')
            return
        end if
#endif
        if (n < 0 .or. size(lp) /= n+1) return
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
        allocate(mark(n), starts(n), ends(n), counts(n), rowlist(n), stat=ierr)
        if (ierr /= 0) then
            call allocation_error(self, status)
            return
        end if
        mark = 0
        first = 1
        do while (first <= n)
            self%node_count = self%node_count+1
            s = self%node_count
            last = first
            m = lp(first+1)-lp(first)
            original_count = m
            mark(li(lp(first):lp(first+1)-1)) = s
            do while (last < n)
                ! A merge follows the elimination-tree parent chain.
                if (lp(last+1)-lp(last) < 2) exit
                if (li(lp(last)+1) /= last+1) exit
                candidate_w = last-first+2
                if (candidate_w > EXACT_WIDTH) exit
                candidate_m = m
                do p = lp(last+1), lp(last+2)-1
                    row = li(p)
                    if (mark(row) == s) cycle
                    mark(row) = s
                    candidate_m = candidate_m+1
                end do
                candidate_count = original_count+lp(last+2)-lp(last+1)
                padded = int(candidate_m,i8)*candidate_w &
                    -int(candidate_w,i8)*(candidate_w-1)/2
                if (padded /= candidate_count) then
                    if (candidate_w > RELAXED_WIDTH) exit
                    if (4_i8*padded > 5_i8*candidate_count) exit
                end if
                last = last+1
                m = candidate_m
                original_count = candidate_count
            end do
            starts(s) = first
            ends(s) = last
            counts(s) = m
            first = last+1
        end do
        allocate(self%panels(self%node_count), self%colnode(n), &
            self%tokens(self%node_count), self%failures(self%node_count), stat=ierr)
        if (ierr /= 0) then
            call allocation_error(self, status)
            return
        end if
        mark = 0
        self%tokens = 0
        self%failures = 0
        maxrows = 0
        do s = 1, self%node_count
            first = starts(s)
            last = ends(s)
            m = counts(s)
            q = 0
            do j = first, last
                self%colnode(j) = s
                do p = lp(j), lp(j+1)-1
                    row = li(p)
                    if (mark(row) == s) cycle
                    mark(row) = s
                    q = q+1
                    rowlist(q) = row
                end do
            end do
            call sort_rows(rowlist, 1, m)
            self%panels(s)%first = first
            self%panels(s)%width = last-first+1
            self%panels(s)%m = m
            allocate(self%panels(s)%rows(m), stat=ierr)
            if (ierr /= 0) then
                call allocation_error(self, status)
                return
            end if
            self%panels(s)%rows = rowlist(1:m)
            self%factor_nnz = self%factor_nnz + int(m,i8)*(last-first+1) &
                -int(last-first+1,i8)*(last-first)/2
            self%panels(s)%offset = self%numeric_count
            self%numeric_count = self%numeric_count+int(m,c_size_t)*(last-first+1)
            self%flops = self%flops+real(m,dp)*(last-first+1)**2
            maxrows = max(maxrows, m)
        end do
        self%factor_storage_bytes = 8_i8*int(self%numeric_count,i8)
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
                dest = self%colnode(self%panels(s)%rows(first))
                group = first+1
                do while (group <= m)
                    if (self%colnode(self%panels(s)%rows(group)) /= dest) exit
                    group = group+1
                end do
                q = group-first
                update_size = int(m-first+1,i8)*q
                if (update_size > huge(maxwork)) then
                    call self%free()
                    call status_set(status, FORTSPARSE_INDEX_OVERFLOW, 'supernode update too large')
                    return
                end if
                maxwork = max(maxwork, int(update_size))
                self%flops = self%flops+2.0_dp*real(update_size,dp)*self%panels(s)%width
                first = group
            end do
        end do
        self%active_threads = 1
        if (nthreads == 0) then
            !$ self%active_threads = omp_get_max_threads()
        end if
        if (nthreads > 0) then
            !$ self%active_threads = nthreads
        end if
        if (self%flops < 1.0e6_dp) self%active_threads = 1
        if (self%device >= 0) self%active_threads = 1
        if (self%device < 0) then
            call setup_incoming(self, status)
            if (status%code == FORTSPARSE_OK) call setup_groups(self, status)
            if (status%code /= FORTSPARSE_OK) then
                call self%free()
                return
            end if
            self%active_threads = min(self%active_threads, max(1,self%task_groups/4))
        end if
        allocate(self%rowmap(n,self%active_threads), &
            self%solve_work(maxrows), stat=ierr)
        if (ierr /= 0) then
            call allocation_error(self, status)
            return
        end if
        self%rowmap = 0
#ifdef FORTSPARSE_HAVE_NATIVE_CUDA
        if (self%device >= 0) then
            call setup_cuda(self, maxwork, status)
            if (status%code /= FORTSPARSE_OK) call self%free()
            self%workspace_bytes = self%memory_bytes()-self%factor_storage_bytes
            return
        end if
#endif
        allocate(self%update_work(maxwork,self%active_threads), stat=ierr)
        if (ierr /= 0) then
            call allocation_error(self, status)
            return
        end if
        self%workspace_bytes = self%memory_bytes()-self%factor_storage_bytes
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine build

    subroutine factor(self, a, status)
        class(supernodes_t), target, intent(inout) :: self
        type(csc_t), intent(in) :: a
        type(fortsparse_status_t), intent(out) :: status
        integer, pointer :: token(:)
#ifdef FORTSPARSE_HAVE_OMP_DEPEND_ITERATOR
        integer, pointer :: sources(:)
        integer :: incoming_start, incoming_end
#endif
        integer :: s, j, p, row, position
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
        started = seconds()
        self%failures = 0
        self%rowmap = 0
        do s = 1, self%node_count
            self%panels(s)%values = 0.0_dp
            do p = 1, self%panels(s)%m
                self%rowmap(self%panels(s)%rows(p),1) = p
            end do
            do j = 1, self%panels(s)%width
                row = self%panels(s)%first+j-1
                do p = a%col_ptr(row), a%col_ptr(row+1)-1
                    if (a%row_idx(p) < row) cycle
                    position = self%rowmap(a%row_idx(p),1)
                    if (position == 0) return
                    self%panels(s)%values(position,j) = a%val(p)
                end do
            end do
            do p = 1, self%panels(s)%m
                self%rowmap(self%panels(s)%rows(p),1) = 0
            end do
        end do
#ifdef FORTSPARSE_HAVE_NATIVE_CUDA
        if (self%device >= 0) then
            call factor_cuda(self, status)
            self%numeric_s = seconds()-started
            return
        end if
#endif
        if (self%active_threads == 1) then
            call factor_serial(self)
        else
#ifdef FORTSPARSE_HAVE_OMP_DEPEND_ITERATOR
            token => self%group_tokens
            sources => self%group_sources
            !$omp parallel num_threads(self%active_threads) default(shared)
            !$omp single
            !$ self%active_threads = omp_get_num_threads()
            do s = 1, self%task_groups
                incoming_start = self%group_ptr(s)
                ! OpenMP iterator end is exclusive, including in Fortran.
                incoming_end = self%group_ptr(s+1)
                !$omp task firstprivate(s) &
                !$omp depend(iterator(incoming=incoming_start:incoming_end), &
                !$omp in:token(sources(incoming))) depend(out:token(s))
                call factor_group(self, s)
                !$omp end task
            end do
            !$omp end single
            !$omp end parallel
#else
            token => self%tokens
            !$omp parallel num_threads(self%active_threads) default(shared)
            !$omp single
            !$ self%active_threads = omp_get_num_threads()
            do s = 1, self%node_count
                !$omp task firstprivate(s) depend(inout:token(s))
                call factor_panel(self, s)
                !$omp end task
                first = self%panels(s)%width+1
                do while (first <= self%panels(s)%m)
                    call update_group(self, s, first, group, dest)
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
        integer :: s
        do s = self%group_first(group), self%group_last(group)
            call factor_destination(self, s)
        end do
    end subroutine factor_group
#endif

    subroutine factor_destination(self, dest)
        class(supernodes_t), intent(inout) :: self
        integer, intent(in) :: dest
        integer :: p, thread
        thread = workspace_thread(self)
        do p = 1, self%panels(dest)%m
            self%rowmap(self%panels(dest)%rows(p),thread) = p
        end do
        do p = self%incoming_ptr(dest), self%incoming_ptr(dest+1)-1
            call accumulate_update(self, self%incoming_source(p), dest, &
                self%incoming_first(p), self%incoming_q(p), thread)
        end do
        call factor_panel(self, dest)
        do p = 1, self%panels(dest)%m
            self%rowmap(self%panels(dest)%rows(p),thread) = 0
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
        if (m > w) call spd_panel_solve(values, w, m, tail, m-w, m)
        if (.not. all(ieee_is_finite(self%panels(s)%values))) self%failures(s) = NOT_SPD
    end subroutine factor_panel

#ifndef FORTSPARSE_HAVE_OMP_DEPEND_ITERATOR
    subroutine update_panel(self, source, dest, first, q)
        class(supernodes_t), intent(inout) :: self
        integer, intent(in) :: source, dest, first, q
        integer :: thread, p
        thread = workspace_thread(self)
        do p = 1, self%panels(dest)%m
            self%rowmap(self%panels(dest)%rows(p),thread) = p
        end do
        call accumulate_update(self, source, dest, first, q, thread)
        do p = 1, self%panels(dest)%m
            self%rowmap(self%panels(dest)%rows(p),thread) = 0
        end do
    end subroutine update_panel
#endif

    integer function workspace_thread(self) result(thread)
        class(supernodes_t), intent(in) :: self
        thread = 1
        if (self%active_threads > 1) then
            !$ thread = 1+omp_get_thread_num()
        end if
    end function workspace_thread

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
        ! Exact dense intersections can update their owning panel directly.
        ! SYRK touches only the diagonal block's lower triangle; GEMM handles
        ! its rectangular tail. Irregular intersections retain mapped scatter.
        position = self%rowmap(self%panels(source)%rows(first),thread)
        contiguous_rows = position > 0
        if (contiguous_rows) then
            do p = first+1, m
                if (self%rowmap(self%panels(source)%rows(p),thread) &
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
            0.0_dp, self%update_work(1,thread), m2)
        do j = 1, q
            column = self%panels(source)%rows(first+j-1)-self%panels(dest)%first+1
            offset = (j-1)*m2+j
            dest_values => self%panels(dest)%values(:,column)
            call spd_scatter_mapped_update(dest_values, &
                self%panels(source)%rows(first+j-1:), self%rowmap(:,thread), &
                self%update_work(offset:,thread), m2-j+1, -1.0_dp)
        end do
    end subroutine accumulate_update

    subroutine solve(self, b, status)
        class(supernodes_t), intent(inout) :: self
        real(dp), contiguous, intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: s, m, w, first, p
        real(dp), pointer, contiguous :: values(:,:), tail(:)
        call status_set(status, FORTSPARSE_NOT_FACTORED, 'supernode factor unavailable')
        if (.not. self%ready) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'incompatible supernode RHS')
        if (size(b) /= self%n) return
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

    subroutine supernodes_free(self)
        class(supernodes_t), intent(inout) :: self
#ifdef FORTSPARSE_HAVE_NATIVE_CUDA
        integer(c_int) :: code
        if (c_associated(self%cuda)) code = spd_cuda_destroy(self%cuda)
#endif
        self%cuda = c_null_ptr
        if (allocated(self%panels)) deallocate(self%panels)
        if (associated(self%host_values)) deallocate(self%host_values)
        if (allocated(self%updates)) deallocate(self%updates)
        if (allocated(self%colnode)) deallocate(self%colnode)
        if (allocated(self%tokens)) deallocate(self%tokens)
        if (allocated(self%failures)) deallocate(self%failures)
        if (allocated(self%rowmap)) deallocate(self%rowmap)
        if (allocated(self%incoming_ptr)) deallocate(self%incoming_ptr)
        if (allocated(self%incoming_source)) deallocate(self%incoming_source)
        if (allocated(self%incoming_first)) deallocate(self%incoming_first)
        if (allocated(self%incoming_q)) deallocate(self%incoming_q)
        if (allocated(self%group_first)) deallocate(self%group_first)
        if (allocated(self%group_last)) deallocate(self%group_last)
        if (allocated(self%group_tokens)) deallocate(self%group_tokens)
        if (allocated(self%group_ptr)) deallocate(self%group_ptr)
        if (allocated(self%group_sources)) deallocate(self%group_sources)
        if (allocated(self%update_work)) deallocate(self%update_work)
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
        self%device = -1
        self%flops = 0.0_dp
        self%numeric_s = 0.0_dp
        self%ready = .false.
    end subroutine supernodes_free

    subroutine finalize_supernodes(self)
        type(supernodes_t), intent(inout) :: self
        call self%free()
    end subroutine finalize_supernodes

    integer(i8) function memory_bytes(self) result(bytes)
        class(supernodes_t), intent(in) :: self
        integer :: s
        bytes = 0_i8
        if (.not. allocated(self%panels)) return
        do s = 1, self%node_count
            if (allocated(self%panels(s)%rows)) &
                bytes = bytes+4_i8*size(self%panels(s)%rows,kind=i8)
            if (associated(self%panels(s)%values)) &
                bytes = bytes+8_i8*size(self%panels(s)%values,kind=i8)
        end do
        if (allocated(self%colnode)) bytes = bytes+4_i8*size(self%colnode,kind=i8)
        if (allocated(self%tokens)) bytes = bytes+4_i8*size(self%tokens,kind=i8)
        if (allocated(self%failures)) bytes = bytes+4_i8*size(self%failures,kind=i8)
        if (allocated(self%rowmap)) bytes = bytes+4_i8*size(self%rowmap,kind=i8)
        if (allocated(self%incoming_ptr)) bytes = bytes+4_i8*size(self%incoming_ptr,kind=i8)
        if (allocated(self%incoming_source)) bytes = bytes+4_i8*size(self%incoming_source,kind=i8)
        if (allocated(self%incoming_first)) bytes = bytes+4_i8*size(self%incoming_first,kind=i8)
        if (allocated(self%incoming_q)) bytes = bytes+4_i8*size(self%incoming_q,kind=i8)
        if (allocated(self%group_first)) bytes = bytes+4_i8*size(self%group_first,kind=i8)
        if (allocated(self%group_last)) bytes = bytes+4_i8*size(self%group_last,kind=i8)
        if (allocated(self%group_tokens)) bytes = bytes+4_i8*size(self%group_tokens,kind=i8)
        if (allocated(self%group_ptr)) bytes = bytes+4_i8*size(self%group_ptr,kind=i8)
        if (allocated(self%group_sources)) bytes = bytes+4_i8*size(self%group_sources,kind=i8)
        if (allocated(self%update_work)) bytes = bytes+8_i8*size(self%update_work,kind=i8)
        if (allocated(self%solve_work)) bytes = bytes+8_i8*size(self%solve_work,kind=i8)
        if (allocated(self%updates)) &
            bytes = bytes+int(storage_size(self%updates)/8,i8)*size(self%updates,kind=i8)
    end function memory_bytes

    subroutine setup_incoming(self, status)
        class(supernodes_t), intent(inout) :: self
        type(fortsparse_status_t), intent(out) :: status
        integer, allocatable :: next(:)
        integer :: s, first, group, dest, p, ierr
        integer(i8) :: entries
        allocate(self%incoming_ptr(self%node_count+1), next(self%node_count), stat=ierr)
        if (ierr /= 0) then
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, 'incoming task map allocation failed')
            return
        end if
        self%incoming_ptr = 0
        entries = 0_i8
        do s = 1, self%node_count
            first = self%panels(s)%width+1
            do while (first <= self%panels(s)%m)
                call update_group(self, s, first, group, dest)
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
        allocate(self%incoming_source(int(entries)), self%incoming_first(int(entries)), &
            self%incoming_q(int(entries)), stat=ierr)
        if (ierr /= 0) then
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, 'incoming task entries allocation failed')
            return
        end if
        next = self%incoming_ptr(:self%node_count)
        do s = 1, self%node_count
            first = self%panels(s)%width+1
            do while (first <= self%panels(s)%m)
                call update_group(self, s, first, group, dest)
                p = next(dest)
                self%incoming_source(p) = s
                self%incoming_first(p) = first
                self%incoming_q(p) = group-first
                next(dest) = p+1
                first = group
            end do
        end do
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine setup_incoming

    subroutine setup_groups(self, status)
        class(supernodes_t), intent(inout) :: self
        type(fortsparse_status_t), intent(out) :: status
        integer, allocatable :: node_group(:), first_node(:), last_node(:), mark(:)
        integer :: s, p, source, g, predecessor, ierr, position
        integer(i8) :: entries
        real(dp) :: work, accumulated

        ! About twelve groups per requested worker, with enough useful work
        ! to amortize task creation. A single expensive panel remains atomic.
        self%group_target_flops = max(1.0e6_dp, &
            self%flops/(12.0_dp*self%active_threads))
        allocate(node_group(self%node_count), first_node(self%node_count), &
            last_node(self%node_count), stat=ierr)
        if (ierr /= 0) then
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, 'task group allocation failed')
            return
        end if
        accumulated = self%group_target_flops
        do s = 1, self%node_count
            if (accumulated >= self%group_target_flops) then
                self%task_groups = self%task_groups+1
                first_node(self%task_groups) = s
                accumulated = 0.0_dp
            end if
            node_group(s) = self%task_groups
            last_node(self%task_groups) = s
            work = real(self%panels(s)%m,dp)*self%panels(s)%width**2
            do p = self%incoming_ptr(s), self%incoming_ptr(s+1)-1
                source = self%incoming_source(p)
                work = work+2.0_dp*(self%panels(source)%m-self%incoming_first(p)+1) &
                    *self%incoming_q(p)*self%panels(source)%width
            end do
            accumulated = accumulated+work
        end do
        allocate(self%group_first(self%task_groups), self%group_last(self%task_groups), &
            self%group_tokens(self%task_groups), self%group_ptr(self%task_groups+1), &
            mark(self%task_groups), stat=ierr)
        if (ierr /= 0) then
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, 'task group headers allocation failed')
            return
        end if
        self%group_first = first_node(:self%task_groups)
        self%group_last = last_node(:self%task_groups)
        self%group_tokens = 0
        self%group_ptr = 0
        mark = 0
        entries = 0_i8
        do g = 1, self%task_groups
            do s = self%group_first(g), self%group_last(g)
                do p = self%incoming_ptr(s), self%incoming_ptr(s+1)-1
                    predecessor = node_group(self%incoming_source(p))
                    if (predecessor == g) cycle
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
        allocate(self%group_sources(int(entries)), stat=ierr)
        if (ierr /= 0) then
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, 'task group edges allocation failed')
            return
        end if
        mark = 0
        do g = 1, self%task_groups
            position = self%group_ptr(g)
            do s = self%group_first(g), self%group_last(g)
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
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine setup_groups

    subroutine update_group(self, s, first, group, dest)
        class(supernodes_t), intent(in) :: self
        integer, intent(in) :: s, first
        integer, intent(out) :: group, dest
        dest = self%colnode(self%panels(s)%rows(first))
        group = first+1
        do while (group <= self%panels(s)%m)
            if (self%colnode(self%panels(s)%rows(group)) /= dest) exit
            group = group+1
        end do
    end subroutine update_group

    recursive subroutine sort_rows(rows, lo, hi)
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
        if (lo < j) call sort_rows(rows, lo, j)
        if (i < hi) call sort_rows(rows, i, hi)
    end subroutine sort_rows

#ifdef FORTSPARSE_HAVE_NATIVE_CUDA
    subroutine setup_cuda(self, maxwork, status)
        class(supernodes_t), intent(inout) :: self
        integer, intent(in) :: maxwork
        type(fortsparse_status_t), intent(out) :: status
        integer(c_int), allocatable, target :: maps(:)
        integer :: s, first, group, dest, q, m2, u, count_updates, base, p, ierr
        integer(i8) :: count_maps
        integer(c_int) :: code

        count_updates = 0
        count_maps = 0_i8
        do s = 1, self%node_count
            first = self%panels(s)%width+1
            do while (first <= self%panels(s)%m)
                call update_group(self, s, first, group, dest)
                count_updates = count_updates+1
                count_maps = count_maps+self%panels(s)%m-first+1+group-first
                first = group
            end do
        end do
        if (count_maps > huge(base)) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, 'CUDA scatter maps too large')
            return
        end if
        allocate(self%updates(count_updates), maps(int(count_maps)), stat=ierr)
        if (ierr /= 0) then
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, 'CUDA scatter map allocation failed')
            return
        end if
        base = 0
        u = 0
        do s = 1, self%node_count
            self%panels(s)%update_first = u+1
            first = self%panels(s)%width+1
            do while (first <= self%panels(s)%m)
                call update_group(self, s, first, group, dest)
                q = group-first
                m2 = self%panels(s)%m-first+1
                u = u+1
                self%updates(u)%dest = dest
                self%updates(u)%first = first
                self%updates(u)%q = q
                self%updates(u)%row_map = int(base,c_size_t)
                do p = 1, self%panels(dest)%m
                    self%rowmap(self%panels(dest)%rows(p),1) = p
                end do
                do p = first, self%panels(s)%m
                    maps(base+p-first+1) = self%rowmap(self%panels(s)%rows(p),1)
                end do
                base = base+m2
                self%updates(u)%column_map = int(base,c_size_t)
                maps(base+1:base+q) = self%panels(s)%rows(first:group-1) &
                    -self%panels(dest)%first+1
                base = base+q
                do p = 1, self%panels(dest)%m
                    self%rowmap(self%panels(dest)%rows(p),1) = 0
                end do
                first = group
            end do
            self%panels(s)%update_last = u
        end do
        code = spd_cuda_create(self%cuda, int(self%device,c_int), &
            self%numeric_count+int(maxwork,c_size_t))
        if (code == 0) then
            if (count_maps > 0) &
                code = spd_cuda_set_maps(self%cuda, c_loc(maps), int(count_maps,c_size_t))
        end if
        call cuda_status(code, status)
        if (code == 0) self%device_bytes = int(spd_cuda_allocated_bytes(self%cuda),i8)
    end subroutine setup_cuda

    subroutine factor_cuda(self, status)
        class(supernodes_t), intent(inout) :: self
        type(fortsparse_status_t), intent(out) :: status
        integer :: s, u, m, w, m2, dest
        integer(c_int) :: code, info, ignored
        integer(c_size_t) :: offset, source_offset
        real(dp), pointer, contiguous :: host(:)

        host => self%host_values
        code = spd_cuda_begin_factor(self%cuda)
        if (code == 0 .and. self%numeric_count > 0) &
            code = spd_cuda_upload(self%cuda, 0_c_size_t, c_loc(host), &
                self%numeric_count)
        do s = 1, self%node_count
            if (code /= 0) exit
            m = self%panels(s)%m
            w = self%panels(s)%width
            offset = self%panels(s)%offset
            code = spd_cuda_panel_factor(self%cuda, offset, w, m)
            if (code == 0 .and. m > w) &
                code = spd_cuda_panel_solve(self%cuda, offset, w, m, offset+w, m-w, m)
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
        if (self%numeric_count > 0) &
            code = spd_cuda_download(self%cuda, 0_c_size_t, c_loc(host), &
                self%numeric_count)
        if (code == 0) code = spd_cuda_synchronize(self%cuda)
        call cuda_status(code, status)
        self%ready = code == 0
    end subroutine factor_cuda

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
        if (code == 2 .or. code == 1003) status_code = FORTSPARSE_OUT_OF_MEMORY
        if (code == 10) status_code = FORTSPARSE_BACKEND_UNAVAILABLE
        write(message, '(a,i0)') 'native CUDA operation failed: ', code
        call status_set(status, status_code, message)
    end subroutine cuda_status
#endif

    subroutine allocation_error(self, status)
        class(supernodes_t), intent(inout) :: self
        type(fortsparse_status_t), intent(out) :: status
        call self%free()
        call status_set(status, FORTSPARSE_OUT_OF_MEMORY, 'supernode allocation failed')
    end subroutine allocation_error

    real(dp) function seconds()
        integer(i8) :: clock, rate
        call system_clock(clock, rate)
        seconds = real(clock,dp)/real(rate,dp)
    end function seconds
end module fortsparse_native_supernodes
