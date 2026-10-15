module fortsparse_native_supernodes
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csc, only: csc_t, csc_is_valid
    use fortsparse_status, only: fortsparse_status_t, status_set, FORTSPARSE_OK, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OUT_OF_MEMORY, &
        FORTSPARSE_NOT_FACTORED, FORTSPARSE_INDEX_OVERFLOW
    use fortsparse_spd_leaves, only: spd_panel_factor, spd_panel_solve, &
        spd_scatter_mapped_update, spd_scatter_update
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    !$ use omp_lib, only: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
    implicit none
    private
    integer, parameter :: WIDTH = 16, NOT_SPD = 8

    type :: panel_t
        integer :: first = 0, width = 0, m = 0
        integer, allocatable :: rows(:)
        real(dp), allocatable :: values(:,:)
    end type panel_t

    type, public :: supernodes_t
        integer :: n = 0, node_count = 0, active_threads = 1
        integer(i8) :: factor_nnz = 0_i8
        real(dp) :: flops = 0.0_dp, numeric_s = 0.0_dp
        logical, private :: ready = .false.
        type(panel_t), allocatable, private :: panels(:)
        integer, allocatable, private :: colnode(:), tokens(:), failures(:), rowmap(:,:)
        real(dp), allocatable, private :: update_work(:,:), solve_work(:)
    contains
        procedure :: build
        procedure :: factor
        procedure :: solve
        procedure :: free => supernodes_free
        procedure :: memory_bytes
    end type supernodes_t

    interface
        subroutine dgemm(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import dp
            character, intent(in) :: ta, tb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            real(dp), intent(in) :: alpha, beta, a(lda,*), b(ldb,*)
            real(dp), intent(inout) :: c(ldc,*)
        end subroutine dgemm
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
    subroutine build(self, n, lp, li, nthreads, status)
        class(supernodes_t), intent(inout) :: self
        integer, intent(in) :: n, lp(:), li(:), nthreads
        type(fortsparse_status_t), intent(out) :: status
        integer, allocatable :: mark(:)
        integer :: s, j, p, row, m, first, last, group, q, dest, ierr, maxwork, maxrows
        integer(i8) :: update_size

        call self%free()
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'invalid symbolic Cholesky pattern')
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
        self%node_count = n/WIDTH
        if (mod(n,WIDTH) /= 0) self%node_count = self%node_count+1
        allocate(self%panels(self%node_count), self%colnode(n), &
            self%tokens(self%node_count), self%failures(self%node_count), mark(n), stat=ierr)
        if (ierr /= 0) then
            call allocation_error(self, status)
            return
        end if
        mark = 0
        self%tokens = 0
        self%failures = 0
        do j = 1, n
            self%colnode(j) = 1+(j-1)/WIDTH
        end do
        maxrows = 0
        do s = 1, self%node_count
            first = 1+(s-1)*WIDTH
            last = min(n, first+WIDTH-1)
            mark(first:last) = s
            do j = first, last
                mark(li(lp(j):lp(j+1)-1)) = s
            end do
            m = count(mark(first:n) == s)
            self%panels(s)%first = first
            self%panels(s)%width = last-first+1
            self%panels(s)%m = m
            allocate(self%panels(s)%rows(m), &
                self%panels(s)%values(m,last-first+1), stat=ierr)
            if (ierr /= 0) then
                call allocation_error(self, status)
                return
            end if
            p = 0
            do row = first, n
                if (mark(row) /= s) cycle
                p = p+1
                self%panels(s)%rows(p) = row
            end do
            self%factor_nnz = self%factor_nnz + int(m,i8)*(last-first+1) &
                -int(last-first+1,i8)*(last-first)/2
            self%flops = self%flops+real(m,dp)*(last-first+1)**2
            maxrows = max(maxrows, m)
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
        allocate(self%rowmap(n,self%active_threads), &
            self%update_work(maxwork,self%active_threads), self%solve_work(maxrows), stat=ierr)
        if (ierr /= 0) then
            call allocation_error(self, status)
            return
        end if
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine build

    subroutine factor(self, a, status)
        class(supernodes_t), target, intent(inout) :: self
        type(csc_t), intent(in) :: a
        type(fortsparse_status_t), intent(out) :: status
        integer, pointer :: token(:)
        integer :: s, j, p, row, position, group, first, q, dest
        real(dp) :: started

        self%ready = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'incompatible supernode matrix')
        if (.not. allocated(self%panels)) return
        if (.not. csc_is_valid(a)) return
        if (a%nrow /= self%n .or. a%ncol /= self%n) return
        if (.not. all(ieee_is_finite(a%val))) return
        started = seconds()
        self%failures = 0
        do s = 1, self%node_count
            self%panels(s)%values = 0.0_dp
            self%rowmap(:,1) = 0
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
        end do
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
                dest = self%colnode(self%panels(s)%rows(first))
                group = first+1
                do while (group <= self%panels(s)%m)
                    if (self%colnode(self%panels(s)%rows(group)) /= dest) exit
                    group = group+1
                end do
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
        self%numeric_s = seconds()-started
        if (any(self%failures /= 0)) then
            call status_set(status, NOT_SPD, 'supernode: nonpositive or nonfinite pivot')
            return
        end if
        self%ready = .true.
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine factor

    subroutine factor_panel(self, s)
        class(supernodes_t), intent(inout) :: self
        integer, intent(in) :: s
        integer :: info, w, m
        if (self%failures(s) /= 0) return
        w = self%panels(s)%width
        m = self%panels(s)%m
        call spd_panel_factor(self%panels(s)%values, w, m, info)
        if (info /= 0) then
            self%failures(s) = NOT_SPD
            return
        end if
        if (m > w) call spd_panel_solve(self%panels(s)%values, w, m, &
            self%panels(s)%values(w+1,1), m-w, m)
        if (.not. all(ieee_is_finite(self%panels(s)%values))) self%failures(s) = NOT_SPD
    end subroutine factor_panel

    subroutine update_panel(self, source, dest, first, q)
        class(supernodes_t), intent(inout) :: self
        integer, intent(in) :: source, dest, first, q
        integer :: thread, m, w, m2, j, p, column, offset
        if (self%failures(source) /= 0) self%failures(dest) = self%failures(source)
        if (self%failures(dest) /= 0) return
        thread = 1
        !$ thread = 1+omp_get_thread_num()
        m = self%panels(source)%m
        w = self%panels(source)%width
        m2 = m-first+1
        call dgemm('N', 'T', m2, q, w, -1.0_dp, &
            self%panels(source)%values(first,1), m, &
            self%panels(source)%values(first,1), m, 0.0_dp, self%update_work(1,thread), m2)
        self%rowmap(:,thread) = 0
        do p = 1, self%panels(dest)%m
            self%rowmap(self%panels(dest)%rows(p),thread) = p
        end do
        do j = 1, q
            column = self%panels(source)%rows(first+j-1)-self%panels(dest)%first+1
            offset = (j-1)*m2+j
            call spd_scatter_mapped_update(self%panels(dest)%values(1,column), &
                self%panels(source)%rows(first+j-1:), self%rowmap(:,thread), &
                self%update_work(offset:,thread), m2-j+1, -1.0_dp)
        end do
    end subroutine update_panel

    subroutine solve(self, b, status)
        class(supernodes_t), intent(inout) :: self
        real(dp), contiguous, intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: s, m, w, first, p
        call status_set(status, FORTSPARSE_NOT_FACTORED, 'supernode factor unavailable')
        if (.not. self%ready) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, 'incompatible supernode RHS')
        if (size(b) /= self%n) return
        do s = 1, self%node_count
            m = self%panels(s)%m
            w = self%panels(s)%width
            first = self%panels(s)%first
            call dtrsv('L', 'N', 'N', w, self%panels(s)%values, m, b(first:), 1)
            if (m == w) cycle
            call dgemv('N', m-w, w, 1.0_dp, self%panels(s)%values(w+1,1), m, &
                b(first:), 1, 0.0_dp, self%solve_work, 1)
            call spd_scatter_update(b, self%panels(s)%rows(w+1:), self%solve_work, m-w, 1.0_dp)
        end do
        do s = self%node_count, 1, -1
            m = self%panels(s)%m
            w = self%panels(s)%width
            first = self%panels(s)%first
            if (m > w) then
                do p = w+1, m
                    self%solve_work(p-w) = b(self%panels(s)%rows(p))
                end do
                call dgemv('T', m-w, w, -1.0_dp, self%panels(s)%values(w+1,1), m, &
                    self%solve_work, 1, 1.0_dp, b(first:), 1)
            end if
            call dtrsv('L', 'T', 'N', w, self%panels(s)%values, m, b(first:), 1)
        end do
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine solve

    subroutine supernodes_free(self)
        class(supernodes_t), intent(inout) :: self
        if (allocated(self%panels)) deallocate(self%panels)
        if (allocated(self%colnode)) deallocate(self%colnode)
        if (allocated(self%tokens)) deallocate(self%tokens)
        if (allocated(self%failures)) deallocate(self%failures)
        if (allocated(self%rowmap)) deallocate(self%rowmap)
        if (allocated(self%update_work)) deallocate(self%update_work)
        if (allocated(self%solve_work)) deallocate(self%solve_work)
        self%n = 0
        self%node_count = 0
        self%active_threads = 1
        self%factor_nnz = 0_i8
        self%flops = 0.0_dp
        self%numeric_s = 0.0_dp
        self%ready = .false.
    end subroutine supernodes_free

    integer(i8) function memory_bytes(self) result(bytes)
        class(supernodes_t), intent(in) :: self
        integer :: s
        bytes = 0_i8
        if (.not. allocated(self%panels)) return
        do s = 1, self%node_count
            if (allocated(self%panels(s)%rows)) &
                bytes = bytes+4_i8*size(self%panels(s)%rows,kind=i8)
            if (allocated(self%panels(s)%values)) &
                bytes = bytes+8_i8*size(self%panels(s)%values,kind=i8)
        end do
        if (allocated(self%colnode)) bytes = bytes+4_i8*size(self%colnode,kind=i8)
        if (allocated(self%tokens)) bytes = bytes+4_i8*size(self%tokens,kind=i8)
        if (allocated(self%failures)) bytes = bytes+4_i8*size(self%failures,kind=i8)
        if (allocated(self%rowmap)) bytes = bytes+4_i8*size(self%rowmap,kind=i8)
        if (allocated(self%update_work)) bytes = bytes+8_i8*size(self%update_work,kind=i8)
        if (allocated(self%solve_work)) bytes = bytes+8_i8*size(self%solve_work,kind=i8)
    end function memory_bytes

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
