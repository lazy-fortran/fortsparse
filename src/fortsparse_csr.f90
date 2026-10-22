module fortsparse_csr
    !! Owned, canonical, one-based real CSR and setup-only sparse algebra.
    use fortsparse_kinds, only: dp, i8
    use fortsparse_status, only: fortsparse_status_t, status_set, &
                  FORTSPARSE_OK, FORTSPARSE_INVALID_MATRIX, FORTSPARSE_RESOURCE_LIMIT, &
                                 FORTSPARSE_OUT_OF_MEMORY, FORTSPARSE_INDEX_OVERFLOW
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
!$  use omp_lib, only: omp_get_max_threads, omp_get_thread_num
    implicit none
    private
    public :: csr_t, csr_valid, csr_matvec, csr_transpose, csr_multiply
    public :: csr_budget_t, csr_allocate, csr_release, csr_copy, csr_find

    type :: csr_t
        integer :: nrow = 0, ncol = 0
        integer, allocatable :: row_ptr(:), col(:)
        real(dp), allocatable :: val(:)
    end type

    type :: csr_budget_t
        integer(i8) :: limit = 0, live = 0, peak = 0
    contains
        procedure :: reserve
        procedure :: release
    end type
contains
    logical function reserve(self, bytes, info) result(ok)
        class(csr_budget_t), intent(inout) :: self
        integer(i8), intent(in) :: bytes
        type(fortsparse_status_t), intent(out) :: info
        ok = .false.
        call status_set(info, FORTSPARSE_RESOURCE_LIMIT, 'CSR/AMG payload limit')
        if (bytes < 0 .or. bytes > huge(bytes) - self%live) return
        if (self%limit > 0) then
            if (bytes > self%limit - self%live) return
        end if
        self%live = self%live + bytes
        self%peak = max(self%peak, self%live)
        call status_set(info, FORTSPARSE_OK, '')
        ok = .true.
    end function

    subroutine release(self, bytes)
        class(csr_budget_t), intent(inout) :: self
        integer(i8), intent(in) :: bytes
        self%live = self%live - bytes
    end subroutine

    subroutine csr_allocate(a, nr, nc, nz, budget, info)
        type(csr_t), intent(out) :: a
        integer, intent(in) :: nr, nc
        integer(i8), intent(in) :: nz
        type(csr_budget_t), intent(inout) :: budget
        type(fortsparse_status_t), intent(out) :: info
        integer :: io
        integer(i8) :: bytes
        call status_set(info, FORTSPARSE_INDEX_OVERFLOW, 'CSR index overflow')
        if (nr < 0 .or. nc < 0 .or. nr == huge(nr)) return
        if (nz < 0 .or. nz >= int(huge(nr), i8)) return
        bytes = 4_i8*(int(nr, i8) + 1) + 12_i8*nz
        if (.not. budget%reserve(bytes, info)) return
        allocate (a%row_ptr(nr + 1), a%col(nz), a%val(nz), stat=io)
        if (io /= 0) then
            call csr_release(a)
            call budget%release(bytes)
            call status_set(info, FORTSPARSE_OUT_OF_MEMORY, 'CSR allocation')
            return
        end if
        a%nrow = nr
        a%ncol = nc
    end subroutine

    subroutine csr_release(a, budget)
        type(csr_t), intent(inout) :: a
        type(csr_budget_t), intent(inout), optional :: budget
        if (allocated(a%row_ptr)) then
            if (present(budget)) call budget%release(4_i8*size(a%row_ptr, kind=i8))
            deallocate (a%row_ptr)
        end if
        if (allocated(a%col)) then
            if (present(budget)) call budget%release(4_i8*size(a%col, kind=i8))
            deallocate (a%col)
        end if
        if (allocated(a%val)) then
            if (present(budget)) call budget%release(8_i8*size(a%val, kind=i8))
            deallocate (a%val)
        end if
        a%nrow = 0
        a%ncol = 0
    end subroutine

    subroutine csr_copy(a, b, budget, info)
        type(csr_t), intent(in) :: a
        type(csr_t), intent(out) :: b
        type(csr_budget_t), intent(inout) :: budget
        type(fortsparse_status_t), intent(out) :: info
        call csr_allocate(b, a%nrow, a%ncol, size(a%val, kind=i8), budget, info)
        if (info%code /= 0) return
        b%row_ptr = a%row_ptr
        b%col = a%col
        b%val = a%val
    end subroutine

    logical function csr_valid(a) result(ok)
        type(csr_t), intent(in) :: a
        integer :: i, k, last
        ok = .false.
        if (a%nrow < 1 .or. a%ncol < 1) return
        if (.not. allocated(a%row_ptr)) return
        if (.not. allocated(a%col)) return
        if (.not. allocated(a%val)) return
        if (size(a%row_ptr, kind=i8) /= int(a%nrow, i8) + 1) return
        if (size(a%col) /= size(a%val)) return
        if (size(a%col) == huge(0)) return
        if (a%row_ptr(1) /= 1) return
        if (a%row_ptr(a%nrow + 1) /= size(a%col) + 1) return
        if (any(a%row_ptr < 1) .or. any(a%row_ptr > size(a%col) + 1)) return
        if (any(a%row_ptr(2:) < a%row_ptr(:a%nrow))) return
        if (.not. all(ieee_is_finite(a%val))) return
        do i = 1, a%nrow
            last = 0
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                if (a%col(k) <= last .or. a%col(k) > a%ncol) return
                last = a%col(k)
            end do
        end do
        ok = .true.
    end function

    integer function csr_find(a, row, col) result(k)
        type(csr_t), intent(in) :: a
        integer, intent(in) :: row, col
        integer :: lo, hi, mid
        lo = a%row_ptr(row)
        hi = a%row_ptr(row + 1) - 1
        k = 0
        do while (lo <= hi)
            mid = lo + (hi - lo)/2
            if (a%col(mid) == col) then
                k = mid
                return
            else if (a%col(mid) < col) then
                lo = mid + 1
            else
                hi = mid - 1
            end if
        end do
    end function

    subroutine csr_matvec(a, x, y)
        !! Shapes and canonical storage are checked by setup, not in hot loops.
        type(csr_t), intent(in) :: a
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: s
        integer :: i, k
        !$omp parallel do private(k,s) if(a%nrow >= 512)
        do i = 1, a%nrow
            s = 0
            !$omp simd reduction(+:s)
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                s = s + a%val(k)*x(a%col(k))
            end do
            y(i) = s
        end do
    end subroutine

    subroutine sort_row(cols, vals)
        integer, intent(inout) :: cols(:)
        real(dp), intent(inout) :: vals(:)
        integer :: j, k, c
        real(dp) :: v
        do j = 2, size(cols)
            c = cols(j)
            v = vals(j)
            k = j - 1
            do while (k >= 1)
                if (cols(k) <= c) exit
                cols(k + 1) = cols(k)
                vals(k + 1) = vals(k)
                k = k - 1
            end do
            cols(k + 1) = c
            vals(k + 1) = v
        end do
    end subroutine

    subroutine csr_transpose(a, b, budget, info)
        type(csr_t), intent(in) :: a
        type(csr_t), intent(out) :: b
        type(csr_budget_t), intent(inout) :: budget
        type(fortsparse_status_t), intent(out) :: info
        integer, allocatable :: next(:)
        integer :: i, k, j, p, io
        call status_set(info, FORTSPARSE_INVALID_MATRIX, 'CSR transpose storage')
        if (.not. csr_valid(a)) return
        call csr_allocate(b, a%ncol, a%nrow, size(a%val, kind=i8), budget, info)
        if (info%code /= 0) return
        if (.not. budget%reserve(4_i8*a%ncol, info)) return
        allocate (next(a%ncol), stat=io)
        if (io /= 0) then
            call budget%release(4_i8*a%ncol)
            call status_set(info, FORTSPARSE_OUT_OF_MEMORY, 'CSR transpose scratch')
            return
        end if
        next = 0
        !$omp parallel do private(j) if(a%nrow >= 512)
        do k = 1, size(a%col)
            j = a%col(k)
            !$omp atomic update
            next(j) = next(j) + 1
        end do
        b%row_ptr(1) = 1
        do i = 1, b%nrow
            b%row_ptr(i + 1) = b%row_ptr(i) + next(i)
        end do
        next = b%row_ptr(:b%nrow)
        !$omp parallel do private(k,j,p) if(a%nrow >= 512)
        do i = 1, a%nrow
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                j = a%col(k)
                !$omp atomic capture
                p = next(j)
                next(j) = next(j) + 1
                !$omp end atomic
                b%col(p) = i
                b%val(p) = a%val(k)
            end do
        end do
        !$omp parallel do if(b%nrow >= 512)
        do i = 1, b%nrow
            call sort_row(b%col(b%row_ptr(i):b%row_ptr(i + 1) - 1), &
                          b%val(b%row_ptr(i):b%row_ptr(i + 1) - 1))
        end do
        call budget%release(4_i8*a%ncol)
    end subroutine

    subroutine csr_multiply(a, b, c, budget, info)
        !! Two-pass row accumulation, O(ncol * OpenMP team) setup scratch.
        !! Keeps structural zeros so Galerkin products retain symmetric patterns.
        type(csr_t), intent(in) :: a, b
        type(csr_t), intent(out) :: c
        type(csr_budget_t), intent(inout) :: budget
        type(fortsparse_status_t), intent(out) :: info
        integer, allocatable :: mark(:, :), count(:)
        real(dp), allocatable :: work(:, :)
        integer :: nt, tid, i, j, k, q, col, nz, p, pass, io
        integer(i8) :: scratch, total
        call status_set(info, FORTSPARSE_INVALID_MATRIX, 'CSR product shape/storage')
        if (.not. csr_valid(a)) return
        if (.not. csr_valid(b)) return
        if (a%ncol /= b%nrow) return
        nt = 1
!$      nt = omp_get_max_threads()
        scratch = 12_i8*b%ncol*nt + 4_i8*a%nrow
        if (.not. budget%reserve(scratch, info)) return
        allocate (mark(b%ncol, nt), work(b%ncol, nt), count(a%nrow), stat=io)
        if (io /= 0) then
            call budget%release(scratch)
            call status_set(info, FORTSPARSE_OUT_OF_MEMORY, 'CSR product scratch')
            return
        end if
        do pass = 1, 2
            mark = 0
            !$omp parallel private(tid,i,j,k,q,col,nz,p)
            tid = 1
!$          tid = omp_get_thread_num() + 1
            !$omp do schedule(static)
            do i = 1, a%nrow
                nz = 0
                if (pass == 2) nz = c%row_ptr(i) - 1
                do j = a%row_ptr(i), a%row_ptr(i + 1) - 1
                    k = a%col(j)
                    do q = b%row_ptr(k), b%row_ptr(k + 1) - 1
                        col = b%col(q)
                        if (mark(col, tid) /= i) then
                            mark(col, tid) = i
                            nz = nz + 1
                            work(col, tid) = 0
                            if (pass == 2) c%col(nz) = col
                        end if
                        if (pass == 2) work(col, tid) = work(col, tid) &
                                                        + a%val(j)*b%val(q)
                    end do
                end do
                if (pass == 1) then
                    count(i) = nz
                else
                    do p = c%row_ptr(i), c%row_ptr(i + 1) - 1
                        c%val(p) = work(c%col(p), tid)
                    end do
                    call sort_row(c%col(c%row_ptr(i):c%row_ptr(i + 1) - 1), &
                                  c%val(c%row_ptr(i):c%row_ptr(i + 1) - 1))
                end if
            end do
            !$omp end do
            !$omp end parallel
            if (pass == 1) then
                total = sum(int(count, i8))
                call csr_allocate(c, a%nrow, b%ncol, total, budget, info)
                if (info%code /= 0) then
                    call budget%release(scratch)
                    return
                end if
                c%row_ptr(1) = 1
                do i = 1, c%nrow
                    c%row_ptr(i + 1) = c%row_ptr(i) + count(i)
                end do
            end if
        end do
        call budget%release(scratch)
    end subroutine
end module
