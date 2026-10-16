module fortsparse_native_hpd
    ! Natural-order simplicial complex Cholesky, A = L L^H.
    ! Sparse reach owns topology; retained FortSym leaves own arithmetic.
    ! Both policies use the same minimal retained workspace in this reference;
    ! LOW_MEMORY therefore changes no layout, precision or mathematical result.
    ! Payload counters exclude object descriptors, caller buffers and allocator/RSS.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csc, only: csc_z_t
    use fortsparse_memory, only: memory_budget_t
    use fortsparse_status, only: fortsparse_status_t, status_set, &
        FORTSPARSE_OK, FORTSPARSE_INVALID_MATRIX, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_OUT_OF_MEMORY, FORTSPARSE_INDEX_OVERFLOW, &
        FORTSPARSE_RESOURCE_LIMIT, FORTSPARSE_NOT_POSITIVE_DEFINITE, &
        FORTSPARSE_INTERNAL_ERROR
    use fortsparse_hpd_leaves, only: hpd_root, hpd_divide, hpd_update, &
        hpd_conj_update, hpd_norm_update
    implicit none
    private
    public :: native_hpd_t
    integer(i8), parameter :: ib = storage_size(0)/8
    integer(i8), parameter :: zb = storage_size(cmplx(0.0_dp, 0.0_dp, dp))/8

    type :: native_hpd_t
        private
        integer :: n = 0
        logical :: ready = .false., low_memory = .false.
        integer, allocatable :: ap(:), ai(:), lp(:), li(:)
        integer, allocatable :: parent(:), count(:), mark(:), reach(:), path(:)
        complex(dp), allocatable :: lv(:), work(:)
        type(memory_budget_t) :: budget
    contains
        procedure :: factor => hpd_factor
        procedure :: refactor => hpd_refactor
        procedure, private :: solve_vector => hpd_solve_vector
        procedure, private :: solve_block => hpd_solve_block
        generic :: solve => solve_vector, solve_block
        procedure :: solve_inplace => hpd_solve_inplace
        procedure :: free => hpd_free
        procedure :: dimension => hpd_dimension
        procedure :: is_factored => hpd_is_factored
        procedure :: factor_nnz => hpd_factor_nnz
        procedure :: memory_current_bytes => hpd_memory_current
        procedure :: memory_peak_bytes => hpd_memory_peak
        procedure :: memory_payloads => hpd_memory_payloads
        procedure :: memory_failure_bytes => hpd_memory_failure
        final :: hpd_finalize
    end type native_hpd_t
contains
    subroutine hpd_factor(self, A, status, memory_limit_bytes, low_memory)
        class(native_hpd_t), intent(inout) :: self
        type(csc_z_t), intent(in) :: A
        type(fortsparse_status_t), intent(out) :: status
        integer(i8), intent(in), optional :: memory_limit_bytes
        logical, intent(in), optional :: low_memory
        integer :: j, p, i, q, top, pos
        integer(i8) :: nz

        call self%free()
        self%budget%peak = 0_i8
        self%budget%failed_required = 0_i8
        self%budget%failed_available = 0_i8
        if (present(memory_limit_bytes)) self%budget%limit = memory_limit_bytes
        if (present(low_memory)) self%low_memory = low_memory
        call self%budget%reserve(0_i8, status, 'HPD memory configuration')
        if (status%code /= FORTSPARSE_OK) return
        call validate_hpd(A, status)
        if (status%code /= FORTSPARSE_OK) return
        self%n = A%ncol
        ! Every allocated payload is reserved separately before allocation.
        call allocate_integer(self%budget, self%ap, self%n + 1, status)
        if (status%code /= FORTSPARSE_OK) goto 900
        call allocate_integer(self%budget, self%ai, A%nnz, status)
        if (status%code /= FORTSPARSE_OK) goto 900
        call allocate_integer(self%budget, self%lp, self%n + 1, status)
        if (status%code /= FORTSPARSE_OK) goto 900
        call allocate_integer(self%budget, self%parent, self%n, status)
        if (status%code /= FORTSPARSE_OK) goto 900
        call allocate_integer(self%budget, self%count, self%n, status)
        if (status%code /= FORTSPARSE_OK) goto 900
        call allocate_integer(self%budget, self%mark, self%n, status)
        if (status%code /= FORTSPARSE_OK) goto 900
        call allocate_integer(self%budget, self%reach, self%n, status)
        if (status%code /= FORTSPARSE_OK) goto 900
        call allocate_integer(self%budget, self%path, self%n, status)
        if (status%code /= FORTSPARSE_OK) goto 900
        self%ap = A%col_ptr
        self%ai = A%row_idx
        self%parent = 0
        self%count = 0
        self%mark = 0
        ! Structural fill count, independent of all numerical values.
        do j = 1, self%n
            do p = self%ap(j), self%ap(j + 1) - 1
                i = self%ai(p)
                if (i >= j) exit
                do while (i /= 0)
                    if (self%mark(i) == j) exit
                    self%mark(i) = j
                    self%count(i) = self%count(i) + 1
                    if (self%parent(i) == 0) self%parent(i) = j
                    i = self%parent(i)
                    if (i == j) exit
                end do
            end do
        end do
        nz = int(self%n, i8)
        do i = 1, self%n
            nz = nz + int(self%count(i), i8)
        end do
        if (nz >= int(huge(0), i8)) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                            'native HPD: factor exceeds supported CSC index range')
            goto 900
        end if
        self%lp(1) = 1
        do j = 1, self%n
            self%lp(j + 1) = self%lp(j) + self%count(j) + 1
        end do
        call allocate_integer(self%budget, self%li, int(nz), status)
        if (status%code /= FORTSPARSE_OK) goto 900
        self%mark = 0
        self%count = 1
        do j = 1, self%n
            self%li(self%lp(j)) = j
            call row_reach(self, j, top)
            do q = top, self%n
                i = self%reach(q)
                pos = self%lp(i) + self%count(i)
                self%li(pos) = j
                self%count(i) = self%count(i) + 1
            end do
        end do
        call allocate_complex(self%budget, self%lv, int(nz), status)
        if (status%code /= FORTSPARSE_OK) goto 900
        call allocate_complex(self%budget, self%work, self%n, status)
        if (status%code /= FORTSPARSE_OK) goto 900
        call numeric_factor(self, A, status)
        return
900     call self%free()
    end subroutine hpd_factor

    subroutine hpd_refactor(self, A, status, reused, memory_limit_bytes, low_memory)
        class(native_hpd_t), intent(inout) :: self
        type(csc_z_t), intent(in) :: A
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(out), optional :: reused
        integer(i8), intent(in), optional :: memory_limit_bytes
        logical, intent(in), optional :: low_memory
        logical :: same, retained_low_memory
        integer :: p
        integer(i8) :: retained_limit

        if (present(reused)) reused = .false.
        self%ready = .false.
        if (present(memory_limit_bytes)) then
            if (memory_limit_bytes > 0_i8) then
                if (self%budget%current > memory_limit_bytes) then
                    self%budget%failed_required = self%budget%current
                    self%budget%failed_available = memory_limit_bytes
                    call status_set(status, FORTSPARSE_RESOURCE_LIMIT, &
                                    'native HPD: new cap is below retained payload')
                    return
                end if
            end if
            self%budget%limit = memory_limit_bytes
        end if
        if (present(low_memory)) self%low_memory = low_memory
        call self%budget%reserve(0_i8, status, 'HPD retained factor payload')
        if (status%code /= FORTSPARSE_OK) return
        call validate_hpd(A, status)
        if (status%code /= FORTSPARSE_OK) return
        same = .false.
        if (allocated(self%ap)) then
            if (self%n == A%ncol .and. size(self%ai) == A%nnz) then
                same = .true.
                do p = 1, self%n + 1
                    if (self%ap(p) /= A%col_ptr(p)) same = .false.
                end do
                if (same) then
                    do p = 1, A%nnz
                        if (self%ai(p) /= A%row_idx(p)) same = .false.
                    end do
                end if
            end if
        end if
        if (same) then
            if (present(reused)) reused = .true.
            call numeric_factor(self, A, status)
        else
            ! Preserve the selected policy explicitly across a new analysis;
            ! this call must not depend on free's configuration-reset semantics.
            retained_limit = self%budget%limit
            retained_low_memory = self%low_memory
            call self%factor(A, status, memory_limit_bytes=retained_limit, &
                             low_memory=retained_low_memory)
        end if
    end subroutine hpd_refactor

    subroutine validate_hpd(A, status)
        type(csc_z_t), intent(in) :: A
        type(fortsparse_status_t), intent(out) :: status
        integer :: n, j, p, i, lo, hi, mid, previous, paired

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                       'native HPD: require finite full Hermitian canonical square CSC')
        if (A%nrow < 0 .or. A%nrow /= A%ncol .or. A%nnz < 0) return
        n = A%ncol
        if (n == huge(n) .or. A%nnz == huge(n)) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                            'native HPD: CSC endpoint overflow')
            return
        end if
        if (.not. allocated(A%col_ptr)) return
        if (.not. allocated(A%row_idx)) return
        if (.not. allocated(A%val)) return
        if (lbound(A%col_ptr, 1) /= 1) return
        if (lbound(A%row_idx, 1) /= 1 .or. lbound(A%val, 1) /= 1) return
        if (size(A%col_ptr) /= n + 1) return
        if (size(A%row_idx) /= A%nnz .or. size(A%val) /= A%nnz) return
        if (A%col_ptr(1) /= 1 .or. A%col_ptr(n + 1) /= A%nnz + 1) return
        do j = 1, n
            if (A%col_ptr(j) < 1) return
            if (A%col_ptr(j + 1) < A%col_ptr(j)) return
            if (A%col_ptr(j + 1) > A%nnz + 1) return
            previous = 0
            do p = A%col_ptr(j), A%col_ptr(j + 1) - 1
                i = A%row_idx(p)
                if (i <= previous .or. i > n) return
                previous = i
                if (.not. finite_complex(A%val(p))) return
                if (i == j) then
                    if (aimag(A%val(p)) /= 0.0_dp) return
                end if
            end do
        end do
        do j = 1, n
            do p = A%col_ptr(j), A%col_ptr(j + 1) - 1
                i = A%row_idx(p)
                lo = A%col_ptr(i)
                hi = A%col_ptr(i + 1) - 1
                paired = 0
                do while (lo <= hi)
                    mid = lo + (hi - lo)/2
                    if (A%row_idx(mid) < j) then
                        lo = mid + 1
                    else if (A%row_idx(mid) > j) then
                        hi = mid - 1
                    else
                        paired = mid
                        exit
                    end if
                end do
                if (paired == 0) return
                if (A%val(p) /= conjg(A%val(paired))) return
            end do
        end do
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine validate_hpd

    subroutine row_reach(self, j, top)
        class(native_hpd_t), intent(inout) :: self
        integer, intent(in) :: j
        integer, intent(out) :: top
        integer :: p, i, length

        top = self%n + 1
        self%mark(j) = j
        do p = self%ap(j), self%ap(j + 1) - 1
            i = self%ai(p)
            if (i >= j) exit
            length = 0
            do while (self%mark(i) /= j)
                length = length + 1
                self%path(length) = i
                self%mark(i) = j
                i = self%parent(i)
            end do
            do while (length > 0)
                top = top - 1
                self%reach(top) = self%path(length)
                length = length - 1
            end do
        end do
    end subroutine row_reach

    subroutine numeric_factor(self, A, status)
        class(native_hpd_t), intent(inout) :: self
        type(csc_z_t), intent(in) :: A
        type(fortsparse_status_t), intent(out) :: status
        integer :: j, p, q, i, top, k
        real(dp) :: d
        complex(dp) :: v

        self%ready = .false.
        self%work = cmplx(0.0_dp, 0.0_dp, dp)
        self%mark = 0
        self%count = 1
        do j = 1, self%n
            call row_reach(self, j, top)
            do p = self%ap(j), self%ap(j + 1) - 1
                i = self%ai(p)
                if (i > j) exit
                self%work(i) = A%val(p)
            end do
            d = real(self%work(j), dp)
            self%work(j) = cmplx(0.0_dp, 0.0_dp, dp)
            do q = top, self%n
                i = self%reach(q)
                v = hpd_divide(self%work(i), real(self%lv(self%lp(i)), dp))
                self%work(i) = cmplx(0.0_dp, 0.0_dp, dp)
                !$omp simd private(k)
                do p = self%lp(i) + 1, self%lp(i) + self%count(i) - 1
                    k = self%li(p)
                    self%work(k) = hpd_update(self%work(k), self%lv(p), v)
                end do
                d = hpd_norm_update(d, v)
                self%lv(self%lp(i) + self%count(i)) = conjg(v)
                self%count(i) = self%count(i) + 1
            end do
            if (.not. ieee_is_finite(d) .or. d <= 0.0_dp) then
                call status_set(status, FORTSPARSE_NOT_POSITIVE_DEFINITE, &
                                'native HPD: nonpositive or nonfinite Cholesky pivot')
                return
            end if
            self%lv(self%lp(j)) = cmplx(hpd_root(d), 0.0_dp, dp)
        end do
        do p = 1, size(self%lv)
            if (.not. finite_complex(self%lv(p))) then
                call status_set(status, FORTSPARSE_NOT_POSITIVE_DEFINITE, &
                                'native HPD: nonfinite Cholesky factor')
                return
            end if
        end do
        self%ready = .true.
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine numeric_factor

    subroutine hpd_solve_vector(self, b, x, status, operation)
        class(native_hpd_t), intent(in) :: self
        complex(dp), intent(in) :: b(:)
        complex(dp), intent(out) :: x(:)
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: operation
        integer :: op

        op = 0
        if (present(operation)) op = operation
        call solve_shape(self, size(b), size(x), op, status)
        if (status%code /= FORTSPARSE_OK) return
        x = b
        call self%solve_inplace(x, status, op)
    end subroutine hpd_solve_vector

    subroutine hpd_solve_inplace(self, b, status, operation)
        class(native_hpd_t), intent(in) :: self
        complex(dp), intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: operation
        integer :: op, j

        op = 0
        if (present(operation)) op = operation
        call solve_shape(self, size(b), size(b), op, status)
        if (status%code /= FORTSPARSE_OK) return
        do j = 1, self%n
            if (.not. finite_complex(b(j))) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                                'native HPD: nonfinite right-hand side')
                return
            end if
        end do
        call triangular_solve(self, b, op)
        do j = 1, self%n
            if (.not. finite_complex(b(j))) then
                call status_set(status, FORTSPARSE_INTERNAL_ERROR, &
                                'native HPD: nonfinite solution')
                return
            end if
        end do
    end subroutine hpd_solve_inplace

    subroutine hpd_solve_block(self, B, X, status, operation)
        class(native_hpd_t), intent(in) :: self
        complex(dp), intent(in) :: B(:, :)
        complex(dp), intent(out) :: X(:, :)
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: operation
        integer :: op, j, rhs

        op = 0
        if (present(operation)) op = operation
        call solve_shape(self, size(B, 1), size(X, 1), op, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(B, 2) /= size(X, 2)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                            'native HPD: block RHS count mismatch')
            return
        end if
        do rhs = 1, size(B, 2)
            do j = 1, self%n
                if (.not. finite_complex(B(j, rhs))) then
                    call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                                    'native HPD: nonfinite block right-hand side')
                    return
                end if
            end do
        end do
        ! Only caller-owned independent output columns are modified.
        !$omp parallel do default(none) shared(self,B,X,op) private(rhs) &
        !$omp if(size(B,2) >= 8)
        do rhs = 1, size(B, 2)
            X(:, rhs) = B(:, rhs)
            call triangular_solve(self, X(:, rhs), op)
        end do
        !$omp end parallel do
        do rhs = 1, size(X, 2)
            do j = 1, self%n
                if (.not. finite_complex(X(j, rhs))) then
                    call status_set(status, FORTSPARSE_INTERNAL_ERROR, &
                                    'native HPD: nonfinite block solution')
                    return
                end if
            end do
        end do
    end subroutine hpd_solve_block

    subroutine solve_shape(self, nb, nx, operation, status)
        class(native_hpd_t), intent(in) :: self
        integer, intent(in) :: nb, nx, operation
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_NOT_FACTORED, &
                        'native HPD: no successful factorization')
        if (.not. self%ready) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        'native HPD: incompatible RHS shape or operation')
        if (nb /= self%n .or. nx /= self%n) return
        if (operation < 0 .or. operation > 2) return
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine solve_shape

    subroutine triangular_solve(self, x, operation)
        class(native_hpd_t), intent(in) :: self
        complex(dp), intent(inout) :: x(:)
        integer, intent(in) :: operation
        integer :: j, p, i
        complex(dp) :: v

        ! A^T = conjg(A), hence x = conjg(A^{-1} conjg(b)). A^H = A.
        if (operation == 1) x = conjg(x)
        do j = 1, self%n
            v = hpd_divide(x(j), real(self%lv(self%lp(j)), dp))
            x(j) = v
            !$omp simd private(i)
            do p = self%lp(j) + 1, self%lp(j + 1) - 1
                i = self%li(p)
                x(i) = hpd_update(x(i), self%lv(p), v)
            end do
        end do
        do j = self%n, 1, -1
            v = x(j)
            do p = self%lp(j) + 1, self%lp(j + 1) - 1
                v = hpd_conj_update(v, x(self%li(p)), self%lv(p))
            end do
            x(j) = hpd_divide(v, real(self%lv(self%lp(j)), dp))
        end do
        if (operation == 1) x = conjg(x)
    end subroutine triangular_solve

    pure logical function finite_complex(z) result(ok)
        complex(dp), intent(in) :: z
        ok = ieee_is_finite(real(z, dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex

    subroutine allocate_integer(budget, array, n, status)
        type(memory_budget_t), intent(inout) :: budget
        integer, allocatable, intent(out) :: array(:)
        integer, intent(in) :: n
        type(fortsparse_status_t), intent(out) :: status
        type(memory_budget_t) :: before
        integer :: ierr

        before = budget
        call budget%reserve(ib*int(n, i8), status, 'HPD integer payload')
        if (status%code /= FORTSPARSE_OK) return
        allocate (array(n), stat=ierr)
        if (ierr /= 0) then
            budget = before
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, &
                            'native HPD: integer allocation failed')
        end if
    end subroutine allocate_integer

    subroutine allocate_complex(budget, array, n, status)
        type(memory_budget_t), intent(inout) :: budget
        complex(dp), allocatable, intent(out) :: array(:)
        integer, intent(in) :: n
        type(fortsparse_status_t), intent(out) :: status
        type(memory_budget_t) :: before
        integer :: ierr

        before = budget
        call budget%reserve(zb*int(n, i8), status, 'HPD complex payload')
        if (status%code /= FORTSPARSE_OK) return
        allocate (array(n), stat=ierr)
        if (ierr /= 0) then
            budget = before
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, &
                            'native HPD: complex allocation failed')
        end if
    end subroutine allocate_complex

    subroutine hpd_memory_payloads(self, factor_bytes, persistent_bytes, &
                                   workspace_bytes)
        class(native_hpd_t), intent(in) :: self
        integer(i8), intent(out) :: factor_bytes, persistent_bytes, workspace_bytes
        factor_bytes = 0_i8
        persistent_bytes = 0_i8
        workspace_bytes = 0_i8
        if (allocated(self%lv)) factor_bytes = zb*size(self%lv, kind=i8)
        if (allocated(self%li)) persistent_bytes = ib*size(self%li, kind=i8)
        if (allocated(self%ap)) persistent_bytes = persistent_bytes + &
            ib*size(self%ap, kind=i8)
        if (allocated(self%ai)) persistent_bytes = persistent_bytes + &
            ib*size(self%ai, kind=i8)
        if (allocated(self%lp)) persistent_bytes = persistent_bytes + &
            ib*size(self%lp, kind=i8)
        if (allocated(self%parent)) persistent_bytes = persistent_bytes + &
            ib*size(self%parent, kind=i8)
        if (allocated(self%count)) workspace_bytes = workspace_bytes + &
            ib*size(self%count, kind=i8)
        if (allocated(self%mark)) workspace_bytes = workspace_bytes + &
            ib*size(self%mark, kind=i8)
        if (allocated(self%reach)) workspace_bytes = workspace_bytes + &
            ib*size(self%reach, kind=i8)
        if (allocated(self%path)) workspace_bytes = workspace_bytes + &
            ib*size(self%path, kind=i8)
        if (allocated(self%work)) workspace_bytes = workspace_bytes + &
            zb*size(self%work, kind=i8)
    end subroutine hpd_memory_payloads

    subroutine hpd_memory_failure(self, required_bytes, available_bytes)
        ! Last denied reservation: increment versus remaining capacity. For a
        ! cap tightened below live payload: retained total versus proposed cap.
        class(native_hpd_t), intent(in) :: self
        integer(i8), intent(out) :: required_bytes, available_bytes
        required_bytes = self%budget%failed_required
        available_bytes = self%budget%failed_available
    end subroutine hpd_memory_failure

    integer(i8) function hpd_memory_current(self) result(bytes)
        class(native_hpd_t), intent(in) :: self
        bytes = self%budget%current
    end function hpd_memory_current

    integer(i8) function hpd_memory_peak(self) result(bytes)
        class(native_hpd_t), intent(in) :: self
        bytes = self%budget%peak
    end function hpd_memory_peak

    integer function hpd_dimension(self) result(n)
        class(native_hpd_t), intent(in) :: self
        n = self%n
    end function hpd_dimension

    logical function hpd_is_factored(self) result(ready)
        class(native_hpd_t), intent(in) :: self
        ready = self%ready
    end function hpd_is_factored

    integer(i8) function hpd_factor_nnz(self) result(nz)
        class(native_hpd_t), intent(in) :: self
        nz = 0_i8
        if (allocated(self%lv)) nz = size(self%lv, kind=i8)
    end function hpd_factor_nnz

    subroutine hpd_free(self)
        class(native_hpd_t), intent(inout) :: self
        if (allocated(self%ap)) deallocate (self%ap)
        if (allocated(self%ai)) deallocate (self%ai)
        if (allocated(self%lp)) deallocate (self%lp)
        if (allocated(self%li)) deallocate (self%li)
        if (allocated(self%parent)) deallocate (self%parent)
        if (allocated(self%count)) deallocate (self%count)
        if (allocated(self%mark)) deallocate (self%mark)
        if (allocated(self%reach)) deallocate (self%reach)
        if (allocated(self%path)) deallocate (self%path)
        if (allocated(self%lv)) deallocate (self%lv)
        if (allocated(self%work)) deallocate (self%work)
        self%n = 0
        self%ready = .false.
        self%budget%current = 0_i8
    end subroutine hpd_free

    subroutine hpd_finalize(self)
        type(native_hpd_t), intent(inout) :: self
        call self%free()
    end subroutine hpd_finalize
end module fortsparse_native_hpd
