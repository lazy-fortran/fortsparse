module fortsparse_native_hermitian_embedding
    ! Complex Hermitian A = Ar + i Ai through the real symmetric embedding
    ! E = [[Ar, -Ai], [Ai, Ar]] (stored interleaved) and the native supernodal
    ! LDL^T. E is the
    ! realification of A: every eigenvalue of A appears twice in E, so the
    ! inertia of E is twice that of A, and E [xr; xi] = [br; bi] solves
    ! A x = b. Used when complex Hermitian or Hermitian positive definite
    ! factors request ldlt_supernodal, until native complex supernodal
    ! factors exist; a positive definite A gives a positive definite E and the
    ! supernodal Cholesky.
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csc, only: csc_t, csc_z_t
    use fortsparse_status, only: fortsparse_status_t, status_set, status_ok, &
        FORTSPARSE_OK, FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_INVALID_MATRIX, &
        FORTSPARSE_NOT_FACTORED
    use fortsparse_backend, only: sparse_backend_t, sparse_factor_stats_t, &
        FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_SYM_SPD
    use fortsparse_inertia, only: FORTSPARSE_INERTIA_UNAVAILABLE
    use fortsparse_native, only: native_backend_t
    implicit none
    private

    type, extends(sparse_backend_t), public :: native_hermitian_embedding_t
        type(native_backend_t) :: inner
        integer :: n = 0
        logical :: ready = .false.
        real(dp), allocatable :: rhs(:), sol(:)
        real(dp), pointer :: scratch(:) => null()
    contains
        procedure :: factor_real => embedding_factor_real
        procedure :: factor_complex => embedding_factor_complex
        procedure :: solve_real => embedding_solve_real
        procedure :: solve_complex => embedding_solve_complex
        procedure :: solve_real_inplace => embedding_solve_real_inplace
        procedure :: solve_complex_inplace => embedding_solve_complex_inplace
        procedure :: factor_stats => embedding_factor_stats
        procedure :: vector => embedding_vector
        procedure :: free => embedding_free
    end type native_hermitian_embedding_t

contains

    subroutine embedding_factor_real(self, A, refine, status)
        class(native_hermitian_embedding_t), intent(inout) :: self
        type(csc_t), intent(in) :: A
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
            'Hermitian embedding: complex matrices only')
    end subroutine embedding_factor_real

    subroutine embedding_factor_complex(self, A, refine, status)
        class(native_hermitian_embedding_t), intent(inout) :: self
        type(csc_z_t), intent(in) :: A
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        type(csc_t) :: e
        integer :: n, j, k, at
        self%ready = .false.
        if (A%nrow /= A%ncol .or. A%ncol < 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'Hermitian embedding: square matrix required')
            return
        end if
        n = A%ncol
        e%nrow = 2*n; e%ncol = 2*n; e%nnz = 4*A%nnz
        allocate(e%col_ptr(2*n + 1), e%row_idx(e%nnz), e%val(e%nnz))
        ! Interleaved unknowns (Re x_1, Im x_1, Re x_2, ...): each entry a_ij is
        ! the 2x2 block [[Re a, -Im a], [Im a, Re a]], so E keeps the pattern
        ! and the band of A.
        at = 0
        e%col_ptr(1) = 1
        do j = 1, n
            do k = A%col_ptr(j), A%col_ptr(j + 1) - 1
                at = at + 1
                e%row_idx(at) = 2*A%row_idx(k) - 1; e%val(at) = real(A%val(k), dp)
                at = at + 1
                e%row_idx(at) = 2*A%row_idx(k); e%val(at) = aimag(A%val(k))
            end do
            e%col_ptr(2*j) = at + 1
            do k = A%col_ptr(j), A%col_ptr(j + 1) - 1
                at = at + 1
                e%row_idx(at) = 2*A%row_idx(k) - 1; e%val(at) = -aimag(A%val(k))
                at = at + 1
                e%row_idx(at) = 2*A%row_idx(k); e%val(at) = real(A%val(k), dp)
            end do
            e%col_ptr(2*j + 1) = at + 1
        end do
        self%inner%options = self%options
        if (self%options%symmetry /= FORTSPARSE_SYM_SPD) then
            self%inner%options%symmetry = FORTSPARSE_SYM_INDEFINITE
            self%inner%options%ldlt_supernodal = .true.
        end if
        call self%inner%factor_real(e, refine, status)
        if (.not. status_ok(status)) return
        self%n = n
        if (allocated(self%rhs)) deallocate(self%rhs, self%sol)
        allocate(self%rhs(2*n), self%sol(2*n))
        self%ready = .true.
    end subroutine embedding_factor_complex

    subroutine embedding_solve_real(self, b, x, status)
        class(native_hermitian_embedding_t), intent(inout) :: self
        real(dp), target, contiguous, intent(in) :: b(:)
        real(dp), target, contiguous, intent(out) :: x(:)
        type(fortsparse_status_t), intent(out) :: status
        x = 0
        call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
            'Hermitian embedding: complex right-hand sides only')
    end subroutine embedding_solve_real

    subroutine embedding_solve_complex(self, b, x, status)
        class(native_hermitian_embedding_t), intent(inout) :: self
        complex(dp), intent(in) :: b(:)
        complex(dp), intent(out) :: x(:)
        type(fortsparse_status_t), intent(out) :: status
        x = 0
        if (.not. self%ready) then
            call status_set(status, FORTSPARSE_NOT_FACTORED, &
                'Hermitian embedding: not factored')
            return
        end if
        if (size(b) /= self%n .or. size(x) /= self%n) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'Hermitian embedding: right-hand side size')
            return
        end if
        self%rhs(1::2) = real(b, dp)
        self%rhs(2::2) = aimag(b)
        call self%inner%solve_real(self%rhs, self%sol, status)
        if (.not. status_ok(status)) return
        x = cmplx(self%sol(1::2), self%sol(2::2), dp)
    end subroutine embedding_solve_complex

    subroutine embedding_solve_real_inplace(self, b, status)
        class(native_hermitian_embedding_t), intent(inout) :: self
        real(dp), intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out) :: status
        call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
            'Hermitian embedding: complex right-hand sides only')
    end subroutine embedding_solve_real_inplace

    subroutine embedding_solve_complex_inplace(self, b, status)
        class(native_hermitian_embedding_t), intent(inout) :: self
        complex(dp), intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out) :: status
        complex(dp) :: x(size(b))
        call embedding_solve_complex(self, b, x, status)
        if (status_ok(status)) b = x
    end subroutine embedding_solve_complex_inplace

    subroutine embedding_factor_stats(self, stats, status)
        ! Statistics of the real factor of E, with the inertia of A.
        class(native_hermitian_embedding_t), intent(inout) :: self
        type(sparse_factor_stats_t), intent(out) :: stats
        type(fortsparse_status_t), intent(out) :: status
        call self%inner%factor_stats(stats, status)
        if (.not. status_ok(status)) return
        if (stats%inertia_status == FORTSPARSE_INERTIA_UNAVAILABLE) return
        if (mod(stats%inertia_negative, 2_i8) /= 0_i8 .or. &
            mod(stats%inertia_zero, 2_i8) /= 0_i8) then
            ! A rounding split of a double eigenvalue: no Hermitian count.
            stats%inertia_status = FORTSPARSE_INERTIA_UNAVAILABLE
            stats%inertia_negative = -1_i8; stats%inertia_zero = -1_i8
            stats%inertia_positive = -1_i8
            return
        end if
        stats%inertia_negative = stats%inertia_negative/2_i8
        stats%inertia_zero = stats%inertia_zero/2_i8
        stats%inertia_positive = stats%inertia_positive/2_i8
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine embedding_factor_stats

    function embedding_vector(self, n) result(p)
        class(native_hermitian_embedding_t), intent(inout) :: self
        integer, intent(in) :: n
        real(dp), pointer :: p(:)
        if (associated(self%scratch)) deallocate(self%scratch)
        allocate(self%scratch(n))
        p => self%scratch
    end function embedding_vector

    subroutine embedding_free(self)
        class(native_hermitian_embedding_t), intent(inout) :: self
        call self%inner%free()
        if (allocated(self%rhs)) deallocate(self%rhs, self%sol)
        if (associated(self%scratch)) deallocate(self%scratch)
        self%ready = .false.
        self%n = 0
    end subroutine embedding_free
end module fortsparse_native_hermitian_embedding
