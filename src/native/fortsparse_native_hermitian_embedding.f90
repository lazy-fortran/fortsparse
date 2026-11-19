module fortsparse_native_hermitian_embedding
    ! Complex Hermitian A = Ar + i Ai through the real symmetric embedding
    ! E = [[Ar, -Ai], [Ai, Ar]] (stored interleaved). E is the realification
    ! of A: every eigenvalue of A appears twice in E, so the inertia of E is
    ! twice that of A, and E [xr; xi] = [br; bi] solves A x = b. A Hermitian
    ! A is factored by the real multifrontal LDL^T; a positive definite A
    ! gives a positive definite E and the supernodal Cholesky.
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csc, only: csc_t, csc_z_t
    use fortsparse_status, only: fortsparse_status_t, status_set, status_ok, &
        FORTSPARSE_OK, FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_INVALID_MATRIX, &
        FORTSPARSE_NOT_FACTORED
    use fortsparse_backend, only: sparse_backend_t, sparse_factor_stats_t, &
        FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_SYM_SPD, FORTSPARSE_OP_NORMAL, &
        FORTSPARSE_OP_TRANSPOSE, FORTSPARSE_OP_ADJOINT
    use fortsparse_inertia, only: FORTSPARSE_INERTIA_UNAVAILABLE
    use fortsparse_native, only: native_backend_t
    use fortsparse_native_multifrontal, only: native_multifrontal_backend_t
    implicit none
    private

    type, extends(sparse_backend_t), public :: native_hermitian_embedding_t
        class(sparse_backend_t), allocatable :: inner
        integer :: n = 0
        logical :: ready = .false.
        real(dp), allocatable :: rhs(:), sol(:)
        real(dp), pointer :: scratch(:) => null()
    contains
        procedure :: factor_real => embedding_factor_real
        procedure :: factor_complex => embedding_factor_complex
        procedure :: refactor_complex => embedding_refactor_complex
        procedure :: solve_complex_op => embedding_solve_complex_op
        procedure :: copy_from => embedding_copy_from
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
        self%ready = .false.
        if (.not. embed(A, e, status)) return
        if (allocated(self%inner)) deallocate (self%inner)
        if (self%options%symmetry == FORTSPARSE_SYM_SPD) then
            allocate (native_backend_t :: self%inner)
        else
            allocate (native_multifrontal_backend_t :: self%inner)
        end if
        self%inner%options = self%options
        if (self%options%symmetry /= FORTSPARSE_SYM_SPD) &
            self%inner%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        call self%inner%factor_real(e, refine, status)
        call finish(self, A%ncol, status)
    end subroutine embedding_factor_complex

    subroutine embedding_refactor_complex(self, A, refine, status, reused)
        ! Same pattern: the inner factor reuses its order and analysis.
        class(native_hermitian_embedding_t), intent(inout) :: self
        type(csc_z_t), intent(in) :: A
        logical, intent(in) :: refine
        type(fortsparse_status_t), intent(out) :: status
        logical, intent(out) :: reused
        type(csc_t) :: e
        reused = .false.
        if (.not. allocated(self%inner)) then
            call self%factor_complex(A, refine, status)
            return
        end if
        if ((self%options%symmetry == FORTSPARSE_SYM_SPD) .neqv. &
            (self%inner%options%symmetry == FORTSPARSE_SYM_SPD)) then
            call self%factor_complex(A, refine, status)
            return
        end if
        self%ready = .false.
        if (.not. embed(A, e, status)) return
        call self%inner%refactor_real(e, refine, status, reused)
        call finish(self, A%ncol, status)
    end subroutine embedding_refactor_complex

    subroutine finish(self, n, status)
        class(native_hermitian_embedding_t), intent(inout) :: self
        integer, intent(in) :: n
        type(fortsparse_status_t), intent(in) :: status
        if (.not. status_ok(status)) return
        self%n = n
        if (allocated(self%rhs)) deallocate(self%rhs, self%sol)
        allocate(self%rhs(2*n), self%sol(2*n))
        self%ready = .true.
    end subroutine finish

    logical function embed(A, e, status) result(ok)
        ! Interleaved unknowns (Re x_1, Im x_1, Re x_2, ...): each entry a_ij is
        ! the 2x2 block [[Re a, -Im a], [Im a, Re a]], so E keeps the pattern
        ! and the band of A.
        type(csc_z_t), intent(in) :: A
        type(csc_t), intent(out) :: e
        type(fortsparse_status_t), intent(out) :: status
        integer :: n, j, k, at
        ok = .false.
        if (A%nrow /= A%ncol .or. A%ncol < 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'Hermitian embedding: square matrix required')
            return
        end if
        n = A%ncol
        e%nrow = 2*n; e%ncol = 2*n; e%nnz = 4*A%nnz
        allocate(e%col_ptr(2*n + 1), e%row_idx(e%nnz), e%val(e%nnz))
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
        call status_set(status, FORTSPARSE_OK, '')
        ok = .true.
    end function embed

    subroutine embedding_solve_complex_op(self, b, x, operation, status)
        ! A^H = A; A^T = conj(A), so A^T x = b is conj(A^-1 conj(b)).
        class(native_hermitian_embedding_t), intent(inout) :: self
        complex(dp), intent(in) :: b(:)
        complex(dp), intent(out) :: x(:)
        integer, intent(in) :: operation
        type(fortsparse_status_t), intent(out) :: status
        select case (operation)
        case (FORTSPARSE_OP_NORMAL, FORTSPARSE_OP_ADJOINT)
            call embedding_solve_complex(self, b, x, status)
        case (FORTSPARSE_OP_TRANSPOSE)
            call embedding_solve_complex(self, conjg(b), x, status)
            x = conjg(x)
        case default
            x = 0
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'Hermitian embedding: unknown solve operation')
        end select
    end subroutine embedding_solve_complex_op

    subroutine embedding_copy_from(self, src, status)
        class(native_hermitian_embedding_t), intent(inout) :: self
        class(sparse_backend_t), intent(in) :: src
        type(fortsparse_status_t), intent(out) :: status
        select type (src)
        type is (native_hermitian_embedding_t)
            if (allocated(self%inner)) deallocate (self%inner)
            call status_set(status, FORTSPARSE_NOT_FACTORED, &
                'Hermitian embedding: source not factored')
            if (.not. allocated(src%inner)) return
            allocate (self%inner, mold=src%inner)
            self%inner%options = src%inner%options
            self%inner%options%host_memory_limit_bytes = self%options%host_memory_limit_bytes
            call self%inner%copy_from(src%inner, status)
            if (.not. status_ok(status)) return
            self%n = src%n
            self%ready = src%ready
            if (allocated(src%rhs)) then
                self%rhs = src%rhs
                self%sol = src%sol
            end if
        class default
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'backend copy: source and destination types differ')
        end select
    end subroutine embedding_copy_from

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
        if (.not. allocated(self%inner)) then
            call status_set(status, FORTSPARSE_NOT_FACTORED, &
                'Hermitian embedding: not factored')
            return
        end if
        call self%inner%factor_stats(stats, status)
        if (.not. status_ok(status)) return
        stats%native_factor_family = 'LDLH_EMBEDDED'
        if (self%options%symmetry == FORTSPARSE_SYM_SPD) stats%native_factor_family = 'HPD_EMBEDDED'
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
        if (allocated(self%inner)) then
            call self%inner%free()
            deallocate (self%inner)
        end if
        if (allocated(self%rhs)) deallocate(self%rhs, self%sol)
        if (associated(self%scratch)) deallocate(self%scratch)
        self%ready = .false.
        self%n = 0
    end subroutine embedding_free
end module fortsparse_native_hermitian_embedding
