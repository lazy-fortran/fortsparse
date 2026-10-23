module fortsparse_precision_csr
    !! Reentrant compressed CSR storage, independent of factorization backends.
    use, intrinsic :: iso_fortran_env, only: real32, real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_precision, only: FP64, FP32, precision_policy_t, &
        precision_pair_t, STAGE_PRECONDITION
    use fortnum_precision_vector, only: precision_vector_t
    use fortnum_refinement, only: mixed_refine, refinement_solver_t
    use fortnum_krylov, only: real_fgmres, real_linear_operator_t
    implicit none
    private
    type, public :: precision_csr_t
        private
        integer :: n=0, storage_format=FP64
        real(real64) :: value_scale=1
        real(real64), allocatable :: wide(:)
        real(real32), allocatable :: single(:)
        integer, allocatable :: row(:), column(:)
        type(precision_vector_t) :: values
    contains
        procedure :: store => csr_store
        procedure :: matvec => csr_matvec
        procedure :: bytes => csr_bytes
    end type
    public :: csr_gmres_ir
    type, extends(real_linear_operator_t) :: exact_csr_operator_t
        type(precision_csr_t), pointer :: matrix => null()
    contains
        procedure :: apply => exact_csr_apply
    end type
    type, extends(real_linear_operator_t) :: low_csr_operator_t
        type(precision_csr_t) :: matrix
        integer :: compute = FP64
        logical :: failed = .false.
    contains
        procedure :: apply => low_csr_apply
    end type
    type, extends(real_linear_operator_t) :: identity_operator_t
    contains
        procedure :: apply => identity_apply
    end type
    type, extends(refinement_solver_t) :: csr_correction_t
        type(precision_csr_t), pointer :: exact => null()
        type(low_csr_operator_t) :: low
        real(real64), allocatable :: original(:)
        integer :: active_storage = -1, inner_iterations, restart
    contains
        procedure :: solve => csr_correction_solve
    end type

contains
    subroutine csr_store(self, n, row, column, values, storage, ok)
        class(precision_csr_t), intent(inout) :: self
        integer, intent(in) :: n, row(:), column(:), storage
        real(real64), intent(in) :: values(:)
        logical, intent(out) :: ok
        type(precision_vector_t) :: next, empty
        ok = .false.
        if (n < 1 .or. size(row) /= n+1) return
        if (size(column) /= size(values)) return
        if (row(1) /= 1 .or. row(n+1) /= size(values)+1) return
        if (any(row(2:) < row(:n))) return
        if (any(row < 1) .or. any(row > size(values)+1)) return
        if (any(column < 1) .or. any(column > n)) return
        call next%store(values, storage, ok)
        if (.not. ok) return
        self%n = n
        self%row = row
        self%column = column
        self%storage_format = storage
        self%value_scale = next%scale_factor()
        if (allocated(self%wide)) deallocate(self%wide)
        if (allocated(self%single)) deallocate(self%single)
        self%values = empty
        select case (storage)
        case (FP64)
            self%wide = values
        case (FP32)
            self%single = real(values/self%value_scale, real32)
        case default
            self%values = next
        end select
    end subroutine

    subroutine csr_matvec(self, x, y, compute, ok)
        class(precision_csr_t), intent(in) :: self
        real(real64), intent(in) :: x(:)
        real(real64), intent(out) :: y(:)
        integer, intent(in) :: compute
        logical, intent(out) :: ok
        real(real64) :: sum64, xs, vs
        real(real32) :: sum32
        integer :: i, j
        ok = .false.
        if (self%n == 0 .or. size(x) /= self%n .or. size(y) /= self%n) return
        if (compute /= FP32 .and. compute /= FP64) return
        if (.not. all(ieee_is_finite(x))) return
        xs = maxval(abs(x))
        if (xs == 0) xs = 1
        vs = self%value_scale
        if (compute == FP32 .and. self%storage_format == FP64) then
            do j=1,size(self%wide)
                if (abs(self%wide(j)) > real(huge(0.0_real32),real64)) return
            end do
        end if
        ! Dispatch storage outside the indexed hot loop: no virtual accessor
        ! or format branch per nonzero. SIMD may reassociate the reduction.
        if (compute == FP64 .and. self%storage_format == FP64) then
            !$omp parallel do private(j,sum64) if(self%n >= 1024)
            do i=1,self%n
                sum64 = 0
                !$omp simd reduction(+:sum64)
                do j=self%row(i), self%row(i+1)-1
                    sum64 = sum64 + self%wide(j)*x(self%column(j))
                end do
                y(i) = sum64
            end do
        else if (compute == FP64 .and. self%storage_format == FP32) then
            !$omp parallel do private(j,sum64) if(self%n >= 1024)
            do i=1,self%n
                sum64 = 0
                !$omp simd reduction(+:sum64)
                do j=self%row(i), self%row(i+1)-1
                    sum64 = sum64 + real(self%single(j),real64)*x(self%column(j))
                end do
                y(i) = sum64*vs
            end do
        else if (compute == FP32 .and. self%storage_format == FP32) then
            !$omp parallel do private(j,sum32) if(self%n >= 1024)
            do i=1,self%n
                sum32 = 0
                !$omp simd reduction(+:sum32)
                do j=self%row(i), self%row(i+1)-1
                    sum32 = sum32 + self%single(j)*real(x(self%column(j))/xs,real32)
                end do
                y(i) = (real(sum32,real64)*vs)*xs
            end do
        else if (compute == FP32 .and. self%storage_format == FP64) then
            !$omp parallel do private(j,sum32) if(self%n >= 1024)
            do i=1,self%n
                sum32 = 0
                !$omp simd reduction(+:sum32)
                do j=self%row(i), self%row(i+1)-1
                    sum32 = sum32 + real(self%wide(j),real32)* &
                        real(x(self%column(j))/xs,real32)
                end do
                y(i) = real(sum32,real64)*xs
            end do
        else if (compute == FP64) then
            !$omp parallel do private(j,sum64) if(self%n >= 1024)
            do i=1,self%n
                sum64 = 0
                do j=self%row(i), self%row(i+1)-1
                    sum64 = sum64 + self%values%at(j)*x(self%column(j))
                end do
                y(i) = sum64
            end do
        else
            !$omp parallel do private(j,sum32) if(self%n >= 1024)
            do i=1,self%n
                sum32 = 0
                do j=self%row(i), self%row(i+1)-1
                    sum32 = sum32 + real(self%values%at_scaled(j), real32) * &
                        real(x(self%column(j))/xs, real32)
                end do
                y(i) = (real(sum32, real64)*vs)*xs
            end do
        end if
        ok = all(ieee_is_finite(y))
    end subroutine

    pure integer(int64) function csr_bytes(self) result(n)
        class(precision_csr_t), intent(in) :: self
        n = self%values%bytes()
        if (allocated(self%wide)) n = n + 8_int64*size(self%wide,kind=int64)
        if (allocated(self%single)) n = n + 4_int64*size(self%single,kind=int64)
        if (allocated(self%row)) n = n + int(size(self%row),int64)* &
            storage_size(0)/8
        if (allocated(self%column)) n = n + int(size(self%column),int64)* &
            storage_size(0)/8
    end function

    subroutine csr_gmres_ir(a, b, x, policy, tolerance, max_steps, &
            inner_iterations, restart, info, steps, residual)
        !! Exact matrix must be fp64. Rounded copies are rebuilt after escalation.
        !! Reuses FortNum FGMRES: basis/orthogonalization fp64, operator selectable.
        type(precision_csr_t), intent(in), target :: a
        real(real64), intent(in) :: b(:), tolerance
        real(real64), intent(inout) :: x(:)
        type(precision_policy_t), intent(inout) :: policy
        integer, intent(in) :: max_steps, inner_iterations, restart
        integer, intent(out) :: info, steps
        real(real64), intent(out) :: residual
        type(exact_csr_operator_t) :: exact
        type(csr_correction_t) :: correction
        info = -1
        steps = 0
        residual = huge(1.0_real64)
        if (a%storage_format /= FP64) return
        if (a%n < 1 .or. inner_iterations < 1) return
        if (restart < 1 .or. restart > inner_iterations) return
        exact%matrix => a
        correction%exact => a
        correction%original = a%wide
        correction%active_storage = -1
        correction%inner_iterations = inner_iterations
        correction%restart = restart
        call mixed_refine(exact, correction, b, x, policy, tolerance, &
            max_steps, info, steps, residual)
    end subroutine

    subroutine exact_csr_apply(self, input, output)
        class(exact_csr_operator_t), intent(inout) :: self
        real(real64), intent(in) :: input(:)
        real(real64), intent(out) :: output(:)
        logical :: ok
        call self%matrix%matvec(input, output, FP64, ok)
    end subroutine

    subroutine low_csr_apply(self, input, output)
        class(low_csr_operator_t), intent(inout) :: self
        real(real64), intent(in) :: input(:)
        real(real64), intent(out) :: output(:)
        logical :: ok
        call self%matrix%matvec(input, output, self%compute, ok)
        if (.not. ok) self%failed = .true.
    end subroutine

    subroutine identity_apply(self, input, output)
        class(identity_operator_t), intent(inout) :: self
        real(real64), intent(in) :: input(:)
        real(real64), intent(out) :: output(:)
        output = input
    end subroutine

    subroutine csr_correction_solve(self, rhs, delta, precision, info)
        class(csr_correction_t), intent(inout) :: self
        real(real64), intent(in) :: rhs(:)
        real(real64), intent(out) :: delta(:)
        type(precision_pair_t), intent(in) :: precision
        integer, intent(out) :: info
        type(identity_operator_t) :: identity
        integer :: iterations
        real(real64) :: norm
        logical :: ok
        if (precision%storage /= self%active_storage) then
            call self%low%matrix%store(self%exact%n, self%exact%row, &
                self%exact%column, self%original, precision%storage, ok)
            if (.not. ok) then
                info = -1
                return
            end if
            self%active_storage = precision%storage
        end if
        self%low%compute = precision%compute
        self%low%failed = .false.
        delta = 0
        call real_fgmres(self%low, identity, rhs, delta, 0.05_real64, &
            self%inner_iterations, self%restart, info, iterations, norm)
        ! The independently recomputed outer defect decides acceptance.
        if (info >= 0 .and. .not. self%low%failed) info = 0
        if (self%low%failed) info = -1
    end subroutine
end module fortsparse_precision_csr
