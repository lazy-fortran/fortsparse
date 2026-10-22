module fortsparse_smoothers
    !! Fixed symmetric polynomial corrections for real SPD CSR operators.
    use fortsparse_kinds, only: dp
    use fortsparse_csr, only: csr_t, csr_valid, csr_find, csr_matvec
    use fortsparse_status, only: fortsparse_status_t, status_set, &
                       FORTSPARSE_OK, FORTSPARSE_INVALID_MATRIX, FORTSPARSE_NOT_FACTORED
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    integer, parameter, public :: SMOOTHER_L1_JACOBI = 1, SMOOTHER_CHEBYSHEV = 2
    public :: csr_smoother_t, smoother_setup, smoother_apply, smoother_destroy
    type :: csr_smoother_t
        private
        integer :: method = 0, degree = 2
        real(dp) :: upper = 0, estimate = 0, ratio = 20
        real(dp), allocatable :: inverse(:), residual(:), direction(:)
    contains
        procedure :: eigenvalue_estimate
        procedure :: eigenvalue_bound
    end type
contains
    subroutine smoother_setup(a, method, handle, info, degree, ratio)
        type(csr_t), intent(in) :: a
        integer, intent(in) :: method
        type(csr_smoother_t), intent(inout) :: handle
        type(fortsparse_status_t), intent(out) :: info
        integer, intent(in), optional :: degree
        real(dp), intent(in), optional :: ratio
        real(dp), allocatable :: x(:), y(:)
        real(dp) :: diagonal, rowsum, bound, norm, rayleigh
        integer :: i, k, j, iter
        call smoother_destroy(handle)
        call status_set(info, FORTSPARSE_INVALID_MATRIX, 'smoother: CSR/parameters')
        if (.not. csr_valid(a)) return
        if (a%nrow /= a%ncol) return
        if (method /= SMOOTHER_L1_JACOBI .and. method /= SMOOTHER_CHEBYSHEV) return
        if (present(degree)) then
            if (degree < 1) return
            handle%degree = degree
        end if
        if (present(ratio)) then
            if (.not. ieee_is_finite(ratio)) return
            if (ratio <= 1) return
            handle%ratio = ratio
        end if
        allocate (handle%inverse(a%nrow), handle%residual(a%nrow), &
                  handle%direction(a%nrow), x(a%nrow), y(a%nrow))
        bound = 0
        !$omp parallel do private(k,j,diagonal,rowsum) reduction(max:bound)
        do i = 1, a%nrow
            k = csr_find(a, i, i)
            diagonal = 0
            if (k /= 0) diagonal = a%val(k)
            rowsum = 0
            do j = a%row_ptr(i), a%row_ptr(i + 1) - 1
                rowsum = rowsum + abs(a%val(j))
            end do
            handle%inverse(i) = 0
            if (diagonal > 0) then
                handle%inverse(i) = 1/diagonal
                bound = max(bound, rowsum/diagonal)
                if (method == SMOOTHER_L1_JACOBI) handle%inverse(i) = 1/rowsum
            end if
        end do
        if (any(handle%inverse <= 0)) then
            call smoother_destroy(handle)
            call status_set(info, FORTSPARSE_INVALID_MATRIX, &
                            'smoother: positive diagonal required')
            return
        end if
        handle%upper = bound*(1 + 32*epsilon(1.0_dp))
        ! Power iteration on D^(-1/2) A D^(-1/2), with a nonconstant seed.
        ! The Gershgorin upper bound, rather than an uncertified Ritz value,
        ! controls the polynomial interval and preserves a positive correction.
        !$omp parallel do
        do i = 1, a%nrow
            x(i) = sin(real(i, dp)*1.23456789_dp) + 0.5_dp
        end do
        rayleigh = 0
        do iter = 1, 12
            norm = sqrt(dot_product(x, x))
            x = x/norm
            y = sqrt(handle%inverse)*x
            call csr_matvec(a, y, handle%residual)
            y = sqrt(handle%inverse)*handle%residual
            rayleigh = dot_product(x, y)
            x = y
            if (dot_product(x, x) <= tiny(1.0_dp)) exit
        end do
        handle%estimate = rayleigh
        handle%method = method
        call status_set(info, FORTSPARSE_OK, '')
    end subroutine

    subroutine smoother_apply(a, handle, b, x, info)
        type(csr_t), intent(in) :: a
        type(csr_smoother_t), intent(inout) :: handle
        real(dp), intent(in) :: b(:)
        real(dp), intent(inout) :: x(:)
        type(fortsparse_status_t), intent(out) :: info
        real(dp) :: center, radius, rho, oldrho, alpha, beta
        integer :: step, i
        call status_set(info, FORTSPARSE_NOT_FACTORED, 'smoother: setup required')
        if (handle%method == 0) return
        call status_set(info, FORTSPARSE_INVALID_MATRIX, 'smoother: vector shape')
        if (size(b) /= size(handle%inverse) .or. size(x) /= size(b)) return
        if (a%nrow /= size(b) .or. a%ncol /= size(b)) return
        center = 0.5_dp*handle%upper*(1 + 1/handle%ratio)
        radius = 0.5_dp*handle%upper*(1 - 1/handle%ratio)
        oldrho = radius/center
        do step = 1, handle%degree
            call csr_matvec(a, x, handle%residual)
            alpha = 1
            beta = 0
            if (handle%method == SMOOTHER_CHEBYSHEV) then
                if (step == 1) then
                    alpha = 1/center
                else
                    rho = 1/(2*center/radius - oldrho)
                    beta = rho*oldrho
                    alpha = 2*rho/radius
                    oldrho = rho
                end if
            end if
            !$omp parallel do if(size(b) >= 512)
            do i = 1, size(b)
                handle%residual(i) = handle%inverse(i)*(b(i) - handle%residual(i))
                if (step == 1) handle%direction(i) = 0
                handle%direction(i) = alpha*handle%residual(i) &
                                      + beta*handle%direction(i)
                x(i) = x(i) + handle%direction(i)
            end do
        end do
        call status_set(info, FORTSPARSE_OK, '')
    end subroutine

    real(dp) function eigenvalue_estimate(self) result(value)
        class(csr_smoother_t), intent(in) :: self
        value = self%estimate
    end function

    real(dp) function eigenvalue_bound(self) result(value)
        class(csr_smoother_t), intent(in) :: self
        value = self%upper
    end function

    subroutine smoother_destroy(handle)
        type(csr_smoother_t), intent(inout) :: handle
        type(csr_smoother_t) :: empty
        handle = empty
    end subroutine
end module
