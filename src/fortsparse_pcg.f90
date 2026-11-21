module fortsparse_pcg
    !! Preconditioned conjugate gradients for real symmetric positive definite
    !! CSR matrices, and for semidefinite ones with a consistent right-hand side
    !! (b in the range of A). Preconditioners: Jacobi, IC(0) and the smoothed
    !! aggregation AMG of fortsparse_amg. Hestenes and Stiefel (1952); IC(0)
    !! with Manteuffel's diagonal shift (Math. Comp. 34, 473-497, 1980).
    use fortsparse_kinds, only: dp
    use fortsparse_csr, only: csr_t, csr_matvec
    use fortsparse_status, only: fortsparse_status_t, status_set, FORTSPARSE_OK, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_NOT_POSITIVE_DEFINITE, FORTSPARSE_NOT_CONVERGED
    use fortsparse_amg, only: amg_params_t, amg_handle_t, amg_setup, amg_apply, &
        amg_destroy
    use fortsparse_ams, only: ams_params_t, ams_handle_t, ams_setup, ams_apply, &
        ams_destroy
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: pcg_params_t, pcg_handle_t, pcg_setup_csr, pcg_setup, pcg_solve, &
        pcg_destroy, pcg_setup_ams
    integer, parameter, public :: PCG_NONE = 0, PCG_JACOBI = 1, PCG_IC0 = 2, &
        PCG_AMG = 3, PCG_AMS = 4

    type :: pcg_params_t
        integer :: preconditioner = PCG_JACOBI
        !! Stop when ||b - A x||_2 <= rtol ||b||_2.
        real(dp) :: rtol = 1.0e-8_dp
        integer :: max_iterations = 1000
        type(amg_params_t) :: amg
        !! PCG_AMS, set up by pcg_setup_ams with the gradient and nodal maps.
        type(ams_params_t) :: ams
    end type

    type :: pcg_handle_t
        private
        type(csr_t) :: a
        integer :: preconditioner = PCG_NONE
        real(dp), allocatable :: inverse_diagonal(:)
        !! IC(0): rows of L (lower triangle with diagonal last), and the
        !! shift that made it exist.
        type(csr_t) :: l
        real(dp) :: shift = 0
        ! Allocatable, so a handle that goes out of scope finalizes only the
        ! preconditioner it holds; gfortran mis-finalizes nested handles
        ! embedded directly.
        type(amg_handle_t), allocatable :: amg
        type(ams_handle_t), allocatable :: ams
        real(dp), allocatable :: r(:), z(:), p(:), q(:)
        logical :: ready = .false.
    end type

contains

    subroutine pcg_setup_csr(n, row_ptr, columns, values, params, handle, info)
        !! Raw-array bridge: one-based CSR with sorted columns, both triangles.
        integer, intent(in) :: n, row_ptr(:), columns(:)
        real(dp), intent(in) :: values(:)
        type(pcg_params_t), intent(in) :: params
        type(pcg_handle_t), intent(inout) :: handle
        type(fortsparse_status_t), intent(out) :: info
        type(csr_t) :: a

        a%nrow = n
        a%ncol = n
        a%row_ptr = row_ptr
        a%col = columns
        a%val = values
        call pcg_setup(a, params, handle, info)
    end subroutine pcg_setup_csr

    subroutine pcg_setup(a, params, handle, info)
        type(csr_t), intent(in) :: a
        type(pcg_params_t), intent(in) :: params
        type(pcg_handle_t), intent(inout) :: handle
        type(fortsparse_status_t), intent(out) :: info
        integer :: n

        call pcg_destroy(handle)
        call status_set(info, FORTSPARSE_INVALID_MATRIX, 'PCG: invalid CSR matrix')
        if (.not. valid(a)) return
        n = a%nrow
        handle%a = a
        handle%preconditioner = params%preconditioner
        select case (params%preconditioner)
        case (PCG_NONE)
        case (PCG_JACOBI)
            call jacobi_setup(handle, info)
        case (PCG_IC0)
            call ic0_setup(handle, info)
        case (PCG_AMG)
            allocate (handle%amg)
            call amg_setup(a, params%amg, handle%amg, info)
        case (PCG_AMS)
            call status_set(info, FORTSPARSE_INVALID_MATRIX, &
                'PCG: AMS needs pcg_setup_ams with G and Pi')
        case default
            call status_set(info, FORTSPARSE_INVALID_MATRIX, 'PCG: unknown preconditioner')
        end select
        if (params%preconditioner /= PCG_NONE .and. info%code /= FORTSPARSE_OK) then
            call pcg_destroy(handle)
            return
        end if
        allocate (handle%r(n), handle%z(n), handle%p(n), handle%q(n))
        handle%ready = .true.
        call status_set(info, FORTSPARSE_OK, '')
    end subroutine pcg_setup

    subroutine pcg_setup_ams(a, g, pi, params, handle, info)
        !! PCG with the AMS preconditioner for an H(curl) matrix A, its
        !! discrete gradient G (n1 x n0) and nodal vector interpolation Pi
        !! (n1 x 3 n0, x, y, z column blocks).
        type(csr_t), intent(in) :: a, g, pi
        type(pcg_params_t), intent(in) :: params
        type(pcg_handle_t), intent(inout) :: handle
        type(fortsparse_status_t), intent(out) :: info
        type(pcg_params_t) :: plain
        integer :: n

        plain = params
        plain%preconditioner = PCG_NONE
        call pcg_setup(a, plain, handle, info)
        if (info%code /= FORTSPARSE_OK) return
        allocate (handle%ams)
        call ams_setup(a, g, pi, params%ams, handle%ams, info)
        if (info%code /= FORTSPARSE_OK) then
            call pcg_destroy(handle)
            return
        end if
        handle%preconditioner = PCG_AMS
        n = a%nrow
    end subroutine pcg_setup_ams

    subroutine pcg_solve(handle, params, b, x, info, iterations, relative_residual)
        !! x holds the initial guess on entry and the iterate on return.
        type(pcg_handle_t), intent(inout) :: handle
        type(pcg_params_t), intent(in) :: params
        real(dp), intent(in) :: b(:)
        real(dp), intent(inout) :: x(:)
        type(fortsparse_status_t), intent(out) :: info
        integer, intent(out), optional :: iterations
        real(dp), intent(out), optional :: relative_residual
        real(dp) :: bnorm, rnorm, target, rz, rz_old, pq, alpha, beta
        integer :: k, n, i

        if (present(iterations)) iterations = 0
        if (present(relative_residual)) relative_residual = 0
        call status_set(info, FORTSPARSE_NOT_FACTORED, 'PCG: setup required')
        if (.not. handle%ready) return
        n = handle%a%nrow
        call status_set(info, FORTSPARSE_INVALID_MATRIX, 'PCG: vector shape')
        if (size(b) /= n .or. size(x) /= n) return
        bnorm = norm(b)
        if (bnorm == 0) then
            x = 0
            call status_set(info, FORTSPARSE_OK, '')
            return
        end if
        target = params%rtol*bnorm
        call csr_matvec(handle%a, x, handle%q)
        handle%r = b - handle%q
        rnorm = norm(handle%r)
        k = 0
        if (rnorm > target) then
            call precondition(handle, handle%r, handle%z, info)
            if (info%code /= FORTSPARSE_OK) return
            handle%p = handle%z
            rz = dot(handle%r, handle%z)
            do k = 1, params%max_iterations
                call csr_matvec(handle%a, handle%p, handle%q)
                pq = dot(handle%p, handle%q)
                ! p in the kernel of A: only a converged residual may end here.
                if (.not. (pq > 0) .or. .not. (rz > 0)) then
                    call status_set(info, FORTSPARSE_NOT_POSITIVE_DEFINITE, &
                        'PCG: curvature p^T A p or r^T z is not positive')
                    call report()
                    return
                end if
                alpha = rz/pq
                call axpy(alpha, handle%p, x)
                call axpy(-alpha, handle%q, handle%r)
                rnorm = norm(handle%r)
                if (.not. ieee_is_finite(rnorm)) then
                    call status_set(info, FORTSPARSE_NOT_POSITIVE_DEFINITE, &
                        'PCG: the residual is not finite')
                    call report()
                    return
                end if
                if (rnorm <= target) exit
                call precondition(handle, handle%r, handle%z, info)
                if (info%code /= FORTSPARSE_OK) return
                rz_old = rz
                rz = dot(handle%r, handle%z)
                beta = rz/rz_old
                !$omp parallel do if(n >= 4096)
                do i = 1, n
                    handle%p(i) = handle%z(i) + beta*handle%p(i)
                end do
            end do
            k = min(k, params%max_iterations)
        end if
        call report()
        if (rnorm <= target) then
            call status_set(info, FORTSPARSE_OK, '')
        else
            call status_set(info, FORTSPARSE_NOT_CONVERGED, &
                'PCG: iteration limit reached above the tolerance')
        end if
    contains
        subroutine report()
            if (present(iterations)) iterations = k
            if (present(relative_residual)) relative_residual = rnorm/bnorm
        end subroutine report
    end subroutine pcg_solve

    subroutine pcg_destroy(handle)
        type(pcg_handle_t), intent(inout) :: handle

        ! Released component by component: assigning a fresh handle runs the
        ! finalizers of a constructor temporary, which gfortran gets wrong
        ! for the nested solver handles.
        if (allocated(handle%amg)) then
            call amg_destroy(handle%amg)
            deallocate (handle%amg)
        end if
        if (allocated(handle%ams)) then
            call ams_destroy(handle%ams)
            deallocate (handle%ams)
        end if
        call release(handle%a)
        call release(handle%l)
        if (allocated(handle%inverse_diagonal)) deallocate (handle%inverse_diagonal)
        if (allocated(handle%r)) deallocate (handle%r, handle%z, handle%p, handle%q)
        handle%preconditioner = PCG_NONE
        handle%shift = 0
        handle%ready = .false.
    end subroutine pcg_destroy

    subroutine release(a)
        type(csr_t), intent(inout) :: a
        if (allocated(a%row_ptr)) deallocate (a%row_ptr)
        if (allocated(a%col)) deallocate (a%col)
        if (allocated(a%val)) deallocate (a%val)
        a%nrow = 0
        a%ncol = 0
    end subroutine release

    subroutine precondition(handle, r, z, info)
        type(pcg_handle_t), intent(inout) :: handle
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        type(fortsparse_status_t), intent(out) :: info
        integer :: i

        call status_set(info, FORTSPARSE_OK, '')
        select case (handle%preconditioner)
        case (PCG_JACOBI)
            !$omp parallel do if(size(r) >= 4096)
            do i = 1, size(r)
                z(i) = handle%inverse_diagonal(i)*r(i)
            end do
        case (PCG_IC0)
            call ic0_apply(handle%l, r, z)
        case (PCG_AMG)
            call amg_apply(handle%amg, r, z, info)
        case (PCG_AMS)
            call ams_apply(handle%ams, r, z, info)
        case default
            z = r
        end select
    end subroutine precondition

    subroutine jacobi_setup(handle, info)
        type(pcg_handle_t), intent(inout) :: handle
        type(fortsparse_status_t), intent(out) :: info
        integer :: i, k
        real(dp) :: d

        allocate (handle%inverse_diagonal(handle%a%nrow))
        do i = 1, handle%a%nrow
            d = 0
            do k = handle%a%row_ptr(i), handle%a%row_ptr(i + 1) - 1
                if (handle%a%col(k) == i) d = handle%a%val(k)
            end do
            if (.not. (d > 0)) then
                call status_set(info, FORTSPARSE_NOT_POSITIVE_DEFINITE, &
                    'PCG: Jacobi needs a positive diagonal')
                return
            end if
            handle%inverse_diagonal(i) = 1/d
        end do
        call status_set(info, FORTSPARSE_OK, '')
    end subroutine jacobi_setup

    subroutine ic0_setup(handle, info)
        !! IC(0) on the pattern of the lower triangle; a pivot that is not
        !! positive restarts with diag(A) scaled by 1 + shift, shift doubling
        !! from 1e-3, as semidefinite and indefinite-ish matrices need.
        type(pcg_handle_t), intent(inout) :: handle
        type(fortsparse_status_t), intent(out) :: info
        real(dp) :: shift
        logical :: ok
        integer :: attempt

        call lower_pattern(handle%a, handle%l)
        shift = 0
        do attempt = 1, 30
            call ic0_factor(handle%a, shift, handle%l, ok)
            if (ok) then
                handle%shift = shift
                call status_set(info, FORTSPARSE_OK, '')
                return
            end if
            if (shift == 0) then
                shift = 1.0e-3_dp
            else
                shift = 2*shift
            end if
        end do
        call status_set(info, FORTSPARSE_NOT_POSITIVE_DEFINITE, &
            'PCG: IC(0) has no positive pivots for any tried shift')
    end subroutine ic0_setup

    subroutine lower_pattern(a, l)
        !! The lower triangle of a with each row's diagonal entry last; a row
        !! without a stored diagonal gets one.
        type(csr_t), intent(in) :: a
        type(csr_t), intent(out) :: l
        integer :: i, k, m

        l%nrow = a%nrow
        l%ncol = a%ncol
        allocate (l%row_ptr(a%nrow + 1))
        l%row_ptr(1) = 1
        do i = 1, a%nrow
            m = 1
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                if (a%col(k) < i) m = m + 1
            end do
            l%row_ptr(i + 1) = l%row_ptr(i) + m
        end do
        allocate (l%col(l%row_ptr(a%nrow + 1) - 1), l%val(l%row_ptr(a%nrow + 1) - 1))
        do i = 1, a%nrow
            m = l%row_ptr(i)
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                if (a%col(k) < i) then
                    l%col(m) = a%col(k)
                    m = m + 1
                end if
            end do
            l%col(m) = i
        end do
    end subroutine lower_pattern

    subroutine ic0_factor(a, shift, l, ok)
        type(csr_t), intent(in) :: a
        real(dp), intent(in) :: shift
        type(csr_t), intent(inout) :: l
        logical, intent(out) :: ok
        integer :: i, k, j, p, q, pend, qend
        real(dp) :: s

        ok = .false.
        ! Values of A on L's pattern; the diagonal scaled by 1 + shift.
        l%val = 0
        do i = 1, a%nrow
            p = l%row_ptr(i)
            pend = l%row_ptr(i + 1) - 1
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                j = a%col(k)
                if (j < i) then
                    do while (l%col(p) < j)
                        p = p + 1
                    end do
                    l%val(p) = a%val(k)
                else if (j == i) then
                    l%val(pend) = a%val(k)*(1 + shift)
                end if
            end do
        end do
        do i = 1, a%nrow
            pend = l%row_ptr(i + 1) - 1
            do p = l%row_ptr(i), pend - 1
                k = l%col(p)
                ! L(i,k) = (A(i,k) - sum_{j<k} L(i,j) L(k,j)) / L(k,k)
                s = l%val(p)
                q = l%row_ptr(k)
                qend = l%row_ptr(k + 1) - 1
                do j = l%row_ptr(i), p - 1
                    do while (q < qend .and. l%col(q) < l%col(j))
                        q = q + 1
                    end do
                    if (q < qend .and. l%col(q) == l%col(j)) s = s - l%val(j)*l%val(q)
                end do
                l%val(p) = s/l%val(qend)
            end do
            s = l%val(pend)
            do p = l%row_ptr(i), pend - 1
                s = s - l%val(p)**2
            end do
            if (.not. (s > 0) .or. .not. ieee_is_finite(s)) return
            l%val(pend) = sqrt(s)
        end do
        ok = .true.
    end subroutine ic0_factor

    subroutine ic0_apply(l, r, z)
        !! z = (L L^T)^{-1} r.
        type(csr_t), intent(in) :: l
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        integer :: i, p, pend
        real(dp) :: s

        do i = 1, l%nrow
            pend = l%row_ptr(i + 1) - 1
            s = r(i)
            do p = l%row_ptr(i), pend - 1
                s = s - l%val(p)*z(l%col(p))
            end do
            z(i) = s/l%val(pend)
        end do
        do i = l%nrow, 1, -1
            pend = l%row_ptr(i + 1) - 1
            z(i) = z(i)/l%val(pend)
            do p = l%row_ptr(i), pend - 1
                z(l%col(p)) = z(l%col(p)) - l%val(p)*z(i)
            end do
        end do
    end subroutine ic0_apply

    logical function valid(a)
        type(csr_t), intent(in) :: a
        integer :: i, k

        valid = .false.
        if (a%nrow <= 0 .or. a%ncol /= a%nrow) return
        if (.not. allocated(a%row_ptr) .or. .not. allocated(a%col) .or. &
            .not. allocated(a%val)) return
        if (size(a%row_ptr) /= a%nrow + 1) return
        if (a%row_ptr(1) /= 1 .or. a%row_ptr(a%nrow + 1) - 1 /= size(a%col)) return
        if (size(a%val) /= size(a%col)) return
        do i = 1, a%nrow
            if (a%row_ptr(i + 1) < a%row_ptr(i)) return
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                if (a%col(k) < 1 .or. a%col(k) > a%ncol) return
                if (k > a%row_ptr(i)) then
                    if (a%col(k) <= a%col(k - 1)) return
                end if
            end do
        end do
        valid = .true.
    end function valid

    real(dp) function dot(x, y)
        real(dp), intent(in) :: x(:), y(:)
        integer :: i

        dot = 0
        !$omp parallel do reduction(+:dot) if(size(x) >= 4096)
        do i = 1, size(x)
            dot = dot + x(i)*y(i)
        end do
    end function dot

    real(dp) function norm(x)
        real(dp), intent(in) :: x(:)

        norm = sqrt(dot(x, x))
    end function norm

    subroutine axpy(alpha, x, y)
        real(dp), intent(in) :: alpha, x(:)
        real(dp), intent(inout) :: y(:)
        integer :: i

        !$omp parallel do if(size(x) >= 4096)
        do i = 1, size(x)
            y(i) = y(i) + alpha*x(i)
        end do
    end subroutine axpy
end module fortsparse_pcg
