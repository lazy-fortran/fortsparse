module amg_test_support
    use fortsparse, only: dp, csr_t, csr_matvec, amg_handle_t, amg_apply, &
                          sparse_ilutp_t, fortsparse_status_t
    implicit none
    private
    public :: poisson_grid, pcg, check
contains
    subroutine check(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (ok) return
        print *, message
        error stop 1
    end subroutine

    subroutine poisson_grid(side, dim, q1, anisotropy, a)
        !! Synthetic Dirichlet interior grid. P1-like uses the 5/7-point graph.
        !! Q1 uses tensor products of exact 1D stiffness and consistent mass.
        integer, intent(in) :: side, dim
        logical, intent(in) :: q1
        real(dp), intent(in) :: anisotropy
        type(csr_t), intent(out) :: a
        integer :: n, ix, iy, iz, dx, dy, dz, jx, jy, jz, row, nz, col
        real(dp) :: value, wx, wy, wz
        real(dp), parameter :: mass(-1:1) = [1/6.0_dp, 2/3.0_dp, 1/6.0_dp]
        real(dp), parameter :: stiff(-1:1) = [-1.0_dp, 2.0_dp, -1.0_dp]
        n = side**dim
        wx = 1
        wy = anisotropy
        wz = anisotropy
        a%nrow = n
        a%ncol = n
        allocate (a%row_ptr(n + 1), a%col(n*3**dim), a%val(n*3**dim))
        nz = 0
        row = 0
        do iz = 1, merge(side, 1, dim == 3)
            do iy = 1, side
                do ix = 1, side
                    row = row + 1
                    a%row_ptr(row) = nz + 1
                    do dz = merge(-1, 0, dim == 3), merge(1, 0, dim == 3)
                        jz = iz + dz
                        if (jz < 1 .or. jz > merge(side, 1, dim == 3)) cycle
                        do dy = -1, 1
                            jy = iy + dy
                            if (jy < 1 .or. jy > side) cycle
                            do dx = -1, 1
                                jx = ix + dx
                                if (jx < 1 .or. jx > side) cycle
                                if (q1) then
                                    value = wx*stiff(dx)*mass(dy) &
                                            + wy*mass(dx)*stiff(dy)
                                    if (dim == 3) value = value*mass(dz) &
                                                        + wz*mass(dx)*mass(dy)*stiff(dz)
                                else
                                    if (abs(dx) + abs(dy) + abs(dz) > 1) cycle
                                    value = 0
                                    if (dx /= 0) value = -wx
                                    if (dy /= 0) value = -wy
                                    if (dz /= 0) value = -wz
                                    if (dx == 0 .and. dy == 0 .and. dz == 0) then
                                        value = 2*(wx + wy)
                                        if (dim == 3) value = value + 2*wz
                                    end if
                                end if
                                if (abs(value) < 1.e-16_dp) cycle
                                col = jx + side*(jy - 1) + side**2*(jz - 1)
                                nz = nz + 1
                                a%col(nz) = col
                                a%val(nz) = value
                            end do
                        end do
                    end do
                end do
            end do
        end do
        a%row_ptr(n + 1) = nz + 1
        a%col = a%col(:nz)
        a%val = a%val(:nz)
    end subroutine

    subroutine pcg(a, b, x, iterations, residual, amg, ilu, tolerance, maxiter)
        type(csr_t), intent(in) :: a
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: x(:), residual
        integer, intent(out) :: iterations
        type(amg_handle_t), intent(inout), optional :: amg
        type(sparse_ilutp_t), intent(inout), optional :: ilu
        real(dp), intent(in), optional :: tolerance
        integer, intent(in), optional :: maxiter
        type(fortsparse_status_t) :: info
        real(dp), allocatable :: r(:), z(:), p(:), ap(:)
        real(dp) :: rho, previous, alpha, denom, normb, tol
        integer :: k, limit, i
        allocate (r(size(b)), z(size(b)), p(size(b)), ap(size(b)))
        tol = 1.e-12_dp
        if (present(tolerance)) tol = tolerance
        limit = 500
        if (present(maxiter)) limit = maxiter
        normb = sqrt(parallel_dot(b, b))
        x = 0
        r = b
        previous = 1
        p = 0
        iterations = 0
        residual = sqrt(parallel_dot(r, r))/max(normb, tiny(1.0_dp))
        if (normb == 0) return
        do k = 1, limit
            if (present(amg)) then
                call amg_apply(amg, r, z, info)
                call check(info%code == 0, 'AMG apply: '//trim(info%msg))
            else if (present(ilu)) then
                call ilu%apply(r, z, info)
                call check(info%code == 0, 'ILUTP apply: '//trim(info%msg))
            else
                z = r
            end if
            rho = parallel_dot(r, z)
            if (rho <= 0) exit
            if (k == 1) then
                !$omp parallel do if(size(b) >= 512)
                do i = 1, size(b)
                    p(i) = z(i)
                end do
            else
                !$omp parallel do if(size(b) >= 512)
                do i = 1, size(b)
                    p(i) = z(i) + (rho/previous)*p(i)
                end do
            end if
            call csr_matvec(a, p, ap)
            denom = parallel_dot(p, ap)
            if (denom <= 0) exit
            alpha = rho/denom
            !$omp parallel do if(size(b) >= 512)
            do i = 1, size(b)
                x(i) = x(i) + alpha*p(i)
                r(i) = r(i) - alpha*ap(i)
            end do
            iterations = k
            previous = rho
            if (sqrt(parallel_dot(r, r)) > tol*normb) cycle
            ! Independent true residual decides convergence, not CG recurrence.
            call csr_matvec(a, x, ap)
            !$omp parallel do if(size(b) >= 512)
            do i = 1, size(b)
                r(i) = b(i) - ap(i)
            end do
            if (sqrt(parallel_dot(r, r)) <= tol*normb) exit
            ! Reliable restart if recursively updated residual drifted.
            previous = 1
            p = 0
            call check(k < limit, 'PCG true residual did not converge')
        end do
        call csr_matvec(a, x, ap)
        !$omp parallel do if(size(b) >= 512)
        do i = 1, size(b)
            r(i) = b(i) - ap(i)
        end do
        residual = sqrt(parallel_dot(r, r))/normb
    end subroutine
    real(dp) function parallel_dot(x, y) result(value)
        real(dp), intent(in) :: x(:), y(:)
        integer :: i
        value = 0
        !$omp parallel do reduction(+:value) if(size(x) >= 512)
        do i = 1, size(x)
            value = value + x(i)*y(i)
        end do
    end function
end module
