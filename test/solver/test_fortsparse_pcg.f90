program test_fortsparse_pcg
    ! PCG on Dirichlet Poisson grids (SPD) and on a Neumann graph Laplacian
    ! (semidefinite, constants in the kernel) with a consistent right-hand
    ! side. Oracles: the true residual ||b - A x|| recomputed here, the error
    ! against the manufactured solution, and iteration counts between the
    ! preconditioners.
    use fortsparse, only: dp, csr_t, csr_matvec, fortsparse_status_t, &
        pcg_params_t, pcg_handle_t, pcg_setup, pcg_setup_csr, pcg_solve, &
        pcg_destroy, PCG_NONE, PCG_JACOBI, PCG_IC0, PCG_AMG, FORTSPARSE_OK, &
        FORTSPARSE_NOT_CONVERGED, FORTSPARSE_INVALID_MATRIX, FORTSPARSE_NOT_FACTORED
    use amg_test_support, only: poisson_grid, check
    implicit none
    type(csr_t) :: a
    type(pcg_params_t) :: params
    type(pcg_handle_t) :: handle
    type(fortsparse_status_t) :: info
    real(dp), allocatable :: b(:), x(:), exact(:), r(:)
    integer :: dim, kind, iterations, counts(0:3), i, n
    integer, parameter :: kinds(4) = [PCG_NONE, PCG_JACOBI, PCG_IC0, PCG_AMG]
    real(dp) :: relative

    params%rtol = 1.0e-10_dp
    params%amg%coarse_size = 16
    do dim = 2, 3
        call poisson_grid(merge(32, 10, dim == 2), dim, .false., 1.0_dp, a)
        call manufactured(a)
        do kind = 1, 4
            params%preconditioner = kinds(kind)
            call pcg_setup(a, params, handle, info)
            call check(info%code == FORTSPARSE_OK, 'setup: '//trim(info%msg))
            x = 0
            call pcg_solve(handle, params, b, x, info, iterations, relative)
            call check(info%code == FORTSPARSE_OK, 'SPD solve: '//trim(info%msg))
            call check(true_residual() <= 1.0e-10_dp*(1 + 1.0e-6_dp), &
                'the true residual meets the tolerance')
            call check(relative <= 1.0e-10_dp, 'the reported residual meets the tolerance')
            call check(maxval(abs(x - exact)) < 1.0e-6_dp, 'SPD solution error')
            counts(kinds(kind)) = iterations
            call pcg_destroy(handle)
        end do
        call check(counts(PCG_IC0) < counts(PCG_JACOBI), 'IC(0) needs fewer iterations than Jacobi')
        call check(counts(PCG_AMG) <= 30, 'AMG iteration bound')
        print '(a,i0,a,4(1x,i0))', 'dim ', dim, ' iterations none/Jacobi/IC0/AMG:', counts
    end do

    ! Semidefinite: the Neumann graph Laplacian; b = A exact lies in range(A).
    call neumann_laplacian(24, a)
    call manufactured(a)
    do kind = 1, 3
        params%preconditioner = kinds(kind)
        call pcg_setup(a, params, handle, info)
        call check(info%code == FORTSPARSE_OK, 'semidefinite setup: '//trim(info%msg))
        x = 0
        call pcg_solve(handle, params, b, x, info, iterations, relative)
        call check(info%code == FORTSPARSE_OK, 'consistent semidefinite solve: '//trim(info%msg))
        call check(true_residual() <= 1.0e-10_dp*(1 + 1.0e-6_dp), &
            'semidefinite true residual')
        ! x - exact is constant: the kernel component is free.
        call check(maxval(x - exact) - minval(x - exact) < 1.0e-6_dp, &
            'semidefinite solution up to the kernel')
        call pcg_destroy(handle)
    end do
    ! An inconsistent right-hand side has no solution and must not pass.
    params%preconditioner = PCG_JACOBI
    params%max_iterations = 500
    call pcg_setup(a, params, handle, info)
    b = 1
    x = 0
    call pcg_solve(handle, params, b, x, info)
    call check(info%code /= FORTSPARSE_OK, 'an inconsistent system is not reported solved')
    call pcg_destroy(handle)

    ! The iteration limit.
    call poisson_grid(32, 2, .false., 1.0_dp, a)
    call manufactured(a)
    params%preconditioner = PCG_NONE
    params%max_iterations = 3
    call pcg_setup(a, params, handle, info)
    x = 0
    call pcg_solve(handle, params, b, x, info, iterations)
    call check(info%code == FORTSPARSE_NOT_CONVERGED .and. iterations == 3, &
        'the iteration limit reports not converged')
    call pcg_destroy(handle)
    params%max_iterations = 1000

    ! Raw arrays, solving before setup and an unsorted row.
    params%preconditioner = PCG_IC0
    call pcg_solve(handle, params, b, x, info)
    call check(info%code == FORTSPARSE_NOT_FACTORED, 'solve before setup')
    call pcg_setup_csr(a%nrow, a%row_ptr, a%col, a%val, params, handle, info)
    x = 0
    call pcg_solve(handle, params, b, x, info)
    call check(info%code == FORTSPARSE_OK .and. maxval(abs(x - exact)) < 1.0e-6_dp, &
        'raw-array setup')
    call pcg_destroy(handle)
    i = a%row_ptr(2)
    a%col([i, i + 1]) = a%col([i + 1, i])
    call pcg_setup(a, params, handle, info)
    call check(info%code == FORTSPARSE_INVALID_MATRIX, 'unsorted columns are refused')
    ! A handle local to a routine is released when the routine returns.
    call poisson_grid(16, 2, .false., 1.0_dp, a)
    call manufactured(a)
    do i = 1, 3
        do kind = 1, 4
            call scoped_solve(kinds(kind))
        end do
    end do
    print '(a)', 'PASS'

contains

    subroutine scoped_solve(preconditioner)
        integer, intent(in) :: preconditioner
        type(pcg_handle_t) :: local
        type(pcg_params_t) :: p
        p%preconditioner = preconditioner
        p%amg%coarse_size = 16
        call pcg_setup(a, p, local, info)
        call check(info%code == FORTSPARSE_OK, 'scoped setup')
        x = 0
        call pcg_solve(local, p, b, x, info)
        call check(info%code == FORTSPARSE_OK, 'scoped solve')
        if (preconditioner /= PCG_JACOBI) call pcg_destroy(local)
    end subroutine scoped_solve

    subroutine manufactured(m)
        type(csr_t), intent(in) :: m
        n = m%nrow
        if (allocated(b)) deallocate (b, x, exact, r)
        allocate (b(n), x(n), exact(n), r(n))
        do i = 1, n
            exact(i) = sin(0.271_dp*i) + cos(0.173_dp*i)
        end do
        call csr_matvec(m, exact, b)
    end subroutine manufactured

    real(dp) function true_residual()
        call csr_matvec(a, x, r)
        true_residual = norm2(b - r)/norm2(b)
    end function true_residual

    subroutine neumann_laplacian(side, m)
        ! 5-point graph Laplacian on a side x side grid: row sums are zero.
        integer, intent(in) :: side
        type(csr_t), intent(out) :: m
        integer :: ix, iy, row, nz, k, jx, jy
        integer, parameter :: ox(5) = [0, -1, 0, 1, 0], oy(5) = [-1, 0, 0, 0, 1]
        m%nrow = side*side
        m%ncol = m%nrow
        allocate (m%row_ptr(m%nrow + 1), m%col(5*m%nrow), m%val(5*m%nrow))
        nz = 0
        row = 0
        do iy = 1, side
            do ix = 1, side
                row = row + 1
                m%row_ptr(row) = nz + 1
                do k = 1, 5
                    jx = ix + ox(k)
                    jy = iy + oy(k)
                    if (jx < 1 .or. jx > side .or. jy < 1 .or. jy > side) cycle
                    nz = nz + 1
                    m%col(nz) = jx + (jy - 1)*side
                    m%val(nz) = -1
                end do
            end do
        end do
        m%row_ptr(m%nrow + 1) = nz + 1
        m%col = m%col(:nz)
        m%val = m%val(:nz)
        ! The diagonal makes every row sum zero.
        do row = 1, m%nrow
            do k = m%row_ptr(row), m%row_ptr(row + 1) - 1
                if (m%col(k) == row) m%val(k) = -(sum(m%val(m%row_ptr(row):m%row_ptr(row + 1) - 1)) + 1)
            end do
        end do
    end subroutine neumann_laplacian
end program test_fortsparse_pcg
