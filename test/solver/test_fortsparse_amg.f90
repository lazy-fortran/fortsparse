program test_fortsparse_amg
    use fortsparse, only: dp, csr_t, csr_matvec, amg_params_t, amg_handle_t, &
                amg_setup, amg_apply, amg_destroy, amg_setup_csr, fortsparse_status_t, &
               SMOOTHER_L1_JACOBI, FORTSPARSE_INVALID_MATRIX, FORTSPARSE_NOT_FACTORED, &
                          FORTSPARSE_NOT_POSITIVE_DEFINITE, FORTSPARSE_RESOURCE_LIMIT
    use fortsparse_csc, only: csc_t
    use fortsparse_native_ldlt, only: native_ldlt_t
    use amg_test_support, only: poisson_grid, pcg, check
    implicit none
    type(csr_t) :: a
    type(csc_t) :: csc
    type(amg_params_t) :: params
    type(amg_handle_t) :: handle
    type(native_ldlt_t) :: direct
    type(fortsparse_status_t) :: info
    real(dp), allocatable :: b(:), x(:), exact(:), reference(:), y(:), z(:)
    real(dp) :: residual, dot1, dot2, eps
    integer :: dim, family, refinement, side, n, i, iterations, cycle, exponent
    params%coarse_size = 12
    do family = 1, 2
        do dim = 2, 3
            do refinement = 1, 2
                side = 4*refinement
                if (dim == 3) side = 2 + 2*refinement
                call poisson_grid(side, dim, family == 2, 1.0_dp, a)
                n = a%nrow
                allocate (b(n), x(n), exact(n), reference(n), y(n), z(n))
                do i = 1, n
                    exact(i) = sin(real(i, dp)*.271_dp) + cos(real(i, dp)*.173_dp)
                    y(i) = cos(real(i, dp))
                end do
                call csr_matvec(a, exact, b)
                do cycle = 1, 2
                    params%cycle = cycle
                    call amg_setup(a, params, handle, info)
                    call check(info%code == 0, 'setup: '//trim(info%msg))
                    call check(handle%levels() > 1, 'must exercise coarsening')
                    call pcg(a, b, x, iterations, residual, amg=handle)
                    call check(iterations <= 40, 'isotropic iteration bound')
                    call check(residual < 1.e-11_dp, 'true residual')
                    call check(maxval(abs(x - exact)) < 1.e-10_dp, 'manufactured error')
                    call amg_apply(handle, b, z, info)
                    call check(info%code == 0, 'apply')
                    dot1 = dot_product(y, z)
                    call amg_apply(handle, y, z, info)
                    call check(info%code == 0, 'apply second RHS')
                    dot2 = dot_product(b, z)
                    call check(abs(dot1 - dot2) < 1.e-11_dp*max(1.0_dp, abs(dot1)), &
                               'preconditioner symmetry')
                    call check(dot_product(y, z) > 0, 'preconditioner positivity')
                end do
                csc%nrow = n
                csc%ncol = n
                csc%nnz = size(a%val)
                csc%col_ptr = a%row_ptr
                csc%row_idx = a%col
                csc%val = a%val
                call direct%factor(csc, info)
                call check(info%code == 0, 'LDLT reference')
                call direct%solve(b, reference, info)
                call check(info%code == 0, 'LDLT solve')
               call check(maxval(abs(x - reference)) < 1.e-10_dp, 'AMG-PCG versus LDLT')
                call direct%free()
                deallocate (b, x, exact, reference, y, z)
            end do
        end do
    end do
    ! Anisotropic solution agreement is independent of the large iteration study.
    params%cycle = 1
    do dim = 2, 3
        do exponent = 2, 6
            eps = 10.0_dp**(-exponent)
            call poisson_grid(4, dim, .false., eps, a)
            n = a%nrow
            allocate(b(n), x(n), exact(n), reference(n))
            do i = 1, n
                exact(i) = sin(real(i, dp)*.271_dp)
            end do
            call csr_matvec(a, exact, b)
            call amg_setup(a, params, handle, info)
            call check(info%code == 0, 'anisotropic setup')
            call pcg(a, b, x, iterations, residual, amg=handle)
            csc%nrow = n
            csc%ncol = n
            csc%nnz = size(a%val)
            csc%col_ptr = a%row_ptr
            csc%row_idx = a%col
            csc%val = a%val
            call direct%factor(csc, info)
            call check(info%code == 0, 'anisotropic LDLT')
            call direct%solve(b, reference, info)
            call check(info%code == 0, 'anisotropic reference solve')
            call check(maxval(abs(x-reference)) < 1.e-10_dp, &
                       'anisotropic AMG-PCG versus LDLT')
            call direct%free()
            deallocate(b, x, exact, reference)
        end do
    end do
    ! Multi-vector caller basis and l1-Jacobi through the raw-array bridge.
    call poisson_grid(8, 2, .false., 1.0_dp, a)
    n = a%nrow
    allocate (params%near_nullspace(n, 2), b(n), x(n), exact(n))
    params%near_nullspace(:, 1) = 1
    do i = 1, n
        params%near_nullspace(i, 2) = real(mod(i - 1, 8), dp)
        exact(i) = sin(real(i, dp))
    end do
    params%smoother = SMOOTHER_L1_JACOBI
    call amg_setup_csr(n, a%row_ptr, a%col, a%val, params, handle, info)
    call check(info%code == 0, 'raw CSR and two-vector basis')
    call csr_matvec(a, exact, b)
    call pcg(a, b, x, iterations, residual, amg=handle)
    call check(maxval(abs(x - exact)) < 1.e-10_dp, 'two-vector solve')
    call amg_destroy(handle)
    call amg_destroy(handle)
    call amg_apply(handle, b, x, info)
    call check(info%code == FORTSPARSE_NOT_FACTORED, 'destroy invalidates handle')
    params%memory_limit_bytes = 1
    call amg_setup(a, params, handle, info)
    call check(info%code == FORTSPARSE_RESOURCE_LIMIT, 'strict cap unsupported')
    params%memory_limit_bytes = 0
    a%val(2) = 1.1_dp*a%val(2)
    call amg_setup(a, params, handle, info)
    call check(info%code == FORTSPARSE_INVALID_MATRIX, 'reject asymmetric CSR')
    a%nrow = 2
    a%ncol = 2
    a%row_ptr = [1, 3, 5]
    a%col = [1, 2, 1, 2]
    a%val = [1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp]
    deallocate (params%near_nullspace)
    call amg_setup(a, params, handle, info)
    call check(info%code == FORTSPARSE_NOT_POSITIVE_DEFINITE, 'indefinite coarse')
end program
