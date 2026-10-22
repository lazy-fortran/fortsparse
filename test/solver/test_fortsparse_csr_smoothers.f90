program test_fortsparse_csr_smoothers
    use fortsparse, only: dp, csr_t, csr_matvec, csr_smoother_t, smoother_setup, &
             smoother_apply, smoother_destroy, SMOOTHER_L1_JACOBI, SMOOTHER_CHEBYSHEV, &
                          fortsparse_status_t
    use fortsparse_csr, only: csr_budget_t, csr_multiply, csr_transpose, csr_valid
    use fortsparse_amg, only: amg_tentative_prolongator
    use amg_test_support, only: poisson_grid, check
    implicit none
    type(csr_t) :: a, t, at, product, rap
    type(csr_budget_t) :: budget
    type(csr_smoother_t) :: smoother
    type(fortsparse_status_t) :: info
    real(dp), allocatable :: b(:), x(:), ax(:), basis(:, :), bc(:, :)
    real(dp), allocatable :: dense(:, :), tdense(:, :), pdense(:, :), coarse(:, :)
    integer :: method, i, k, n, j
    call poisson_grid(4, 2, .false., 1.0_dp, a)
    n = a%nrow
    allocate (b(n), x(n), ax(n), basis(n, 3))
    b = 1
    do method = SMOOTHER_L1_JACOBI, SMOOTHER_CHEBYSHEV
        call smoother_setup(a, method, smoother, info, degree=4)
        call check(info%code == 0, 'smoother setup')
        ! Exact upper spectral edge of this separable Dirichlet operator.
        call check(smoother%eigenvalue_estimate() <= smoother%eigenvalue_bound(), &
                   'power estimate below certified bound')
        if (method == SMOOTHER_CHEBYSHEV) then
            call check(smoother%eigenvalue_bound() >= 1 + cos(acos(-1.0_dp)/5), &
                       'known eigenvalue below bound')
        end if
        x = 0
        call smoother_apply(a, smoother, b, x, info)
        call check(info%code == 0, 'smoother apply')
        call csr_matvec(a, x, ax)
        call check(sum((b - ax)**2) < sum(b*b), 'standalone residual contraction')
        b = 0
        x = 0
        call smoother_apply(a, smoother, b, x, info)
        call check(info%code == 0 .and. maxval(abs(x)) == 0, 'zero RHS')
        b = 1
        call smoother_destroy(smoother)
    end do
    basis(:, 1) = 1
    do i = 1, n
        basis(i, 2) = real(mod(i - 1, 4), dp)
    end do
    basis(:, 3) = 2*basis(:, 1) ! locally rank-deficient caller basis
    call amg_tentative_prolongator(a, .08_dp, basis, t, bc, info)
    call check(info%code == 0, 'tentative QR')
    do j = 1, 3
        call csr_matvec(t, bc(:, j), ax)
        call check(maxval(abs(ax - basis(:, j))) < 1.e-12_dp, 'T Bc reproduces B')
    end do
    call csr_transpose(t, at, budget, info)
    call check(info%code == 0 .and. csr_valid(at), 'CSR transpose')
    call csr_multiply(a, t, product, budget, info)
    call check(info%code == 0 .and. csr_valid(product), 'CSR product')
    allocate (dense(n, n), tdense(n, t%ncol), pdense(n, t%ncol))
    dense = 0
    tdense = 0
    pdense = 0
    do i = 1, n
        do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
            dense(i, a%col(k)) = a%val(k)
        end do
        do k = t%row_ptr(i), t%row_ptr(i + 1) - 1
            tdense(i, t%col(k)) = t%val(k)
        end do
        do k = product%row_ptr(i), product%row_ptr(i + 1) - 1
            pdense(i, product%col(k)) = product%val(k)
        end do
    end do
    call check(maxval(abs(pdense - matmul(dense, tdense))) < 1.e-12_dp, &
               'sparse product versus dense oracle')
    call csr_multiply(at, product, rap, budget, info)
    call check(info%code == 0, 'Galerkin sparse product')
    allocate (coarse(t%ncol, t%ncol))
    coarse = 0
    do i = 1, rap%nrow
        do k = rap%row_ptr(i), rap%row_ptr(i + 1) - 1
            coarse(i, rap%col(k)) = rap%val(k)
        end do
    end do
    call check(maxval(abs(coarse - matmul(transpose(tdense), &
                        matmul(dense, tdense)))) < 1.e-12_dp, 'RAP versus dense oracle')
    do j = 1, at%nrow
        do k = at%row_ptr(j), at%row_ptr(j + 1) - 1
            call check(abs(at%val(k) - tdense(at%col(k), j)) < 1.e-14_dp, &
                       'transpose versus dense oracle')
        end do
    end do
end program
