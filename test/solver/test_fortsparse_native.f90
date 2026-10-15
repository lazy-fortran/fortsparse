program test_fortsparse_native
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortsparse
    implicit none
    type(csc_t) :: A, D, changed, invalid
    type(csc_z_t) :: Z
    type(sparse_solver_t) :: solver, other
    type(fortsparse_status_t) :: status
    type(sparse_factor_stats_t) :: stats
    integer, parameter :: n = 19
    integer :: i, j, p, k, ordering, rows(n*n), cols(n*n)
    real(dp) :: dense(n, n), vals(n*n), exact(n), b(n), x(n), y(n)
    real(dp) :: Bblock(n, 8), Xblock(n, 8), dx(n), seed(n), bbar(n)
    real(dp), allocatable :: abar(:)
    real(dp) :: zero(n), lhs, rhs
    logical :: reused

    ! A sparse SPD graph with fill: weighted graph Laplacian + I.
    ! This independent fixture is also invariant under the reversed numbering.
    dense = 0.0_dp
    do i = 1, n
        dense(i, i) = 1.0_dp
    end do
    do i = 1, n - 1
        call edge(i, i + 1, 0.5_dp)
        if (i + 4 <= n) call edge(i, i + 4, 0.3_dp)
    end do
    k = 0
    do j = 1, n
        do i = 1, n
            if (dense(i, j) == 0.0_dp) cycle
            k = k + 1
            rows(k) = i
            cols(k) = j
            vals(k) = dense(i, j)
        end do
    end do
    call csc_from_triplet(n, n, rows(:k), cols(:k), vals(:k), A, status)
    call require(status_ok(status), 'construct')
    solver%backend_id = FORTSPARSE_BACKEND_NATIVE
    solver%options%symmetry = FORTSPARSE_SYM_SPD
    other%backend_id = FORTSPARSE_BACKEND_NATIVE
    other%options%symmetry = FORTSPARSE_SYM_SPD
    exact = [(sin(real(i, dp)), i=1, n)]
    b = matmul(dense, exact)
    call sparse_factor(solver, A, status)
    call require(status_ok(status), 'factor sparse graph')
    call sparse_solve(solver, b, x, status)
    call require(status_ok(status), 'solve')
    call close(x, exact, 'known independent solution')
    call sparse_solve(solver, b(:n - 1), y, status)
    call require(status%code == FORTSPARSE_INVALID_MATRIX, 'RHS dimensions')
    y = b
    call sparse_solve(solver, y, status)
    call close(y, exact, 'in-place')
    do j = 1, 8
        Bblock(:, j) = j*b
    end do
    call sparse_solve_block(solver, Bblock, Xblock, status)
    call require(status_ok(status), 'block solve')
    do j = 1, 8
        call close(Xblock(:, j), j*exact, 'block known solution')
    end do
    call sparse_factor_stats(solver, stats, status)
    call require(status_ok(status), 'stats')
    call require(stats%factor_nnz >= (A%nnz + n)/2, 'factor fill')
    call require(stats%factor_mem_mb > 0, 'memory accounting')
    D = A
    D%val = 0.125_dp*A%val
    zero = 0.0_dp
    call sparse_solve_jvp(solver, D, x, zero, dx, status)
    call require(status_ok(status), 'implicit JVP')
    call close(dx, -0.125_dp*exact, 'analytic tangent')
    seed = [(cos(real(i, dp)), i=1, n)]
    allocate (abar(A%nnz))
    call sparse_solve_vjp(solver, A, x, seed, bbar, abar, status)
    call require(status_ok(status), 'implicit VJP')
    lhs = dot_product(seed, dx)
    rhs = dot_product(abar, D%val)
    call require(abs(lhs - rhs) < 2.0e-12_dp, 'tangent/adjoint pairing')
    changed = A
    changed%val = 2*A%val
    call sparse_refactor(solver, changed, status, reused)
    call require(status_ok(status) .and. reused, 'reuse analysis')
    call sparse_solve(solver, 2*b, x, status)
    call close(x, exact, 'refactor values')
    ! A second handle must be independent, including after first is destroyed.
    call sparse_factor(other, A, status)
    call sparse_destroy(solver)
    call sparse_solve(other, b, x, status)
    call close(x, exact, 'independent handles')
    ! Explicit ordering exercises permutation and changed-value reuse, or
    ! returns the optional dependency's documented unavailable status.
    do ordering = FORTSPARSE_ORDER_AMD, FORTSPARSE_ORDER_ND
        other%options%ordering = ordering
        call sparse_factor(other, A, status)
        if (status%code == FORTSPARSE_BACKEND_UNAVAILABLE) cycle
        call require(status_ok(status), 'explicit ordering')
        call sparse_factor_stats(other, stats, status)
        call require(stats%ordering_used == ordering, 'reported ordering')
        call sparse_solve(other, b, x, status)
        call close(x, exact, 'permuted solve')
        call sparse_refactor(other, changed, status, reused)
        call require(status_ok(status) .and. reused, 'permuted refactor')
        call sparse_solve(other, 2*b, x, status)
        call close(x, exact, 'permuted refactor solve')
    end do
    other%options%ordering = FORTSPARSE_ORDER_DEFAULT
    other%backend_id = FORTSPARSE_BACKEND_NATIVE_REFERENCE
    call sparse_factor(other, A, status)
    call require(status_ok(status), 'explicit reference selection')
    call sparse_solve(other, b, x, status)
    call close(x, exact, 'reference factory switch')
    other%backend_id = FORTSPARSE_BACKEND_NATIVE
    call sparse_factor(other, A, status)
    invalid = A
    invalid%val(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call sparse_factor(other, invalid, status)
    call require(status%code == FORTSPARSE_INVALID_MATRIX, 'NaN rejected')
    call sparse_solve(other, b, x, status)
    call require(status%code == FORTSPARSE_NOT_FACTORED, 'failed factor invalidates')
    invalid = A
    invalid%val(2) = invalid%val(2) + 1.0_dp
    call sparse_factor(other, invalid, status)
    call require(status%code == FORTSPARSE_INVALID_MATRIX, 'asymmetry rejected')
    invalid = A
    invalid%val = -A%val
    call sparse_factor(other, invalid, status)
    call require(status%code == FORTSPARSE_NOT_POSITIVE_DEFINITE, &
        'nonpositive pivot rejected')
    call sparse_factor(other, Z, status)
    call require(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
        'complex explicitly unsupported')
    ! Changed graph (diagonal-only) cannot reuse the old graph analysis.
    call csc_from_triplet(n, n, [(i, i=1, n)], [(i, i=1, n)], &
        [(2.0_dp, i=1, n)], changed, status)
    call sparse_refactor(other, changed, status, reused)
    call require(status_ok(status) .and. .not. reused, 'changed pattern')
    call sparse_solve(other, 2*exact, x, status)
    call close(x, exact, 'disconnected diagonal blocks')
    call sparse_destroy(other)
    print '(a)', 'native real SPD value, reuse, block, derivative and failure tests PASS'
contains
    subroutine edge(a, b, weight)
        integer, intent(in) :: a, b
        real(dp), intent(in) :: weight
        dense(a, a) = dense(a, a) + weight
        dense(b, b) = dense(b, b) + weight
        dense(a, b) = -weight
        dense(b, a) = -weight
    end subroutine edge
    subroutine require(ok, name)
        logical, intent(in) :: ok
        character(*), intent(in) :: name
        if (.not. ok) then
            print *, 'FAIL ', name, ': ', trim(status%msg)
            error stop 1
        end if
    end subroutine require
    subroutine close(actual, expected, name)
        real(dp), intent(in) :: actual(:), expected(:)
        character(*), intent(in) :: name
        call require(maxval(abs(actual - expected)) < 2.0e-12_dp, name)
    end subroutine close
end program test_fortsparse_native
