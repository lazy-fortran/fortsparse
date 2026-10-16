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
    integer :: i, j, p, k, ordering, rows(n*n), cols(n*n), device
    character(16) :: argument
    real(dp) :: dense(n, n), vals(n*n), exact(n), b(n), x(n), y(n)
    real(dp) :: Bblock(n, 8), Xblock(n, 8), dx(n), seed(n), bbar(n)
    real(dp), allocatable :: abar(:)
    real(dp) :: zero(n), lhs, rhs
    logical :: reused

    call get_command_argument(1, argument)
    device = -1
    if (trim(argument) == '--cuda') device = 0

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
    solver%options%device = device
    other%backend_id = FORTSPARSE_BACKEND_NATIVE
    other%options%symmetry = FORTSPARSE_SYM_SPD
    other%options%device = device
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
    other%options%device = -1
    other%backend_id = FORTSPARSE_BACKEND_NATIVE_REFERENCE
    call sparse_factor(other, A, status)
    call require(status_ok(status), 'explicit reference selection')
    call sparse_solve(other, b, x, status)
    call close(x, exact, 'reference factory switch')
    other%backend_id = FORTSPARSE_BACKEND_NATIVE
    other%options%device = device
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
    call blocked_ordering_regression()
    print '(a)', 'native real SPD value, reuse, block, derivative and failure tests PASS'
contains
    subroutine blocked_ordering_regression()
        ! The bandwidth-16 graph contains a 17-clique, so every ordering has
        ! a column with more than 12 factor entries. With n >= 64 this uses
        ! blocked panels, including the original-CSC permutation loader.
        integer, parameter :: nb = 96, band = 16
        type(csc_t) :: matrix, tangent, scaled
        type(sparse_solver_t) :: ordered_solver
        type(sparse_factor_stats_t) :: ordered_stats
        integer :: ii, jj, nz, order, rr(nb*nb), cc(nb*nb), block
        real(dp) :: dd(nb,nb), vv(nb*nb), weight, solution(nb), rhs0(nb)
        real(dp) :: result(nb), rhs_block(nb,4), result_block(nb,4)
        real(dp) :: direction(nb), cotangent(nb), rhs_cotangent(nb), zeros(nb)
        real(dp), allocatable :: matrix_cotangent(:)
        logical :: analysis_reused

        dd = 0.0_dp
        do ii = 1, nb
            dd(ii,ii) = 1.0_dp
        end do
        do jj = 1, nb
            do ii = jj + 1, min(nb,jj + band)
                weight = 0.01_dp*real(1 + mod(ii + 3*jj,7),dp)
                dd(ii,jj) = -weight
                dd(jj,ii) = -weight
                dd(ii,ii) = dd(ii,ii) + weight
                dd(jj,jj) = dd(jj,jj) + weight
            end do
        end do
        nz = 0
        do jj = 1, nb
            do ii = 1, nb
                if (dd(ii,jj) == 0.0_dp) cycle
                nz = nz + 1
                rr(nz) = ii
                cc(nz) = jj
                vv(nz) = dd(ii,jj)
            end do
        end do
        call csc_from_triplet(nb, nb, rr(:nz), cc(:nz), vv(:nz), matrix, status)
        call require(status_ok(status), 'blocked ordering construct')
        solution = [(sin(real(ii,dp)/7.0_dp), ii=1,nb)]
        rhs0 = matmul(dd, solution)
        tangent = matrix
        tangent%val = 0.125_dp*matrix%val
        scaled = matrix
        zeros = 0.0_dp
        cotangent = [(cos(real(ii,dp)/9.0_dp), ii=1,nb)]
        allocate (matrix_cotangent(matrix%nnz))
        ordered_solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        ordered_solver%options%symmetry = FORTSPARSE_SYM_SPD
        ordered_solver%options%device = device
        do order = FORTSPARSE_ORDER_AMD, FORTSPARSE_ORDER_ND
            ordered_solver%options%ordering = order
            call sparse_factor(ordered_solver, matrix, status)
            if (status%code == FORTSPARSE_BACKEND_UNAVAILABLE) cycle
            call require(status_ok(status), 'blocked ordered factor')
            call sparse_factor_stats(ordered_solver, ordered_stats, status)
            call require(status_ok(status), 'blocked ordered stats')
            call require(ordered_stats%ordering_used == order, 'blocked ordering reported')
            call sparse_solve(ordered_solver, rhs0, result, status)
            call require(status_ok(status), 'blocked ordered solve')
            call close(result, solution, 'blocked ordered independent solution')
            do block = 1, 4
                rhs_block(:,block) = block*rhs0
            end do
            call sparse_solve_block(ordered_solver, rhs_block, result_block, status)
            call require(status_ok(status), 'blocked ordered block solve')
            do block = 1, 4
                call close(result_block(:,block), block*solution, 'blocked ordered block values')
            end do
            call sparse_solve_jvp(ordered_solver, tangent, result, zeros, direction, status)
            call require(status_ok(status), 'blocked ordered JVP')
            call close(direction, -0.125_dp*solution, 'blocked ordered analytic tangent')
            call sparse_solve_vjp(ordered_solver, matrix, result, cotangent, &
                rhs_cotangent, matrix_cotangent, status)
            call require(status_ok(status), 'blocked ordered VJP')
            call require(abs(dot_product(cotangent,direction) - &
                dot_product(matrix_cotangent,tangent%val)) < 2.0e-12_dp, &
                'blocked ordered tangent/adjoint pairing')
            do block = 1, 2
                scaled%val = real(block + 1,dp)*matrix%val
                call sparse_refactor(ordered_solver, scaled, status, analysis_reused)
                call require(status_ok(status) .and. analysis_reused, 'blocked ordered refactor')
                call sparse_solve(ordered_solver, real(block + 1,dp)*rhs0, result, status)
                call require(status_ok(status), 'blocked ordered refactor solve')
                call close(result, solution, 'blocked ordered refactor values')
            end do
        end do
        call sparse_destroy(ordered_solver)
    end subroutine blocked_ordering_regression
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
