program test_fortsparse_ldlt_blocked
    use fortsparse
    implicit none
    integer, parameter :: n = 96
    type(csc_t) :: matrix, tangent_matrix
    type(sparse_solver_t) :: solver
    type(fortsparse_status_t) :: status
    type(sparse_factor_stats_t) :: stats
    real(dp) :: l(n, n), d(n, n), a(n, n), x(n), b(n), answer(n)
    integer :: i, j, k, rows(n*n), columns(n*n), negative, state
    real(dp) :: values(n*n)
    real(dp) :: tangent(n), bdot(n), bbar(n), rhs_block(n, 2), solution_block(n, 2)
    real(dp), allocatable :: abar(:)
    integer :: trial
    integer(i8) :: nn, nz, np
    logical :: reused
    do trial = 1, 4
        l = 0.0_dp
        d = 0.0_dp
        do i = 1, n
            l(i, i) = 1.0_dp
            d(i, i) = real((-1)**i, dp)*(2.0_dp + real(i, dp)/n)
        end do
        negative = n/2
        if (trial == 2) then
            ! Adjacent 2x2 pivots, with exactly one sign of each kind per pair.
            do i = 1, n - 1, 2
                d(i, i) = 0.0_dp
                d(i + 1, i + 1) = 0.0_dp
                d(i + 1, i) = 2.0_dp
                d(i, i + 1) = 2.0_dp
            end do
        end if
        do i = 1, n - 4
            l(i + 2, i) = 0.125_dp
            l(i + 4, i) = -0.0625_dp
        end do
        a = matmul(matmul(l, d), transpose(l))
        if (trial == 3) then
            ! First pivot must interchange with a distant candidate. The
            ! private dense trial must be discarded, then unrestricted BK used.
            ! Use a pure permutation of a block-diagonal signed spectrum.
            a = 0.0_dp
            do i = 2, n - 1
                a(i, i) = real((-1)**i, dp)*3.0_dp
            end do
            a(1, n) = 8.0_dp
            a(n, 1) = 8.0_dp
        end if
        k = 0
        do j = 1, n
            do i = 1, n
                if (trial /= 4 .and. abs(a(i, j)) <= 0.0_dp) cycle
                k = k + 1
                rows(k) = i; columns(k) = j; values(k) = a(i, j)
            end do
        end do
        call csc_from_triplet(n, n, rows(:k), columns(:k), values(:k), matrix, status)
        call check()
        x = [(sin(real(i, dp)), i=1, n)]
        b = matmul(a, x)
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        solver%options%ldlt_blocked = .true.
        solver%refine = .false.
        call sparse_factor(solver, matrix, status)
        call check()
        call sparse_factor_stats(solver, stats, status)
        call check()
        call expect(stats%ldlt_blocked_panels > 0, 'no blocked panel exercised')
        call sparse_inertia(solver, nn, nz, np, state, status)
        call check()
        call expect(nn == negative .and. nz == 0 .and. np == n - negative, 'congruence inertia')
        call expect(state == FORTSPARSE_INERTIA_RELIABLE, 'well separated inertia uncertain')
        call sparse_solve(solver, b, answer, status)
        call check()
        call expect(maxval(abs(answer - x)) < 1.0e-12_dp, 'manufactured solve')
        call expect(maxval(abs(matmul(a, answer) - b)) < 1.0e-12_dp, 'original matrix residual')
        rhs_block(:, 1) = b
        rhs_block(:, 2) = -2.0_dp*b
        call sparse_solve_block(solver, rhs_block, solution_block, status)
        call check()
        call expect(maxval(abs(solution_block(:, 2) + 2.0_dp*x)) < 1e-12_dp, 'block solve')
        tangent = cos(x)
        tangent_matrix = matrix
        tangent_matrix%val = 0.125_dp*matrix%val
        bdot = matmul(a, tangent) + 0.125_dp*b
        call sparse_solve_jvp(solver, tangent_matrix, x, bdot, answer, status)
        call check()
        call expect(maxval(abs(answer - tangent)) < 1e-12_dp, 'implicit tangent')
        if (allocated(abar)) deallocate(abar)
        allocate(abar(matrix%nnz))
        bdot = matmul(a, tangent)
        call sparse_solve_vjp(solver, matrix, x, bdot, bbar, abar, status, use_primal_factors=.true.)
        call check()
        call expect(maxval(abs(bbar - tangent)) < 1e-12_dp, 'implicit adjoint')
        do j = 1, n
            do i = matrix%col_ptr(j), matrix%col_ptr(j + 1) - 1
                call expect(abs(abar(i) + tangent(matrix%row_idx(i))*x(j)) < 1e-12_dp, 'sparse cotangent')
            end do
        end do
        matrix%val = -matrix%val
        call sparse_refactor(solver, matrix, status, reused)
        call check()
        call expect(reused, 'same-pattern blocked refactor')
        call sparse_solve(solver, -b, answer, status)
        call check()
        call expect(maxval(abs(answer - x)) < 1.0e-12_dp, 'refactored solve')
        call sparse_destroy(solver)
    end do
    call wide_recovery()
    print '(a)', 'blocked BK known congruences, 2x2 pivots and global recovery pass'
contains
    subroutine wide_recovery()
        ! A far 2x2 pivot has a wide, populated trailing support. This exercises
        ! the row-merge recovery with multiple workers, against independently
        ! formed L D L^T, its exact congruence inertia and a known solution.
        integer, parameter :: wide_n = 160
        real(dp) :: ll(wide_n, wide_n), dd(wide_n, wide_n)
        real(dp) :: aa(wide_n, wide_n), known(wide_n), rhs(wide_n), got(wide_n)
        real(dp) :: entries(wide_n*wide_n)
        integer :: rr(wide_n*wide_n), cc(wide_n*wide_n), row, col, entry
        type(csc_t) :: wide
        type(sparse_solver_t) :: factor
        ll = 0.0_dp
        dd = 0.0_dp
        do row = 1, wide_n
            ll(row, row) = 1.0_dp
            dd(row, row) = real((-1)**row, dp)*3.0_dp
        end do
        ll(2:, 1) = 0.0078125_dp
        dd(1, 1) = 0.25_dp
        dd(wide_n, wide_n) = -0.25_dp
        dd(1, wide_n) = 8.0_dp
        dd(wide_n, 1) = 8.0_dp
        aa = matmul(matmul(ll, dd), transpose(ll))
        entry = 0
        do col = 1, wide_n
            do row = 1, wide_n
                entry = entry + 1
                rr(entry) = row
                cc(entry) = col
                entries(entry) = aa(row, col)
            end do
        end do
        call csc_from_triplet(wide_n, wide_n, rr, cc, entries, wide, status)
        call check()
        known = [(sin(real(row, dp)), row=1, wide_n)]
        rhs = matmul(aa, known)
        factor%backend_id = FORTSPARSE_BACKEND_NATIVE
        factor%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        factor%options%ldlt_blocked = .true.
        factor%refine = .false.
        call sparse_factor(factor, wide, status)
        call check()
        call sparse_inertia(factor, nn, nz, np, state, status)
        call check()
        call expect(nn == wide_n/2 .and. nz == 0 .and. np == wide_n/2, &
            'wide recovery congruence inertia')
        call expect(state == FORTSPARSE_INERTIA_RELIABLE, 'wide recovery reliability')
        call sparse_solve(factor, rhs, got, status)
        call check()
        call expect(maxval(abs(got - known)) < 1e-12_dp, 'wide recovery solution')
        call expect(maxval(abs(matmul(aa, got) - rhs)) < 1e-12_dp, &
            'wide recovery original matrix residual')
        call sparse_destroy(factor)
    end subroutine wide_recovery

    subroutine check()
        if (status%code /= FORTSPARSE_OK) then
            print *, status%code, trim(status%msg)
            error stop 'blocked LDLT status'
        end if
    end subroutine check
    subroutine expect(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print *, message
            error stop 'blocked LDLT oracle'
        end if
    end subroutine expect
end program test_fortsparse_ldlt_blocked
