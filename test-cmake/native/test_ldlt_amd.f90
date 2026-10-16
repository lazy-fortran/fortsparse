program test_ldlt_amd
    use fortsparse
    use, intrinsic :: iso_c_binding, only: c_int64_t
    use fortsparse_native_amd, only: bounded_amd_available
    use fortsparse_native_ldlt, only: native_ldlt_t
    implicit none
    integer, parameter :: n = 7
    type(sparse_solver_t) :: solver
    type(native_ldlt_t) :: direct, natural_layout
    type(csc_t) :: a, adot, diagonal_a
    type(sparse_factor_stats_t) :: stats
    type(fortsparse_status_t) :: status
    real(dp) :: dense(n, n), exact(n), rhs(n), answer(n), block(n, 2), result(n, 2)
    real(dp) :: tangent(n), bdot(n), seed(n), bbar(n)
    real(dp), allocatable :: abar(:)
    integer :: i, j, k, rr(n*n), cc(n*n), q(n), badq(n), copied(n)
    real(dp) :: vv(n*n)
    integer(i8) :: cap
    integer :: calls_before
    logical :: reused, matched
    interface
        integer function calls() bind(C, name='fsparse_spd_probe_amd_calls')
        end function
        subroutine probe_begin(nth) bind(C, name='fsparse_spd_probe_begin')
            import c_int64_t
            integer(c_int64_t), value :: nth
        end subroutine
        subroutine probe_end(count, current, maximum) bind(C, name='fsparse_spd_probe_end')
            import c_int64_t
            integer(c_int64_t), intent(out) :: count, current, maximum
        end subroutine
    end interface
    dense = 0
    dense(1, 4) = 2; dense(4, &
 1) = 2 ! A zero-diagonal pivot pair in original coordinates.
    do i = 2, n
        if (i /= 4) dense(i, i) = (-1._dp)**i*(3 + i)
    end do
    do i = 2, n
        if (i /= 4) then
            dense(i, 4) = 0.1_dp*i; dense(4, i) = dense(i, 4)
        end if
    end do
    exact = [(sin(real(i, dp)), i=1, n)]
    call matrix()
    q = [7, 6, 5, 4, 3, 2, 1]
    diagonal_a%nrow = n; diagonal_a%ncol = n; diagonal_a%nnz = n
    diagonal_a%col_ptr = [(i, i=1, n + 1)]
    diagonal_a%row_idx = [(i, i=1, n)]
    diagonal_a%val = [(real((-1)**i*(i + 2), dp), i=1, n)]
    call natural_layout%factor(diagonal_a, status)
    call checked()
    call direct%factor(diagonal_a, status, column_permutation=q)
    call checked()
    call expect(direct%retained_bytes() == natural_layout%retained_bytes() + 4_i8*n, &
                'initial order payload not counted exactly')
    call direct%free()
    call natural_layout%free()
    call direct%factor(a, status, column_permutation=q)
    call checked()
    call direct%copy_column_order(n, a%col_ptr, a%row_idx, copied, matched)
    call expect(matched .and. all(copied == q), 'initial order copy')
    rhs = matmul(dense, exact)
    call direct%solve(rhs, answer, status)
    call checked()
    call expect(maxval(abs(answer - exact)) < 1e-12_dp, 'direct ordered LDLT solution')
    badq = q; badq(1) = badq(2)
    call direct%factor(a, status, column_permutation=badq)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX, 'duplicate initial order accepted')
    call direct%factor(a, status, column_permutation=q)
    call checked()
    call direct%refactor(a, status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX, 'changed initial ordering reused')
    call direct%free()

    solver%backend_id = FORTSPARSE_BACKEND_NATIVE
    solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
    solver%options%ordering = FORTSPARSE_ORDER_AMD
    solver%options%memory_policy = FORTSPARSE_MEMORY_LOW
    solver%refine = .false.
    call sparse_factor(solver, a, status)
    if (bounded_amd_available() /= 1) then
        call expect(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, 'AMD ABI unavailable status')
        stop
    end if
    call checked()
    call operations()
    call sparse_factor_stats(solver, stats, status)
    call checked()
    call expect(stats%native_factor_family == 'LDLT_REAL', 'ordered family changed')
    call expect(stats%ordering_used == FORTSPARSE_ORDER_AMD, 'actual order missing')
    call expect(stats%strict_memory_supported .and. stats%host_peak_bytes > 0, &
                'strict order peak')
    cap = stats%host_peak_bytes
    call sparse_destroy(solver)
    solver%options%host_memory_limit_bytes = cap - 1_i8
    call sparse_factor(solver, a, status)
    call expect(status%code == FORTSPARSE_RESOURCE_LIMIT, 'one byte below exact fresh peak')
    call sparse_destroy(solver)
    solver%options%host_memory_limit_bytes = cap
    call sparse_factor(solver, a, status)
    call checked()
    call operations()
    call allocation_checks()
    ! Reserve the original-matrix checker before exercising ordered refinement.
    solver%refine = .true.
    ! Changed pivot values keep coordinates, while numerical pivoting is fresh.
    dense(1, 1) = 8; dense(4, 4) = -7
    call matrix()
    ! Adding diagonal entries changes the pattern: reuse must be false.
    solver%options%host_memory_limit_bytes = 1000000_i8
    call sparse_refactor(solver, a, status, reused)
    call checked()
    call expect(.not. reused, 'changed original pattern reused')
    call operations()
    dense(1, 1) = -2; dense(4, 4) = 9
    call matrix()
    calls_before = calls()
    call sparse_refactor(solver, a, status, reused)
    call checked()
    call expect(calls() == calls_before, 'same pattern/order reran AMD')
    call expect(reused, 'same pattern/order not reused')
    call operations()
    calls_before = calls()
    a%val(a%col_ptr(1) + 1) = a%val(a%col_ptr(1) + 1) + 0.5_dp
    call sparse_factor(solver, a, status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX, 'asymmetric declared LDLT accepted')
    call expect(calls() == calls_before, 'invalid symmetry reached AMD')
    call matrix()
    call sparse_factor(solver, a, status)
    call checked()
    a%val = 0
    call sparse_refactor(solver, a, status, reused)
    call expect(status%code == FORTSPARSE_SINGULAR .and. .not. reused, &
                'singular ordered refactor')
    call sparse_solve(solver, rhs, answer, status)
    call expect(status%code == FORTSPARSE_NOT_FACTORED, &
                'failed ordered factors remained usable')
    call matrix()
    call sparse_refactor(solver, a, status, reused)
    call checked()
    call operations()
    solver%options%ordering = FORTSPARSE_ORDER_DEFAULT
    call sparse_refactor(solver, a, status, reused)
    call checked()
    call expect(.not. reused, 'changed ordering reused')
    call operations()
    call sparse_destroy(solver)
    solver%options%ordering = FORTSPARSE_ORDER_AMD
    solver%options%host_memory_limit_bytes = 0_i8
    solver%options%memory_policy = FORTSPARSE_MEMORY_AUTO
    call sparse_factor(solver, a, status)
    call checked()
    call operations()
    call sparse_factor_stats(solver, stats, status)
    call checked()
    call expect(stats%host_peak_bytes < 0 .and. .not. stats%strict_memory_supported, &
                'opaque uncapped AMD claimed a complete peak')
    call sparse_destroy(solver)
    solver%backend_id = FORTSPARSE_BACKEND_NATIVE_REFERENCE
    solver%options%ordering = FORTSPARSE_ORDER_AMD
    call sparse_factor(solver, a, status)
    call expect(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, 'reference accepted AMD')
    call sparse_destroy(solver)
    print '(a)', 'Ordered LDLT original-coordinate solves, derivatives, reuse and caps pass'
contains
    subroutine allocation_checks()
        integer(c_int64_t) :: allocations, live, high, nth, ignored
        call sparse_free(solver)
        solver%options%host_memory_limit_bytes = 1000000_i8
        call probe_begin(0_c_int64_t)
        call sparse_factor(solver, a, status)
        call checked()
        call sparse_factor_stats(solver, stats, status)
        call checked()
        call probe_end(allocations, live, high)
        call expect(live == stats%host_current_bytes, 'strict LDLT retained payload differs from allocation oracle')
        call expect(high == stats%host_peak_bytes, 'strict LDLT peak differs from allocation oracle')
        print '(a,i0,a,i0,a,i0)', 'LDLT allocation oracle: ', allocations, &
            ' allocations; retained=', live, '; peak=', high
        do nth = 1, allocations
            call sparse_free(solver)
            call probe_begin(nth)
            call sparse_factor(solver, a, status)
            call sparse_free(solver)
            call probe_end(ignored, live, high)
            call expect(live == 0_i8, 'failed LDLT allocation leaked owned payload')
            call expect(status%code == FORTSPARSE_OUT_OF_MEMORY, 'injected LDLT allocation failure status')
            call expect(.not. solver%factored, 'failed allocation retained a public factor')
        end do
        call sparse_free(solver)
        call sparse_factor(solver, a, status)
        call checked()
        call operations()
    end subroutine
    subroutine matrix()
        k = 0
        do j = 1, n
            do i = 1, n
                if (dense(i, j) == 0) cycle
                k = k + 1; rr(k) = i; cc(k) = j; vv(k) = dense(i, j)
            end do
        end do
        call csc_from_triplet(n, n, rr(:k), cc(:k), vv(:k), a, status)
        call checked()
    end subroutine
    subroutine operations()
        integer :: op
        rhs = matmul(dense, exact)
        do op = FORTSPARSE_OP_NORMAL, FORTSPARSE_OP_ADJOINT
            call sparse_solve_op(solver, rhs, answer, op, status)
            call checked()
            call expect(maxval(abs(answer - exact)) < 1e-11_dp, 'ordered N/T/H solve')
            block(:, 1) = rhs; block(:, 2) = -2*rhs
            call sparse_solve_block_op(solver, block, result, op, status)
            call checked()
            call expect(maxval(abs(result(:, 1) - exact)) < 1e-11_dp, 'ordered block solve')
            call expect(maxval(abs(result(:, 2) + 2*exact)) < 1e-11_dp, 'ordered block sign')
        end do
        if (allocated(abar)) deallocate (abar)
        allocate (abar(a%nnz))
        adot = a; adot%val = 0.125_dp*a%val
        tangent = cos(exact); bdot = matmul(dense, tangent) + 0.125_dp*rhs
        call sparse_solve_jvp(solver, adot, exact, bdot, answer, status)
        call checked()
        call expect(maxval(abs(answer - tangent)) < 1e-11_dp, 'ordered implicit tangent')
        seed = matmul(transpose(dense), tangent)
        call sparse_solve_vjp(solver, a, exact, seed, bbar, abar, status, use_primal_factors=.true.)
        call checked()
        call expect(maxval(abs(bbar - tangent)) < 1e-11_dp, 'ordered implicit adjoint')
        do j = 1, n
            do i = a%col_ptr(j), a%col_ptr(j + 1) - 1
                call expect(abs(abar(i) + tangent(a%row_idx(i))*exact(j)) < 1e-11_dp, &
                            'ordered sparse cotangent')
            end do
        end do
    end subroutine
    subroutine checked()
        if (.not. status_ok(status)) then
            print *, status%code, trim(status%msg)
            error stop 'unexpected LDLT AMD status'
        end if
    end subroutine
    subroutine expect(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, message
            error stop 'LDLT AMD oracle'
        end if
    end subroutine
end program
