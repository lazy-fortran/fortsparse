program test_spd_amd_bounded
    use fortsparse
    use fortsparse_native, only: native_backend_t
    use fortsparse_native_amd, only: bounded_amd_available
    use, intrinsic :: iso_c_binding, only: c_int, c_int64_t
    implicit none
    type(csc_t) :: a, changed, adot
    type(native_backend_t) :: factor
    type(sparse_solver_t) :: solver
    type(sparse_factor_stats_t) :: stats, initial
    type(sparse_memory_plan_t) :: plan
    type(fortsparse_status_t) :: status
    integer, parameter :: n = 9
    real(dp) :: dense(n, n), exact(n), rhs(n), x(n), block(n, 3), result(n, &
                                                                         3)
    real(dp) :: zero(n), dx(n), seed(n), bbar(n), abar(3*n - 2)
    real(dp) :: pairing
    integer :: i, j, op, before
    integer(c_int64_t) :: calls, live, peak, allocations, nth
    integer(i8) :: cap, quota
    logical :: reused
    interface
        subroutine begin(nth) bind(C, name='fsparse_spd_probe_begin')
            import c_int64_t
            integer(c_int64_t), value :: nth
        end subroutine
        subroutine finish(count, current, maximum) bind(C, &
                                                        name='fsparse_spd_probe_end')
            import c_int64_t
            integer(c_int64_t), intent(out) :: count, current, maximum
        end subroutine
        integer(c_int) function amd_calls() bind(C, &
                                                 name='fsparse_spd_probe_amd_calls')
            import c_int
        end function
        integer(c_int) function opaque_calls() &
            bind(C, name='fsparse_spd_probe_opaque_calls')
            import c_int
        end function
    end interface
    dense = 0.0_dp
    do i = 1, n
        dense(i, i) = 4.0_dp + 0.125_dp*i
        if (i > 1) then
            dense(i, 1) = -0.25_dp; dense(1, i) = -0.25_dp
        end if
        exact(i) = sin(real(i, dp))
    end do
    call from_dense(dense, a)
    rhs = matmul(dense, exact)
    factor%options%symmetry = FORTSPARSE_SYM_SPD
    factor%options%ordering = FORTSPARSE_ORDER_AMD
    factor%options%memory_policy = FORTSPARSE_MEMORY_LOW
    call begin(0_c_int64_t)
    call factor%factor_real(a, .false., status)
    call finish(calls, live, peak)
    if (bounded_amd_available() /= 1) then
        call expect(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
                    'AMD-off status')
        call expect(amd_calls() == 0 .and. opaque_calls() == 0, &
                    'unavailable AMD invoked')
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options = factor%options
        call sparse_analyze(solver, a, plan, status)
        call expect(status%code == FORTSPARSE_BACKEND_UNAVAILABLE .and. &
                    .not. plan%fits, &
                    'unavailable AMD planning did not fail closed')
        stop
    end if
    call checked()
    call expect(any(factor%perm /= [(i, i=1, n)]), &
                'fixture needs a nontrivial permutation')
    call factor%factor_stats(initial, status)
    call checked()
    call expect(initial%host_current_bytes == live, &
                'independent retained payload')
    call expect(initial%host_peak_bytes == peak, 'independent peak payload')
    call expect(initial%strict_memory_supported, 'strict AMD capability')
    call expect(initial%ordering_used == FORTSPARSE_ORDER_AMD, &
                'actual AMD ordering')
    call expect(opaque_calls() == 0, 'strict path used opaque AMD')
    call planning_oracle()
    cap = initial%host_peak_bytes; allocations = calls
    do op = 0, 2
        call factor%solve_real_op(rhs, x, op, status)
        call checked()
        call expect(maxval(abs(x - exact)) < 1e-12_dp, &
                    'original-A N/T/H oracle')
    end do
    call factor%free()
    factor%options%host_memory_limit_bytes = cap - 1
    call factor%factor_real(a, .false., status)
    call expect(status%code == FORTSPARSE_RESOURCE_LIMIT, &
                'one byte below compact peak')
    call factor%factor_stats(stats, status)
    call checked()
    call expect(stats%host_peak_bytes <= cap - 1, 'refused cap exceeded')
    factor%options%host_memory_limit_bytes = cap
    call factor%factor_real(a, .false., status)
    call checked()
    changed = a; changed%val = 2*a%val
    before = amd_calls()
    call begin(0_c_int64_t)
    call factor%refactor_real(changed, .false., status, reused)
    call finish(calls, live, peak)
    call checked()
    call expect(reused .and. amd_calls() == before, &
                'same-pattern AMD was not reused')
    call expect(calls == 0, 'same-pattern numeric refactor allocated')
    call factor%solve_real(2*rhs, x, status)
    call checked()
    call expect(maxval(abs(x - exact)) < 1e-12_dp, &
                'refactor original-A oracle')
    ! Each controlled allocator failure must leave an honest live ledger and
    ! unusable factor, then permit a fresh successful recovery.
    factor%options%host_memory_limit_bytes = 0_i8
    do nth = 1, allocations
        call factor%free()
        call begin(nth)
        call factor%factor_real(a, .false., status)
        call finish(calls, live, peak)
        call expect(status%code == FORTSPARSE_OUT_OF_MEMORY, &
                    'injected allocation status')
        call factor%factor_stats(stats, status)
        call checked()
        if (stats%host_current_bytes /= live) print *, &
            'allocation index/calls/live/ledger', &
            nth, calls, live, stats%host_current_bytes
        call expect(stats%host_current_bytes == live, &
                    'partial allocation live ledger')
        call factor%solve_real(rhs, x, status)
        call expect(status%code == FORTSPARSE_NOT_FACTORED, &
                    'failed factor remained usable')
    end do
    call factor%factor_real(a, .false., status)
    call checked()
    call factor%solve_real(rhs, x, status)
    call checked()
    call expect(maxval(abs(x - exact)) < 1e-12_dp, &
                'allocation-failure recovery')
    call factor%free()
    ! Public ordered planning reuses actual AMD and preserves retained factors.
    solver%backend_id = FORTSPARSE_BACKEND_NATIVE
    solver%options%symmetry = FORTSPARSE_SYM_SPD
    solver%options%ordering = FORTSPARSE_ORDER_AMD
    solver%options%memory_policy = FORTSPARSE_MEMORY_LOW
    solver%refine = .true.
    solver%options%host_memory_limit_bytes = 1000000_i8
    call sparse_analyze(solver, a, plan, status)
    call checked()
    call expect(plan%fits .and. plan%prediction_quality == MEMORY_PLAN_EXACT_LAYOUT, &
                'ordered exact plan unavailable')
    call sparse_factor(solver, a, status)
    call checked()
    call sparse_analyze(solver, a, plan, status)
    call checked()
    call expect(plan%fits .and. solver%factored .and. &
                plan%prediction_quality == MEMORY_PLAN_BOUND, &
                'ordered planning changed retained factors')
    call sparse_solve(solver, rhs, x, status)
    call checked()
    call expect(maxval(abs(x - exact)) < 1e-12_dp, &
                'public refined AMD solve')
    zero = 0.0_dp; seed = cos(exact); adot = a; adot%val = 0.125_dp*a%val
    call sparse_solve_jvp(solver, adot, x, zero, dx, status)
    call checked()
    call expect(maxval(abs(dx + 0.125_dp*x)) < 1e-12_dp, 'AMD implicit JVP')
    call sparse_solve_vjp(solver, a, x, seed, bbar, abar, status, &
                          use_primal_factors=.true.)
    call checked()
    pairing = dot_product(seed, dx) - dot_product(abar, adot%val)
    call expect(abs(pairing) < 1e-12_dp, 'AMD VJP pairing')
    do j = 1, 3
        block(:, j) = j*rhs
    end do
    call sparse_solve_block(solver, block, result, status)
    call checked()
    call expect(maxval(abs(matmul(dense, result) - block)) < 1e-12_dp, &
                'public block original-A residual')
    call sparse_destroy(solver)
    call refinement_oracle()
    call blocked_oracle()
    call large_count_oracle()
    call default_order_oracle()
    print '(a)', &
        'Strict AMD SPD allocation, recovery, permutation, reuse and refinement pass'
contains
    subroutine planning_oracle()
        type(sparse_solver_t) :: request
        type(sparse_memory_plan_t) :: predicted
        type(sparse_factor_stats_t) :: retained_stats, unchanged
        integer :: backend, refinement, row, col, pivot
        logical :: filled(n, n)
        integer(i8) :: expected_peak, retained_bytes
        integer(c_int64_t) :: measured_calls, measured_live, measured_peak, &
            failures, failed_index
        do backend = FORTSPARSE_BACKEND_NATIVE, FORTSPARSE_BACKEND_NATIVE_REFERENCE
            if (backend /= FORTSPARSE_BACKEND_NATIVE .and. &
                backend /= FORTSPARSE_BACKEND_NATIVE_REFERENCE) cycle
            do refinement = 0, 1
                request%backend_id = backend
                request%options%symmetry = FORTSPARSE_SYM_SPD
                request%options%ordering = FORTSPARSE_ORDER_AMD
                request%options%memory_policy = FORTSPARSE_MEMORY_LOW
                request%refine = refinement == 1
                call sparse_analyze(request, a, predicted, status, nrhs=7)
                call checked()
                call expect(predicted%fits .and. predicted%max_rhs_batch == 7, &
                    'AMD RHS plan')
                call expect(predicted%prediction_quality == &
                    MEMORY_PLAN_EXACT_LAYOUT, 'AMD layout quality')
                expected_peak = predicted%peak_host_bytes
                retained_bytes = predicted%factor_host_bytes + &
                    predicted%persistent_host_bytes
                ! Use concrete production backend: excludes the polymorphic object
                ! descriptor allocated by the public factory, outside payload scope.
                call factor%free()
                factor%reference_only = backend == FORTSPARSE_BACKEND_NATIVE_REFERENCE
                factor%options%host_memory_limit_bytes = expected_peak
                call begin(0_c_int64_t)
                call factor%factor_real(a, request%refine, status)
                call finish(measured_calls, measured_live, measured_peak)
                call checked()
                call expect(measured_peak == expected_peak, &
                    'predicted vs allocator peak')
                call expect(measured_live == retained_bytes, &
                    'predicted vs allocator retained')
                filled = .false.
                do col = 1, n
                    do row = 1, n
                        filled(row, col) = dense(factor%perm(row), &
                            factor%perm(col)) /= 0.0_dp
                    end do
                end do
                do pivot = 1, n
                    do col = pivot + 1, n
                        do row = pivot + 1, n
                            if (filled(row, pivot) .and. filled(col, &
                                pivot)) filled(row, col) = .true.
                        end do
                    end do
                end do
                row = 0
                do col = 1, n
                    row = row + count(filled(col:n, col))
                end do
                call expect(predicted%symbolic_nnz == row, &
                    'independent permuted dense fill count')
                call factor%free()
                factor%options%host_memory_limit_bytes = expected_peak - 1_i8
                call factor%factor_real(a, request%refine, status)
                call expect(status%code == FORTSPARSE_RESOURCE_LIMIT, &
                    'planned cap minus one factor')
                call factor%free()
                request%options%host_memory_limit_bytes = expected_peak - 1_i8
                call sparse_analyze(request, a, predicted, status)
                call expect(status%code == FORTSPARSE_RESOURCE_LIMIT .and. &
                    .not. predicted%fits, &
                            'planned cap minus one analysis')
                request%options%host_memory_limit_bytes = expected_peak
                call sparse_factor(request, a, status)
                call checked()
                call sparse_factor_stats(request, retained_stats, status)
                call checked()
                call sparse_analyze(request, a, predicted, status)
                call expect(status%code == FORTSPARSE_RESOURCE_LIMIT, &
                    'retained overlap gate')
                call sparse_factor_stats(request, unchanged, status)
                call checked()
                call expect(unchanged%host_current_bytes == &
                    retained_stats%host_current_bytes .and. unchanged%host_peak_bytes == &
                    retained_stats%host_peak_bytes, 'analysis mutated factor ledger')
                call sparse_solve(request, rhs, x, status)
                call checked()
                call expect(maxval(abs(x - exact)) < 1e-12_dp, &
                    'failed analysis preserved factor')
                request%options%host_memory_limit_bytes = expected_peak + retained_bytes
                call begin(0_c_int64_t)
                call sparse_analyze(request, a, predicted, status)
                call finish(measured_calls, measured_live, measured_peak)
                call checked()
                call expect(measured_live == 0_i8, 'analysis leaked ordered temporary')
                call expect(measured_peak + retained_bytes <= &
                            request%options%host_memory_limit_bytes, &
                                'old plus actual analysis peak exceeds cap')
                call expect(predicted%fits .and. &
                    predicted%prediction_quality == MEMORY_PLAN_BOUND, &
                            'retained conservative bound')
                ! The first 13 calls are explicit payload arrays: AMD outputs
                ! and arena (3), ordered CSC/pattern/cursor (7), count scratch
                ! (3). GNU finalization then allocates runtime descriptors
                ! without STAT; these are outside solver payload accounting.
                failures = 13_c_int64_t
                call expect(measured_calls >= failures, &
                    'missing payload allocation observations')
                do failed_index = 1, failures
                    call begin(failed_index)
                    call sparse_analyze(request, a, predicted, status)
                    call finish(measured_calls, measured_live, measured_peak)
                    call expect(status%code == FORTSPARSE_OUT_OF_MEMORY .and. &
                        .not. predicted%fits, &
                                'planning allocation failure status')
                    call expect(measured_live == 0_i8, &
                        'failed analysis leaked temporary')
                end do
                call sparse_solve(request, rhs, x, status)
                call checked()
                call expect(maxval(abs(x - exact)) < 1e-12_dp, &
                    'OOM analysis preserved old factor')
                call sparse_destroy(request)
            end do
        end do
        request%backend_id = FORTSPARSE_BACKEND_NATIVE_REFERENCE
        request%options = factor%options
        request%options%memory_policy = FORTSPARSE_MEMORY_AUTO
        call sparse_analyze(request, a, predicted, status)
        call expect(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
            'opaque reference AMD plan')
        request%options%memory_policy = FORTSPARSE_MEMORY_LOW
        request%options%host_memory_limit_bytes = 0_i8
        call sparse_analyze(request, a, predicted, status, nrhs=-1)
        call expect(status%code == FORTSPARSE_INVALID_MATRIX, &
            'negative ordered RHS count')
        request%options%host_memory_limit_bytes = -1_i8
        call sparse_analyze(request, a, predicted, status)
        call expect(status%code == FORTSPARSE_INVALID_MATRIX, 'negative ordered cap')
        request%options%host_memory_limit_bytes = 0_i8
        request%options%ordering = FORTSPARSE_ORDER_ND
        call sparse_analyze(request, a, predicted, status)
        call expect(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
            'ND planning unsupported')
        request%options%ordering = FORTSPARSE_ORDER_AMD
        request%options%device = 0
        call sparse_analyze(request, a, predicted, status)
        call expect(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
            'GPU planning unsupported')
        factor%reference_only = .false.
        factor%options%host_memory_limit_bytes = 0_i8
        call factor%factor_real(a, .false., status)
        call checked()
    end subroutine planning_oracle

    subroutine checked()
        if (.not. status_ok(status)) then
            print *, status%code, trim(status%msg)
            error stop 'unexpected bounded SPD status'
        end if
    end subroutine
    subroutine expect(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, message
            error stop 'bounded SPD oracle'
        end if
    end subroutine
    subroutine from_dense(matrix, csc)
        real(dp), intent(in) :: matrix(:, :)
        type(csc_t), intent(out) :: csc
        integer :: m, i, j, p
        m = size(matrix, 1)
        csc%nrow = m; csc%ncol = m; csc%nnz = count(matrix /= 0.0_dp)
        allocate (csc%col_ptr(m + 1), csc%row_idx(csc%nnz), csc%val(csc%nnz))
        p = 0
        do j = 1, m
            csc%col_ptr(j) = p + 1
            do i = 1, m
                if (matrix(i, j) == 0.0_dp) cycle
                p = p + 1; csc%row_idx(p) = i; csc%val(p) = matrix(i, j)
            end do
        end do
        csc%col_ptr(m + 1) = p + 1
    end subroutine
    subroutine refinement_oracle()
        type(csc_t) :: updated
        type(sparse_factor_stats_t) :: refined
        integer(c_int64_t) :: refinement_allocations
        quota = 12_i8*a%nnz + 20_i8*n + 4_i8
        factor%options%host_memory_limit_bytes = quota - 1
        before = amd_calls()
        call begin(0_c_int64_t)
        call factor%factor_real(a, .true., status)
        call finish(calls, live, peak)
        call expect(status%code == FORTSPARSE_RESOURCE_LIMIT, &
                    'refinement quota gate')
        call expect(calls == 0 .and. amd_calls() == before, &
                    'quota gate allocated or ordered')
        factor%options%host_memory_limit_bytes = 0_i8
        call begin(0_c_int64_t)
        call factor%factor_real(a, .true., status)
        call finish(calls, live, peak)
        call checked()
        refinement_allocations = calls
        call factor%factor_stats(refined, status)
        call checked()
        call expect(refined%host_current_bytes == live, &
                    'refinement live payload')
        call expect(refined%host_peak_bytes == peak, &
                    'refinement peak payload')
        call expect(live == initial%host_current_bytes + quota, &
                    'exact refinement quota')
        factor%lv(factor%lp(1)) = factor%lv(factor%lp(1))*(1 + 1e-4_dp)
        do op = 0, 2
            call begin(0_c_int64_t)
            call factor%solve_real_op(rhs, x, op, status)
            call finish(calls, live, peak)
            call checked()
            call expect(calls == 0, 'permuted refinement correction allocated')
            call expect(factor%stats%refinement_steps > 0, &
                        'correction oracle did not correct')
            call expect(maxval(abs(x - exact)) < 1e-12_dp, &
                        'permuted correction original-A oracle')
        end do
        updated = a; updated%val = 2*a%val
        before = amd_calls()
        call factor%refactor_real(updated, .true., status, reused)
        call checked()
        call expect(reused .and. amd_calls() == before, &
                    'refined pattern reuse')
        factor%lv(factor%lp(1)) = factor%lv(factor%lp(1))*(1 + 1e-4_dp)
        call factor%solve_real(2*rhs, x, status)
        call checked()
        call expect(maxval(abs(x - exact)) < 1e-12_dp, &
                    'refined refactor kept stale original values')
        call factor%free()
        do nth = 1, refinement_allocations
            call begin(nth)
            call factor%factor_real(a, .true., status)
            call finish(calls, live, peak)
            call expect(status%code == FORTSPARSE_OUT_OF_MEMORY, &
                        'refined injected allocation status')
            call factor%factor_stats(stats, status)
            call checked()
            call expect(stats%host_current_bytes == live, &
                        'refined partial allocation live ledger')
            call factor%free()
        end do
        factor%options%host_memory_limit_bytes = 0_i8
        factor%options%memory_policy = FORTSPARSE_MEMORY_LOW
        call factor%factor_real(a, .false., status)
        call checked()
        block
            real(dp) :: altered(n, n)
            type(csc_t) :: pattern
            altered = dense; altered(2, 3) = 0.125_dp; altered(3, 2) = 0.125_dp
            call from_dense(altered, pattern)
            before = amd_calls()
            call factor%refactor_real(pattern, .false., status, reused)
            call checked()
            call expect(.not. reused .and. amd_calls() == before + 1, &
                        'changed graph reused AMD')
        end block
        updated = a
        do j = 1, n
            do i = updated%col_ptr(j), updated%col_ptr(j + 1) - 1
                if (updated%row_idx(i) == j) updated%val(i) = -4.0_dp
            end do
        end do
        call factor%refactor_real(updated, .false., status, reused)
        call expect(status%code == FORTSPARSE_NOT_POSITIVE_DEFINITE, &
                    'AMD nonpositive pivot status')
        call factor%solve_real(rhs, x, status)
        call expect(status%code == FORTSPARSE_NOT_FACTORED, &
                    'nonpositive factor remained usable')
        call factor%refactor_real(a, .false., status, reused)
        call checked()
        call factor%solve_real(rhs, x, status)
        call checked()
        call expect(maxval(abs(x - exact)) < 1e-12_dp, &
                    'AMD pivot-failure recovery')
        factor%options%ordering = FORTSPARSE_ORDER_DEFAULT
        call factor%refactor_real(a, .false., status, reused)
        call checked()
        call expect(.not. reused, 'changed ordering reused old analysis')
        factor%options%ordering = FORTSPARSE_ORDER_ND
        call factor%factor_real(a, .false., status)
        call expect(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
                    'strict ND accepted')
        factor%options%ordering = FORTSPARSE_ORDER_AMD
        factor%options%device = 0
        call factor%factor_real(a, .false., status)
        call expect(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
                    'strict CUDA accepted')
        factor%options%device = -1
        call factor%free()
    end subroutine
    subroutine blocked_oracle()
        integer, parameter :: m = 64
        type(native_backend_t) :: blocked
        type(csc_t) :: full
        real(dp) :: matrix(m, m), truth(m), b(m), result(m)
        matrix = 0.01_dp
        do i = 1, m
            matrix(i, i) = 2.0_dp; truth(i) = cos(real(i, dp))
        end do
        call from_dense(matrix, full)
        b = matmul(matrix, truth)
        blocked%options%symmetry = FORTSPARSE_SYM_SPD
        blocked%options%ordering = FORTSPARSE_ORDER_AMD
        blocked%options%host_memory_limit_bytes = 10000000_i8
        blocked%options%nthreads = 4
        call begin(0_c_int64_t)
        call blocked%factor_real(full, .false., status)
        call finish(calls, live, peak)
        call checked()
        call blocked%factor_stats(stats, status)
        call checked()
        call expect(.not. stats%compact_memory_used, &
                    'dense fixture missed blocked schedule')
        call expect(stats%host_current_bytes == live, &
                    'blocked independent retained payload')
        call expect(stats%host_peak_bytes == peak, &
                    'blocked independent peak payload')
        call expect(.not. allocated(blocked%source), &
                    'blocked retained source map')
        call expect(.not. allocated(blocked%ordered%val), &
                    'blocked retained ordered values')
        call blocked%solve_real(b, result, status)
        call checked()
        call expect(maxval(abs(result - truth)) < 1e-12_dp, &
                    'blocked original-A solve')
        call blocked%free()
        ! Dense symbolic fill is analytically m*(m+1)/2. This cap fits the
        ! compact ordered layout but cannot coexist with the panel build.
        cap = 12_i8*(m*(m + 1_i8)/2_i8) + 20_i8*m*m + 48_i8*m + 12_i8
        blocked%options%host_memory_limit_bytes = cap
        call begin(0_c_int64_t)
        call blocked%factor_real(full, .false., status)
        call finish(calls, live, peak)
        call checked()
        call blocked%factor_stats(stats, status)
        call checked()
        call expect(stats%compact_memory_used, 'AMD cap did not use compact fallback')
        call expect(stats%host_current_bytes == cap, 'AMD fallback exact layout')
        call expect(stats%host_current_bytes == live, 'fallback actual live payload')
        call expect(stats%host_peak_bytes == peak .and. peak <= cap, &
                    'fallback actual peak payload')
        call blocked%solve_real(b, result, status)
        call checked()
        call expect(maxval(abs(result - truth)) < 1e-12_dp, 'fallback original-A solve')
        call blocked%free()
    end subroutine
    subroutine large_count_oracle()
        integer, parameter :: m = 4096
        type(native_backend_t) :: large
        type(sparse_solver_t) :: request
        type(sparse_memory_plan_t) :: forecast
        type(csc_t) :: tridiagonal
        real(dp) :: b(m), truth(m), result(m)
        integer :: p, i, j
        tridiagonal%nrow = m; tridiagonal%ncol = m; tridiagonal%nnz = 3*m - 2
        allocate (tridiagonal%col_ptr(m + 1), tridiagonal%row_idx(3*m - 2), &
                  tridiagonal%val(3*m - 2))
        p = 0
        do j = 1, m
            tridiagonal%col_ptr(j) = p + 1
            do i = max(1, j - 1), min(m, j + 1)
                p = p + 1; tridiagonal%row_idx(p) = i; tridiagonal%val(p) = -1.0_dp
                if (i == j) tridiagonal%val(p) = 4.0_dp
            end do
            truth(j) = sin(real(j, dp))
        end do
        tridiagonal%col_ptr(m + 1) = p + 1
        b = 4*truth
        b(2:) = b(2:) - truth(:m - 1); b(:m - 1) = b(:m - 1) - truth(2:)
        large%options%symmetry = FORTSPARSE_SYM_SPD
        large%options%ordering = FORTSPARSE_ORDER_AMD
        large%options%memory_policy = FORTSPARSE_MEMORY_LOW
        request%backend_id = FORTSPARSE_BACKEND_NATIVE_REFERENCE
        request%options = large%options
        request%refine = .false.
        call sparse_analyze(request, tridiagonal, forecast, status)
        call checked()
        call expect(forecast%workspace_host_bytes > 0_i8, &
            'large count scratch not forecast')
        large%options%host_memory_limit_bytes = forecast%peak_host_bytes
        call begin(0_c_int64_t)
        call large%factor_real(tridiagonal, .false., status)
        call finish(calls, live, peak)
        call checked()
        call large%factor_stats(stats, status)
        call checked()
        call expect(stats%host_current_bytes == live, &
                    'large independent retained payload')
        call expect(stats%host_peak_bytes == peak, &
                    'large graph-count workspace peak')
        call expect(forecast%peak_host_bytes == peak, &
            'large predicted actual allocator peak')
        call expect(forecast%symbolic_nnz == 2_i8*m - 1_i8, &
            'large analytic tridiagonal fill')
        call expect(forecast%factor_host_bytes + &
                    forecast%persistent_host_bytes == live, &
                    'large predicted actual retained payload')
        call large%solve_real(b, result, status)
        call checked()
        call expect(maxval(abs(result - truth)) < 1e-12_dp, &
                    'large original tridiagonal oracle')
        call large%free()
        large%options%host_memory_limit_bytes = forecast%peak_host_bytes - 1_i8
        call large%factor_real(tridiagonal, .false., status)
        call expect(status%code == FORTSPARSE_RESOURCE_LIMIT, &
            'large peak minus one cap')
        call large%free()
        request%options%host_memory_limit_bytes = forecast%peak_host_bytes - 1_i8
        call sparse_analyze(request, tridiagonal, forecast, status)
        call expect(status%code == FORTSPARSE_RESOURCE_LIMIT, &
            'large planning peak minus one')
    end subroutine

    subroutine default_order_oracle()
        ! DEFAULT on a row-numbered 65 x 65 grid (n = 4225, bandwidth 65)
        ! selects AMD under LOW_MEMORY, refinement and a plain cap, and the
        ! exact plan and ledger follow that order.
        integer, parameter :: mx = 65, my = 65, m = mx*my
        type(native_backend_t) :: strip
        type(sparse_solver_t) :: request
        type(sparse_memory_plan_t) :: natural, ordered, default
        type(csc_t) :: lap
        real(dp) :: b(m), truth(m), result(m)
        integer :: p, i, j, refinement
        lap%nrow = m; lap%ncol = m; lap%nnz = 5*m - 2*mx - 2*my
        allocate (lap%col_ptr(m + 1), lap%row_idx(lap%nnz), lap%val(lap%nnz))
        p = 0
        do j = 1, m
            lap%col_ptr(j) = p + 1
            do i = j - mx, j + mx
                if (i < 1 .or. i > m) cycle
                if (abs(i - j) == 1 .and. &
                    (min(i, j) - 1)/mx /= (max(i, j) - 1)/mx) cycle
                if (abs(i - j) /= 0 .and. abs(i - j) /= 1 .and. abs(i - j) /= mx) cycle
                p = p + 1; lap%row_idx(p) = i; lap%val(p) = -1.0_dp
                if (i == j) lap%val(p) = 4.0_dp
            end do
            truth(j) = sin(real(j, dp))
        end do
        lap%col_ptr(m + 1) = p + 1
        call expect(p == lap%nnz, 'grid Laplacian entry count')
        b = csc_matvec(lap, truth)
        do refinement = 0, 1
            request%backend_id = FORTSPARSE_BACKEND_NATIVE_REFERENCE
            request%options%symmetry = FORTSPARSE_SYM_SPD
            request%options%memory_policy = FORTSPARSE_MEMORY_LOW
            request%options%ordering = FORTSPARSE_ORDER_DEFAULT
            request%refine = refinement == 1
            call sparse_analyze(request, lap, natural, status)
            call checked()
            request%backend_id = FORTSPARSE_BACKEND_NATIVE
            request%options%ordering = FORTSPARSE_ORDER_AMD
            call sparse_analyze(request, lap, ordered, status)
            call checked()
            request%options%ordering = FORTSPARSE_ORDER_DEFAULT
            call sparse_analyze(request, lap, default, status)
            call checked()
            call expect(default%symbolic_nnz == ordered%symbolic_nnz .and. &
                default%peak_host_bytes == ordered%peak_host_bytes .and. &
                default%factor_host_bytes == ordered%factor_host_bytes .and. &
                default%persistent_host_bytes == ordered%persistent_host_bytes, &
                'DEFAULT plan is not the AMD plan')
            call expect(default%symbolic_nnz < natural%symbolic_nnz/2, &
                'DEFAULT plan fill not below natural order')
            strip%options = request%options
            strip%options%host_memory_limit_bytes = default%peak_host_bytes
            before = amd_calls()
            call begin(0_c_int64_t)
            call strip%factor_real(lap, request%refine, status)
            call finish(calls, live, peak)
            call checked()
            call strip%factor_stats(stats, status)
            call checked()
            call expect(stats%ordering_used == FORTSPARSE_ORDER_AMD .and. &
                amd_calls() == before + 1, 'bounded DEFAULT did not order by AMD')
            call expect(stats%host_peak_bytes == peak .and. &
                default%peak_host_bytes == peak, 'DEFAULT predicted vs allocator peak')
            call expect(stats%host_current_bytes == live .and. &
                default%factor_host_bytes + default%persistent_host_bytes == live, &
                'DEFAULT predicted vs allocator retained')
            call strip%solve_real(b, result, status)
            call checked()
            call expect(maxval(abs(result - truth)) < 1e-12_dp, &
                'bounded DEFAULT original-A oracle')
            call strip%free()
            strip%options%host_memory_limit_bytes = default%peak_host_bytes - 1_i8
            call strip%factor_real(lap, request%refine, status)
            call expect(status%code == FORTSPARSE_RESOURCE_LIMIT, &
                'bounded DEFAULT peak minus one cap')
            call strip%free()
        end do
        ! A plain cap keeps the blocked schedule and its exact ledger.
        strip%options = request%options
        strip%options%memory_policy = FORTSPARSE_MEMORY_AUTO
        strip%options%host_memory_limit_bytes = 1000000000_i8
        call begin(0_c_int64_t)
        call strip%factor_real(lap, .false., status)
        call finish(calls, live, peak)
        call checked()
        call strip%factor_stats(stats, status)
        call checked()
        call expect(stats%ordering_used == FORTSPARSE_ORDER_AMD, &
            'capped DEFAULT did not order by AMD')
        call expect(stats%host_current_bytes == live .and. &
            stats%host_peak_bytes == peak, 'capped DEFAULT ledger')
        call strip%solve_real(b, result, status)
        call checked()
        call expect(maxval(abs(result - truth)) < 1e-12_dp, &
            'capped DEFAULT original-A oracle')
        call strip%free()
        call sparse_destroy(request)
    end subroutine default_order_oracle

end program test_spd_amd_bounded
