program test_matching_api
    use fortsparse
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    integer, parameter::n = 5
    type(sparse_solver_t)::solver
    type(csc_t)::a, adot, bad
    type(csc_z_t)::z, zdot
    type(fortsparse_status_t)::status
    type(sparse_factor_stats_t)::stats, fresh
    real(dp)::dense(n, n), exact(n), rhs(n, 2), answer(n, 2), dx(n), zero(n), seed(n), &
               bbar(n)
    real(dp), allocatable::abar(:)
    complex(dp)::zdense(n, n), zexact(n), zrhs(n, 2), zanswer(n, 2), zdx(n), &
                  zzero(n), zbbar(n)
    complex(dp), allocatable::zabar(:)
    integer::row(n*n), col(n*n), i, j, k, op, cycle, ord
    integer(i8)::cap, old_live
    logical::reused
    interface
        integer(c_int) function amd_available() &
            bind(C, name='fsparse_native_order_amd_available')
            import c_int
        end function
        integer(c_int) function bounded_available() &
            bind(C, name='fsparse_native_amd_bounded_available')
            import c_int
        end function
    end interface
    k = 0
    do j = 1, n
        do i = 1, n
            k = k + 1; row(k) = i; col(k) = j
        end do
    end do
    allocate (abar(n*n), zabar(n*n))
    exact = [(sin(real(i, dp)), i=1, n)]; zexact = cmplx(exact, &
 0.25_dp*exact + 0.1_dp, dp)
    solver%backend_id = FORTSPARSE_BACKEND_NATIVE
    solver%options%structural_matching = .true.
    solver%options%pivot_tolerance = 0.1_dp
    do ord = FORTSPARSE_ORDER_DEFAULT, FORTSPARSE_ORDER_ND
        if (ord == 1 .and. amd_available() /= 1) cycle
        if (.not. native_ordering_available(ord)) cycle
        solver%options%ordering = ord
        do cycle = 1, 2
            dense = 0
            do i = 1, n
                j = mod(i + cycle - 1, n) + 1
                dense(i, j) = 4 + i
            end do
            call csc_from_triplet(n, n, row, col, reshape(dense, [n*n]), a, status)
            call checked('real canonical input')
            solver%options%structural_matching = .false.
            call sparse_factor(solver, a, status)
            call checked('unmatched general ordering')
            call real_operations()
            if (ord == FORTSPARSE_ORDER_ND) then
                call sparse_factor_stats(solver, stats, status)
                call checked('ND statistics')
                call check(stats%ordering_used == FORTSPARSE_ORDER_ND, &
                           'actual ND ordering is reported')
                call check(stats%host_peak_bytes == -1 .and. &
                           .not. stats%strict_memory_supported, &
                           'opaque METIS peak remains unknown')
            end if
            solver%options%structural_matching = .true.
            if (cycle == 1) then
                call sparse_factor(solver, a, status)
            else
                call sparse_refactor(solver, a, status, reused)
            end if
            call checked('real matching factor/refactor')
            if (cycle == 2) then
                call check(.not. reused, 'matching/order recomputed after zero changes')
            end if
            call real_operations()
            ! Original stored pattern stays full; nonzero graph changes again
            ! for complex. A new match is required and adjoint must conjugate.
            zdense = cmplx(dense, 0.5_dp*dense, dp)
            call csc_from_triplet(n, n, row, col, reshape(zdense, [n*n]), z, &
                                  status)
            solver%options%structural_matching = .false.
            call sparse_factor(solver, z, status)
            call checked('unmatched complex ordering')
            zrhs(:, 1) = matmul(zdense, zexact)
            zrhs(:, 2) = 2*zrhs(:, 1)
            call sparse_solve_block_op(solver, zrhs, zanswer, 0, status)
            call checked('unmatched complex block solve')
            call check(maxval(abs(zanswer(:, 1) - zexact)) < 1e-12_dp, &
                       'unmatched complex original matrix oracle')
            solver%options%structural_matching = .true.
            call sparse_factor(solver, z, status)
            call checked('complex matching')
            do op = 0, 2
                select case (op)
                case (0)
                    zrhs(:, 1) = matmul(zdense, zexact)
                case (1)
                    zrhs(:, 1) = matmul(transpose(zdense), zexact)
                case (2)
                    zrhs(:, 1) = matmul(transpose(conjg(zdense)), zexact)
                end select
                zrhs(:, 2) = 2*zrhs(:, 1)
                call sparse_solve_block_op(solver, zrhs, zanswer, op, status)
                call checked('complex N/T/H block solve')
                call check(maxval(abs(zanswer(:, 1) - zexact)) < 1e-12_dp, &
                           'complex independent known vector')
                call check(maxval(abs(zanswer(:, 2) - 2*zexact)) < 1e-12_dp, &
                           'complex second RHS')
            end do
            ! Raw complex path borrows original input, with the same policy.
            call sparse_factor(solver, n, n, z%nnz, z%col_ptr, z%row_idx, &
                               z%val, status)
            call checked('raw complex matching')
            zdot = z; zdot%val = 0.125_dp*z%val; zzero = cmplx(0.0_dp, 0.0_dp, dp)
            call sparse_solve_jvp(solver, zdot, zexact, zzero, zdx, status)
            call checked('complex matching implicit JVP')
            call check(maxval(abs(zdx + 0.125_dp*zexact)) < 1e-12_dp, &
                       'complex JVP equation oracle')
            call sparse_solve_vjp(solver, z, zexact, zexact, zbbar, zabar, status, &
                                  use_primal_factors=.true.)
            call checked('complex matching retained-factor VJP')
            call check(abs(real(dot_product(zexact, zdx) - dot_product(zabar, &
                                                                       zdot%val), dp)) &
                       < 1e-12_dp, 'complex real-pairing oracle')
        end do
        call sparse_destroy(solver)
    end do
    solver%options%ordering = FORTSPARSE_ORDER_DEFAULT
    solver%options%memory_policy = FORTSPARSE_MEMORY_LOW
    call sparse_factor(solver, a, status)
    call checked('LOW natural matching')
    call sparse_factor_stats(solver, fresh, status)
    cap = fresh%host_peak_bytes
    call check(fresh%strict_memory_supported .and. cap > 0, &
               'matching has exact natural payload statistics')
    call sparse_destroy(solver)
    solver%options%host_memory_limit_bytes = cap
    call sparse_factor(solver, a, status)
    call checked('exact fresh peak cap')
    call sparse_factor_stats(solver, stats, status)
    call check(stats%host_peak_bytes <= cap, 'fresh cap respected')
    old_live = stats%host_current_bytes
    ! This cap must fail before allocating the complete matching workspace,
    ! even though the previous factor remains allocated during analysis.
    solver%options%host_memory_limit_bytes = old_live + 28_i8*n - 1
    call sparse_refactor(solver, a, status, reused)
    call check(status%code == FORTSPARSE_RESOURCE_LIMIT, &
               'old factors plus matching workspace cap')
    call sparse_factor_stats(solver, stats, status)
    call check(stats%host_current_bytes == old_live, &
               'analysis cap failure releases transient arrays')
    call check(stats%budget_required_bytes == 28_i8*n, &
               'matching reservation diagnostic')
    solver%options%host_memory_limit_bytes = 0
    call sparse_factor(solver, a, status)
    call checked('recovery after cap refusal')
    call real_operations()
    ! Hall deficient with all stored columns present; zeros are not edges.
    bad = a; bad%val = 0
    do j = 1, n
        do k = bad%col_ptr(j), bad%col_ptr(j + 1) - 1
            if (bad%row_idx(k) == 1) bad%val(k) = 1
        end do
    end do
    call sparse_factor(solver, bad, status)
    call check(status%code == FORTSPARSE_SINGULAR .and. .not. solver%factored, &
               'Hall failure invalidates solve')
    call sparse_factor(solver, a, status)
    call checked('recovery after Hall failure')
    call real_operations()
    solver%refine = .true.
    call sparse_factor(solver, n, n, a%nnz, a%col_ptr, a%row_idx, a%val, status)
    call checked('raw matching with retained original-matrix refinement')
    call real_operations()
    adot = a; adot%val = 0.125_dp*a%val; zero = 0; seed = exact
    call sparse_solve_jvp(solver, adot, exact, zero, dx, status)
    call checked('matched implicit JVP')
    call check(maxval(abs(dx + 0.125_dp*exact)) < 1e-12_dp, &
               'JVP equation oracle')
    call sparse_solve_vjp(solver, a, exact, seed, bbar, abar, status, &
                          use_primal_factors=.true.)
    call checked('matched implicit VJP')
    call check(abs(dot_product(seed, dx) - dot_product(abar, adot%val)) &
               < 1e-12_dp, 'adjoint pairing')
    call sparse_destroy(solver)
    if (amd_available() == 1 .and. bounded_available() == 1) then
        solver%refine = .false.
        solver%options%ordering = FORTSPARSE_ORDER_AMD
        call sparse_factor(solver, a, status)
        call checked('LOW matched AMD with bounded workspace')
        call sparse_factor_stats(solver, fresh, status)
        call check(fresh%strict_memory_supported, &
                   'bounded matching AMD statistics')
        cap = fresh%host_peak_bytes
        call sparse_destroy(solver)
        solver%options%host_memory_limit_bytes = cap
        call sparse_factor(solver, a, status)
        call checked('bounded AMD exact cap includes transformed graph')
        call real_operations()
        call sparse_destroy(solver)
    end if
    solver%options%ordering = FORTSPARSE_ORDER_ND
    solver%options%host_memory_limit_bytes = 1024_i8*1024_i8
    solver%options%memory_policy = FORTSPARSE_MEMORY_AUTO
    call sparse_factor(solver, a, status)
    call check(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
               'ND rejects a strict cap before opaque ordering')
    solver%options%host_memory_limit_bytes = 0
    solver%options%memory_policy = FORTSPARSE_MEMORY_LOW
    call sparse_factor(solver, a, status)
    call check(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
               'ND rejects LOW before opaque ordering')
    solver%options%memory_policy = FORTSPARSE_MEMORY_AUTO
    if (.not. native_ordering_available(FORTSPARSE_ORDER_ND)) then
        call sparse_factor(solver, a, status)
        call check(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
                   'METIS-off explicit ordering status')
    end if
    solver%options%host_memory_limit_bytes = 0
    solver%options%ordering = FORTSPARSE_ORDER_DEFAULT
    solver%options%structural_matching = .false.
    call sparse_factor(solver,a,status)
    call checked('switch back to default matching-disabled policy')
    call real_operations()
    solver%options%structural_matching = .true.
    solver%options%symmetry = FORTSPARSE_SYM_SPD
    call sparse_factor(solver, a, status)
    call check(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
               'SPD matching explicitly unsupported')
    solver%options%symmetry = FORTSPARSE_SYM_GENERAL
    solver%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC
    call sparse_factor(solver, a, status)
    call check(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
               'external matching explicitly unsupported')
contains
    subroutine real_operations()
        integer::op
        do op = 0, 2
            if (op == 0) then
                rhs(:, 1) = matmul(dense, exact)
            else
                rhs(:, 1) = matmul(transpose(dense), exact)
            end if
            rhs(:, 2) = 2*rhs(:, 1)
            call sparse_solve_block_op(solver, rhs, answer, op, status)
            call checked('real N/T/H block solve')
            call check(maxval(abs(answer(:, 1) - exact)) < 1e-12_dp, &
                       'real original dense oracle')
            call check(maxval(abs(answer(:, 2) - 2*exact)) < 1e-12_dp, &
                       'real second RHS')
        end do
    end subroutine
    subroutine checked(message)
        character(*), intent(in)::message
        if (.not. status_ok(status)) print *, status%msg
        call check(status_ok(status), message)
    end subroutine
    subroutine check(ok, message)
        logical, intent(in)::ok
        character(*), intent(in)::message
        if (.not. ok) then
            print *, message
            error stop 1
        end if
    end subroutine
end program
