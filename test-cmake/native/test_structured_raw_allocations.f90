module structured_raw_test_support
    use fortsparse
    use fortsparse_backend, only: sparse_backend_t
    use fortsparse_native, only: native_backend_t
    use fortsparse_native_general, only: native_general_backend_t
    use fortsparse_native_amd, only: bounded_amd_available
    use, intrinsic :: iso_c_binding, only: c_int64_t, c_int
    implicit none
    interface
        subroutine begin(nth) bind(C, name='fsparse_spd_probe_begin')
            import c_int64_t
            integer(c_int64_t), value :: nth
        end subroutine
        subroutine finish(count, current, maximum) bind(C, name='fsparse_spd_probe_end')
            import c_int64_t
            integer(c_int64_t), intent(out) :: count, current, maximum
        end subroutine
        integer(c_int) function amd_calls() bind(C, &
            name='fsparse_spd_probe_amd_calls')
            import c_int
        end function
        integer(c_int) function opaque_calls() bind(C, &
            name='fsparse_spd_probe_opaque_calls')
            import c_int
        end function
    end interface
contains
    subroutine exercise(backend, family)
        class(sparse_backend_t), intent(inout) :: backend
        integer, intent(in) :: family
        integer :: ptr(3), rows(4), policy, op, refining
        real(dp) :: a(2, 2), values(4), exact(2), b(2), x(2)
        real(dp), pointer :: vector(:)
        complex(dp) :: z(2, 2), zvalues(4), ze(2), zb(2), zx(2)
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: stats
        integer(c_int64_t) :: count, current, maximum, total, nth, cap
        logical :: reused, refine
        ptr = [1, 3, 5]; rows = [1, 2, 1, 2]
        a = reshape([4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2, 2])
        z = cmplx(a, 0.0_dp, dp)
        if (family == 2) a = reshape([0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp], [2, 2])
        if (family >= 3) then
            z(1, 2) = cmplx(1.0_dp, 1.0_dp, dp)
            z(2, 1) = conjg(z(1, 2))
        end if
        if (family == 4) then
            z(1, 1) = 0; z(2, 2) = 0
        else if (family == 5) then
            z(2, 1) = z(1, 2); z(1, 1) = cmplx(0.0_dp, 1.0_dp, dp)
            z(2, 2) = cmplx(-2.0_dp, 0.5_dp, dp)
        end if
        values = reshape(a, [4]); zvalues = reshape(z, [4])
        exact = [2.0_dp, -3.0_dp]; ze = cmplx(exact, [1.0_dp, 2.0_dp], dp)
        backend%options%symmetry = FORTSPARSE_SYM_SPD
        if (family == 2 .or. &
            family == 5) backend%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        if (family == 4) backend%options%symmetry = FORTSPARSE_SYM_HERMITIAN
        if (family == 1) then
            refine = .false.
            backend%options%memory_policy = FORTSPARSE_MEMORY_AUTO
            backend%options%host_memory_limit_bytes = 0_i8
            call begin(0_c_int64_t)
            call factor()
            call ok(status_ok(status), 'uncapped SPD factor')
            call backend%factor_stats(stats, status)
            call finish(count, current, maximum)
            call ok(current == stats%host_current_bytes, &
                    'uncapped SPD retained payload')
            call backend%free()
        end if
        do refining = 0, 1
            refine = refining == 1
            do policy = 0, 1
                ! AUTO is bounded by a generous cap; LOW also works unlimited.
                backend%options%memory_policy = policy
                backend%options%host_memory_limit_bytes = 1000000_i8
                if (policy == 1) backend%options%host_memory_limit_bytes = 0_i8
                call begin(0_c_int64_t)
                call factor()
                call ok(status_ok(status), 'factor')
                call backend%factor_stats(stats, status)
                call finish(count, current, maximum)
                if (current /= stats%host_current_bytes .or. &
                    maximum /= stats%host_peak_bytes) then
                    print *, family, policy, refining, current, &
                        stats%host_current_bytes, maximum, stats%host_peak_bytes
                end if
                call ok(current == stats%host_current_bytes, &
                        'retained payload equals allocator')
                call ok(maximum == stats%host_peak_bytes, &
                        'simultaneous peak equals allocator')
                call ok(stats%strict_memory_supported, 'strict supported')
                total = count; cap = maximum
                call solve_all()
                call backend%free()
                ! Same fixed schedule must refuse a byte below observed peak.
                backend%options%host_memory_limit_bytes = cap - 1
                call begin(0_c_int64_t)
                call factor()
                call ok(status%code == FORTSPARSE_RESOURCE_LIMIT, &
                        'one-byte-below peak refuses')
                call backend%factor_stats(stats, status)
                call ok(stats%budget_required_bytes > stats%budget_available_bytes, &
                        'failure diagnostics')
                call backend%free()
                call finish(count, current, maximum)
                call ok(current == 0 .and. maximum <= cap - 1, &
                        'failure releases and respects cap')
                backend%options%host_memory_limit_bytes = 1000000_i8
                ! One trace spans old factors + conversion + retained checker.
                call begin(0_c_int64_t)
                call factor()
                if (family <= 2) then
                    call backend%refactor_real_raw(2, 2, 4, ptr, rows, values, &
                                                   refine, status, reused)
                else
                    call backend%refactor_complex_raw(2, 2, 4, ptr, rows, &
                                                      zvalues, refine, status, reused)
                end if
                if (.not. status_ok(status) .or. .not. reused) print *, family, &
                    policy, refining, status%code, reused, trim(status%msg)
                call ok(status_ok(status) .and. reused, 'same pattern reuse')
                call backend%factor_stats(stats, status)
                call finish(count, current, maximum)
                call ok(current == stats%host_current_bytes, &
                        'refactor retained allocator')
                call ok(maximum == stats%host_peak_bytes, &
                        'refactor conversion overlap peak')
                call solve_all()
                call backend%free()
                ! A retained factor must share the initial conversion admission.
                call factor()
                call backend%factor_stats(stats, status)
                backend%options%host_memory_limit_bytes = stats%host_current_bytes &
                                                          + 28_i8 + merge(32_i8, &
                                                              64_i8, family <= 2) - 1_i8
                call begin(0_c_int64_t)
                if (family <= 2) then
                    call backend%refactor_real_raw(2, 2, 4, ptr, rows, values, &
                                                   refine, status, reused)
                else
                    call backend%refactor_complex_raw(2, 2, 4, ptr, rows, &
                                                      zvalues, refine, status, reused)
                end if
                call ok(status%code == FORTSPARSE_RESOURCE_LIMIT, &
                        'old factor plus copy admission')
                call finish(count, current, maximum)
                call ok(count == 0, &
                        'old factor overlap rejected before conversion allocation')
                call backend%factor_stats(stats, status)
                call ok(stats%host_current_bytes == 0, &
                        'failed direct refactor invalidated and freed')
                ! After conversion release, vector allocation uses the full cap.
                backend%options%host_memory_limit_bytes = 1000000_i8
                call factor()
                call backend%factor_stats(stats, status)
                backend%options%host_memory_limit_bytes = stats%host_current_bytes + 80_i8
                vector => backend%vector(10)
                call ok(associated(vector), &
                        'released RAW budget available to user vector')
                call backend%free()
                nullify (vector)
                backend%options%host_memory_limit_bytes = 1000000_i8
                do nth = 1, total
                    call begin(nth)
                    call factor()
                    call ok(status%code == FORTSPARSE_OUT_OF_MEMORY, &
                            'allocation fault status')
                    call backend%free()
                    call finish(count, current, maximum)
                    call ok(current == 0, 'allocation fault cleanup')
                    call factor()
                    call ok(status_ok(status), 'allocation fault recovery')
                    call backend%free()
                end do
            end do
        end do
    contains
        subroutine factor()
            if (family <= 2) then
                call backend%factor_real_raw(2, 2, 4, ptr, rows, values, refine, status)
            else
                call backend%factor_complex_raw(2, 2, 4, ptr, rows, zvalues, &
                                                refine, status)
            end if
        end subroutine
        subroutine solve_all()
            do op = 0, 2
                if (family <= 2) then
                    b = matmul(a, exact)
                    call backend%solve_real_op(b, x, op, status)
                    call ok(status_ok(status) .and. maxval(abs(x - exact)) < 1e-12_dp, &
                        'real N/T/H')
                else
                    select case (op)
                    case (0)
                        zb = matmul(z, ze)
                    case (1)
                        zb = matmul(transpose(z), ze)
                    case (2)
                        zb = matmul(conjg(transpose(z)), ze)
                    end select
                    call backend%solve_complex_op(zb, zx, op, status)
                    call ok(status_ok(status) .and. maxval(abs(zx - ze)) < 1e-12_dp, &
                        'complex N/T/H')
                end if
            end do
        end subroutine
    end subroutine
    subroutine exercise_raw_amd()
        type(native_backend_t) :: backend
        type(sparse_solver_t) :: request
        type(csc_t) :: a, changed
        type(sparse_memory_plan_t) :: plan, retained_plan
        type(sparse_factor_stats_t) :: stats, unchanged
        type(fortsparse_status_t) :: status
        integer, parameter :: n = 9, nz = 3*n - 2
        integer :: rr(nz), cc(nz), i, k, refinement, before, vectors
        real(dp) :: dense(n, n), vv(nz), exact(n), rhs(n), x(n)
        real(dp) :: block(n, 2), answer(n, 2)
        real(dp), pointer :: vector(:)
        integer(c_int64_t) :: count, live, peak
        integer(i8) :: conversion, cap, retained
        logical :: reused

        if (bounded_amd_available() /= 1) return
        dense = 0.0_dp
        do i = 1, n
            dense(i, i) = 4.0_dp + 0.125_dp*i
            if (i > 1) then
                dense(i, 1) = -0.25_dp
                dense(1, i) = -0.25_dp
            end if
            exact(i) = sin(real(i, dp))
        end do
        k = 0
        do i = 1, n
            k = k + 1; rr(k) = i; cc(k) = i; vv(k) = dense(i, i)
            if (i == 1) cycle
            k = k + 1; rr(k) = i; cc(k) = 1; vv(k) = dense(i, 1)
            k = k + 1; rr(k) = 1; cc(k) = i; vv(k) = dense(1, i)
        end do
        call csc_from_triplet(n, n, rr, cc, vv, a, status)
        call ok(status_ok(status), 'AMD seam input')
        changed = a; changed%val = 2*a%val
        rhs = matmul(dense, exact)
        block(:, 1) = rhs; block(:, 2) = -2*rhs
        conversion = 4_i8*(n + 1 + nz) + 8_i8*nz
        request%backend_id = FORTSPARSE_BACKEND_NATIVE
        request%options%symmetry = FORTSPARSE_SYM_SPD
        request%options%ordering = FORTSPARSE_ORDER_AMD
        request%options%memory_policy = FORTSPARSE_MEMORY_LOW
        do refinement = 0, 1
            request%refine = refinement == 1
            request%options%host_memory_limit_bytes = 0_i8
            call sparse_analyze(request, a, plan, status, nrhs=2)
            call ok(status_ok(status) .and. plan%fits, 'AMD CSC plan')
            call ok(plan%prediction_quality == MEMORY_PLAN_EXACT_LAYOUT, &
                    'fresh AMD CSC exact layout')
            backend%options = request%options
            before = amd_calls()
            call begin(0_c_int64_t)
            call factor()
            call ok(status_ok(status), 'AMD LOW RAW factor')
            call backend%factor_stats(stats, status)
            call finish(count, live, peak)
            call ok(amd_calls() == before + 1 .and. opaque_calls() == 0, &
                    'RAW uses audited AMD wrapper')
            call ok(any(backend%perm /= [(i, i=1, n)]), 'nontrivial RAW permutation')
            call ok(live == stats%host_current_bytes .and. &
                peak == stats%host_peak_bytes, &
                    'AMD RAW allocator equality')
            call ok(peak == plan%peak_host_bytes + conversion, &
                    'exact CSC plan plus RAW temporary quota')
            call solve_original()
            cap = peak
            call backend%free()
            backend%options%host_memory_limit_bytes = cap - 1_i8
            call begin(0_c_int64_t)
            call factor()
            call ok(status%code == FORTSPARSE_RESOURCE_LIMIT, 'AMD RAW byte below peak')
            call backend%free()
            call finish(count, live, peak)
            call ok(live == 0 .and. peak <= cap - 1_i8, &
                'AMD RAW refusal cap and cleanup')
            backend%options%host_memory_limit_bytes = cap
            call begin(0_c_int64_t)
            call factor()
            call ok(status_ok(status), 'AMD RAW exact peak cap')
            call backend%factor_stats(stats, status)
            retained = stats%host_current_bytes
            before = amd_calls()
            call backend%refactor_real_raw(n, n, nz, a%col_ptr, a%row_idx, &
                                           changed%val, request%refine, status, reused)
            call ok(status_ok(status) .and. reused, 'AMD RAW retained numeric reuse')
            call ok(amd_calls() == before, 'retained RAW reuses AMD ordering')
            call backend%factor_stats(stats, status)
            call finish(count, live, peak)
            call ok(live == stats%host_current_bytes .and. &
                peak == stats%host_peak_bytes, &
                    'AMD RAW retained conversion overlap equality')
            call ok(peak == max(cap, retained + conversion), 'retained RAW exact peak')
            call backend%solve_real_block(2*block, answer, status)
            call ok(status_ok(status), 'retained RAW block status')
            call ok(maxval(abs(matmul(dense, answer) - block)) < 1e-12_dp, &
                    'retained RAW block original-matrix oracle')
            call ok(backend%options%host_memory_limit_bytes == cap, &
                    'RAW restores original cap after refactor')
            vectors = int((cap - live)/8_i8)
            call ok(vectors > 0, 'released conversion provides vector headroom')
            vector => backend%vector(vectors)
            call ok(associated(vector), 'vector uses full restored RAW cap')
            call backend%free()
            nullify (vector)
            call factor()
            call ok(status_ok(status), 'retained AMD RAW admission setup')
            call backend%factor_stats(stats, status)
            backend%options%host_memory_limit_bytes = &
                stats%host_current_bytes + conversion - 1_i8
            call begin(0_c_int64_t)
            call backend%refactor_real_raw(n, n, nz, a%col_ptr, a%row_idx, &
                                           changed%val, request%refine, status, reused)
            call ok(status%code == FORTSPARSE_RESOURCE_LIMIT .and. .not. reused, &
                    'retained AMD RAW conversion quota refusal')
            call finish(count, live, peak)
            call ok(count == 0 .and. live == 0, &
                'retained AMD refuses before copy allocation')
            ! Public planning preserves the retained RAW factor and its history;
            ! its BOUND quality explicitly accounts a second factor's overlap.
            call sparse_factor(request, n, n, nz, a%col_ptr, a%row_idx, a%val, status)
            call ok(status_ok(status), 'public AMD LOW RAW factor')
            call sparse_factor_stats(request, stats, status)
            request%options%host_memory_limit_bytes = &
                stats%host_current_bytes + plan%peak_host_bytes - 1_i8
            call sparse_analyze(request, a, retained_plan, status, nrhs=2)
            call ok(status%code == FORTSPARSE_RESOURCE_LIMIT .and. &
                .not. retained_plan%fits, &
                    'retained RAW plus planned factor byte-below refusal')
            call ok(request%factored, 'refused planning preserves retained RAW factor')
            request%options%host_memory_limit_bytes = &
                stats%host_current_bytes + plan%peak_host_bytes
            call sparse_analyze(request, a, retained_plan, status, nrhs=2)
            call ok(status_ok(status) .and. retained_plan%fits, &
                    'planning beside retained RAW factor')
            call ok(retained_plan%prediction_quality == MEMORY_PLAN_BOUND .and. &
                    retained_plan%peak_host_bytes == plan%peak_host_bytes, &
                    'retained RAW plan quality and layout')
            call sparse_factor_stats(request, unchanged, status)
            call ok(unchanged%host_current_bytes == stats%host_current_bytes .and. &
                    unchanged%host_peak_bytes == stats%host_peak_bytes, &
                    'planning preserves retained RAW allocation history')
            call sparse_solve(request, rhs, x, status)
            call ok(status_ok(status) .and. maxval(abs(x - exact)) < 1e-12_dp, &
                    'retained RAW survives public planning')
            call sparse_destroy(request)
        end do
    contains
        subroutine factor()
            call backend%factor_real_raw(n, n, nz, a%col_ptr, a%row_idx, a%val, &
                                         request%refine, status)
        end subroutine
        subroutine solve_original()
            integer :: operation
            do operation = 0, 2
                call backend%solve_real_op(rhs, x, operation, status)
                call ok(status_ok(status) .and. maxval(abs(x - exact)) < 1e-12_dp, &
                        'AMD RAW original scalar N/T/H oracle')
            end do
            call backend%solve_real_block(block, answer, status)
            call ok(status_ok(status), 'AMD RAW original block status')
            call ok(maxval(abs(matmul(dense, answer) - block)) < 1e-12_dp, &
                    'AMD RAW original block oracle')
        end subroutine
    end subroutine exercise_raw_amd

    subroutine ok(condition, message)
        logical, intent(in)::condition
        character(*), intent(in)::message
        if (.not. condition) then
            print *, message
            error stop 'structured RAW allocation oracle'
        end if
    end subroutine
end module
program test_structured_raw_allocations
    use structured_raw_test_support
    implicit none
    type(native_backend_t) :: spd
    type(native_general_backend_t) :: general
    integer :: family
    call exercise_raw_amd()
    call exercise(spd, 1)
    do family = 2, 5
        call exercise(general, family)
    end do
end program
