module structured_raw_test_support
    use fortsparse
    use fortsparse_backend, only: sparse_backend_t
    use fortsparse_native, only: native_backend_t
    use fortsparse_native_general, only: native_general_backend_t
    use, intrinsic :: iso_c_binding, only: c_int64_t
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
    call exercise(spd, 1)
    do family = 2, 5
        call exercise(general, family)
    end do
end program
