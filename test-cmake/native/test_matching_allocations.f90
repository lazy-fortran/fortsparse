program test_matching_allocations
    use fortsparse
    use fortsparse_native_general, only: native_general_backend_t
    use, intrinsic :: iso_c_binding, only: c_int64_t, c_int
    implicit none
    type(native_general_backend_t)::backend
    type(csc_t)::a
    type(csc_z_t)::z
    type(fortsparse_status_t)::status
    type(sparse_factor_stats_t)::stats
    integer, parameter::n = 5
    integer::i, j, k, ord, rr(n*n), cc(n*n)
    real(dp)::vv(n*n), dense(n, n), exact(n), rhs(n), x(n)
    complex(dp) :: zrhs(n),zx(n)
    integer(c_int64_t)::count, current, maximum, nth, total, cap
    interface
        integer(c_int) function bounded_available() &
            bind(C, name='fsparse_native_amd_bounded_available')
            import c_int
        end function
        subroutine begin(nth) bind(C, name='fsparse_spd_probe_begin')
            import c_int64_t
            integer(c_int64_t), value::nth
        end subroutine
        subroutine finish(count, current, maximum) bind(C, &
                                                        name='fsparse_spd_probe_end')
            import c_int64_t
            integer(c_int64_t), intent(out)::count, current, maximum
        end subroutine
    end interface
    dense = 0
    do i = 1, n
        dense(i, mod(i, n) + 1) = 4 + i
        dense(i, mod(i + 1, n) + 1) = 0.25_dp
    end do
    k = 0
    do j = 1, n
        do i = 1, n
            k = k + 1; rr(k) = i; cc(k) = j; vv(k) = dense(i, j)
        end do
    end do
    call csc_from_triplet(n, n, rr, cc, vv, a, status)
    call csc_from_triplet(n, n, rr, cc, cmplx(vv, 0.5_dp*vv, dp), z, status)
    exact = [(sin(real(i, dp)), i=1, n)]; rhs = matmul(dense, exact)
    zrhs=cmplx(rhs,0.5_dp*rhs,dp)
    backend%options%structural_matching = .true.
    backend%options%pivot_tolerance = 0.1_dp
    backend%options%memory_policy = FORTSPARSE_MEMORY_LOW
    do ord = 0, 1
        if (ord == 1 .and. bounded_available() /= 1) cycle
        backend%options%ordering = merge(FORTSPARSE_ORDER_AMD, &
                                         FORTSPARSE_ORDER_DEFAULT, ord == 1)
        call begin(0_c_int64_t)
        call backend%factor_real(a, .false., status)
        call check(status_ok(status), 'traced matching factor')
        call backend%factor_stats(stats, status)
        call finish(count, current, maximum)
        call check(current == stats%host_current_bytes, &
                   'external retained array payload equality')
        call check(maximum == stats%host_peak_bytes, &
                   'external peak including matching overlap equality')
        total = count; cap = maximum
        call backend%free()
        backend%options%host_memory_limit_bytes = cap
        call begin(0_c_int64_t)
        call backend%factor_real(a, .false., status)
        call check(status_ok(status), 'external exact peak cap factor')
        call backend%free()
        call finish(count, current, maximum)
        call check(current == 0, 'all tracked arrays released')
        call check(maximum <= cap, 'actual allocation trace respects cap')
        backend%options%host_memory_limit_bytes = 28_c_int64_t*n - 1
        call begin(0_c_int64_t)
        call backend%factor_real(a, .false., status)
        call check(status%code == FORTSPARSE_RESOURCE_LIMIT, &
                   'matching refused before allocation')
        call finish(count, current, maximum)
        call check(count == 0 .and. current == 0, &
                   'no allocation on initial matching cap refusal')
        backend%options%host_memory_limit_bytes = 0
        ! Fail every reached allocator call, including partial multi-allocation,
        ! then independently check cleanup and successful reuse of the same object.
        do nth = 1, total
            call begin(nth)
            call backend%factor_real(a, .false., status)
            call check(status%code == FORTSPARSE_OUT_OF_MEMORY, &
                       'injected matching/factor allocation failure')
            call backend%free()
            call finish(count, current, maximum)
            call check(current == 0, &
                       'injected failure releases every owned allocation')
            call backend%factor_real(a, .false., status)
            call check(status_ok(status), 'recovery after allocation failure')
            call backend%solve_real(rhs, x, status)
            call check(status_ok(status) .and. maxval(abs(x - exact)) < 1e-12_dp, &
                       'recovered original-matrix solve')
            call backend%free()
        end do
        ! Complex original-matrix refinement arrays and matched AMD workspace
        ! are covered by the same external lifetime/payload oracle.
        call begin(0_c_int64_t)
        call backend%factor_complex(z, .true., status)
        call check(status_ok(status), &
                   'complex matching with original checker storage')
        call backend%factor_stats(stats, status)
        call finish(count, current, maximum)
        total = count
        call check(current == stats%host_current_bytes, &
                   'complex retained payload equality')
        call check(maximum == stats%host_peak_bytes, &
                   'complex actual peak equality')
        ! Existing refinement admission reserves its future checker quota while
        ! factoring, even when the checker arrays will be allocated afterward.
        cap = maximum + 4_c_int64_t*(n + 1 + n*n) + 16_c_int64_t*(n*n + 2*n)
        call backend%free()
        backend%options%host_memory_limit_bytes = cap
        call begin(0_c_int64_t)
        call backend%factor_complex(z, .true., status)
        call check(status_ok(status), &
                   'complex cap includes conservative future checker quota')
        call backend%free()
        call finish(count, current, maximum)
        call check(current == 0 .and. maximum <= cap, 'complex cap and lifetime')
        backend%options%host_memory_limit_bytes = 0
        do nth = 1, total
            call begin(nth)
            call backend%factor_complex(z, .true., status)
            call check(status%code == FORTSPARSE_OUT_OF_MEMORY, &
                       'complex injected allocation failure')
            call backend%free()
            call finish(count, current, maximum)
            call check(current == 0, 'complex failure releases all allocations')
            call backend%factor_complex(z, .true., status)
            call check(status_ok(status), 'complex allocation failure recovery')
            call backend%solve_complex(zrhs, zx, status)
            call check(status_ok(status) .and. maxval(abs(zx - exact)) < 1e-12_dp, &
                       'complex recovered solve')
            call backend%free()
        end do
    end do
contains
    subroutine check(ok, message)
        logical, intent(in)::ok
        character(*), intent(in)::message
        if (.not. ok) then
            print *, message, count, current, maximum, stats%host_current_bytes, &
                stats%host_peak_bytes
            print *, status%code, status%msg
            error stop 1
        end if
    end subroutine
end program
