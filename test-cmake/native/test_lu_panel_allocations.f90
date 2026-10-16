program test_lu_panel_allocations
    use, intrinsic :: iso_c_binding, only: c_int, c_int64_t
    use fortsparse_kinds, only: dp, i8
    use fortsparse_status
    use fortsparse_native_lu, only: lu_factor_t
    implicit none
    interface
        subroutine arm(after) bind(C, name='fs_panel_arm')
            import c_int
            integer(c_int), value :: after
        end subroutine
        integer(c_int) function stop() bind(C, name='fs_panel_stop')
            import c_int
        end function
        integer(c_int) function failed() bind(C, name='fs_panel_failed')
            import c_int
        end function
        integer(c_int64_t) function current() bind(C, name='fs_panel_current')
            import c_int64_t
        end function
        integer(c_int64_t) function peak() bind(C, name='fs_panel_peak')
            import c_int64_t
        end function
    end interface
    integer, parameter :: n = 17
    type(lu_factor_t) :: factor
    type(fortsparse_status_t) :: status
    integer :: ptr(n + 1), rows(n*n), i, j, k, allocs, ignored, point
    real(dp) :: values(n*n), rhs(n)
    integer(i8) :: maxlive, retained
    k = 1
    do j = 1, n
        ptr(j) = k
        do i = 1, n
            rows(k) = i
            values(k) = real(mod(i + 3*j, 7) - 3, dp)/64
            if (i == j) values(k) = 20
            k = k + 1
        end do
    end do
    ptr(n + 1) = k
    call arm(-1)
    call factor%factor_raw(n, ptr, rows, values, status, pruning=.true., panel_width=8)
    allocs = stop()
    call checked(status)
    if (current() /= factor%memory_bytes()) &
        error stop 'reached panel live payload differs from allocator'
    if (peak() /= factor%peak_memory_bytes) &
        error stop 'panel peak/old-new overlap differs from allocator'
    maxlive = peak(); retained = current()
    ! Same-pattern refactor reuses every owned buffer: no numerical allocation.
    call arm(-1)
    call factor%factor_raw(n, ptr, rows, values, status, pruning=.true., panel_width=8)
    ignored = stop()
    call checked(status)
    if (ignored /= 0 .or. current() /= retained) error stop 'panel refactor allocated'
    if (peak() /= factor%peak_memory_bytes) error stop 'refactor observed peak'
    rhs = 1
    call arm(0)
    call factor%solve(rhs, status)
    ignored = stop()
    call checked(status)
    if (ignored /= 0 .or. failed() /= 0) error stop 'panel solve allocated'
    call factor%free()
    if (current() /= 0) error stop 'panel payload leaked after free'
    ! Every reached allocation including all three panel setup arrays.
    do point = 0, allocs - 1
        call arm(point)
        call factor%factor_raw(n, ptr, rows, values, status, pruning=.true., &
                               panel_width=8)
        ignored = stop()
        if (failed() /= 1) error stop 'allocation failure point not reached'
        if (status%code /= FORTSPARSE_OUT_OF_MEMORY) &
            error stop 'panel allocation failure status'
        if (current() /= factor%memory_bytes()) &
            error stop 'failure retained payload differs'
        if (peak() /= factor%peak_memory_bytes) &
            error stop 'failure reached peak differs'
        rhs = 1
        call factor%solve(rhs, status)
        if (status%code /= FORTSPARSE_NOT_FACTORED) error stop 'partial factor ready'
        ! Recover directly from partially allocated retained setup/factors.
        call arm(-1)
        call factor%factor_raw(n, ptr, rows, values, status, pruning=.true., &
                               panel_width=8)
        ignored = stop(); call checked(status)
        call factor%free()
        if (current() /= 0) error stop 'allocation recovery leaked'
    end do
    ! Fresh cap boundary is independently observed. Growth may compact spare
    ! capacity to fit a smaller cap, so do not infer a one-byte threshold from
    ! an unlimited peak. Find the exact minimal successful cap for this fixture.
    block
        integer(i8) :: lower, upper, mid, minimum
        lower = 0; upper = maxlive
        do while (upper - lower > 1)
            mid = lower + (upper - lower)/2
            call arm(-1)
            call factor%factor_raw(n, ptr, rows, values, status, pruning=.true., &
                                   memory_limit_bytes=mid, panel_width=8)
            ignored = stop()
            if (current() > mid .or. peak() > mid) &
                error stop 'strict panel cap exceeded'
            if (status%code == FORTSPARSE_OK) then
                upper = mid
            else
                if (status%code /= FORTSPARSE_RESOURCE_LIMIT) &
                    error stop 'cap search status'
                lower = mid
            end if
            call factor%free()
        end do
        minimum = upper
        call arm(-1)
        call factor%factor_raw(n, ptr, rows, values, status, pruning=.true., &
                               memory_limit_bytes=minimum - 1, panel_width=8)
        ignored = stop()
        if (status%code /= FORTSPARSE_RESOURCE_LIMIT) &
            error stop 'one-byte-below panel cap'
        if (factor%failed_required_bytes <= factor%failed_available_bytes) &
            error stop 'cap diagnostic'
        call factor%free()
        call arm(-1)
        call factor%factor_raw(n, ptr, rows, values, status, pruning=.true., &
                               memory_limit_bytes=minimum, panel_width=8)
        ignored = stop(); call checked(status)
        if (peak() > minimum .or. current() /= factor%memory_bytes()) &
            error stop 'exact cap accounting'
        call factor%free()
    end block
    call arm(-1)
    call factor%factor_raw(n, ptr, rows, values, status, pruning=.true.)
    ignored = stop(); call checked(status)
    retained = current()
    call arm(-1)
    call factor%factor_raw(n, ptr, rows, values, status, pruning=.true., &
                           memory_limit_bytes=retained + 68_i8*n + 63_i8, panel_width=8)
    ignored = stop()
    if (status%code /= FORTSPARSE_RESOURCE_LIMIT .or. ignored /= 0) &
        error stop 'panel reuse cap must reject before allocation'
    if (current() /= retained .or. peak() /= retained) &
        error stop 'panel cap changed retained payload'
    if (factor%failed_required_bytes /= 68_i8*n + 64_i8) &
        error stop 'panel reuse required workspace diagnostic'
    call arm(-1)
    call factor%factor_raw(n, ptr, rows, values, status, pruning=.true., &
                           memory_limit_bytes=retained + 68_i8*n + 64_i8, panel_width=8)
    ignored = stop(); call checked(status)
    if (current() /= retained + 68_i8*n + 64_i8 .or. ignored /= 3) &
        error stop 'exact panel reuse admission/layout'
    if (current() /= factor%memory_bytes() .or. peak() /= factor%peak_memory_bytes) &
        error stop 'exact panel reuse allocator observation'
    call factor%free()
    if (current() /= 0) error stop 'final panel allocation leak'
    print '(a,i0,a)', 'Panel allocation payload, overlap, cap and ', allocs, &
        ' fault points pass'
contains
    subroutine checked(s)
        type(fortsparse_status_t), intent(in) :: s
        if (s%code /= FORTSPARSE_OK) then
            print *, s%code, trim(s%msg)
            error stop 'unexpected panel allocation status'
        end if
    end subroutine checked
end program test_lu_panel_allocations
