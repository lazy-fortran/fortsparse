program test_ilutp_allocations
    use, intrinsic :: iso_c_binding, only: c_int
    use fortsparse, only: dp, sparse_ilutp_t, sparse_ilutp_stats_t, &
                          fortsparse_status_t, FORTSPARSE_OUT_OF_MEMORY
    implicit none
    interface
        subroutine arm(after) bind(C, name='fs_ilutp_alloc_arm')
            import :: c_int
            integer(c_int), value :: after
        end subroutine
        integer(c_int) function stop_probe() bind(C, &
                                                  name='fs_ilutp_alloc_stop')
            import :: c_int
        end function
    end interface
    type(sparse_ilutp_t) :: factor
    type(sparse_ilutp_stats_t) :: info
    type(fortsparse_status_t) :: status
    integer(c_int) :: calls, step
    real(dp) :: x(2)
    do step = 0, 12
        call arm(step)
        call setup()
        calls = stop_probe()
        if (status%code /= FORTSPARSE_OUT_OF_MEMORY) &
            error stop 'injected OOM not returned'
        info = factor%stats()
        if (info%ready .or. info%current_bytes /= 0) &
            error stop 'partial OOM arrays retained'
        call setup()
        if (status%code /= 0) &
            error stop 'setup did not recover after injected OOM'
    end do
    call arm(-1_c_int)
    do step = 1, 5
        call factor%apply([1.0_dp, 2.0_dp], x, status)
        if (status%code /= 0) error stop 'retained apply failed'
    end do
    calls = stop_probe()
    if (calls /= 0) error stop 'ILUTP apply allocated'
    if (maxval(abs(x - [.25_dp, .5_dp])) > epsilon(1.0_dp)) &
        error stop 'apply oracle'
    call factor%free()
    print '(a)', 'PASS setup allocation recovery and allocation-free apply'
contains
    subroutine setup()
        call factor%setup(2, [1, 2, 4], [2, 1, 2], &
                          [2.0_dp, 2.0_dp, 3.0_dp], 1, 1, 0.0_dp, &
                          .1_dp, status)
    end subroutine
end program
