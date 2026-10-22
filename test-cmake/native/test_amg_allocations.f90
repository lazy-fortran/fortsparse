program test_amg_allocations
    use, intrinsic :: iso_c_binding, only: c_int
    use fortsparse, only: dp, csr_t, amg_params_t, amg_handle_t, amg_setup, &
                          amg_apply, amg_destroy, fortsparse_status_t
    use amg_test_support, only: poisson_grid, check
    implicit none
    interface
        subroutine arm(after) bind(C, name='fs_ilutp_alloc_arm')
            import c_int
            integer(c_int), value :: after
        end subroutine
        integer(c_int) function stop_probe() bind(C, name='fs_ilutp_alloc_stop')
            import c_int
        end function
    end interface
    type(csr_t) :: a
    type(amg_params_t) :: params
    type(amg_handle_t) :: handle
    type(fortsparse_status_t) :: info
    real(dp) :: b(1024), z(1024)
    integer :: repeat, cycle
    integer(c_int) :: calls
    call poisson_grid(32, 2, .false., 1.0_dp, a)
    b = 1
    do cycle = 1, 2
        params%cycle = cycle
        call amg_setup(a, params, handle, info)
        call check(info%code == 0, 'allocation probe setup')
        call amg_apply(handle, b, z, info) ! warm the OpenMP runtime
        call check(info%code == 0, 'allocation probe warmup')
        call arm(-1_c_int)
        do repeat = 1, 5
            call amg_apply(handle, b, z, info)
            if (info%code /= 0) error stop 'allocation probe apply'
        end do
        calls = stop_probe()
        call check(calls == 0, 'AMG apply allocated')
        call check(dot_product(b, z) > 0, 'allocation probe positive correction')
        call amg_destroy(handle)
    end do
end program
