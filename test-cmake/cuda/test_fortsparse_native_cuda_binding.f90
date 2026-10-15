program test_fortsparse_native_cuda_binding
    use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_size_t, c_int, &
        c_double, c_loc
    use fortsparse_native_cuda, only: native_cuda_available, spd_cuda_create, &
        spd_cuda_destroy, spd_cuda_upload, spd_cuda_download, spd_cuda_panel_factor, &
        spd_cuda_factor_info, spd_cuda_synchronize, spd_cuda_allocated_bytes
    implicit none
    type(c_ptr) :: context
    real(c_double), target :: panel(2,2)
    integer(c_int) :: status, info
    if (.not. native_cuda_available) error stop 'CUDA test needs CUDA build'
    context = c_null_ptr
    status = spd_cuda_create(context, 0_c_int, 4_c_size_t)
    call checked(status)
    if (spd_cuda_allocated_bytes(context) < 32_c_size_t) error stop 'CUDA allocation size'
    panel = reshape([4.0_c_double, 2.0_c_double, 2.0_c_double, 10.0_c_double], [2,2])
    call checked(spd_cuda_upload(context, 0_c_size_t, c_loc(panel), 4_c_size_t))
    call checked(spd_cuda_panel_factor(context, 0_c_size_t, 2_c_int, 2_c_int))
    call checked(spd_cuda_factor_info(context, info))
    if (info /= 0) error stop 'Fortran CUDA pivot status'
    call checked(spd_cuda_download(context, 0_c_size_t, c_loc(panel), 4_c_size_t))
    call checked(spd_cuda_synchronize(context))
    if (abs(panel(1,1)-2.0_c_double) > 1.0e-14_c_double) error stop 'CUDA factor diagonal'
    if (abs(panel(2,1)-1.0_c_double) > 1.0e-14_c_double) error stop 'CUDA factor offdiagonal'
    if (abs(panel(2,2)-3.0_c_double) > 1.0e-14_c_double) error stop 'CUDA factor diagonal'
    call checked(spd_cuda_destroy(context))
    print '(a)', 'Fortran ISO_C_BINDING retained CUDA context and panel factor pass'
contains
    subroutine checked(code)
        integer(c_int), intent(in) :: code
        if (code /= 0) error stop 'Fortran CUDA C ABI failed'
    end subroutine checked
end program test_fortsparse_native_cuda_binding
