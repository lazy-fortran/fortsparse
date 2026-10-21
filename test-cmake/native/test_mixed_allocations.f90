program test_mixed_allocations
    use fortsparse
    use fortnum_precision, only: FP32
    use iso_c_binding, only: c_int
    implicit none
    interface
        subroutine arm(after) bind(C,name='fs_ilutp_alloc_arm')
            import c_int
            integer(c_int), value :: after
        end subroutine
        integer(c_int) function stop_probe() bind(C,name='fs_ilutp_alloc_stop')
            import c_int
        end function
    end interface
    type(csc_t) :: a
    type(csc_z_t) :: z
    type(mixed_lu_t) :: factor
    type(mixed_lu_z_t) :: zfactor
    type(mixed_solve_report_t) :: report
    type(fortsparse_status_t) :: status
    real(dp) :: b(3),x(3)
    complex(dp) :: zb(3),zx(3)
    integer :: k,calls
    call csc_from_triplet(3,3,[1,2,2,3],[1,1,2,3], &
        [1.1_dp,0.2_dp,2.3_dp,3.7_dp],a,status)
    call csc_from_triplet(3,3,[1,2,2,3],[1,1,2,3], &
        cmplx(a%val,0.13_dp,dp),z,status)
    call factor%factor(a,status,FP32)
    if (status%code/=FORTSPARSE_OK) error stop 'real factor'
    call zfactor%factor(z,status,FP32)
    if (status%code/=FORTSPARSE_OK) error stop 'complex factor'
    b=[1.0_dp,2.0_dp,3.0_dp]
    zb=cmplx(b,0.5_dp,dp)
    call arm(-1_c_int)
    do k=1,18
        call factor%solve(b,x,report,status,mod(k,3))
        if (status%code/=FORTSPARSE_OK) error stop 'real solve'
        call zfactor%solve(zb,zx,report,status,mod(k,3))
        if (status%code/=FORTSPARSE_OK) error stop 'complex solve'
        x=b
        call factor%apply(x,status)
        zx=zb
        call zfactor%apply(zx,status)
    end do
    calls=stop_probe()
    if (calls/=0) then
        print *, 'unexpected retained-solve allocations: ',calls
        error stop 1
    end if
    print *, 'real/complex mixed N/T/H and preconditioner apply allocate nothing PASS'
end program
