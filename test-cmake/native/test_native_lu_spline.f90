program test_native_lu_spline
    use fortsparse
    implicit none
    type(csc_t) :: a
    type(sparse_solver_t) :: solver
    type(fortsparse_status_t) :: status
    integer, allocatable :: ir(:), jc(:)
    real(dp), allocatable :: av(:), b(:), x(:)
    real(dp) :: nodes(8), expected(8), evaluated(8), error, h
    integer :: n, m, nz, i, u, nrhs
    character(512) :: directory
    character(120) :: header

    call get_command_argument(1,directory)
    open(newunit=u,file=trim(directory)//'/spline54.mtx',status='old')
    read(u,'(a)')header
    read(u,*)n,m,nz
    allocate(ir(nz),jc(nz),av(nz),b(n),x(n))
    do i=1,nz
        read(u,*)ir(i),jc(i),av(i)
    end do
    close(u)
    call csc_from_triplet(n,m,ir,jc,av,a,status)
    call check_status(status)
    open(newunit=u,file=trim(directory)//'/spline54-rhs.mtx',status='old')
    read(u,'(a)')header
    read(u,*)n,nrhs
    do i=1,n
        read(u,*)b(i)
    end do
    close(u)
    nodes=[0.0_dp,0.08_dp,0.21_dp,0.37_dp,0.52_dp,0.70_dp,0.86_dp,1.0_dp]
    expected=sin(2*nodes)+0.25_dp*nodes*nodes
    solver%backend_id=FORTSPARSE_BACKEND_NATIVE
    solver%refine=.false.
    call sparse_factor(solver,a,status)
    call check_status(status)
    call sparse_solve(solver,b,x,status)
    call check_status(status)
    evaluated(:7)=x(1:43:7)
    h=nodes(8)-nodes(7)
    evaluated(8)=x(43)+h*(x(44)+h*(x(45)+h*x(46)))
    error=maxval(abs(evaluated-expected))
    if(error>1.0e-12_dp) then
        print *, 'FAIL: original spline node observable',error
        error stop 1
    end if
    call sparse_free(solver)
contains
    subroutine check_status(status)
        type(fortsparse_status_t), intent(in) :: status
        if(.not.status_ok(status)) then
            print *,status%code,status%msg
            error stop 2
        end if
    end subroutine
end program
