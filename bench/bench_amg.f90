program bench_amg
    use fortsparse, only: dp, i8, csr_t, csr_matvec, amg_params_t, amg_handle_t, &
                          amg_setup, amg_destroy, sparse_ilutp_t, fortsparse_status_t
    use fortsparse_native_ldlt, only: native_ldlt_t
    use fortsparse_csc, only: csc_t
    use amg_test_support, only: poisson_grid, pcg, check
    use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options
    use, intrinsic :: iso_c_binding, only: c_double
!$  use omp_lib, only: omp_get_num_threads
    implicit none
    interface
        subroutine metrics(rss, childrss, childcpu) &
            bind(C, name='fsbench_process_metrics')
            import c_double
            real(c_double), intent(out) :: rss, childrss, childcpu
        end subroutine
        real(c_double) function process_cpu() bind(C, name='fsbench_parent_cpu_seconds')
            import c_double
        end function
    end interface
    type(csr_t) :: a
    type(csc_t) :: csc
    type(amg_params_t) :: params
    type(amg_handle_t) :: amg
    type(sparse_ilutp_t) :: ilu
    type(native_ldlt_t) :: direct
    type(fortsparse_status_t) :: info
    real(dp), allocatable :: b(:), x(:), exact(:), work(:)
    real(dp) :: eps, setup_time, solve_time, residual, error, start
    real(c_double) :: rss, childrss, childcpu
    integer :: side, dim, family, n, i, iterations, threads, levels, code
    integer(i8) :: ticks, rate
    character(64) :: arg, solver, revision, host, job
    if (command_argument_count() /= 7) then
        print *, 'bench_amg SIDE DIM FAMILY(1=P1,2=Q1) EPS THETA SOLVER MAXITER'
        stop 2
    end if
    call get_command_argument(1, arg)
    read (arg, *) side
    call get_command_argument(2, arg)
    read (arg, *) dim
    call get_command_argument(3, arg)
    read (arg, *) family
    call get_command_argument(4, arg)
    read (arg, *) eps
    call get_command_argument(5, arg)
    read (arg, *) params%strength_threshold
    call get_command_argument(6, solver)
    call get_command_argument(7, arg)
    read (arg, *) code
    call get_environment_variable('AMG_REVISION', revision)
    call get_environment_variable('HOSTNAME', host)
    call get_environment_variable('SLURM_JOB_ID', job)
    threads = 1
    !$omp parallel
    !$omp single
!$  threads = omp_get_num_threads()
    !$omp end single
    !$omp end parallel
    call poisson_grid(side, dim, family == 2, eps, a)
    n = a%nrow
    allocate (b(n), x(n), exact(n), work(n))
    !$omp parallel do
    do i = 1, n
        exact(i) = sin(real(i, dp)*.271_dp) + cos(real(i, dp)*.173_dp)
    end do
    call csr_matvec(a, exact, b)
    start = wall()
    select case (trim(solver))
    case ('amg')
        call amg_setup(a, params, amg, info)
    case ('ilutp')
        call ilu%setup(n, a%row_ptr, a%col, a%val, 12, 12, &
                       1.e-3_dp, 0.0_dp, info)
    case ('ldlt')
        csc%nrow = n
        csc%ncol = n
        csc%nnz = size(a%val)
        csc%col_ptr = a%row_ptr
        csc%row_idx = a%col
        csc%val = a%val
        call direct%factor(csc, info)
    case default
        error stop 'unknown solver'
    end select
    setup_time = wall() - start
    call check(info%code == 0, 'setup failed: '//trim(info%msg))
    start = wall()
    select case (trim(solver))
    case ('amg')
        call pcg(a, b, x, iterations, residual, amg=amg, maxiter=code)
    case ('ilutp')
        call pcg(a, b, x, iterations, residual, ilu=ilu, maxiter=code)
    case ('ldlt')
        call direct%solve(b, x, info)
        call check(info%code == 0, 'direct solve')
        call csr_matvec(a, x, work)
        residual = sqrt(sum((b - work)**2)/sum(b*b))
        iterations = 1
    end select
    solve_time = wall() - start
    error = maxval(abs(x - exact))
    levels = amg%levels()
    call metrics(rss, childrss, childcpu)
    print '(a)', '# '//compiler_version()//' '//compiler_options()
    print '(a)', '# revision='//trim(revision)//' host='//trim(host)//' job='//trim(job)
    print '(a)', 'solver,side,dim,family,eps,theta,rows,nnz,threads,levels,'// &
        'complexity,iterations,residual,max_error,setup_s,solve_s,'// &
        'process_cpu_to_sample_s,peak_rss_kib,converged'
    print '(*(g0,:,","))', trim(solver), side, dim, family, eps, &
        params%strength_threshold, n, size(a%val), threads, levels, &
        amg%operator_complexity(), iterations, residual, error, setup_time, &
        solve_time, process_cpu(), rss, residual < 1.e-11_dp
    call amg_destroy(amg)
    call ilu%free()
    call direct%free()
contains
    real(dp) function wall() result(t)
        call system_clock(ticks, rate)
        t = real(ticks, dp)/real(rate, dp)
    end function
end program
