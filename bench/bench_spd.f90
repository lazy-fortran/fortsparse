program bench_spd
    use fortsparse, only: dp, csc_t, sparse_solver_t, fortsparse_status_t, &
        status_ok, sparse_factor, sparse_solve, sparse_solve_block, &
        sparse_refactor, sparse_free, sparse_factor_stats_t, &
        sparse_factor_stats, FORTSPARSE_SYM_SPD
    use spd_fixtures, only: fixture_t, make_fixture, fixture_product, &
        backward_error, wall_time
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    interface
        subroutine metrics(rss, child_rss, child_cpu) &
                bind(C, name="fsbench_process_metrics")
            import :: c_double
            real(c_double), intent(out) :: rss, child_rss, child_cpu
        end subroutine
    end interface
    type(fixture_t) :: fixture
    type(csc_t) :: a
    type(sparse_solver_t) :: solver
    type(fortsparse_status_t) :: status
    type(sparse_factor_stats_t) :: stats
    character(32) :: backend, family, argument, workload
    integer :: level, threads, nrhs, j
    logical :: held_out, reused
    real(dp), allocatable :: b(:, :), x(:, :), reference(:)
    real(dp) :: start, cpu_start, factor_s, solve_s, refactor_s, cpu_s
    real(dp) :: error, rss, child_rss, child_cpu, total_start, total_s

    total_start = wall_time()
    call cpu_time(cpu_start)
    call get_command_argument(1, backend)
    call get_command_argument(2, family)
    call get_command_argument(3, argument)
    read (argument, *) level
    call get_command_argument(4, argument)
    read (argument, *) threads
    call get_command_argument(5, workload)
    call get_command_argument(6, argument)
    held_out = trim(argument) == "heldout"
    nrhs = 1
    if (trim(workload) == "rhs8") nrhs = 8
    select case (trim(backend))
    case ("superlu")
        solver%backend_id = 1
    case ("umfpack")
        solver%backend_id = 2
    case ("mumps")
        solver%backend_id = 3
    case ("native")
        solver%backend_id = 4
    case default
        error stop "unknown backend"
    end select
    if (trim(workload) /= "rhs1" .and. trim(workload) /= "rhs8" &
        .and. trim(workload) /= "refactor") error stop "unknown workload"
    call make_fixture(trim(family), level, held_out, fixture)
    a%nrow = fixture%n
    a%ncol = fixture%n
    a%nnz = size(fixture%val)
    a%col_ptr = fixture%ptr
    a%row_idx = fixture%row
    a%val = fixture%val
    allocate (b(a%nrow, nrhs), x(a%nrow, nrhs), reference(a%nrow))
    do j = 1, nrhs
        reference = 1.0_dp + real(j, dp)/8.0_dp
        reference(1:a%nrow:2) = -0.5_dp
        call fixture_product(fixture, reference, b(:, j))
    end do
    solver%options%symmetry = FORTSPARSE_SYM_SPD
    solver%options%nthreads = threads
    start = wall_time()
    call sparse_factor(solver, a, status)
    factor_s = wall_time() - start
    if (.not. status_ok(status)) then
        write (*, '(a,i0)') "backend status: ", status%code
        stop 2
    end if
    call sparse_factor_stats(solver, stats, status)
    refactor_s = -1.0_dp
    reused = .false.
    if (trim(workload) == "refactor") then
        a%val = a%val*1.125_dp
        fixture%val = a%val
        b = b*1.125_dp
        start = wall_time()
        call sparse_refactor(solver, a, status, reused)
        refactor_s = wall_time() - start
        if (.not. status_ok(status)) error stop "refactor failed"
    end if
    start = wall_time()
    if (nrhs == 1) then
        call sparse_solve(solver, b(:, 1), x(:, 1), status)
    else
        call sparse_solve_block(solver, b, x, status)
    end if
    solve_s = wall_time() - start
    if (.not. status_ok(status)) error stop "solve failed"
    error = 0.0_dp
    do j = 1, nrhs
        error = max(error, backward_error(fixture, x(:, j), b(:, j)))
    end do
    if (.not. (error < 1.0e-12_dp)) error stop "backward error failed"
    call sparse_free(solver)
    call cpu_time(cpu_s)
    cpu_s = cpu_s - cpu_start
    total_s = wall_time() - total_start
    call metrics(rss, child_rss, child_cpu)
    write (*, '(a)') "backend,family,level,n,nnz,threads,nrhs,factor_s,"// &
        "solve_s,refactor_s,reused,parent_cpu_s,total_s,parent_rss_kib,"// &
        "helper_rss_kib,helper_cpu_s,factor_nnz,factor_mem_mb,"// &
        "analysis_s,numeric_s,backward_error"
    write (*, '(*(g0,:,","))') trim(backend), trim(family), level, a%nrow, &
        a%nnz, threads, nrhs, factor_s, solve_s, refactor_s, reused, cpu_s, &
        total_s, rss, child_rss, child_cpu, stats%factor_nnz, &
        stats%factor_mem_mb, stats%analysis_s, stats%factor_s, error
end program
