program bench_spd
    use fortsparse, only: dp, csc_t, sparse_solver_t, fortsparse_status_t, &
        status_ok, sparse_factor, sparse_solve, sparse_solve_block, &
        sparse_refactor, sparse_free, sparse_factor_stats_t, &
        sparse_factor_stats, sparse_solve_jvp, sparse_solve_vjp, &
        FORTSPARSE_SYM_SPD
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
    type(csc_t) :: a, a_dot
    type(sparse_solver_t) :: solver
    type(fortsparse_status_t) :: status
    type(sparse_factor_stats_t) :: stats
    character(32) :: backend, family, argument, workload
    integer :: level, threads, nrhs, j, profile_iterations
    logical :: held_out, reused, reused_one
    real(dp), allocatable :: b(:, :), x(:, :), reference(:)
    real(dp), allocatable :: tangent(:), adjoint(:), matrix_bar(:), b_dot(:)
    real(dp), allocatable :: original_b(:, :)
    real(dp) :: start, cpu_start, factor_s, solve_s, refactor_s, cpu_s
    real(dp) :: error, rss, child_rss, child_cpu, total_start, total_s
    real(dp) :: jvp_s, vjp_s, lhs, rhs, cold_factor_s
    real(dp) :: profile_start, profile_s, factor_sum, solve_sum, refactor_sum
    real(dp) :: iteration_start, predicted_s, scale

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
    case ("native-cuda")
        solver%backend_id = 4
        solver%options%device = 0
    case ("native-reference")
        solver%backend_id = 5
    case default
        error stop "unknown backend"
    end select
    if (trim(workload) /= "rhs1" .and. trim(workload) /= "rhs8" &
        .and. trim(workload) /= "refactor" .and. &
        trim(workload) /= "derivatives" .and. trim(workload) /= "warm" &
        .and. trim(workload) /= "profile-factor" .and. &
        trim(workload) /= "profile-refactor") &
        error stop "unknown workload"
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
    cold_factor_s = factor_s
    if (trim(workload) == "warm") then
        ! Retain loaded numerical runtimes but redo analysis and factorization.
        call sparse_free(solver)
        start = wall_time()
        call sparse_factor(solver, a, status)
        factor_s = wall_time() - start
        if (.not. status_ok(status)) error stop "warm factor failed"
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
    profile_iterations = 0
    profile_s = -1.0_dp
    if (trim(workload) == "profile-factor" .or. &
        trim(workload) == "profile-refactor") then
        predicted_s = factor_s + solve_s
        if (wall_time() - total_start + predicted_s > 0.9_dp) &
            error stop "profiling case exceeds one-second replay budget"
        allocate (original_b(a%nrow, nrhs))
        original_b = b
        factor_sum = 0.0_dp
        solve_sum = 0.0_dp
        refactor_sum = 0.0_dp
        reused = .true.
        profile_start = wall_time()
        do
            iteration_start = wall_time()
            profile_iterations = profile_iterations + 1
            if (trim(workload) == "profile-factor") then
                call sparse_free(solver)
                start = wall_time()
                call sparse_factor(solver, a, status)
                factor_sum = factor_sum + wall_time() - start
            else
                scale = 1.0_dp
                if (mod(profile_iterations, 2) == 1) scale = 1.125_dp
                a%val = fixture%val*scale
                b = original_b*scale
                start = wall_time()
                call sparse_refactor(solver, a, status, reused_one)
                refactor_sum = refactor_sum + wall_time() - start
                reused = reused .and. reused_one
            end if
            if (.not. status_ok(status)) error stop "profiling factor failed"
            start = wall_time()
            call sparse_solve(solver, b(:, 1), x(:, 1), status)
            solve_sum = solve_sum + wall_time() - start
            if (.not. status_ok(status)) error stop "profiling solve failed"
            predicted_s = wall_time() - iteration_start
            if (wall_time() - profile_start >= 0.2_dp) exit
            if (wall_time() - total_start + predicted_s > 0.9_dp) exit
        end do
        profile_s = wall_time() - profile_start
        solve_s = solve_sum/real(profile_iterations, dp)
        if (trim(workload) == "profile-factor") then
            factor_s = factor_sum/real(profile_iterations, dp)
            reused = .false.
        else
            refactor_s = refactor_sum/real(profile_iterations, dp)
            fixture%val = a%val
        end if
        call sparse_factor_stats(solver, stats, status)
    end if
    error = 0.0_dp
    do j = 1, nrhs
        error = max(error, backward_error(fixture, x(:, j), b(:, j)))
    end do
    if (.not. (error < 1.0e-12_dp)) error stop "backward error failed"
    jvp_s = -1.0_dp
    vjp_s = -1.0_dp
    if (trim(workload) == "derivatives") then
        a_dot = a
        a_dot%val = 0.1_dp*a%val
        allocate (tangent(a%nrow), adjoint(a%nrow), matrix_bar(a%nnz), &
            b_dot(a%nrow))
        b_dot = 0.2_dp*b(:, 1)
        reference = 1.0_dp
        start = wall_time()
        call sparse_solve_jvp(solver, a_dot, x(:, 1), b_dot, tangent, status)
        jvp_s = wall_time() - start
        if (.not. status_ok(status)) error stop "JVP failed"
        start = wall_time()
        ! A is symmetric, so its factors also solve the adjoint equation.
        call sparse_solve_vjp(solver, a, x(:, 1), reference, adjoint, &
            matrix_bar, status)
        vjp_s = wall_time() - start
        if (.not. status_ok(status)) error stop "VJP failed"
        if (maxval(abs(tangent - 0.1_dp*x(:, 1))) > 1.0e-10_dp) &
            error stop "JVP known-solution oracle failed"
        lhs = dot_product(reference, tangent)
        rhs = dot_product(adjoint, b_dot) + dot_product(matrix_bar, a_dot%val)
        if (abs(lhs - rhs) > 1.0e-10_dp*max(1.0_dp, abs(lhs))) &
            error stop "derivative adjoint identity failed"
    end if
    call sparse_free(solver)
    call cpu_time(cpu_s)
    cpu_s = cpu_s - cpu_start
    total_s = wall_time() - total_start
    call metrics(rss, child_rss, child_cpu)
    write (*, '(a)') "backend,family,level,n,nnz,threads,nrhs,factor_s,"// &
        "solve_s,refactor_s,reused,parent_cpu_s,total_s,parent_rss_kib,"// &
        "helper_rss_kib,helper_cpu_s,factor_nnz,factor_mem_mb,"// &
        "analysis_s,numeric_s,backward_error,jvp_s,vjp_s,cold_factor_s,"// &
        "requested_threads,active_threads,requested_device,actual_device,"// &
        "profile_iterations,profile_s"
    write (*, '(*(g0,:,","))') trim(backend), trim(family), level, a%nrow, &
        a%nnz, threads, nrhs, factor_s, solve_s, refactor_s, reused, cpu_s, &
        total_s, rss, child_rss, child_cpu, stats%factor_nnz, &
        stats%factor_mem_mb, stats%analysis_s, stats%factor_s, error, jvp_s, &
        vjp_s, cold_factor_s, threads, stats%active_threads, &
        solver%options%device, -1, profile_iterations, profile_s
end program
