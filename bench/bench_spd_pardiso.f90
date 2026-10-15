program bench_spd_pardiso
    ! This executable links only oneMKL, never libfortsparse/OpenBLAS.
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_int64_t
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spd_fixtures, only: fixture_t, make_fixture, fixture_product, &
        backward_error, wall_time
    implicit none
    interface
        subroutine pardisoinit(pt, mtype, iparm) bind(C, name="pardisoinit")
            import :: c_int, c_int64_t
            integer(c_int64_t), intent(inout) :: pt(*)
            integer(c_int), intent(in) :: mtype
            integer(c_int), intent(out) :: iparm(*)
        end subroutine
        subroutine pardiso(pt, maxfct, mnum, mtype, phase, n, av, ia, ja, &
                perm, nrhs, iparm, msglvl, b, x, error) bind(C, name="pardiso")
            import :: c_int, c_int64_t, c_double
            integer(c_int64_t), intent(inout) :: pt(*)
            integer(c_int), intent(in) :: maxfct, mnum, mtype, phase, n
            real(c_double), intent(in) :: av(*)
            integer(c_int), intent(in) :: ia(*), ja(*), nrhs, msglvl
            integer(c_int), intent(inout) :: perm(*), iparm(*)
            real(c_double), intent(inout) :: b(*), x(*)
            integer(c_int), intent(out) :: error
        end subroutine
        subroutine metrics(rss, child_rss, child_cpu) &
                bind(C, name="fsbench_process_metrics")
            import :: c_double
            real(c_double), intent(out) :: rss, child_rss, child_cpu
        end subroutine
    end interface
    type(fixture_t) :: a
    character(32) :: family, argument, workload
    integer :: level, threads, j, col, p, position, profile_iterations
    integer(c_int64_t) :: pt(64)
    integer(c_int) :: iparm(64), nrhs, code
    integer(c_int), allocatable :: ia(:), ja(:), perm(:)
    real(dp), allocatable :: av(:), b(:, :), x(:, :), reference(:)
    real(dp), allocatable :: primal(:), tangent(:), b_dot(:), matrix_bar(:)
    real(dp), allocatable :: original_b(:, :), original_av(:)
    real(dp) :: start, total_start, cpu_start, cpu_s, total_s, factor_s
    real(dp) :: analysis_s, numeric_s, solve_s, refactor_s, error
    real(dp) :: rss, child_rss, child_cpu, factor_mem
    real(dp) :: jvp_s, vjp_s, cold_factor_s, lhs, rhs
    real(dp) :: profile_start, profile_s, factor_sum, solve_sum, refactor_sum
    real(dp) :: iteration_start, predicted_s, scale
    logical :: held_out, reused

    total_start = wall_time()
    call cpu_time(cpu_start)
    call get_command_argument(1, family)
    call get_command_argument(2, argument)
    read (argument, *) level
    call get_command_argument(3, argument)
    read (argument, *) threads
    call get_command_argument(4, workload)
    call get_command_argument(5, argument)
    held_out = trim(argument) == "heldout"
    nrhs = 1
    if (trim(workload) == "rhs8") nrhs = 8
    if (trim(workload) /= "rhs1" .and. trim(workload) /= "rhs8" &
        .and. trim(workload) /= "refactor" .and. &
        trim(workload) /= "derivatives" .and. trim(workload) /= "warm" &
        .and. trim(workload) /= "profile-factor" .and. &
        trim(workload) /= "profile-refactor") &
        error stop "unknown workload"
    call make_fixture(trim(family), level, held_out, a)
    allocate (ia(a%n + 1), ja(size(a%val)), av(size(a%val)), perm(a%n))
    av = 0.0_dp
    allocate (b(a%n, nrhs), x(a%n, nrhs), reference(a%n))
    position = 1
    do col = 1, a%n
        ia(col) = position
        do p = a%ptr(col), a%ptr(col + 1) - 1
            if (a%row(p) < col) cycle
            ! Symmetric CSC's lower triangle is the upper triangle in CSR.
            ja(position) = a%row(p)
            av(position) = a%val(p)
            position = position + 1
        end do
    end do
    ia(a%n + 1) = position
    perm = 0
    do j = 1, nrhs
        reference = 1.0_dp + real(j, dp)/8.0_dp
        reference(1:a%n:2) = -0.5_dp
        call fixture_product(a, reference, b(:, j))
    end do
    pt = 0
    call pardisoinit(pt, 2_c_int, iparm)
    iparm(2) = 2 ! METIS ordering; other settings retain oneMKL defaults.
    start = wall_time()
    call phase(11)
    analysis_s = wall_time() - start
    start = wall_time()
    call phase(22)
    numeric_s = wall_time() - start
    factor_s = analysis_s + numeric_s
    cold_factor_s = factor_s
    if (trim(workload) == "warm") then
        call phase(-1)
        pt = 0
        call pardisoinit(pt, 2_c_int, iparm)
        iparm(2) = 2
        start = wall_time()
        call phase(11)
        analysis_s = wall_time() - start
        start = wall_time()
        call phase(22)
        numeric_s = wall_time() - start
        factor_s = analysis_s + numeric_s
    end if
    factor_mem = real(max(iparm(15), iparm(16) + iparm(17)), dp)/1024.0_dp
    refactor_s = -1.0_dp
    reused = .false.
    if (trim(workload) == "refactor") then
        av = av*1.125_dp
        a%val = a%val*1.125_dp
        b = b*1.125_dp
        start = wall_time()
        call phase(22)
        refactor_s = wall_time() - start
        reused = .true.
    end if
    start = wall_time()
    call phase(33)
    solve_s = wall_time() - start
    profile_iterations = 0
    profile_s = -1.0_dp
    if (trim(workload) == "profile-factor" .or. &
        trim(workload) == "profile-refactor") then
        predicted_s = factor_s + solve_s
        if (wall_time() - total_start + predicted_s > 0.9_dp) &
            error stop "profiling case exceeds one-second replay budget"
        allocate (original_b(a%n, nrhs), original_av(size(av)))
        original_b = b
        original_av = av
        factor_sum = 0.0_dp
        solve_sum = 0.0_dp
        refactor_sum = 0.0_dp
        profile_start = wall_time()
        do
            iteration_start = wall_time()
            profile_iterations = profile_iterations + 1
            if (trim(workload) == "profile-factor") then
                call phase(-1)
                pt = 0
                call pardisoinit(pt, 2_c_int, iparm)
                iparm(2) = 2
                start = wall_time()
                call phase(11)
                call phase(22)
                factor_sum = factor_sum + wall_time() - start
            else
                scale = 1.0_dp
                if (mod(profile_iterations, 2) == 1) scale = 1.125_dp
                av = original_av*scale
                b = original_b*scale
                start = wall_time()
                call phase(22)
                refactor_sum = refactor_sum + wall_time() - start
            end if
            start = wall_time()
            call phase(33)
            solve_sum = solve_sum + wall_time() - start
            predicted_s = wall_time() - iteration_start
            if (wall_time() - profile_start >= 0.2_dp) exit
            if (wall_time() - total_start + predicted_s > 0.9_dp) exit
        end do
        profile_s = wall_time() - profile_start
        solve_s = solve_sum/real(profile_iterations, dp)
        if (trim(workload) == "profile-factor") then
            factor_s = factor_sum/real(profile_iterations, dp)
        else
            refactor_s = refactor_sum/real(profile_iterations, dp)
            a%val = a%val*scale
            reused = .true.
        end if
    end if
    error = 0.0_dp
    do j = 1, nrhs
        error = max(error, backward_error(a, x(:, j), b(:, j)))
    end do
    if (.not. (error < 1.0e-12_dp)) error stop "backward error failed"
    jvp_s = -1.0_dp
    vjp_s = -1.0_dp
    if (trim(workload) == "derivatives") then
        allocate (primal(a%n), tangent(a%n), b_dot(a%n), &
            matrix_bar(size(a%val)))
        primal = x(:, 1)
        b_dot = 0.2_dp*b(:, 1)
        start = wall_time()
        call fixture_product(a, primal, reference)
        b(:, 1) = b_dot - 0.1_dp*reference
        call phase(33)
        tangent = x(:, 1)
        jvp_s = wall_time() - start
        start = wall_time()
        b(:, 1) = 1.0_dp
        call phase(33)
        do col = 1, a%n
            do p = a%ptr(col), a%ptr(col + 1) - 1
                matrix_bar(p) = -x(a%row(p), 1)*primal(col)
            end do
        end do
        vjp_s = wall_time() - start
        if (maxval(abs(tangent - 0.1_dp*primal)) > 1.0e-10_dp) &
            error stop "JVP known-solution oracle failed"
        lhs = sum(tangent)
        rhs = dot_product(x(:, 1), b_dot) + &
            0.1_dp*dot_product(matrix_bar, a%val)
        if (abs(lhs - rhs) > 1.0e-10_dp*max(1.0_dp, abs(lhs))) &
            error stop "derivative adjoint identity failed"
    end if
    call phase(-1)
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
    write (*, '(*(g0,:,","))') "mkl-pardiso", trim(family), level, a%n, &
        size(a%val), threads, nrhs, factor_s, solve_s, refactor_s, reused, &
        cpu_s, total_s, rss, child_rss, child_cpu, iparm(18), factor_mem, &
        analysis_s, numeric_s, error, jvp_s, vjp_s, cold_factor_s, &
        threads, -1, -1, -1, profile_iterations, profile_s

contains
    subroutine phase(number)
        integer(c_int), intent(in) :: number
        call pardiso(pt, 1_c_int, 1_c_int, 2_c_int, number, int(a%n, c_int), &
            av, ia, ja, perm, nrhs, iparm, 0_c_int, b, x, code)
        if (code /= 0) then
            write (*, '(a,i0,a,i0)') "PARDISO phase ", number, " error ", code
            stop 2
        end if
    end subroutine
end program
