program bench_spd_cudss
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_int64_t
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use spd_fixtures, only: fixture_t, make_fixture, fixture_known_rhs, &
        fixture_spd_update, backward_error, wall_time
    implicit none
    interface
        integer(c_int) function run(n, ptr, row, values, updated, rhs, solution, &
                nrhs, threads, device, workload, resident, ordering, permutation, &
                repeats, times, counts) bind(C, name="fsbench_cudss_run")
            import :: c_int, c_double, c_int64_t
            integer(c_int), value :: n, nrhs, threads, device, workload, resident
            integer(c_int), value :: ordering, repeats
            integer(c_int), intent(in) :: ptr(*), row(*), permutation(*)
            real(c_double), intent(in) :: values(*), updated(*), rhs(*)
            real(c_double), intent(out) :: solution(*), times(13)
            integer(c_int64_t), intent(out) :: counts(21)
        end function run
        subroutine metrics(rss, child_rss, child_cpu) &
                bind(C, name="fsbench_process_metrics")
            import :: c_double
            real(c_double), intent(out) :: rss, child_rss, child_cpu
        end subroutine metrics
    end interface
    type(fixture_t) :: a, changed
    character(32) :: family, argument, workload, variant, io_mode, ordering
    character(1024) :: matrix_path, permutation_path
    character(16) :: ordering_hash
    integer :: level, threads, device, mode, resident, order, nrhs, repeats
    integer :: r, j, unit, ios, permutation_size, expected_minor
    integer(c_int) :: code
    integer(c_int), allocatable :: permutation(:)
    logical, allocatable :: seen(:)
    real(dp), allocatable :: rhs(:, :, :), reference(:, :, :), solution(:, :, :)
    real(dp) :: start, cpu_start, cpu_s, total_s, check_start, check_s
    real(dp) :: error, forward_error, rss, child_rss, child_cpu
    real(c_double) :: times(13)
    integer(c_int64_t) :: counts(21)

    start = wall_time()
    call cpu_time(cpu_start)
    call get_command_argument(1, family)
    call get_command_argument(2, argument)
    read (argument, *) level
    call get_command_argument(3, argument)
    read (argument, *) threads
    call get_command_argument(4, workload)
    call get_command_argument(5, variant)
    call get_command_argument(6, argument)
    read (argument, *) device
    call get_command_argument(7, io_mode)
    call get_command_argument(8, ordering)
    call get_command_argument(9, matrix_path)
    if (trim(matrix_path) == '-') matrix_path = ''
    call get_command_argument(10, argument)
    repeats = 1
    if (trim(workload) == 'repeated-rhs') repeats = 8
    if (len_trim(argument) > 0) read (argument, *) repeats
    call get_command_argument(11, permutation_path)
    ! Position 12 is only a registered rejection oracle, never a benchmark knob.
    call get_command_argument(12, argument)
    expected_minor = 0
    if (len_trim(argument) > 0) read (argument, *) expected_minor
    if (threads < 1 .or. device < 0 .or. repeats < 1 .or. repeats > 32) &
        error stop 'invalid thread/device/repeat argument'
    if (trim(variant) /= 'training' .and. trim(variant) /= 'heldout') &
        error stop 'unknown fixture variant'
    resident = 0
    if (trim(io_mode) == 'resident') resident = 1
    if (trim(io_mode) /= 'host' .and. trim(io_mode) /= 'resident') &
        error stop 'unknown IO boundary'
    nrhs = 1
    mode = 0
    select case (trim(workload))
    case ('rhs1')
    case ('rhs8')
        nrhs = 8
    case ('warm')
        mode = 2
    case ('refactor')
        mode = 3
    case ('repeated-rhs')
        mode = 4
    case default
        error stop 'unknown GPU workload'
    end select
    if (mode /= 4 .and. repeats /= 1) error stop 'repeats need repeated-rhs'
    select case (trim(ordering))
    case ('tuned')
        order = 0
    case ('amd')
        order = 1
    case ('nd')
        order = 2
    case ('common')
        order = 3
    case default
        error stop 'unknown ordering'
    end select
    call make_fixture(trim(family), level, trim(variant) == 'heldout', a, &
        trim(matrix_path))
    allocate (permutation(a%n), seen(0:a%n - 1))
    permutation = [(j - 1, j=1, a%n)]
    if (order == 3) then
        open (newunit=unit, file=trim(permutation_path), status='old', &
            action='read', iostat=ios)
        if (ios /= 0) error stop 'common permutation file missing'
        read (unit, *, iostat=ios) permutation_size
        if (ios /= 0) error stop 'invalid common permutation size'
        if (permutation_size /= a%n) error stop 'common permutation size mismatch'
        read (unit, *, iostat=ios) permutation
        if (ios /= 0) error stop 'incomplete common permutation'
        close (unit)
        seen = .false.
        do j = 1, a%n
            if (permutation(j) < 0 .or. permutation(j) >= a%n) &
                error stop 'common permutation outside graph'
            if (seen(permutation(j))) error stop 'common permutation is not bijective'
            seen(permutation(j)) = .true.
        end do
    end if
    call fixture_spd_update(a, changed)
    allocate (rhs(a%n, nrhs, repeats), reference(a%n, nrhs, repeats), &
        solution(a%n, nrhs, repeats))
    do r = 1, repeats
        if (mode == 3) then
            call fixture_known_rhs(changed, nrhs, r, rhs(:, :, r), reference(:, :, r))
        else
            call fixture_known_rhs(a, nrhs, r, rhs(:, :, r), reference(:, :, r))
        end if
    end do
    solution = 0.0_dp
    code = run(a%n, a%ptr, a%row, a%val, changed%val, rhs, solution, nrhs, &
        threads, device, mode, resident, order, permutation, repeats, times, counts)
    if (code == 0 .and. mode == 2) then
        if (counts(18) /= 1 .or. counts(19) /= 2 .or. &
            counts(20) /= 2 .or. counts(21) /= 2) &
            error stop 'warm workflow did not retain runtime and redo analysis/factors'
    end if
    check_start = wall_time()
    error = -1.0_dp
    forward_error = -1.0_dp
    if (code == 0) then
        error = 0.0_dp
        forward_error = 0.0_dp
        do r = 1, repeats
            do j = 1, nrhs
                if (.not. all(ieee_is_finite(solution(:, j, r)))) then
                    code = 5
                    cycle
                end if
                if (mode == 3) then
                    error = max(error, backward_error(changed, &
                        solution(:, j, r), rhs(:, j, r)))
                else
                    error = max(error, backward_error(a, &
                        solution(:, j, r), rhs(:, j, r)))
                end if
                forward_error = max(forward_error, &
                    maxval(abs(solution(:, j, r) - reference(:, j, r)))/ &
                    maxval(abs(reference(:, j, r))))
            end do
        end do
        if (error > 1.0e-12_dp .or. forward_error > 1.0e-8_dp) code = 5
    end if
    check_s = wall_time() - check_start
    call cpu_time(cpu_s)
    cpu_s = cpu_s - cpu_start
    total_s = wall_time() - start
    call metrics(rss, child_rss, child_cpu)
    write (ordering_hash, '(z16.16)') counts(11)
    write (*, '(a)') 'backend,family,level,n,nnz,threads,nrhs,workload,variant,'// &
        'actual_device,io_mode,ordering,repeats,vendor_status,backward_error,'// &
        'factor_s,cold_factor_s,analysis_s,numeric_s,refactor_s,solve_s,total_s,'// &
        'parent_cpu_s,parent_rss_kib,helper_rss_kib,helper_cpu_s,runtime_setup_s,'// &
        'conversion_s,h2d_s,d2h_s,synchronization_s,cleanup_s,'// &
        'device_requested_peak_bytes,device_used_peak_bytes,'// &
        'device_reserved_peak_bytes,wrapper_device_bytes,'// &
        'device_process_sampled_max_bytes,device_total_peak_bytes,'// &
        'pinned_host_bytes,factor_nnz,requested_threads,active_threads,'// &
        'ordering_hash,solver_version,execution_mode,gpu_kernel_calls,'// &
        'gpu_offload_verified,device_memory_estimate_peak_bytes,'// &
        'host_memory_estimate_peak_bytes,check_s,forward_error,'// &
        'pool_limit_bytes,memory_sampling_s,device_process_sample_count,'// &
        'runtime_handle_creations,data_object_creations,'// &
        'symbolic_analysis_calls,numeric_factor_calls'
    write (*, '(*(g0,:,","))') 'cudss', trim(family), level, a%n, size(a%val), &
        threads, nrhs, trim(workload), trim(variant), counts(9), trim(io_mode), &
        trim(ordering), repeats, code, error, times(11), times(12), times(3), &
        times(4), times(5), times(6), total_s, cpu_s, rss, child_rss, child_cpu, &
        times(1), times(2), times(7), times(8), times(9), times(10), &
        counts(1:8), threads, counts(10), ordering_hash, '0.8.0.10-cuda13', &
        'gpu-only-factor-solve', -1, 0, counts(14:15), check_s, forward_error, &
        counts(16), times(13), counts(17:21)
    if (expected_minor > 0) then
        if (code /= 1000 + expected_minor) &
            error stop 'expected exact vendor nonpositive minor was not reported'
        print '(a)', 'PASS cuDSS rejected expected nonpositive minor'
        stop
    end if
    if (code /= 0) stop 2
end program bench_spd_cudss
