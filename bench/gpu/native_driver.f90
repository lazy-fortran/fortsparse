program native_gpu_comparison
    use fortsparse, only: dp, csc_t, sparse_solver_t, fortsparse_status_t, &
        sparse_factor_stats_t, sparse_factor, sparse_refactor, sparse_solve_block, &
        sparse_free, sparse_factor_stats, status_ok, FORTSPARSE_SYM_SPD, &
        FORTSPARSE_ORDER_AMD, FORTSPARSE_ORDER_ND
    use spd_fixtures, only: fixture_t, make_fixture, fixture_spd_update, &
        fixture_known_rhs, backward_error, wall_time
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_int64_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    interface
        subroutine metrics(rss, child_rss, child_cpu) &
                bind(C, name='fsbench_process_metrics')
            import :: c_double
            real(c_double), intent(out) :: rss, child_rss, child_cpu
        end subroutine
        integer(c_int64_t) function device_memory(device) &
                bind(C, name='fsbench_gpu_process_memory')
            import :: c_int, c_int64_t
            integer(c_int), value :: device
        end function
    end interface
    type(fixture_t) :: fixture, changed
    type(csc_t) :: a
    type(sparse_solver_t) :: solver
    type(fortsparse_status_t) :: status
    type(sparse_factor_stats_t) :: stats
    character(32) :: family, workload, variant, io_mode, ordering, argument
    character(32) :: backend, execution_mode
    character(1024) :: matrix_path
    integer :: level, threads, device, nrhs, repeats, j, repetition, n, nnz
    integer :: offload, device_unknown
    integer :: sample_count
    integer(c_int64_t) :: sampled_device_bytes, sample_bytes
    real(dp), allocatable :: rhs(:, :), reference(:, :), x(:, :)
    real(dp) :: total_start, cpu_start, start, factor_s, cold_factor_s, solve_s
    real(dp) :: refactor_s, total_s, cpu_s, cleanup_s, error, known_error
    real(dp) :: rss, child_rss, child_cpu, conversion_s, h2d_s, d2h_s
    real(dp) :: sampling_s
    logical :: reused

    total_start = wall_time()
    call cpu_time(cpu_start)
    sampled_device_bytes = -1_c_int64_t
    sample_count = 0
    sampling_s = 0.0_dp
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
    if (matrix_path == '-') matrix_path = ''
    repeats = 1
    if (workload == 'repeated-rhs') repeats = 8
    call get_command_argument(10, argument)
    if (len_trim(argument) > 0) read (argument, *) repeats
    if (repeats < 1) error stop 'positive repeat count required'
    if (workload /= 'repeated-rhs' .and. repeats /= 1) &
        error stop 'repeats apply only to retained-factor RHS workflow'
    if (io_mode /= 'host') error stop 'native resident workflow unavailable'
    if (variant /= 'training' .and. variant /= 'heldout') &
        error stop 'unknown coefficient variant'
    select case (workload)
    case ('rhs1', 'rhs8', 'warm', 'refactor', 'repeated-rhs')
    case default
        error stop 'unsupported native workload'
    end select
    select case (ordering)
    case ('tuned')
    case ('amd')
        solver%options%ordering = FORTSPARSE_ORDER_AMD
    case ('nd')
        solver%options%ordering = FORTSPARSE_ORDER_ND
    case default
        error stop 'native supplied permutation unavailable'
    end select
    backend = 'native'
    execution_mode = 'cpu'
    if (device >= 0) then
        backend = 'native-cuda'
        execution_mode = 'gpu-factor-host-solve'
    end if
    solver%backend_id = 4
    solver%options%symmetry = FORTSPARSE_SYM_SPD
    solver%options%nthreads = threads
    solver%options%device = device
    call make_fixture(trim(family), level, variant == 'heldout', fixture, &
        trim(matrix_path))
    n = fixture%n
    nnz = size(fixture%val)
    start = wall_time()
    a%nrow = n
    a%ncol = n
    a%nnz = nnz
    a%col_ptr = fixture%ptr
    a%row_idx = fixture%row
    a%val = fixture%val
    conversion_s = wall_time() - start
    nrhs = 1
    if (workload == 'rhs8') nrhs = 8
    allocate (rhs(n, nrhs), reference(n, nrhs), x(n, nrhs))
    call fixture_known_rhs(fixture, nrhs, 1, rhs, reference)
    start = wall_time()
    call sparse_factor(solver, a, status)
    cold_factor_s = wall_time() - start + conversion_s
    call require_ok('factor')
    call sample_device()
    factor_s = cold_factor_s
    if (workload == 'warm') then
        call sparse_free(solver)
        start = wall_time()
        call sparse_factor(solver, a, status)
        factor_s = wall_time() - start
        call require_ok('warm factor')
        call sample_device()
    end if
    refactor_s = -1.0_dp
    if (workload == 'refactor') then
        call fixture_spd_update(fixture, changed)
        call move_alloc(changed%val, fixture%val)
        call fixture_known_rhs(fixture, nrhs, 1, rhs, reference)
        start = wall_time()
        a%val = fixture%val
        call sparse_refactor(solver, a, status, reused)
        refactor_s = wall_time() - start
        call require_ok('refactor')
        call sample_device()
        if (.not. reused) error stop 'native symbolic reuse failed'
    end if
    solve_s = 0.0_dp
    error = 0.0_dp
    do repetition = 1, repeats
        call fixture_known_rhs(fixture, nrhs, repetition, rhs, reference)
        start = wall_time()
        call sparse_solve_block(solver, rhs, x, status)
        solve_s = solve_s + wall_time() - start
        call require_ok('solve')
        call sample_device()
        do j = 1, nrhs
            error = max(error, backward_error(fixture, x(:, j), rhs(:, j)))
            known_error = maxval(abs(x(:, j) - reference(:, j))) / &
                max(1.0_dp, maxval(abs(reference(:, j))))
            if (.not. ieee_is_finite(known_error)) &
                error stop 'native nonfinite solution'
            if (known_error > 1.0e-8_dp) error stop 'native known solution failed'
        end do
    end do
    if (.not. ieee_is_finite(error)) error stop 'native nonfinite residual'
    if (error > 1.0e-12_dp) error stop 'native backward error failed'
    call sparse_factor_stats(solver, stats, status)
    call require_ok('statistics')
    offload = 0
    device_unknown = 0
    h2d_s = 0.0_dp
    d2h_s = 0.0_dp
    if (stats%active_device >= 0) then
        ! Successful native CUDA factorization synchronizes queued device work.
        ! The companion profile supplies kernel-level evidence separately.
        offload = 1
        device_unknown = -1
        h2d_s = -1.0_dp
        d2h_s = -1.0_dp
    end if
    if (device >= 0 .and. offload == 0) error stop 'native CUDA offload absent'
    start = wall_time()
    call sparse_free(solver)
    deallocate (rhs, reference, x, a%col_ptr, a%row_idx, a%val)
    deallocate (fixture%ptr, fixture%row, fixture%val)
    if (allocated(changed%ptr)) deallocate (changed%ptr, changed%row)
    cleanup_s = wall_time() - start
    call cpu_time(cpu_s)
    cpu_s = cpu_s - cpu_start
    total_s = wall_time() - total_start
    call metrics(rss, child_rss, child_cpu)
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
        'gpu_offload_verified,device_process_sample_count,memory_sampling_s'
    write (*, '(*(g0,:,","))') trim(backend), trim(family), level, n, nnz, &
        threads, nrhs, trim(workload), trim(variant), stats%active_device, &
        trim(io_mode), trim(ordering), repeats, 0, error, factor_s, cold_factor_s, &
        stats%analysis_s, stats%factor_s, refactor_s, solve_s, total_s, cpu_s, &
        rss, child_rss, child_cpu, -1.0_dp, conversion_s, h2d_s, d2h_s, &
        -1.0_dp, cleanup_s, device_unknown, device_unknown, device_unknown, &
        stats%device_bytes, sampled_device_bytes, device_unknown, device_unknown, &
        stats%factor_nnz, threads, stats%active_threads, '-1', 'fortsparse-native', &
        trim(execution_mode), -1, offload, sample_count, sampling_s
contains
    subroutine sample_device()
        real(dp) :: sample_start
        if (device < 0) then
            sampled_device_bytes = 0_c_int64_t
            return
        end if
        sample_start = wall_time()
        sample_bytes = device_memory(int(device, c_int))
        sampling_s = sampling_s + wall_time() - sample_start
        if (sample_bytes < 0) return
        sampled_device_bytes = max(sampled_device_bytes, sample_bytes)
        sample_count = sample_count + 1
    end subroutine

    subroutine require_ok(phase)
        character(*), intent(in) :: phase
        if (status_ok(status)) return
        write (*, '(a,a,a,i0)') 'native ', phase, ' status: ', status%code
        stop 2
    end subroutine
end program
