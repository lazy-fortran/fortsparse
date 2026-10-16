program bench_spd_cholmod
    ! Standalone GPL comparison: fixture/oracle code only, no FortSparse link.
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_int64_t, c_double, &
        c_associated
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use spd_fixtures, only: fixture_t, make_fixture, fixture_product, &
        backward_error, wall_time, fixture_spd_update, fixture_known_rhs
    implicit none
    interface
        function create(n, nnz, ptr, row, val, threads, device, ordering, &
                perm, pool, times, status) bind(C, name="fsbench_cholmod_create") &
                result(context)
            import c_ptr, c_int, c_int64_t, c_double
            integer(c_int), value :: n, nnz, threads, device, ordering
            integer(c_int), intent(in) :: ptr(*), row(*)
            real(c_double), intent(in) :: val(*)
            integer(c_int64_t), intent(in) :: perm(*)
            integer(c_int64_t), value :: pool
            real(c_double), intent(out) :: times(*)
            integer(c_int), intent(out) :: status
            type(c_ptr) :: context
        end function
        function analyze(context) bind(C, name="fsbench_cholmod_analyze") &
                result(code)
            import c_ptr, c_int
            type(c_ptr), value :: context
            integer(c_int) :: code
        end function
        function factor(context, values) &
                bind(C, name="fsbench_cholmod_factor") result(code)
            import c_ptr, c_int, c_double
            type(c_ptr), value :: context
            real(c_double), intent(in) :: values(*)
            integer(c_int) :: code
        end function
        function solve(context, nrhs, b, x) &
                bind(C, name="fsbench_cholmod_solve") result(code)
            import c_ptr, c_int, c_double
            type(c_ptr), value :: context
            integer(c_int), value :: nrhs
            real(c_double), intent(in) :: b(*)
            real(c_double), intent(out) :: x(*)
            integer(c_int) :: code
        end function
        subroutine stats(context, values, perm) bind(C, name="fsbench_cholmod_stats")
            import c_ptr, c_double, c_int64_t
            type(c_ptr), value :: context
            real(c_double), intent(out) :: values(*)
            integer(c_int64_t), intent(out) :: perm(*)
        end subroutine
        function ordering_hash(context) bind(C, name="fsbench_cholmod_ordering_hash") &
                result(hash)
            import c_ptr, c_int64_t
            type(c_ptr), value :: context
            integer(c_int64_t) :: hash
        end function
        subroutine version(numbers) bind(C, name="fsbench_cholmod_version")
            import c_int
            integer(c_int), intent(out) :: numbers(3)
        end subroutine
        subroutine destroy(context) bind(C, name="fsbench_cholmod_free")
            import c_ptr
            type(c_ptr), value :: context
        end subroutine
        subroutine metrics(rss, child_rss, child_cpu) &
                bind(C, name="fsbench_process_metrics")
            import c_double
            real(c_double), intent(out) :: rss, child_rss, child_cpu
        end subroutine
    end interface
    type(fixture_t) :: a
    type(c_ptr) :: context
    character(32) :: family, arg, workload, variant, io_mode, ordering
    character(32) :: order_hash, version_string, execution_mode
    character(1024) :: matrix_path, permutation_path, permutation_out
    integer :: level, threads, device, nrhs, repeats, j, k, unit, ios, m
    integer :: offload_verified
    integer(c_int) :: code, order_id, version_numbers(3)
    integer(c_int64_t) :: pool
    integer(c_int64_t), allocatable :: permutation(:), actual_permutation(:)
    logical, allocatable :: seen(:)
    real(dp), allocatable :: b(:, :), x(:, :), reference(:, :)
    real(dp) :: total_start, cpu_start, total_s, cpu_s, start
    real(dp) :: analysis_s, numeric_s, factor_s, cold_factor_s, refactor_s, solve_s
    real(dp) :: cleanup_s, gpu_calls
    real(dp) :: setup_times(2), info(10), rss, child_rss, child_cpu
    real(dp) :: error, forward_error, fixture_bytes, host_work_bytes
    logical :: reused

    total_start = wall_time()
    call cpu_time(cpu_start)
    call get_command_argument(1, family)
    call get_command_argument(2, arg)
    read (arg, *) level
    call get_command_argument(3, arg)
    read (arg, *) threads
    call get_command_argument(4, workload)
    call get_command_argument(5, variant)
    call get_command_argument(6, arg)
    read (arg, *) device
    call get_command_argument(7, io_mode)
    call get_command_argument(8, ordering)
    call get_command_argument(9, matrix_path)
    call get_command_argument(10, arg)
    repeats = 8
    if (len_trim(arg) > 0) read (arg, *) repeats
    call get_command_argument(11, permutation_path)
    if (io_mode /= "host") then
        write (*, '(a)') "CHOLMOD resident IO unavailable: "// &
            "factors and solves reside on host"
        stop 77
    end if
    if (threads < 1 .or. device < 0 .or. repeats < 1) error stop "invalid arguments"
    select case (trim(workload))
    case ("rhs1", "rhs8", "warm", "refactor", "repeated-rhs")
    case default
        error stop "unknown workload"
    end select
    call make_fixture(trim(family), level, variant == "heldout", a, trim(matrix_path))
    nrhs = 1
    if (workload == "rhs8") nrhs = 8
    allocate (b(a%n, nrhs), x(a%n, nrhs), reference(a%n, nrhs))
    allocate (permutation(a%n), actual_permutation(a%n), seen(a%n))
    permutation = 0
    order_id = 0
    select case (trim(ordering))
    case ("tuned")
    case ("amd")
        order_id = 1
    case ("nd")
        order_id = 2
    case ("common")
        order_id = 3
        open (newunit=unit, file=trim(permutation_path), status="old", iostat=ios)
        if (ios /= 0) error stop "common permutation file unavailable"
        read (unit, *, iostat=ios) m
        if (ios /= 0) error stop "invalid permutation dimension"
        if (m /= a%n) error stop "common permutation dimension mismatch"
        read (unit, *, iostat=ios) permutation
        close (unit)
        if (ios /= 0) error stop "invalid common permutation"
        seen = .false.
        do j = 1, a%n
            k = int(permutation(j)) + 1
            if (k < 1 .or. k > a%n) error stop "permutation index out of range"
            if (seen(k)) error stop "permutation is not bijective"
            seen(k) = .true.
        end do
    case default
        error stop "unknown ordering"
    end select
    call fixture_known_rhs(a, nrhs, 1, b, reference)
    pool = 268435456_c_int64_t
    call get_environment_variable("CHOLMOD_GPU_POOL_BYTES", arg)
    if (len_trim(arg) > 0) read (arg, *) pool
    if (pool < 67108864_c_int64_t) error stop "GPU pool must be at least 64 MiB"
    context = create(a%n, size(a%val), a%ptr, a%row, a%val, threads, device, &
        order_id, permutation, pool, setup_times, code)
    if (.not. c_associated(context)) then
        write (*, '(a,i0)') "CHOLMOD create error ", code
        stop 2
    end if
    call run_factor()
    factor_s = factor_s + sum(setup_times)
    cold_factor_s = factor_s
    if (workload == "warm") call run_factor()
    refactor_s = -1.0_dp
    reused = .false.
    if (workload == "refactor") then
        block
            type(fixture_t) :: changed
            call fixture_spd_update(a, changed)
            call move_alloc(changed%val, a%val)
        end block
        call fixture_known_rhs(a, nrhs, 1, b, reference)
        start = wall_time()
        code = factor(context, a%val)
        call require_success("refactor")
        refactor_s = wall_time() - start
        reused = .true.
    end if
    solve_s = 0.0_dp
    if (workload /= "repeated-rhs") repeats = 1
    do j = 1, repeats
        call fixture_known_rhs(a, nrhs, j, b, reference)
        start = wall_time()
        code = solve(context, nrhs, b, x)
        call require_success("solve")
        solve_s = solve_s + wall_time() - start
    end do
    if (any(.not. ieee_is_finite(x))) error stop "nonfinite solution"
    error = 0.0_dp
    forward_error = 0.0_dp
    do j = 1, nrhs
        error = max(error, backward_error(a, x(:, j), b(:, j)))
        forward_error = max(forward_error, maxval(abs(x(:, j) - reference(:, j)))/ &
            max(1.0_dp, maxval(abs(reference(:, j)))))
    end do
    if (.not. ieee_is_finite(error)) error stop "nonfinite backward error"
    if (.not. (error < 1.0e-12_dp)) error stop "independent residual failed"
    if (.not. ieee_is_finite(forward_error)) error stop "nonfinite forward error"
    if (forward_error > 1.0e-8_dp) error stop "known-solution oracle failed"
    call stats(context, info, actual_permutation)
    write (order_hash, '(z16.16)') ordering_hash(context)
    call version(version_numbers)
    write (version_string, '(i0,a,i0,a,i0)') version_numbers(1), ".", &
        version_numbers(2), ".", version_numbers(3)
    gpu_calls = sum(info(5:8))
    ! The public API does not expose whether BLAS_TIMER counters are built.
    ! Positive counters are evidence; zero is inconclusive in ordinary builds.
    execution_mode = "offload-unknown"
    offload_verified = -1
    if (gpu_calls > 0) then
        execution_mode = "gpu-factor-host-solve"
        offload_verified = 1
    else
        gpu_calls = -1
    end if
    if (order_id == 3) then
        if (any(actual_permutation /= permutation)) &
            error stop "common permutation changed"
    end if
    call get_environment_variable("CHOLMOD_PERMUTATION_OUT", permutation_out)
    if (len_trim(permutation_out) > 0) then
        open (newunit=unit, file=trim(permutation_out), status="replace")
        write (unit, '(i0)') a%n
        write (unit, '(i0)') actual_permutation
        close (unit)
    end if
    fixture_bytes = real(size(a%ptr), dp)*4 + real(size(a%row), dp)*4 + &
        real(size(a%val), dp)*8
    host_work_bytes = real(size(b) + size(x) + size(reference), dp)*8 + &
        real(size(permutation) + size(actual_permutation), dp)*8 + &
        real(size(seen), dp)*storage_size(.true.)/8
    start = wall_time()
    call destroy(context)
    cleanup_s = wall_time() - start
    call cpu_time(cpu_s)
    cpu_s = cpu_s - cpu_start
    total_s = wall_time() - total_start
    call metrics(rss, child_rss, child_cpu)
    write (*, '(a)') "backend,family,level,n,nnz,threads,nrhs,workload,io_mode,"// &
        "ordering,factor_s,solve_s,refactor_s,reused,parent_cpu_s,total_s,"// &
        "parent_rss_kib,analysis_s,numeric_s,backward_error,forward_error,"// &
        "cold_factor_s,runtime_setup_s,conversion_s,h2d_s,d2h_s,synchronization_s,"// &
        "requested_threads,active_threads,requested_device,actual_device,"// &
        "repeats,wrapper_device_bytes,device_requested_peak_bytes,"// &
        "device_used_peak_bytes,device_reserved_peak_bytes,"// &
        "device_process_sampled_max_bytes,cholmod_device_pool_peak_bytes,"// &
        "cholmod_host_peak_bytes,wrapper_host_bytes,"// &
        "pinned_host_bytes,fixture_bytes,host_work_bytes,factor_nnz,"// &
        "gpu_gemm_calls,gpu_syrk_calls,gpu_trsm_calls,gpu_potrf_calls,"// &
        "gpu_pool_budget_bytes,"// &
        "variant,cleanup_s,ordering_hash,solver_version,execution_mode,"// &
        "gpu_operation_calls,"// &
        "gpu_offload_verified,fixture_update_peak_extra_bytes"
    write (*, '(*(g0,:,","))') "cholmod-cuda", trim(family), level, a%n, &
        size(a%val), threads, nrhs, trim(workload), trim(io_mode), trim(ordering), &
        factor_s, solve_s, refactor_s, reused, cpu_s, total_s, rss, &
        analysis_s, numeric_s, error, forward_error, cold_factor_s, &
        setup_times, -1.0_dp, -1.0_dp, -1.0_dp, threads, -1, device, int(info(10)), &
        repeats, 0, -1, -1, -1, -1, info(3), info(1), info(2), info(4), &
        fixture_bytes, host_work_bytes, info(9), info(5:8), pool, &
        trim(variant), cleanup_s, trim(order_hash), trim(version_string), &
        trim(execution_mode), gpu_calls, offload_verified, &
        merge(fixture_bytes, 0.0_dp, workload == "refactor")

contains
    subroutine require_success(phase)
        character(*), intent(in) :: phase
        if (code == 0) return
        write (*, '(a,a,a,i0)') "CHOLMOD ", phase, " error ", code
        call destroy(context)
        stop 2
    end subroutine

    subroutine run_factor()
        start = wall_time()
        code = analyze(context)
        call require_success("analyze")
        analysis_s = wall_time() - start
        start = wall_time()
        code = factor(context, a%val)
        call require_success("factor")
        numeric_s = wall_time() - start
        factor_s = analysis_s + numeric_s
    end subroutine

end program
