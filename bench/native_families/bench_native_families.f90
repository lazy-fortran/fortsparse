program bench_native_families
    use fortsparse, only: csc_t, csc_z_t, sparse_solver_t, fortsparse_status_t, &
                   sparse_factor, sparse_solve_block_op, sparse_refactor, sparse_free, &
                   sparse_factor_stats, sparse_factor_stats_t, FORTSPARSE_SYM_GENERAL, &
                   FORTSPARSE_SYM_SPD, FORTSPARSE_SYM_INDEFINITE, &
                   FORTSPARSE_ORDER_AMD
    use native_family_fixtures, only: dp, family_fixture_t, &
                                      fixture_backward_error, wall_time
    use native_family_real_io, only: load_benchmark_fixture, load_benchmark_rhs, &
                                    export_real_solution
    use native_family_records, only: family_record_t, process_metrics, &
        parent_cpu_seconds
    implicit none
    type(family_fixture_t) :: fixture
    type(family_record_t) :: r
    type(csc_t) :: ar
    type(csc_z_t) :: az
    type(sparse_solver_t) :: solver
    type(fortsparse_status_t) :: status, stats_status
    type(sparse_factor_stats_t) :: stats
    complex(dp), allocatable :: bz(:, :), xz(:, :), truth(:, :)
    real(dp), allocatable :: br(:, :), xr(:, :)
    character(64) :: argument
    character(1024) :: mini5_path
    real(dp) :: total_start, cpu_start, start, cpu_end
    logical :: heldout
    integer :: operation, j

    total_start = wall_time()
    call cpu_time(cpu_start)
    call get_command_argument(1, r%backend)
    call get_command_argument(2, r%family)
    call get_command_argument(3, argument)
    read (argument, *) r%level
    call get_command_argument(4, argument)
    read (argument, *) r%threads
    call get_command_argument(5, r%workload)
    call get_command_argument(6, r%operation)
    call get_command_argument(7, r%variant)
    call get_command_argument(8, mini5_path)
    if (r%threads < 1) error stop 'positive thread budget required'
    heldout = trim(r%variant) == 'heldout'
    operation = 0
    select case (trim(r%operation))
    case ('N')
        operation = 0
    case ('T')
        operation = 1
    case ('H')
        operation = 2
    case default
        error stop 'operation must be N/T/H'
    end select
    select case (trim(r%backend))
    case ('superlu')
        solver%backend_id = 1
    case ('umfpack')
        solver%backend_id = 2
    case ('mumps')
        solver%backend_id = 3
    case ('native','native-amd')
        solver%backend_id = 4
        if (trim(r%backend) == 'native-amd') &
            solver%options%ordering = FORTSPARSE_ORDER_AMD
    case default
        error stop 'unknown backend'
    end select
    r%nrhs = 1
    select case (trim(r%workload))
    case ('rhs1', 'refactor', 'warm')
    case ('rhs4')
        r%nrhs = 4
    case ('rhs8')
        r%nrhs = 8
    case default
        error stop 'workload must be rhs1/rhs4/rhs8/refactor/warm'
    end select
    call load_benchmark_fixture(trim(r%family), r%level, heldout, fixture, &
                             trim(mini5_path))
    r%n = fixture%n
    r%nnz = size(fixture%val)
    allocate (bz(r%n, r%nrhs), xz(r%n, r%nrhs), truth(r%n, r%nrhs))
    call load_benchmark_rhs(fixture, trim(r%family), r%nrhs, heldout, operation, &
                            bz, truth, trim(mini5_path))
    solver%options%nthreads = r%threads
    solver%options%symmetry = FORTSPARSE_SYM_GENERAL
    if (fixture%mtype == 4) solver%options%symmetry = FORTSPARSE_SYM_SPD
    if (fixture%mtype == -2) solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
    if (fixture%is_complex) then
        az%nrow = r%n
        az%ncol = r%n
        az%nnz = r%nnz
        az%col_ptr = fixture%ptr
        az%row_idx = fixture%row
        az%val = fixture%val
    else
        ar%nrow = r%n
        ar%ncol = r%n
        ar%nnz = r%nnz
        ar%col_ptr = fixture%ptr
        ar%row_idx = fixture%row
        ar%val = real(fixture%val, dp)
        allocate (br(r%n, r%nrhs), xr(r%n, r%nrhs))
        br = real(bz, dp)
    end if
    start = wall_time()
    call factor()
    r%factor_s = wall_time() - start
    r%status_code = status%code
    if (status%code /= 0) goto 900
    if (trim(r%workload) == 'warm') then
        r%prime_factor_s = r%factor_s
        call sparse_free(solver)
        start = wall_time()
        call factor()
        r%factor_s = wall_time() - start
        r%status_code = status%code
        if (status%code /= 0) goto 900
    end if
    call sparse_factor_stats(solver, stats, stats_status)
    r%analysis_s = stats%analysis_s
    r%numeric_s = stats%factor_s
    if (trim(r%workload) == 'refactor') then
        fixture%val = 1.125_dp*fixture%val
        bz = 1.125_dp*bz
        if (fixture%is_complex) then
            az%val = fixture%val
            start = wall_time()
            call sparse_refactor(solver, az, status, r%reused)
        else
            ar%val = real(fixture%val, dp)
            br = real(bz, dp)
            start = wall_time()
            call sparse_refactor(solver, ar, status, r%reused)
        end if
        r%refactor_s = wall_time() - start
        r%status_code = status%code
        if (status%code /= 0) goto 900
    end if
    start = wall_time()
    if (fixture%is_complex) then
        call sparse_solve_block_op(solver, bz, xz, operation, status)
    else
        call sparse_solve_block_op(solver, br, xr, operation, status)
        xz = cmplx(xr, 0.0_dp, dp)
    end if
    r%solve_s = wall_time() - start
    r%status_code = status%code
    if (status%code /= 0) goto 900
    r%backward_error = 0.0_dp
    do j = 1, r%nrhs
        r%backward_error = max(r%backward_error, &
                         fixture_backward_error(fixture, xz(:, j), bz(:, j), operation))
    end do
    r%forward_error = maxval(abs(xz - truth))/max(1.0_dp, maxval(abs(truth)))
    if (.not. (r%backward_error < 1.e-12_dp)) r%status_code = 5
    if (.not. (r%forward_error < 1.e-10_dp)) r%status_code = 5
    call export_real_solution(trim(r%family), xz)
900 continue
    call sparse_factor_stats(solver, stats, stats_status)
    r%factor_nnz = stats%factor_nnz
    r%factor_mem_mb = stats%factor_mem_mb
    r%host_current_bytes = stats%host_current_bytes
    r%host_peak_bytes = stats%host_peak_bytes
    r%active_threads = stats%active_threads
    r%native_factor_family = stats%native_factor_family
    if (index(trim(r%backend),'native') == 1 .and. stats%active_threads == 1) &
        r%native_factor_serial = 1
    call sparse_free(solver)
    r%total_s = wall_time() - total_start
    call cpu_time(cpu_end)
    r%parent_main_cpu_s = cpu_end - cpu_start
    call process_metrics(r%parent_rss_kib, r%helper_rss_kib, r%helper_cpu_s)
    r%parent_cpu_s = parent_cpu_seconds()
    call r%write_csv()
    if (r%status_code /= 0) stop 2
contains
    subroutine factor()
        if (fixture%is_complex) then
            call sparse_factor(solver, az, status)
        else
            call sparse_factor(solver, ar, status)
        end if
    end subroutine
end program bench_native_families
