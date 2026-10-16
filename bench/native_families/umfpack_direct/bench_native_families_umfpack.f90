! SPDX-License-Identifier: GPL-2.0-or-later
! Separate benchmark executable linking GPL UMFPACK; never part of FortSparse.
program bench_native_families_umfpack
    use, intrinsic :: iso_c_binding, only: c_int, c_int32_t, c_int64_t, &
                      c_double, c_double_complex, c_ptr, c_null_ptr, c_associated, c_loc
    use native_family_fixtures, only: dp, family_fixture_t, &
                                      fixture_backward_error, wall_time
    use native_family_real_io, only: load_benchmark_fixture, load_benchmark_rhs, &
                                    export_real_solution
    use native_family_records, only: family_record_t, process_metrics, &
                                     parent_cpu_seconds
    implicit none
    interface
        function umf_create(is_complex) result(owner) bind(C, name='fsbench_umf_create')
            import :: c_int, c_ptr
            integer(c_int), value :: is_complex
            type(c_ptr) :: owner
        end function
        integer(c_int) function umf_symbolic(owner, n, ap, ai, ax) &
            bind(C, name='fsbench_umf_symbolic')
            import :: c_int, c_int32_t, c_ptr
            type(c_ptr), value :: owner, ax
            integer(c_int32_t), value :: n
            integer(c_int32_t), intent(in) :: ap(*), ai(*)
        end function
        integer(c_int) function umf_numeric(owner, ap, ai, ax) &
            bind(C, name='fsbench_umf_numeric')
            import :: c_int, c_int32_t, c_ptr
            type(c_ptr), value :: owner, ax
            integer(c_int32_t), intent(in) :: ap(*), ai(*)
        end function
        integer(c_int) function umf_solve(owner, operation, ap, ai, ax, b, x) &
            bind(C, name='fsbench_umf_solve')
            import :: c_int, c_int32_t, c_ptr
            type(c_ptr), value :: owner, ax, b, x
            integer(c_int), value :: operation
            integer(c_int32_t), intent(in) :: ap(*), ai(*)
        end function
        subroutine umf_stats(owner, nnz, peak_mb) bind(C, name='fsbench_umf_stats')
            import :: c_ptr, c_int64_t, c_double
            type(c_ptr), value :: owner
            integer(c_int64_t), intent(out) :: nnz
            real(c_double), intent(out) :: peak_mb
        end subroutine
        subroutine umf_destroy(owner) bind(C, name='fsbench_umf_destroy')
            import :: c_ptr
            type(c_ptr), value :: owner
        end subroutine
    end interface
    type(family_fixture_t) :: a
    type(family_record_t) :: r
    type(c_ptr) :: owner = c_null_ptr, ax
    integer(c_int32_t), allocatable :: ap(:), ai(:)
    complex(c_double_complex), allocatable, target :: az(:), bz(:, :), xz(:, :)
    real(c_double), allocatable, target :: ar(:), br(:, :), xr(:, :)
    complex(dp), allocatable :: truth(:, :)
    character(64) :: argument
    character(1024) :: mini5_path
    integer :: operation, j
    integer(c_int) :: code
    real(dp) :: total_start, cpu_start, cpu_end, start
    logical :: heldout

    total_start = wall_time()
    call cpu_time(cpu_start)
    if (dp /= c_double .or. dp /= c_double_complex) &
        error stop 'fixture real64 and C double complex kinds differ'
    r%backend = 'umfpack-direct'
    call get_command_argument(1, r%family)
    call get_command_argument(2, argument)
    read (argument, *) r%level
    call get_command_argument(3, argument)
    read (argument, *) r%threads
    call get_command_argument(4, r%workload)
    call get_command_argument(5, r%operation)
    call get_command_argument(6, r%variant)
    call get_command_argument(7, mini5_path)
    if (r%threads < 1) error stop 'positive thread budget required'
    heldout = trim(r%variant) == 'heldout'
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
    call load_benchmark_fixture(trim(r%family), r%level, heldout, a, trim(mini5_path))
    r%n = a%n
    r%nnz = size(a%val)
    allocate (ap(a%n + 1), ai(r%nnz), bz(a%n, r%nrhs), xz(a%n, r%nrhs), &
              truth(a%n, r%nrhs))
    ap = int(a%ptr - 1, c_int32_t)
    ai = int(a%row - 1, c_int32_t)
    if (a%is_complex) then
        allocate (az(r%nnz))
        az = a%val
        ax = c_loc(az)
    else
        allocate (ar(r%nnz), br(a%n, r%nrhs), xr(a%n, r%nrhs))
        ar = real(a%val, dp)
        ax = c_loc(ar)
    end if
    call load_benchmark_rhs(a, trim(r%family), r%nrhs, heldout, operation, &
                            bz, truth, trim(mini5_path))
    if (.not. a%is_complex) br = real(bz, dp)
    owner = umf_create(int(merge(1, 0, a%is_complex), c_int))
    code = 6_c_int
    if (.not. c_associated(owner)) goto 900
    call factor()
    if (code /= 0) goto 900
    if (trim(r%workload) == 'warm') then
        r%prime_factor_s = r%factor_s
        call umf_destroy(owner)
        owner = umf_create(int(merge(1, 0, a%is_complex), c_int))
        code = 6_c_int
        if (.not. c_associated(owner)) goto 900
        call factor()
        if (code /= 0) goto 900
    end if
    if (trim(r%workload) == 'refactor') then
        a%val = 1.125_dp*a%val
        bz = 1.125_dp*bz
        if (a%is_complex) then
            az = 1.125_dp*az
        else
            ar = 1.125_dp*ar
            br = real(bz, dp)
        end if
        start = wall_time()
        code = umf_numeric(owner, ap, ai, ax)
        r%refactor_s = wall_time() - start
        if (code /= 0) goto 900
        r%reused = .true.
    end if
    start = wall_time()
    do j = 1, r%nrhs
        if (a%is_complex) then
            code = umf_solve(owner, int(operation, c_int), ap, ai, ax, &
                             c_loc(bz(1, j)), c_loc(xz(1, j)))
        else
            code = umf_solve(owner, int(operation, c_int), ap, ai, ax, &
                             c_loc(br(1, j)), c_loc(xr(1, j)))
        end if
        if (code /= 0) exit
    end do
    r%solve_s = wall_time() - start
    if (code /= 0) goto 900
    if (.not. a%is_complex) xz = cmplx(xr, 0.0_dp, dp)
    r%backward_error = 0.0_dp
    do j = 1, r%nrhs
        r%backward_error = max(r%backward_error, &
                               fixture_backward_error(a, xz(:, j), bz(:, j), operation))
    end do
    r%forward_error = maxval(abs(xz - truth))/max(1.0_dp, maxval(abs(truth)))
    if (.not. (r%backward_error < 1.e-12_dp)) code = 5
    if (.not. (r%forward_error < 1.e-10_dp)) code = 5
    call export_real_solution(trim(r%family), xz)
900 continue
    r%status_code = int(code)
    if (c_associated(owner)) then
        call umf_stats(owner, r%factor_nnz, r%factor_mem_mb)
        call umf_destroy(owner)
    end if
    r%total_s = wall_time() - total_start
    call cpu_time(cpu_end)
    r%parent_main_cpu_s = cpu_end - cpu_start
    call process_metrics(r%parent_rss_kib, r%helper_rss_kib, r%helper_cpu_s)
    r%parent_cpu_s = parent_cpu_seconds()
    call r%write_csv()
    if (r%status_code /= 0) stop 2
contains
    subroutine factor()
        start = wall_time()
        code = umf_symbolic(owner, int(a%n, c_int32_t), ap, ai, ax)
        r%analysis_s = wall_time() - start
        if (code /= 0) return
        start = wall_time()
        code = umf_numeric(owner, ap, ai, ax)
        r%numeric_s = wall_time() - start
        r%factor_s = r%analysis_s + r%numeric_s
    end subroutine factor
end program bench_native_families_umfpack
