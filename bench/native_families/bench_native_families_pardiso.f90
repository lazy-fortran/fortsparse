program bench_native_families_pardiso
    ! Intel oneMKL PARDISO only. No FortSparse/OpenBLAS objects or interfaces.
    use, intrinsic :: iso_c_binding, only: c_int, c_int64_t, c_ptr, c_loc
    use native_family_fixtures, only: dp, family_fixture_t, &
                                      fixture_backward_error, wall_time
    use native_family_real_io, only: load_benchmark_fixture, load_benchmark_rhs, &
                                    export_real_solution
    use native_family_records, only: family_record_t, process_metrics, &
        parent_cpu_seconds
    implicit none
    interface
        subroutine pardisoinit(pt, mtype, iparm) bind(C, name='pardisoinit')
            import :: c_int, c_int64_t
            integer(c_int64_t), intent(inout) :: pt(*)
            integer(c_int), intent(in) :: mtype
            integer(c_int), intent(out) :: iparm(*)
        end subroutine
        subroutine pardiso(pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, &
                           iparm, msglvl, b, x, error) bind(C, name='pardiso')
            import :: c_int, c_int64_t, c_ptr
            integer(c_int64_t), intent(inout) :: pt(*)
            integer(c_int), intent(in) :: maxfct, mnum, mtype, phase, n, nrhs, msglvl
            type(c_ptr), value :: a, b, x
            integer(c_int), intent(in) :: ia(*), ja(*)
            integer(c_int), intent(inout) :: perm(*), iparm(*)
            integer(c_int), intent(out) :: error
        end subroutine
        integer(c_int) function set_threads(n) bind(C, name='MKL_Set_Num_Threads_Local')
            import :: c_int
            integer(c_int), value :: n
        end function
        subroutine set_dynamic(n) bind(C, name='MKL_Set_Dynamic')
            import :: c_int
            integer(c_int), value :: n
        end subroutine
    end interface
    type(family_fixture_t) :: a
    type(family_record_t) :: r
    integer(c_int64_t) :: pt(64)
    integer(c_int) :: iparm(64), code, mtype, ignored, n, nrhs
    integer(c_int), allocatable :: ia(:), ja(:), perm(:), cursor(:)
    complex(dp), allocatable, target :: avz(:), bz(:, :), xz(:, :)
    real(dp), allocatable, target :: avr(:), br(:, :), xr(:, :)
    complex(dp), allocatable :: truth(:, :)
    character(64) :: argument
    character(1024) :: mini5_path
    integer :: operation, j, p, i, pos, nnz
    real(dp) :: total_start, cpu_start, cpu_end, start
    logical :: heldout, initialized

    total_start = wall_time()
    call cpu_time(cpu_start)
    r%backend = 'intel-pardiso'
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
    n = int(r%n, c_int)
    nrhs = int(r%nrhs, c_int)
    mtype = int(a%mtype, c_int)
    initialized = .false.
    allocate (ia(n + 1), perm(n), cursor(n))
    ia = 0
    nnz = 0
    do j = 1, n
        do p = a%ptr(j), a%ptr(j + 1) - 1
            i = a%row(p)
            if (mtype == 4 .or. mtype == -2) then
                if (i > j) cycle
            end if
            nnz = nnz + 1
            ia(i + 1) = ia(i + 1) + 1
        end do
    end do
    ia(1) = 1
    do i = 1, n
        ia(i + 1) = ia(i + 1) + ia(i)
    end do
    cursor = ia(:n)
    allocate (ja(nnz), truth(n, nrhs), bz(n, nrhs), xz(n, nrhs))
    if (a%is_complex) then
        allocate (avz(nnz))
    else
        allocate (avr(nnz), br(n, nrhs), xr(n, nrhs))
    end if
    ! True CSR of original A, with sorted columns. Hermitian/symmetric uses
    ! the upper triangle; no mistaken CSC-as-CSR transpose shortcut.
    do j = 1, n
        do p = a%ptr(j), a%ptr(j + 1) - 1
            i = a%row(p)
            if (mtype == 4 .or. mtype == -2) then
                if (i > j) cycle
            end if
            pos = cursor(i)
            ja(pos) = int(j, c_int)
            if (a%is_complex) then
                avz(pos) = a%val(p)
            else
                avr(pos) = real(a%val(p), dp)
            end if
            cursor(i) = pos + 1
        end do
    end do
    call load_benchmark_rhs(a, trim(r%family), r%nrhs, heldout, operation, &
                            bz, truth, trim(mini5_path))
    if (.not. a%is_complex) br = real(bz, dp)
    perm = 0
    pt = 0_c_int64_t
    ignored = set_threads(int(r%threads, c_int))
    call set_dynamic(0_c_int)
    call pardisoinit(pt, mtype, iparm)
    iparm(2) = 2_c_int ! METIS; default numeric stability settings retained.
    iparm(18) = -1_c_int
    iparm(19) = -1_c_int
    initialized = .true.
    start = wall_time()
    call phase(11)
    r%analysis_s = wall_time() - start
    if (code /= 0) goto 900
    start = wall_time()
    call phase(22)
    r%numeric_s = wall_time() - start
    r%factor_s = r%analysis_s + r%numeric_s
    if (code /= 0) goto 900
    if (trim(r%workload) == 'warm') then
        r%prime_factor_s = r%factor_s
        call phase(-1)
        pt = 0_c_int64_t
        call pardisoinit(pt, mtype, iparm)
        iparm(2) = 2_c_int
        iparm(18) = -1_c_int
        iparm(19) = -1_c_int
        start = wall_time()
        call phase(11)
        r%analysis_s = wall_time() - start
        if (code /= 0) goto 900
        start = wall_time()
        call phase(22)
        r%numeric_s = wall_time() - start
        r%factor_s = r%analysis_s + r%numeric_s
        if (code /= 0) goto 900
    end if
    if (trim(r%workload) == 'refactor') then
        a%val = 1.125_dp*a%val
        bz = 1.125_dp*bz
        if (a%is_complex) then
            avz = 1.125_dp*avz
        else
            avr = 1.125_dp*avr
            br = real(bz, dp)
        end if
        start = wall_time()
        call phase(22)
        r%refactor_s = wall_time() - start
        if (code /= 0) goto 900
        r%reused = .true.
    end if
    start = wall_time()
    if (mtype == 4 .and. operation == 1) bz = conjg(bz)
    iparm(12) = 0_c_int
    if (mtype == 11 .or. mtype == 13) then
        if (operation == 1) iparm(12) = 2_c_int
        if (operation == 2) iparm(12) = 1_c_int
    end if
    call phase(33)
    if (mtype == 4 .and. operation == 1) then
        bz = conjg(bz)
        xz = conjg(xz)
    end if
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
    ! Vendor-reported factor memory differs from owned payload and process RSS.
    r%factor_nnz = int(iparm(18), kind(r%factor_nnz))
    r%factor_mem_mb = real(max(iparm(15), iparm(16) + iparm(17)), dp)/1024
    if (initialized) call phase(-1)
    r%total_s = wall_time() - total_start
    call cpu_time(cpu_end)
    r%parent_main_cpu_s = cpu_end - cpu_start
    call process_metrics(r%parent_rss_kib, r%helper_rss_kib, r%helper_cpu_s)
    r%parent_cpu_s = parent_cpu_seconds()
    call r%write_csv()
    if (r%status_code /= 0) stop 2
contains
    subroutine phase(action)
        integer, intent(in) :: action
        if (a%is_complex) then
            call pardiso(pt, 1_c_int, 1_c_int, mtype, int(action, c_int), n, &
             c_loc(avz), ia, ja, perm, nrhs, iparm, 0_c_int, c_loc(bz), c_loc(xz), code)
        else
            call pardiso(pt, 1_c_int, 1_c_int, mtype, int(action, c_int), n, &
             c_loc(avr), ia, ja, perm, nrhs, iparm, 0_c_int, c_loc(br), c_loc(xr), code)
        end if
    end subroutine
end program bench_native_families_pardiso
