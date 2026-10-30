program bench_symmetric_ordering
    ! One fresh process per row: build, factor, solve, check against A.
    ! Usage: bench_symmetric_ordering FAMILY ORDER M SHIFT [MODE]
    !   FAMILY ldlt  real SYM_INDEFINITE, five-point Laplacian with diagonal 4-s
    !          kkt   real SYM_INDEFINITE, [[L, B^T], [B, 0]], B couples vertex pairs
    !          hpd   complex SYM_SPD, unitary gauge D L D^H of the Laplacian
    !          ldlh  complex SYM_HERMITIAN, D (L - s I) D^H
    !          zsym  complex SYM_INDEFINITE (transpose-symmetric), (1+0.2i)(L - s I)
    !          spd   real SYM_SPD Laplacian; MODE low | cap | refine
    !   ORDER  natural | amd | nd | default;  M grid side;  SHIFT s.
    use, intrinsic :: iso_c_binding, only: c_double
    use fortsparse
    implicit none
    interface
        subroutine process_metrics(rss, child_rss, child_cpu) &
                bind(C, name='fsbench_process_metrics')
            import :: c_double
            real(c_double), intent(out) :: rss, child_rss, child_cpu
        end subroutine
        real(c_double) function parent_cpu() bind(C, name='fsbench_parent_cpu_seconds')
            import :: c_double
        end function
    end interface
    character(16) :: family, order_name, mode, arg
    integer :: m, n, nnz, i, exact, reliability, code
    integer(8) :: t0, t1, t2, rate
    integer(i8) :: neg, zero, pos, plan_peak
    real(dp) :: shift, berr, anorm, rss, crss, ccpu
    type(csc_t) :: A
    type(csc_z_t) :: Z
    type(sparse_solver_t) :: solver
    type(fortsparse_status_t) :: status, fstatus
    type(sparse_factor_stats_t) :: stats
    type(sparse_memory_plan_t) :: plan
    real(dp), allocatable :: b(:), x(:)
    complex(dp), allocatable :: zb(:), zx(:)

    call get_command_argument(1, family)
    call get_command_argument(2, order_name)
    call get_command_argument(3, arg); read (arg, *) m
    call get_command_argument(4, arg); read (arg, *) shift
    mode = 'auto'
    if (command_argument_count() >= 5) call get_command_argument(5, mode)
    solver%backend_id = FORTSPARSE_BACKEND_NATIVE
    solver%options%ordering = FORTSPARSE_ORDER_DEFAULT
    if (order_name == 'natural') solver%options%ordering = FORTSPARSE_ORDER_NATURAL
    if (order_name == 'amd') solver%options%ordering = FORTSPARSE_ORDER_AMD
    if (order_name == 'nd') solver%options%ordering = FORTSPARSE_ORDER_ND
    call system_clock(count_rate=rate)
    exact = -1
    plan_peak = -1
    berr = -1.0_dp
    neg = -1; zero = -1; pos = -1
    select case (family)
    case ('ldlt', 'kkt', 'spd')
        if (family == 'kkt') then
            call saddle(m, A)
            exact = (m*m)/2
        else
            call real_laplacian(m, shift, A)
            exact = negatives(m, shift)
        end if
        solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        if (family == 'spd') then
            solver%options%symmetry = FORTSPARSE_SYM_SPD
            if (mode == 'low' .or. mode == 'refine') &
                solver%options%memory_policy = FORTSPARSE_MEMORY_LOW
            if (mode == 'cap') solver%options%host_memory_limit_bytes = 8000000000_i8
            solver%refine = mode == 'refine'
            call sparse_analyze(solver, A, plan, status)
            if (status_ok(status)) plan_peak = plan%peak_host_bytes
        end if
        n = A%nrow; nnz = A%nnz
        allocate (b(n), x(n))
        b = [(sin(0.37_dp*real(i, dp)) + 1.0_dp, i=1, n)]
        call system_clock(t0)
        call sparse_factor(solver, A, fstatus)
        call system_clock(t1)
        if (status_ok(fstatus)) then
            call sparse_solve(solver, b, x, status)
            anorm = csc_norm(A%col_ptr, abs(A%val))
            berr = maxval(abs(csc_matvec(A, x) - b))/(anorm*maxval(abs(x)) + maxval(abs(b)))
        end if
        call system_clock(t2)
    case default
        select case (family)
        case ('hpd')
            call complex_laplacian(m, 0.0_dp, .true., Z)
            solver%options%symmetry = FORTSPARSE_SYM_SPD
            exact = 0
        case ('ldlh')
            call complex_laplacian(m, shift, .true., Z)
            solver%options%symmetry = FORTSPARSE_SYM_HERMITIAN
            exact = negatives(m, shift)
        case default
            call complex_laplacian(m, shift, .false., Z)
            solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        end select
        n = Z%nrow; nnz = Z%nnz
        allocate (zb(n), zx(n))
        zb = [(cmplx(cos(0.21_dp*real(i, dp)), 1.0_dp/real(i, dp), dp), i=1, n)]
        call system_clock(t0)
        call sparse_factor(solver, Z, fstatus)
        call system_clock(t1)
        if (status_ok(fstatus)) then
            call sparse_solve(solver, zb, zx, status)
            anorm = csc_norm(Z%col_ptr, abs(Z%val))
            berr = maxval(abs(csc_matvec(Z, zx) - zb))/(anorm*maxval(abs(zx)) + maxval(abs(zb)))
        end if
        call system_clock(t2)
    end select
    call sparse_factor_stats(solver, stats, status)
    call sparse_inertia(solver, neg, zero, pos, reliability, status)
    code = fstatus%code
    call process_metrics(rss, crss, ccpu)
    print '(a,",",a,",",a,",",i0,",",i0,",",f0.4,4(",",i0),4(",",es10.3),3(",",i0),",",i0,",",f12.0,",",f10.3,2(",",i0),",",i0)', &
        trim(family), trim(order_name), trim(mode), m, n, shift, nnz, code, &
        stats%ordering_used, stats%factor_nnz, stats%analysis_s, stats%factor_s, &
        real(t1 - t0, dp)/real(rate, dp), berr, neg, zero, pos, exact, rss, parent_cpu(), &
        stats%host_peak_bytes, plan_peak, stats%active_threads
contains
    pure real(dp) function csc_norm(ptr, mag) result(norm)
        ! Max column sum = infinity norm for the symmetric/Hermitian matrices here.
        integer, intent(in) :: ptr(:)
        real(dp), intent(in) :: mag(:)
        integer :: j
        norm = 0.0_dp
        do j = 1, size(ptr) - 1
            norm = max(norm, sum(mag(ptr(j):ptr(j + 1) - 1)))
        end do
    end function
    integer function negatives(m, s) result(count)
        integer, intent(in) :: m
        real(dp), intent(in) :: s
        real(dp) :: h
        integer :: i, j
        h = acos(-1.0_dp)/real(m + 1, dp)
        count = 0
        do j = 1, m
            do i = 1, m
                if (4.0_dp - 2.0_dp*cos(i*h) - 2.0_dp*cos(j*h) < s) count = count + 1
            end do
        end do
    end function
    subroutine stencil(m, r, c, k)
        integer, intent(in) :: m
        integer, allocatable, intent(out) :: r(:), c(:)
        integer, intent(out) :: k
        integer :: ix, iy, p
        allocate (r(5*m*m), c(5*m*m))
        k = 0
        do iy = 1, m
            do ix = 1, m
                p = (iy - 1)*m + ix
                k = k + 1; r(k) = p; c(k) = p
                if (ix > 1) then; k = k + 1; r(k) = p; c(k) = p - 1; end if
                if (ix < m) then; k = k + 1; r(k) = p; c(k) = p + 1; end if
                if (iy > 1) then; k = k + 1; r(k) = p; c(k) = p - m; end if
                if (iy < m) then; k = k + 1; r(k) = p; c(k) = p + m; end if
            end do
        end do
    end subroutine
    subroutine real_laplacian(m, s, A)
        integer, intent(in) :: m
        real(dp), intent(in) :: s
        type(csc_t), intent(out) :: A
        integer, allocatable :: r(:), c(:)
        real(dp), allocatable :: v(:)
        integer :: k
        call stencil(m, r, c, k)
        allocate (v(k))
        where (r(:k) == c(:k))
            v = 4.0_dp - s
        elsewhere
            v = -1.0_dp
        end where
        call csc_from_triplet(m*m, m*m, r(:k), c(:k), v, A, status)
        if (.not. status_ok(status)) error stop 'laplacian'
    end subroutine
    subroutine complex_laplacian(m, s, hermitian, Z)
        integer, intent(in) :: m
        real(dp), intent(in) :: s
        logical, intent(in) :: hermitian
        type(csc_z_t), intent(out) :: Z
        integer, allocatable :: r(:), c(:)
        complex(dp), allocatable :: v(:)
        integer :: k, p
        call stencil(m, r, c, k)
        allocate (v(k))
        do p = 1, k
            if (hermitian) then
                v(p) = -exp(cmplx(0.0_dp, 0.3_dp*real(r(p) - c(p), dp), dp))
                if (r(p) == c(p)) v(p) = 4.0_dp - s
            else
                v(p) = -(1.0_dp, 0.2_dp)
                if (r(p) == c(p)) v(p) = (1.0_dp, 0.2_dp)*(4.0_dp - s)
            end if
        end do
        call csc_from_triplet(m*m, m*m, r(:k), c(:k), v, Z, status)
        if (.not. status_ok(status)) error stop 'complex laplacian'
    end subroutine
    subroutine saddle(m, A)
        integer, intent(in) :: m
        type(csc_t), intent(out) :: A
        integer, allocatable :: r(:), c(:), rr(:), cc(:)
        real(dp), allocatable :: v(:)
        integer :: k, nh, i
        call stencil(m, r, c, k)
        nh = m*m
        allocate (rr(k + 2*nh), cc(k + 2*nh), v(k + 2*nh))
        rr(:k) = r(:k); cc(:k) = c(:k)
        where (r(:k) == c(:k))
            v(:k) = 4.0_dp
        elsewhere
            v(:k) = -1.0_dp
        end where
        do i = 1, nh/2
            rr(k + 1:k + 4) = [nh + i, 2*i - 1, nh + i, 2*i]
            cc(k + 1:k + 4) = [2*i - 1, nh + i, 2*i, nh + i]
            v(k + 1:k + 4) = [1.0_dp, 1.0_dp, -2.0_dp, -2.0_dp]
            k = k + 4
        end do
        call csc_from_triplet(nh + nh/2, nh + nh/2, rr(:k), cc(:k), v(:k), A, status)
        if (.not. status_ok(status)) error stop 'saddle'
    end subroutine
end program bench_symmetric_ordering
