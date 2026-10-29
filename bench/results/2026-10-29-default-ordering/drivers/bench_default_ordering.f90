program bench_default_ordering
    ! One fresh process per (case, ordering): build, factor, solve, residual.
    ! Usage: bench_default_ordering laplace|elastic natural|default size
    use, intrinsic :: iso_c_binding, only: c_double
    use fortsparse
    implicit none
    interface
        subroutine process_metrics(rss, child_rss, child_cpu) bind(C, name='fsbench_process_metrics')
            import :: c_double
            real(c_double), intent(out) :: rss, child_rss, child_cpu
        end subroutine
        real(c_double) function parent_cpu() bind(C, name='fsbench_parent_cpu_seconds')
            import :: c_double
        end function
    end interface
    character(16) :: kind, order_name, arg
    integer :: m, ordering, i
    integer(8) :: t0, t1, t2, rate
    type(csc_t) :: A
    type(csc_z_t) :: Z
    type(sparse_solver_t) :: solver
    type(fortsparse_status_t) :: status
    type(sparse_factor_stats_t) :: stats
    real(dp), allocatable :: b(:), x(:)
    complex(dp), allocatable :: zb(:), zx(:)
    real(dp) :: rel, rss, crss, ccpu
    integer :: n, nnz

    call get_command_argument(1, kind)
    call get_command_argument(2, order_name)
    call get_command_argument(3, arg)
    read (arg, *) m
    ordering = FORTSPARSE_ORDER_DEFAULT
    if (order_name == 'natural') ordering = FORTSPARSE_ORDER_NATURAL
    solver%backend_id = FORTSPARSE_BACKEND_NATIVE
    solver%options%ordering = ordering
    call system_clock(count_rate=rate)
    if (kind == 'laplace') then
        call laplacian(m, A)
        n = A%nrow; nnz = A%nnz
        allocate (b(n), x(n))
        b = [(sin(0.37_dp*real(i, dp)) + 1.0_dp, i=1, n)]
        call system_clock(t0)
        call sparse_factor(solver, A, status)
        call system_clock(t1)
        if (.not. status_ok(status)) error stop 'factor failed'
        call sparse_solve(solver, b, x, status)
        call system_clock(t2)
        rel = maxval(abs(csc_matvec(A, x) - b))/(8.0_dp*maxval(abs(x)) + maxval(abs(b)))
    else
        call elasticity(m, Z)
        n = Z%nrow; nnz = Z%nnz
        allocate (zb(n), zx(n))
        zb = [(cmplx(cos(0.21_dp*real(i, dp)), 1.0_dp/real(i, dp), dp), i=1, n)]
        call system_clock(t0)
        call sparse_factor(solver, Z, status)
        call system_clock(t1)
        if (.not. status_ok(status)) error stop 'factor failed'
        call sparse_solve(solver, zb, zx, status)
        call system_clock(t2)
        rel = maxval(abs(csc_matvec(Z, zx) - zb))/(20.0_dp*maxval(abs(zx)) + maxval(abs(zb)))
    end if
    if (.not. status_ok(status)) error stop 'solve failed'
    call sparse_factor_stats(solver, stats, status)
    call process_metrics(rss, crss, ccpu)
    print '(a,",",a,",",i0,",",i0,",",i0,",",i0,",",es10.3,",",es10.3,",",f10.3,",",f8.3,",",es9.2,",",f12.0,",",f10.3,",",i0)', &
        trim(kind), trim(order_name), n, nnz, stats%ordering_used, stats%factor_nnz, &
        stats%analysis_s, stats%factor_s, real(t1 - t0, dp)/real(rate, dp), &
        real(t2 - t1, dp)/real(rate, dp), rel, rss, parent_cpu(), stats%active_threads
contains
    subroutine laplacian(m, A)
        integer, intent(in) :: m
        type(csc_t), intent(out) :: A
        integer, allocatable :: r(:), c(:)
        real(dp), allocatable :: v(:)
        integer :: ix, iy, k, p
        allocate (r(5*m*m), c(5*m*m), v(5*m*m))
        k = 0
        do iy = 1, m
            do ix = 1, m
                p = (iy - 1)*m + ix
                k = k + 1; r(k) = p; c(k) = p; v(k) = 4.0_dp
                if (ix > 1) then; k = k + 1; r(k) = p; c(k) = p - 1; v(k) = -1.0_dp; end if
                if (ix < m) then; k = k + 1; r(k) = p; c(k) = p + 1; v(k) = -1.0_dp; end if
                if (iy > 1) then; k = k + 1; r(k) = p; c(k) = p - m; v(k) = -1.0_dp; end if
                if (iy < m) then; k = k + 1; r(k) = p; c(k) = p + m; v(k) = -1.0_dp; end if
            end do
        end do
        call csc_from_triplet(m*m, m*m, r(:k), c(:k), v(:k), A, status)
        if (.not. status_ok(status)) error stop 'laplacian'
    end subroutine
    subroutine elasticity(m, A)
        ! Same element as test_fortsparse_native_default_ordering.
        integer, intent(in) :: m
        type(csc_z_t), intent(out) :: A
        integer, allocatable :: r(:), c(:)
        complex(dp), allocatable :: v(:)
        real(dp), parameter :: lap(3, 3) = reshape([1.1_dp, -0.9_dp, 0.1_dp, &
            -0.9_dp, 2.1_dp, -0.9_dp, 0.1_dp, -0.9_dp, 1.1_dp], [3, 3])
        real(dp), parameter :: cp(2, 2) = reshape([2.0_dp, 0.5_dp, 0.5_dp, 2.0_dp], [2, 2])
        integer :: ix, iy, k, nn, i, j, q, bb, t, v00, vert(3)
        nn = 2*m*m
        allocate (r(72*(m - 1)**2 + nn), c(72*(m - 1)**2 + nn), v(72*(m - 1)**2 + nn))
        k = 0
        do iy = 1, m - 1
            do ix = 1, m - 1
                v00 = (iy - 1)*m + ix
                do t = 1, 2
                    if (t == 1) vert = [v00, v00 + 1, v00 + m + 1]
                    if (t == 2) vert = [v00, v00 + m, v00 + m + 1]
                    do j = 1, 3
                        do bb = 1, 2
                            do i = 1, 3
                                do q = 1, 2
                                    k = k + 1
                                    r(k) = 2*(vert(i) - 1) + q
                                    c(k) = 2*(vert(j) - 1) + bb
                                    v(k) = (1.0_dp, 0.2_dp)*lap(i, j)*cp(q, bb)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
        do i = 1, nn
            k = k + 1; r(k) = i; c(k) = i; v(k) = (0.05_dp, 0.3_dp)
        end do
        call csc_from_triplet(nn, nn, r(:k), c(:k), v(:k), A, status)
        if (.not. status_ok(status)) error stop 'elasticity'
    end subroutine
end program bench_default_ordering
