program bench_csc_matvec
    ! Strong-scaling benchmark of csc_matvec (real and complex): 3D 7-point
    ! Laplacian pattern with random values. Run with OMP_NUM_THREADS=t.
    ! Also times a row-partition variant (each thread owns a row range and
    ! bisects the sorted rows of every column) for comparison.
    use fortsparse, only: dp, csc_t, csc_z_t, csc_from_triplet, csc_matvec, &
        fortsparse_status_t
    !$ use omp_lib, only: omp_get_max_threads, omp_get_wtime, omp_get_thread_num
    implicit none
    integer :: m, n, nt, rep, nreps
    character(len=32) :: arg
    type(csc_t) :: A
    type(csc_z_t) :: Z
    real(dp), allocatable :: x(:), y(:), y2(:)
    complex(dp), allocatable :: zx(:), zy(:)
    real(dp) :: t0, t1, tr, tz, tp
    integer :: i
    real(dp) :: err

    m = 100
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg); read (arg, *) m
    end if
    n = m**3
    call build(m, A, Z)
    allocate (x(n), y(n), y2(n), zx(n), zy(n))
    call random_number(x)
    zx = cmplx(x, x(n:1:-1), dp)
    nt = 1
    !$ nt = omp_get_max_threads()
    nreps = max(5, 2000000*5/A%nnz*10)
    y = csc_matvec(A, x)
    t0 = wtime()
    do rep = 1, nreps
        y = csc_matvec(A, x)
    end do
    tr = (wtime() - t0)/nreps
    zy = csc_matvec(Z, zx)
    t0 = wtime()
    do rep = 1, nreps
        zy = csc_matvec(Z, zx)
    end do
    tz = (wtime() - t0)/nreps
    call rowpart(A, x, y2)
    t0 = wtime()
    do rep = 1, nreps
        call rowpart(A, x, y2)
    end do
    tp = (wtime() - t0)/nreps
    err = maxval(abs(y - y2))/maxval(abs(y))
    write (*, "(a,i0,a,i0,a,i0,a,3(es11.3),a,es9.2)") "n=", n, " nnz=", A%nnz, &
        " threads=", nt, " t_real,t_cplx,t_rowpart[s]=", tr, tz, tp, " relerr=", err
    i = 0

contains

    function wtime() result(t)
        real(dp) :: t
        integer(8) :: c, r
        t = 0.0_dp
        !$ t = omp_get_wtime()
        !$ return
        call system_clock(c, r)
        t = real(c, dp)/real(r, dp)
    end function wtime

    subroutine rowpart(M, x, y)
        type(csc_t), intent(in) :: M
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer :: t, nth, lo, hi, j, p, q, a0, b0, mid

        nth = 1
        !$ nth = omp_get_max_threads()
        !$omp parallel num_threads(nth) private(t, lo, hi, j, p, q, a0, b0, mid)
        t = 0
        !$ t = omp_get_thread_num()
        lo = (int(M%nrow, 8)*t)/nth + 1
        hi = (int(M%nrow, 8)*(t + 1))/nth
        y(lo:hi) = 0.0_dp
        do j = 1, M%ncol
            a0 = M%col_ptr(j); b0 = M%col_ptr(j + 1)
            do while (a0 < b0)
                mid = (a0 + b0)/2
                if (M%row_idx(mid) < lo) then
                    a0 = mid + 1
                else
                    b0 = mid
                end if
            end do
            do p = a0, M%col_ptr(j + 1) - 1
                if (M%row_idx(p) > hi) exit
                y(M%row_idx(p)) = y(M%row_idx(p)) + M%val(p)*x(j)
            end do
        end do
        !$omp end parallel
    end subroutine rowpart

    subroutine build(m, A, Z)
        integer, intent(in) :: m
        type(csc_t), intent(out) :: A
        type(csc_z_t), intent(out) :: Z
        integer :: ix, iy, iz, k, c, nn, d, dx(7), dy(7), dz(7), jx, jy, jz
        integer, allocatable :: r(:), cc(:)
        real(dp), allocatable :: v(:)
        type(fortsparse_status_t) :: st
        real(dp), allocatable :: vi(:)

        dx = [0, 1, -1, 0, 0, 0, 0]; dy = [0, 0, 0, 1, -1, 0, 0]
        dz = [0, 0, 0, 0, 0, 1, -1]
        nn = m**3
        allocate (r(7*nn), cc(7*nn), v(7*nn))
        k = 0
        do iz = 1, m; do iy = 1, m; do ix = 1, m
            c = ix + m*(iy - 1) + m*m*(iz - 1)
            do d = 1, 7
                jx = ix + dx(d); jy = iy + dy(d); jz = iz + dz(d)
                if (jx < 1 .or. jx > m .or. jy < 1 .or. jy > m) cycle
                if (jz < 1 .or. jz > m) cycle
                k = k + 1
                r(k) = jx + m*(jy - 1) + m*m*(jz - 1)
                cc(k) = c
            end do
        end do; end do; end do
        allocate (vi(k))
        call random_number(v(1:k)); call random_number(vi)
        call csc_from_triplet(nn, nn, r(1:k), cc(1:k), v(1:k), A, st)
        call csc_from_triplet(nn, nn, r(1:k), cc(1:k), cmplx(v(1:k), vi, dp), Z, st)
    end subroutine build

end program bench_csc_matvec
