program test_mpm_allocations
    ! External malloc/calloc/realloc oracle: weighted matching and scaling run
    ! entirely in caller-owned storage, for real and complex, cached and
    ! recomputed costs, including the singular early-return paths.
    use fortsparse, only: dp, fortsparse_status_t, status_ok, FORTSPARSE_SINGULAR
    use fortsparse_native_mpm, only: mpm_stats_t, mpm_iwork_size, &
                                     mpm_rwork_size, mpm_workspace_bytes, &
                                     max_product_matching, mpm_scaling
    use, intrinsic :: iso_c_binding, only: c_int64_t
    implicit none
    integer, parameter :: n = 300, per = 5
    integer :: ptr(n + 1), rows(per*n), cmatch(n), j, k, e
    integer, allocatable :: iwork(:)
    real(dp) :: values(per*n), u(n), v(n), dr(n), dc(n)
    real(dp), allocatable :: rcache(:), rplain(:)
    complex(dp) :: z(per*n)
    type(fortsparse_status_t) :: status
    type(mpm_stats_t) :: stats
    integer(c_int64_t) :: count, current, maximum
    integer, parameter :: hall_ptr(4) = [1, 2, 3, 6], hall_rows(5) = [1, 1, 1, 2, 3]
    real(dp), parameter :: hall_values(5) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    interface
        subroutine begin(nth) bind(C, name='fsparse_spd_probe_begin')
            import c_int64_t
            integer(c_int64_t), value :: nth
        end subroutine
        subroutine finish(count, current, maximum) &
            bind(C, name='fsparse_spd_probe_end')
            import c_int64_t
            integer(c_int64_t), intent(out) :: count, current, maximum
        end subroutine
    end interface
    ! Sorted strided pattern containing the anti-diagonal transversal.
    e = 0; ptr(1) = 1
    do j = 1, n
        do k = 0, per - 1
            e = e + 1
            rows(e) = 1 + mod(n - j + k*7, n)
            values(e) = real(1 + mod(13*e, 17), dp)*10.0_dp**mod(e, 9)
        end do
        call sort(rows(ptr(j):e), values(ptr(j):e))
        ptr(j + 1) = e + 1
    end do
    z = cmplx(values, -0.5_dp*values, dp)
    allocate (iwork(mpm_iwork_size(n)), rcache(mpm_rwork_size(n, e, .true.)), &
              rplain(mpm_rwork_size(n, e, .false.)))
    if (4*size(iwork) + 8*size(rcache) /= mpm_workspace_bytes(n, e, .true.)) &
        error stop 'cached workspace identity'
    if (4*size(iwork) + 8*size(rplain) /= mpm_workspace_bytes(n, e, .false.)) &
        error stop 'plain workspace identity'
    call begin(0_c_int64_t)
    call max_product_matching(n, ptr, rows, values, cmatch, u, v, iwork, rcache, &
                              status, stats)
    if (.not. status_ok(status)) error stop 'real cached matching'
    call mpm_scaling(n, ptr, rows, values, u, v, dr, dc, status)
    if (.not. status_ok(status)) error stop 'real scaling'
    call max_product_matching(n, ptr, rows, z, cmatch, u, v, iwork, rplain, &
                              status, stats)
    if (.not. status_ok(status)) error stop 'complex plain matching'
    call mpm_scaling(n, ptr, rows, z, u, v, dr, dc, status)
    if (.not. status_ok(status)) error stop 'complex scaling'
    values(ptr(7):ptr(8) - 1) = 0
    call max_product_matching(n, ptr, rows, values, cmatch, u, v, iwork, rcache, &
                              status, stats)
    if (status%code /= FORTSPARSE_SINGULAR) error stop 'zero column singular'
    ! Empty-row and failed-augmentation messages are formatted without I/O.
    call max_product_matching(3, hall_ptr, hall_rows, hall_values, cmatch(1:3), &
                              u(1:3), v(1:3), iwork(1:18), rplain(1:6), status, &
                              stats)
    if (status%code /= FORTSPARSE_SINGULAR) error stop 'Hall singular'
    call max_product_matching(3, [1, 3, 5, 7], [1, 3, 1, 3, 1, 3], &
                              [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], &
                              cmatch(1:3), u(1:3), v(1:3), iwork(1:18), &
                              rplain(1:6), status, stats)
    if (status%code /= FORTSPARSE_SINGULAR) error stop 'empty row singular'
    call finish(count, current, maximum)
    if (count /= 0 .or. maximum /= 0) then
        print *, count, current, maximum
        error stop 'mpm allocated'
    end if
    if (stats%greedy_matched < 0) error stop 'stats'
contains
    subroutine sort(r, x)
        integer, intent(inout) :: r(:)
        real(dp), intent(inout) :: x(:)
        integer :: i, m, t
        real(dp) :: y
        do i = 2, size(r)
            t = r(i); y = x(i); m = i - 1
            do while (m >= 1)
                if (r(m) <= t) exit
                r(m + 1) = r(m); x(m + 1) = x(m); m = m - 1
            end do
            r(m + 1) = t; x(m + 1) = y
        end do
    end subroutine
end program
