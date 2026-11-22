program test_fortsparse_symmetric_triangle
    ! One-triangle symmetric input (sparse_options_t%triangle). Oracles: the
    ! solution and inertia of the same matrix given in full storage, the exact
    ! inertia of a saddle-point matrix [[H, B^T], [B, 0]] (n_H positive, m_B
    ! negative eigenvalues) whose zero block forces 2x2 pivots, independence
    ! of the unreferenced triangle (a full matrix that is symmetric only to
    ! round-off), retained-analysis refactor, and the error paths.
    use fortsparse
    implicit none
    integer, parameter :: nh = 60, mb = 17, n = nh + mb
    type(csc_t) :: full, lower, upper, skew, shifted_full, shifted_lower
    type(sparse_solver_t) :: s_full, s_low, s_up, s_ref
    type(fortsparse_status_t) :: status
    real(dp) :: b(n), x_full(n), x_low(n), x_up(n), x_skew(n), x_new(n), x_ref(n)
    integer(i8) :: neg, zero, pos
    integer :: failures, reliability
    logical :: reused
    integer :: k, nl, nf
    integer :: lr(1000), lc(1000), fr(1000), fc(1000)
    real(dp) :: lv(1000), fv(1000)

    failures = 0
    call saddle_point(1.0_dp, full, lower, upper)
    b = [(sin(0.7_dp*real(k, dp)) + 0.25_dp, k=1, n)]

    ! Reference: full storage, exactly symmetric.
    call configure(s_full, FORTSPARSE_TRIANGLE_FULL)
    call sparse_factor(s_full, full, status)
    call expect(status_ok(status), 'full factor')
    call sparse_solve(s_full, b, x_full, status)
    call expect(status_ok(status), 'full solve')
    call sparse_inertia(s_full, neg, zero, pos, reliability, status)
    call expect(status_ok(status), 'full inertia available')
    call expect(neg == mb .and. zero == 0 .and. pos == nh, 'full inertia is (m, 0, n)')

    ! Lower triangle only.
    call configure(s_low, FORTSPARSE_TRIANGLE_LOWER)
    call sparse_factor(s_low, lower, status)
    call expect(status_ok(status), 'lower-triangle factor')
    call sparse_solve(s_low, b, x_low, status)
    call expect(status_ok(status), 'lower-triangle solve')
    call expect(all(x_low == x_full), 'lower triangle gives the full-storage solution')
    call sparse_inertia(s_low, neg, zero, pos, reliability, status)
    call expect(status_ok(status), 'lower-triangle inertia available')
    call expect(neg == mb .and. zero == 0 .and. pos == nh, 'lower inertia is (m, 0, n)')
    call expect(residual(full, x_low, b) < 1.0e-12_dp, 'lower-triangle residual')

    ! Upper triangle only.
    call configure(s_up, FORTSPARSE_TRIANGLE_UPPER)
    call sparse_factor(s_up, upper, status)
    call expect(status_ok(status), 'upper-triangle factor')
    call sparse_solve(s_up, b, x_up, status)
    call expect(all(x_up == x_full), 'upper triangle gives the full-storage solution')
    call sparse_inertia(s_up, neg, zero, pos, reliability, status)
    call expect(neg == mb .and. zero == 0 .and. pos == nh, 'upper inertia is (m, 0, n)')

    ! A full matrix symmetric only to round-off: rejected as full storage,
    ! accepted when its lower triangle is selected (the upper one is not read).
    skew = full
    call perturb_upper(skew)
    call configure(s_ref, FORTSPARSE_TRIANGLE_FULL)
    call sparse_factor(s_ref, skew, status)
    call expect(.not. status_ok(status), 'full storage still requires exact symmetry')
    call sparse_destroy(s_ref)
    call configure(s_ref, FORTSPARSE_TRIANGLE_LOWER)
    call sparse_factor(s_ref, skew, status)
    call expect(status_ok(status), 'lower triangle of a round-off asymmetric matrix')
    call sparse_solve(s_ref, b, x_skew, status)
    call expect(all(x_skew == x_full), 'unreferenced triangle does not matter')
    call sparse_destroy(s_ref)

    ! Refactor with new values: analysis kept, result equal to a fresh factor.
    call saddle_point(1.6_dp, shifted_full, shifted_lower, upper)
    call sparse_refactor(s_low, shifted_lower, status, reused)
    call expect(status_ok(status), 'lower-triangle refactor')
    call expect(reused, 'lower-triangle refactor keeps the analysis')
    call sparse_solve(s_low, b, x_new, status)
    call sparse_inertia(s_low, neg, zero, pos, reliability, status)
    call expect(status_ok(status), 'refactor inertia available')
    call expect(neg == mb .and. zero == 0 .and. pos == nh, 'refactor inertia')
    call configure(s_ref, FORTSPARSE_TRIANGLE_FULL)
    call sparse_factor(s_ref, shifted_full, status)
    call sparse_solve(s_ref, b, x_ref, status)
    call expect(all(x_new == x_ref), 'refactor equals fresh full-storage factor')
    call expect(residual(shifted_full, x_new, b) < 1.0e-12_dp, 'refactor residual')
    call sparse_destroy(s_ref)

    ! Scalar Bunch--Kaufman reference path with one triangle.
    call configure(s_ref, FORTSPARSE_TRIANGLE_LOWER)
    s_ref%backend_id = FORTSPARSE_BACKEND_NATIVE_REFERENCE
    s_ref%options%ordering = FORTSPARSE_ORDER_NATURAL
    call sparse_factor(s_ref, lower, status)
    call expect(status_ok(status), 'reference lower-triangle factor')
    if (.not. status_ok(status)) print '(a)', trim(status%msg)
    if (status_ok(status)) then
        call sparse_solve(s_ref, b, x_ref, status)
        call expect(residual(full, x_ref, b) < 1.0e-10_dp, 'reference residual')
        call sparse_inertia(s_ref, neg, zero, pos, reliability, status)
        call expect(status_ok(status) .and. neg == mb .and. pos == nh, &
            'reference lower-triangle inertia')
    end if
    call sparse_destroy(s_ref)

    ! Symmetric positive definite, lower triangle.
    call spd_check()

    ! Raw CSC arrays, one triangle.
    call configure(s_ref, FORTSPARSE_TRIANGLE_LOWER)
    call sparse_factor(s_ref, n, n, lower%nnz, lower%col_ptr, lower%row_idx, &
        lower%val, status)
    call expect(status_ok(status), 'raw lower-triangle factor')
    call sparse_solve(s_ref, b, x_ref, status)
    call expect(all(x_ref == x_full), 'raw lower-triangle solution')
    call sparse_refactor(s_ref, n, n, shifted_lower%nnz, shifted_lower%col_ptr, &
        shifted_lower%row_idx, shifted_lower%val, status, reused)
    call expect(status_ok(status) .and. reused, 'raw lower-triangle refactor reuse')
    call sparse_destroy(s_ref)

    ! A triangle needs a symmetric structure hint.
    call configure(s_ref, FORTSPARSE_TRIANGLE_LOWER)
    s_ref%options%symmetry = FORTSPARSE_SYM_GENERAL
    call sparse_factor(s_ref, lower, status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX, &
        'GENERAL with a triangle is rejected')
    call sparse_destroy(s_ref)

    call sparse_destroy(s_full)
    call sparse_destroy(s_low)
    call sparse_destroy(s_up)
    if (failures /= 0) then
        print '(a,i0,a)', 'test_fortsparse_symmetric_triangle: ', failures, ' failures'
        error stop 1
    end if
    print '(a)', 'test_fortsparse_symmetric_triangle: ok'

contains

    subroutine configure(s, triangle)
        type(sparse_solver_t), intent(inout) :: s
        integer, intent(in) :: triangle
        s%backend_id = FORTSPARSE_BACKEND_NATIVE
        s%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        s%options%ordering = FORTSPARSE_ORDER_AMD
        s%options%nthreads = 1
        s%options%triangle = triangle
    end subroutine configure

    subroutine expect(condition, what)
        logical, intent(in) :: condition
        character(*), intent(in) :: what
        if (condition) return
        failures = failures + 1
        print '(2a)', 'FAIL: ', what
    end subroutine expect

    real(dp) function residual(m, x, rhs) result(berr)
        ! Normwise backward error of x for the full matrix m.
        type(csc_t), intent(in) :: m
        real(dp), intent(in) :: x(:), rhs(:)
        real(dp), allocatable :: r(:), rowsum(:)
        integer :: j, p
        allocate (r(m%ncol), rowsum(m%ncol))
        r = rhs
        rowsum = 0
        do j = 1, m%ncol
            do p = m%col_ptr(j), m%col_ptr(j + 1) - 1
                r(m%row_idx(p)) = r(m%row_idx(p)) - m%val(p)*x(j)
                rowsum(m%row_idx(p)) = rowsum(m%row_idx(p)) + abs(m%val(p))
            end do
        end do
        berr = maxval(abs(r))/(maxval(rowsum)*maxval(abs(x)) + maxval(abs(rhs)))
    end function residual

    subroutine add(r, c, v)
        ! Append a stored (lower) entry; its mirror goes to the full list.
        integer, intent(in) :: r, c
        real(dp), intent(in) :: v
        nl = nl + 1
        lr(nl) = r; lc(nl) = c; lv(nl) = v
        nf = nf + 1
        fr(nf) = r; fc(nf) = c; fv(nf) = v
        if (r /= c) then
            nf = nf + 1
            fr(nf) = c; fc(nf) = r; fv(nf) = v
        end if
    end subroutine add

    subroutine finish(dim, a_full, a_low, a_up)
        integer, intent(in) :: dim
        type(csc_t), intent(out) :: a_full, a_low, a_up
        call csc_from_triplet(dim, dim, lr(1:nl), lc(1:nl), lv(1:nl), a_low, status)
        call csc_from_triplet(dim, dim, lc(1:nl), lr(1:nl), lv(1:nl), a_up, status)
        call csc_from_triplet(dim, dim, fr(1:nf), fc(1:nf), fv(1:nf), a_full, status)
    end subroutine finish

    subroutine saddle_point(scale, a_full, a_low, a_up)
        ! H: tridiagonal SPD (diagonal 2 scale + 1, off-diagonal -1) plus a
        ! second off-diagonal; B: mb sparse constraint rows; zero (2,2) block.
        real(dp), intent(in) :: scale
        type(csc_t), intent(out) :: a_full, a_low, a_up
        integer :: i, c
        nl = 0; nf = 0
        do i = 1, nh
            call add(i, i, 2.0_dp*scale + 1.0_dp)
            if (i < nh) call add(i + 1, i, -1.0_dp)
            if (i + 3 <= nh) call add(i + 3, i, 0.2_dp*scale)
        end do
        do c = 1, mb
            call add(nh + c, 1 + mod(5*c, nh), 1.0_dp + 0.1_dp*scale*real(c, dp))
            call add(nh + c, 1 + mod(11*c + 3, nh), -0.5_dp)
        end do
        call finish(n, a_full, a_low, a_up)
    end subroutine saddle_point

    subroutine perturb_upper(m)
        ! Round-off relative change of every strictly upper entry.
        type(csc_t), intent(inout) :: m
        integer :: j, p
        do j = 1, m%ncol
            do p = m%col_ptr(j), m%col_ptr(j + 1) - 1
                if (m%row_idx(p) < j) m%val(p) = m%val(p)*(1.0_dp + 3.0e-16_dp)
            end do
        end do
    end subroutine perturb_upper

    subroutine spd_check()
        integer, parameter :: g = 12, m = g*g
        integer :: i, j
        real(dp) :: rhs(m), xf(m), xl(m)
        type(csc_t) :: af, al, au
        type(sparse_solver_t) :: sf, sl
        nl = 0; nf = 0
        do j = 1, g
            do i = 1, g
                call add(i + (j - 1)*g, i + (j - 1)*g, 4.5_dp)
                if (i < g) call add(i + 1 + (j - 1)*g, i + (j - 1)*g, -1.0_dp)
                if (j < g) call add(i + j*g, i + (j - 1)*g, -1.0_dp)
            end do
        end do
        call finish(m, af, al, au)
        rhs = [(1.0_dp + cos(0.3_dp*real(i, dp)), i=1, m)]
        sf%backend_id = FORTSPARSE_BACKEND_NATIVE
        sf%options%symmetry = FORTSPARSE_SYM_SPD
        sl%backend_id = FORTSPARSE_BACKEND_NATIVE
        sl%options%symmetry = FORTSPARSE_SYM_SPD
        sl%options%triangle = FORTSPARSE_TRIANGLE_LOWER
        call sparse_factor(sf, af, status)
        call expect(status_ok(status), 'SPD full factor')
        call sparse_solve(sf, rhs, xf, status)
        call sparse_factor(sl, al, status)
        call expect(status_ok(status), 'SPD lower-triangle factor')
        call sparse_solve(sl, rhs, xl, status)
        call expect(all(xl == xf), 'SPD lower triangle gives the full-storage solution')
        call sparse_destroy(sf)
        call sparse_destroy(sl)
    end subroutine spd_check

end program test_fortsparse_symmetric_triangle
