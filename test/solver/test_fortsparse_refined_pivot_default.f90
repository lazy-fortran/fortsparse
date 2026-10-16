program test_fortsparse_refined_pivot_default
    ! Refined native LU default pivot threshold on a transport-like system with
    ! four interleaved fields (solution scales 1e13, 1e3, 1e-9, 1e-9). Rows of
    ! the first field couple to the small fields with coefficients up to 1e24,
    ! which dominate their row maxima but not |A||x|. Row-scaled maximum
    ! pivoting (threshold 1) then leaves a large componentwise backward error,
    ! and refinement stagnates for the transpose. The refined default (0.01)
    ! keeps the diagonal pivots. Independent oracle: real128 dense residuals.
    use fortsparse
    implicit none
    integer, parameter :: qp = selected_real_kind(30)
    integer, parameter :: nodes = 50, n = 4*nodes
    real(dp), parameter :: eps = epsilon(1.0_dp)
    real(dp) :: a(n,n), x(n), b(n)

    call build()
    call real_default()
    call real_explicit()
    call complex_default()
    print '(a)', 'refined native LU pivot default passed'
contains
    subroutine build()
        ! Five-point second (S) and first (D) derivative stencils with radially
        ! varying coefficients; a 1e-3 relative diagonal shift as time term.
        real(dp), parameter :: s(-2:2) = [-1.0_dp, 16.0_dp, -30.0_dp, 16.0_dp, -1.0_dp]/12.0_dp
        real(dp), parameter :: d(-2:2) = [1.0_dp, -8.0_dp, 0.0_dp, 8.0_dp, -1.0_dp]/12.0_dp
        real(dp), parameter :: scales(0:3) = [1.0e13_dp, 1.0e3_dp, 1.0e-9_dp, 1.0e-9_dp]
        integer :: i, k, r
        a = 0.0_dp
        do i = 0, nodes-1
            call couple(i, 0, 0, s, -1.8e7_dp, 0.3_dp, 1)
            call couple(i, 0, 1, d, 10.0_dp, 0.5_dp, 2)
            call couple(i, 0, 2, s, -2.2e23_dp, 1.0_dp, 3)
            call couple(i, 0, 3, s, -1.4e19_dp, 1.0_dp, 4)
            call couple(i, 1, 0, d, 1.0e-15_dp, 0.5_dp, 5)
            call couple(i, 1, 1, s, -1.8e7_dp, 0.3_dp, 6)
            call couple(i, 1, 2, d, 5.0e6_dp, 0.5_dp, 7)
            call couple(i, 1, 3, d, -1.4e7_dp, 0.5_dp, 8)
            call couple(i, 2, 0, s, -1.0e-21_dp, 0.5_dp, 9)
            call couple(i, 2, 2, s, -1.8e7_dp, 0.3_dp, 10)
            call couple(i, 2, 3, s, 3.0e-4_dp, 0.5_dp, 11)
            call couple(i, 3, 0, s, -4.0e-22_dp, 0.5_dp, 12)
            call couple(i, 3, 2, s, -6.0_dp, 0.5_dp, 13)
            call couple(i, 3, 3, s, -2.2e7_dp, 0.3_dp, 14)
            do k = 0, 3
                r = 4*i+k+1
                a(r,r) = 1.001_dp*a(r,r)
            end do
        end do
        do r = 1, n
            k = mod(r-1, 4)
            x(r) = scales(k)*(1.5_dp+sin(3.0_dp*real((r-1)/4,dp)/real(nodes-1,dp)+real(k,dp)))
        end do
        b = product(a, x)
    end subroutine build

    subroutine couple(i, row_field, col_field, stencil, base, ramp, phase)
        integer, intent(in) :: i, row_field, col_field, phase
        real(dp), intent(in) :: stencil(-2:2), base, ramp
        real(dp) :: t, c
        integer :: k
        t = real(i,dp)/real(nodes-1,dp)
        c = base*10.0_dp**(ramp*t+0.05_dp*sin(7.3_dp*real(i,dp)+real(phase,dp)))
        do k = -2, 2
            if (i+k < 0 .or. i+k >= nodes .or. stencil(k) == 0.0_dp) cycle
            a(4*i+row_field+1, 4*(i+k)+col_field+1) = &
                a(4*i+row_field+1, 4*(i+k)+col_field+1)+c*stencil(k)
        end do
    end subroutine couple

    subroutine real_default()
        type(sparse_solver_t) :: solver
        type(sparse_factor_stats_t) :: stats
        type(fortsparse_status_t) :: status
        type(csc_t) :: c
        real(dp) :: y(n), rhs(n,2), sol(n,2), nw, cw
        integer :: j
        call to_csc(a, c)
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        call sparse_factor(solver, c, status)
        call checked(status)
        call sparse_factor_stats(solver, stats, status)
        call check(stats%effective_pivot_tolerance == 1.0_dp, 'unrefined default threshold 1')
        solver%refine = .true.
        call sparse_factor(solver, c, status)
        call checked(status)
        call sparse_factor_stats(solver, stats, status)
        call check(stats%effective_pivot_tolerance == 0.01_dp, 'refined default threshold 0.01')
        call sparse_solve(solver, b, y, status)
        call checked(status)
        call errors(a, b, y, nw, cw)
        call check(cw <= 4.0_dp*eps .and. nw <= 32.0_dp*eps, 'refined default normal solve')
        call sparse_solve_op(solver, b, y, FORTSPARSE_OP_TRANSPOSE, status)
        call checked(status)
        call errors(transpose(a), b, y, nw, cw)
        call check(cw <= 4.0_dp*eps .and. nw <= 32.0_dp*eps, 'refined default transpose solve')
        rhs(:,1) = b
        rhs(:,2) = product(a, x(n:1:-1))
        call sparse_solve_block_op(solver, rhs, sol, FORTSPARSE_OP_TRANSPOSE, status)
        call checked(status)
        do j = 1, 2
            call errors(transpose(a), rhs(:,j), sol(:,j), nw, cw)
            call check(cw <= 4.0_dp*eps, 'refined default block transpose solve')
        end do
        call sparse_solve_block(solver, rhs, sol, status)
        call checked(status)
        do j = 1, 2
            call errors(a, rhs(:,j), sol(:,j), nw, cw)
            call check(cw <= 4.0_dp*eps, 'refined default block solve')
        end do
        y = b
        call sparse_solve(solver, y, status)
        call checked(status)
        call errors(a, b, y, nw, cw)
        call check(cw <= 4.0_dp*eps, 'refined default in-place solve')
        call sparse_destroy(solver)
    end subroutine real_default

    subroutine real_explicit()
        ! The fixture discriminates: explicit maximum pivoting leaves a large
        ! componentwise error, and an explicit threshold is honored when refined.
        type(sparse_solver_t) :: solver
        type(sparse_factor_stats_t) :: stats
        type(fortsparse_status_t) :: status
        type(csc_t) :: c
        real(dp) :: y(n), nw, cw
        call to_csc(a, c)
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%pivot_tolerance = 1.0_dp
        call sparse_factor(solver, c, status)
        call checked(status)
        call sparse_solve(solver, b, y, status)
        call checked(status)
        call errors(a, b, y, nw, cw)
        call check(nw <= 32.0_dp*eps .and. cw >= 1.0e-8_dp, &
            'maximum pivoting passes normwise with a large componentwise error')
        solver%refine = .true.
        call sparse_factor(solver, c, status)
        call checked(status)
        call sparse_factor_stats(solver, stats, status)
        call check(stats%effective_pivot_tolerance == 1.0_dp, 'explicit refined threshold')
        call sparse_destroy(solver)
    end subroutine real_explicit

    subroutine complex_default()
        ! Row-dependent complex phases: pivot moduli scale with whole rows, so
        ! the real fixture's pivoting behavior carries over.
        type(sparse_solver_t) :: solver
        type(sparse_factor_stats_t) :: stats
        type(fortsparse_status_t) :: status
        type(csc_z_t) :: c
        complex(dp) :: za(n,n), zx(n), zb(n), zy(n)
        real(dp) :: nw, cw
        integer :: i, j
        do j = 1, n
            do i = 1, n
                za(i,j) = a(i,j)*cmplx(1.0_dp, 0.5_dp*sin(real(i,dp)), dp)
            end do
        end do
        zx = x*cmplx(1.0_dp, 0.25_dp, dp)
        zb = product_z(za, zx)
        call to_csc_z(za, c)
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%pivot_tolerance = 1.0_dp
        call sparse_factor(solver, c, status)
        call checked(status)
        call sparse_solve(solver, zb, zy, status)
        call checked(status)
        call errors_z(za, zb, zy, nw, cw)
        call check(cw >= 1.0e-8_dp, 'complex maximum pivoting has a large componentwise error')
        solver%options%pivot_tolerance = -1.0_dp
        solver%refine = .true.
        call sparse_factor(solver, c, status)
        call checked(status)
        call sparse_factor_stats(solver, stats, status)
        call check(stats%effective_pivot_tolerance == 0.01_dp, 'complex refined default threshold')
        call sparse_solve(solver, zb, zy, status)
        call checked(status)
        call errors_z(za, zb, zy, nw, cw)
        call check(cw <= 4.0_dp*eps .and. nw <= 32.0_dp*eps, 'complex refined default normal solve')
        call sparse_solve_op(solver, zb, zy, FORTSPARSE_OP_TRANSPOSE, status)
        call checked(status)
        call errors_z(transpose(za), zb, zy, nw, cw)
        call check(cw <= 4.0_dp*eps, 'complex refined default transpose solve')
        call sparse_solve_op(solver, zb, zy, FORTSPARSE_OP_ADJOINT, status)
        call checked(status)
        call errors_z(conjg(transpose(za)), zb, zy, nw, cw)
        call check(cw <= 4.0_dp*eps, 'complex refined default adjoint solve')
        call sparse_destroy(solver)
    end subroutine complex_default

    subroutine to_csc(m, c)
        real(dp), intent(in) :: m(:,:)
        type(csc_t), intent(out) :: c
        type(fortsparse_status_t) :: status
        integer :: rows(count(m /= 0.0_dp)), cols(size(rows)), i, j, k
        real(dp) :: values(size(rows))
        k = 0
        do j = 1, size(m,2)
            do i = 1, size(m,1)
                if (m(i,j) /= 0.0_dp) then
                    k = k+1
                    rows(k) = i
                    cols(k) = j
                    values(k) = m(i,j)
                end if
            end do
        end do
        call csc_from_triplet(size(m,1), size(m,2), rows, cols, values, c, status)
        call checked(status)
    end subroutine to_csc

    subroutine to_csc_z(m, c)
        complex(dp), intent(in) :: m(:,:)
        type(csc_z_t), intent(out) :: c
        type(fortsparse_status_t) :: status
        integer :: rows(count(m /= (0.0_dp,0.0_dp))), cols(size(rows)), i, j, k
        complex(dp) :: values(size(rows))
        k = 0
        do j = 1, size(m,2)
            do i = 1, size(m,1)
                if (m(i,j) /= (0.0_dp,0.0_dp)) then
                    k = k+1
                    rows(k) = i
                    cols(k) = j
                    values(k) = m(i,j)
                end if
            end do
        end do
        call csc_from_triplet(size(m,1), size(m,2), rows, cols, values, c, status)
        call checked(status)
    end subroutine to_csc_z

    subroutine errors(m, rhs, y, normwise, componentwise)
        ! Independent real128 normwise and componentwise backward errors of
        ! the dense matrix, accumulating its nonzero entries only.
        real(dp), intent(in) :: m(:,:), rhs(:), y(:)
        real(dp), intent(out) :: normwise, componentwise
        real(qp) :: r(size(y)), d(size(y)), rowsum(size(y)), term
        integer :: i, j
        r = real(rhs,qp)
        d = abs(real(rhs,qp))
        rowsum = 0.0_qp
        do j = 1, size(m,2)
            do i = 1, size(m,1)
                if (m(i,j) == 0.0_dp) cycle
                term = real(m(i,j),qp)*real(y(j),qp)
                r(i) = r(i)-term
                d(i) = d(i)+abs(term)
                rowsum(i) = rowsum(i)+abs(real(m(i,j),qp))
            end do
        end do
        normwise = real(maxval(abs(r))/(maxval(rowsum)*maxval(abs(real(y,qp))) &
            +maxval(abs(real(rhs,qp)))), dp)
        componentwise = real(maxval(abs(r)/merge(d, 1.0_qp, d > 0.0_qp)), dp)
    end subroutine errors

    subroutine errors_z(m, rhs, y, normwise, componentwise)
        complex(dp), intent(in) :: m(:,:), rhs(:), y(:)
        real(dp), intent(out) :: normwise, componentwise
        complex(qp) :: r(size(y)), term
        real(qp) :: d(size(y)), rowsum(size(y))
        integer :: i, j
        r = cmplx(rhs,kind=qp)
        d = abs(r)
        rowsum = 0.0_qp
        do j = 1, size(m,2)
            do i = 1, size(m,1)
                if (m(i,j) == (0.0_dp,0.0_dp)) cycle
                term = cmplx(m(i,j),kind=qp)*cmplx(y(j),kind=qp)
                r(i) = r(i)-term
                d(i) = d(i)+abs(cmplx(m(i,j),kind=qp))*abs(cmplx(y(j),kind=qp))
                rowsum(i) = rowsum(i)+abs(cmplx(m(i,j),kind=qp))
            end do
        end do
        normwise = real(maxval(abs(r))/(maxval(rowsum)*maxval(abs(cmplx(y,kind=qp))) &
            +maxval(abs(cmplx(rhs,kind=qp)))), dp)
        componentwise = real(maxval(abs(r)/merge(d, 1.0_qp, d > 0.0_qp)), dp)
    end subroutine errors_z

    function product(m, v) result(p)
        ! real128 product, rounded once to binary64.
        real(dp), intent(in) :: m(:,:), v(:)
        real(dp) :: p(size(m,1))
        real(qp) :: acc(size(m,1))
        integer :: i, j
        acc = 0.0_qp
        do j = 1, size(m,2)
            do i = 1, size(m,1)
                if (m(i,j) /= 0.0_dp) acc(i) = acc(i)+real(m(i,j),qp)*real(v(j),qp)
            end do
        end do
        p = real(acc, dp)
    end function product

    function product_z(m, v) result(p)
        complex(dp), intent(in) :: m(:,:), v(:)
        complex(dp) :: p(size(m,1))
        complex(qp) :: acc(size(m,1))
        integer :: i, j
        acc = (0.0_qp, 0.0_qp)
        do j = 1, size(m,2)
            do i = 1, size(m,1)
                if (m(i,j) /= (0.0_dp,0.0_dp)) &
                    acc(i) = acc(i)+cmplx(m(i,j),kind=qp)*cmplx(v(j),kind=qp)
            end do
        end do
        p = cmplx(acc, kind=dp)
    end function product_z

    subroutine checked(status)
        type(fortsparse_status_t), intent(in) :: status
        if (.not. status_ok(status)) then
            print '(a)', trim(status%msg)
            error stop 'unexpected status'
        end if
    end subroutine checked

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print '(a)', message
            error stop 1
        end if
    end subroutine check
end program test_fortsparse_refined_pivot_default
