program test_fortsparse_componentwise_refinement
    ! Componentwise (Oettli-Prager) refinement in the native general adapter:
    ! badly row-scaled systems whose unrefined LU solution passes the normwise
    ! 32*eps gate but has a large componentwise backward error, stagnation,
    ! the three-correction cap and the reported errors. Independent oracle:
    ! real128 residuals of the original dense matrix.
    use fortsparse
    implicit none
    integer, parameter :: qp = selected_real_kind(30)
    real(dp), parameter :: eps = epsilon(1.0_dp)
    integer :: seed
    real(dp), allocatable :: a(:,:), x(:), b(:)

    call row_scaled()
    call stagnation()
    call perturbed_stagnation()
    call step_cap()
    print '(a)', 'native componentwise refinement passed'
contains
    real(dp) function uniform()
        seed = int(modulo(1103515245_i8*seed+12345_i8, 2147483648_i8))
        uniform = real(seed, dp)/2147483648.0_dp
    end function uniform

    subroutine build(n, density, diagonal, scale)
        ! Sparse rows scaled by 10**(+-scale), a subdiagonal cycle carrying
        ! the pivots, a small diagonal and a graded signed solution.
        integer, intent(in) :: n, density, scale
        real(dp), intent(in) :: diagonal
        integer :: i, j, k
        real(dp) :: row_scale
        if (allocated(a)) deallocate(a, x, b)
        allocate(a(n,n), x(n), b(n))
        a = 0.0_dp
        do k = 1, density*n
            i = 1+int(uniform()*n)
            j = 1+int(uniform()*n)
            a(i,j) = 2.0_dp*uniform()-1.0_dp
        end do
        do j = 1, n
            a(j,j) = diagonal*uniform()
            a(mod(j,n)+1,j) = 1.0_dp+uniform()
        end do
        do i = 1, n
            row_scale = 10.0_dp**(real(scale,dp)*(2.0_dp*uniform()-1.0_dp))
            a(i,:) = row_scale*a(i,:)
        end do
        do j = 1, n
            x(j) = 10.0_dp**(-real(scale,dp)*uniform())*(2.0_dp*uniform()-1.0_dp)
        end do
        b = real(matmul(real(a,qp), real(x,qp)), dp)
    end subroutine build

    subroutine to_csc(m, c)
        real(dp), intent(in) :: m(:,:)
        type(csc_t), intent(out) :: c
        type(fortsparse_status_t) :: status
        integer :: rows(size(m)), cols(size(m)), i, j, k
        real(dp) :: values(size(m))
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
        call csc_from_triplet(size(m,1), size(m,2), rows(:k), cols(:k), values(:k), &
            c, status)
        call checked(status)
    end subroutine to_csc

    subroutine errors(m, rhs, y, normwise, componentwise)
        ! Independent real128 normwise and componentwise backward errors.
        real(dp), intent(in) :: m(:,:), rhs(:), y(:)
        real(dp), intent(out) :: normwise, componentwise
        real(qp) :: r(size(y)), d(size(y))
        r = real(rhs,qp)-matmul(real(m,qp), real(y,qp))
        d = abs(real(rhs,qp))+matmul(abs(real(m,qp)), abs(real(y,qp)))
        normwise = real(maxval(abs(r))/(maxval(sum(abs(real(m,qp)),2))* &
            maxval(abs(real(y,qp)))+maxval(abs(real(rhs,qp)))), dp)
        componentwise = real(maxval(abs(r)/merge(d, 1.0_qp, d > 0.0_qp)), dp)
    end subroutine errors

    subroutine errors_z(m, rhs, y, normwise, componentwise)
        complex(dp), intent(in) :: m(:,:), rhs(:), y(:)
        real(dp), intent(out) :: normwise, componentwise
        complex(qp) :: r(size(y))
        real(qp) :: d(size(y))
        r = cmplx(rhs,kind=qp)-matmul(cmplx(m,kind=qp), cmplx(y,kind=qp))
        d = abs(cmplx(rhs,kind=qp))+matmul(abs(cmplx(m,kind=qp)), abs(cmplx(y,kind=qp)))
        normwise = real(maxval(abs(r))/(maxval(sum(abs(cmplx(m,kind=qp)),2))* &
            maxval(abs(cmplx(y,kind=qp)))+maxval(abs(cmplx(rhs,kind=qp)))), dp)
        componentwise = real(maxval(abs(r)/merge(d, 1.0_qp, d > 0.0_qp)), dp)
    end subroutine errors_z

    subroutine row_scaled()
        ! Rows scaled over 1e+-5: unrefined LU passes the normwise gate with
        ! a componentwise error far above eps; refinement must reduce it.
        type(sparse_solver_t) :: plain, refined
        type(sparse_factor_stats_t) :: stats
        type(fortsparse_status_t) :: status
        type(csc_t) :: c
        type(csc_z_t) :: z
        complex(dp), allocatable :: za(:,:), zx(:), zb(:), zy(:)
        real(dp) :: y(30), nw, cw, block_b(30,2), block_y(30,2), worst
        integer :: trial, method, i, j
        seed = 4242
        do trial = 1, 5
            call build(30, 3, 0.1_dp, trial)
        end do
        call to_csc(a, c)
        plain%backend_id = FORTSPARSE_BACKEND_NATIVE
        call sparse_factor(plain, c, status)
        call checked(status)
        call sparse_solve(plain, b, y, status)
        call checked(status)
        call errors(a, b, y, nw, cw)
        call check(nw <= 32.0_dp*eps, 'unrefined LU passes the normwise gate')
        call check(cw >= 100.0_dp*eps, 'unrefined componentwise error is large')
        do method = FORTSPARSE_LU_SCALAR, FORTSPARSE_LU_SUPERNODAL
            refined%backend_id = FORTSPARSE_BACKEND_NATIVE
            refined%refine = .true.
            refined%options%lu_method = method
            call sparse_factor(refined, c, status)
            call checked(status)
            call sparse_solve(refined, b, y, status)
            call checked(status)
            call sparse_factor_stats(refined, stats, status)
            call errors(a, b, y, nw, cw)
            call check(stats%refinement_steps >= 1 .and. stats%refinement_steps <= 3, &
                'componentwise criterion takes a correction after the normwise pass')
            call check(cw <= 4.0_dp*eps, 'refined componentwise error at a few eps')
            call check(nw <= 32.0_dp*eps, 'refined normwise error')
            call check(stats%backward_error >= 0.0_dp .and. &
                stats%backward_error <= 32.0_dp*eps, 'reported normwise error')
            call check(stats%componentwise_backward_error >= 0.0_dp .and. &
                stats%componentwise_backward_error <= 4.0_dp*eps, &
                'reported componentwise error')
            call check(abs(stats%componentwise_backward_error-cw) <= 2.0_dp*eps, &
                'reported componentwise error matches the real128 oracle')
            ! Transpose solve: the denominator accumulates per column.
            call sparse_solve_op(refined, b, y, FORTSPARSE_OP_TRANSPOSE, status)
            call checked(status)
            call sparse_factor_stats(refined, stats, status)
            call errors(transpose(a), b, y, nw, cw)
            call check(cw <= 4.0_dp*eps .and. &
                stats%componentwise_backward_error <= 4.0_dp*eps, &
                'refined transpose componentwise error')
            ! Block: the reported errors are the maxima over the columns.
            block_b(:,1) = b
            block_b(:,2) = real(matmul(real(a,qp), real(x(30:1:-1),qp)), dp)
            call sparse_solve_block(refined, block_b, block_y, status)
            call checked(status)
            call sparse_factor_stats(refined, stats, status)
            call errors(a, block_b(:,1), block_y(:,1), nw, worst)
            call errors(a, block_b(:,2), block_y(:,2), nw, cw)
            worst = max(worst, cw)
            call check(worst <= 4.0_dp*eps .and. &
                abs(stats%componentwise_backward_error-worst) <= 2.0_dp*eps, &
                'block componentwise maximum')
            call sparse_destroy(refined)
        end do
        call sparse_destroy(plain)
        ! Complex LU with the same scaling.
        allocate(za(30,30))
        do j = 1, 30
            do i = 1, 30
                za(i,j) = cmplx(a(i,j), real(1-2*mod(i+j,2),dp)*0.5_dp*a(i,j), dp)
            end do
        end do
        zx = cmplx(x, 0.25_dp*x(30:1:-1), dp)
        zb = cmplx(matmul(cmplx(za,kind=qp), cmplx(zx,kind=qp)), kind=dp)
        allocate(zy(30))
        call to_csc_z(za, z)
        refined%backend_id = FORTSPARSE_BACKEND_NATIVE
        refined%refine = .true.
        refined%options%lu_method = FORTSPARSE_LU_SCALAR
        call sparse_factor(refined, z, status)
        call checked(status)
        call sparse_solve(refined, zb, zy, status)
        call checked(status)
        call sparse_factor_stats(refined, stats, status)
        call errors_z(za, zb, zy, nw, cw)
        call check(cw <= 4.0_dp*eps .and. nw <= 32.0_dp*eps, 'complex refined errors')
        call check(stats%componentwise_backward_error >= 0.0_dp .and. &
            abs(stats%componentwise_backward_error-cw) <= 2.0_dp*eps, &
            'complex reported componentwise error')
        call sparse_destroy(refined)
    end subroutine row_scaled

    subroutine to_csc_z(m, c)
        complex(dp), intent(in) :: m(:,:)
        type(csc_z_t), intent(out) :: c
        type(fortsparse_status_t) :: status
        integer :: rows(size(m)), cols(size(m)), i, j, k
        complex(dp) :: values(size(m))
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
        call csc_from_triplet(size(m,1), size(m,2), rows(:k), cols(:k), values(:k), &
            c, status)
        call checked(status)
    end subroutine to_csc_z

    subroutine stagnation()
        ! Dense rows put the rounding floor of the computed componentwise
        ! error near or above eps. A solve that returns with that error above
        ! eps before the cap stopped because the error did not halve; whether
        ! one does depends on rounding, so perturbed_stagnation tests the rule.
        type(sparse_solver_t) :: solver
        type(sparse_factor_stats_t) :: stats
        type(fortsparse_status_t) :: status
        type(csc_t) :: c
        real(dp) :: y(40), nw, cw
        integer :: trial, method, stagnated
        seed = 4242
        stagnated = 0
        do trial = 1, 6
            call build(40, 40, 0.1_dp, trial)
            call to_csc(a, c)
            do method = FORTSPARSE_LU_SCALAR, FORTSPARSE_LU_SUPERNODAL
                solver%backend_id = FORTSPARSE_BACKEND_NATIVE
                solver%refine = .true.
                solver%options%lu_method = method
                call sparse_factor(solver, c, status)
                call checked(status)
                call sparse_solve(solver, b, y, status)
                call checked(status)
                call sparse_factor_stats(solver, stats, status)
                call errors(a, b, y, nw, cw)
                call check(nw <= 32.0_dp*eps .and. cw <= 8.0_dp*eps, &
                    'dense-row refined errors')
                call check(stats%refinement_steps <= 3, 'correction cap')
                if (stats%componentwise_backward_error > eps .and. &
                    stats%refinement_steps < 3) stagnated = stagnated+1
                call sparse_destroy(solver)
            end do
        end do
        print '(a,i0,a)', 'stagnation stops: ', stagnated, ' of 12 solves'
    end subroutine stagnation

    subroutine perturbed_stagnation()
        ! Refinement that contracts by a known factor rho > 1/2 per step.
        ! 2x2 blocks [1 .1; .1 1], except s*[1 c; c 1] (c = .9, s = 2**-50)
        ! in front, n = 1000 so that one perturbed pivot passes the gate.
        ! Matching scaling makes the special block [1 c t; c/t 1] (t > 0,
        ! scaled entries at most 1), whose second pivot p = 1-c**2 = .19 is
        ! below tau*||M||_inf with tau = .4 and .76 <= tau*||M||_inf <= .8;
        ! the unit blocks' pivots .99 are not. The factor is that of A+E,
        ! with E nonzero only in that pivot, raised by delta to
        ! p' = tau*||M||_inf. The refinement error e_k+1 = (A+E)^-1 E e_k:
        ! (A+E)^-1 E has rank one with the eigenvalue rho = delta/p' in
        ! [.75, .7625], and e_0 = -(A+E)^-1 E x is its eigenvector, so every
        ! correction multiplies the residual by exactly rho. It lies in the
        ! block's second row, whose |A||x|+|b| changes by a relative O(eta)
        ! with x = (1, eta) there: the componentwise error contracts by
        ! rho*(1+O(eta)) and cannot halve. The normwise error stays below
        ! s*eta, so the rule stops after one correction with the error far
        ! above eps, and so would every later correction.
        integer, parameter :: n = 1000
        real(dp), parameter :: c = 0.9_dp, s = 2.0_dp**(-50), eta = 2.0_dp**(-10)
        type(sparse_solver_t) :: refined, plain
        type(sparse_factor_stats_t) :: stats
        type(fortsparse_status_t) :: status
        type(csc_t) :: m
        real(dp) :: y(n), x0(n), x2(n), r(n), dx(n), nw, cw0, cw1, cw2
        integer :: j
        if (allocated(a)) deallocate(a, x, b)
        allocate(a(n,n), x(n), b(n))
        a = 0.0_dp
        do j = 1, n, 2
            a(j:j+1,j:j+1) = reshape([1.0_dp, 0.1_dp, 0.1_dp, 1.0_dp], [2,2])
        end do
        a(1:2,1:2) = s*reshape([1.0_dp, c, c, 1.0_dp], [2,2])
        x = 1.0_dp
        x(2) = eta
        b = real(matmul(real(a,qp), real(x,qp)), dp)
        call to_csc(a, m)
        refined%backend_id = FORTSPARSE_BACKEND_NATIVE
        refined%refine = .true.
        refined%options%lu_method = FORTSPARSE_LU_SUPERNODAL
        refined%options%pivot_perturbation = 0.4_dp
        call sparse_factor(refined, m, status)
        call checked(status)
        call sparse_factor_stats(refined, stats, status)
        call check(stats%native_factor_family == 'LU_SUPERNODAL' .and. &
            stats%perturbed_pivots == 1, 'one perturbed pivot in the special block')
        call sparse_solve(refined, b, y, status)
        call checked(status)
        call sparse_factor_stats(refined, stats, status)
        call errors(a, b, y, nw, cw1)
        call check(stats%refinement_steps == 1, 'componentwise stagnation stops after one correction')
        call check(stats%backward_error <= 32.0_dp*eps .and. nw <= 32.0_dp*eps, &
            'normwise gate passes before refinement')
        call check(stats%componentwise_backward_error > 1.0e6_dp*eps .and. &
            abs(stats%componentwise_backward_error-cw1) <= 1.0e-3_dp*cw1, &
            'stagnated componentwise error matches the real128 oracle')
        ! The contraction itself, with the same factor unrefined: x0 is the
        ! unrefined solution, y = x1 and x2 one further correction.
        plain%backend_id = FORTSPARSE_BACKEND_NATIVE
        plain%options%lu_method = FORTSPARSE_LU_SUPERNODAL
        plain%options%pivot_perturbation = 0.4_dp
        plain%options%verify_factor = .false.
        plain%options%verify_perturbed = .false.
        call sparse_factor(plain, m, status)
        call checked(status)
        call sparse_solve(plain, b, x0, status)
        call checked(status)
        r = real(real(b,qp)-matmul(real(a,qp), real(y,qp)), dp)
        call sparse_solve(plain, r, dx, status)
        call checked(status)
        x2 = y+dx
        call errors(a, b, x0, nw, cw0)
        call errors(a, b, x2, nw, cw2)
        call check(cw1 >= 0.74_dp*cw0 .and. cw1 <= 0.77_dp*cw0 .and. &
            cw2 >= 0.74_dp*cw1 .and. cw2 <= 0.77_dp*cw1, &
            'each correction contracts the componentwise error by rho')
        call sparse_destroy(plain)
        call sparse_destroy(refined)
    end subroutine perturbed_stagnation

    subroutine step_cap()
        ! A 20x20 block with diagonals near 1e-5, rows times 1e-12, next to a
        ! unit block. Threshold pivoting (1e-20) accepts the small diagonal
        ! pivots: growth makes refinement contract slowly. The normwise gate
        ! passes early (the block rows are tiny) while the componentwise
        ! error still halves at the third correction: such a solve is
        ! accepted at the cap, and one more independent correction with an
        ! unrefined solver still halves the error (cap, not stagnation).
        integer, parameter :: m = 20, n = 24
        type(sparse_solver_t) :: solver, plain
        type(sparse_factor_stats_t) :: stats
        type(fortsparse_status_t) :: status, solved
        type(csc_t) :: c
        real(dp) :: y(n), dy(n), r(n), nw, cw, cw_next
        integer :: trial, capped, i, j, k
        seed = 99
        capped = 0
        if (allocated(a)) deallocate(a, x, b)
        allocate(a(n,n), x(n), b(n))
        do trial = 1, 8
            a = 0.0_dp
            do k = 1, 3*m
                i = 1+int(uniform()*m)
                j = 1+int(uniform()*m)
                a(i,j) = 2.0_dp*uniform()-1.0_dp
            end do
            do j = 1, m
                a(j,j) = 1.0e-5_dp*(1.0_dp+uniform())
                a(mod(j,m)+1,j) = 1.0_dp+uniform()
            end do
            a(1:m,1:m) = 1.0e-12_dp*a(1:m,1:m)
            do j = m+1, n
                a(j,j) = 1.0_dp
            end do
            do j = 1, n
                x(j) = 2.0_dp*uniform()-1.0_dp
            end do
            b = real(matmul(real(a,qp), real(x,qp)), dp)
            call to_csc(a, c)
            solver%backend_id = FORTSPARSE_BACKEND_NATIVE
            solver%refine = .true.
            solver%options%pivot_tolerance = 1.0e-20_dp
            call sparse_factor(solver, c, status)
            call checked(status)
            call sparse_solve(solver, b, y, solved)
            call sparse_factor_stats(solver, stats, status)
            call checked(status)
            call check(stats%refinement_steps <= 3, 'correction cap respected')
            if (solved%code /= FORTSPARSE_OK) then
                call check(solved%code == FORTSPARSE_SINGULAR .and. &
                    stats%backward_error > 32.0_dp*eps, &
                    'explicit normwise gate failure')
                call sparse_destroy(solver)
                cycle
            end if
            call errors(a, b, y, nw, cw)
            call check(stats%backward_error <= 32.0_dp*eps .and. nw <= 32.0_dp*eps, &
                'accepted normwise gate')
            if (stats%refinement_steps == 3 .and. &
                stats%componentwise_backward_error > 1.0e3_dp*eps) then
                capped = capped+1
                call check(cw > 1.0e2_dp*eps, 'oracle confirms the capped error')
                plain%backend_id = FORTSPARSE_BACKEND_NATIVE
                plain%options%pivot_tolerance = 1.0e-20_dp
                call sparse_factor(plain, c, status)
                call checked(status)
                r = real(real(b,qp)-matmul(real(a,qp), real(y,qp)), dp)
                call sparse_solve(plain, r, dy, status)
                call checked(status)
                call errors(a, b, y+dy, nw, cw_next)
                call check(2.0_dp*cw_next <= cw, &
                    'cap, not stagnation, stopped refinement')
                call sparse_destroy(plain)
            end if
            call sparse_destroy(solver)
        end do
        print '(a,i0,a)', 'componentwise cap stops: ', capped, ' of 8 solves'
        call check(capped >= 2, 'still-converging componentwise errors hit the cap')
    end subroutine step_cap

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
end program test_fortsparse_componentwise_refinement
