program test_fortsparse_slu_z_api
    ! Public API with lu_method = SUPERNODAL on complex (non-Hermitian)
    ! matrices: known solutions for A, A^T and A^H with refinement, block
    ! RHS, same-pattern refactor reuse, implicit derivatives (JVP against a
    ! central difference, VJP through the adjoint solve of the primal
    ! factors against the real adjoint identity), the gates and scalar
    ! fallback, singular input and cap exhaustion with recovery.
    use fortsparse
    implicit none
    integer :: seed, trial
    complex(dp), allocatable :: dense(:,:)

    seed = 90210
    do trial = 1, 6
        call random_dense(8+17*trial, trial, dense)
        call check_public(dense)
    end do
    call random_dense(40, 3, dense)
    call check_derivatives(dense)
    call check_fallback()
    call check_singular()
    call random_dense(70, 6, dense)
    call check_caps(dense)
    print '(a)', 'complex supernodal LU public API passed'
contains
    real(dp) function uniform()
        seed = int(modulo(1103515245_i8*seed+12345_i8, 2147483648_i8))
        uniform = real(seed, dp)/2147483648.0_dp
    end function uniform

    complex(dp) function random_z()
        random_z = cmplx(2.0_dp*uniform()-1.0_dp, 2.0_dp*uniform()-1.0_dp, dp)
    end function random_z

    subroutine random_dense(n, trial, a)
        ! Near-symmetric pattern; every third case has a zero diagonal and
        ! a random matching carries the dominant entries.
        integer, intent(in) :: n, trial
        complex(dp), allocatable, intent(out) :: a(:,:)
        integer :: perm(n), i, j, k, swap
        allocate(a(n,n))
        a = (0.0_dp, 0.0_dp)
        perm = [(i, i=1, n)]
        if (mod(trial, 3) == 0) then
            do i = n, 2, -1
                j = 1+int(uniform()*i)
                swap = perm(i)
                perm(i) = perm(j)
                perm(j) = swap
            end do
        end if
        do k = 1, 3*n
            i = 1+int(uniform()*n)
            j = 1+int(uniform()*n)
            a(i,j) = random_z()
            if (uniform() < 0.85_dp) a(j,i) = random_z()
        end do
        if (mod(trial, 3) == 0) then
            do j = 1, n
                a(j,j) = (0.0_dp, 0.0_dp)
            end do
        end if
        do j = 1, n
            a(perm(j),j) = (2.0_dp+uniform())*exp(cmplx(0.0_dp, 6.0_dp*uniform(), dp))
        end do
    end subroutine random_dense

    subroutine to_csc(a, matrix)
        complex(dp), intent(in) :: a(:,:)
        type(csc_z_t), intent(out) :: matrix
        integer, allocatable :: rows(:), cols(:)
        complex(dp), allocatable :: values(:)
        type(fortsparse_status_t) :: status
        integer :: n, i, j, k
        n = size(a,1)
        k = count(a /= (0.0_dp, 0.0_dp))
        allocate(rows(k), cols(k), values(k))
        k = 0
        do j = 1, n
            do i = 1, n
                if (a(i,j) == (0.0_dp, 0.0_dp)) cycle
                k = k+1
                rows(k) = i
                cols(k) = j
                values(k) = a(i,j)
            end do
        end do
        call csc_from_triplet(n, n, rows, cols, values, matrix, status)
        call check(status_ok(status), 'CSC construction')
    end subroutine to_csc

    subroutine supernodal(solver)
        type(sparse_solver_t), intent(inout) :: solver
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%lu_method = FORTSPARSE_LU_SUPERNODAL
    end subroutine supernodal

    function op_matrix(a, op) result(m)
        complex(dp), intent(in) :: a(:,:)
        integer, intent(in) :: op
        complex(dp), allocatable :: m(:,:)
        m = a
        if (op == FORTSPARSE_OP_TRANSPOSE) m = transpose(a)
        if (op == FORTSPARSE_OP_ADJOINT) m = conjg(transpose(a))
    end function op_matrix

    subroutine check_public(a)
        complex(dp), intent(in) :: a(:,:)
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: stats
        type(csc_z_t) :: matrix
        complex(dp), allocatable :: exact(:), b(:), x(:), block(:,:), result(:,:)
        integer :: n, i, op, j
        logical :: reused
        n = size(a,1)
        call to_csc(a, matrix)
        call supernodal(solver)
        solver%refine = .true.
        call sparse_factor(solver, matrix, status)
        call check(status_ok(status), 'complex supernodal factor: '//trim(status%msg))
        call sparse_factor_stats(solver, stats, status)
        call check(stats%native_factor_family == 'LU_SUPERNODAL', 'family reported')
        call check(stats%lu_method_used == FORTSPARSE_LU_SUPERNODAL, 'method reported')
        call check(.not. stats%lu_fallback, 'no fallback on a regular matrix')
        call check(stats%matching_used == FORTSPARSE_MATCHING_WEIGHTED, 'weighted matching')
        call check(stats%perturbed_pivots == 0, 'no perturbation')
        call check(stats%pivot_growth > 0.0_dp .and. stats%min_pivot > 0.0_dp, &
            'pivot statistics')
        call check(stats%factor_nnz > 0 .and. stats%host_current_bytes > 0, 'sizes reported')
        allocate(exact(n), b(n), x(n), block(n,3), result(n,3))
        exact = [(cmplx(cos(real(i,dp))+0.5_dp, sin(real(3*i,dp)), dp), i=1, n)]
        do op = FORTSPARSE_OP_NORMAL, FORTSPARSE_OP_ADJOINT
            b = matmul(op_matrix(a, op), exact)
            call sparse_solve_op(solver, b, x, op, status)
            call check(status_ok(status), 'solve operation')
            call check(maxval(abs(x-exact)) <= 1.0e-10_dp*maxval(abs(exact)), &
                'known solution recovered (N, T, H)')
            call sparse_factor_stats(solver, stats, status)
            call check(stats%backward_error <= 32.0_dp*epsilon(1.0_dp), &
                'refined backward error gate')
            do j = 1, 3
                block(:,j) = cmplx(real(j,dp), real(1-j,dp), dp)*b
            end do
            call sparse_solve_block_op(solver, block, result, op, status)
            call check(status_ok(status), 'block operation')
            call check(maxval(abs(result(:,3)-(3.0_dp,-2.0_dp)*exact)) <= &
                4.0e-10_dp*maxval(abs(exact)), 'block oracle')
        end do
        ! Same pattern, new complex values: the analysis is reused.
        matrix%val = (1.5_dp, 0.5_dp)*matrix%val
        call sparse_refactor(solver, matrix, status, reused)
        call check(status_ok(status), 'refactor')
        call check(reused, 'same-pattern refactor reuses the analysis')
        b = (1.5_dp, 0.5_dp)*matmul(a, exact)
        call sparse_solve(solver, b, x, status)
        call check(maxval(abs(x-exact)) <= 1.0e-10_dp*maxval(abs(exact)), 'refactor solve')
        call sparse_factor_stats(solver, stats, status)
        call check(stats%native_factor_family == 'LU_SUPERNODAL', 'refactor family')
        call sparse_destroy(solver)
    end subroutine check_public

    subroutine check_derivatives(a)
        ! x(t) solves (A + t Adot) x = b + t bdot. The VJP uses the adjoint
        ! solve on the primal supernodal factors.
        complex(dp), intent(in) :: a(:,:)
        type(sparse_solver_t) :: solver, plus, minus
        type(fortsparse_status_t) :: status
        type(csc_z_t) :: matrix, tangent, shifted
        complex(dp), allocatable :: b(:), bdot(:), x(:), xdot(:), xp(:), xm(:)
        complex(dp), allocatable :: adot(:,:), xbar(:), bbar(:), abar(:)
        real(dp) :: h, pairing, adjoint
        integer :: n, i, j, k
        n = size(a,1)
        call to_csc(a, matrix)
        allocate(adot(n,n), b(n), bdot(n), x(n), xdot(n), xp(n), xm(n), xbar(n), &
            bbar(n), abar(matrix%nnz))
        adot = (0.0_dp, 0.0_dp)
        ! Bounded complex tangent on the pattern of A.
        where (a /= (0.0_dp, 0.0_dp)) adot = cmplx(sin(3.0_dp*real(a,dp)), &
            0.5_dp*cos(2.0_dp*aimag(a)), dp)
        call to_csc(adot, tangent)
        b = [(cmplx(sin(real(2*i,dp)), cos(real(i,dp)), dp), i=1, n)]
        bdot = [(cmplx(cos(real(i,dp)), 0.3_dp, dp), i=1, n)]
        call supernodal(solver)
        call sparse_factor(solver, matrix, status)
        call check(status_ok(status), 'derivative primal factor')
        call sparse_solve(solver, b, x, status)
        call sparse_solve_jvp(solver, tangent, x, bdot, xdot, status)
        call check(status_ok(status), 'JVP')
        h = 1.0e-6_dp
        call to_csc(a+h*adot, shifted)
        call supernodal(plus)
        call sparse_factor(plus, shifted, status)
        call sparse_solve(plus, b+h*bdot, xp, status)
        call to_csc(a-h*adot, shifted)
        call supernodal(minus)
        call sparse_factor(minus, shifted, status)
        call sparse_solve(minus, b-h*bdot, xm, status)
        call check(maxval(abs(xdot-(xp-xm)/(2.0_dp*h))) <= 1.0e-6_dp*maxval(abs(xdot)), &
            'JVP equals the central difference')
        xbar = [(cmplx(1.0_dp/real(i,dp), -0.5_dp/real(i,dp), dp), i=1, n)]
        call sparse_solve_vjp(solver, matrix, x, xbar, bbar, abar, status, &
            use_primal_factors=.true.)
        call check(status_ok(status), 'VJP through the adjoint solve')
        pairing = real(dot_product(xbar, xdot), dp)
        adjoint = real(dot_product(bbar, bdot), dp)
        k = 0
        do j = 1, n
            do i = matrix%col_ptr(j), matrix%col_ptr(j+1)-1
                k = k+1
                adjoint = adjoint+real(conjg(abar(k))*adot(matrix%row_idx(i), j), dp)
            end do
        end do
        call check(abs(pairing-adjoint) <= 1.0e-10_dp*(abs(pairing)+1.0_dp), &
            'adjoint identity Re<xbar,xdot> = Re(<bbar,bdot> + <Abar,Adot>)')
        call sparse_destroy(solver)
        call sparse_destroy(plus)
        call sparse_destroy(minus)
    end subroutine check_derivatives

    subroutine blocks(n, special, a)
        ! Block-diagonal 2x2 blocks z [1 c; conj(c) 1] with |z| = 1, |c| = 0.9
        ! (special) or 0.1: the special static second pivot has modulus 0.19
        ! with ||A||_inf = 1.9, the identity matching and unit scaling are
        ! optimal. The block phases z make the matrix non-Hermitian.
        integer, intent(in) :: n, special
        type(csc_z_t), intent(out) :: a
        integer :: rows(2*n), cols(2*n), k, j
        complex(dp) :: values(2*n), c, z
        type(fortsparse_status_t) :: status
        k = 0
        do j = 1, n, 2
            c = 0.1_dp*exp(cmplx(0.0_dp, 0.3_dp*j, dp))
            z = exp(cmplx(0.0_dp, 0.7_dp*j, dp))
            if (j/2 < special) c = 9.0_dp*c
            rows(k+1:k+4) = [j, j+1, j, j+1]
            cols(k+1:k+4) = [j, j, j+1, j+1]
            values(k+1:k+4) = z*[(1.0_dp, 0.0_dp), conjg(c), c, (1.0_dp, 0.0_dp)]
            k = k+4
        end do
        call csc_from_triplet(n, n, rows, cols, values, a, status)
        call check(status_ok(status), 'block matrix construction')
    end subroutine blocks

    subroutine run(a, tau, limit, status, stats)
        ! Factor and solve a known solution; the returned solution is
        ! checked against an independent original-matrix residual.
        type(csc_z_t), intent(in) :: a
        real(dp), intent(in) :: tau
        integer(i8), intent(in) :: limit
        type(fortsparse_status_t), intent(out) :: status
        type(sparse_factor_stats_t), intent(out) :: stats
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: factor_status
        complex(dp), allocatable :: exact(:), b(:), x(:), r(:)
        real(dp), allocatable :: rowsum(:)
        integer :: j, k
        call supernodal(solver)
        solver%options%pivot_perturbation = tau
        solver%options%host_memory_limit_bytes = limit
        call sparse_factor(solver, a, factor_status)
        call sparse_factor_stats(solver, stats, status)
        if (status_ok(factor_status)) then
            allocate(exact(a%ncol), b(a%ncol), x(a%ncol), r(a%ncol), rowsum(a%ncol))
            exact = [(cmplx(1.0_dp+mod(3*j,11), mod(j,4), dp), j=1, a%ncol)]
            b = (0.0_dp, 0.0_dp)
            rowsum = 0.0_dp
            do j = 1, a%ncol
                do k = a%col_ptr(j), a%col_ptr(j+1)-1
                    b(a%row_idx(k)) = b(a%row_idx(k))+a%val(k)*exact(j)
                    rowsum(a%row_idx(k)) = rowsum(a%row_idx(k))+abs(a%val(k))
                end do
            end do
            call sparse_solve(solver, b, x, status)
            call check(status_ok(status), 'solve')
            r = b
            do j = 1, a%ncol
                do k = a%col_ptr(j), a%col_ptr(j+1)-1
                    r(a%row_idx(k)) = r(a%row_idx(k))-a%val(k)*x(j)
                end do
            end do
            call check(maxval(abs(r)) <= 64.0_dp*epsilon(1.0_dp)*(maxval(abs(b)) &
                +maxval(rowsum)*maxval(abs(x))), 'independent original-matrix residual')
        end if
        status = factor_status
        call sparse_destroy(solver)
    end subroutine run

    subroutine check_fallback()
        ! n = 1200: perturbation limit 1. tight perturbs the special pivot
        ! by a relative 1e-6 (refinement converges); twice tight fails the
        ! validation; two special blocks exceed the limit; 0.5 forces it.
        integer, parameter :: n = 1200
        real(dp), parameter :: tight = 0.1_dp*(1.0_dp+1.0e-6_dp)
        type(csc_z_t) :: a
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: stats
        call blocks(n, 1, a)
        call run(a, tight, 0_i8, status, stats)
        call check(status_ok(status), 'one perturbation within the limit')
        call check(stats%perturbed_pivots == 1 .and. .not. stats%lu_fallback, &
            'validation passes after refinement')
        call check(stats%native_factor_family == 'LU_SUPERNODAL', 'supernodal kept')
        call run(a, 2.0_dp*tight, 0_i8, status, stats)
        call check(status_ok(status), 'refinement fallback succeeds')
        call check(stats%lu_fallback .and. stats%lu_fallback_reason == 'REFINEMENT', &
            'fallback reason REFINEMENT')
        call blocks(n, 2, a)
        call run(a, tight, 0_i8, status, stats)
        call check(status_ok(status), 'perturbation fallback succeeds')
        call check(stats%lu_fallback .and. stats%lu_fallback_reason == 'PERTURBATION', &
            'fallback reason PERTURBATION')
        call check(stats%native_factor_family == 'LU_COMPLEX', 'fallback family')
        call check(stats%lu_method_used == FORTSPARSE_LU_SCALAR, 'fallback method')
        call check(stats%perturbed_pivots == 2, 'attempt counters kept after fallback')
        call run(a, 0.5_dp, 0_i8, status, stats)
        call check(status_ok(status) .and. stats%lu_fallback_reason == 'PERTURBATION', &
            'forced fallback')
    end subroutine check_fallback

    subroutine check_singular()
        ! Column 2 is empty: structurally singular matching; the scalar LU
        ! reports SINGULAR.
        type(csc_z_t) :: a
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: stats
        call csc_from_triplet(3, 3, [1, 2, 3], [1, 1, 3], &
            [(1.0_dp, 1.0_dp), (2.0_dp, 0.0_dp), (0.0_dp, 3.0_dp)], a, status)
        call check(status_ok(status), 'singular construction')
        call supernodal(solver)
        call sparse_factor(solver, a, status)
        call check(status%code == FORTSPARSE_SINGULAR, 'singular input reported')
        call sparse_factor_stats(solver, stats, status)
        call check(stats%lu_fallback_reason == 'SINGULAR', 'singular reason recorded')
        call sparse_destroy(solver)
    end subroutine check_singular

    subroutine check_caps(a)
        ! A non-binding cap measures the exact capped peak; one byte less
        ! is refused with required > available; the measured peak admits
        ! the factor again (recovery in a fresh and in the same context).
        complex(dp), intent(in) :: a(:,:)
        type(csc_z_t) :: matrix
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: stats
        type(sparse_solver_t) :: solver
        integer(i8) :: peak
        call to_csc(a, matrix)
        call run(matrix, 1.0e-13_dp, huge(1_i8)/4, status, stats)
        call check(status_ok(status) .and. .not. stats%lu_fallback, 'capped factor')
        peak = stats%host_peak_bytes
        call check(peak > 0_i8, 'exact capped peak')
        call run(matrix, 1.0e-13_dp, peak-1_i8, status, stats)
        call check(status%code == FORTSPARSE_RESOURCE_LIMIT, 'refused below the peak')
        call check(stats%budget_required_bytes > stats%budget_available_bytes, &
            'required and available reported')
        call run(matrix, 1.0e-13_dp, peak, status, stats)
        call check(status_ok(status) .and. .not. stats%lu_fallback, 'fits its own peak')
        call check(stats%host_peak_bytes <= peak, 'capped peak within the cap')
        call supernodal(solver)
        solver%options%host_memory_limit_bytes = peak-1_i8
        call sparse_factor(solver, matrix, status)
        call check(status%code == FORTSPARSE_RESOURCE_LIMIT, 'context refused')
        solver%options%host_memory_limit_bytes = peak
        call sparse_factor(solver, matrix, status)
        call check(status_ok(status), 'same context recovers')
        call sparse_destroy(solver)
    end subroutine check_caps

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (condition) return
        print '(a)', 'FAIL: '//message
        error stop 1
    end subroutine check
end program test_fortsparse_slu_z_api
