program test_fortsparse_slu_api
    ! Public API with lu_method = SUPERNODAL: dense oracles, transpose and
    ! block solves, same-pattern refactor reuse, implicit derivatives
    ! against finite differences and the adjoint identity, context copy
    ! lifetime, explicit refusals and the default ordering selection.
    ! Complex coverage: test_fortsparse_slu_z_api.
    use fortsparse
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    integer :: seed, trial
    real(dp), allocatable :: dense(:,:)

    seed = 31337
    do trial = 1, 6
        call random_dense(8+17*trial, trial, dense)
        call check_public(dense)
    end do
    call random_dense(40, 3, dense)
    call check_derivatives(dense)
    call check_copy(dense)
    call check_refusals(dense)
    call check_ordering(dense)
    print '(a)', 'supernodal LU public API passed'
contains
    real(dp) function uniform()
        seed = modulo(1103515245_i8*seed+12345_i8, 2147483648_i8)
        uniform = real(seed, dp)/2147483648.0_dp
    end function uniform

    subroutine random_dense(n, trial, a)
        integer, intent(in) :: n, trial
        real(dp), allocatable, intent(out) :: a(:,:)
        integer :: perm(n), i, j, k, swap
        allocate(a(n,n))
        a = 0.0_dp
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
            a(i,j) = 2.0_dp*uniform()-1.0_dp
            if (uniform() < 0.85_dp) a(j,i) = 2.0_dp*uniform()-1.0_dp
        end do
        if (mod(trial, 3) == 0) then
            do j = 1, n
                a(j,j) = 0.0_dp
            end do
        end if
        do j = 1, n
            a(perm(j),j) = 2.0_dp+uniform()
        end do
    end subroutine random_dense

    subroutine to_csc(a, matrix)
        real(dp), intent(in) :: a(:,:)
        type(csc_t), intent(out) :: matrix
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: values(:)
        type(fortsparse_status_t) :: status
        integer :: n, i, j, k
        n = size(a,1)
        k = count(a /= 0.0_dp)
        allocate(rows(k), cols(k), values(k))
        k = 0
        do j = 1, n
            do i = 1, n
                if (a(i,j) == 0.0_dp) cycle
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

    subroutine check_public(a)
        real(dp), intent(in) :: a(:,:)
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: stats
        type(csc_t) :: matrix
        real(dp), allocatable :: exact(:), b(:), x(:), block(:,:), result(:,:)
        integer :: n, i, op, j
        logical :: reused
        n = size(a,1)
        call to_csc(a, matrix)
        call supernodal(solver)
        call sparse_factor(solver, matrix, status)
        call check(status_ok(status), 'supernodal factor: '//trim(status%msg))
        call sparse_factor_stats(solver, stats, status)
        call check(stats%native_factor_family == 'LU_SUPERNODAL', 'family reported')
        call check(stats%lu_method_used == FORTSPARSE_LU_SUPERNODAL, 'method reported')
        call check(.not. stats%lu_fallback, 'no fallback on a regular matrix')
        call check(stats%matching_used == FORTSPARSE_MATCHING_WEIGHTED, 'weighted matching')
        call check(stats%perturbed_pivots == 0, 'no perturbation')
        call check(stats%block_pivot_swaps >= 0 .and. stats%pivot_growth > 0.0_dp, &
            'pivot statistics')
        call check(stats%min_pivot > 0.0_dp .and. stats%scaling_log10_range >= 0.0_dp, &
            'pivot and scaling statistics')
        call check(stats%factor_nnz > 0 .and. stats%host_current_bytes > 0, 'sizes reported')
        allocate(exact(n), b(n), x(n), block(n,4), result(n,4))
        exact = [(cos(real(i,dp))+0.5_dp, i=1, n)]
        do op = 0, 2
            if (op == 0) then
                b = matmul(a, exact)
            else
                b = matmul(transpose(a), exact)
            end if
            call sparse_solve_op(solver, b, x, op, status)
            call check(status_ok(status), 'solve operation')
            call check(maxval(abs(x-exact)) <= 1.0e-10_dp*maxval(abs(exact)), &
                'known solution recovered')
            call sparse_factor_stats(solver, stats, status)
            call check(stats%backward_error <= 32.0_dp*epsilon(1.0_dp), &
                'refined backward error gate')
            do j = 1, 4
                block(:,j) = real(j,dp)*b
            end do
            call sparse_solve_block_op(solver, block, result, op, status)
            call check(status_ok(status), 'block operation')
            call check(maxval(abs(result(:,4)-4.0_dp*exact)) <= 4.0e-10_dp*maxval(abs(exact)), &
                'block oracle')
        end do
        ! Same pattern, new values: the analysis is reused.
        matrix%val = 1.5_dp*matrix%val
        call sparse_refactor(solver, matrix, status, reused)
        call check(status_ok(status), 'refactor')
        call check(reused, 'same-pattern refactor reuses the analysis')
        b = 1.5_dp*matmul(a, exact)
        call sparse_solve(solver, b, x, status)
        call check(maxval(abs(x-exact)) <= 1.0e-10_dp*maxval(abs(exact)), 'refactor solve')
        call sparse_factor_stats(solver, stats, status)
        call check(stats%native_factor_family == 'LU_SUPERNODAL', 'refactor family')
        call sparse_destroy(solver)
    end subroutine check_public

    subroutine check_derivatives(a)
        ! x(t) solves (A + t Adot) x = b + t bdot. JVP against a central
        ! difference; VJP against the adjoint pairing <xbar, xdot>.
        real(dp), intent(in) :: a(:,:)
        type(sparse_solver_t) :: solver, plus, minus
        type(fortsparse_status_t) :: status
        type(csc_t) :: matrix, tangent, shifted
        real(dp), allocatable :: b(:), bdot(:), x(:), xdot(:), xp(:), xm(:), adot(:,:)
        real(dp), allocatable :: xbar(:), bbar(:), abar(:)
        real(dp) :: h, pairing, adjoint
        integer :: n, i, j, k
        n = size(a,1)
        call to_csc(a, matrix)
        allocate(adot(n,n), b(n), bdot(n), x(n), xdot(n), xp(n), xm(n), xbar(n), &
            bbar(n), abar(matrix%nnz))
        adot = 0.0_dp
        where (a /= 0.0_dp) adot = sin(3.0_dp*a)
        call to_csc(adot, tangent)
        b = [(sin(real(2*i,dp)), i=1, n)]
        bdot = [(cos(real(i,dp)), i=1, n)]
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
        xbar = [(1.0_dp/real(i,dp), i=1, n)]
        call sparse_solve_vjp(solver, matrix, x, xbar, bbar, abar, status, &
            use_primal_factors=.true.)
        call check(status_ok(status), 'VJP')
        pairing = dot_product(xbar, xdot)
        adjoint = dot_product(bbar, bdot)
        k = 0
        do j = 1, n
            do i = matrix%col_ptr(j), matrix%col_ptr(j+1)-1
                k = k+1
                adjoint = adjoint+abar(k)*adot(matrix%row_idx(i), j)
            end do
        end do
        call check(abs(pairing-adjoint) <= 1.0e-10_dp*(abs(pairing)+1.0_dp), &
            'adjoint identity <xbar,xdot> = <bbar,bdot> + <Abar,Adot>')
        call sparse_destroy(solver)
        call sparse_destroy(plus)
        call sparse_destroy(minus)
    end subroutine check_derivatives

    subroutine check_copy(a)
        real(dp), intent(in) :: a(:,:)
        type(sparse_solver_t), allocatable :: source
        type(sparse_solver_t) :: copy
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: stats
        type(csc_t) :: matrix
        real(dp), allocatable :: exact(:), b(:), x(:)
        integer :: n, i
        n = size(a,1)
        call to_csc(a, matrix)
        allocate(source)
        call supernodal(source)
        call sparse_factor(source, matrix, status)
        call check(status_ok(status), 'copy source factor')
        call sparse_copy(copy, source, status)
        call check(status_ok(status), 'context copy')
        call sparse_destroy(source)
        deallocate(source)
        allocate(exact(n), b(n), x(n))
        exact = [(real(i,dp)/real(n,dp), i=1, n)]
        b = matmul(transpose(a), exact)
        call sparse_solve_op(copy, b, x, FORTSPARSE_OP_TRANSPOSE, status)
        call check(status_ok(status), 'copied factors solve after the source is gone')
        call check(maxval(abs(x-exact)) <= 1.0e-10_dp, 'copied factor oracle')
        call sparse_factor_stats(copy, stats, status)
        call check(stats%native_factor_family == 'LU_SUPERNODAL', 'copied family')
        call sparse_destroy(copy)
    end subroutine check_copy

    subroutine check_refusals(a)
        real(dp), intent(in) :: a(:,:)
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(csc_t) :: matrix
        type(csc_z_t) :: zmatrix
        integer :: n
        n = size(a,1)
        call to_csc(a, matrix)
        call supernodal(solver)
        solver%options%lu_method = 7
        call sparse_factor(solver, matrix, status)
        call check(status%code == FORTSPARSE_INVALID_MATRIX, 'invalid LU method refused')
        call supernodal(solver)
        solver%options%symmetry = FORTSPARSE_SYM_SPD
        call sparse_factor(solver, matrix, status)
        call check(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
            'supernodal LU is GENERAL only')
        solver%options%symmetry = FORTSPARSE_SYM_GENERAL
        zmatrix%nrow = n
        zmatrix%ncol = n
        zmatrix%nnz = matrix%nnz
        zmatrix%col_ptr = matrix%col_ptr
        zmatrix%row_idx = matrix%row_idx
        zmatrix%val = cmplx(matrix%val, 0.5_dp, dp)
        call sparse_factor(solver, zmatrix, status)
        call check(status_ok(status), 'complex supernodal LU factors (stage 8)')
        matrix%val(1) = ieee_value(1.0_dp, ieee_quiet_nan)
        call sparse_factor(solver, matrix, status)
        call check(status%code == FORTSPARSE_INVALID_MATRIX, 'non-finite input refused')
        call sparse_destroy(solver)
    end subroutine check_refusals

    subroutine check_ordering(a)
        ! DEFAULT: METIS ND without a host cap or LOW_MEMORY when METIS is
        ! built, else bounded AMD; an explicit request wins. Each factor
        ! must solve against the dense oracle and report the ordering used.
        real(dp), intent(in) :: a(:,:)
        integer(i8), parameter :: cap = 1000000000_i8
        integer, parameter :: auto = FORTSPARSE_MEMORY_AUTO, low = FORTSPARSE_MEMORY_LOW
        integer, parameter :: default = FORTSPARSE_ORDER_DEFAULT
        integer, parameter :: amd = FORTSPARSE_ORDER_AMD, nd_order = FORTSPARSE_ORDER_ND
        type(csc_t) :: matrix
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        integer :: nd
        nd = amd
        if (native_ordering_available(nd_order)) nd = nd_order
        call to_csc(a, matrix)
        call ordered(a, matrix, default, 0_i8, auto, nd, 'default without cap')
        call ordered(a, matrix, default, cap, auto, amd, 'default under cap')
        call ordered(a, matrix, default, 0_i8, low, amd, 'default with LOW_MEMORY')
        call ordered(a, matrix, amd, 0_i8, auto, amd, 'explicit AMD without cap')
        call ordered(a, matrix, amd, cap, low, amd, 'explicit AMD, cap and LOW')
        if (nd /= nd_order) return
        call ordered(a, matrix, nd_order, 0_i8, auto, nd_order, 'explicit ND')
        ! An explicit ND request under a host cap is still refused.
        call supernodal(solver)
        solver%options%ordering = nd_order
        solver%options%host_memory_limit_bytes = cap
        call sparse_factor(solver, matrix, status)
        call check(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
            'explicit ND under a cap is refused')
        call sparse_destroy(solver)
    end subroutine check_ordering

    subroutine ordered(a, matrix, ordering, cap, policy, expected, label)
        real(dp), intent(in) :: a(:,:)
        type(csc_t), intent(in) :: matrix
        integer, intent(in) :: ordering, policy, expected
        integer(i8), intent(in) :: cap
        character(*), intent(in) :: label
        type(sparse_solver_t) :: solver
        type(sparse_factor_stats_t) :: stats
        type(fortsparse_status_t) :: status
        real(dp) :: known(size(a,1)), b(size(a,1)), x(size(a,1))
        integer :: i
        call supernodal(solver)
        solver%options%ordering = ordering
        solver%options%host_memory_limit_bytes = cap
        solver%options%memory_policy = policy
        solver%refine = .true.
        call sparse_factor(solver, matrix, status)
        call check(status_ok(status), label//': factor')
        call sparse_factor_stats(solver, stats, status)
        call check(stats%native_factor_family == 'LU_SUPERNODAL', label//': family')
        call check(stats%ordering_used == expected, label//': ordering used')
        known = [(sin(real(i,dp)), i = 1, size(a,1))]
        b = matmul(a, known)
        call sparse_solve(solver, b, x, status)
        call check(status_ok(status), label//': solve')
        call check(maxval(abs(x-known)) <= 1.0e-10_dp*maxval(abs(known)), &
            label//': solution')
        call sparse_destroy(solver)
    end subroutine ordered

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (condition) return
        print '(a)', 'FAIL: '//message
        error stop 1
    end subroutine check
end program test_fortsparse_slu_api
