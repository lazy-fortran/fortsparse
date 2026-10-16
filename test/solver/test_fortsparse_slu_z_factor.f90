program test_fortsparse_slu_z_factor
    ! Complex supernodal LU numerics against independent oracles: LAPACK
    ! zgesv, original-matrix residuals of known solutions for A, A^T and
    ! A^H, a dense complex restricted-pivot emulation written here (exchange
    ! and perturbation counts, solutions), the bilinear transpose identity
    ! c^T A^-1 b = (A^-T c)^T b, the sesquilinear adjoint identity
    ! <c, A^-1 b> = <A^-H c, b>, block RHS, zero diagonals that need the
    ! weighted matching, value-type refusals and the exact 16-byte plan.
    use fortsparse_kinds, only: dp, i8
    use fortsparse_memory, only: memory_budget_t
    use fortsparse_status, only: fortsparse_status_t, FORTSPARSE_OK, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_RESOURCE_LIMIT, FORTSPARSE_SINGULAR
    use fortsparse_native_slu, only: slu_factor_t, SLU_PHASE_COUNT
    implicit none
    interface
        subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import dp
            integer, intent(in) :: n, nrhs, lda, ldb
            complex(dp), intent(inout) :: a(lda,*), b(ldb,*)
            integer, intent(out) :: ipiv(*), info
        end subroutine zgesv
    end interface
    integer :: seed, trial, total_swaps
    complex(dp), allocatable :: dense(:,:)

    seed = 7717
    total_swaps = 0
    do trial = 1, 12
        call random_dense(5+mod(29*trial, 120), trial, dense)
        call check_matrix(dense, .false., total_swaps)
        call check_matrix(dense, .true., total_swaps)
    end do
    call check(total_swaps > 0, 'random cases exercise in-block exchanges')
    ! n = 300 gives 128-wide dense supernodes: the ZTRSM panel path.
    call dense_case(300, dense)
    call check_matrix(dense, .true., total_swaps)
    call check_zero_pivot()
    call check_refusals()
    call random_dense(60, 3, dense)
    call check_caps(dense)
    print '(a,i0,a)', 'complex supernodal LU factor: ', total_swaps, ' block exchanges'
contains
    real(dp) function uniform()
        seed = int(modulo(1103515245_i8*seed+12345_i8, 2147483648_i8))
        uniform = real(seed, dp)/2147483648.0_dp
    end function uniform

    complex(dp) function random_z()
        random_z = cmplx(2.0_dp*uniform()-1.0_dp, 2.0_dp*uniform()-1.0_dp, dp)
    end function random_z

    subroutine random_dense(n, trial, a)
        ! Near-symmetric random pattern of genuinely complex (non-Hermitian)
        ! entries; a random perfect matching carries the larger entries and
        ! every third case has no diagonal at all. Wide column scales.
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
            a(perm(j),j) = (1.0_dp+uniform())*exp(cmplx(0.0_dp, 6.0_dp*uniform(), dp))
        end do
        do j = 1, n
            a(:,j) = a(:,j)*10.0_dp**(int(6*uniform())-3)
        end do
    end subroutine random_dense

    subroutine dense_case(n, a)
        integer, intent(in) :: n
        complex(dp), allocatable, intent(out) :: a(:,:)
        integer :: i, j
        allocate(a(n,n))
        do j = 1, n
            do i = 1, n
                a(i,j) = random_z()
            end do
            a(j,j) = a(j,j)+(4.0_dp, 2.0_dp)
        end do
    end subroutine dense_case

    subroutine to_csc(a, ptr, rows, values)
        complex(dp), intent(in) :: a(:,:)
        integer, allocatable, intent(out) :: ptr(:), rows(:)
        complex(dp), allocatable, intent(out) :: values(:)
        integer :: n, i, j, k
        n = size(a,1)
        k = count(a /= (0.0_dp, 0.0_dp))
        allocate(ptr(n+1), rows(k), values(k))
        k = 0
        ptr(1) = 1
        do j = 1, n
            do i = 1, n
                if (a(i,j) == (0.0_dp, 0.0_dp)) cycle
                k = k+1
                rows(k) = i
                values(k) = a(i,j)
            end do
            ptr(j+1) = k+1
        end do
    end subroutine to_csc

    real(dp) function backward(a, x, b) result(error)
        complex(dp), intent(in) :: a(:,:), x(:), b(:)
        error = maxval(abs(b-matmul(a,x)))/(maxval(sum(abs(a),2))*maxval(abs(x)) &
            +maxval(abs(b)))
    end function backward

    subroutine check_matrix(a, blas, total_swaps)
        complex(dp), intent(in) :: a(:,:)
        logical, intent(in) :: blas
        integer, intent(inout) :: total_swaps
        type(slu_factor_t) :: lu
        type(fortsparse_status_t) :: status
        integer, allocatable :: ptr(:), rows(:), ipiv(:)
        complex(dp), allocatable :: values(:), exact(:), b(:), x(:), ref(:,:), rhs(:,:)
        complex(dp), allocatable :: emulated(:), lapack(:)
        integer :: n, info, i, op, swaps, perturbed
        n = size(a,1)
        call to_csc(a, ptr, rows, values)
        call lu%analyze(n, ptr, rows, values, status, force_blas=blas)
        call check(status%code == FORTSPARSE_OK, 'analysis: '//trim(status%msg))
        call check(lu%is_complex(), 'complex analysis')
        call check(lu%plan_bytes == lu%memory_bytes(), 'exact 16-byte plan')
        call lu%factor(ptr, rows, values, status)
        call check(status%code == FORTSPARSE_OK, 'factor: '//trim(status%msg))
        allocate(exact(n), b(n), x(n), ipiv(n))
        exact = [(cmplx(sin(real(i,dp))+0.25_dp, cos(real(2*i,dp)), dp), i=1, n)]
        do op = 0, 2
            b = matmul(operator_matrix(a, op), exact)
            x = b
            call lu%solve(x, status, op)
            call check(status%code == FORTSPARSE_OK, 'solve')
            call check(backward(operator_matrix(a, op), x, b) < 1.0e-13_dp, &
                'original-matrix backward error (N, T, H)')
            ref = operator_matrix(a, op)
            rhs = reshape(b, [n,1])
            call zgesv(n, 1, ref, n, ipiv, rhs, n, info)
            call check(info == 0, 'LAPACK zgesv reference')
            call check(maxval(abs(x-rhs(:,1))) <= 1.0e-8_dp*maxval(abs(rhs(:,1))), &
                'solution agrees with zgesv')
            if (op == 0) lapack = rhs(:,1)
        end do
        call check_identities(lu, n)
        call emulate(lu, a, swaps, perturbed, emulated)
        call check(swaps == lu%block_swaps, 'block exchanges equal the emulation')
        call check(perturbed == lu%perturbed_pivots, 'perturbations equal the emulation')
        b = matmul(a, exact)
        x = b
        call lu%solve(x, status)
        ! Rounding differs (left- vs right-looking sums); the zgesv spread
        ! measures the conditioning of the case.
        call check(maxval(abs(x-emulated)) <= 1.0e-10_dp*maxval(abs(emulated)) &
            +100.0_dp*maxval(abs(lapack-emulated)), 'solution equals the emulated solve')
        call check(lu%max_growth > 0.0_dp .and. lu%min_pivot > 0.0_dp, 'growth statistics')
        if (n >= 300 .and. blas) call check(wide_panel(lu), 'ZTRSM-sized panel present')
        total_swaps = total_swaps+lu%block_swaps
        call lu%free()
    end subroutine check_matrix

    logical function wide_panel(lu)
        ! A supernode with width >= 96 and >= 32 off-diagonal rows.
        type(slu_factor_t), intent(in) :: lu
        integer, allocatable :: rows(:)
        integer :: s, first, width, m
        wide_panel = .false.
        do s = 1, lu%node_count
            call lu%node_structure(s, first, width, m, rows)
            if (width >= 96 .and. m-width >= 32) wide_panel = .true.
        end do
    end function wide_panel

    function operator_matrix(a, op) result(m)
        complex(dp), intent(in) :: a(:,:)
        integer, intent(in) :: op
        complex(dp), allocatable :: m(:,:)
        select case (op)
        case (0)
            m = a
        case (1)
            m = transpose(a)
        case default
            m = conjg(transpose(a))
        end select
    end function operator_matrix

    subroutine check_identities(lu, n)
        ! Bilinear transpose and sesquilinear adjoint identities with
        ! independent right-hand sides, and block columns equal to vector
        ! solves for every operation.
        type(slu_factor_t), intent(inout) :: lu
        integer, intent(in) :: n
        type(fortsparse_status_t) :: status
        complex(dp) :: b(n), c(n), x(n), y(n), block(n,3)
        integer :: i, op
        b = [(cmplx(cos(real(3*i,dp)), sin(real(i,dp)), dp), i=1, n)]
        c = [(cmplx(sin(real(5*i,dp)), 0.5_dp*cos(real(7*i,dp)), dp), i=1, n)]
        x = b
        call lu%solve(x, status, 0)
        y = c
        call lu%solve(y, status, 1)
        call check(abs(sum(c*x)-sum(y*b)) <= 1.0e-10_dp*(sum(abs(c*x))+sum(abs(y*b))+1.0_dp), &
            'transpose identity c^T A^-1 b = (A^-T c)^T b')
        y = c
        call lu%solve(y, status, 2)
        call check(abs(dot_product(c,x)-dot_product(y,b)) <= &
            1.0e-10_dp*(sum(abs(c*x))+sum(abs(y*b))+1.0_dp), &
            'adjoint identity <c, A^-1 b> = <A^-H c, b>')
        do op = 0, 2
            block(:,1) = b
            block(:,2) = c
            block(:,3) = b-(0.0_dp, 2.0_dp)*c
            call lu%solve_block(block, status, op)
            call check(status%code == FORTSPARSE_OK, 'block solve')
            y = b-(0.0_dp, 2.0_dp)*c
            call lu%solve(y, status, op)
            call check(all(block(:,3) == y), 'block column equals vector solve')
        end do
    end subroutine check_identities

    subroutine emulate(lu, a, swaps, perturbed, x)
        ! Right-looking dense LU of M = Dr Pm A Dc in the final order with
        ! partial pivoting restricted to each supernode's rows and the same
        ! perturbation rule tau*p/|p|; solves A x = A*exact with it.
        type(slu_factor_t), intent(in) :: lu
        complex(dp), intent(in) :: a(:,:)
        integer, intent(out) :: swaps, perturbed
        complex(dp), allocatable, intent(out) :: x(:)
        integer, allocatable :: row_old(:), col_old(:), node_end(:), rows(:), rowperm(:)
        real(dp), allocatable :: rs(:), cs(:)
        complex(dp), allocatable :: f(:,:), t(:), saved(:), exact(:)
        integer :: n, k, i, r, s, first, width, m, swap_row
        real(dp) :: tau_abs, best
        n = size(a,1)
        call lu%ordering(row_old, col_old, rs, cs)
        allocate(node_end(n), f(n,n), t(n), saved(n), rowperm(n), x(n), exact(n))
        do s = 1, lu%node_count
            call lu%node_structure(s, first, width, m, rows)
            node_end(first:first+width-1) = first+width-1
        end do
        do k = 1, n
            do i = 1, n
                f(i,k) = rs(i)*a(row_old(i), col_old(k))*cs(k)
            end do
        end do
        tau_abs = 1.0e-13_dp*maxval(sum(abs(f),2))
        rowperm = [(i, i=1, n)]
        swaps = 0
        perturbed = 0
        do k = 1, n
            r = k
            best = abs(f(k,k))
            do i = k+1, node_end(k)
                if (abs(f(i,k)) > best) then
                    best = abs(f(i,k))
                    r = i
                end if
            end do
            if (r /= k) then
                saved = f(k,:)
                f(k,:) = f(r,:)
                f(r,:) = saved
                swap_row = rowperm(k)
                rowperm(k) = rowperm(r)
                rowperm(r) = swap_row
                swaps = swaps+1
            end if
            if (abs(f(k,k)) < tau_abs) then
                if (abs(f(k,k)) > 0.0_dp) then
                    f(k,k) = tau_abs*f(k,k)/abs(f(k,k))
                else
                    f(k,k) = cmplx(tau_abs, 0.0_dp, dp)
                end if
                perturbed = perturbed+1
            end if
            f(k+1:,k) = f(k+1:,k)/f(k,k)
            do i = k+1, n
                f(k+1:,i) = f(k+1:,i)-f(k+1:,k)*f(k,i)
            end do
        end do
        exact = [(cmplx(sin(real(i,dp))+0.25_dp, cos(real(2*i,dp)), dp), i=1, n)]
        saved = matmul(a, exact)
        do i = 1, n
            t(i) = rs(rowperm(i))*saved(row_old(rowperm(i)))
        end do
        do k = 1, n
            t(k+1:) = t(k+1:)-f(k+1:,k)*t(k)
        end do
        do k = n, 1, -1
            t(k) = t(k)/f(k,k)
            t(:k-1) = t(:k-1)-f(:k-1,k)*t(k)
        end do
        do k = 1, n
            x(col_old(k)) = cs(k)*t(k)
        end do
    end subroutine emulate

    subroutine check_zero_pivot()
        ! One 2x2 supernode, identity matching optimal, second static pivot
        ! exactly zero: perturbed once to +tau*||M||.
        type(slu_factor_t) :: lu
        type(fortsparse_status_t) :: status
        integer :: ptr(3), rows(4)
        complex(dp) :: values(4), x(2)
        ptr = [1, 3, 5]
        rows = [1, 2, 1, 2]
        values = [(0.0_dp, 1.0_dp), (0.0_dp, 1.0_dp), (1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp)]
        call lu%analyze(2, ptr, rows, values, status)
        call check(status%code == FORTSPARSE_OK, 'zero-pivot analysis')
        call check(lu%node_count == 1, 'zero-pivot matrix is one supernode')
        call lu%factor(ptr, rows, values, status)
        call check(status%code == FORTSPARSE_OK, 'perturbed factor succeeds')
        call check(lu%perturbed_pivots == 1 .and. lu%block_swaps == 0, &
            'exactly one perturbed pivot, no exchange')
        call check(abs(lu%min_pivot-2.0e-13_dp) <= 1.0e-15_dp, &
            'perturbed pivot magnitude tau*||M||')
        x = [(0.0_dp, 0.0_dp), (1.0_dp, 0.0_dp)]
        call lu%solve(x, status)
        call check(abs(x(2)) >= 1.0e12_dp, 'solve uses the perturbed pivot')
    end subroutine check_zero_pivot

    subroutine check_refusals()
        ! Value-type mismatches, a pattern outside the structure and a
        ! structurally singular matrix.
        type(slu_factor_t) :: lu
        type(fortsparse_status_t) :: status
        integer :: ptr(4), rows(3), ptr2(4), rows2(4), ptr3(4), rows3(3)
        complex(dp) :: values(3), values2(4), x(3)
        real(dp) :: rvalues(3), rx(3)
        ptr = [1, 2, 3, 4]
        rows = [1, 2, 3]
        values = [(2.0_dp, 1.0_dp), (0.0_dp, 3.0_dp), (4.0_dp, -1.0_dp)]
        rvalues = [2.0_dp, 3.0_dp, 4.0_dp]
        call lu%analyze(3, ptr, rows, values, status)
        call check(status%code == FORTSPARSE_OK, 'diagonal analysis')
        call lu%factor(ptr, rows, rvalues, status)
        call check(status%code == FORTSPARSE_INVALID_MATRIX, 'real values refused')
        ptr2 = [1, 3, 4, 5]
        rows2 = [1, 3, 2, 3]
        values2 = [(2.0_dp, 0.0_dp), (1.0_dp, 1.0_dp), (3.0_dp, 0.0_dp), (4.0_dp, 0.0_dp)]
        call lu%factor(ptr2, rows2, values2, status)
        call check(status%code == FORTSPARSE_INVALID_MATRIX, 'entry outside structure refused')
        call check(.not. lu%is_factored(), 'refused factor is not usable')
        call lu%factor(ptr, rows, values, status)
        call check(status%code == FORTSPARSE_OK, 'same structure refactors')
        rx = 1.0_dp
        call lu%solve(rx, status)
        call check(status%code == FORTSPARSE_INVALID_MATRIX, 'real RHS refused')
        x = [(2.0_dp, 1.0_dp), (0.0_dp, 3.0_dp), (4.0_dp, -1.0_dp)]
        call lu%solve(x, status, 2)
        call check(maxval(abs(x-[(0.6_dp,0.8_dp), (-1.0_dp,0.0_dp), &
            cmplx(15.0_dp/17.0_dp,-8.0_dp/17.0_dp,dp)])) <= 1.0e-15_dp, 'diagonal adjoint solve')
        ptr3 = [1, 3, 3, 4]
        rows3 = [1, 2, 3]
        call lu%analyze(3, ptr3, rows3, values, status)
        call check(status%code == FORTSPARSE_SINGULAR, 'structurally singular matching')
        call check(.not. lu%is_analyzed(), 'singular analysis leaves no factor')
    end subroutine check_refusals

    subroutine check_caps(a)
        ! Exact ledger: an unlimited analysis defines the peak; a cap one
        ! byte short is refused with required > available and a restored
        ! ledger; each injected phase refusal is recoverable.
        complex(dp), intent(in) :: a(:,:)
        type(slu_factor_t) :: lu, real_lu
        type(fortsparse_status_t) :: status
        type(memory_budget_t) :: ledger
        integer, allocatable :: ptr(:), rows(:)
        complex(dp), allocatable :: values(:)
        integer(i8) :: peak
        integer :: phase
        call to_csc(a, ptr, rows, values)
        ledger%limit = huge(1_i8)/4
        call lu%analyze(size(a,1), ptr, rows, values, status, budget=ledger)
        call check(status%code == FORTSPARSE_OK, 'capped analysis')
        call check(ledger%current == lu%memory_bytes(), 'ledger holds the retained layout')
        peak = lu%peak_memory_bytes
        call real_lu%analyze(size(a,1), ptr, rows, real(values), status)
        call check(lu%memory_bytes() > real_lu%memory_bytes(), 'complex layout is larger')
        call lu%free(ledger, status)
        call check(ledger%current == 0_i8, 'free releases the ledger')
        ledger = memory_budget_t(limit=peak-1_i8)
        call lu%analyze(size(a,1), ptr, rows, values, status, budget=ledger)
        call check(status%code == FORTSPARSE_RESOURCE_LIMIT, 'cap one byte short refused')
        call check(lu%failed_required_bytes > lu%failed_available_bytes, &
            'required and available reported')
        call check(ledger%current == 0_i8, 'refusal restores the ledger')
        do phase = 1, SLU_PHASE_COUNT
            ledger = memory_budget_t(limit=huge(1_i8)/4)
            call lu%analyze(size(a,1), ptr, rows, values, status, budget=ledger, &
                inject_phase=phase)
            call check(status%code == FORTSPARSE_RESOURCE_LIMIT, 'injected phase refused')
            call check(lu%failed_phase == phase, 'failing phase named')
            call check(ledger%current == 0_i8, 'injected refusal restores the ledger')
        end do
        ledger = memory_budget_t(limit=peak)
        call lu%analyze(size(a,1), ptr, rows, values, status, budget=ledger)
        call check(status%code == FORTSPARSE_OK, 'recovers at its measured peak')
        call lu%factor(ptr, rows, values, status)
        call check(status%code == FORTSPARSE_OK, 'factor after recovery')
    end subroutine check_caps

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (condition) return
        print '(a)', 'FAIL: '//message
        error stop 1
    end subroutine check
end program test_fortsparse_slu_z_factor
