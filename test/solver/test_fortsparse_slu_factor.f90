program test_fortsparse_slu_factor
    ! Serial supernodal LU numerics against independent oracles: LAPACK dgesv,
    ! original-matrix residuals of known solutions, a dense restricted-pivot
    ! emulation written here (exchange and perturbation counts, solutions),
    ! the transpose identity c^T A^-1 b = (A^-T c)^T b and block RHS.
    use fortsparse_kinds, only: dp, i8
    use fortsparse_status, only: fortsparse_status_t, FORTSPARSE_OK, &
        FORTSPARSE_INVALID_MATRIX
    use fortsparse_native_slu, only: slu_factor_t, SLU_ORDER_AMD, SLU_ORDER_ND
    implicit none
    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda,*), b(ldb,*)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface
    integer :: seed, trial, total_swaps
    real(dp), allocatable :: dense(:,:)

    seed = 977
    total_swaps = 0
    do trial = 1, 24
        call random_dense(5+mod(37*trial, 196), trial, dense)
        call check_matrix(dense, .false., total_swaps)
        call check_matrix(dense, .true., total_swaps)
    end do
    call check(total_swaps > 0, 'random cases exercise in-block exchanges')
    call block_diagonal(dense)
    call check_matrix(dense, .false., total_swaps)
    call dense_case(260, dense)
    call check_matrix(dense, .true., total_swaps)
    call check_zero_pivot([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], 1.0_dp)
    call check_zero_pivot([1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp-2.0_dp**(-50)], -1.0_dp)
    call check_outside_pattern()
    print '(a,i0,a)', 'supernodal LU factor: ', total_swaps, ' block exchanges verified'
contains
    real(dp) function uniform()
        seed = modulo(1103515245_i8*seed+12345_i8, 2147483648_i8)
        uniform = real(seed, dp)/2147483648.0_dp
    end function uniform

    subroutine random_dense(n, trial, a)
        ! Near-symmetric random pattern; a random perfect matching carries
        ! the larger entries, and every third case has no diagonal at all.
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
            a(perm(j),j) = (1.0_dp+uniform())*merge(1.0_dp, -1.0_dp, uniform() < 0.5_dp)
        end do
        ! Wide column scales exercise the matching scaling.
        do j = 1, n
            a(:,j) = a(:,j)*10.0_dp**(int(6*uniform())-3)
        end do
    end subroutine random_dense

    subroutine block_diagonal(a)
        ! Three disconnected random blocks under a random symmetric relabelling.
        real(dp), allocatable, intent(out) :: a(:,:)
        real(dp), allocatable :: b(:,:), c(:,:)
        integer :: perm(75), i, j, swap, offset, k
        allocate(c(75,75))
        c = 0.0_dp
        offset = 0
        do k = 1, 3
            call random_dense(15+5*k, 3*k, b)
            c(offset+1:offset+size(b,1), offset+1:offset+size(b,1)) = b
            offset = offset+size(b,1)
        end do
        perm = [(i, i=1, 75)]
        do i = 75, 2, -1
            j = 1+int(uniform()*i)
            swap = perm(i)
            perm(i) = perm(j)
            perm(j) = swap
        end do
        allocate(a(75,75))
        a = c(perm, perm)
    end subroutine block_diagonal

    subroutine dense_case(n, a)
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: a(:,:)
        integer :: i, j
        allocate(a(n,n))
        do j = 1, n
            do i = 1, n
                a(i,j) = 2.0_dp*uniform()-1.0_dp
            end do
            a(j,j) = a(j,j)+4.0_dp
        end do
    end subroutine dense_case

    subroutine to_csc(a, ptr, rows, values)
        real(dp), intent(in) :: a(:,:)
        integer, allocatable, intent(out) :: ptr(:), rows(:)
        real(dp), allocatable, intent(out) :: values(:)
        integer :: n, i, j, k
        n = size(a,1)
        allocate(ptr(n+1), rows(count(a /= 0.0_dp)), values(count(a /= 0.0_dp)))
        k = 0
        ptr(1) = 1
        do j = 1, n
            do i = 1, n
                if (a(i,j) == 0.0_dp) cycle
                k = k+1
                rows(k) = i
                values(k) = a(i,j)
            end do
            ptr(j+1) = k+1
        end do
    end subroutine to_csc

    real(dp) function backward(a, x, b) result(error)
        real(dp), intent(in) :: a(:,:), x(:), b(:)
        error = maxval(abs(b-matmul(a,x)))/(maxval(sum(abs(a),2))*maxval(abs(x)) &
            +maxval(abs(b)))
    end function backward

    subroutine check_matrix(a, blas, total_swaps)
        real(dp), intent(in) :: a(:,:)
        logical, intent(in) :: blas
        integer, intent(inout) :: total_swaps
        type(slu_factor_t) :: lu
        type(fortsparse_status_t) :: status
        integer, allocatable :: ptr(:), rows(:), ipiv(:)
        real(dp), allocatable :: values(:), exact(:), b(:), x(:), ref(:,:), lapack(:,:)
        real(dp), allocatable :: c(:), y(:), block(:,:), emulated(:)
        integer :: n, info, i, swaps, perturbed
        n = size(a,1)
        call to_csc(a, ptr, rows, values)
        call lu%analyze(n, ptr, rows, values, status, force_blas=blas)
        call check(status%code == FORTSPARSE_OK, 'analysis: '//trim(status%msg))
        call lu%factor(ptr, rows, values, status)
        call check(status%code == FORTSPARSE_OK, 'factor: '//trim(status%msg))
        allocate(exact(n), b(n), x(n), c(n), y(n), block(n,3), ipiv(n))
        exact = [(sin(real(i,dp))+0.25_dp, i=1, n)]
        b = matmul(a, exact)
        x = b
        call lu%solve(x, status)
        call check(status%code == FORTSPARSE_OK, 'normal solve')
        call check(backward(a, x, b) < 1.0e-13_dp, 'normal original-matrix backward error')
        ref = a
        lapack = reshape(b, [n,1])
        call dgesv(n, 1, ref, n, ipiv, lapack, n, info)
        call check(info == 0, 'LAPACK reference')
        call check(maxval(abs(x-lapack(:,1))) <= 1.0e-8_dp*maxval(abs(lapack(:,1))), &
            'solution agrees with dgesv')
        b = matmul(transpose(a), exact)
        x = b
        call lu%solve(x, status, 1)
        call check(status%code == FORTSPARSE_OK, 'transpose solve')
        call check(backward(transpose(a), x, b) < 1.0e-13_dp, &
            'transpose original-matrix backward error')
        ! Transpose identity with independent right-hand sides.
        b = [(cos(real(3*i,dp)), i=1, n)]
        c = [(sin(real(5*i,dp)), i=1, n)]
        x = b
        call lu%solve(x, status, 0)
        y = c
        call lu%solve(y, status, 2)
        call check(abs(dot_product(c,x)-dot_product(y,b)) <= &
            1.0e-10_dp*(sum(abs(c*x))+sum(abs(y*b))+1.0_dp), 'transpose identity')
        block(:,1) = b
        block(:,2) = c
        block(:,3) = b-c
        call lu%solve_block(block, status, 1)
        call check(status%code == FORTSPARSE_OK, 'block solve')
        y = b
        call lu%solve(y, status, 1)
        call check(all(block(:,1) == y), 'block column equals vector solve')
        y = b-c
        call lu%solve(y, status, 1)
        call check(all(block(:,3) == y), 'third block column equals vector solve')
        ! Dense restricted-pivot emulation in the same order and blocks.
        call emulate(lu, a, swaps, perturbed, emulated)
        call check(swaps == lu%block_swaps, 'block exchanges equal the emulation')
        call check(perturbed == lu%perturbed_pivots, 'perturbations equal the emulation')
        b = matmul(a, exact)
        x = b
        call lu%solve(x, status)
        ! Rounding differs (left- vs right-looking sums); the LAPACK spread
        ! measures the conditioning of the case.
        call check(maxval(abs(x-emulated)) <= 1.0e-10_dp*maxval(abs(emulated)) &
            +100.0_dp*maxval(abs(lapack(:,1)-emulated)), &
            'solution equals the emulated factor solve')
        call check(lu%max_growth > 0.0_dp .and. lu%min_pivot > 0.0_dp, 'growth statistics')
        total_swaps = total_swaps+lu%block_swaps
        call lu%free()
    end subroutine check_matrix

    subroutine emulate(lu, a, swaps, perturbed, x)
        ! Right-looking dense LU of M = Dr Pm A Dc in the final order with
        ! partial pivoting restricted to each supernode's rows and the same
        ! perturbation rule; solves A x = A*exact with it.
        type(slu_factor_t), intent(in) :: lu
        real(dp), intent(in) :: a(:,:)
        integer, intent(out) :: swaps, perturbed
        real(dp), allocatable, intent(out) :: x(:)
        integer, allocatable :: row_old(:), col_old(:), node_end(:), rows(:), rowperm(:)
        real(dp), allocatable :: rs(:), cs(:), f(:,:), t(:), saved(:), exact(:)
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
                f(k,k) = merge(tau_abs, -tau_abs, f(k,k) >= 0.0_dp)
                perturbed = perturbed+1
            end if
            f(k+1:,k) = f(k+1:,k)/f(k,k)
            do i = k+1, n
                f(k+1:,i) = f(k+1:,i)-f(k+1:,k)*f(k,i)
            end do
        end do
        exact = [(sin(real(i,dp))+0.25_dp, i=1, n)]
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

    subroutine check_zero_pivot(values, sign_expected)
        ! A 2x2 supernode whose identity matching is optimal and whose second
        ! pivot is 0 (all ones) or -2^-50 (columns [1 1], [-1 -1-2^-50]):
        ! below tau*||M||, it becomes sign(p)*tau*||M|| with sign(0) = +.
        real(dp), intent(in) :: values(4), sign_expected
        type(slu_factor_t) :: lu
        type(fortsparse_status_t) :: status
        integer :: ptr(3), rows(4)
        real(dp) :: x(2)
        ptr = [1, 3, 5]
        rows = [1, 2, 1, 2]
        call lu%analyze(2, ptr, rows, values, status)
        call check(status%code == FORTSPARSE_OK, 'zero-pivot analysis')
        call check(lu%node_count == 1, 'zero-pivot matrix is one supernode')
        call lu%factor(ptr, rows, values, status)
        call check(status%code == FORTSPARSE_OK, 'perturbed factor succeeds')
        call check(lu%perturbed_pivots == 1, 'exactly one perturbed pivot')
        call check(lu%block_swaps == 0, 'no exchange on equal magnitudes')
        call check(abs(lu%min_pivot-2.0e-13_dp) <= 1.0e-15_dp, &
            'perturbed pivot magnitude tau*||M||')
        x = [0.0_dp, 1.0_dp]
        call lu%solve(x, status)
        call check(sign(1.0_dp, x(2)) == sign_expected, 'perturbation keeps the sign, sign(0)=+')
        call check(abs(x(2)) >= 1.0e12_dp, 'solve uses the perturbed pivot')
    end subroutine check_zero_pivot

    subroutine check_outside_pattern()
        type(slu_factor_t) :: lu
        type(fortsparse_status_t) :: status
        integer :: ptr(4), rows(3), ptr2(4), rows2(4)
        real(dp) :: values(3), values2(4)
        ptr = [1, 2, 3, 4]
        rows = [1, 2, 3]
        values = [2.0_dp, 3.0_dp, 4.0_dp]
        call lu%analyze(3, ptr, rows, values, status)
        call check(status%code == FORTSPARSE_OK, 'diagonal analysis')
        ptr2 = [1, 3, 4, 5]
        rows2 = [1, 3, 2, 3]
        values2 = [2.0_dp, 1.0_dp, 3.0_dp, 4.0_dp]
        call lu%factor(ptr2, rows2, values2, status)
        call check(status%code == FORTSPARSE_INVALID_MATRIX, 'entry outside structure refused')
        call check(.not. lu%is_factored(), 'refused factor is not usable')
        call lu%factor(ptr, rows, values, status)
        call check(status%code == FORTSPARSE_OK, 'same structure refactors')
    end subroutine check_outside_pattern

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (condition) return
        print '(a)', 'FAIL: '//message
        error stop 1
    end subroutine check
end program test_fortsparse_slu_factor
