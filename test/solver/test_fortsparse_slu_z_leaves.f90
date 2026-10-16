program test_fortsparse_slu_z_leaves
    ! Compiled complex supernodal LU kernels against independent dense
    ! references: intrinsic complex division (including scales where |b|**2
    ! overflows), P A = L U reconstruction with textbook partial pivoting
    ! restricted to the block, matmul products for the triangular solves and
    ! the Schur update (plain transposes), the block solves, exchanges and
    ! the scatter/gather leaves, and the perturbation count and phase.
    use fortsparse_kinds, only: dp, i8
    use fortsparse_slu_z_leaves, only: slu_z_divide, slu_z_update, slu_z_scale, &
        slu_z_block_factor, slu_z_trsm_upper, slu_z_trsm_unit_lower_t, &
        slu_z_schur_update, slu_z_apply_swaps, slu_z_block_lower, slu_z_block_upper, &
        slu_z_block_upper_t, slu_z_block_lower_t, slu_z_scatter_mapped_update, &
        slu_z_scatter_update, slu_z_gather_update
    implicit none
    real(dp), parameter :: eps = epsilon(1.0_dp)
    integer :: seed

    seed = 4242
    call check_scalar_leaves()
    call check_block_factor(7, 11)
    call check_block_factor(33, 33)
    call check_perturbation()
    call check_triangular(9, 13)
    call check_schur(17, 6, 5)
    call check_block_solves(12)
    call check_scatter(15)
    print '(a)', 'complex supernodal LU leaves passed'
contains
    real(dp) function uniform()
        seed = int(modulo(1103515245_i8*seed+12345_i8, 2147483648_i8))
        uniform = real(seed, dp)/2147483648.0_dp
    end function uniform

    complex(dp) function random_z()
        random_z = cmplx(2.0_dp*uniform()-1.0_dp, 2.0_dp*uniform()-1.0_dp, dp)
    end function random_z

    subroutine fill(a)
        complex(dp), intent(out) :: a(:,:)
        integer :: i, j
        do j = 1, size(a,2)
            do i = 1, size(a,1)
                a(i,j) = random_z()
            end do
        end do
    end subroutine fill

    subroutine check_scalar_leaves()
        complex(dp) :: a, b, c, q, a0, b0
        integer :: k
        real(dp) :: scales(4)
        ! Common power-of-two scales: the quotient is unchanged, while a naive
        ! |b|**2 would overflow (2**1000) or underflow (2**-1000).
        scales = [1.0_dp, 2.0_dp**1000, 2.0_dp**(-1000), 2.0_dp**900]
        do k = 1, 200
            a0 = random_z()
            b0 = random_z()+(0.1_dp, 0.0_dp)
            a = a0*scales(1+mod(k,4))
            b = b0*scales(1+mod(k,4))
            q = slu_z_divide(a, b)
            call check(abs(q-a0/b0) <= 8.0_dp*eps*abs(a0/b0), 'complex quotient')
        end do
        a = (1.0e300_dp, 3.0e300_dp)
        b = (2.0e300_dp, -1.0e300_dp)
        q = slu_z_divide(a, b)
        call check(abs(q-(-0.2_dp, 1.4_dp)) <= 8.0_dp*eps, 'quotient without overflow')
        call check(slu_z_divide((0.0_dp,0.0_dp), (0.0_dp,2.0_dp)) == (0.0_dp,0.0_dp), &
            'zero numerator')
        a = (1.0_dp, 2.0_dp)
        b = (3.0_dp, -1.0_dp)
        c = (0.5_dp, 4.0_dp)
        call check(abs(slu_z_update(a, b, c)-(a-b*c)) <= 4.0_dp*eps*abs(b*c), &
            'complex update a - b c (not conjugated)')
        call check(abs(slu_z_scale(2.0_dp, b, 0.25_dp)-(1.5_dp, -0.5_dp)) <= eps, &
            'real scaling of a complex value')
    end subroutine check_scalar_leaves

    subroutine check_block_factor(w, lda)
        ! Rows w+1..lda below the block must stay untouched.
        integer, intent(in) :: w, lda
        complex(dp) :: a(lda,w), original(lda,w), ref(w,w), l(w,w), u(w,w)
        integer :: ipiv(w), ref_piv(w), i, k, r, perturbed, swaps, ref_swaps
        real(dp) :: min_pivot, max_u, best
        complex(dp) :: row(w)
        call fill(a)
        original = a
        perturbed = 0
        swaps = 0
        min_pivot = huge(1.0_dp)
        max_u = 0.0_dp
        call slu_z_block_factor(a, lda, w, 0.0_dp, ipiv, perturbed, swaps, &
            min_pivot, max_u)
        call check(all(a(w+1:,:) == original(w+1:,:)), 'rows below the block untouched')
        ! Independent textbook partial pivoting on the block.
        ref = original(1:w,1:w)
        ref_swaps = 0
        do k = 1, w
            r = k
            best = abs(ref(k,k))
            do i = k+1, w
                if (abs(ref(i,k)) > best) then
                    best = abs(ref(i,k))
                    r = i
                end if
            end do
            ref_piv(k) = r
            if (r /= k) then
                row = ref(k,:)
                ref(k,:) = ref(r,:)
                ref(r,:) = row
                ref_swaps = ref_swaps+1
            end if
            ref(k+1:,k) = ref(k+1:,k)/ref(k,k)
            do i = k+1, w
                ref(k+1:,i) = ref(k+1:,i)-ref(k+1:,k)*ref(k,i)
            end do
        end do
        call check(all(ipiv == ref_piv), 'pivot rows equal partial pivoting')
        call check(swaps == ref_swaps .and. perturbed == 0, 'exchange count')
        l = (0.0_dp, 0.0_dp)
        u = (0.0_dp, 0.0_dp)
        do k = 1, w
            l(k,k) = (1.0_dp, 0.0_dp)
            l(k+1:,k) = a(k+1:w,k)
            u(1:k,k) = a(1:k,k)
        end do
        ref = original(1:w,1:w)
        do k = 1, w
            row = ref(k,:)
            ref(k,:) = ref(ipiv(k),:)
            ref(ipiv(k),:) = row
        end do
        call check(maxval(abs(matmul(l,u)-ref)) <= 64.0_dp*w*eps*maxval(abs(ref)), &
            'P A = L U reconstruction')
        call check(abs(min_pivot-minval([(abs(u(k,k)), k=1, w)])) == 0.0_dp, &
            'minimum pivot magnitude')
        call check(max_u == maxval(abs(u)), 'maximum |U|')
    end subroutine check_block_factor

    subroutine check_perturbation()
        ! Second pivot exactly zero (perturbed to +tau) or tiny with phase
        ! (3+4i)/5 (perturbed to tau*(3+4i)/5); exactly one each.
        complex(dp) :: a(2,2)
        integer :: ipiv(2), perturbed, swaps
        real(dp) :: min_pivot, max_u, tau
        tau = 1.0e-10_dp
        a = reshape([(1.0_dp,1.0_dp), (1.0_dp,1.0_dp), (2.0_dp,0.0_dp), &
            (2.0_dp,0.0_dp)], [2,2])
        perturbed = 0
        swaps = 0
        min_pivot = huge(1.0_dp)
        max_u = 0.0_dp
        call slu_z_block_factor(a, 2, 2, tau, ipiv, perturbed, swaps, min_pivot, max_u)
        call check(perturbed == 1 .and. swaps == 0, 'one perturbation, no exchange')
        call check(a(2,2) == cmplx(tau, 0.0_dp, dp), 'zero pivot becomes +tau')
        a = reshape([(1.0_dp,0.0_dp), (1.0_dp,0.0_dp), (1.0_dp,0.0_dp), &
            cmplx(1.0_dp+3.0e-17_dp, 4.0e-17_dp, dp)], [2,2])
        a(2,2) = (1.0_dp,0.0_dp)+cmplx(3.0e-12_dp, 4.0e-12_dp, dp)
        perturbed = 0
        call slu_z_block_factor(a, 2, 2, tau, ipiv, perturbed, swaps, min_pivot, max_u)
        call check(perturbed == 1, 'tiny pivot perturbed')
        call check(abs(a(2,2)-tau*(0.6_dp,0.8_dp)) <= 1.0e-4_dp*tau, &
            'perturbation keeps the pivot phase')
        call check(abs(min_pivot-tau) <= 1.0e-12_dp*tau, 'min pivot is tau')
    end subroutine check_perturbation

    subroutine check_triangular(w, rows)
        integer, intent(in) :: w, rows
        complex(dp) :: f(w,w), x(rows,w), b(rows,w), l(w,w), u(w,w)
        integer :: k
        call fill(f)
        call fill(x)
        do k = 1, w
            f(k,k) = f(k,k)+(3.0_dp, 0.0_dp)
        end do
        u = (0.0_dp, 0.0_dp)
        l = (0.0_dp, 0.0_dp)
        do k = 1, w
            u(1:k,k) = f(1:k,k)
            l(k,k) = (1.0_dp, 0.0_dp)
            l(k+1:,k) = f(k+1:,k)
        end do
        b = matmul(x, u)
        call slu_z_trsm_upper(f, w, w, b, rows, rows)
        call check(maxval(abs(b-x)) <= 1.0e3_dp*eps, 'B U^-1')
        b = matmul(x, transpose(l))
        call slu_z_trsm_unit_lower_t(f, w, w, b, rows, rows)
        call check(maxval(abs(b-x)) <= 1.0e3_dp*eps, 'B L^-T (plain transpose)')
    end subroutine check_triangular

    subroutine check_schur(m, n, k)
        integer, intent(in) :: m, n, k
        complex(dp) :: c(m,n), c0(m,n), left(m,k), right(n,k)
        call fill(c)
        call fill(left)
        call fill(right)
        c0 = c
        call slu_z_schur_update(c, m, left, m, right, n, m, n, k)
        call check(maxval(abs(c-(c0-matmul(left, transpose(right))))) <= 64.0_dp*eps, &
            'Schur update C - left right^T')
    end subroutine check_schur

    subroutine check_block_solves(w)
        integer, intent(in) :: w
        complex(dp) :: f(w,w), l(w,w), u(w,w), x(w), t(w)
        integer :: k, ipiv(w)
        call fill(f)
        do k = 1, w
            f(k,k) = f(k,k)+(4.0_dp, 1.0_dp)
        end do
        l = (0.0_dp, 0.0_dp)
        u = (0.0_dp, 0.0_dp)
        do k = 1, w
            l(k,k) = (1.0_dp, 0.0_dp)
            l(k+1:,k) = f(k+1:,k)
            u(1:k,k) = f(1:k,k)
        end do
        x = [(random_z(), k=1, w)]
        t = matmul(l, x)
        call slu_z_block_lower(f, w, w, t)
        call check(maxval(abs(t-x)) <= 1.0e3_dp*eps, 'L^-1 t')
        t = matmul(u, x)
        call slu_z_block_upper(f, w, w, t)
        call check(maxval(abs(t-x)) <= 1.0e3_dp*eps, 'U^-1 t')
        t = matmul(transpose(u), x)
        call slu_z_block_upper_t(f, w, w, t)
        call check(maxval(abs(t-x)) <= 1.0e3_dp*eps, 'U^-T t')
        t = matmul(transpose(l), x)
        call slu_z_block_lower_t(f, w, w, t)
        call check(maxval(abs(t-x)) <= 1.0e3_dp*eps, 'L^-T t')
        do k = 1, w
            ipiv(k) = k+int(uniform()*(w-k+1))
        end do
        t = x
        call slu_z_apply_swaps(t, ipiv, w, .false.)
        call slu_z_apply_swaps(t, ipiv, w, .true.)
        call check(all(t == x), 'exchanges and their transpose are inverse')
    end subroutine check_block_solves

    subroutine check_scatter(n)
        integer, intent(in) :: n
        complex(dp) :: x(2*n), x0(2*n), values(n), v, g, g0
        integer :: indices(n), rowmap(2*n), k
        do k = 1, n
            indices(k) = 2*k-mod(k,2)
            values(k) = random_z()
        end do
        x = [(random_z(), k=1, 2*n)]
        x0 = x
        v = (0.3_dp, -0.7_dp)
        call slu_z_scatter_update(x, indices, values, n, v)
        call check(maxval(abs(x(indices)-(x0(indices)-values*v))) <= 8.0_dp*eps, &
            'scatter update')
        g0 = (1.0_dp, 2.0_dp)
        g = g0
        call slu_z_gather_update(g, indices, values, x, n)
        call check(abs(g-(g0-sum(values*x(indices)))) <= 64.0_dp*eps, 'gather update')
        rowmap = 0
        do k = 1, n, 2
            rowmap(indices(k)) = k
        end do
        x = (0.0_dp, 0.0_dp)
        call slu_z_scatter_mapped_update(x, indices, rowmap, values, n, v)
        do k = 1, n
            if (mod(k,2) == 1) then
                call check(abs(x(k)+values(k)*v) <= 8.0_dp*eps, 'mapped entry')
            else
                call check(x(k) == (0.0_dp, 0.0_dp), 'unmapped entry skipped')
            end if
        end do
    end subroutine check_scatter

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (condition) return
        print '(a)', 'FAIL: '//message
        error stop 1
    end subroutine check
end program test_fortsparse_slu_z_leaves
