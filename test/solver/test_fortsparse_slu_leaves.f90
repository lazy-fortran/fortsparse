program test_fortsparse_slu_leaves
    ! Each compiled generated supernodal-LU kernel against an independent
    ! dense reference: textbook partial pivoting restricted to the block,
    ! factor reconstruction, matmul residuals and exact counts.
    use fortsparse_kinds, only: dp, i8
    use fortsparse_slu_leaves, only: slu_update, slu_divide, slu_scale, &
        slu_block_factor, slu_trsm_upper, slu_trsm_unit_lower_t, slu_schur_update, &
        slu_apply_swaps, slu_block_lower, slu_block_upper, slu_block_upper_t, &
        slu_block_lower_t
    implicit none
    integer :: seed, w

    seed = 4242
    call check(slu_update(5.0_dp, 2.0_dp, 3.0_dp) == -1.0_dp, 'update leaf')
    call check(slu_divide(3.0_dp, 4.0_dp) == 0.75_dp, 'divide leaf')
    call check(slu_scale(2.0_dp, 3.0_dp, 0.5_dp) == 3.0_dp, 'scale leaf')
    do w = 1, 17
        call check_block_factor(w)
        call check_triangular(w)
    end do
    call check_perturbation()
    print '(a)', 'supernodal LU generated leaves passed'
contains
    real(dp) function uniform()
        seed = modulo(1103515245_i8*seed+12345_i8, 2147483648_i8)
        uniform = 2.0_dp*real(seed, dp)/2147483648.0_dp-1.0_dp
    end function uniform

    subroutine reference_lu(a, ipiv)
        ! Textbook right-looking partial pivoting (first maximum), no
        ! perturbation; whole rows are exchanged.
        real(dp), intent(inout) :: a(:,:)
        integer, intent(out) :: ipiv(:)
        integer :: n, k, i, r
        real(dp) :: row(size(a,2))
        n = size(a,1)
        do k = 1, n
            r = k
            do i = k+1, n
                if (abs(a(i,k)) > abs(a(r,k))) r = i
            end do
            ipiv(k) = r
            row = a(k,:)
            a(k,:) = a(r,:)
            a(r,:) = row
            a(k+1:,k) = a(k+1:,k)/a(k,k)
            do i = k+1, n
                a(k+1:,i) = a(k+1:,i)-a(k+1:,k)*a(k,i)
            end do
        end do
    end subroutine reference_lu

    subroutine check_block_factor(w)
        integer, intent(in) :: w
        integer, parameter :: extra = 3
        real(dp) :: a(w+extra,w), original(w+extra,w), ref(w,w), l(w,w), u(w,w)
        real(dp) :: permuted(w,w), min_pivot, max_u, row(w)
        integer :: ipiv(w), ref_ipiv(w), perturbed, swaps, i, j, k
        do j = 1, w
            do i = 1, w+extra
                a(i,j) = uniform()
            end do
        end do
        original = a
        perturbed = 0
        swaps = 0
        min_pivot = huge(1.0_dp)
        max_u = 0.0_dp
        call slu_block_factor(a, w+extra, w, 1.0e-300_dp, ipiv, perturbed, swaps, &
            min_pivot, max_u)
        call check(all(a(w+1:,:) == original(w+1:,:)), 'rows below the block untouched')
        ref = original(1:w,:)
        call reference_lu(ref, ref_ipiv)
        call check(all(ipiv == ref_ipiv), 'pivot rows equal textbook partial pivoting')
        call check(maxval(abs(a(1:w,:)-ref)) <= 1.0e-12_dp*max(1.0_dp, maxval(abs(ref))), &
            'packed factors equal the reference')
        call check(perturbed == 0, 'no perturbation above the threshold')
        call check(swaps == count(ipiv /= [(k, k=1, w)]), 'exchange count')
        l = 0.0_dp
        u = 0.0_dp
        do j = 1, w
            l(j,j) = 1.0_dp
            l(j+1:,j) = a(j+1:w,j)
            u(1:j,j) = a(1:j,j)
        end do
        permuted = original(1:w,:)
        do k = 1, w
            row = permuted(k,:)
            permuted(k,:) = permuted(ipiv(k),:)
            permuted(ipiv(k),:) = row
        end do
        call check(maxval(abs(matmul(l,u)-permuted)) <= 1.0e-12_dp*w, 'P A = L U')
        call check(abs(min_pivot-minval([(abs(u(k,k)), k=1, w)])) == 0.0_dp, 'min pivot')
        call check(max_u == maxval(abs(u)), 'max |U|')
    end subroutine check_block_factor

    subroutine check_triangular(w)
        integer, intent(in) :: w
        integer, parameter :: rows = 5
        real(dp) :: lu(w+2,w), b(rows+1,w), x(rows+1,w), l(w,w), u(w,w)
        real(dp) :: left(rows,4), right(w,4), c(rows,w), t(w), s(w)
        integer :: ipiv(w), i, j, k
        do j = 1, w
            do i = 1, w+2
                lu(i,j) = uniform()
            end do
            lu(j,j) = 2.0_dp+abs(lu(j,j))
            do i = 1, rows+1
                b(i,j) = uniform()
            end do
        end do
        l = 0.0_dp
        u = 0.0_dp
        do j = 1, w
            l(j,j) = 1.0_dp
            l(j+1:,j) = lu(j+1:w,j)
            u(1:j,j) = lu(1:j,j)
        end do
        x = b
        call slu_trsm_upper(lu, w+2, w, x, rows+1, rows)
        call check(maxval(abs(matmul(x(1:rows,:),u)-b(1:rows,:))) <= 1.0e-13_dp*w, 'X U = B')
        call check(all(x(rows+1,:) == b(rows+1,:)), 'trsm leaves rows beyond the block')
        x = b
        call slu_trsm_unit_lower_t(lu, w+2, w, x, rows+1, rows)
        call check(maxval(abs(matmul(x(1:rows,:),transpose(l))-b(1:rows,:))) &
            <= 1.0e-13_dp*w, 'X L^T = B')
        do k = 1, 4
            do i = 1, rows
                left(i,k) = uniform()
            end do
            do j = 1, w
                right(j,k) = uniform()
            end do
        end do
        c = b(1:rows,:)
        call slu_schur_update(c, rows, left, rows, right, w, rows, w, 4)
        call check(maxval(abs(c-(b(1:rows,:)-matmul(left,transpose(right))))) &
            <= 1.0e-14_dp*4, 'Schur update')
        s = b(1,:)
        t = s
        call slu_block_lower(lu, w+2, w, t)
        call check(maxval(abs(matmul(l,t)-s)) <= 1.0e-13_dp*w, 'L t = s')
        t = s
        call slu_block_upper(lu, w+2, w, t)
        call check(maxval(abs(matmul(u,t)-s)) <= 1.0e-13_dp*w, 'U t = s')
        t = s
        call slu_block_upper_t(lu, w+2, w, t)
        call check(maxval(abs(matmul(transpose(u),t)-s)) <= 1.0e-13_dp*w, 'U^T t = s')
        t = s
        call slu_block_lower_t(lu, w+2, w, t)
        call check(maxval(abs(matmul(transpose(l),t)-s)) <= 1.0e-13_dp*w, 'L^T t = s')
        do k = 1, w
            ipiv(k) = k+int((w-k+1)*(0.5_dp+0.5_dp*uniform()))
            ipiv(k) = min(w, max(k, ipiv(k)))
        end do
        t = s
        do k = 1, w
            x(1,1) = t(k)
            t(k) = t(ipiv(k))
            t(ipiv(k)) = x(1,1)
        end do
        x(1,:) = s
        call slu_apply_swaps(x(1,:), ipiv, w, .false.)
        call check(all(x(1,:) == t), 'forward exchanges')
        call slu_apply_swaps(t, ipiv, w, .true.)
        call check(all(t == s), 'transposed exchanges invert the forward ones')
    end subroutine check_triangular

    subroutine check_perturbation()
        ! Exact zero second pivot is perturbed to +tau; a tiny negative one
        ! to -tau; counts and min pivot are exact.
        real(dp) :: a(2,2), min_pivot, max_u
        integer :: ipiv(2), perturbed, swaps
        a = reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [2,2])
        perturbed = 0
        swaps = 0
        min_pivot = huge(1.0_dp)
        max_u = 0.0_dp
        call slu_block_factor(a, 2, 2, 1.0e-10_dp, ipiv, perturbed, swaps, min_pivot, max_u)
        call check(perturbed == 1 .and. a(2,2) == 1.0e-10_dp, 'zero pivot becomes +tau')
        call check(min_pivot == 1.0e-10_dp .and. swaps == 0, 'zero pivot statistics')
        a = reshape([1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp-2.0_dp**(-40)], [2,2])
        call slu_block_factor(a, 2, 2, 1.0e-10_dp, ipiv, perturbed, swaps, min_pivot, max_u)
        call check(perturbed == 2 .and. a(2,2) == -1.0e-10_dp, 'negative pivot becomes -tau')
        a = reshape([1.0e-12_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2,2])
        call slu_block_factor(a, 2, 2, 1.0e-10_dp, ipiv, perturbed, swaps, min_pivot, max_u)
        call check(perturbed == 3 .and. a(1,1) == 1.0e-10_dp, 'small first pivot perturbed')
    end subroutine check_perturbation

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (condition) return
        print '(a)', 'FAIL: '//message
        error stop 1
    end subroutine check
end program test_fortsparse_slu_leaves
