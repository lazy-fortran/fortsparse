    ! Complex execution schedules only: every arithmetic operation is one of
    ! the generated Cartesian leaves above. Pivot selection, exchanges, the
    ! perturbation decision and power-of-two normalization are comparisons,
    ! data movement and exact exponent shifts. Magnitudes are the modulus
    ! |z| (the maximum-product matching weight). Dummies are assumed-size as
    ! in the real schedule. Nothing is allocated.

    pure elemental function slu_z_update(a, b, c) result(value)
        complex(dp), intent(in) :: a, b, c
        complex(dp) :: value
        real(dp) :: vr, vi
        call slu_z_update_value(real(a,dp), aimag(a), real(b,dp), aimag(b), &
            real(c,dp), aimag(c), vr, vi)
        value = cmplx(vr, vi, dp)
    end function slu_z_update

    ! a/b for b /= 0. Both operands are normalized by exact powers of two
    ! first, so |b|**2 is never formed at the original scale.
    pure elemental function slu_z_divide(a, b) result(value)
        complex(dp), intent(in) :: a, b
        complex(dp) :: value
        real(dp) :: am, bm, vr, vi
        integer :: ea, eb
        am = max(abs(real(a,dp)), abs(aimag(a)))
        bm = max(abs(real(b,dp)), abs(aimag(b)))
        if (am == 0.0_dp .and. bm > 0.0_dp) then
            value = (0.0_dp, 0.0_dp)
            return
        end if
        ea = exponent(am)
        eb = exponent(bm)
        call slu_z_quotient_value(scale(real(a,dp),-ea), scale(aimag(a),-ea), &
            scale(real(b,dp),-eb), scale(aimag(b),-eb), vr, vi)
        value = cmplx(scale(vr,ea-eb), scale(vi,ea-eb), dp)
    end function slu_z_divide

    ! a*b*c with real scales a and c.
    pure elemental function slu_z_scale(a, b, c) result(value)
        real(dp), intent(in) :: a, c
        complex(dp), intent(in) :: b
        complex(dp) :: value
        real(dp) :: vr, vi
        call slu_z_scale_value(a, real(b,dp), aimag(b), c, vr, vi)
        value = cmplx(vr, vi, dp)
    end function slu_z_scale

    ! Restricted-pivot LU of the w-by-w diagonal block (see the real
    ! slu_block_factor). A pivot with |p| < tau_abs becomes tau_abs*p/|p|,
    ! or +tau_abs for p = 0, and is counted.
    pure subroutine slu_z_block_factor(a, lda, w, tau_abs, ipiv, perturbed, &
            swaps, min_pivot, max_u)
        integer, intent(in) :: lda, w
        complex(dp), intent(inout) :: a(lda,*)
        real(dp), intent(in) :: tau_abs
        integer, intent(out) :: ipiv(*)
        integer, intent(inout) :: perturbed, swaps
        real(dp), intent(inout) :: min_pivot, max_u
        integer :: i, j, k, r
        real(dp) :: best
        complex(dp) :: saved
        do k = 1, w
            r = k
            best = abs(a(k,k))
            do i = k+1, w
                if (abs(a(i,k)) > best) then
                    best = abs(a(i,k))
                    r = i
                end if
            end do
            ipiv(k) = r
            if (r /= k) then
                swaps = swaps+1
                do j = 1, w
                    saved = a(k,j)
                    a(k,j) = a(r,j)
                    a(r,j) = saved
                end do
            end if
            if (best < tau_abs) then
                if (best > 0.0_dp) then
                    a(k,k) = slu_z_scale(tau_abs, slu_z_divide(a(k,k), &
                        cmplx(best, 0.0_dp, dp)), 1.0_dp)
                else
                    a(k,k) = cmplx(tau_abs, 0.0_dp, dp)
                end if
                perturbed = perturbed+1
            end if
            min_pivot = min(min_pivot, abs(a(k,k)))
            do i = k+1, w
                a(i,k) = slu_z_divide(a(i,k), a(k,k))
            end do
            do j = k+1, w
                !$omp simd
                do i = k+1, w
                    a(i,j) = slu_z_update(a(i,j), a(i,k), a(k,j))
                end do
            end do
        end do
        do j = 1, w
            do i = 1, j
                max_u = max(max_u, abs(a(i,j)))
            end do
        end do
    end subroutine slu_z_block_factor

    ! B := B U^(-1); U upper non-unit in u(1:w,1:w).
    pure subroutine slu_z_trsm_upper(u, ldu, w, b, ldb, rows)
        integer, intent(in) :: ldu, w, ldb, rows
        complex(dp), intent(in) :: u(ldu,*)
        complex(dp), intent(inout) :: b(ldb,*)
        integer :: i, j, k
        do j = 1, w
            do k = 1, j-1
                !$omp simd
                do i = 1, rows
                    b(i,j) = slu_z_update(b(i,j), b(i,k), u(k,j))
                end do
            end do
            do i = 1, rows
                b(i,j) = slu_z_divide(b(i,j), u(j,j))
            end do
        end do
    end subroutine slu_z_trsm_upper

    ! B := B L^(-T) (plain transpose, no conjugation); L unit lower.
    pure subroutine slu_z_trsm_unit_lower_t(l, ldl, w, b, ldb, rows)
        integer, intent(in) :: ldl, w, ldb, rows
        complex(dp), intent(in) :: l(ldl,*)
        complex(dp), intent(inout) :: b(ldb,*)
        integer :: i, j, k
        do j = 1, w
            do k = 1, j-1
                !$omp simd
                do i = 1, rows
                    b(i,j) = slu_z_update(b(i,j), b(i,k), l(j,k))
                end do
            end do
        end do
    end subroutine slu_z_trsm_unit_lower_t

    ! C(1:m,1:n) -= left(1:m,1:k) right(1:n,1:k)^T (plain transpose).
    pure subroutine slu_z_schur_update(c, ldc, left, ldl, right, ldr, m, n, k)
        integer, intent(in) :: ldc, ldl, ldr, m, n, k
        complex(dp), intent(inout) :: c(ldc,*)
        complex(dp), intent(in) :: left(ldl,*), right(ldr,*)
        integer :: i, j, p
        do j = 1, n
            do p = 1, k
                !$omp simd
                do i = 1, m
                    c(i,j) = slu_z_update(c(i,j), left(i,p), right(j,p))
                end do
            end do
        end do
    end subroutine slu_z_schur_update

    ! Block-local exchanges: k = 1..w for P, k = w..1 for P^T.
    pure subroutine slu_z_apply_swaps(t, ipiv, w, transpose)
        integer, intent(in) :: w, ipiv(*)
        complex(dp), intent(inout) :: t(*)
        logical, intent(in) :: transpose
        integer :: k, first, last, step
        complex(dp) :: saved
        first = 1
        last = w
        step = 1
        if (transpose) then
            first = w
            last = 1
            step = -1
        end if
        do k = first, last, step
            if (ipiv(k) == k) cycle
            saved = t(k)
            t(k) = t(ipiv(k))
            t(ipiv(k)) = saved
        end do
    end subroutine slu_z_apply_swaps

    ! t := L^(-1) t, unit lower triangle of the packed block.
    pure subroutine slu_z_block_lower(a, lda, w, t)
        integer, intent(in) :: lda, w
        complex(dp), intent(in) :: a(lda,*)
        complex(dp), intent(inout) :: t(*)
        integer :: i, j
        do j = 1, w
            !$omp simd
            do i = j+1, w
                t(i) = slu_z_update(t(i), a(i,j), t(j))
            end do
        end do
    end subroutine slu_z_block_lower

    ! t := U^(-1) t.
    pure subroutine slu_z_block_upper(a, lda, w, t)
        integer, intent(in) :: lda, w
        complex(dp), intent(in) :: a(lda,*)
        complex(dp), intent(inout) :: t(*)
        integer :: i, j
        do j = w, 1, -1
            t(j) = slu_z_divide(t(j), a(j,j))
            !$omp simd
            do i = 1, j-1
                t(i) = slu_z_update(t(i), a(i,j), t(j))
            end do
        end do
    end subroutine slu_z_block_upper

    ! t := U^(-T) t (plain transpose).
    pure subroutine slu_z_block_upper_t(a, lda, w, t)
        integer, intent(in) :: lda, w
        complex(dp), intent(in) :: a(lda,*)
        complex(dp), intent(inout) :: t(*)
        integer :: i, j
        complex(dp) :: value
        do j = 1, w
            value = t(j)
            do i = 1, j-1
                value = slu_z_update(value, a(i,j), t(i))
            end do
            t(j) = slu_z_divide(value, a(j,j))
        end do
    end subroutine slu_z_block_upper_t

    ! t := L^(-T) t, unit lower triangle.
    pure subroutine slu_z_block_lower_t(a, lda, w, t)
        integer, intent(in) :: lda, w
        complex(dp), intent(in) :: a(lda,*)
        complex(dp), intent(inout) :: t(*)
        integer :: i, j
        complex(dp) :: value
        do j = w, 1, -1
            value = t(j)
            do i = j+1, w
                value = slu_z_update(value, a(i,j), t(i))
            end do
            t(j) = value
        end do
    end subroutine slu_z_block_lower_t

    ! x(rowmap(rows(i))) -= values(i)*v for mapped (nonzero) positions.
    pure subroutine slu_z_scatter_mapped_update(x, rows, rowmap, values, n, v)
        integer, intent(in) :: n, rows(*), rowmap(*)
        complex(dp), intent(inout) :: x(*)
        complex(dp), intent(in) :: values(*), v
        integer :: i, position
        do i = 1, n
            position = rowmap(rows(i))
            if (position > 0) x(position) = slu_z_update(x(position), values(i), v)
        end do
    end subroutine slu_z_scatter_mapped_update

    ! x(indices(i)) -= values(i)*v.
    pure subroutine slu_z_scatter_update(x, indices, values, n, v)
        integer, intent(in) :: n, indices(*)
        complex(dp), intent(inout) :: x(*)
        complex(dp), intent(in) :: values(*), v
        integer :: i
        do i = 1, n
            x(indices(i)) = slu_z_update(x(indices(i)), values(i), v)
        end do
    end subroutine slu_z_scatter_update

    ! value -= sum_i values(i)*x(indices(i)).
    pure subroutine slu_z_gather_update(value, indices, values, x, n)
        integer, intent(in) :: n, indices(*)
        complex(dp), intent(inout) :: value
        complex(dp), intent(in) :: values(*), x(*)
        integer :: i
        do i = 1, n
            value = slu_z_update(value, values(i), x(indices(i)))
        end do
    end subroutine slu_z_gather_update
