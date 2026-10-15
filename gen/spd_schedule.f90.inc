    ! Execution schedules only: all numerical leaves come from FortSym above.
    ! Factor lower triangle in place. Upper triangle is never read or changed.
    subroutine spd_panel_factor(a, n, lda, info)
        integer, intent(in) :: n, lda
        real(dp), intent(inout) :: a(lda, n)
        integer, intent(out) :: info
        integer :: i, j, k
        info = 0
        do j = 1, n
            do k = 1, j-1
                a(j,j) = spd_update(a(j,j), a(j,k), a(j,k))
            end do
            if (.not. a(j,j) > 0.0_dp) then
                info = j
                return
            end if
            a(j,j) = spd_root(a(j,j))
            do k = 1, j-1
                !$omp simd
                do i = j+1, n
                    a(i,j) = spd_update(a(i,j), a(i,k), a(j,k))
                end do
            end do
            !$omp simd
            do i = j+1, n
                a(i,j) = spd_divide(a(i,j), a(j,j))
            end do
        end do
    end subroutine spd_panel_factor

    ! Right-side triangular solve B := B L^(-T), B is m by n.
    subroutine spd_panel_solve(l, n, ldl, b, m, ldb)
        integer, intent(in) :: n, ldl, m, ldb
        real(dp), intent(in) :: l(ldl, n)
        real(dp), intent(inout) :: b(ldb, n)
        integer :: i, j, k
        do j = 1, n
            do k = 1, j-1
                !$omp simd
                do i = 1, m
                    b(i,j) = spd_update(b(i,j), b(i,k), l(j,k))
                end do
            end do
            !$omp simd
            do i = 1, m
                b(i,j) = spd_divide(b(i,j), l(j,j))
            end do
        end do
    end subroutine spd_panel_solve

    ! Full rectangular Schur update C := C - left right^T.
    ! The caller owns parallel tasks; SIMD avoids nested thread teams here.
    subroutine spd_schur_update(c, ldc, left, ldl, right, ldr, m, n, k)
        integer, intent(in) :: ldc, ldl, ldr, m, n, k
        real(dp), intent(inout) :: c(ldc, n)
        real(dp), intent(in) :: left(ldl, k), right(ldr, k)
        integer :: i, j, p
        do j = 1, n
            do p = 1, k
                !$omp simd
                do i = 1, m
                    c(i,j) = spd_update(c(i,j), left(i,p), right(j,p))
                end do
            end do
        end do
    end subroutine spd_schur_update
