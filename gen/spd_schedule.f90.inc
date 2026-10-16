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
        interface
            subroutine dtrsm(side, uplo, trans, diag, rows, cols, alpha, &
                    factor, ldf, values, ldv)
                import dp
                character, intent(in) :: side, uplo, trans, diag
                integer, intent(in) :: rows, cols, ldf, ldv
                real(dp), intent(in) :: alpha, factor(ldf,*)
                real(dp), intent(inout) :: values(ldv,*)
            end subroutine dtrsm
        end interface
        if (n >= 96 .and. m >= 32) then
            call dtrsm('R', 'L', 'T', 'N', m, n, 1.0_dp, l, ldl, b, ldb)
            return
        end if
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
    subroutine spd_schur_update(c, ldc, left, ldl, right, ldr, m, n, k, first)
        integer, intent(in) :: ldc, ldl, ldr, m, n, k
        real(dp), intent(inout) :: c(ldc, n)
        real(dp), intent(in) :: left(ldl, k), right(ldr, k)
        integer, optional, intent(in) :: first
        integer :: i, j, p, offset
        offset = 0
        if (present(first)) offset = first-1
        do j = 1, n
            do p = 1, k
                !$omp simd
                do i = 1, m
                    c(i,j) = spd_update(c(i,j), left(i+offset,p), right(j+offset,p))
                end do
            end do
        end do
    end subroutine spd_schur_update

    ! Unique, valid indices are owned and checked by sparse symbolic analysis.
    ! Keep scalar leaves inside one compiled loop rather than cross-module
    ! calls per sparse entry. These kernels allocate no scratch.
    pure subroutine spd_scatter_update(x, indices, values, n, v)
        integer, intent(in) :: n, indices(*)
        real(dp), intent(inout) :: x(*)
        real(dp), intent(in) :: values(*), v
        integer :: i
        !$omp simd
        do i = 1, n
            x(indices(i)) = spd_update(x(indices(i)), values(i), v)
        end do
    end subroutine spd_scatter_update

    pure subroutine spd_gather_update(value, indices, values, x, n)
        integer, intent(in) :: n, indices(*)
        real(dp), intent(inout) :: value
        real(dp), intent(in) :: values(*), x(*)
        integer :: i
        !$omp simd reduction(+:value)
        do i = 1, n
            value = spd_update(value, values(i), x(indices(i)))
        end do
    end subroutine spd_gather_update

    ! The row map encodes zero for structurally absent destination entries.
    ! Source rows are distinct; mapped nonzero destination indices are unique.
    pure subroutine spd_scatter_mapped_update(x, rows, rowmap, values, n, v)
        integer, intent(in) :: n, rows(*), rowmap(*)
        real(dp), intent(inout) :: x(*)
        real(dp), intent(in) :: values(*), v
        integer :: i, position
        !$omp simd private(position)
        do i = 1, n
            position = rowmap(rows(i))
            if (position > 0) &
                x(position) = spd_update(x(position), values(i), v)
        end do
    end subroutine spd_scatter_mapped_update
