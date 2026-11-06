module fortsparse_native_sytrf
    ! Dense symmetric indefinite block factor Q^T A Q = L D L^T (lower) in the
    ! DSYTRF_RK storage: unit L in the strict lower triangle of A with zero
    ! entries inside 2x2 blocks, diag(D) on the diagonal, subdiag(D) in E and
    ! interchanges k <-> |ipiv(k)| applied in order k = 1..n (ipiv(k) < 0 marks
    ! a 2x2 block at k, k+1). LAPACK < 3.7 (e.g. Apple Accelerate) has no
    ! DSYTRF_RK; the build then factors with DSYTRF and converts its storage.
    use fortsparse_kinds, only: dp
    implicit none
    private
    public :: sytrf_rk_lower, sytrf_to_rk_lower

contains

    subroutine sytrf_rk_lower(n, a, lda, e, ipiv, work, lwork, info)
        integer, intent(in) :: n, lda, lwork
        real(dp), intent(inout) :: a(lda,*)
        real(dp), intent(out) :: e(*), work(*)
        integer, intent(out) :: ipiv(*), info
        interface
#ifdef FORTSPARSE_LAPACK_NO_DSYTRF_RK
            subroutine dsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
                import dp
                character, intent(in) :: uplo
                integer, intent(in) :: n, lda, lwork
                real(dp), intent(inout) :: a(lda,*)
                real(dp), intent(out) :: work(*)
                integer, intent(out) :: ipiv(*), info
            end subroutine dsytrf
#else
            subroutine dsytrf_rk(uplo, n, a, lda, e, ipiv, work, lwork, info)
                import dp
                character, intent(in) :: uplo
                integer, intent(in) :: n, lda, lwork
                real(dp), intent(inout) :: a(lda,*)
                real(dp), intent(out) :: e(*), work(*)
                integer, intent(out) :: ipiv(*), info
            end subroutine dsytrf_rk
#endif
        end interface
#ifdef FORTSPARSE_LAPACK_NO_DSYTRF_RK
        call dsytrf('L', n, a, lda, ipiv, work, lwork, info)
        if (info >= 0) call sytrf_to_rk_lower(n, a, lda, e, ipiv)
#else
        call dsytrf_rk('L', n, a, lda, e, ipiv, work, lwork, info)
#endif
    end subroutine sytrf_rk_lower

    subroutine sytrf_to_rk_lower(n, a, lda, e, ipiv)
        ! DSYTRF stores L = P(1) L(1) P(2) L(2) ... with D(k+1,k) in A and one
        ! interchange k+1 <-> -ipiv(k) per 2x2 block; LAPACK 3.7 DSYCONVF
        ! ('L', 'C') defines this conversion to the DSYTRF_RK storage.
        integer, intent(in) :: n, lda
        real(dp), intent(inout) :: a(lda,*)
        real(dp), intent(out) :: e(*)
        integer, intent(inout) :: ipiv(*)
        integer :: i, ip, j
        real(dp) :: saved
        if (n < 1) return
        e(n) = 0.0_dp
        i = 1
        do while (i <= n)
            if (ipiv(i) < 0 .and. i < n) then
                e(i) = a(i+1,i)
                e(i+1) = 0.0_dp
                a(i+1,i) = 0.0_dp
                i = i+2
            else
                e(i) = 0.0_dp
                i = i+1
            end if
        end do
        i = 1
        do while (i <= n)
            if (ipiv(i) > 0) then
                ip = ipiv(i)
                if (ip /= i) then
                    do j = 1, i-1
                        saved = a(i,j)
                        a(i,j) = a(ip,j)
                        a(ip,j) = saved
                    end do
                end if
                i = i+1
            else
                ip = -ipiv(i)
                if (ip /= i+1) then
                    do j = 1, i-1
                        saved = a(i+1,j)
                        a(i+1,j) = a(ip,j)
                        a(ip,j) = saved
                    end do
                end if
                ipiv(i) = -i
                i = i+2
            end if
        end do
    end subroutine sytrf_to_rk_lower
end module fortsparse_native_sytrf
