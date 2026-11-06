program test_fortsparse_sytrf_rk
    ! Dense indefinite block factor in DSYTRF_RK storage. Oracle: the
    ! reconstruction Q L D L^T Q^T of the input matrix. The DSYTRF storage
    ! conversion is checked directly, so it is covered on every LAPACK, and
    ! the factor the build selected (DSYTRF_RK or DSYTRF + conversion) too.
    ! Zero diagonals force 2x2 pivots; graded entries force interchanges.
    use fortsparse_kinds, only: dp
    use fortsparse_native_sytrf, only: sytrf_rk_lower, sytrf_to_rk_lower
    implicit none
    interface
        subroutine dsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
            import dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda,*)
            real(dp), intent(out) :: work(*)
            integer, intent(out) :: ipiv(*), info
        end subroutine dsytrf
    end interface
    integer, parameter :: n = 12, lda = n+3, lwork = n*64
    real(dp) :: a(n,n), f(lda,n), e(n), work(lwork)
    integer :: ipiv(n), info, failures, i, j, trial
    failures = 0
    do trial = 1, 3
        call sample(trial, a)
        f = 0.0_dp
        f(1:n,1:n) = a
        call dsytrf('L', n, f, lda, ipiv, work, lwork, info)
        call expect(info == 0, 'DSYTRF factors the sample')
        call expect(any(ipiv < 0), 'sample needs 2x2 pivots')
        call expect(any([(abs(ipiv(i)) > i+1, i = 1, n)]), 'sample needs interchanges')
        call sytrf_to_rk_lower(n, f, lda, e, ipiv)
        call check(a, f, e, ipiv, 'converted DSYTRF storage')
        f = 0.0_dp
        f(1:n,1:n) = a
        call sytrf_rk_lower(n, f, lda, e, ipiv, work, lwork, info)
        call expect(info == 0, 'selected factor succeeds')
        call check(a, f, e, ipiv, 'selected factor storage')
    end do
    if (failures /= 0) then
        write(*,'(a,i0,a)') 'FAIL: ', failures, ' check(s)'
        error stop 1
    end if
    write(*,'(a)') 'PASS'

contains

    subroutine sample(seed, m)
        integer, intent(in) :: seed
        real(dp), intent(out) :: m(n,n)
        integer :: p, q
        do q = 1, n
            do p = q, n
                m(p,q) = sin(real(seed*37+p*n+q, dp))*real(1+mod(p*q+seed, 7), dp)
                m(q,p) = m(p,q)
            end do
            m(q,q) = 0.0_dp
        end do
        if (seed == 3) then
            do q = 1, n, 3
                m(q,q) = 10.0_dp*real(q, dp)
            end do
        end if
    end subroutine sample

    subroutine check(m, factor, offdiag, pivots, what)
        ! A = Q L D L^T Q^T with Q = P(1) ... P(n), P(k) swapping k, |ipiv(k)|.
        real(dp), intent(in) :: m(n,n), factor(lda,n), offdiag(n)
        integer, intent(in) :: pivots(n)
        character(*), intent(in) :: what
        real(dp) :: l(n,n), d(n,n), r(n,n), row(n)
        integer :: k, p
        l = 0.0_dp
        d = 0.0_dp
        do k = 1, n
            l(k,k) = 1.0_dp
            l(k+1:n,k) = factor(k+1:n,k)
            d(k,k) = factor(k,k)
        end do
        k = 1
        do while (k <= n)
            if (pivots(k) < 0 .and. k < n) then
                if (factor(k+1,k) /= 0.0_dp) call expect(.false., what//': 2x2 L block')
                l(k+1,k) = 0.0_dp
                d(k+1,k) = offdiag(k)
                d(k,k+1) = offdiag(k)
                k = k+2
            else
                k = k+1
            end if
        end do
        r = matmul(l, matmul(d, transpose(l)))
        do k = n, 1, -1
            p = abs(pivots(k))
            if (p == k) cycle
            row = r(k,:)
            r(k,:) = r(p,:)
            r(p,:) = row
            row = r(:,k)
            r(:,k) = r(:,p)
            r(:,p) = row
        end do
        call expect(maxval(abs(r-m)) <= 1.0e-12_dp*maxval(abs(m)), &
            what//': Q L D L^T Q^T reproduces A')
    end subroutine check

    subroutine expect(condition, what)
        logical, intent(in) :: condition
        character(*), intent(in) :: what
        if (condition) return
        failures = failures+1
        write(*,'(a)') 'FAIL: '//what
    end subroutine expect
end program test_fortsparse_sytrf_rk
