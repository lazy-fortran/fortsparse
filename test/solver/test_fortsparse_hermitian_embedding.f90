program test_fortsparse_hermitian_embedding
    ! Complex Hermitian block-tridiagonal A = L D L^H (unit lower block
    ! bidiagonal L, real diagonal D with prescribed signs): Sylvester's law
    ! gives the exact inertia. The supernodal real embedding must return it,
    ! agree with the scalar complex LDL^H and solve A x = b.
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csc, only: csc_z_t, csc_from_triplet, csc_matvec
    use fortsparse_status, only: fortsparse_status_t, status_ok
    use fortsparse, only: sparse_solver_t, sparse_factor, sparse_solve, &
        sparse_inertia, sparse_destroy, FORTSPARSE_BACKEND_NATIVE, &
        FORTSPARSE_SYM_HERMITIAN, FORTSPARSE_SYM_SPD
    implicit none
    integer, parameter :: nb = 12, w = 9, n = nb*w
    complex(dp), allocatable :: dense(:, :), lower(:, :), b(:), x(:), vals(:)
    real(dp), allocatable :: d(:)
    integer, allocatable :: rows(:), cols(:)
    type(csc_z_t) :: a
    type(fortsparse_status_t) :: status
    type(sparse_solver_t) :: embedded, scalar
    integer(i8) :: neg(2), zer(2), pos(2)
    integer :: i, j, k, rel, expected
    allocate(lower(n, n), d(n), dense(n, n))
    lower = 0
    do i = 1, n
        lower(i, i) = 1
        d(i) = (1 + mod(7*i, 5))*merge(-1.0_dp, 1.0_dp, mod(i, 7) == 3)
    end do
    do j = 1, n
        do i = j + 1, min(n, ((j - 1)/w + 2)*w)
            lower(i, j) = cmplx(sin(0.3_dp*i + 0.7_dp*j), cos(1.1_dp*i - 0.2_dp*j), dp)/3
        end do
    end do
    expected = count(d < 0)
    dense = matmul(lower, matmul(diag(d), conjg(transpose(lower))))
    dense = (dense + conjg(transpose(dense)))/2
    allocate(rows(n*n), cols(n*n), vals(n*n))
    k = 0
    do j = 1, n
        do i = 1, n
            if (abs(dense(i, j)) <= 0) cycle
            k = k + 1; rows(k) = i; cols(k) = j; vals(k) = dense(i, j)
            if (i == j) vals(k) = cmplx(real(dense(i, j), dp), 0, dp)
        end do
    end do
    call csc_from_triplet(n, n, rows(:k), cols(:k), vals(:k), a, status)
    if (.not. status_ok(status)) error stop 'matrix'
    embedded%backend_id = FORTSPARSE_BACKEND_NATIVE
    embedded%options%symmetry = FORTSPARSE_SYM_HERMITIAN
    embedded%options%ldlt_supernodal = .true.
    scalar%backend_id = FORTSPARSE_BACKEND_NATIVE
    scalar%options%symmetry = FORTSPARSE_SYM_HERMITIAN
    call sparse_factor(embedded, a, status)
    if (.not. status_ok(status)) then
        print '(a,i0,1x,a)', 'status ', status%code, trim(status%msg)
        error stop 'embedded factor'
    end if
    call sparse_factor(scalar, a, status)
    if (.not. status_ok(status)) error stop 'scalar factor'
    call sparse_inertia(embedded, neg(1), zer(1), pos(1), rel, status)
    if (.not. status_ok(status)) error stop 'embedded inertia'
    call sparse_inertia(scalar, neg(2), zer(2), pos(2), rel, status)
    if (.not. status_ok(status)) error stop 'scalar inertia'
    if (neg(1) /= expected .or. zer(1) /= 0 .or. pos(1) /= n - expected) &
        error stop 'embedded inertia differs from Sylvester'
    if (any([neg(2), zer(2), pos(2)] /= [neg(1), zer(1), pos(1)])) &
        error stop 'embedded and scalar inertia differ'
    allocate(b(n), x(n))
    do i = 1, n
        b(i) = cmplx(cos(0.5_dp*i), sin(0.9_dp*i), dp)
    end do
    call sparse_solve(embedded, b, x, status)
    if (.not. status_ok(status)) error stop 'embedded solve'
    ! Normwise backward error, independent of the conditioning of A.
    if (maxval(abs(csc_matvec(a, x) - b)) > 1e-12_dp*(maxval(abs(a%val)) &
        *sum(abs(x)) + maxval(abs(b)))) error stop 'embedded solve backward error'
    call sparse_destroy(embedded); call sparse_destroy(scalar)
    ! Positive definite: |D| gives the supernodal Cholesky of the embedding.
    dense = matmul(lower, matmul(diag(abs(d)), conjg(transpose(lower))))
    dense = (dense + conjg(transpose(dense)))/2
    k = 0
    do j = 1, n
        do i = 1, n
            if (abs(dense(i, j)) <= 0) cycle
            k = k + 1; rows(k) = i; cols(k) = j; vals(k) = dense(i, j)
            if (i == j) vals(k) = cmplx(real(dense(i, j), dp), 0, dp)
        end do
    end do
    call csc_from_triplet(n, n, rows(:k), cols(:k), vals(:k), a, status)
    embedded%backend_id = FORTSPARSE_BACKEND_NATIVE
    embedded%options%symmetry = FORTSPARSE_SYM_SPD
    embedded%options%ldlt_supernodal = .true.
    call sparse_factor(embedded, a, status)
    if (.not. status_ok(status)) error stop 'embedded positive definite factor'
    call sparse_solve(embedded, b, x, status)
    if (.not. status_ok(status)) error stop 'embedded positive definite solve'
    if (maxval(abs(csc_matvec(a, x) - b)) > 1e-12_dp*(maxval(abs(a%val)) &
        *sum(abs(x)) + maxval(abs(b)))) error stop 'positive definite backward error'
    call sparse_destroy(embedded)
    print '(a,i0,a)', 'PASS Hermitian embedding: ', expected, ' negative pivots'
contains
    pure function diag(v) result(m)
        real(dp), intent(in) :: v(:)
        complex(dp) :: m(size(v), size(v))
        integer :: q
        m = 0
        do q = 1, size(v)
            m(q, q) = v(q)
        end do
    end function diag
end program test_fortsparse_hermitian_embedding
