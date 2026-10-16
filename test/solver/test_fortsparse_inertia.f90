program test_fortsparse_inertia
    ! Sylvester inertia of native LDL^T/LDL^H/Cholesky factors against
    ! independent dense LAPACK eigenvalue counts (dsyev, dsygv, zheev).
    use fortsparse
    implicit none
    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character, intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine dsyev
        subroutine dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info)
            import :: dp
            integer, intent(in) :: itype, n, lda, ldb, lwork
            character, intent(in) :: jobz, uplo
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            real(dp), intent(out) :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine dsygv
        subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
            import :: dp
            character, intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            complex(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*), rwork(*)
            complex(dp), intent(out) :: work(*)
            integer, intent(out) :: info
        end subroutine zheev
    end interface
    integer :: ordering
    integer, parameter :: orderings(2) = [FORTSPARSE_ORDER_DEFAULT, &
                                          FORTSPARSE_ORDER_AMD]

    call random_seed_fixed()
    do ordering = 1, 2
        if (.not. ordering_available(orderings(ordering))) then
            print '(a)', 'SKIP: AMD ordering unavailable in this build (fpm)'
            cycle
        end if
        call random_shifts(37, 0.15_dp, orderings(ordering))
        call random_shifts(150, 0.03_dp, orderings(ordering))
        call block_pencil([3, 5, 4, 6, 2, 7, 3, 5, 4, 6, 3, 2, 5, 4], &
                          orderings(ordering), 24)
        call saddle_point(12, 5, orderings(ordering))
    end do
    call exact_singular()
    call near_singular()
    call scaled_blocks()
    call extreme_pairs()
    call hermitian_shifts(30)
    call unsupported_families()
    print '(a)', 'native inertia: all checks passed'
contains
    ! The fpm build has no AMD ordering; CMake builds must provide it.
    logical function ordering_available(order) result(available)
        integer, intent(in) :: order
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(csc_t) :: a
        call csc_from_triplet(2, 2, [1, 2, 1, 2], [1, 1, 2, 2], &
                              [1.0_dp, 2.0_dp, 2.0_dp, -1.0_dp], a, status)
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        solver%options%ordering = order
        call sparse_factor(solver, a, status)
        available = status%code /= FORTSPARSE_BACKEND_UNAVAILABLE
        call sparse_destroy(solver)
    end function ordering_available

    subroutine random_seed_fixed()
        integer :: n
        integer, allocatable :: seed(:)
        call random_seed(size=n)
        allocate (seed(n))
        seed = 20261016
        call random_seed(put=seed)
    end subroutine random_seed_fixed

    subroutine expect(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine expect

    subroutine dense_to_csc(dense, mask, a)
        real(dp), intent(in) :: dense(:, :)
        logical, intent(in) :: mask(:, :)
        type(csc_t), intent(out) :: a
        type(fortsparse_status_t) :: status
        integer, allocatable :: rows(:), columns(:)
        real(dp), allocatable :: values(:)
        integer :: i, j, k, n
        n = size(dense, 1)
        k = count(mask)
        allocate (rows(k), columns(k), values(k))
        k = 0
        do j = 1, n
            do i = 1, n
                if (.not. mask(i, j)) cycle
                k = k + 1
                rows(k) = i
                columns(k) = j
                values(k) = dense(i, j)
            end do
        end do
        call csc_from_triplet(n, n, rows, columns, values, a, status)
        call expect(status%code == FORTSPARSE_OK, 'triplet assembly')
    end subroutine dense_to_csc

    subroutine symmetric_eigenvalues(dense, w)
        real(dp), intent(in) :: dense(:, :)
        real(dp), intent(out) :: w(:)
        real(dp), allocatable :: copy(:, :), work(:)
        integer :: n, info
        n = size(dense, 1)
        copy = dense
        allocate (work(max(1, 3*n)))
        call dsyev('N', 'L', n, copy, n, w, work, size(work), info)
        call expect(info == 0, 'dsyev oracle')
    end subroutine symmetric_eigenvalues

    ! Shift halfway between two neighbouring oracle eigenvalues: the exact
    ! inertia then has a gap of half the spacing on both sides.
    real(dp) function gap_shift(w, k) result(sigma)
        real(dp), intent(in) :: w(:)
        integer, intent(in) :: k
        sigma = 0.5_dp*(w(k) + w(k + 1))
    end function gap_shift

    subroutine query(solver, negative, zero, positive, reliability, status)
        type(sparse_solver_t), intent(inout) :: solver
        integer(i8), intent(out) :: negative, zero, positive
        integer, intent(out) :: reliability
        type(fortsparse_status_t), intent(out) :: status
        call sparse_inertia(solver, negative, zero, positive, reliability, status)
    end subroutine query

    subroutine check_counts(solver, expected_negative, n, label)
        type(sparse_solver_t), intent(inout) :: solver
        integer, intent(in) :: expected_negative, n
        character(*), intent(in) :: label
        type(fortsparse_status_t) :: status
        integer(i8) :: negative, zero, positive
        integer :: reliability
        call query(solver, negative, zero, positive, reliability, status)
        call expect(status%code == FORTSPARSE_OK, label//': inertia status')
        call expect(reliability == FORTSPARSE_INERTIA_RELIABLE, label//': reliability')
        call expect(negative == expected_negative, label//': negative count')
        call expect(zero == 0_i8 .and. positive == int(n - expected_negative, i8), &
                    label//': zero/positive counts')
    end subroutine check_counts

    subroutine random_shifts(n, density, order)
        integer, intent(in) :: n, order
        real(dp), intent(in) :: density
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(csc_t) :: a
        real(dp), allocatable :: base(:, :), shifted(:, :), w(:)
        logical, allocatable :: mask(:, :)
        real(dp) :: r
        integer :: i, j, k, trial
        logical :: reused
        allocate (base(n, n), mask(n, n), w(n))
        base = 0.0_dp
        mask = .false.
        do j = 1, n
            mask(j, j) = .true.
            call random_number(r)
            base(j, j) = 4.0_dp*r - 2.0_dp
            do i = j + 1, n
                call random_number(r)
                if (r > density) cycle
                call random_number(r)
                base(i, j) = 2.0_dp*r - 1.0_dp
                base(j, i) = base(i, j)
                mask(i, j) = .true.
                mask(j, i) = .true.
            end do
        end do
        ! Zero part of the diagonal to force 2x2 Bunch--Kaufman pivots.
        do j = 1, n, 3
            base(j, j) = 0.0_dp
        end do
        call symmetric_eigenvalues(base, w)
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        solver%options%ordering = order
        do trial = 1, 7
            k = max(1, min(n - 1, (trial*n)/8))
            shifted = base
            do j = 1, n
                shifted(j, j) = shifted(j, j) - gap_shift(w, k)
            end do
            call dense_to_csc(shifted, mask, a)
            if (trial == 1) then
                call sparse_factor(solver, a, status)
            else
                call sparse_refactor(solver, a, status, reused)
                call expect(reused, 'random shift: pattern reuse')
            end if
            call expect(status%code == FORTSPARSE_OK, 'random shift factor')
            call check_counts(solver, k, n, 'random symmetric shift')
        end do
        ! Unshifted matrix (zero diagonal entries): compare with oracle signs.
        call dense_to_csc(base, mask, a)
        call sparse_factor(solver, a, status)
        call expect(status%code == FORTSPARSE_OK, 'unshifted factor')
        if (minval(abs(w)) > 1.0e-8_dp) &
            call check_counts(solver, count(w < 0.0_dp), n, 'unshifted indefinite')
        call sparse_destroy(solver)
    end subroutine random_shifts

    ! GLISS-like pencil K - sigma M: variable-width dense blocks on a
    ! block-tridiagonal pattern, symmetric indefinite K and SPD M. The number
    ! of generalized eigenvalues below sigma equals the negative inertia.
    subroutine block_pencil(widths, order, nshift)
        integer, intent(in) :: widths(:), order, nshift
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(csc_t) :: a
        real(dp), allocatable :: k(:, :), m(:, :), shifted(:, :), w(:), work(:)
        real(dp), allocatable :: kc(:, :), mc(:, :)
        logical, allocatable :: mask(:, :)
        integer, allocatable :: start(:)
        real(dp) :: r
        integer :: n, b, i, j, info, s, below
        logical :: reused
        n = sum(widths)
        allocate (k(n, n), m(n, n), mask(n, n), w(n), start(size(widths) + 1))
        start(1) = 1
        do b = 1, size(widths)
            start(b + 1) = start(b) + widths(b)
        end do
        mask = .false.
        do b = 1, size(widths)
            mask(start(b):start(b + 1) - 1, start(b):start(b + 1) - 1) = .true.
            if (b < size(widths)) then
                mask(start(b + 1):start(b + 2) - 1, start(b):start(b + 1) - 1) = .true.
                mask(start(b):start(b + 1) - 1, start(b + 1):start(b + 2) - 1) = .true.
            end if
        end do
        k = 0.0_dp
        m = 0.0_dp
        do j = 1, n
            do i = j, n
                if (.not. mask(i, j)) cycle
                call random_number(r)
                k(i, j) = 2.0_dp*r - 1.0_dp
                k(j, i) = k(i, j)
                call random_number(r)
                m(i, j) = 0.1_dp*r
                m(j, i) = m(i, j)
            end do
        end do
        do j = 1, n
            ! Gershgorin: M is SPD by strict diagonal dominance.
            m(j, j) = sum(abs(m(:, j))) + 1.0_dp
        end do
        kc = k
        mc = m
        allocate (work(3*n))
        call dsygv(1, 'N', 'L', n, kc, n, mc, n, w, work, size(work), info)
        call expect(info == 0, 'dsygv oracle')
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        solver%options%ordering = order
        do s = 1, nshift
            below = 1 + ((s - 1)*(n - 2))/(nshift - 1)
            shifted = k - gap_shift(w, below)*m
            call dense_to_csc(shifted, mask, a)
            if (s == 1) then
                call sparse_factor(solver, a, status)
            else
                call sparse_refactor(solver, a, status, reused)
                call expect(reused, 'pencil shift: analysis reuse')
            end if
            call expect(status%code == FORTSPARSE_OK, 'pencil shift factor')
            call check_counts(solver, below, n, 'block-tridiagonal pencil')
        end do
        call sparse_destroy(solver)
    end subroutine block_pencil

    ! KKT matrix [[H, B^T], [B, 0]] with SPD H and full-rank B has inertia
    ! (p, 0, m); the zero block forces 2x2 pivots.
    subroutine saddle_point(m, p, order)
        integer, intent(in) :: m, p, order
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(csc_t) :: a
        real(dp) :: dense(m + p, m + p), r, w(m + p)
        logical :: mask(m + p, m + p)
        integer :: i, j
        dense = 0.0_dp
        mask = .false.
        do i = 1, m
            dense(i, i) = 1.0e-3_dp*i
            mask(i, i) = .true.
        end do
        do i = 1, p
            mask(m + i, m + i) = .true.
            do j = 1, m
                if (j /= i .and. j /= i + p .and. mod(i + j, 4) /= 0) cycle
                call random_number(r)
                dense(m + i, j) = 1.0_dp + r
                dense(j, m + i) = dense(m + i, j)
                mask(m + i, j) = .true.
                mask(j, m + i) = .true.
            end do
        end do
        call symmetric_eigenvalues(dense, w)
        call expect(count(w < 0.0_dp) == p .and. minval(abs(w)) > 1.0e-6_dp, &
                    'saddle-point oracle inertia')
        call dense_to_csc(dense, mask, a)
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        solver%options%ordering = order
        call sparse_factor(solver, a, status)
        call expect(status%code == FORTSPARSE_OK, 'saddle-point factor')
        call check_counts(solver, p, m + p, 'saddle-point 2x2 pivots')
        call sparse_destroy(solver)
    end subroutine saddle_point

    ! [[1, 2], [2, 4]] has an exactly representable zero Schur complement;
    ! the rest of the matrix is a diagonal with known signs.
    subroutine exact_singular()
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(csc_t) :: a
        real(dp) :: dense(6, 6), x(6), b(6)
        logical :: mask(6, 6)
        integer(i8) :: negative, zero, positive
        integer :: reliability, order
        dense = 0.0_dp
        dense(1:2, 1:2) = reshape([1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp], [2, 2])
        dense(3, 3) = -3.0_dp
        dense(4, 4) = 5.0_dp
        dense(5, 5) = -1.0_dp
        dense(6, 6) = 0.0_dp
        dense(5, 6) = 0.0_dp
        mask = dense /= 0.0_dp
        mask(6, 6) = .true.
        do order = 1, 2
            if (.not. ordering_available(orderings(order))) cycle
            call dense_to_csc(dense, mask, a)
            solver%backend_id = FORTSPARSE_BACKEND_NATIVE
            solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
            solver%options%ordering = orderings(order)
            call sparse_factor(solver, a, status)
            call expect(status%code == FORTSPARSE_SINGULAR, 'exact singular status')
            call query(solver, negative, zero, positive, reliability, status)
            call expect(status%code == FORTSPARSE_OK, 'singular inertia available')
            call expect(reliability == FORTSPARSE_INERTIA_SINGULAR, &
                        'singular reliability')
            ! Oracle eigenvalues: 0, 5 (pair), -3, 5, -1, 0.
            call expect(negative == 2_i8 .and. zero == 2_i8 .and. positive == 2_i8, &
                        'exact singular counts')
            b = 1.0_dp
            call sparse_solve(solver, b, x, status)
            call expect(status%code == FORTSPARSE_NOT_FACTORED, &
                        'singular solve refused')
            call sparse_destroy(solver)
        end do
    end subroutine exact_singular

    ! A Schur pivot of 2^-50 is far below n eps growth for n = 120: the sign
    ! is reported but flagged UNCERTAIN; the remaining counts stay exact.
    subroutine near_singular()
        integer, parameter :: n = 120
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(csc_t) :: a
        real(dp) :: dense(n, n)
        logical :: mask(n, n)
        integer(i8) :: negative, zero, positive
        integer :: reliability, j
        dense = 0.0_dp
        do j = 3, n
            dense(j, j) = merge(2.0_dp, -3.0_dp, mod(j, 3) == 0)
        end do
        dense(1:2, 1:2) = reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp + 2.0_dp**(-50)], &
                                  [2, 2])
        mask = dense /= 0.0_dp
        call dense_to_csc(dense, mask, a)
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        call sparse_factor(solver, a, status)
        call expect(status%code == FORTSPARSE_OK, 'near-singular factor')
        call query(solver, negative, zero, positive, reliability, status)
        call expect(status%code == FORTSPARSE_OK, 'near-singular inertia')
        call expect(reliability == FORTSPARSE_INERTIA_UNCERTAIN, 'near-singular flag')
        call expect(zero == 0_i8 .and. negative == int(count([(mod(j, 3) /= 0, &
                    j=3, n)]), i8) .and. negative + positive == int(n, i8), &
                    'near-singular counts')
        call sparse_destroy(solver)
    end subroutine near_singular

    ! Two decoupled blocks at scales 1e12 and 1e-8 (as for a pencil whose K
    ! has zero rows): the thresholds are local, so the small block's pivots
    ! stay RELIABLE. Oracle: separate dsyev counts of each block.
    subroutine scaled_blocks()
        integer, parameter :: nb = 30, ns = 10, n = nb + ns
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(csc_t) :: a
        real(dp) :: dense(n, n), big(nb, nb), small(ns, ns), wb(nb), ws(ns), r
        logical :: mask(n, n)
        integer :: i, j
        big = 0.0_dp
        small = 0.0_dp
        do j = 1, nb
            do i = j, min(nb, j + 2)
                call random_number(r)
                big(i, j) = 1.0e12_dp*(2.0_dp*r - 1.0_dp)
                big(j, i) = big(i, j)
            end do
        end do
        do j = 1, ns
            small(j, j) = merge(1.0e-8_dp, -2.0e-8_dp, mod(j, 3) == 0)
            if (j < ns) then
                small(j + 1, j) = 3.0e-9_dp
                small(j, j + 1) = 3.0e-9_dp
            end if
        end do
        call symmetric_eigenvalues(big, wb)
        call symmetric_eigenvalues(small, ws)
        call expect(minval(abs(wb)) > 1.0e-6_dp*maxval(abs(wb)) .and. &
                    minval(abs(ws)) > 1.0e-3_dp*maxval(abs(ws)), 'scaled oracle gaps')
        dense = 0.0_dp
        dense(1:nb, 1:nb) = big
        dense(nb + 1:, nb + 1:) = small
        mask = dense /= 0.0_dp
        call dense_to_csc(dense, mask, a)
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        call sparse_factor(solver, a, status)
        call expect(status%code == FORTSPARSE_OK, 'scaled blocks factor')
        call check_counts(solver, count(wb < 0.0_dp) + count(ws < 0.0_dp), n, &
                          'decoupled 1e12/1e-8 blocks')
        call sparse_destroy(solver)
    end subroutine scaled_blocks

    ! Extreme-scale 2x2 pivots: inertia must not overflow or underflow.
    subroutine extreme_pairs()
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(csc_t) :: a
        real(dp) :: dense(4, 4), scale
        logical :: mask(4, 4)
        integer :: s
        do s = 1, 2
            scale = merge(1.0e300_dp, 1.0e-300_dp, s == 1)
            dense = 0.0_dp
            dense(1, 2) = scale
            dense(2, 1) = scale
            dense(3, 4) = 2.0_dp*scale
            dense(4, 3) = 2.0_dp*scale
            dense(3, 3) = 1.0e-3_dp*scale
            mask = dense /= 0.0_dp
            call dense_to_csc(dense, mask, a)
            solver%backend_id = FORTSPARSE_BACKEND_NATIVE
            solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
            call sparse_factor(solver, a, status)
            call expect(status%code == FORTSPARSE_OK, 'extreme pair factor')
            call check_counts(solver, 2, 4, 'extreme-scale 2x2 pivots')
            call sparse_destroy(solver)
        end do
    end subroutine extreme_pairs

    subroutine hermitian_shifts(n)
        integer, intent(in) :: n
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(csc_z_t) :: a
        complex(dp), allocatable :: dense(:, :), copy(:, :), work(:), values(:)
        real(dp), allocatable :: w(:), rwork(:)
        integer, allocatable :: rows(:), columns(:)
        real(dp) :: r, q
        integer :: i, j, k, info, trial
        logical :: reused
        allocate (dense(n, n), w(n), rwork(3*n), work(4*n))
        dense = (0.0_dp, 0.0_dp)
        do j = 1, n
            call random_number(r)
            dense(j, j) = cmplx(merge(0.0_dp, 2.0_dp*r - 1.0_dp, mod(j, 4) == 0), &
                                0.0_dp, dp)
            do i = j + 1, n
                call random_number(r)
                if (r > 0.2_dp .and. i /= j + 1) cycle
                call random_number(r)
                call random_number(q)
                dense(i, j) = cmplx(2.0_dp*r - 1.0_dp, 2.0_dp*q - 1.0_dp, dp)
                dense(j, i) = conjg(dense(i, j))
            end do
        end do
        copy = dense
        call zheev('N', 'L', n, copy, n, w, work, size(work), rwork, info)
        call expect(info == 0, 'zheev oracle')
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry = FORTSPARSE_SYM_HERMITIAN
        do trial = 1, 5
            k = (trial*n)/6
            copy = dense
            do j = 1, n
                copy(j, j) = copy(j, j) - gap_shift(w, k)
            end do
            ! Keep the pattern (including the shifted diagonal) fixed.
            allocate (rows(0), columns(0), values(0))
            do j = 1, n
                do i = 1, n
                    if (i /= j .and. abs(copy(i, j)) <= 0.0_dp) cycle
                    rows = [rows, i]
                    columns = [columns, j]
                    values = [values, copy(i, j)]
                end do
            end do
            call csc_from_triplet(n, n, rows, columns, values, a, status)
            call expect(status%code == FORTSPARSE_OK, 'Hermitian assembly')
            deallocate (rows, columns, values)
            if (trial == 1) then
                call sparse_factor(solver, a, status)
            else
                call sparse_refactor(solver, a, status, reused)
            end if
            call expect(status%code == FORTSPARSE_OK, 'Hermitian factor')
            call check_counts(solver, k, n, 'Hermitian LDL^H shift')
        end do
        call sparse_destroy(solver)
    end subroutine hermitian_shifts

    subroutine unsupported_families()
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(csc_t) :: a
        type(csc_z_t) :: z
        real(dp) :: dense(3, 3)
        logical :: mask(3, 3)
        integer(i8) :: negative, zero, positive
        integer :: reliability
        dense = reshape([4.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, &
                         0.0_dp, 1.0_dp, 2.0_dp], [3, 3])
        mask = dense /= 0.0_dp
        call dense_to_csc(dense, mask, a)
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        call query(solver, negative, zero, positive, reliability, status)
        call expect(status%code == FORTSPARSE_NOT_FACTORED, 'inertia before factor')
        ! Real SPD Cholesky: all pivots positive on success.
        solver%options%symmetry = FORTSPARSE_SYM_SPD
        call sparse_factor(solver, a, status)
        call expect(status%code == FORTSPARSE_OK, 'SPD factor')
        call check_counts(solver, 0, 3, 'SPD Cholesky')
        call sparse_destroy(solver)
        ! General LU (scalar and statically perturbed supernodal): no inertia.
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry = FORTSPARSE_SYM_GENERAL
        call sparse_factor(solver, a, status)
        call expect(status%code == FORTSPARSE_OK, 'LU factor')
        call query(solver, negative, zero, positive, reliability, status)
        call expect(status%code == FORTSPARSE_BACKEND_UNAVAILABLE .and. &
                    reliability == FORTSPARSE_INERTIA_UNAVAILABLE .and. &
                    negative == -1_i8, 'LU inertia must be unsupported')
        call sparse_destroy(solver)
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry = FORTSPARSE_SYM_GENERAL
        solver%options%lu_method = FORTSPARSE_LU_SUPERNODAL
        call sparse_factor(solver, a, status)
        ! The fpm build has no supernodal LU; either way no inertia results.
        call expect(status%code == FORTSPARSE_OK .or. &
                    status%code == FORTSPARSE_BACKEND_UNAVAILABLE, 'supernodal LU')
        call query(solver, negative, zero, positive, reliability, status)
        call expect(status%code == FORTSPARSE_BACKEND_UNAVAILABLE .and. &
                    reliability == FORTSPARSE_INERTIA_UNAVAILABLE, &
                    'perturbed supernodal LU inertia must be unsupported')
        call sparse_destroy(solver)
        ! Complex symmetric (non-Hermitian) LDL^T has no inertia.
        call csc_from_triplet(2, 2, [1, 2, 1, 2], [1, 1, 2, 2], &
            [(1.0_dp, 1.0_dp), (2.0_dp, 0.0_dp), (2.0_dp, 0.0_dp), (0.0_dp, 1.0_dp)], &
            z, status)
        call expect(status%code == FORTSPARSE_OK, 'complex symmetric assembly')
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options = sparse_options_t()
        solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        call sparse_factor(solver, z, status)
        call expect(status%code == FORTSPARSE_OK, 'complex symmetric factor')
        call query(solver, negative, zero, positive, reliability, status)
        call expect(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
                    'complex symmetric inertia must be unsupported')
        call sparse_destroy(solver)
    end subroutine unsupported_families
end program test_fortsparse_inertia
