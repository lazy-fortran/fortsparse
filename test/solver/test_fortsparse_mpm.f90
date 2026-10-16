program test_fortsparse_mpm
    ! Maximum-product matching: exhaustive optimality for n <= 8, an independent
    ! dense Hungarian oracle for n = 200, LP certificates (dual feasibility,
    ! complementary slackness, zero duality gap), scaling bounds, invalid and
    ! singular inputs, and the caller workspace contract.
    use fortsparse, only: dp, fortsparse_status_t, status_ok, &
                          FORTSPARSE_SINGULAR, FORTSPARSE_INVALID_MATRIX
    use fortsparse_kinds, only: i8
    use fortsparse_native_mpm, only: mpm_stats_t, mpm_iwork_size, &
                                     mpm_rwork_size, mpm_workspace_bytes, &
                                     max_product_matching, mpm_scaling
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
                                             ieee_positive_inf
    implicit none
    integer, parameter :: NMAX = 8
    integer(i8) :: seed = 20261016_i8
    integer :: trial, n, singular_seen, optimal_seen
    real(dp) :: dense(NMAX, NMAX)
    logical :: stored(NMAX, NMAX)

    call check_workspace_sizes()
    singular_seen = 0; optimal_seen = 0
    do trial = 1, 600
        n = 1 + mod(trial - 1, NMAX)
        call random_dense(n, trial, dense, stored)
        call exhaustive_case(n, dense(1:n, 1:n), stored(1:n, 1:n), &
                             mod(trial, 3) == 0)
    end do
    call check(singular_seen > 20 .and. optimal_seen > 300, &
               'random corpus covers singular and regular matrices')
    call edge_cases()
    call invalid_and_singular_cases()
    call disconnected_blocks()
    call hungarian_regression(200)
    call large_certificate(2000)
    print '(a)', 'test_fortsparse_mpm: all checks passed'
contains
    subroutine check(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine

    real(dp) function uniform()
        ! Park--Miller minimal standard generator; deterministic everywhere.
        seed = mod(16807_i8*seed, 2147483647_i8)
        uniform = real(seed, dp)/2147483647.0_dp
    end function

    integer function pick(m)
        integer, intent(in) :: m
        pick = min(m, 1 + int(uniform()*m))
    end function

    subroutine random_dense(m, kind, a, mask)
        ! Kinds rotate through wide magnitudes, heavy ties and explicit zeros.
        integer, intent(in) :: m, kind
        real(dp), intent(out) :: a(:, :)
        logical, intent(out) :: mask(:, :)
        real(dp) :: density, sgn
        integer :: i, j
        density = 0.25_dp + 0.6_dp*uniform()
        a = 0; mask = .false.
        do j = 1, m
            do i = 1, m
                if (uniform() > density) cycle
                mask(i, j) = .true.
                sgn = merge(1.0_dp, -1.0_dp, uniform() < 0.5_dp)
                select case (mod(kind, 4))
                case (0)
                    a(i, j) = sgn*10.0_dp**(24*uniform() - 12)
                case (1)
                    a(i, j) = sgn*2.0_dp**pick(3)
                case (2)
                    a(i, j) = 1.0_dp
                case default
                    a(i, j) = sgn*10.0_dp**(24*uniform() - 12)
                    if (uniform() < 0.2_dp) a(i, j) = 0
                end select
            end do
        end do
    end subroutine

    subroutine to_csc(m, a, mask, ptr, rows, values)
        integer, intent(in) :: m
        real(dp), intent(in) :: a(:, :)
        logical, intent(in) :: mask(:, :)
        integer, allocatable, intent(out) :: ptr(:), rows(:)
        real(dp), allocatable, intent(out) :: values(:)
        integer :: i, j, k
        allocate (ptr(m + 1), rows(count(mask)), values(count(mask)))
        k = 0; ptr(1) = 1
        do j = 1, m
            do i = 1, m
                if (.not. mask(i, j)) cycle
                k = k + 1; rows(k) = i; values(k) = a(i, j)
            end do
            ptr(j + 1) = k + 1
        end do
    end subroutine

    recursive subroutine search(m, logs, ok, col, used, partial, best)
        ! Exhaustive enumeration of every nonzero perfect matching.
        integer, intent(in) :: m, col
        real(dp), intent(in) :: logs(:, :), partial
        logical, intent(in) :: ok(:, :)
        logical, intent(inout) :: used(:)
        real(dp), intent(inout) :: best
        integer :: i
        if (col > m) then
            best = max(best, partial)
            return
        end if
        do i = 1, m
            if (used(i) .or. .not. ok(i, col)) cycle
            used(i) = .true.
            call search(m, logs, ok, col + 1, used, partial + logs(i, col), best)
            used(i) = .false.
        end do
    end subroutine

    real(dp) function brute_force(m, a, mask) result(best)
        ! Maximum sum of log|a| over all permutations; -huge when none exists.
        integer, intent(in) :: m
        real(dp), intent(in) :: a(:, :)
        logical, intent(in) :: mask(:, :)
        real(dp) :: logs(m, m)
        logical :: ok(m, m), used(m)
        ok = mask(1:m, 1:m) .and. a(1:m, 1:m) /= 0
        logs = 0
        where (ok) logs = log(abs(a(1:m, 1:m)))
        used = .false.; best = -huge(1.0_dp)
        call search(m, logs, ok, 1, used, 0.0_dp, best)
    end function

    subroutine run_real(m, ptr, rows, values, cache, cmatch, u, v, status, stats)
        integer, intent(in) :: m, ptr(:), rows(:)
        real(dp), intent(in) :: values(:)
        logical, intent(in) :: cache
        integer, intent(out) :: cmatch(:)
        real(dp), intent(out) :: u(:), v(:)
        type(fortsparse_status_t), intent(out) :: status
        type(mpm_stats_t), intent(out) :: stats
        integer, allocatable :: iwork(:)
        real(dp), allocatable :: rwork(:)
        integer :: k
        allocate (iwork(mpm_iwork_size(m)))
        allocate (rwork(mpm_rwork_size(m, size(rows), cache)))
        call check(4_i8*size(iwork, kind=i8) + 8_i8*size(rwork, kind=i8) == &
                   mpm_workspace_bytes(m, size(rows), cache), &
                   'workspace byte identity')
        call max_product_matching(m, ptr, rows, values, cmatch, u, v, iwork, &
                                  rwork, status, stats)
        if (status_ok(status)) then
            call check(all(iwork(1:m) > 0), 'row inverse complete')
            call check(all(iwork(cmatch) == [(k, k=1, m)]), 'row inverse')
        end if
    end subroutine

    subroutine certificate(m, ptr, rows, mags, cmatch, u, v, objective)
        ! Independent LP certificate of optimality and the scaling property.
        integer, intent(in) :: m, ptr(:), rows(:), cmatch(:)
        real(dp), intent(in) :: mags(:), u(:), v(:)
        real(dp), intent(out) :: objective
        real(dp) :: colmax, c, red, tol, dual, scaled, logprod
        logical :: hit(m), found
        integer :: i, j, e
        hit = .false.
        do j = 1, m
            call check(cmatch(j) >= 1 .and. cmatch(j) <= m, 'matched row range')
            call check(.not. hit(cmatch(j)), 'matching is a permutation')
            hit(cmatch(j)) = .true.
        end do
        objective = 0; logprod = 0; dual = sum(u) + sum(v)
        tol = 1e-12_dp*(1 + maxval(abs(u)) + maxval(abs(v)))*max(1, m)
        do j = 1, m
            colmax = maxval(mags(ptr(j):ptr(j + 1) - 1))
            found = .false.
            do e = ptr(j), ptr(j + 1) - 1
                if (mags(e) == 0) cycle
                i = rows(e)
                c = log(colmax) - log(mags(e))
                red = c - u(i) - v(j)
                call check(red >= -tol, 'dual feasibility')
                scaled = exp(u(i))*mags(e)*(exp(v(j))/colmax)
                call check(scaled <= 1 + 1e-12_dp*max(1, m), 'scaled bound')
                if (i == cmatch(j)) then
                    found = .true.
                    call check(abs(red) <= tol, 'complementary slackness')
                    call check(abs(scaled - 1) <= 1e-11_dp*max(1, m), &
                               'matched scaled entry is one')
                    objective = objective + c
                    logprod = logprod + log(mags(e))
                end if
            end do
            call check(found, 'matched entry is numerically nonzero')
        end do
        call check(abs(objective - dual) <= tol*m, 'zero duality gap')
        objective = logprod
    end subroutine

    subroutine check_scaling(m, ptr, rows, values, cmatch, u, v)
        integer, intent(in) :: m, ptr(:), rows(:), cmatch(:)
        real(dp), intent(in) :: values(:), u(:), v(:)
        real(dp) :: dr(m), dc(m), s
        type(fortsparse_status_t) :: status
        integer :: j, e
        call mpm_scaling(m, ptr, rows, values, u, v, dr, dc, status)
        call check(status_ok(status), 'scaling status')
        do j = 1, m
            do e = ptr(j), ptr(j + 1) - 1
                s = abs(dr(rows(e))*values(e)*dc(j))
                call check(s <= 1 + 1e-12_dp*max(1, m), 'mpm_scaling bound')
                if (rows(e) == cmatch(j)) then
                    call check(abs(s - 1) <= 1e-11_dp*max(1, m), &
                               'mpm_scaling matched entry')
                end if
            end do
        end do
    end subroutine

    subroutine exhaustive_case(m, a, mask, with_complex)
        integer, intent(in) :: m
        real(dp), intent(in) :: a(:, :)
        logical, intent(in) :: mask(:, :), with_complex
        integer, allocatable :: ptr(:), rows(:), iwork(:)
        real(dp), allocatable :: values(:), rwork(:)
        complex(dp), allocatable :: z(:)
        integer :: cmatch(m), cm2(m), e, j
        real(dp) :: u(m), v(m), u2(m), v2(m), best, objective, dr(m), dc(m)
        type(fortsparse_status_t) :: status, status2
        type(mpm_stats_t) :: stats, stats2
        call to_csc(m, a, mask, ptr, rows, values)
        best = brute_force(m, a, mask)
        call run_real(m, ptr, rows, values, .true., cmatch, u, v, status, stats)
        call run_real(m, ptr, rows, values, .false., cm2, u2, v2, status2, stats2)
        call check(status%code == status2%code, 'cached/uncached status agree')
        if (best == -huge(1.0_dp)) then
            call check(status%code == FORTSPARSE_SINGULAR, &
                       'no nonzero perfect matching gives SINGULAR')
            singular_seen = singular_seen + 1
            return
        end if
        call check(status_ok(status), 'nonsingular matrix matched')
        optimal_seen = optimal_seen + 1
        call check(all(cmatch == cm2) .and. all(u == u2) .and. all(v == v2), &
                   'cached and recomputed costs are bitwise identical')
        call check(stats%greedy_matched + stats%augmentations == m, &
                   'statistics account for every column')
        call certificate(m, ptr, rows, abs(values), cmatch, u, v, objective)
        call check(abs(objective - best) <= 1e-9_dp*(1 + abs(best)), &
                   'exhaustive optimality (real)')
        call check_scaling(m, ptr, rows, values, cmatch, u, v)
        if (.not. with_complex) return
        ! Complex entries with the same moduli and arbitrary phases.
        allocate (z(size(values)))
        do e = 1, size(values)
            z(e) = values(e)*cmplx(cos(7.0_dp*e), sin(7.0_dp*e), dp)
        end do
        allocate (iwork(mpm_iwork_size(m)), rwork(mpm_rwork_size(m, size(z), &
                                                                 .false.)))
        call max_product_matching(m, ptr, rows, z, cm2, u2, v2, iwork, rwork, &
                                  status, stats)
        call check(status_ok(status), 'complex matched')
        call certificate(m, ptr, rows, abs(z), cm2, u2, v2, objective)
        call check(abs(objective - best) <= 1e-9_dp*(1 + abs(best)), &
                   'exhaustive optimality (complex modulus)')
        call mpm_scaling(m, ptr, rows, z, u2, v2, dr, dc, status)
        call check(status_ok(status), 'complex scaling status')
        do j = 1, m
            do e = ptr(j), ptr(j + 1) - 1
                call check(abs(dr(rows(e))*z(e)*dc(j)) <= 1 + 1e-12_dp*m, &
                           'complex scaling bound')
                if (rows(e) == cm2(j)) then
                    call check(abs(abs(dr(rows(e))*z(e)*dc(j)) - 1) <= &
                               1e-11_dp*m, 'complex matched scaled entry')
                end if
            end do
        end do
    end subroutine

    subroutine check_workspace_sizes()
        call check(mpm_workspace_bytes(10, 50, .true.) == 4*60 + 8*(20 + 50), &
                   'cached workspace bytes')
        call check(mpm_workspace_bytes(10, 50, .false.) == 4*60 + 8*20, &
                   'uncached workspace bytes')
        call check(mpm_workspace_bytes(0, 0, .true.) == 0, 'empty workspace')
        call check(mpm_workspace_bytes(-1, 0, .false.) == -1, 'negative n')
        call check(mpm_workspace_bytes(10, -1, .false.) == -1, 'negative nnz')
        call check(mpm_workspace_bytes(huge(0), 0, .false.) == -1, &
                   'iwork index overflow')
        call check(mpm_workspace_bytes(huge(0)/6, huge(0), .true.) == -1, &
                   'rwork index overflow')
        call check(mpm_workspace_bytes(300000000, 0, .false.) == &
                   4_i8*1800000000_i8 + 8_i8*600000000_i8, 'large exact bytes')
    end subroutine

    subroutine edge_cases()
        integer :: ptr0(1), rows0(0), cm0(0), iw0(0), ptr1(2), rows1(1), cm1(1)
        integer :: iw1(6)
        real(dp) :: v0(0), u0(0), w0(0), rw0(0), u1(1), v1(1), rw1(3), dr1(1), dc1(1)
        type(fortsparse_status_t) :: status
        type(mpm_stats_t) :: stats
        ptr0 = 1
        call max_product_matching(0, ptr0, rows0, v0, cm0, u0, w0, iw0, rw0, &
                                  status, stats)
        call check(status_ok(status), 'n = 0 is valid')
        call mpm_scaling(0, ptr0, rows0, v0, u0, w0, rw0, rw0, status)
        call check(status_ok(status), 'n = 0 scaling')
        ptr1 = [1, 2]; rows1 = 1
        call max_product_matching(1, ptr1, rows1, [-3.0_dp], cm1, u1, v1, iw1, &
                                  rw1, status, stats)
        call check(status_ok(status) .and. cm1(1) == 1, 'n = 1 match')
        call check(abs(u1(1) + v1(1)) <= 0, 'n = 1 zero cost')
        call check(stats%greedy_matched == 1, 'n = 1 greedy')
        call mpm_scaling(1, ptr1, rows1, [-3.0_dp], u1, v1, dr1, dc1, status)
        call check(abs(abs(-3*dr1(1)*dc1(1)) - 1) < 1e-15_dp, 'n = 1 scaling')
        call max_product_matching(1, ptr1, rows1, [0.0_dp], cm1, u1, v1, iw1, &
                                  rw1, status, stats)
        call check(status%code == FORTSPARSE_SINGULAR, 'n = 1 explicit zero')
    end subroutine

    subroutine run3(ptr, rows, values, expected, message, text)
        integer, intent(in) :: ptr(:), rows(:), expected
        real(dp), intent(in) :: values(:)
        character(*), intent(in) :: message
        character(*), intent(in), optional :: text
        integer :: cm(3)
        real(dp) :: u(3), v(3)
        type(fortsparse_status_t) :: status
        type(mpm_stats_t) :: stats
        call run_real(3, ptr, rows, values, .false., cm, u, v, status, stats)
        if (status%code /= expected) print '(a)', trim(status%msg)
        call check(status%code == expected, message)
        if (present(text)) call check(index(status%msg, text) > 0, &
                                      message//': explicit message')
    end subroutine

    subroutine invalid_and_singular_cases()
        real(dp) :: nan, inf
        integer :: cm(3), iw(18)
        real(dp) :: u(3), v(3), rw(6), rwbad(7)
        complex(dp) :: z(5)
        type(fortsparse_status_t) :: status
        type(mpm_stats_t) :: stats
        nan = ieee_value(nan, ieee_quiet_nan)
        inf = ieee_value(inf, ieee_positive_inf)
        ! Every row and column is nonzero, but rows 2 and 3 share one column.
        call run3([1, 2, 3, 6], [1, 1, 1, 2, 3], [1.0_dp, 2.0_dp, 3.0_dp, &
                  4.0_dp, 5.0_dp], FORTSPARSE_SINGULAR, 'Hall violation', &
                  'augmenting path from column 2 (2 columns matched before)')
        call run3([1, 3, 5, 7], [1, 2, 2, 3, 1, 3], [1.0_dp, 0.0_dp, 0.0_dp, &
                  0.0_dp, 1.0_dp, 1.0_dp], FORTSPARSE_SINGULAR, &
                  'column of explicit zeros', 'column 2 has no numerically')
        call run3([1, 3, 5, 7], [1, 3, 1, 3, 1, 3], [1.0_dp, 2.0_dp, 3.0_dp, &
                  4.0_dp, 5.0_dp, 6.0_dp], FORTSPARSE_SINGULAR, 'empty row', &
                  'row 2 has no numerically')
        ! Zeros that would complete a perfect matching are excluded.
        call run3([1, 3, 5, 7], [1, 2, 1, 2, 2, 3], [1.0_dp, 1.0_dp, 1.0_dp, &
                  1.0_dp, 1.0_dp, 0.0_dp], FORTSPARSE_SINGULAR, &
                  'numerically zero matching edge')
        call run3([1, 2, 3, 4], [1, 2, 3], [1.0_dp, nan, 1.0_dp], &
                  FORTSPARSE_INVALID_MATRIX, 'NaN')
        call run3([1, 2, 3, 4], [1, 2, 3], [1.0_dp, -inf, 1.0_dp], &
                  FORTSPARSE_INVALID_MATRIX, 'Inf')
        call run3([1, 3, 4, 5], [2, 1, 2, 3], [1.0_dp, 1.0_dp, 1.0_dp, &
                  1.0_dp], FORTSPARSE_INVALID_MATRIX, 'unsorted rows')
        call run3([1, 2, 3, 4], [1, 4, 3], [1.0_dp, 1.0_dp, 1.0_dp], &
                  FORTSPARSE_INVALID_MATRIX, 'row out of range')
        z = [(1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), (1.0_dp, 1.0_dp), &
             (2.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)]
        z(3) = cmplx(1.0_dp, nan, dp)
        call max_product_matching(3, [1, 2, 4, 6], [1, 1, 2, 2, 3], z, cm, u, v, &
                                  iw, rw, status, stats)
        call check(status%code == FORTSPARSE_INVALID_MATRIX, 'complex NaN')
        z(3) = cmplx(huge(1.0_dp), huge(1.0_dp), dp)
        call max_product_matching(3, [1, 2, 4, 6], [1, 1, 2, 2, 3], z, cm, u, v, &
                                  iw, rw, status, stats)
        call check(status%code == FORTSPARSE_INVALID_MATRIX, 'complex overflow')
        z(3) = (1.0_dp, 1.0_dp)
        call max_product_matching(3, [1, 2, 4, 6], [1, 1, 2, 2, 3], z, cm, u, v, &
                                  iw, rw, status, stats)
        call check(status%code == FORTSPARSE_SINGULAR, 'complex zero excluded')
        call max_product_matching(3, [1, 2, 3, 4], [1, 2, 3], [1.0_dp, 1.0_dp, &
                                  1.0_dp], cm, u, v, iw(1:17), rw, status, stats)
        call check(status%code == FORTSPARSE_INVALID_MATRIX, 'short iwork')
        call max_product_matching(3, [1, 2, 3, 4], [1, 2, 3], [1.0_dp, 1.0_dp, &
                                  1.0_dp], cm, u, v, iw, rwbad, status, stats)
        call check(status%code == FORTSPARSE_INVALID_MATRIX, 'odd rwork size')
        call max_product_matching(3, [1, 2, 3, 4], [1, 2, 3], [2.0_dp, 1.0_dp, &
                                  1.0_dp], cm, u, v, iw, rwbad(1:6), status, stats)
        call check(status_ok(status) .and. all(cm == [1, 2, 3]), 'diagonal')
    end subroutine

    subroutine disconnected_blocks()
        ! Two blocks scattered by a symmetric permutation: the optimum is the
        ! sum of the block optima and no edge crosses the blocks.
        integer, parameter :: m1 = 4, m2 = 5, m = m1 + m2
        integer :: p(m), i, j, attempt
        real(dp) :: a1(NMAX, NMAX), a2(NMAX, NMAX), a(m, m), best1, best2, obj
        logical :: s1(NMAX, NMAX), s2(NMAX, NMAX), s(m, m)
        integer, allocatable :: ptr(:), rows(:)
        real(dp), allocatable :: values(:)
        integer :: cm(m)
        real(dp) :: u(m), v(m)
        type(fortsparse_status_t) :: status
        type(mpm_stats_t) :: stats
        p = [3, 7, 1, 9, 5, 2, 8, 4, 6]
        do attempt = 1, 20
            best1 = -huge(1.0_dp); best2 = best1
            do while (best1 == -huge(1.0_dp))
                call random_dense(m1, 0, a1, s1)
                best1 = brute_force(m1, a1, s1)
            end do
            do while (best2 == -huge(1.0_dp))
                call random_dense(m2, 3, a2, s2)
                best2 = brute_force(m2, a2, s2)
            end do
            a = 0; s = .false.
            do j = 1, m1
                do i = 1, m1
                    a(p(i), p(j)) = a1(i, j); s(p(i), p(j)) = s1(i, j)
                end do
            end do
            do j = 1, m2
                do i = 1, m2
                    a(p(m1 + i), p(m1 + j)) = a2(i, j)
                    s(p(m1 + i), p(m1 + j)) = s2(i, j)
                end do
            end do
            call to_csc(m, a, s, ptr, rows, values)
            call run_real(m, ptr, rows, values, .true., cm, u, v, status, stats)
            call check(status_ok(status), 'block matrix matched')
            call certificate(m, ptr, rows, abs(values), cm, u, v, obj)
            call check(abs(obj - best1 - best2) <= 1e-9_dp*(1 + abs(obj)), &
                       'block optimum is the sum of block optima')
            do j = 1, m
                call check((findloc(p(1:m1), j, 1) > 0) .eqv. &
                           (findloc(p(1:m1), cm(j), 1) > 0), 'no cross-block edge')
            end do
        end do
    end subroutine

    subroutine hungarian(m, cost, assign)
        ! Dense O(m^3) shortest-augmenting-path Hungarian method on potentials,
        ! row by row (an independent, deliberately simple oracle).
        integer, intent(in) :: m
        real(dp), intent(in) :: cost(:, :)
        integer, intent(out) :: assign(:)
        real(dp) :: pu(0:m), pv(0:m), minv(0:m), delta, cur
        integer :: p(0:m), way(0:m), i, j, j0, j1, i0
        logical :: used(0:m)
        pu = 0; pv = 0; p = 0; way = 0
        do i = 1, m
            p(0) = i; j0 = 0; minv = huge(1.0_dp); used = .false.
            do
                used(j0) = .true.; i0 = p(j0); delta = huge(1.0_dp); j1 = 0
                do j = 1, m
                    if (used(j)) cycle
                    cur = cost(i0, j) - pu(i0) - pv(j)
                    if (cur < minv(j)) then
                        minv(j) = cur; way(j) = j0
                    end if
                    if (minv(j) < delta) then
                        delta = minv(j); j1 = j
                    end if
                end do
                do j = 0, m
                    if (used(j)) then
                        pu(p(j)) = pu(p(j)) + delta; pv(j) = pv(j) - delta
                    else
                        minv(j) = minv(j) - delta
                    end if
                end do
                j0 = j1
                if (p(j0) == 0) exit
            end do
            do
                j1 = way(j0); p(j0) = p(j1); j0 = j1
                if (j0 == 0) exit
            end do
        end do
        assign = p(1:m)
    end subroutine

    subroutine sparse_random(m, per_column, a, mask)
        ! Guaranteed nonzero transversal on a random permutation plus random
        ! off-transversal entries; magnitudes 1e-12..1e12.
        integer, intent(in) :: m, per_column
        real(dp), intent(out) :: a(:, :)
        logical, intent(out) :: mask(:, :)
        integer :: perm(m), i, j, k, t
        perm = [(i, i=1, m)]
        do i = m, 2, -1
            k = pick(i); t = perm(i); perm(i) = perm(k); perm(k) = t
        end do
        a = 0; mask = .false.
        do j = 1, m
            mask(perm(j), j) = .true.
            a(perm(j), j) = 10.0_dp**(24*uniform() - 12)
            do k = 1, per_column
                i = pick(m)
                mask(i, j) = .true.
                a(i, j) = merge(1, -1, uniform() < 0.5_dp)* &
                          10.0_dp**(24*uniform() - 12)
            end do
        end do
    end subroutine

    subroutine hungarian_regression(m)
        integer, intent(in) :: m
        real(dp), allocatable :: a(:, :), cost(:, :), values(:)
        logical, allocatable :: mask(:, :)
        integer, allocatable :: ptr(:), rows(:)
        integer :: cm(m), assign(m), i, j
        real(dp) :: u(m), v(m), obj, oracle
        real(dp), parameter :: PENALTY = 1e6_dp
        type(fortsparse_status_t) :: status
        type(mpm_stats_t) :: stats
        allocate (a(m, m), cost(m, m), mask(m, m))
        call sparse_random(m, 6, a, mask)
        cost = PENALTY
        do j = 1, m
            do i = 1, m
                if (mask(i, j) .and. a(i, j) /= 0) cost(i, j) = -log(abs(a(i, j)))
            end do
        end do
        call hungarian(m, cost, assign)
        oracle = 0
        do j = 1, m
            call check(cost(assign(j), j) < PENALTY, 'oracle uses nonzero edges')
            oracle = oracle - cost(assign(j), j)
        end do
        call to_csc(m, a, mask, ptr, rows, values)
        call run_real(m, ptr, rows, values, .true., cm, u, v, status, stats)
        call check(status_ok(status), 'n = 200 matched')
        call certificate(m, ptr, rows, abs(values), cm, u, v, obj)
        call check(abs(obj - oracle) <= 1e-9_dp*(1 + abs(oracle)), &
                   'Hungarian oracle objective')
        call check_scaling(m, ptr, rows, values, cm, u, v)
    end subroutine

    subroutine large_certificate(m)
        ! Deterministic 2000 x 2000 unsymmetric case; the independent LP
        ! certificate proves optimality without enumerating matchings.
        integer, intent(in) :: m
        integer, allocatable :: ptr(:), rows(:), perm(:)
        real(dp), allocatable :: values(:)
        integer :: cm(m), cm2(m), i, j, k, t, e
        real(dp) :: u(m), v(m), u2(m), v2(m), obj
        type(fortsparse_status_t) :: status
        type(mpm_stats_t) :: stats, stats2
        allocate (ptr(m + 1), rows(6*m), values(6*m), perm(m))
        perm = [(i, i=1, m)]
        do i = m, 2, -1
            k = pick(i); t = perm(i); perm(i) = perm(k); perm(k) = t
        end do
        e = 0; ptr(1) = 1
        do j = 1, m
            ! Sorted rows: the transversal row plus five strided neighbours.
            do k = 0, 5
                i = 1 + mod(perm(j) - 1 + k*(1 + mod(7*j, 97)), m)
                if (any(rows(ptr(j):e) == i)) cycle
                e = e + 1; rows(e) = i
                values(e) = merge(1, -1, uniform() < 0.5_dp)* &
                            10.0_dp**(24*uniform() - 12)
            end do
            call sort_column(rows(ptr(j):e), values(ptr(j):e))
            ptr(j + 1) = e + 1
        end do
        call run_real(m, ptr, rows(1:e), values(1:e), .true., cm, u, v, status, &
                      stats)
        call check(status_ok(status), 'n = 2000 matched')
        call run_real(m, ptr, rows(1:e), values(1:e), .false., cm2, u2, v2, &
                      status, stats2)
        call check(all(cm == cm2) .and. all(u == u2) .and. all(v == v2), &
                   'n = 2000 cached/uncached identical')
        call check(stats%heap_pops == stats2%heap_pops .and. &
                   stats%relaxations == stats2%relaxations, 'identical work')
        call certificate(m, ptr, rows(1:e), abs(values(1:e)), cm, u, v, obj)
        call check_scaling(m, ptr, rows(1:e), values(1:e), cm, u, v)
        call check(stats%log_range_dr > 0 .and. stats%log_range_dc > 0, &
                   'scaling ranges reported')
    end subroutine

    subroutine sort_column(r, x)
        integer, intent(inout) :: r(:)
        real(dp), intent(inout) :: x(:)
        integer :: i, k, t
        real(dp) :: y
        do i = 2, size(r)
            t = r(i); y = x(i); k = i - 1
            do while (k >= 1)
                if (r(k) <= t) exit
                r(k + 1) = r(k); x(k + 1) = x(k); k = k - 1
            end do
            r(k + 1) = t; x(k + 1) = y
        end do
    end subroutine
end program
