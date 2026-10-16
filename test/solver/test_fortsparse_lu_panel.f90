program test_fortsparse_lu_panel
    use fortsparse_kinds, only: dp, i8
    use fortsparse_status
    use fortsparse_csc, only: csc_t
    use fortsparse_native_lu, only: lu_factor_t, lu_panel_prefix, lu_prune_column
    use fortsparse_native_lu_z, only: lu_factor_z_t
    use fortsparse_spd_leaves, only: spd_panel_scatter, spd_panel_scatter_jvp, &
                                     spd_panel_scatter_vjp
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    call kernel_oracle()
    call closure_oracle()
    call factor_oracles()
    call memory_and_failure_oracles()
    print '(a)', &
        'Real LU panels: kernels, closure, pivots, reuse and caps pass'
contains
    subroutine checked(s)
        type(fortsparse_status_t), intent(in) :: s
        if (s%code /= FORTSPARSE_OK) then
            print *, s%code, trim(s%msg)
            error stop 'unexpected panel status'
        end if
    end subroutine checked

    subroutine kernel_oracle()
        integer, parameter :: n = 11, m = 4, rows(m) = [9, 2, 7, 4]
        real(dp) :: x(8, n), original(8, n), y(8, n), dx(8, n), dx0(8, n)
        real(dp) :: seed(8, n), xb(8, n), v(m), dv(m), c(8), dc(8), vb(m), cb(8)
        real(dp) :: left, right, expected, tangent
        integer :: i, j, mask, active, marks(n), checkmarks(n), lanes
        do j = 1, n
            do i = 1, 8
                original(i, j) = real(3*i - j, dp)/17
                dx0(i, j) = real(i + 2*j, dp)/31
                seed(i, j) = real(2*i - j, dp)/19
            end do
        end do
        v = [0.5_dp, -2.0_dp, 0.0_dp, 1.25_dp]
        dv = [0.25_dp, 0.3_dp, -0.2_dp, 0.1_dp]
        do i = 1, 8
            c(i) = real(i - 3, dp)/7
            dc(i) = real(i + 1, dp)/13
        end do
        do lanes = 1, 8
            do mask = 0, 255
                active = mask; marks = 0; checkmarks = 0
                x = original; dx = dx0; y = original
                call spd_panel_scatter(x, marks, rows, v, m, c, lanes, active)
                call spd_panel_scatter_jvp(y, dx, checkmarks, rows, v, dv, &
                    m, c, dc, lanes, active)
                if (any(marks /= checkmarks)) &
                    error stop 'panel primal/JVP graph mismatch'
                do j = 1, n
                    do i = 1, 8
                        expected = original(i, j); tangent = dx0(i, j)
                        if (any(rows == j) .and. i <= lanes) then
                            if (btest(mask, i - 1)) then
                                block
                                    integer :: k
                                    do k = 1, m
                                        if (rows(k) /= j) cycle
                                        expected = original(i, j) - v(k)*c(i)
                                        tangent = dx0(i, j) - dv(k)*c(i) - v(k)*dc(i)
                                    end do
                                end block
                            end if
                        end if
                        if (abs(x(i, j) - expected) > 2e-15_dp .or. &
                            abs(y(i, j) - expected) > 2e-15_dp) &
                            error stop 'panel masked value/padding oracle'
                        if (abs(dx(i, j) - tangent) > 2e-15_dp) &
                            error stop 'panel tangent oracle'
                    end do
                end do
                xb = seed; vb = 0; cb = 0
                call spd_panel_scatter_vjp(original, xb, rows, v, m, c, &
                    lanes, active, vb, cb)
                left = sum(seed*dx)
                right = sum(xb*dx0) + sum(vb*dv) + sum(cb*dc)
                if (abs(left - right) > 2e-13_dp) error stop 'panel VJP duality'
            end do
        end do
        marks = 4; x = original
        call spd_panel_scatter(x, marks, rows, v, 0, c, 8, 255)
        if (any(x /= original) .or. any(marks /= 4)) error stop 'empty panel scatter'
    end subroutine kernel_oracle

    function closure(graph) result(r)
        logical, intent(in) :: graph(:, :)
        logical :: r(size(graph, 1), size(graph, 2))
        integer :: i, j, k
        r = graph
        do i = 1, size(r, 1)
            r(i, i) = .true.
        end do
        do k = 1, size(r, 1)
            do j = 1, size(r, 1)
                do i = 1, size(r, 1)
                    r(i, j) = r(i, j) .or. (r(i, k) .and. r(k, j))
                end do
            end do
        end do
    end function closure

    recursive subroutine graph_visit(row, graph, seen, reach, count)
        integer, intent(in) :: row
        logical, intent(in) :: graph(:, :)
        logical, intent(inout) :: seen(:)
        integer, intent(inout) :: reach(:), count
        integer :: child
        if (seen(row)) return
        seen(row) = .true.
        do child = 1, size(seen)
            if (graph(row, child)) call graph_visit(child, graph, seen, reach, count)
        end do
        count = count + 1
        reach(count) = row
    end subroutine graph_visit

    subroutine closure_oracle()
        integer, parameter :: n = 5, m = n*(n - 1)/2, p(n) = [3, 5, 1, 4, 2]
        integer :: lp(n + 1), rows(m), ends(n), marks(n), inv(n), reach(n)
        integer :: mask, bit, i, j, k, e, next, count, lane, root, expected(n)
        real(dp) :: values(m), dense(n, n), f(8, n), reference(8, n), coeff(8)
        logical :: full(n, n), reduced(n, n), actual(n, n), oracle(n, n), seen(n)
        do i = 1, n
            inv(p(i)) = i
        end do
        ! All 1024 DAGs, nontrivial row labels, and actual paired pruning swaps.
        do mask = 0, 2**m - 1
            full = .false.; dense = 0; next = 1; bit = 0
            do j = 1, n
                lp(j) = next
                do i = j + 1, n
                    if (btest(mask, bit)) then
                        rows(next) = p(i)
                        values(next) = real(2*j - i, dp)/16
                        ! Keep graph entries nonzero, including cancellation later.
                        if (values(next) == 0) values(next) = 0.125_dp
                        dense(p(i), p(j)) = values(next)
                        full(p(j), p(i)) = .true.
                        next = next + 1
                    end if
                    bit = bit + 1
                end do
            end do
            lp(n + 1) = next; ends = lp(2:); marks = 0
            do j = 1, n
                do e = lp(j), lp(j + 1) - 1
                    marks(rows(e)) = -j
                end do
                do k = 1, j - 1
                    call lu_prune_column(rows, values, lp(k), lp(k + 1), p(j), marks, &
                                         j, ends(k))
                end do
            end do
            reduced = .false.
            do j = 1, n
                do e = lp(j), abs(ends(j)) - 1
                    reduced(p(j), rows(e)) = .true.
                end do
                do e = lp(j), lp(j + 1) - 1
                    if (values(e) /= dense(rows(e), p(j))) &
                        error stop 'pruned L tag misalignment'
                end do
            end do
            actual = closure(reduced); oracle = closure(full)
            if (any(actual .neqv. oracle)) error stop 'pruned graph closure differs'
            f = 0; marks = 0; expected = 0; seen = .false.; count = 0
            do lane = 1, 3
                root = p(2*lane - 1)
                f(lane, root) = real(lane, dp)/3
                marks(root) = ibset(marks(root), lane - 1)
                call graph_visit(root, reduced, seen, reach, count)
                do i = 1, n
                    if (oracle(root, i)) expected(i) = ibset(expected(i), lane - 1)
                end do
            end do
            reference = f
            ! Independent dense unit-lower forward elimination in pivot order.
            do j = 1, n
                do i = j + 1, n
                    reference(:, p(i)) = reference(:, p(i)) &
                        - dense(p(i), p(j))*reference(:, p(j))
                end do
            end do
            call lu_panel_prefix(f, marks, coeff, 3, reach, count, inv, lp, &
                rows, values)
            if (any(marks /= expected)) &
                error stop 'production panel lane reach vs Boolean closure'
            if (maxval(abs(f - reference)) > 3e-15_dp) &
                error stop 'production prefix vs dense oracle'
        end do
    end subroutine closure_oracle

    subroutine dense_csc(a, ptr, rows, values)
        real(dp), intent(in) :: a(:, :)
        integer, allocatable, intent(out) :: ptr(:), rows(:)
        real(dp), allocatable, intent(out) :: values(:)
        integer :: i, j, k, n
        n = size(a, 1)
        allocate (ptr(n + 1), rows(count(a /= 0)), values(count(a /= 0)))
        k = 1
        do j = 1, n
            ptr(j) = k
            do i = 1, n
                if (a(i, j) == 0) cycle
                rows(k) = i; values(k) = a(i, j); k = k + 1
            end do
        end do
        ptr(n + 1) = k
    end subroutine dense_csc

    subroutine factor_oracles()
        integer, parameter :: sizes(*) = [1, 7, 8, 9, 15, 16, 17, 25]
        type(lu_factor_t) :: scalar, panel, plus, minus
        type(fortsparse_status_t) :: s
        real(dp), allocatable :: a(:, :), x(:, :), b(:, :), b0(:, :)
        real(dp), allocatable :: d(:, :), da(:, :), db(:, :)
        real(dp), allocatable :: values(:), vp(:), vm(:), bp(:, :), bm(:, :)
        integer, allocatable :: ptr(:), rows(:), q(:), pref(:)
        integer :: n, caseid, variant, i, j, k, mode, pass
        real(dp), parameter :: h = 1e-5_dp
        do caseid = 1, size(sizes)
            n = sizes(caseid)
            allocate (a(n, n), x(n, 3), b(n, 3), b0(n, 3), d(n, 3), &
                da(n, n), db(n, 3), q(n), pref(n))
            do i = 1, n
                q(i) = n + 1 - i
                pref(i) = mod(i + 3, n) + 1
                do j = 1, 3
                    x(i, j) = real(i - 2*j, dp)/11
                    d(i, j) = real(2*i + j, dp)/23
                end do
            end do
            do variant = 1, 3
                a = 0
                do j = 1, n
                    ! Variant1 dense, variant2 sparse, variant3 disconnected blocks.
                    do i = 1, n
                        if (variant == 2 .and. &
                            mod(3*i + 2*j, 7) /= 0 .and. i /= j) cycle
                        if (variant == 3 .and. (i - 1)/3 /= (j - 1)/3) cycle
                        a(pref(i), j) = real(mod(3*i + 5*j, 9) - 4, dp)/32
                    end do
                    a(pref(j), j) = real(n + 3, dp)
                end do
                call dense_csc(a, ptr, rows, values)
                do pass = 0, 1
                    call scalar%factor_raw(n, ptr, rows, values, s, q, 0.1_dp, &
                                           pruning=pass == 1, &
                                           row_preference=pref)
                    call checked(s)
                    call panel%factor_raw(n, ptr, rows, values, s, q, 0.1_dp, &
                                          pruning=pass == 1, &
                                          row_preference=pref, panel_width=8)
                    call checked(s)
                    if (panel%panel_width /= 8) error stop 'panel width metadata'
                    do mode = 0, 2
                        b0 = matmul(a, x)
                        if (mode /= 0) b0 = matmul(transpose(a), x)
                        b = b0
                        call panel%solve(b, s, mode); call checked(s)
                        if (maxval(abs(b - x)) > 2e-12_dp) &
                            error stop 'panel known-x N/T/H block'
                        b = b0
                        call scalar%solve(b, s, mode); call checked(s)
                        if (maxval(abs(b - x)) > 2e-12_dp) &
                            error stop 'reference known-x'
                    end do
                    ! Same-pattern numeric refactor changes every value.
                    call panel%factor_raw(n, ptr, rows, 1.125_dp*values, s, q, 0.1_dp, &
                                          pruning=pass == 1, row_preference=pref, &
                                          panel_width=8)
                    call checked(s)
                    if (.not. panel%last_reused_pattern) &
                        error stop 'panel refactor reuse'
                    b = 1.125_dp*matmul(a, x)
                    call panel%solve(b, s); call checked(s)
                    if (maxval(abs(b - x)) > 2e-12_dp) error stop 'panel changed values'
                end do
                if (variant == 1) then
                    da = 0.03125_dp*a
                    db = matmul(da, x) + matmul(a, d)
                    b = db - matmul(da, x)
                    call panel%factor_raw(n, ptr, rows, values, s, q, pruning=.true., &
                                          panel_width=8)
                    call checked(s)
                    call panel%solve(b, s); call checked(s)
                    if (maxval(abs(b - d)) > 2e-12_dp) error stop 'panel implicit JVP'
                    b0 = d
                    call panel%solve(b0, s, 1); call checked(s)
                    if (abs(sum(d*d) - sum(b0*(db - matmul(da, x)))) > 2e-10_dp) &
                        error stop 'panel implicit VJP duality'
                    vp = values*(1 + h*0.03125_dp); vm = values*(1 - h*0.03125_dp)
                    bp = matmul(a, x) + h*db; bm = matmul(a, x) - h*db
                    call plus%factor_raw(n, ptr, rows, vp, s, q, &
                                         panel_width=8); call checked(s)
                    call minus%factor_raw(n, ptr, rows, vm, s, q, &
                                          panel_width=8); call checked(s)
                    call plus%solve(bp, s); call checked(s)
                    call minus%solve(bm, s); call checked(s)
                    if (maxval(abs((bp - bm)/(2*h) - d)) > 2e-8_dp) &
                        error stop 'panel finite-difference derivative'
                end if
                do j = 1, n
                    do k = ptr(j), ptr(j + 1) - 1
                        if (values(k) /= a(rows(k), j)) &
                            error stop 'panel changed borrowed values'
                    end do
                    if (q(j) /= n + 1 - j .or. pref(j) /= mod(j + 3, n) + 1) &
                        error stop 'panel changed borrowed permutation'
                end do
                deallocate (ptr, rows, values)
            end do
            call panel%free(); call scalar%free(); call plus%free(); call minus%free()
            deallocate (a, x, b, b0, d, da, db, q, pref)
        end do
        ! Exact cancellation: dependent-looking prefix but nonsingular remainder.
        block
            real(dp) :: a(9, 9), rhs(9, 1)
            integer :: j
            a = 0
            do j = 1, 9
                a(j, j) = 2
            end do
            a(1:3, 1) = [1.0_dp, 1.0_dp, 1.0_dp]
            a(1:3, 2) = [1.0_dp, 1.0_dp, 0.0_dp]
            a(1:3, 3) = [0.0_dp, 1.0_dp, 1.0_dp]
            call dense_csc(a, ptr, rows, values)
            call panel%factor_raw(9, ptr, rows, values, s, pruning=.true., &
                panel_width=8)
            call checked(s)
            rhs = matmul(a, reshape([(real(j, dp), j=1, 9)], [9, 1]))
            call panel%solve(rhs, s); call checked(s)
            if (maxval(abs(rhs(:, 1) - [(real(j, dp), j=1, 9)])) > 1e-13_dp) &
                error stop 'panel exact cancellation'
        end block
    end subroutine factor_oracles

    subroutine memory_and_failure_oracles()
        type(lu_factor_t) :: scalar, panel
        type(lu_factor_z_t) :: z
        type(csc_t) :: a
        type(fortsparse_status_t) :: s
        integer :: ptr(18), rows(17), i
        real(dp) :: values(17), rhs(17)
        integer(i8) :: delta, peak
        do i = 1, 17
            ptr(i) = i; rows(i) = i; values(i) = 2
        end do
        ptr(18) = 18
        call scalar%factor_raw(17, ptr, rows, values, s); call checked(s)
        call panel%factor_raw(17, ptr, rows, values, s, panel_width=8); call checked(s)
        delta = 68_i8*17 + 64
        if (panel%memory_bytes() - scalar%memory_bytes() /= delta) &
            error stop 'panel exact payload delta'
        peak = panel%peak_memory_bytes
        call panel%free()
        call panel%factor_raw(17, ptr, rows, values, s, memory_limit_bytes=peak, &
                              panel_width=8)
        call checked(s)
        if (panel%peak_memory_bytes > peak) error stop 'panel exceeded cap'
        call panel%free()
        call panel%factor_raw(17, ptr, rows, values, s, memory_limit_bytes=1_i8, &
                              panel_width=8)
        if (s%code /= FORTSPARSE_RESOURCE_LIMIT) error stop 'tiny panel cap'
        if (panel%memory_bytes() > 1_i8) error stop 'tiny cap allocated payload'
        call panel%factor_raw(17, ptr, rows, values, s, panel_width=8, &
            low_memory=.true.)
        if (s%code /= FORTSPARSE_BACKEND_UNAVAILABLE) error stop 'panel LOW combination'
        call panel%factor_raw(17, ptr, rows, values, s, panel_width=3)
        if (s%code /= FORTSPARSE_INVALID_MATRIX) error stop 'invalid panel width'
        call z%factor_raw(17, ptr, rows, cmplx(values, 0.0_dp, dp), s, panel_width=8)
        if (s%code /= FORTSPARSE_BACKEND_UNAVAILABLE) &
            error stop 'complex panel unavailable'
        ! The container wrapper forwards the same final optional argument.
        a%nrow = 17; a%ncol = 17; a%nnz = 17
        a%col_ptr = ptr; a%row_idx = rows; a%val = values
        call panel%factor(a, s, panel_width=8); call checked(s)
        if (any(a%val /= values) .or. &
            any(a%row_idx /= rows) .or. any(a%col_ptr /= ptr)) &
            error stop 'panel mutated borrowed CSC'
        values(2) = 0
        call panel%factor_raw(17, ptr, rows, values, s, panel_width=8)
        if (s%code /= FORTSPARSE_SINGULAR) error stop 'singular panel'
        rhs = 1
        call panel%solve(rhs, s)
        if (s%code /= FORTSPARSE_NOT_FACTORED) error stop 'failed panel ready'
        values = 2
        call panel%factor_raw(17, ptr, rows, values, s, panel_width=8); call checked(s)
        rhs = values
        call panel%solve(rhs, s); call checked(s)
        if (any(rhs /= 1)) error stop 'failed panel recovery'
        values(3) = ieee_value(1.0_dp, ieee_quiet_nan)
        call panel%factor_raw(17, ptr, rows, values, s, panel_width=8)
        if (s%code /= FORTSPARSE_INVALID_MATRIX) error stop 'panel nonfinite guard'
        call panel%free(); call panel%free()
        block
            real(dp) :: empty(0)
            integer :: empty_rows(0)
            call scalar%free()
            call scalar%factor_raw(0, [1], empty_rows, empty, s)
            call checked(s)
            call panel%factor_raw(0, [1], empty_rows, empty, s, panel_width=8)
            call checked(s)
            if (panel%memory_bytes() - scalar%memory_bytes() /= 64_i8) &
                error stop 'zero-dimensional panel payload'
            call panel%solve(empty, s); call checked(s)
            call panel%free()
        end block
        if (panel%memory_bytes() /= 0) error stop 'panel idempotent free'
    end subroutine memory_and_failure_oracles
end program test_fortsparse_lu_panel
