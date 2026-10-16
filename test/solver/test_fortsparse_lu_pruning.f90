program test_fortsparse_lu_pruning
    use fortsparse
    use fortsparse_native_lu, only: lu_factor_t, prune_real => lu_prune_column
    use fortsparse_native_lu_z, only: lu_factor_z_t, prune_complex => lu_prune_column
    implicit none
    call graph_oracle()
    call numerical_oracles()
    call memory_oracle()
    print '(a)', &
        'Witnessed LU graph closure, cancellation, pivots, derivatives and caps pass'
contains
    subroutine checked(status)
        type(fortsparse_status_t), intent(in) :: status
        if (.not. status_ok(status)) then
            print *, status%code, trim(status%msg)
            error stop 'unexpected pruning status'
        end if
    end subroutine checked

    function closure(graph) result(reach)
        logical, intent(in) :: graph(:, :)
        logical :: reach(size(graph, 1), size(graph, 2))
        integer :: i, j, k
        reach = graph
        do i = 1, size(graph, 1)
            reach(i, i) = .true.
        end do
        ! Independent Boolean Floyd-Warshall; no DFS or pruning rule here.
        do k = 1, size(graph, 1)
            do j = 1, size(graph, 1)
                do i = 1, size(graph, 1)
                    reach(i, j) = reach(i, j) .or. (reach(i, k) .and. reach(k, j))
                end do
            end do
        end do
    end function closure

    subroutine graph_oracle()
        integer, parameter :: n = 5, m = n*(n - 1)/2, p(n) = [3, 5, 1, 4, 2]
        integer :: lp(n + 1), rows(m), zrows(m), ends(n), zends(n), marks(n)
        real(dp) :: values(m), tag
        complex(dp) :: zvalues(m)
        logical :: full(n, n), reduced(n, n), actual(n, n), expected(n, n)
        integer :: mask, bit, i, j, k, e, next, removed
        removed = 0
        ! Exhaust all 1024 DAGs for a nontrivial ordering of original row IDs.
        do mask = 0, 2**m - 1
            next = 1; bit = 0; full = .false.; rows = 0; values = 0.0_dp
            do j = 1, n
                lp(j) = next
                do i = j + 1, n
                    if (btest(mask, bit)) then
                        rows(next) = p(i)
                        values(next) = real(16*p(j) + p(i), dp)
                        next = next + 1
                        full(p(j), p(i)) = .true.
                    end if
                    bit = bit + 1
                end do
            end do
            lp(n + 1) = next
            zrows = rows
            zvalues = cmplx(values, -values, dp)
            ends = lp(2:); zends = ends; marks = 0
            do j = 1, n
                do e = lp(j), lp(j + 1) - 1
                    marks(rows(e)) = -j
                end do
                do k = 1, j - 1
                    call prune_real(rows, values, lp(k), lp(k + 1), p(j), marks, j, &
                                    ends(k))
                    call prune_complex(zrows, zvalues, lp(k), lp(k + 1), p(j), &
                                       marks, j, zends(k))
                    if (any(ends /= zends) &
                        .or. any(rows(:next - 1) /= zrows(:next - 1))) &
                        error stop 'real/complex graph partition differs'
                    reduced = .false.
                    do i = 1, j
                        do e = lp(i), abs(ends(i)) - 1
                            reduced(p(i), rows(e)) = .true.
                        end do
                    end do
                    actual = closure(reduced)
                    expected = full
                    do i = j + 1, n
                        expected(p(i), :) = .false.
                    end do
                    expected = closure(expected)
                    if (any(actual .neqv. expected)) &
                        error stop 'pruning changed reach closure'
                    do i = 1, n
                        do e = lp(i), lp(i + 1) - 1
                            tag = real(16*p(i) + rows(e), dp)
                            if (values(e) /= tag .or. zvalues(e) /= cmplx(tag, &
                                                                          -tag, dp)) &
                                error stop 'pruning detached a numerical coefficient'
                        end do
                    end do
                end do
            end do
            removed = removed + sum(lp(2:) - abs(ends))
        end do
        if (removed == 0) error stop 'graph oracle exercised no removed edges'
    end subroutine graph_oracle

    subroutine from_dense(matrix, a, z)
        real(dp), intent(in) :: matrix(:, :)
        type(csc_t), intent(out) :: a
        type(csc_z_t), intent(out) :: z
        integer :: i, j, k, n, nz
        n = size(matrix, 1); nz = count(matrix /= 0.0_dp)
        a%nrow = n; a%ncol = n; a%nnz = nz
        allocate (a%col_ptr(n + 1), a%row_idx(nz), a%val(nz))
        k = 1
        do j = 1, n
            a%col_ptr(j) = k
            do i = 1, n
                if (matrix(i, j) == 0.0_dp) cycle
                a%row_idx(k) = i; a%val(k) = matrix(i, j); k = k + 1
            end do
        end do
        a%col_ptr(n + 1) = k
        z%nrow = n; z%ncol = n; z%nnz = nz
        z%col_ptr = a%col_ptr; z%row_idx = a%row_idx
        z%val = cmplx(a%val, a%val, dp)
    end subroutine from_dense

    subroutine numerical_oracles()
        integer, parameter :: n = 4, q(n) = [3, 1, 4, 2]
        type(csc_t) :: a, adot, bad
        type(csc_z_t) :: z, zbad, zdot
        type(lu_factor_t) :: candidate, reference
        type(lu_factor_z_t) :: zcandidate, zreference
        type(sparse_solver_t) :: solver
        type(sparse_factor_stats_t) :: production_stats,reference_stats
        type(fortsparse_status_t) :: status
        real(dp) :: matrix(n, n), original(n, n), truth(n, 3), rhs(n, 3), answer(n, 3)
        real(dp) :: badmatrix(n, n)
        real(dp) :: b(n), x(n), dx(n), zero(n), seed(n), bbar(n), values_bar(n*n)
        complex(dp) :: zmatrix(n, n), ztruth(n, 3), zrhs(n, 3), zanswer(n, 3)
        complex(dp) :: zx(n),zdx(n),zseed(n),zbbar(n),zvalues_bar(n*n),zzero(n)
        integer :: i, j, op, repeat, fixture
        ! L32 cancels exactly; deleting1->3 using only a symmetric pair is wrong.
        original = reshape([2._dp, 1._dp, 1._dp, 0._dp, &
                            1._dp, 2._dp, .5_dp, 0._dp, &
                            1._dp, .5_dp, 2._dp, 0._dp, &
                            1._dp, .5_dp, 0._dp, 2._dp], [n, n])
        do fixture = 1, 2
            matrix = original
            if (fixture == 2) then
                matrix = 1.0_dp
                do i = 1, n
                    matrix(i, i) = 4.0_dp
                end do
            end if
            do j = 1, 3
                do i = 1, n
                    truth(i, j) = real(i - 2*j, dp)/8
                end do
            end do
            ztruth = cmplx(truth, -.5_dp*truth, dp)
            do repeat = 1, 2
                if (repeat == 2) matrix = matrix([3, 1, 4, 2], :)
                call from_dense(matrix, a, z)
                zmatrix = cmplx(matrix, matrix, dp)
                if (repeat == 1) then
                    call candidate%factor(a, status, pruning=.true.)
                    call checked(status)
                    call reference%factor(a, status, pruning=.false.)
                    call checked(status)
                    call zcandidate%factor(z, status, pruning=.true.)
                    call checked(status)
                    call zreference%factor(z, status, pruning=.false.)
                else
                    call candidate%factor(a, status, q, pruning=.true.)
                    call checked(status)
                    call reference%factor(a, status, q, pruning=.false.)
                    call checked(status)
                    call zcandidate%factor(z, status, q, pruning=.true.)
                    call checked(status)
                    call zreference%factor(z, status, q, pruning=.false.)
                end if
                call checked(status)
                do op = 0, 2
                    if (op == 0) then
                        rhs = matmul(matrix, truth)
                        zrhs = matmul(zmatrix, ztruth)
                    else
                        rhs = matmul(transpose(matrix), truth)
                        if (op == 1) zrhs = matmul(transpose(zmatrix), ztruth)
                        if (op == 2) zrhs = matmul(conjg(transpose(zmatrix)), ztruth)
                    end if
                    answer = rhs; zanswer = zrhs
                    call candidate%solve(answer, status, op)
                    call checked(status)
                    call zcandidate%solve(zanswer, status, op)
                    call checked(status)
                    if (maxval(abs(answer - truth)) > 1e-13_dp .or. &
                        maxval(abs(zanswer - ztruth)) > 1e-13_dp) &
                        error stop 'pruned original-matrix N/T/H block oracle'
                    call reference%solve(rhs, status, op)
                    call checked(status)
                    call zreference%solve(zrhs, status, op)
                    call checked(status)
                    if (maxval(abs(rhs - answer)) > 1e-13_dp .or. &
                        maxval(abs(zrhs - zanswer)) > 1e-13_dp) &
                        error stop 'pruned/full-graph factors disagree'
                end do
            end do
            solver%backend_id = FORTSPARSE_BACKEND_NATIVE
            solver%refine = .true.
            call sparse_factor(solver, a, status)
            call checked(status)
            b = matmul(matrix, truth(:, 1)); zero = 0; seed = truth(:, 2)
            call sparse_solve(solver, b, x, status)
            call checked(status)
            adot = a; adot%val = .125_dp*a%val
            call sparse_solve_jvp(solver, adot, x, zero, dx, status)
            call checked(status)
            if (maxval(abs(dx + .125_dp*truth(:, 1))) > 1e-13_dp) &
                error stop 'pruned implicit JVP oracle'
            call sparse_solve_vjp(solver, a, x, seed, bbar, values_bar(:a%nnz), &
                                  status, use_primal_factors=.true.)
            call checked(status)
            if (abs(dot_product(seed, dx) - dot_product(values_bar(:a%nnz), adot%val)) &
                > 1e-13_dp) error stop 'pruned implicit VJP pairing'
            call sparse_factor_stats(solver,production_stats,status)
            call checked(status)
            solver%backend_id=FORTSPARSE_BACKEND_NATIVE_REFERENCE
            call sparse_factor(solver,a,status)
            call checked(status)
            call sparse_factor_stats(solver,reference_stats,status)
            call checked(status)
            if(production_stats%host_current_bytes-reference_stats%host_current_bytes &
                /=int(storage_size(0)/8,i8)*n) error stop 'public reference kept prune metadata'
            solver%backend_id=FORTSPARSE_BACKEND_NATIVE
            call sparse_factor(solver,z,status)
            call checked(status)
            zrhs=matmul(zmatrix,ztruth)
            call sparse_solve(solver,zrhs(:,1),zx,status)
            call checked(status)
            zdot=z;zdot%val=.125_dp*z%val
            zzero=cmplx(0.0_dp,0.0_dp,dp);zseed=ztruth(:,2)
            call sparse_solve_jvp(solver,zdot,zx,zzero,zdx,status)
            call checked(status)
            if(maxval(abs(zdx+.125_dp*ztruth(:,1)))>1e-13_dp) &
                error stop 'complex pruned implicit JVP'
            call sparse_solve_vjp(solver,z,zx,zseed,zbbar,zvalues_bar(:z%nnz), &
                status,use_primal_factors=.true.)
            call checked(status)
            if(abs(real(dot_product(zseed,zdx)- &
                dot_product(zvalues_bar(:z%nnz),zdot%val),dp))>1e-13_dp) &
                error stop 'complex pruned implicit VJP'
            call sparse_factor_stats(solver,production_stats,status)
            call checked(status)
            solver%backend_id=FORTSPARSE_BACKEND_NATIVE_REFERENCE
            call sparse_factor(solver,z,status)
            call checked(status)
            call sparse_factor_stats(solver,reference_stats,status)
            call checked(status)
            if(production_stats%host_current_bytes-reference_stats%host_current_bytes &
                /=int(storage_size(0)/8,i8)*n) error stop 'complex reference prune metadata'
            call sparse_destroy(solver)
        end do
        ! Fail only after earlier columns have performed witnessed pruning.
        ! All entries remain present, so the retry reuses the same allocation.
        badmatrix = matrix
        badmatrix(n, :) = badmatrix(n - 1, :)
        call from_dense(badmatrix, bad, zbad)
        call candidate%factor(bad, status, pruning=.true.)
        if (status%code /= FORTSPARSE_SINGULAR) error stop 'pruned late singular status'
        call zcandidate%factor(zbad, status, pruning=.true.)
        if (status%code /= FORTSPARSE_SINGULAR) &
            error stop 'complex late singular status'
        call candidate%factor(a, status, pruning=.true.)
        call checked(status)
        call zcandidate%factor(z, status, pruning=.true.)
        call checked(status)
        if (.not. candidate%last_reused_pattern &
            .or. .not. zcandidate%last_reused_pattern) &
            error stop 'pruned failed-pattern reuse'
        b = matmul(matrix, truth(:, 1))
        zrhs = matmul(zmatrix, ztruth)
        call candidate%solve(b, status)
        call checked(status)
        call zcandidate%solve(zrhs, status)
        call checked(status)
        if (maxval(abs(b - truth(:, 1))) > 1e-13_dp .or. &
            maxval(abs(zrhs - ztruth)) > 1e-13_dp) &
                error stop 'pruning failure/reuse reset'
    end subroutine numerical_oracles

    subroutine memory_oracle()
        type(lu_factor_t) :: factor
        type(fortsparse_status_t) :: status
        integer(i8) :: base, extra, current, ib, vb
        real(dp) :: values(3) = [2._dp, 3._dp, 4._dp]
        ib = int(storage_size(0)/8, i8)
        vb = int(storage_size(0.0_dp)/8, i8)
        base = ib*(10_i8*3 + 3 + 3) + 2_i8*vb*3
        extra = ib*3
        call factor%factor_raw(3, [1, 2, 3, 4], [1, 2, 3], values, status, &
                               memory_limit_bytes=base + extra - 1_i8, pruning=.true.)
        if (status%code /= FORTSPARSE_RESOURCE_LIMIT) &
            error stop 'prune fixed cap boundary'
        if (factor%memory_bytes() /= 0 &
            .or. factor%failed_required_bytes /= base + extra) &
            error stop 'prune fixed quota includes metadata'
        call factor%factor_raw(3, [1, 2, 3, 4], [1, 2, 3], values, status, &
                               pruning=.false.)
        call checked(status)
        current = factor%memory_bytes()
        call factor%factor_raw(3, [1, 2, 3, 4], [1, 2, 3], values, status, &
                               memory_limit_bytes=current + extra - 1_i8, &
                               pruning=.true.)
        if (status%code /= FORTSPARSE_RESOURCE_LIMIT) &
            error stop 'prune reuse cap boundary'
        if (factor%memory_bytes() /= current &
            .or. factor%failed_required_bytes /= extra) &
            error stop 'prune retained-factor overlap quota'
        call factor%factor_raw(3, [1, 2, 3, 4], [1, 2, 3], values, status, &
                               memory_limit_bytes=current + extra, pruning=.true.)
        call checked(status)
        if (factor%memory_bytes() /= current + extra .or. &
            factor%peak_memory_bytes /= current + extra) &
                error stop 'prune exact reuse cap'
        call factor%factor_raw(3, [1, 2, 3, 4], [1, 2, 3], values, status, &
                               pruning=.false.)
        call checked(status)
        if (factor%memory_bytes() /= current) &
            error stop 'unpruned path retained metadata'
        call factor%free()
    end subroutine memory_oracle
end program test_fortsparse_lu_pruning
