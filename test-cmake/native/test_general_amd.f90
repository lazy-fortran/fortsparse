program test_general_amd
    use fortsparse
    use fortsparse_native_general, only: native_general_backend_t
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    integer, parameter :: n = 9
    type(sparse_solver_t) :: solver
    type(csc_t) :: A, Adot, changed, bad
    type(csc_z_t) :: Z, Zdot
    type(fortsparse_status_t) :: status
    type(sparse_factor_stats_t) :: stats
    real(dp) :: dense(n, n), exact(n), rhs(n), x(n), dx(n), zeros(n), seed(n), bbar(n)
    real(dp), allocatable :: abar(:)
    complex(dp) :: zdense(n, n), zexact(n), zrhs(n), zx(n), zdx(n), zzeros(n), &
                   zseed(n), zbbar(n)
    complex(dp), allocatable :: zbar(:)
    integer :: i, j, p, op, before
    logical :: reused, available
    interface
        integer(c_int) function amd_calls() bind(C, &
                                                 name='fsparse_test_general_amd_calls')
            import c_int
        end function
        integer(c_int) function amd_nonidentity() &
            bind(C, name='fsparse_test_general_amd_nonidentity')
            import c_int
        end function
        integer(c_int) function amd_available() &
            bind(C, name='fsparse_test_general_amd_available')
            import c_int
        end function
    end interface
    available = amd_available() /= 0
    dense = 0.0_dp
    do i = 1, n
        dense(i, i) = 4.0_dp + .1_dp*i
        dense(i, 1) = .2_dp*i
        exact(i) = sin(real(i, dp))
        seed(i) = cos(real(i, dp))
    end do
    dense(1, 1) = 0.0_dp ! Requires numerical row pivoting.
    dense(1, 2) = 1.0_dp
    dense(2, 1) = 2.0_dp
    call matrix_from_dense(dense, A)
    call require(A%nnz < n*n, 'genuinely asymmetric sparse fixture')
    zeros = 0.0_dp
    solver%backend_id = FORTSPARSE_BACKEND_NATIVE
    solver%options%ordering = FORTSPARSE_ORDER_AMD
    before = amd_calls()
    solver%options%host_memory_limit_bytes = 1024_i8*1024_i8
    call sparse_factor(solver, A, status)
    call require(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
                 'opaque AMD cap rejected')
    call require(amd_calls() == before, 'cap rejects before AMD')
    solver%options%host_memory_limit_bytes = 0_i8
    solver%options%memory_policy = FORTSPARSE_MEMORY_LOW
    call sparse_factor(solver, A, status)
    call require(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
                 'LOW_MEMORY AMD rejected')
    call require(amd_calls() == before, 'LOW_MEMORY rejects before AMD')
    solver%options%memory_policy = FORTSPARSE_MEMORY_AUTO
    solver%refine = .true.
    call sparse_factor(solver, A, status)
    if (.not. available) then
        call require(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
                     'AMD-off explicit status')
    else
        call checked(status)
        call require(amd_calls() == before + 1, 'initial AMD call')
        call require(amd_nonidentity() /= 0, 'ordering actually permutes columns')
        call real_operations()
        call sparse_factor_stats(solver, stats, status)
        call checked(status)
        call require(stats%ordering_used == FORTSPARSE_ORDER_AMD, &
                     'actual AMD statistic')
        call require(stats%host_current_bytes > 0, &
                     'exact retained payload remains reported')
        call require(stats%host_peak_bytes == -1 .and. stats%factor_mem_mb < 0, &
                     'opaque complete peak explicitly unknown')
        call require(.not. stats%strict_memory_supported, &
                     'no false strict AMD guarantee')
        before = amd_calls()
        A%val = 1.125_dp*A%val
        dense = 1.125_dp*dense
        call sparse_refactor(solver, A, status, reused)
        call checked(status)
        call require(reused .and. amd_calls() == before, 'value-only refactor reuses q')
        call real_operations()
        call sparse_factor_stats(solver, stats, status)
        call checked(status)
        call require(stats%host_peak_bytes == -1, &
                     'AMD refactor complete peak stays unknown')
        allocate (abar(A%nnz))
        Adot = A
        Adot%val = .125_dp*A%val
        rhs = matmul(dense, exact)
        call sparse_solve(solver, rhs, x, status)
        call checked(status)
        call sparse_solve_jvp(solver, Adot, x, zeros, dx, status)
        call checked(status)
        call require(maxval(abs(dx + .125_dp*exact)) < 2.e-12_dp, 'AMD real JVP oracle')
        call sparse_solve_vjp(solver, A, x, seed, bbar, abar, status, &
                              use_primal_factors=.true.)
        call checked(status)
        call require(abs(dot_product(seed, dx) - dot_product(abar, Adot%val)) &
                     < 2.e-12_dp, 'AMD real adjoint pairing')
        call require(amd_calls() == before, 'derivatives retain AMD factors')
        ! Change a pivot candidate without changing the graph or column order.
        do p = A%col_ptr(1), A%col_ptr(2) - 1
            if (A%row_idx(p) == 3) A%val(p) = 8.0_dp
        end do
        dense(3, 1) = 8.0_dp
        call sparse_refactor(solver, A, status, reused)
        call checked(status)
        call require(reused .and. amd_calls() == before, 'numeric pivots retain q')
        call real_operations()
        dense(n, 2) = .03125_dp
        call matrix_from_dense(dense, changed)
        call sparse_refactor(solver, changed, status, reused)
        call checked(status)
        call require(.not. reused .and. amd_calls() == before + 1, &
                     'changed pattern recomputes q')
        A = changed
        call real_operations()
        solver%options%ordering = FORTSPARSE_ORDER_DEFAULT
        call sparse_refactor(solver, A, status, reused)
        call checked(status)
        call require(.not. reused, 'changed ordering strategy invalidates reuse')
        call sparse_factor_stats(solver, stats, status)
        call require(stats%ordering_used == FORTSPARSE_ORDER_DEFAULT, &
                     'natural order restored')
        call require(stats%host_peak_bytes >= 0 .and. stats%strict_memory_supported, &
                     'natural payload peak and cap contract preserved')
    end if
    solver%options%ordering = FORTSPARSE_ORDER_DEFAULT
    call sparse_factor(solver, A, status)
    call checked(status)
    call real_operations()
    solver%options%ordering = FORTSPARSE_ORDER_AMD
    solver%backend_id = FORTSPARSE_BACKEND_NATIVE_REFERENCE
    before = amd_calls()
    call sparse_factor(solver, A, status)
    call require(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
                 'retained handle switches to reference')
    call require(amd_calls() == before, 'reference rejects before AMD')
    solver%backend_id = FORTSPARSE_BACKEND_NATIVE
    solver%options%ordering = FORTSPARSE_ORDER_DEFAULT
    call sparse_factor(solver, A, status)
    call checked(status)
    call sparse_destroy(solver)
    solver%backend_id = FORTSPARSE_BACKEND_NATIVE
    solver%options%ordering = FORTSPARSE_ORDER_AMD
    bad = A
    bad%col_ptr(2) = bad%nnz + 2
    call sparse_factor(solver, bad, status)
    if (available) call require(status%code == FORTSPARSE_INVALID_MATRIX, &
                                'malformed CSC guard')
    if (available) then
        Z%nrow = n
        Z%ncol = n
        Z%nnz = A%nnz
        Z%col_ptr = A%col_ptr
        Z%row_idx = A%row_idx
        allocate (Z%val(Z%nnz))
        zdense = cmplx(0.0_dp, 0.0_dp, dp)
        do j = 1, n
            do p = A%col_ptr(j), A%col_ptr(j + 1) - 1
                i = A%row_idx(p)
                Z%val(p) = cmplx(A%val(p), .03125_dp*(i - j), dp)
                zdense(i, j) = Z%val(p)
            end do
        end do
        zexact = cmplx(exact, seed, dp)
        zseed = cmplx(seed, -exact, dp)
        zzeros = cmplx(0.0_dp, 0.0_dp, dp)
        before = amd_calls()
        call sparse_factor(solver, Z, status)
        call checked(status)
        call require(amd_calls() == before + 1, 'complex AMD initial analysis')
        call complex_operations()
        before = amd_calls()
        Z%val = cmplx(1.125_dp, .0625_dp, dp)*Z%val
        zdense = cmplx(1.125_dp, .0625_dp, dp)*zdense
        call sparse_refactor(solver, Z, status, reused)
        call checked(status)
        call require(reused .and. amd_calls() == before, 'complex q retained')
        call complex_operations()
        allocate (zbar(Z%nnz))
        Zdot = Z
        Zdot%val = .125_dp*Z%val
        zrhs = matmul(zdense, zexact)
        call sparse_solve(solver, zrhs, zx, status)
        call checked(status)
        call sparse_solve_jvp(solver, Zdot, zx, zzeros, zdx, status)
        call checked(status)
        call require(maxval(abs(zdx + .125_dp*zexact)) < 2.e-12_dp, &
                     'AMD complex JVP oracle')
        call sparse_solve_vjp(solver, Z, zx, zseed, zbbar, zbar, status, &
                              use_primal_factors=.true.)
        call checked(status)
        call require(abs(real(dot_product(zseed, zdx) - dot_product(zbar, &
                                                                    Zdot%val), dp)) &
                     < 2.e-12_dp, 'AMD complex adjoint pairing')
        call require(amd_calls() == before, 'complex derivatives retain factors')
    end if
    call sparse_destroy(solver)
    call raw_operations()
    print '(a)', &
        'General AMD ordering, N/T/H, reuse, derivatives and cap contracts pass'
contains
    subroutine raw_operations()
        type(native_general_backend_t) :: raw
        real(dp) :: values(A%nnz)
        complex(dp) :: zvalues(A%nnz)
        integer :: calls
        raw%options%ordering = FORTSPARSE_ORDER_AMD
        calls = amd_calls()
        call raw%factor_real_raw(n, n, A%nnz, A%col_ptr, A%row_idx, A%val, .true., &
                                 status)
        if (.not. available) then
            call require(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
                         'raw AMD-off status')
            call raw%free()
            return
        end if
        call checked(status)
        call require(amd_calls() == calls + 1, 'raw real initial ordering')
        values = 1.25_dp*A%val
        call raw%refactor_real_raw(n, n, A%nnz, A%col_ptr, A%row_idx, values, &
                                   .true., status, reused)
        call checked(status)
        call require(reused .and. amd_calls() == calls + 1, 'raw real q reuse')
        do op = 0, 2
            if (op == 0) then
                rhs = 1.25_dp*matmul(dense, exact)
            else
                rhs = 1.25_dp*matmul(transpose(dense), exact)
            end if
            call raw%solve_real_op(rhs, x, op, status)
            call checked(status)
            call require(maxval(abs(x - exact)) < 2.e-12_dp, 'raw real N/T/H oracle')
        end do
        call raw%factor_real_raw(n, n, A%nnz, A%col_ptr(:n), A%row_idx, A%val, &
                                 .false., status)
        call require(status%code == FORTSPARSE_INVALID_MATRIX, &
                     'raw pointer extent guard')
        calls = amd_calls()
        call raw%factor_complex_raw(n, n, Z%nnz, Z%col_ptr, Z%row_idx, Z%val, &
                                    .true., status)
        call checked(status)
        call require(amd_calls() == calls + 1, 'raw complex initial ordering')
        zvalues = cmplx(.875_dp, .03125_dp, dp)*Z%val
        call raw%refactor_complex_raw(n, n, Z%nnz, Z%col_ptr, Z%row_idx, zvalues, &
                                      .true., status, reused)
        call checked(status)
        call require(reused .and. amd_calls() == calls + 1, 'raw complex q reuse')
        do op = 0, 2
            select case (op)
            case (0)
                zrhs = cmplx(.875_dp, .03125_dp, dp)*matmul(zdense, zexact)
            case (1)
                zrhs = cmplx(.875_dp, .03125_dp, dp)* &
                       matmul(transpose(zdense), zexact)
            case (2)
                zrhs = cmplx(.875_dp, -.03125_dp, dp)* &
                       matmul(conjg(transpose(zdense)), zexact)
            end select
            call raw%solve_complex_op(zrhs, zx, op, status)
            call checked(status)
            call require(maxval(abs(zx - zexact)) < 2.e-12_dp, &
                         'raw complex N/T/H oracle')
        end do
        call raw%free()
    end subroutine raw_operations
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print '(a)', message
            error stop 1
        end if
    end subroutine
    subroutine checked(s)
        type(fortsparse_status_t), intent(in) :: s
        if (.not. status_ok(s)) then
            print '(i0,2a)', s%code, ': ', s%msg
            error stop 1
        end if
    end subroutine
    subroutine matrix_from_dense(values, result)
        real(dp), intent(in) :: values(:, :)
        type(csc_t), intent(out) :: result
        integer :: col, row, k, nz
        nz = count(abs(values) > 0.0_dp) + 1 ! Stored zero first diagonal.
        result%nrow = n
        result%ncol = n
        result%nnz = nz
        allocate (result%col_ptr(n + 1), result%row_idx(nz), result%val(nz))
        k = 1
        do col = 1, n
            result%col_ptr(col) = k
            do row = 1, n
                if (values(row, col) == 0.0_dp .and. (row /= 1 .or. col /= 1)) cycle
                result%row_idx(k) = row
                result%val(k) = values(row, col)
                k = k + 1
            end do
        end do
        result%col_ptr(n + 1) = k
    end subroutine
    subroutine real_operations()
        real(dp) :: rhs_block(n, 2), answer(n, 2)
        do op = 0, 2
            if (op == 0) then
                rhs = matmul(dense, exact)
            else
                rhs = matmul(transpose(dense), exact)
            end if
            call sparse_solve_op(solver, rhs, x, op, status)
            call checked(status)
            call require(maxval(abs(x - exact)) < 2.e-12_dp, &
                         'AMD real N/T/H with refinement')
            rhs_block(:, 1) = rhs
            rhs_block(:, 2) = -2.0_dp*rhs
            call sparse_solve_block_op(solver, rhs_block, answer, op, status)
            call checked(status)
            call require(maxval(abs(answer(:, 1) - exact)) < 2.e-12_dp .and. &
                         maxval(abs(answer(:, 2) + 2.0_dp*exact)) < 4.e-12_dp, &
                         'AMD real block N/T/H oracle')
        end do
    end subroutine
    subroutine complex_operations()
        complex(dp) :: rhs_block(n, 2), answer(n, 2)
        do op = 0, 2
            select case (op)
            case (0)
                zrhs = matmul(zdense, zexact)
            case (1)
                zrhs = matmul(transpose(zdense), zexact)
            case (2)
                zrhs = matmul(conjg(transpose(zdense)), zexact)
            end select
            call sparse_solve_op(solver, zrhs, zx, op, status)
            call checked(status)
            call require(maxval(abs(zx - zexact)) < 2.e-12_dp, &
                         'AMD complex N/T/H refinement')
            rhs_block(:, 1) = zrhs
            rhs_block(:, 2) = cmplx(-1.0_dp, .5_dp, dp)*zrhs
            call sparse_solve_block_op(solver, rhs_block, answer, op, status)
            call checked(status)
            call require(maxval(abs(answer(:, 1) - zexact)) < 2.e-12_dp .and. &
                         maxval(abs(answer(:, 2) - &
                                    cmplx(-1.0_dp, .5_dp, dp)*zexact)) < 4.e-12_dp, &
                         'AMD complex block N/T/H oracle')
        end do
    end subroutine
end program test_general_amd
