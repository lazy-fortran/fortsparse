program test_lu_row_preference
    use fortsparse
    use fortsparse_native_lu, only: lu_factor_t
    use fortsparse_native_lu_z, only: lu_factor_z_t
    implicit none
    integer, parameter :: n = 4, q(n) = [3, 1, 4, 2], identity(n) = [1, 2, 3, 4]
    type(lu_factor_t) :: f, baseline
    type(lu_factor_z_t) :: zf, zbaseline
    type(csc_t) :: a
    type(csc_z_t) :: za
    type(fortsparse_status_t) :: status
    real(dp) :: dense(n, n), values(n*n), exact(n, 2), rhs(n, 2), short_values(4), small_rhs(2)
    complex(dp) :: zdense(n, n), zvalues(n*n), zexact(n, 2), zrhs(n, 2), zsmall(4), zsmall_rhs(2)
    integer :: i, j, k, op, mode, ptr(n + 1), rows(n*n), preference(n), saved(n), bad(n)
    integer(i8) :: cap, zcap, bytes, zbytes
    do j = 1, n
        do i = 1, n
            dense(i, j) = 0.1_dp*(2*i + j)
            if (i == j) dense(i, j) = 0
            if (i == 1 + mod(j, n)) dense(i, j) = 5 + j
            zdense(i, j) = dense(i, j)*cmplx(1._dp, 0.1_dp*(i - j), dp)
            k = i + (j - 1)*n
            rows(k) = i; values(k) = dense(i, j); zvalues(k) = zdense(i, j)
        end do
    end do
    ptr = [(1 + (j - 1)*n, j=1, n + 1)]
    exact(:, 1) = [1._dp, -2._dp, 3._dp, -4._dp]; exact(:, 2) = -0.5_dp*exact(:, 1)
    zexact = cmplx(exact, 2*exact, dp)
    a%nrow = n; a%ncol = n; a%nnz = n*n
    a%col_ptr = ptr; a%row_idx = rows; a%val = values
    za%nrow = n; za%ncol = n; za%nnz = n*n
    za%col_ptr = ptr; za%row_idx = rows; za%val = zvalues
    do mode = 0, 1
        call baseline%free(); call zbaseline%free(); call f%free(); call zf%free()
        call baseline%factor(a, status, column_permutation=q, pruning=mode == 1)
        call checked()
        call zbaseline%factor(za, status, column_permutation=q, pruning=mode == 1)
        call checked()
        cap = baseline%peak_memory_bytes; bytes = baseline%memory_bytes()
        zcap = zbaseline%peak_memory_bytes; zbytes = zbaseline%memory_bytes()
        ! Identity preference exactly recovers the default diagonal hint and cap.
        call f%factor(a, status, q, memory_limit_bytes=cap, pruning=mode == 1, row_preference=identity)
        call checked()
        call zf%factor(za, status, q, memory_limit_bytes=zcap, pruning=mode == 1, row_preference=identity)
        call checked()
        call expect(f%memory_bytes() == bytes .and. f%peak_memory_bytes == cap, 'real preference added owned payload')
        call expect(zf%memory_bytes() == zbytes .and. zf%peak_memory_bytes == zcap, &
                    'complex preference added owned payload')
        ! A nontrivial original-column map, independent of requested column order.
        preference = [2, 3, 4, 1]; saved = preference
        call f%factor_raw(n, ptr, rows, values, status, q, pivot_threshold=0.1_dp, &
                          pruning=mode == 1, row_preference=preference)
        call checked()
        call zf%factor_raw(n, ptr, rows, zvalues, status, q, pivot_threshold=0.1_dp, &
                           pruning=mode == 1, row_preference=preference)
        call checked()
        call expect(all(preference == saved), 'borrowed preference mutated')
        preference = 0 ! Factor/solve retains no borrowed preference reference.
        call operations()
        do j = 1, 4
            bad = identity
            select case (j)
            case (1)
                bad(1) = 0
            case (2)
                bad(n) = n + 1
            case (3)
                bad(n) = bad(1)
            case (4)
                call f%factor(a, status, row_preference=bad(:n - 1))
                call expect(status%code == FORTSPARSE_INVALID_MATRIX, 'real short preference accepted')
                call zf%factor(za, status, row_preference=bad(:n - 1))
                call expect(status%code == FORTSPARSE_INVALID_MATRIX, 'complex short preference accepted')
                cycle
            end select
            call f%factor_raw(n, ptr, rows, values, status, row_preference=bad)
            call expect(status%code == FORTSPARSE_INVALID_MATRIX, 'real invalid preference accepted')
            call zf%factor_raw(n, ptr, rows, zvalues, status, row_preference=bad)
            call expect(status%code == FORTSPARSE_INVALID_MATRIX, 'complex invalid preference accepted')
        end do
        rhs = 1; zrhs = 1
        call f%solve(rhs, status)
        call expect(status%code == FORTSPARSE_NOT_FACTORED, 'real invalid preference left usable factor')
        call zf%solve(zrhs, status)
        call expect(status%code == FORTSPARSE_NOT_FACTORED, 'complex invalid preference left usable factor')
        call f%factor(a, status, q, row_preference=saved)
        call checked()
        call zf%factor(za, status, q, row_preference=saved)
        call checked()
        call operations()
    end do
    ! A qualifying hinted tiny pivot differs behaviorally from the default.
    short_values = [1._dp, tiny(1._dp)*0.01_dp, 1._dp, 1._dp]
    zsmall = cmplx(short_values, 0._dp, dp)
    call f%factor_raw(2, [1, 3, 5], [1, 2, 1, 2], short_values, status, pivot_threshold=0._dp, row_preference=[2, 1])
    call expect(status%code == FORTSPARSE_SINGULAR, 'real preferred tiny pivot ignored')
    call zf%factor_raw(2, [1, 3, 5], [1, 2, 1, 2], zsmall, status, pivot_threshold=0._dp, row_preference=[2, 1])
    call expect(status%code == FORTSPARSE_SINGULAR, 'complex preferred tiny pivot ignored')
    call f%factor_raw(2, [1, 3, 5], [1, 2, 1, 2], short_values, status, pivot_threshold=1._dp, row_preference=[2, 1])
    call checked()
    call zf%factor_raw(2, [1, 3, 5], [1, 2, 1, 2], zsmall, status, pivot_threshold=1._dp, row_preference=[2, 1])
    call checked()
    ! Recursive structural singletons override a dangerous numerical hint.
    short_values(3) = 0; zsmall = cmplx(short_values, 0._dp, dp)
    call f%factor_raw(2, [1, 3, 5], [1, 2, 1, 2], short_values, status, pivot_threshold=0._dp, row_preference=[2, 1])
    call checked()
    call zf%factor_raw(2, [1, 3, 5], [1, 2, 1, 2], zsmall, status, pivot_threshold=0._dp, row_preference=[2, 1])
    call checked()
    small_rhs = [2._dp, 3._dp]; zsmall_rhs = cmplx(small_rhs, 0._dp, dp)
    call f%solve(small_rhs, status); call checked()
    call zf%solve(zsmall_rhs, status); call checked()
    call expect(maxval(abs(small_rhs - [2._dp, 3._dp])) < 1e-13_dp, 'real singleton hint priority')
    call expect(maxval(abs(zsmall_rhs - cmplx([2._dp, 3._dp], 0._dp, dp))) < 1e-13_dp, &
                'complex singleton hint priority')
    call missing_hint()
    print '(a)', 'Borrowed real/complex LU row preferences, thresholds, singleton priority, caps and recovery pass'
contains
    subroutine missing_hint()
        real(dp) :: v(6), b(3)
        complex(dp) :: zv(6), zb(3)
        ! Every row has degree two; identity hints are absent, not singletons.
        v = [1._dp, 4._dp, 2._dp, 1._dp, 1._dp, 3._dp]; zv = cmplx(v, 0.2_dp*v, dp)
        call f%factor_raw(3, [1, 3, 5, 7], [2, 3, 1, 3, 1, 2], v, status, &
                          column_permutation=[2, 3, 1], pivot_threshold=0._dp, row_preference=[1, 2, 3])
        call checked()
        call zf%factor_raw(3, [1, 3, 5, 7], [2, 3, 1, 3, 1, 2], zv, status, &
                           column_permutation=[2, 3, 1], pivot_threshold=0._dp, row_preference=[1, 2, 3])
        call checked()
        b = [7._dp, 10._dp, 6._dp]; zb = cmplx(b, 0.2_dp*b, dp)
        call f%solve(b, status); call checked()
        call zf%solve(zb, status); call checked()
        call expect(maxval(abs(b - [1._dp, 2._dp, 3._dp])) < 1e-12_dp, 'missing real hint did not fall back')
        call expect(maxval(abs(zb - cmplx([1._dp, 2._dp, 3._dp], 0._dp, dp))) < 1e-12_dp, &
                    'missing complex hint did not fall back')
    end subroutine
    subroutine operations()
        do op = FORTSPARSE_OP_NORMAL, FORTSPARSE_OP_ADJOINT
            select case (op)
            case (FORTSPARSE_OP_NORMAL)
                rhs = matmul(dense, exact); zrhs = matmul(zdense, zexact)
            case (FORTSPARSE_OP_TRANSPOSE)
                rhs = matmul(transpose(dense), exact); zrhs = matmul(transpose(zdense), zexact)
            case (FORTSPARSE_OP_ADJOINT)
                rhs = matmul(transpose(dense), exact); zrhs = matmul(transpose(conjg(zdense)), zexact)
            end select
            call f%solve(rhs, status, op); call checked()
            call zf%solve(zrhs, status, op); call checked()
            call expect(maxval(abs(rhs - exact)) < 1e-11_dp, 'real preferred N/T/H block solve')
            call expect(maxval(abs(zrhs - zexact)) < 1e-11_dp, 'complex preferred N/T/H block solve')
        end do
    end subroutine
    subroutine checked()
        if (.not. status_ok(status)) then
            print *, status%code, trim(status%msg)
            error stop 'unexpected row preference status'
        end if
    end subroutine
    subroutine expect(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, message
            error stop 'row preference oracle'
        end if
    end subroutine
end program
