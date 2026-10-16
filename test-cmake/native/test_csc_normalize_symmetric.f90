program test_csc_normalize_symmetric
    use fortsparse
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_finite
    use, intrinsic :: iso_c_binding, only: c_int64_t
    implicit none
    integer, parameter :: n = 4
    integer :: ptr(n + 1), rows(n*n), i, j, k, broken(n*n), badptr(n + 1)
    real(dp) :: values(n*n), original(n*n), dense(n, n), delta, actual
    real(dp) :: pair(2), saved_pair(2), invalid_tolerance
    integer(c_int64_t) :: calls, bytes
    type(fortsparse_status_t) :: status
    type(sparse_solver_t) :: solver
    interface
        subroutine begin_count() bind(C, name='fsparse_refinement_alloc_begin')
        end subroutine
        subroutine end_count(count, requested) bind(C, name='fsparse_refinement_alloc_end')
            import c_int64_t
            integer(c_int64_t), intent(out) :: count, requested
        end subroutine
    end interface
    ptr = [(1 + (j - 1)*n, j=1, n + 1)]
    do j = 1, n
        do i = 1, n
            k = i + (j - 1)*n
            rows(k) = i
            values(k) = real(i + 3*j, dp)
            if (i < j) values(k) = real(j + 3*i, dp)*(1 + 2*epsilon(1._dp))
        end do
    end do
    original = values
    call begin_count()
    call csc_normalize_symmetric_raw(n, ptr, rows, values, 1e-12_dp, delta, status)
    call end_count(calls, bytes)
    call checked()
    call expect(calls == 0 .and. bytes == 0, 'normalization allocated storage')
    dense = reshape(values, [n, n])
    call expect(all(dense == transpose(dense)), 'output is not exactly symmetric')
    actual = maxval(abs(values - original))
    call expect(delta >= actual .and. delta > 0, 'reported correction does not enclose changes')
    do j = 1, n
        do i = j, n
            k = i + (j - 1)*n
            call expect(values(k) == original(k), 'authoritative lower/diagonal changed')
            call expect(dense(j, i) == original(k), 'upper differs from original lower oracle')
        end do
    end do
    ! Exactly symmetric zero-tolerance input is idempotent.
    call csc_normalize_symmetric_raw(n, ptr, rows, values, 0._dp, delta, status)
    call checked()
    call expect(delta == 0, 'exact symmetry reports correction')
    ! A late offending pair must reject before modifying an earlier tiny mismatch.
    values = original
    values(n + (n - 2)*n) = values(n + (n - 2)*n) + 1
    call rejected(values, 'substantive asymmetry')
    values = original
    call csc_normalize_symmetric_raw(n, ptr, rows, values, 0._dp, delta, status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX .and. all(values == original), &
                'zero tolerance accepted or mutated asymmetric input')
    values = original; values(n*n) = ieee_value(0._dp, ieee_quiet_nan)
    call rejected(values, 'nonfinite last entry')
    values = original; values(n*n) = ieee_value(0._dp, ieee_positive_inf)
    call rejected(values, 'infinite last entry')
    values = original; broken = rows; broken(n*n) = n - 1
    call csc_normalize_symmetric_raw(n, ptr, broken, values, 1e-12_dp, delta, status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX .and. all(values == original), 'duplicate row accepted')
    broken = rows; broken(1) = 0
    call csc_normalize_symmetric_raw(n, ptr, broken, values, 1e-12_dp, delta, status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX .and. all(values == original), 'invalid row accepted')
    badptr = ptr; badptr(n) = size(values) + 2
    call csc_normalize_symmetric_raw(n, badptr, rows, values, 1e-12_dp, delta, status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX .and. all(values == original), 'bad pointer accepted')
    call csc_normalize_symmetric_raw(n, ptr(:n), rows, values, 1e-12_dp, delta, status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX .and. all(values == original), 'short pointer accepted')
    call csc_normalize_symmetric_raw(n, ptr, rows(:n), values, 1e-12_dp, delta, status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX .and. all(values == original), 'short row extent accepted')
    invalid_tolerance = ieee_value(0._dp, ieee_quiet_nan)
    call csc_normalize_symmetric_raw(n, ptr, rows, values, invalid_tolerance, delta, status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX .and. all(values == original), 'NaN tolerance accepted')
    call csc_normalize_symmetric_raw(n, ptr, rows, values, -1._dp, delta, status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX .and. all(values == original), 'negative tolerance accepted')
    ! Missing reflected stored zero is still a malformed full structure.
    pair = 0
    call csc_normalize_symmetric_raw(2, [1, 2, 2], [2], pair(:1), 1e-12_dp, delta, status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX, 'missing reflected zero accepted')
    call csc_normalize_symmetric_raw(0, [1], rows(:0), values(:0), 0._dp, delta, status)
    call checked()
    call expect(delta == 0, 'empty correction')
    call csc_normalize_symmetric_raw(huge(0), ptr, rows, values, 0._dp, delta, status)
    call expect(status%code == FORTSPARSE_INDEX_OVERFLOW, 'pointer extent overflow')
    ! This also runs with GNU invalid/zero/overflow traps enabled.
    pair = [huge(0._dp), -huge(0._dp)]; saved_pair = pair
    call csc_normalize_symmetric_raw(2, [1, 2, 3], [2, 1], pair, 1._dp, delta, status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX .and. all(pair == saved_pair), 'extreme asymmetry accepted')
    call expect(.not. ieee_is_finite(delta), 'unrepresentable correction not reported')
    pair = [tiny(0._dp), nearest(tiny(0._dp), 1._dp)]
    call csc_normalize_symmetric_raw(2, [1, 2, 3], [2, 1], pair, 1e-12_dp, delta, status)
    call checked()
    call expect(pair(1) == pair(2) .and. delta > 0, 'small-scale normalization')
    ! Zero-diagonal indefinite structured factors remain strict and use the normalized pair.
    pair = [1._dp, 1._dp + epsilon(1._dp)]
    solver%backend_id = FORTSPARSE_BACKEND_NATIVE
    solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
    call sparse_factor(solver, 2, 2, 2, [1, 2, 3], [2, 1], pair, status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX, 'factor silently loosened symmetry')
    call csc_normalize_symmetric_raw(2, [1, 2, 3], [2, 1], pair, 1e-12_dp, delta, status)
    call checked()
    call sparse_factor(solver, 2, 2, 2, [1, 2, 3], [2, 1], pair, status)
    call checked()
    call sparse_destroy(solver)
    print '(a)', 'Explicit allocation-free lower-to-upper symmetry normalization passes'
contains
    subroutine rejected(input, message)
        real(dp), intent(inout) :: input(:)
        character(*), intent(in) :: message
        ! Bitwise unchanged includes the NaN payload, without floating NaN comparisons.
        integer(i8) :: saved(size(input)), after(size(input))
        saved = transfer(input, saved)
        call csc_normalize_symmetric_raw(n, ptr, rows, input, 1e-12_dp, delta, status)
        after = transfer(input, after)
        call expect(status%code == FORTSPARSE_INVALID_MATRIX .and. all(saved == after), message)
    end subroutine
    subroutine checked()
        if (.not. status_ok(status)) then
            print *, status%code, trim(status%msg)
            error stop 'unexpected normalization status'
        end if
    end subroutine
    subroutine expect(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, message
            error stop 'symmetry normalization oracle'
        end if
    end subroutine
end program
