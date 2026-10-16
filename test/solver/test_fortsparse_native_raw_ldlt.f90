program test_fortsparse_native_raw_ldlt
    use fortsparse
    use fortsparse_native_general, only: native_general_backend_t
    implicit none
    type(sparse_solver_t) :: solver
    type(native_general_backend_t) :: direct
    type(sparse_factor_stats_t) :: stats
    type(fortsparse_status_t) :: status
    integer :: ptr(3) = [1, 2, 3], rows(2) = [2, 1], saved_ptr(3), saved_rows(2)
    real(dp) :: values(2) = [1.0_dp, 1.0_dp], saved_values(2), exact(2) = [2.0_dp, -3.0_dp]
    real(dp) :: rhs(2) = [-3.0_dp, 2.0_dp], result(2), block(2, 2), answer(2, 2)
    logical :: reused
    integer :: structure, operation
    saved_ptr = ptr; saved_rows = rows
    solver%backend_id = FORTSPARSE_BACKEND_NATIVE
    solver%refine = .true.
    do structure = FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_SYM_HERMITIAN
        solver%options%symmetry = structure
        values = 1.0_dp; saved_values = values
        call sparse_factor(solver, 2, 2, 2, ptr, rows, values, status)
        call checked()
        call expect(all(ptr == saved_ptr) .and. all(rows == saved_rows), 'borrowed indices changed')
        call expect(all(values == saved_values), 'borrowed values changed')
        call sparse_factor_stats(solver, stats, status)
        call checked()
        call expect(stats%native_factor_family == 'LDLT_REAL', 'wrong structured factor family')
        call expect(.not. stats%strict_memory_supported .and. stats%host_peak_bytes < 0, &
                    'raw conversion claimed a complete strict peak')
        values = 17.0_dp ! Neither the factors nor retained checker may borrow these values.
        do operation = 0, 2
            call sparse_solve_op(solver, rhs, result, operation, status)
            call checked()
            call expect(maxval(abs(result - exact)) < 1e-13_dp, 'zero-diagonal pair N/T/H solve')
        end do
        values = 2.0_dp
        call sparse_refactor(solver, 2, 2, 2, ptr, rows, values, status, reused)
        call checked()
        call expect(reused, 'raw structured pattern not reused')
        block(:, 1) = 2*rhs; block(:, 2) = 4*rhs
        call sparse_solve_block(solver, block, answer, status)
        call checked()
        call expect(maxval(abs(answer(:, 1) - exact)) < 1e-13_dp, 'raw refactor values')
        call expect(maxval(abs(answer(:, 2) - 2*exact)) < 1e-13_dp, 'raw block values')
        call sparse_factor(solver, 2, 2, 2, ptr(:2), rows, values, status)
        call expect(status%code == FORTSPARSE_INVALID_MATRIX, 'short pointer extent')
        call sparse_solve(solver, rhs, result, status)
        call expect(status%code == FORTSPARSE_NOT_FACTORED, 'invalid raw factor remained usable')
        values = 0.0_dp
        call sparse_factor(solver, 2, 2, 2, ptr, rows, values, status)
        call expect(status%code == FORTSPARSE_SINGULAR, 'singular structured raw status')
        values = 1.0_dp
        call sparse_factor(solver, 2, 2, 2, ptr, rows, values, status)
        call checked()
        solver%options%host_memory_limit_bytes = 1_i8
        call sparse_refactor(solver, 2, 2, 2, ptr, rows, values, status, reused)
        call expect(status%code == FORTSPARSE_BACKEND_UNAVAILABLE .and. .not. reused, &
                    'capped conversion accepted')
        solver%options%host_memory_limit_bytes = 0_i8
        solver%options%memory_policy = FORTSPARSE_MEMORY_LOW
        call sparse_factor(solver, 2, 2, 2, ptr, rows, values, status)
        call expect(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, 'LOW conversion accepted')
        solver%options%memory_policy = FORTSPARSE_MEMORY_AUTO
        call sparse_factor(solver, 2, 2, 2, ptr, rows, values, status)
        call checked()
        call sparse_solve(solver, rhs, result, status)
        call checked()
        call expect(maxval(abs(result - exact)) < 1e-13_dp, 'raw recovery')
        call sparse_destroy(solver)
    end do
    ! The backend also guards direct raw calls, before their owned conversion.
    direct%options%symmetry = FORTSPARSE_SYM_INDEFINITE
    direct%options%host_memory_limit_bytes = 1_i8
    call direct%factor_real_raw(2, 2, 2, ptr, rows, values, .false., status)
    call expect(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, 'direct capped raw conversion')
    direct%options%host_memory_limit_bytes = 0_i8
    call direct%factor_real_raw(2, 2, 2, ptr(:2), rows, values, .false., status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX, 'direct raw extent guard')
    call direct%factor_real_raw(2, 2, 2, ptr, rows(:1), values, .false., status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX, 'direct row extent guard')
    call direct%factor_real_raw(2, huge(0), 2, ptr, rows, values, .false., status)
    call expect(status%code == FORTSPARSE_INDEX_OVERFLOW, 'direct pointer overflow guard')
    print '(a)', 'Structured raw LDLT pairs, borrowed input, refactor and cap refusal pass'
contains
    subroutine checked()
        if (.not. status_ok(status)) then
            print *, status%code, trim(status%msg)
            error stop 'unexpected structured raw status'
        end if
    end subroutine
    subroutine expect(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, message
            error stop 'structured raw oracle'
        end if
    end subroutine
end program test_fortsparse_native_raw_ldlt
