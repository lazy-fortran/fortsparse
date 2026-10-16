program test_fortsparse_structured_raw
    use fortsparse
    implicit none
    type(sparse_solver_t) :: solver
    type(fortsparse_status_t) :: status
    type(sparse_factor_stats_t) :: stats
    integer :: ptr(3), rows(4), family, operation
    real(dp) :: rv(4), rx(2), rb(2), re(2)
    complex(dp) :: a(2, 2), values(4), exact(2), rhs(2), x(2)
    logical :: reused
    ptr = [1, 3, 5]; rows = [1, 2, 1, 2]
    solver%backend_id = FORTSPARSE_BACKEND_NATIVE
    solver%refine = .true.
    solver%options%memory_policy = FORTSPARSE_MEMORY_LOW
    solver%options%host_memory_limit_bytes = 100000_i8
    solver%options%symmetry = FORTSPARSE_SYM_SPD
    rv = [4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp]; re = [2.0_dp, &
 -3.0_dp]; rb = [5.0_dp, -7.0_dp]
    call sparse_factor(solver, 2, 2, 4, ptr, rows, rv, status)
    call check(status_ok(status), 'real SPD factor')
    call sparse_solve(solver, rb, rx, status)
    call check(status_ok(status) .and. maxval(abs(rx - re)) < 1e-13_dp, &
               'real SPD known solution')
    rv = -rv
    call sparse_factor(solver, 2, 2, 4, ptr, rows, rv, status)
    call check(status%code == FORTSPARSE_NOT_POSITIVE_DEFINITE, &
               'real SPD nonpositive status')
    call sparse_destroy(solver)
    do family = 1, 3
        solver%options%symmetry = FORTSPARSE_SYM_SPD
        a = reshape([cmplx(4.0_dp, 0.0_dp, dp), cmplx(1.0_dp, -1.0_dp, dp), &
                     cmplx(1.0_dp, 1.0_dp, dp), cmplx(3.0_dp, 0.0_dp, dp)], [2, 2])
        if (family == 2) then
            solver%options%symmetry = FORTSPARSE_SYM_HERMITIAN
            a(1, 1) = 0; a(2, 2) = 0
        else if (family == 3) then
            solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
            a(2, 1) = a(1, 2); a(1, 1) = cmplx(0.0_dp, 1.0_dp, dp)
        end if
        values = reshape(a, [4]); exact = [cmplx(2.0_dp, 1.0_dp, dp), &
 cmplx(-3.0_dp, 2.0_dp, dp)]
        call sparse_factor(solver, 2, 2, 4, ptr, rows, values, status)
        call check(status_ok(status), 'complex structured factor')
        call sparse_factor_stats(solver, stats, status)
        call check(stats%strict_memory_supported .and. stats%host_peak_bytes > 0, &
                   'complex raw strict stats')
        values = 99 ! factor/checker must own their retained data
        do operation = 0, 2
            select case (operation)
            case (0)
                rhs = matmul(a, exact)
            case (1)
                rhs = matmul(transpose(a), exact)
            case (2)
                rhs = matmul(conjg(transpose(a)), exact)
            end select
            call sparse_solve_op(solver, rhs, x, operation, status)
            call check(status_ok(status) .and. maxval(abs(x - exact)) < 1e-12_dp, &
                       'complex known N/T/H')
        end do
        values = reshape(2*a, [4])
        call sparse_refactor(solver, 2, 2, 4, ptr, rows, values, status, reused)
        call check(status_ok(status) .and. reused, 'complex structured refactor')
        rhs = matmul(2*a, exact)
        call sparse_solve(solver, rhs, x, status)
        call check(status_ok(status) .and. maxval(abs(x - exact)) < 1e-12_dp, &
                   'complex changed values')
        call sparse_factor(solver, 2, 2, 4, ptr, rows(:3), values, status)
        call check(status%code == FORTSPARSE_INVALID_MATRIX, 'short complex row extent')
        call sparse_solve(solver, rhs, x, status)
        call check(status%code == FORTSPARSE_NOT_FACTORED, &
                   'failed complex raw invalidates factor')
        values = reshape(a, [4]); values(2) = cmplx(17.0_dp, 2.0_dp, dp)
        call sparse_factor(solver, 2, 2, 4, ptr, rows, values, status)
        call check(status%code == FORTSPARSE_INVALID_MATRIX, &
                   'exact declared complex symmetry validation')
        values = reshape(a, [4])
        solver%options%host_memory_limit_bytes = 1_i8
        call sparse_factor(solver, 2, 2, 4, ptr, rows, values, status)
        call check(status%code == FORTSPARSE_RESOURCE_LIMIT, 'complex copy cap refusal')
        solver%options%host_memory_limit_bytes = 100000_i8
        call sparse_factor(solver, 2, 2, 4, ptr, rows, values, status)
        call check(status_ok(status), 'complex cap recovery')
        call sparse_destroy(solver)
    end do
contains
    subroutine check(condition, message)
        logical, intent(in)::condition
        character(*), intent(in)::message
        if (.not. condition) then
            print *, message, status%code, trim(status%msg)
            error stop 'public structured RAW oracle'
        end if
    end subroutine
end program
