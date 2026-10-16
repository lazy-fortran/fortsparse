program test_fortsparse_umfpack_controls
    ! Real UMFPACK, not a mocked helper: reject before numeric, recover on the
    ! same solver, and solve well-conditioned systems with tiny diagonal pivots.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: dp, i8, sparse_solver_t, sparse_factor_stats_t, &
        fortsparse_status_t, status_ok, sparse_factor, sparse_solve, &
        sparse_factor_stats, sparse_free, FORTSPARSE_BACKEND_UMFPACK_IPC, &
        FORTSPARSE_RESOURCE_LIMIT, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_BACKEND_UNAVAILABLE
    implicit none
    integer :: failures = 0, helper_length
    character(16) :: mode
    character(4096) :: helper

    call get_environment_variable('FORTSPARSE_UMFPACK_HELPER', helper, helper_length)
    if (helper_length == 0) then
        print '(a)', 'SKIP: FORTSPARSE_UMFPACK_HELPER not set (MIT-only build)'
        stop 0
    end if
    call get_command_argument(1, mode)
    if (trim(mode) == 'legacy') then
        call check_legacy()
    else
        call check_real()
        call check_complex()
    end if
    if (failures /= 0) then
        write (error_unit, '(i0,a)') failures, ' checks failed'
        stop 1
    end if
    print '(a)', 'PASS: UMFPACK pivot control, symbolic budget and recovery'
contains
    subroutine check_legacy()
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        character(4096) :: receipt
        integer :: unit, ios, requests
        solver%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC
        solver%options%memory_budget_bytes = 1_i8
        call sparse_factor(solver, 1, 1, 1, [1, 2], [1], [1.0_dp], status)
        call require(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
            'legacy helper rejected before factor')
        call sparse_free(solver)
        call get_environment_variable('FORTSPARSE_LEGACY_RECEIPT', receipt)
        open (newunit=unit, file=trim(receipt), status='old', iostat=ios)
        call require(ios == 0, 'legacy request receipt exists')
        if (ios /= 0) return
        read (unit, *, iostat=ios) requests
        close (unit, status='delete')
        call require(ios == 0, 'legacy request count readable')
        call require(requests == 0, 'no factor request reached legacy helper')
    end subroutine check_legacy

    subroutine check_real()
        type(sparse_solver_t) :: solver
        type(sparse_factor_stats_t) :: stats
        type(fortsparse_status_t) :: status
        integer :: cp(4), ri(7)
        real(dp) :: a(7), b(3), x(3), expected(3)
        integer(i8) :: estimate

        ! A = [1,1,0; 1,eps,1; 0,1,1], cond_2 tends to 2 as eps -> 0.
        cp = [1, 3, 6, 8]
        ri = [1, 2, 1, 2, 3, 2, 3]
        a = [1.0_dp, 1.0_dp, 1.0_dp, 1.e-16_dp, 1.0_dp, 1.0_dp, 1.0_dp]
        expected = [1.0_dp, 2.0_dp, 3.0_dp]
        b = [3.0_dp, 4.0_dp + 2.e-16_dp, 5.0_dp]
        solver%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC
        solver%options%memory_budget_bytes = 1_i8
        call sparse_factor(solver, 3, 3, 7, cp, ri, a, status)
        call require(status%code == FORTSPARSE_RESOURCE_LIMIT, 'real reject')
        call sparse_factor_stats(solver, stats, status)
        call require(status_ok(status), 'real rejected stats available')
        call check_rejected(stats, 'real')
        call require(abs(stats%effective_pivot_tolerance - 0.1_dp) < &
            epsilon(1.0_dp), 'native default pivot unchanged')
        estimate = stats%estimated_peak_bytes
        call sparse_solve(solver, b, x, status)
        call require(status%code == FORTSPARSE_NOT_FACTORED, 'reject not factored')
        solver%options%memory_budget_bytes = estimate
        solver%options%pivot_tolerance = 1.0_dp
        call sparse_factor(solver, 3, 3, 7, cp, ri, a, status)
        call require(status_ok(status), 'real accepted estimate boundary')
        call sparse_factor_stats(solver, stats, status)
        call check_accepted(stats, 'real')
        call sparse_solve(solver, b, x, status)
        call require(status_ok(status), 'real recovery solve')
        call require(maxval(abs(x - expected)) < 5.e-14_dp, 'real known solution')
        call sparse_free(solver)
    end subroutine check_real

    subroutine check_complex()
        type(sparse_solver_t) :: solver
        type(sparse_factor_stats_t) :: stats
        type(fortsparse_status_t) :: status
        integer :: cp(4), ri(7), k
        real(dp) :: ar(7)
        complex(dp) :: a(7), b(3), x(3), expected(3), phase
        integer(i8) :: estimate

        cp = [1, 3, 6, 8]
        ri = [1, 2, 1, 2, 3, 2, 3]
        ar = [1.0_dp, 1.0_dp, 1.0_dp, 1.e-16_dp, 1.0_dp, 1.0_dp, 1.0_dp]
        phase = cmplx(1.0_dp, 0.75_dp, dp)
        a = ar*phase
        expected = cmplx([1.0_dp, 2.0_dp, 3.0_dp], &
            [-1.0_dp, 0.5_dp, 2.0_dp], dp)
        ! Independent analytic product of the three rows, no CSC matvec.
        b = phase*[expected(1) + expected(2), &
            expected(1) + 1.e-16_dp*expected(2) + expected(3), &
            expected(2) + expected(3)]
        solver%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC
        solver%options%pivot_tolerance = 1.0_dp
        call sparse_factor(solver, 3, 3, 7, cp, ri, a, status)
        call require(status_ok(status), 'complex factor')
        call sparse_factor_stats(solver, stats, status)
        call check_accepted(stats, 'complex')
        estimate = stats%estimated_peak_bytes
        call sparse_solve(solver, b, x, status)
        call require(status_ok(status), 'complex initial solve')
        call require(maxval(abs(x - expected)) < 5.e-14_dp, 'complex known solution')
        ! Rejecting an attempted refactor must invalidate the old live factor.
        do k = 1, 2
            solver%options%memory_budget_bytes = 1_i8
            call sparse_factor(solver, 3, 3, 7, cp, ri, a, status)
            call require(status%code == FORTSPARSE_RESOURCE_LIMIT, 'complex reject')
            call sparse_factor_stats(solver, stats, status)
            call check_rejected(stats, 'complex')
            call sparse_solve(solver, b, x, status)
            call require(status%code == FORTSPARSE_NOT_FACTORED, 'no stale factor')
            solver%options%memory_budget_bytes = estimate
            call sparse_factor(solver, 3, 3, 7, cp, ri, a, status)
            call require(status_ok(status), 'complex recover factor')
            call sparse_solve(solver, b, x, status)
            call require(status_ok(status), 'complex recover solve')
            call require(maxval(abs(x - expected)) < 5.e-14_dp, 'complex recovery')
        end do
        solver%options%pivot_tolerance = 1.01_dp
        call sparse_factor(solver, 3, 3, 7, cp, ri, a, status)
        call require(status%code == FORTSPARSE_INVALID_MATRIX, 'invalid pivot')
        call sparse_free(solver)
    end subroutine check_complex

    subroutine check_rejected(stats, label)
        type(sparse_factor_stats_t), intent(in) :: stats
        character(*), intent(in) :: label
        call require(stats%symbolic_peak_bytes > 0_i8, label//' symbolic memory')
        call require(stats%estimated_peak_bytes > 1_i8, label//' estimate')
        call require(stats%ipc_bytes > 0_i8, label//' separate IPC bytes')
        call require(.not. stats%numeric_started, label//' no numeric call')
        call require(stats%numeric_peak_bytes == -1_i8, label//' numeric unknown')
    end subroutine check_rejected

    subroutine check_accepted(stats, label)
        type(sparse_factor_stats_t), intent(in) :: stats
        character(*), intent(in) :: label
        call require(stats%numeric_started, label//' numeric called')
        call require(stats%numeric_peak_bytes > 0_i8, label//' numeric memory')
        call require(stats%estimated_peak_bytes > 0_i8, label//' symbolic estimate')
        call require(stats%effective_pivot_tolerance == 1.0_dp, label//' pivot receipt')
        print '(a,3(1x,i0))', label//' UMFPACK estimate/actual/IPC bytes:', &
            stats%estimated_peak_bytes, stats%numeric_peak_bytes, stats%ipc_bytes
    end subroutine check_accepted

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (condition) return
        failures = failures + 1
        write (error_unit, '(a)') 'FAIL: '//message
    end subroutine require
end program test_fortsparse_umfpack_controls
