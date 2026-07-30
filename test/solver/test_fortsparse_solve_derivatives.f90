program test_fortsparse_solve_derivatives
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: &
        csc_from_triplet, csc_t, csc_transpose, dp, fortsparse_status_t, &
        sparse_factor, sparse_free, sparse_solve, sparse_solve_jvp, &
        sparse_solve_once, sparse_solve_vjp, sparse_solver_t, status_ok
    implicit none

    integer, parameter :: rows(7) = [1, 1, 2, 2, 2, 3, 3]
    integer, parameter :: columns(7) = [1, 2, 1, 2, 3, 2, 3]
    real(dp), parameter :: values(7) = [ &
        4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]
    real(dp), parameter :: values_dot(7) = [ &
        0.2_dp, -0.1_dp, 0.05_dp, 0.3_dp, -0.2_dp, 0.15_dp, 0.1_dp]
    real(dp), parameter :: step = 1.0e-6_dp
    real(dp) :: A_bar(7), b(3), b_bar(3), b_dot(3), b_minus(3), b_plus(3)
    real(dp) :: lhs, rhs
    real(dp) :: x(3), x_bar(3), x_dot(3), x_minus(3), x_plus(3)
    type(csc_t) :: A, A_dot, A_minus, A_plus, transpose_A
    type(fortsparse_status_t) :: status
    type(sparse_solver_t) :: solver, transpose_solver
    integer :: failures

    failures = 0
    b = [6.0_dp, 10.0_dp, 8.0_dp]
    b_dot = [0.3_dp, -0.4_dp, 0.2_dp]
    x_bar = [0.7_dp, -0.2_dp, 0.5_dp]
    call build_matrix(values, A)
    call build_matrix(values_dot, A_dot)
    call sparse_factor(solver, A, status)
    call check(status_ok(status), "primal factorization succeeds")
    call sparse_solve(solver, b, x, status)
    call check(status_ok(status), "primal solve succeeds")

    call sparse_solve_jvp(solver, A_dot, x, b_dot, x_dot, status)
    call check(status_ok(status), "implicit solve JVP succeeds")
    call build_matrix(values + step*values_dot, A_plus)
    call build_matrix(values - step*values_dot, A_minus)
    b_plus = b + step*b_dot
    b_minus = b - step*b_dot
    call sparse_solve_once(A_plus, b_plus, x_plus, status)
    call sparse_solve_once(A_minus, b_minus, x_minus, status)
    call check(maxval(abs( &
        x_dot - (x_plus - x_minus)/(2.0_dp*step))) < 2.0e-9_dp, &
        "implicit JVP matches an independent central difference")

    call csc_transpose(A, transpose_A, status)
    call sparse_factor(transpose_solver, transpose_A, status)
    call sparse_solve_vjp( &
        transpose_solver, A, x, x_bar, b_bar, A_bar, status)
    call check(status_ok(status), "implicit solve VJP succeeds")
    lhs = dot_product(x_bar, x_dot)
    rhs = dot_product(b_bar, b_dot) + dot_product(A_bar, A_dot%val)
    call check(abs(lhs - rhs) < 2.0e-12_dp, &
        "solve JVP and VJP satisfy the adjoint identity")
    call check(abs(lhs - dot_product( &
        x_bar, (x_plus - x_minus)/(2.0_dp*step))) < 2.0e-9_dp, &
        "solve VJP matches a scalar central-difference objective")

    call sparse_free(solver)
    call sparse_free(transpose_solver)
    if (failures > 0) then
        write (error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"

contains

    subroutine build_matrix(entries, matrix)
        real(dp), intent(in) :: entries(:)
        type(csc_t), intent(out) :: matrix

        call csc_from_triplet( &
            3, 3, rows, columns, entries, matrix, status)
        if (.not. status_ok(status)) error stop "matrix construction failed"
    end subroutine build_matrix

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write (error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_fortsparse_solve_derivatives
