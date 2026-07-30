program test_fortsparse_complex_solve_derivatives
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: &
        csc_conjugate_transpose, csc_from_triplet, csc_z_t, dp, &
        fortsparse_status_t, sparse_factor, sparse_free, sparse_solve, &
        sparse_solve_jvp, sparse_solve_once, sparse_solve_vjp, &
        sparse_solver_t, status_ok
    implicit none

    integer, parameter :: rows(7) = [1, 1, 2, 2, 2, 3, 3]
    integer, parameter :: columns(7) = [1, 2, 1, 2, 3, 2, 3]
    complex(dp), parameter :: values(7) = [ &
        cmplx(4.0_dp, 0.2_dp, dp), cmplx(1.0_dp, -0.3_dp, dp), &
        cmplx(0.7_dp, 0.4_dp, dp), cmplx(3.0_dp, -0.1_dp, dp), &
        cmplx(1.0_dp, 0.2_dp, dp), cmplx(0.8_dp, -0.5_dp, dp), &
        cmplx(2.0_dp, 0.3_dp, dp)]
    complex(dp), parameter :: values_dot(7) = [ &
        cmplx(0.2_dp, -0.05_dp, dp), cmplx(-0.1_dp, 0.02_dp, dp), &
        cmplx(0.05_dp, 0.03_dp, dp), cmplx(0.3_dp, -0.04_dp, dp), &
        cmplx(-0.2_dp, 0.06_dp, dp), cmplx(0.15_dp, -0.02_dp, dp), &
        cmplx(0.1_dp, 0.05_dp, dp)]
    real(dp), parameter :: step = 1.0e-6_dp
    complex(dp) :: A_bar(7), b(3), b_bar(3), b_dot(3), b_minus(3), b_plus(3)
    complex(dp) :: x(3), x_bar(3), x_dot(3), x_minus(3), x_plus(3)
    real(dp) :: lhs, rhs
    type(csc_z_t) :: A, A_adjoint, A_dot, A_minus, A_plus
    type(fortsparse_status_t) :: status
    type(sparse_solver_t) :: adjoint_solver, solver
    integer :: failures

    failures = 0
    b = [ &
        cmplx(6.0_dp, 0.5_dp, dp), cmplx(10.0_dp, -0.7_dp, dp), &
        cmplx(8.0_dp, 0.2_dp, dp)]
    b_dot = [ &
        cmplx(0.3_dp, -0.1_dp, dp), cmplx(-0.4_dp, 0.2_dp, dp), &
        cmplx(0.2_dp, 0.05_dp, dp)]
    x_bar = [ &
        cmplx(0.7_dp, 0.1_dp, dp), cmplx(-0.2_dp, 0.3_dp, dp), &
        cmplx(0.5_dp, -0.4_dp, dp)]
    call build_matrix(values, A)
    call build_matrix(values_dot, A_dot)
    call sparse_factor(solver, A, status)
    call check(status_ok(status), "complex primal factorization succeeds")
    call sparse_solve(solver, b, x, status)
    call check(status_ok(status), "complex primal solve succeeds")

    call sparse_solve_jvp(solver, A_dot, x, b_dot, x_dot, status)
    call check(status_ok(status), "complex implicit solve JVP succeeds")
    call build_matrix(values + step*values_dot, A_plus)
    call build_matrix(values - step*values_dot, A_minus)
    b_plus = b + step*b_dot
    b_minus = b - step*b_dot
    call sparse_solve_once(A_plus, b_plus, x_plus, status)
    call sparse_solve_once(A_minus, b_minus, x_minus, status)
    call check(maxval(abs( &
        x_dot - (x_plus - x_minus)/(2.0_dp*step))) < 2.0e-9_dp, &
        "complex implicit JVP matches an independent central difference")

    call csc_conjugate_transpose(A, A_adjoint, status)
    call sparse_factor(adjoint_solver, A_adjoint, status)
    call sparse_solve_vjp( &
        adjoint_solver, A, x, x_bar, b_bar, A_bar, status)
    call check(status_ok(status), "complex implicit solve VJP succeeds")
    lhs = real(sum(conjg(x_bar)*x_dot), dp)
    rhs = real(sum(conjg(b_bar)*b_dot) + sum(conjg(A_bar)*A_dot%val), dp)
    call check(abs(lhs - rhs) < 2.0e-12_dp, &
        "complex solve products satisfy the real adjoint identity")
    call check(abs(lhs - real(sum(conjg(x_bar)*( &
        x_plus - x_minus)/(2.0_dp*step)), dp)) < 2.0e-9_dp, &
        "complex solve VJP matches a scalar central-difference objective")

    call sparse_free(solver)
    call sparse_free(adjoint_solver)
    if (failures > 0) then
        write (error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"

contains

    subroutine build_matrix(entries, matrix)
        complex(dp), intent(in) :: entries(:)
        type(csc_z_t), intent(out) :: matrix

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

end program test_fortsparse_complex_solve_derivatives
