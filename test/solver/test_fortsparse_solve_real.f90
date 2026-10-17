program test_fortsparse_solve_real
    ! Real solve on the default backend: factor once, solve two RHS
    ! reusing the factorization, and exercise sparse_solve_once. Checks the
    ! solution against the known answer and the residual A x - b.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortsparse, only: dp, csc_t, csc_z_t, csc_from_triplet, csc_matvec, &
        sparse_solver_t, sparse_factor, sparse_refactor, sparse_solve, sparse_free, &
        sparse_solve_op, sparse_solve_block_op, sparse_solve_jvp, sparse_solve_vjp, &
        sparse_solve_once, fortsparse_status_t, status_ok, &
        FORTSPARSE_BACKEND_SUPERLU_ILU, FORTSPARSE_BACKEND_SUPERLU, &
        FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_OP_NORMAL, FORTSPARSE_OP_ADJOINT
    implicit none

    integer :: nfail
    nfail = 0

    call run(nfail)
    call raw_factor(nfail)
    call ilu_scaling(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    subroutine run(nfail)
        integer, intent(inout) :: nfail

        type(csc_t)               :: A
        type(sparse_solver_t)     :: solver
        type(fortsparse_status_t) :: status
        integer                   :: rows(7), cols(7)
        real(dp)                  :: vals(7)
        real(dp)                  :: b(3), x(3), x2(3), xe(3), xe2(3)

        ! A = [[4,1,0],[1,3,1],[0,1,2]]
        rows = [1, 1, 2, 2, 2, 3, 3]
        cols = [1, 2, 1, 2, 3, 2, 3]
        vals = [4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]
        call csc_from_triplet(3, 3, rows, cols, vals, A, status)
        call check_true("build_ok", status_ok(status), nfail)

        xe = [1.0_dp, 2.0_dp, 3.0_dp]
        b = [6.0_dp, 10.0_dp, 8.0_dp]
        call sparse_factor(solver, A, status)
        call check_true("factor_ok", status_ok(status), nfail)

        call sparse_solve(solver, b, x, status)
        call check_true("solve_ok", status_ok(status), nfail)
        call check_err("solve_x", x, xe, nfail)
        call check_residual("solve_res", A, x, b, nfail)

        ! Second RHS reuses the factorization.
        xe2 = [0.0_dp, 1.0_dp, 0.0_dp]
        call sparse_solve(solver, csc_matvec(A, xe2), x2, status)
        call check_true("solve2_ok", status_ok(status), nfail)
        call check_err("solve2_x", x2, xe2, nfail)
        call sparse_free(solver)

        ! One-shot driver.
        call sparse_solve_once(A, b, x, status)
        call check_true("once_ok", status_ok(status), nfail)
        call check_err("once_x", x, xe, nfail)
    end subroutine run

    subroutine raw_factor(nfail)
        integer, intent(inout) :: nfail

        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        integer :: colptr(4), rowidx(7), pass
        real(dp) :: val(7), b(3), x(3)

        ! Nonsymmetric matrix; the backend must own everything it needs after
        ! factor returns, including the original matrix when refining.
        do pass = 1, 3
            colptr = [1, 3, 6, 8]
            rowidx = [1, 2, 1, 2, 3, 2, 3]
            val = [4.0_dp, -2.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 2.0_dp, 5.0_dp]
            solver%refine = pass == 2
            call sparse_factor(solver, 3, 3, 7, colptr, rowidx, val, status)
            call check_true("raw_factor", status_ok(status), nfail)
            val = 0.0_dp
            rowidx = 1
            colptr = 1
            b = [6.0_dp, 10.0_dp, 17.0_dp]
            call sparse_solve(solver, b, x, status)
            call check_true("raw_solve", status_ok(status), nfail)
            call check_err("raw_x", x, [1.0_dp, 2.0_dp, 3.0_dp], nfail)
            call sparse_solve(solver, 2.0_dp*b, x, status)
            call check_true("raw_solve_again", status_ok(status), nfail)
            call check_err("raw_x_again", x, [2.0_dp, 4.0_dp, 6.0_dp], nfail)
        end do
        call sparse_free(solver)
    end subroutine raw_factor

    subroutine ilu_scaling(nfail)
        integer, intent(inout) :: nfail
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        integer :: colptr(4), rowidx(7), pass, j, k, operation
        logical :: reused
        type(csc_t) :: A
        type(csc_z_t) :: Az
        real(dp) :: val(7), dense(3,3), base(3,3), rs(3), cs(3)
        real(dp) :: want(3), b(3), x(3), scale, wt(3), bt(3), block_b(3,2), block_x(3,2)
        real(dp) :: tangent(3), values_bar(7)

        base = reshape([4.0_dp, -2.0_dp, 0.0_dp, 1.0_dp, 3.0_dp, &
            1.0_dp, 0.0_dp, 2.0_dp, 5.0_dp], [3,3])
        colptr = [1, 3, 6, 8]
        rowidx = [1, 2, 1, 2, 3, 2, 3]
        solver%backend_id = FORTSPARSE_BACKEND_SUPERLU_ILU
        do pass = 1, 4
            rs = 1.0_dp
            cs = 1.0_dp
            if (pass == 2 .or. pass == 4) rs = [1.0e-8_dp, 1.0_dp, 1.0e8_dp]
            if (pass == 3 .or. pass == 4) cs = [1.0e6_dp, 1.0e-3_dp, 1.0e3_dp]
            do j = 1, 3
                dense(:,j) = rs*base(:,j)*cs(j)
                do k = colptr(j), colptr(j + 1) - 1
                    val(k) = dense(rowidx(k),j)
                end do
            end do
            want = [1.0_dp, 2.0_dp, 3.0_dp]/cs
            b = matmul(dense, want)
            call sparse_factor(solver, 3, 3, 7, colptr, rowidx, val, status)
            call check_true("ilu_factor", status_ok(status), nfail)
            call sparse_solve(solver, b, x, status)
            call check_true("ilu_solve", status_ok(status), nfail)
            call check_true("ilu_scaled_solution", &
                maxval(abs((x - want)/want)) < 1.0e-12_dp, nfail)
            scale = maxval(sum(abs(dense), dim=2))*maxval(abs(x)) + maxval(abs(b))
            call check_true("ilu_original_residual", &
                maxval(abs(matmul(dense, x) - b))/scale < 1.0e-14_dp, nfail)
            call sparse_solve(solver, 2.0_dp*b, x, status)
            call check_true("ilu_repeated", status_ok(status), nfail)
            call check_true("ilu_repeated_solution", &
                maxval(abs((x - 2.0_dp*want)/want)) < 1.0e-12_dp, nfail)
            ! A small matrix with no dropped fill supplies a dense independent
            ! control for every operation, including both equilibration sides.
            x = b
            call sparse_solve(solver, x, status)
            call check_true("ilu_inplace_status", status_ok(status), nfail)
            call check_true("ilu_inplace_solution", &
                maxval(abs((x - want)/want)) < 1.0e-12_dp, nfail)
            do operation = FORTSPARSE_OP_NORMAL, FORTSPARSE_OP_ADJOINT
                if (operation == FORTSPARSE_OP_NORMAL) then
                    wt = want
                    bt = b
                else
                    wt = [1.0_dp, 2.0_dp, 3.0_dp]/rs
                    bt = matmul(transpose(dense), wt)
                end if
                call sparse_solve_op(solver, bt, x, operation, status)
                call check_true("ilu_operation_status", status_ok(status), nfail)
                call check_true("ilu_operation_solution", &
                    maxval(abs((x - wt)/wt)) < 1.0e-12_dp, nfail)
                block_b(:,1) = bt
                block_b(:,2) = -2.0_dp*bt
                call sparse_solve_block_op(solver, block_b, block_x, operation, status)
                call check_true("ilu_block_status", status_ok(status), nfail)
                call check_true("ilu_block_solution", &
                    maxval(abs((block_x(:,1) - wt)/wt)) < 1.0e-12_dp .and. &
                    maxval(abs((block_x(:,2) + 2.0_dp*wt)/wt)) < 2.0e-12_dp, nfail)
            end do
        end do

        ! Changed values must build a fresh equilibrated ILU, not exact dgstrf
        ! factors through the direct backend's same-pattern reuse route.
        val = 1.25_dp*val
        dense = 1.25_dp*dense
        b = matmul(dense, want)
        call sparse_refactor(solver, 3, 3, 7, colptr, rowidx, val, status, reused)
        call check_true("ilu_refactor_status", status_ok(status), nfail)
        call check_true("ilu_refactor_not_reused", .not. reused, nfail)
        call sparse_solve(solver, b, x, status)
        call check_true("ilu_refactor_solution", &
            status_ok(status) .and. maxval(abs((x - want)/want)) < 1.0e-12_dp, nfail)

        A%nrow = 3
        A%ncol = 3
        A%nnz = 7
        A%col_ptr = colptr
        A%row_idx = rowidx
        A%val = val
        ! Selection is mutable, but the retained approximate factors must never
        ! be admitted as derivatives of the exact equation A*x=b.
        solver%backend_id = FORTSPARSE_BACKEND_SUPERLU
        call sparse_solve_jvp(solver, A, want, b, tangent, status)
        call check_true("ilu_jvp_rejected_retained_mode", &
            status%code == FORTSPARSE_BACKEND_UNAVAILABLE, nfail)
        call sparse_solve_vjp(solver, A, want, b, tangent, values_bar, status)
        call check_true("ilu_vjp_rejected_retained_mode", &
            status%code == FORTSPARSE_BACKEND_UNAVAILABLE, nfail)
        call sparse_solve_vjp(solver, A, want, b, tangent, values_bar, status, &
            use_primal_factors=.true.)
        call check_true("ilu_primal_vjp_rejected", &
            status%code == FORTSPARSE_BACKEND_UNAVAILABLE, nfail)
        solver%backend_id = FORTSPARSE_BACKEND_SUPERLU_ILU

        ! Unsupported requests invalidate factors and allow later recovery.
        solver%refine = .true.
        call sparse_factor(solver, A, status)
        call check_true("ilu_refinement_rejected", &
            status%code == FORTSPARSE_BACKEND_UNAVAILABLE, nfail)
        call sparse_solve(solver, b, x, status)
        call check_true("ilu_rejected_factor_unavailable", &
            status%code == FORTSPARSE_NOT_FACTORED, nfail)
        solver%refine = .false.
        Az%nrow = A%nrow
        Az%ncol = A%ncol
        Az%nnz = A%nnz
        Az%col_ptr = A%col_ptr
        Az%row_idx = A%row_idx
        Az%val = cmplx(A%val, 0.0_dp, dp)
        call sparse_factor(solver, Az, status)
        call check_true("ilu_complex_rejected", &
            status%code == FORTSPARSE_BACKEND_UNAVAILABLE, nfail)
        call sparse_solve(solver, b, x, status)
        call check_true("ilu_complex_failure_unfactored", &
            status%code == FORTSPARSE_NOT_FACTORED, nfail)
        ! An ILU zero pivot is perturbed by SuperLU, so its approximate inverse
        ! remains applicable even when the original matrix is rank deficient.
        call sparse_factor(solver, 2, 2, 4, [1,3,5], [1,2,1,2], &
            [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], status)
        call check_true("ilu_perturbed_factor", status_ok(status), nfail)
        call sparse_solve(solver, [2.0_dp, 2.0_dp], x(1:2), status)
        call check_true("ilu_perturbed_apply", status_ok(status), nfail)
        call check_true("ilu_perturbed_finite", all(ieee_is_finite(x(1:2))), nfail)
        ! Changing from an incomplete to an exact factor must reset its mode.
        solver%backend_id = FORTSPARSE_BACKEND_SUPERLU
        call sparse_factor(solver, 3, 3, 7, colptr, rowidx, val, status)
        call check_true("ilu_to_lu_factor", status_ok(status), nfail)
        call sparse_solve(solver, b, x, status)
        call check_true("ilu_to_lu_solve", status_ok(status), nfail)
        call check_true("ilu_to_lu_solution", &
            maxval(abs((x - want)/want)) < 1.0e-12_dp, nfail)
        solver%backend_id = FORTSPARSE_BACKEND_SUPERLU_ILU
        A%val = 0.0_dp
        call sparse_solve_jvp(solver, A, want, 0.0_dp*b, tangent, status)
        call check_true("lu_jvp_retained_exact_mode", status_ok(status), nfail)
        call check_true("lu_jvp_zero", all(tangent == 0.0_dp), nfail)
        A%val = val
        call sparse_factor(solver, A, status)
        call check_true("lu_to_ilu_factor", status_ok(status), nfail)
        call sparse_solve(solver, b, x, status)
        call check_true("lu_to_ilu_solution", &
            status_ok(status) .and. maxval(abs((x - want)/want)) < 1.0e-12_dp, nfail)
        call sparse_free(solver)
        call sparse_free(solver)
    end subroutine ilu_scaling

    subroutine check_err(label, got, want, nfail)
        character(*), intent(in)    :: label
        real(dp),     intent(in)    :: got(:), want(:)
        integer,      intent(inout) :: nfail
        if (maxval(abs(got - want)) > 1.0e-9_dp) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,es13.6)") "FAIL [", label, &
                "] max|err| ", maxval(abs(got - want))
        end if
    end subroutine check_err

    subroutine check_residual(label, A, x, b, nfail)
        character(*), intent(in)    :: label
        type(csc_t),  intent(in)    :: A
        real(dp),     intent(in)    :: x(:), b(:)
        integer,      intent(inout) :: nfail
        real(dp)                    :: res

        res = maxval(abs(csc_matvec(A, x) - b))
        if (res > 1.0e-9_dp) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,es13.6)") "FAIL [", label, &
                "] residual ", res
        end if
    end subroutine check_residual

    subroutine check_true(label, cond, nfail)
        character(*), intent(in)    :: label
        logical,      intent(in)    :: cond
        integer,      intent(inout) :: nfail
        if (.not. cond) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a)") "FAIL [", label, "] expected .true."
        end if
    end subroutine check_true

end program test_fortsparse_solve_real
