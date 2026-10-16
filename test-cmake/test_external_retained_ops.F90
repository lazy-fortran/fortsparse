program test_external_retained_ops
    use fortsparse, only: dp, csc_t, csc_z_t, sparse_solver_t, &
        fortsparse_status_t, status_ok, csc_from_triplet, sparse_factor, &
        sparse_solve, sparse_solve_op, sparse_solve_block_op, sparse_free, &
        FORTSPARSE_BACKEND_MUMPS, FORTSPARSE_BACKEND_UNAVAILABLE, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_SYM_GENERAL, FORTSPARSE_SYM_HERMITIAN, FORTSPARSE_SYM_INDEFINITE
#ifdef FS_TEST_MUMPS
    use fs_test_mpi, only: test_mpi_init, test_mpi_finalize
#endif
    implicit none
    integer :: backend, refinement, kind, operation, i
    integer, parameter :: rows(9) = [1, 2, 3, 1, 2, 3, 1, 2, 3]
    integer, parameter :: cols(9) = [1, 1, 1, 2, 2, 2, 3, 3, 3]
    real(dp) :: a(3, 3), op_a(3, 3), rx(3), rb(3), ry(3), rq(3)
    real(dp) :: rknown(3, 2), rblock(3, 2), rresult(3, 2)
    real(dp) :: rempty(3, 0), rempty_result(3, 0), rbadshape(2, 0)
    complex(dp) :: z(3, 3), op_z(3, 3), zx(3), zb(3), zy(3), zq(3)
    complex(dp) :: zknown(3, 2), zblock(3, 2), zresult(3, 2), values(9)
    complex(dp) :: zempty(3, 0), zempty_result(3, 0), zbadshape(2, 0)
    type(csc_t) :: matrix
    type(csc_z_t) :: zmatrix
    type(sparse_solver_t) :: solver
    type(fortsparse_status_t) :: status
    character(16) :: argument
    logical :: stale

    call get_command_argument(1, argument)
    read (argument, *) backend
    call get_command_argument(2, argument)
    stale = trim(argument) == 'stale'
#ifdef FS_TEST_MUMPS
    if (backend == FORTSPARSE_BACKEND_MUMPS) call test_mpi_init()
#endif
    a = reshape([0.0_dp, 3.0_dp, 1.0_dp, 2.0_dp, 4.0_dp, -2.0_dp, &
                 1.0_dp, -1.0_dp, 5.0_dp], [3, 3])
    rknown(:, 1) = [0.7_dp, -1.2_dp, 2.1_dp]
    rknown(:, 2) = [-0.3_dp, 0.2_dp, 1.7_dp]
    zknown = cmplx(rknown, reshape([0.2_dp, -0.6_dp, 0.7_dp, &
                                    -0.4_dp, 0.1_dp, 0.9_dp], [3, 2]), dp)
    call check_unfactored_empty()
    do refinement = 0, 1
        solver%backend_id = backend
        solver%refine = refinement == 1
        solver%options%nthreads = 1
        call csc_from_triplet(3, 3, rows, cols, reshape(a, [9]), matrix, status)
        call require_ok('real matrix')
        call sparse_factor(solver, matrix, status)
        if (status%code == FORTSPARSE_BACKEND_UNAVAILABLE) stop 77
        call require_ok('real factor')
        call check_real_empty()
        do operation = 0, 2
            op_a = a
            if (operation /= 0) op_a = transpose(a)
            rb = matmul(op_a, rknown(:, 1))
            call sparse_solve_op(solver, rb, rx, operation, status)
            if (stale .and. operation /= 0) then
         if (status%code /= FORTSPARSE_BACKEND_UNAVAILABLE) error stop 'stale helper op'
                cycle
            end if
            call require_ok('real retained op')
         if (maxval(abs(rx - rknown(:, 1))) > 2e-12_dp) error stop 'real known solution'
            rblock = matmul(op_a, rknown)
            call sparse_solve_block_op(solver, rblock, rresult, operation, status)
            call require_ok('real block op')
          if (maxval(abs(rresult - rknown)) > 2e-12_dp) error stop 'real block solution'
        end do
        ! Normal legacy entry after transpose/adjoint must not inherit a flag.
        rb = matmul(a, rknown(:, 1))
        call sparse_solve(solver, rb, rx, status)
        call require_ok('real normal after op')
        if (maxval(abs(rx - rknown(:, 1))) > 2e-12_dp) error stop 'real state leak'
        if (.not. stale) then
            rq = [0.1_dp, 0.8_dp, -0.4_dp]
            call sparse_solve_op(solver, rq, ry, 2, status)
            call require_ok('real adjoint pairing')
            if (abs(dot_product(rq, rx) - dot_product(ry, rb)) > 2e-12_dp) &
                error stop 'real inverse adjoint identity'
        end if
        call sparse_free(solver)
        do kind = 0, 2
            z = cmplx(a, reshape([0.3_dp, 0.1_dp, -0.2_dp, 0.4_dp, 0.2_dp, &
                                  0.7_dp, -0.1_dp, 0.5_dp, -0.3_dp], [3, 3]), dp)
            solver%options%symmetry = FORTSPARSE_SYM_GENERAL
            if (kind == 1) then
                z = z + conjg(transpose(z))
                solver%options%symmetry = FORTSPARSE_SYM_HERMITIAN
            else if (kind == 2) then
                z = z + transpose(z)
                solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
            end if
            do i = 1, 3
                z(i, i) = z(i, i) + 10.0_dp
            end do
            values = reshape(z, [9])
            call csc_from_triplet(3, 3, rows, cols, values, zmatrix, status)
            call sparse_factor(solver, zmatrix, status)
            if (status%code == FORTSPARSE_BACKEND_UNAVAILABLE) then
                if (backend /= FORTSPARSE_BACKEND_MUMPS) &
                    error stop 'unexpected complex unavailable'
                print *, 'ZMUMPS unavailable: complex retained operations not built'
                exit
            end if
            call require_ok('complex factor')
            call check_complex_empty()
            do operation = 0, 2
                op_z = z
                if (operation == 1) op_z = transpose(z)
                if (operation == 2) op_z = conjg(transpose(z))
                zb = matmul(op_z, zknown(:, 1))
                call sparse_solve_op(solver, zb, zx, operation, status)
                if (stale .and. operation /= 0) then
        if (status%code /= FORTSPARSE_BACKEND_UNAVAILABLE) error stop 'stale complex op'
                    cycle
                end if
                call require_ok('complex retained op')
      if (maxval(abs(zx - zknown(:, 1))) > 2e-12_dp) error stop 'complex known solution'
                zblock = matmul(op_z, zknown)
                call sparse_solve_block_op(solver, zblock, zresult, operation, status)
                call require_ok('complex block op')
       if (maxval(abs(zresult - zknown)) > 2e-12_dp) error stop 'complex block solution'
            end do
            zb = matmul(z, zknown(:, 1))
            call sparse_solve(solver, zb, zx, status)
            call require_ok('complex normal after op')
          if (maxval(abs(zx - zknown(:, 1))) > 2e-12_dp) error stop 'complex state leak'
            if (.not. stale) then
                zq = cmplx([0.1_dp, 0.8_dp, -0.4_dp], [0.3_dp, -0.2_dp, 0.5_dp], dp)
                call sparse_solve_op(solver, zq, zy, 2, status)
                call require_ok('complex adjoint pairing')
                if (abs(dot_product(zq, zx) - dot_product(zy, zb)) > 2e-12_dp) &
                    error stop 'complex inverse adjoint identity'
            end if
            call sparse_free(solver)
        end do
        solver%options%symmetry = FORTSPARSE_SYM_GENERAL
    end do
    call check_unfactored_empty()
#ifdef FS_TEST_MUMPS
    if (backend == FORTSPARSE_BACKEND_MUMPS) call test_mpi_finalize()
#endif
    print *, 'PASS retained external operations'
contains
    subroutine check_real_empty()
        integer :: op
        do op = 0, 2
            call sparse_solve_block_op(solver, rempty, rempty_result, op, status)
            call require_ok('real empty block')
            call sparse_solve_block_op(solver, rempty, rbadshape, op, status)
            if (status%code /= FORTSPARSE_INVALID_MATRIX) error stop 'empty real shape'
        end do
        do op = -1, 3, 4
            call sparse_solve_block_op(solver, rempty, rempty_result, op, status)
            if (status%code /= FORTSPARSE_INVALID_MATRIX) error stop 'empty real operation'
        end do
    end subroutine check_real_empty

    subroutine check_complex_empty()
        integer :: op
        do op = 0, 2
            call sparse_solve_block_op(solver, zempty, zempty_result, op, status)
            call require_ok('complex empty block')
            call sparse_solve_block_op(solver, zempty, zbadshape, op, status)
            if (status%code /= FORTSPARSE_INVALID_MATRIX) error stop 'empty complex shape'
        end do
        do op = -1, 3, 4
            call sparse_solve_block_op(solver, zempty, zempty_result, op, status)
            if (status%code /= FORTSPARSE_INVALID_MATRIX) error stop 'empty complex operation'
        end do
    end subroutine check_complex_empty

    subroutine check_unfactored_empty()
        integer :: op
        do op = 0, 2
            call sparse_solve_block_op(solver, rempty, rempty_result, op, status)
            if (status%code /= FORTSPARSE_NOT_FACTORED) error stop 'empty real unfactored'
            call sparse_solve_block_op(solver, zempty, zempty_result, op, status)
            if (status%code /= FORTSPARSE_NOT_FACTORED) error stop 'empty complex unfactored'
        end do
    end subroutine check_unfactored_empty

    subroutine require_ok(label)
        character(*), intent(in) :: label
        if (.not. status_ok(status)) then
            print *, label, status%code, status%msg
            error stop 'external retained operation failed'
        end if
    end subroutine
end program
