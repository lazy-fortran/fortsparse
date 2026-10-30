program test_fortsparse_native_default_ordering
    ! Default native ordering is fill-reducing on large row-numbered 2D meshes,
    ! leaves the solution unchanged within roundoff, and NATURAL keeps the
    ! input order. Oracles: explicit natural factors and original residuals.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse
    implicit none
    integer, parameter :: grid = 65, vertices = 46
    real(dp), parameter :: fill_ratio = 0.5_dp
    type(csc_t) :: laplace
    type(csc_z_t) :: elastic
    integer :: failures

    failures = 0
    call build_laplacian(grid, 0.0_dp, laplace)
    call build_elasticity(vertices, elastic)
    call compare_real(laplace, FORTSPARSE_SYM_GENERAL, 'real LU')
    call compare_real(laplace, FORTSPARSE_SYM_SPD, 'real SPD')
    call compare_complex(elastic, 'complex LU')
    if (failures /= 0) then
        write (error_unit, '(i0,a)') failures, ' default-ordering checks failed'
        stop 1
    end if
    print '(a)', 'PASS'
    stop 0

contains

    subroutine compare_real(A, symmetry, label)
        type(csc_t), intent(in) :: A
        integer, intent(in) :: symmetry
        character(*), intent(in) :: label
        real(dp), allocatable :: b(:), x_default(:), x_natural(:)
        integer(i8) :: nnz_default, nnz_natural
        integer :: i

        allocate (b(A%nrow))
        b = [(sin(0.37_dp*real(i, dp)) + 1.0_dp, i=1, A%nrow)]
        call solve_real(A, symmetry, FORTSPARSE_ORDER_DEFAULT, b, x_default, &
            nnz_default, FORTSPARSE_ORDER_AMD, label//' default')
        call solve_real(A, symmetry, FORTSPARSE_ORDER_NATURAL, b, x_natural, &
            nnz_natural, FORTSPARSE_ORDER_DEFAULT, label//' natural')
        if (failures /= 0) return
        call check(real(nnz_default, dp) < fill_ratio*real(nnz_natural, dp), &
            label//': default fill well below natural order')
        call check(maxval(abs(x_default - x_natural)) <= &
            1.0e-10_dp*maxval(abs(x_natural)), label//': same solution')
    end subroutine compare_real

    subroutine solve_real(A, symmetry, ordering, b, x, factor_nnz, used, label)
        type(csc_t), intent(in) :: A
        integer, intent(in) :: symmetry, ordering, used
        real(dp), intent(in) :: b(:)
        real(dp), allocatable, intent(out) :: x(:)
        integer(i8), intent(out) :: factor_nnz
        character(*), intent(in) :: label
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: stats
        real(dp), allocatable :: r(:)

        allocate (x(size(b)), r(size(b)))
        factor_nnz = -1_i8
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry = symmetry
        solver%options%ordering = ordering
        call sparse_factor(solver, A, status)
        call check(status_ok(status), label//': factor')
        if (.not. status_ok(status)) return
        call sparse_solve(solver, b, x, status)
        call check(status_ok(status), label//': solve')
        call sparse_factor_stats(solver, stats, status)
        call check(stats%ordering_used == used, label//': reported ordering')
        factor_nnz = stats%factor_nnz
        call check(factor_nnz > 0_i8, label//': factor entries reported')
        r = csc_matvec(A, x) - b
        call check(maxval(abs(r)) <= 1.0e-12_dp*(8.0_dp*maxval(abs(x)) + &
            maxval(abs(b))), label//': original residual')
        call sparse_free(solver)
    end subroutine solve_real

    subroutine compare_complex(A, label)
        type(csc_z_t), intent(in) :: A
        character(*), intent(in) :: label
        complex(dp), allocatable :: b(:), x_default(:), x_natural(:)
        integer(i8) :: nnz_default, nnz_natural
        integer :: i

        allocate (b(A%nrow))
        b = [(cmplx(cos(0.21_dp*real(i, dp)), 1.0_dp/real(i, dp), dp), i=1, A%nrow)]
        call solve_complex(A, FORTSPARSE_ORDER_DEFAULT, b, x_default, nnz_default, &
            FORTSPARSE_ORDER_AMD, label//' default')
        call solve_complex(A, FORTSPARSE_ORDER_NATURAL, b, x_natural, nnz_natural, &
            FORTSPARSE_ORDER_DEFAULT, label//' natural')
        if (failures /= 0) return
        call check(real(nnz_default, dp) < fill_ratio*real(nnz_natural, dp), &
            label//': default fill well below natural order')
        call check(maxval(abs(x_default - x_natural)) <= &
            1.0e-10_dp*maxval(abs(x_natural)), label//': same solution')
    end subroutine compare_complex

    subroutine solve_complex(A, ordering, b, x, factor_nnz, used, label)
        type(csc_z_t), intent(in) :: A
        integer, intent(in) :: ordering, used
        complex(dp), intent(in) :: b(:)
        complex(dp), allocatable, intent(out) :: x(:)
        integer(i8), intent(out) :: factor_nnz
        character(*), intent(in) :: label
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: stats
        complex(dp), allocatable :: r(:)

        allocate (x(size(b)), r(size(b)))
        factor_nnz = -1_i8
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%ordering = ordering
        call sparse_factor(solver, A, status)
        call check(status_ok(status), label//': factor')
        if (.not. status_ok(status)) return
        call sparse_solve(solver, b, x, status)
        call check(status_ok(status), label//': solve')
        call sparse_factor_stats(solver, stats, status)
        call check(stats%ordering_used == used, label//': reported ordering')
        factor_nnz = stats%factor_nnz
        call check(factor_nnz > 0_i8, label//': factor entries reported')
        r = csc_matvec(A, x) - b
        call check(maxval(abs(r)) <= 1.0e-12_dp*(20.0_dp*maxval(abs(x)) + &
            maxval(abs(b))), label//': original residual')
        call sparse_free(solver)
    end subroutine solve_complex

    subroutine build_laplacian(m, shift, A)
        ! Five-point Laplacian, row-wise numbering: bandwidth m.
        integer, intent(in) :: m
        real(dp), intent(in) :: shift
        type(csc_t), intent(out) :: A
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        type(fortsparse_status_t) :: status
        integer :: ix, iy, k, c

        allocate (rows(5*m*m), cols(5*m*m), vals(5*m*m))
        k = 0
        do iy = 1, m
            do ix = 1, m
                c = (iy - 1)*m + ix
                call push(rows, cols, vals, k, c, c, 4.0_dp - shift)
                if (ix > 1) call push(rows, cols, vals, k, c, c - 1, -1.0_dp)
                if (ix < m) call push(rows, cols, vals, k, c, c + 1, -1.0_dp)
                if (iy > 1) call push(rows, cols, vals, k, c, c - m, -1.0_dp)
                if (iy < m) call push(rows, cols, vals, k, c, c + m, -1.0_dp)
            end do
        end do
        call csc_from_triplet(m*m, m*m, rows(:k), cols(:k), vals(:k), A, status)
        call check(status_ok(status), 'Laplacian construction')
    end subroutine build_laplacian

    subroutine push(rows, cols, vals, k, i, j, v)
        integer, intent(inout) :: rows(:), cols(:), k
        real(dp), intent(inout) :: vals(:)
        integer, intent(in) :: i, j
        real(dp), intent(in) :: v
        k = k + 1
        rows(k) = i
        cols(k) = j
        vals(k) = v
    end subroutine push

    subroutine build_elasticity(m, A)
        ! Two dofs per vertex of a structured P1 triangulation, row-wise vertex
        ! numbering, 6x6 element blocks L (x) C with PSD L and SPD C, scaled by
        ! a complex factor plus a complex shift: complex symmetric, nonsingular.
        integer, intent(in) :: m
        type(csc_z_t), intent(out) :: A
        integer, allocatable :: rows(:), cols(:)
        complex(dp), allocatable :: vals(:)
        type(fortsparse_status_t) :: status
        complex(dp), parameter :: shift = (0.05_dp, 0.3_dp)
        integer :: ix, iy, k, n, i, v00

        n = 2*m*m
        allocate (rows(72*(m - 1)**2 + n), cols(72*(m - 1)**2 + n), &
            vals(72*(m - 1)**2 + n))
        k = 0
        do iy = 1, m - 1
            do ix = 1, m - 1
                v00 = (iy - 1)*m + ix
                call element(rows, cols, vals, k, [v00, v00 + 1, v00 + m + 1])
                call element(rows, cols, vals, k, [v00, v00 + m, v00 + m + 1])
            end do
        end do
        do i = 1, n
            k = k + 1
            rows(k) = i
            cols(k) = i
            vals(k) = shift
        end do
        call csc_from_triplet(n, n, rows(:k), cols(:k), vals(:k), A, status)
        call check(status_ok(status), 'elasticity construction')
    end subroutine build_elasticity

    subroutine element(rows, cols, vals, k, vertex)
        integer, intent(inout) :: rows(:), cols(:), k
        complex(dp), intent(inout) :: vals(:)
        integer, intent(in) :: vertex(3)
        real(dp), parameter :: lap(3, 3) = reshape([1.1_dp, -0.9_dp, 0.1_dp, &
            -0.9_dp, 2.1_dp, -0.9_dp, 0.1_dp, -0.9_dp, 1.1_dp], [3, 3])
        real(dp), parameter :: coupling(2, 2) = reshape([2.0_dp, 0.5_dp, &
            0.5_dp, 2.0_dp], [2, 2])
        complex(dp), parameter :: scale = (1.0_dp, 0.2_dp)
        integer :: i, j, a, b
        do j = 1, 3
            do b = 1, 2
                do i = 1, 3
                    do a = 1, 2
                        k = k + 1
                        rows(k) = 2*(vertex(i) - 1) + a
                        cols(k) = 2*(vertex(j) - 1) + b
                        vals(k) = scale*lap(i, j)*coupling(a, b)
                    end do
                end do
            end do
        end do
    end subroutine element

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (condition) return
        failures = failures + 1
        write (error_unit, '(a)') 'FAIL: '//label
    end subroutine check

end program test_fortsparse_native_default_ordering
