program test_fortsparse_native_symmetric_ordering
    ! Fill-reducing orders for the native symmetric and Hermitian factors:
    ! real LDL^T, complex HPD Cholesky, Hermitian LDL^H and complex symmetric
    ! LDL^T. Oracles: exact Laplacian spectra for inertia, explicit natural
    ! factors for fill, and residuals against the original matrix.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse
    implicit none
    ! Grid shapes (x, y), numbered along x first. large: a 14 x 300 strip,
    ! n = 4200 with bandwidth 14, reaches the default AMD rule; small: 30 x 30.
    integer, parameter :: large(2) = [14, 300], small(2) = [30, 30]
    real(dp), parameter :: strong = 3.5_dp, near = 2.0137_dp
    integer :: failures

    failures = 0
    call real_default_orders()
    call real_fill_survives_pivoting()
    call real_saddle_point()
    call hpd_orders()
    call hermitian_orders()
    call complex_symmetric_orders()
    if (failures /= 0) then
        write (error_unit, '(i0,a)') failures, ' symmetric-ordering checks failed'
        stop 1
    end if
    print '(a)', 'PASS'
    stop 0

contains

    subroutine real_default_orders()
        type(csc_t) :: A
        integer(i8) :: nnz
        call laplacian(large, strong, A)
        call factor_real(A, FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_ORDER_DEFAULT, &
            FORTSPARSE_ORDER_AMD, negatives(large, strong), nnz, 'real LDLT default')
    end subroutine real_default_orders

    subroutine real_fill_survives_pivoting()
        ! Unrestricted Bunch-Kaufman interchanges undid AMD here (3-4x fill).
        type(csc_t) :: A
        integer(i8) :: definite, ordered, natural
        real(dp) :: shifts(2)
        integer :: s
        call laplacian(small, 0.0_dp, A)
        call factor_real(A, FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_ORDER_AMD, &
            FORTSPARSE_ORDER_AMD, 0, definite, 'real LDLT AMD s=0')
        shifts = [near, strong]
        do s = 1, size(shifts)
            call laplacian(small, shifts(s), A)
            call factor_real(A, FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_ORDER_AMD, &
                FORTSPARSE_ORDER_AMD, negatives(small, shifts(s)), ordered, &
                'real LDLT AMD indefinite')
            call factor_real(A, FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_ORDER_NATURAL, &
                FORTSPARSE_ORDER_DEFAULT, negatives(small, shifts(s)), natural, &
                'real LDLT natural indefinite')
            call check(real(ordered, dp) <= 1.5_dp*real(definite, dp), &
                'real LDLT: indefinite AMD fill stays near definite AMD fill')
            call check(real(ordered, dp) < 0.5_dp*real(natural, dp), &
                'real LDLT: indefinite AMD fill well below natural order')
        end do
    end subroutine real_fill_survives_pivoting

    subroutine real_saddle_point()
        ! Zero diagonal block: AMD puts constraints first, so every one of
        ! them needs a 2x2 pivot or a delay. Inertia (n_h positive, n_c negative).
        type(csc_t) :: A
        integer(i8) :: ordered, natural
        integer :: nh
        nh = 16*16
        call saddle_point([16, 16], A)
        call factor_real(A, FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_ORDER_AMD, &
            FORTSPARSE_ORDER_AMD, nh/2, ordered, 'real KKT AMD')
        call factor_real(A, FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_ORDER_NATURAL, &
            FORTSPARSE_ORDER_DEFAULT, nh/2, natural, 'real KKT natural')
        call check(real(ordered, dp) < 0.5_dp*real(natural, dp), &
            'real KKT: AMD fill well below natural order')
    end subroutine real_saddle_point

    subroutine hpd_orders()
        type(csc_z_t) :: A
        integer(i8) :: ordered, natural, nnz
        call hermitian_laplacian(large, 0.0_dp, A)
        call factor_complex(A, FORTSPARSE_SYM_SPD, FORTSPARSE_ORDER_DEFAULT, &
            FORTSPARSE_ORDER_AMD, 0, nnz, 'HPD default', .false.)
        call hermitian_laplacian(small, 0.0_dp, A)
        call factor_complex(A, FORTSPARSE_SYM_SPD, FORTSPARSE_ORDER_AMD, &
            FORTSPARSE_ORDER_AMD, 0, ordered, 'HPD AMD', .true.)
        call factor_complex(A, FORTSPARSE_SYM_SPD, FORTSPARSE_ORDER_NATURAL, &
            FORTSPARSE_ORDER_DEFAULT, 0, natural, 'HPD natural', .true.)
        call check(real(ordered, dp) < 0.5_dp*real(natural, dp), &
            'HPD: AMD fill well below natural order')
    end subroutine hpd_orders

    subroutine hermitian_orders()
        type(csc_z_t) :: A
        integer(i8) :: definite, ordered, natural, nnz
        call hermitian_laplacian(large, strong, A)
        call factor_complex(A, FORTSPARSE_SYM_HERMITIAN, FORTSPARSE_ORDER_DEFAULT, &
            FORTSPARSE_ORDER_AMD, negatives(large, strong), nnz, 'LDLH default', .false.)
        call hermitian_laplacian(small, 0.0_dp, A)
        call factor_complex(A, FORTSPARSE_SYM_HERMITIAN, FORTSPARSE_ORDER_AMD, &
            FORTSPARSE_ORDER_AMD, 0, definite, 'LDLH AMD s=0', .false.)
        call hermitian_laplacian(small, strong, A)
        call factor_complex(A, FORTSPARSE_SYM_HERMITIAN, FORTSPARSE_ORDER_AMD, &
            FORTSPARSE_ORDER_AMD, negatives(small, strong), ordered, 'LDLH AMD', .true.)
        call factor_complex(A, FORTSPARSE_SYM_HERMITIAN, FORTSPARSE_ORDER_NATURAL, &
            FORTSPARSE_ORDER_DEFAULT, negatives(small, strong), natural, &
            'LDLH natural', .false.)
        call check(real(ordered, dp) <= 1.5_dp*real(definite, dp), &
            'LDLH: indefinite AMD fill stays near definite AMD fill')
        call check(real(ordered, dp) < 0.5_dp*real(natural, dp), &
            'LDLH: indefinite AMD fill well below natural order')
    end subroutine hermitian_orders

    subroutine complex_symmetric_orders()
        ! (1+0.2i)(L - s I): complex symmetric, not Hermitian; no inertia.
        type(csc_z_t) :: A
        integer(i8) :: ordered, natural, nnz
        call symmetric_laplacian(large, strong, A)
        call factor_complex(A, FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_ORDER_DEFAULT, &
            FORTSPARSE_ORDER_AMD, -1, nnz, 'complex LDLT default', .false.)
        call symmetric_laplacian(small, strong, A)
        call factor_complex(A, FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_ORDER_AMD, &
            FORTSPARSE_ORDER_AMD, -1, ordered, 'complex LDLT AMD', .true.)
        call factor_complex(A, FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_ORDER_NATURAL, &
            FORTSPARSE_ORDER_DEFAULT, -1, natural, 'complex LDLT natural', .false.)
        call check(real(ordered, dp) < 0.5_dp*real(natural, dp), &
            'complex LDLT: AMD fill well below natural order')
    end subroutine complex_symmetric_orders

    subroutine factor_real(A, symmetry, ordering, used, negative, nnz, label)
        ! negative: exact number of negative eigenvalues.
        type(csc_t), intent(in) :: A
        integer, intent(in) :: symmetry, ordering, used, negative
        integer(i8), intent(out) :: nnz
        character(*), intent(in) :: label
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: stats
        real(dp), allocatable :: b(:), x(:)
        integer(i8) :: nneg, nzero, npos
        integer :: i, reliability

        nnz = huge(nnz)
        allocate (b(A%nrow), x(A%nrow))
        b = [(sin(0.37_dp*real(i, dp)) + 1.0_dp, i=1, A%nrow)]
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry = symmetry
        solver%options%ordering = ordering
        call sparse_factor(solver, A, status)
        call check(status_ok(status), label//': factor')
        if (.not. status_ok(status)) return
        call sparse_factor_stats(solver, stats, status)
        call check(stats%ordering_used == used, label//': reported ordering')
        nnz = stats%factor_nnz
        call sparse_inertia(solver, nneg, nzero, npos, reliability, status)
        call check(status_ok(status), label//': inertia available')
        call check(nneg == negative .and. nzero == 0_i8 .and. &
            npos == int(A%nrow - negative, i8), label//': exact inertia')
        call sparse_solve(solver, b, x, status)
        call check(status_ok(status), label//': solve')
        call check(maxval(abs(csc_matvec(A, x) - b)) <= 1.0e-12_dp* &
            (8.0_dp*maxval(abs(x)) + maxval(abs(b))), label//': original residual')
        call sparse_free(solver)
    end subroutine factor_real

    subroutine factor_complex(A, symmetry, ordering, used, negative, nnz, label, reuse)
        ! negative < 0: no inertia oracle (complex symmetric).
        type(csc_z_t), intent(in) :: A
        integer, intent(in) :: symmetry, ordering, used, negative
        integer(i8), intent(out) :: nnz
        character(*), intent(in) :: label
        logical, intent(in) :: reuse
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: stats
        complex(dp), allocatable :: b(:), x(:)
        integer(i8) :: nneg, nzero, npos
        integer :: i, reliability
        logical :: reused

        nnz = huge(nnz)
        allocate (b(A%nrow), x(A%nrow))
        b = [(cmplx(cos(0.21_dp*real(i, dp)), 1.0_dp/real(i, dp), dp), i=1, A%nrow)]
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry = symmetry
        solver%options%ordering = ordering
        call sparse_factor(solver, A, status)
        call check(status_ok(status), label//': factor')
        if (.not. status_ok(status)) return
        if (reuse) then
            ! Same pattern, new values: the retained order and analysis apply.
            call sparse_refactor(solver, A, status, reused)
            call check(status_ok(status) .and. reused, label//': same-pattern reuse')
        end if
        call sparse_factor_stats(solver, stats, status)
        call check(stats%ordering_used == used, label//': reported ordering')
        nnz = stats%factor_nnz
        if (negative >= 0) then
            call sparse_inertia(solver, nneg, nzero, npos, reliability, status)
            call check(status_ok(status), label//': inertia available')
            call check(nneg == negative .and. nzero == 0_i8 .and. &
                npos == int(A%nrow - negative, i8), label//': exact inertia')
        end if
        call sparse_solve(solver, b, x, status)
        call check(status_ok(status), label//': solve')
        call check(maxval(abs(csc_matvec(A, x) - b)) <= 1.0e-12_dp* &
            (8.0_dp*maxval(abs(x)) + maxval(abs(b))), label//': original residual')
        call sparse_free(solver)
    end subroutine factor_complex

    integer function negatives(g, shift) result(count)
        ! Eigenvalues of the Dirichlet five-point Laplacian minus shift.
        integer, intent(in) :: g(2)
        real(dp), intent(in) :: shift
        real(dp) :: pi, hx, hy
        integer :: i, j
        pi = acos(-1.0_dp)
        hx = pi/real(g(1) + 1, dp)
        hy = pi/real(g(2) + 1, dp)
        count = 0
        do j = 1, g(2)
            do i = 1, g(1)
                if (4.0_dp - 2.0_dp*cos(i*hx) - 2.0_dp*cos(j*hy) < shift) &
                    count = count + 1
            end do
        end do
    end function negatives

    subroutine stencil(g, rows, cols, k)
        ! Five-point stencil numbered along x first; bandwidth g(1).
        integer, intent(in) :: g(2)
        integer, allocatable, intent(out) :: rows(:), cols(:)
        integer, intent(out) :: k
        integer :: ix, iy, c
        allocate (rows(5*g(1)*g(2)), cols(5*g(1)*g(2)))
        k = 0
        do iy = 1, g(2)
            do ix = 1, g(1)
                c = (iy - 1)*g(1) + ix
                call push(rows, cols, k, c, c)
                if (ix > 1) call push(rows, cols, k, c, c - 1)
                if (ix < g(1)) call push(rows, cols, k, c, c + 1)
                if (iy > 1) call push(rows, cols, k, c, c - g(1))
                if (iy < g(2)) call push(rows, cols, k, c, c + g(1))
            end do
        end do
    end subroutine stencil

    subroutine push(rows, cols, k, i, j)
        integer, intent(inout) :: rows(:), cols(:), k
        integer, intent(in) :: i, j
        k = k + 1
        rows(k) = i
        cols(k) = j
    end subroutine push

    subroutine laplacian(g, shift, A)
        integer, intent(in) :: g(2)
        real(dp), intent(in) :: shift
        type(csc_t), intent(out) :: A
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        type(fortsparse_status_t) :: status
        integer :: k
        call stencil(g, rows, cols, k)
        allocate (vals(k))
        where (rows(:k) == cols(:k))
            vals = 4.0_dp - shift
        elsewhere
            vals = -1.0_dp
        end where
        call csc_from_triplet(g(1)*g(2), g(1)*g(2), rows(:k), cols(:k), vals, A, status)
        call check(status_ok(status), 'Laplacian construction')
    end subroutine laplacian

    subroutine hermitian_laplacian(g, shift, A)
        ! Unitary gauge D (L - s I) D^H, D = diag(exp(0.3 i v)): Hermitian
        ! with complex off-diagonals and the spectrum of the real Laplacian.
        integer, intent(in) :: g(2)
        real(dp), intent(in) :: shift
        type(csc_z_t), intent(out) :: A
        integer, allocatable :: rows(:), cols(:)
        complex(dp), allocatable :: vals(:)
        type(fortsparse_status_t) :: status
        integer :: k, p
        call stencil(g, rows, cols, k)
        allocate (vals(k))
        do p = 1, k
            if (rows(p) == cols(p)) then
                vals(p) = cmplx(4.0_dp - shift, 0.0_dp, dp)
            else
                vals(p) = -exp(cmplx(0.0_dp, 0.3_dp*real(rows(p) - cols(p), dp), dp))
            end if
        end do
        call csc_from_triplet(g(1)*g(2), g(1)*g(2), rows(:k), cols(:k), vals, A, status)
        call check(status_ok(status), 'Hermitian Laplacian construction')
    end subroutine hermitian_laplacian

    subroutine symmetric_laplacian(g, shift, A)
        integer, intent(in) :: g(2)
        real(dp), intent(in) :: shift
        type(csc_z_t), intent(out) :: A
        integer, allocatable :: rows(:), cols(:)
        complex(dp), allocatable :: vals(:)
        type(fortsparse_status_t) :: status
        integer :: k
        call stencil(g, rows, cols, k)
        allocate (vals(k))
        where (rows(:k) == cols(:k))
            vals = (1.0_dp, 0.2_dp)*(4.0_dp - shift)
        elsewhere
            vals = -(1.0_dp, 0.2_dp)
        end where
        call csc_from_triplet(g(1)*g(2), g(1)*g(2), rows(:k), cols(:k), vals, A, status)
        call check(status_ok(status), 'complex symmetric Laplacian construction')
    end subroutine symmetric_laplacian

    subroutine saddle_point(g, A)
        ! [[H, B^T], [B, 0]], H the grid Laplacian with n_h vertices, B couples
        ! vertices 2i-1 and 2i with (1, -2): full row rank, so the inertia is
        ! (n_h positive, n_h/2 negative).
        integer, intent(in) :: g(2)
        type(csc_t), intent(out) :: A
        integer, allocatable :: rows(:), cols(:), r(:), c(:)
        real(dp), allocatable :: vals(:)
        type(fortsparse_status_t) :: status
        integer :: k, nh, i
        call stencil(g, rows, cols, k)
        nh = g(1)*g(2)
        allocate (r(k + 2*nh), c(k + 2*nh), vals(k + 2*nh))
        r(:k) = rows(:k)
        c(:k) = cols(:k)
        where (rows(:k) == cols(:k))
            vals(:k) = 4.0_dp
        elsewhere
            vals(:k) = -1.0_dp
        end where
        do i = 1, nh/2
            r(k + 1:k + 4) = [nh + i, 2*i - 1, nh + i, 2*i]
            c(k + 1:k + 4) = [2*i - 1, nh + i, 2*i, nh + i]
            vals(k + 1:k + 4) = [1.0_dp, 1.0_dp, -2.0_dp, -2.0_dp]
            k = k + 4
        end do
        call csc_from_triplet(nh + nh/2, nh + nh/2, r(:k), c(:k), vals(:k), A, status)
        call check(status_ok(status), 'saddle-point construction')
    end subroutine saddle_point

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (condition) return
        failures = failures + 1
        write (error_unit, '(a)') 'FAIL: '//label
    end subroutine check

end program test_fortsparse_native_symmetric_ordering
