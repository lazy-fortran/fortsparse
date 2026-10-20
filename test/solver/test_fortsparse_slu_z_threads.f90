program test_fortsparse_slu_z_threads
    ! OpenMP subtree-task complex supernodal LU against the serial factor of
    ! the same analysis: original-matrix backward errors for A, A^T and A^H,
    ! identical perturbation and exchange counts and gate decisions, growth
    ! and minimum pivot within rounding, the actual team, the public team
    ! report, and (argument "repeat") repeated factor/refactor reuse. The
    ! matrices have no near-threshold or tied pivots, so only the update
    ! summation order may differ between teams.
    use fortsparse
    use fortsparse_native_slu, only: slu_factor_t
    implicit none
    integer :: seed
    integer, allocatable :: ptr(:), rows(:)
    complex(dp), allocatable :: values(:)
    logical :: threaded
    character(16) :: mode

    seed = 5150
    threaded = .false.
    call get_command_argument(1, mode)
    call grid2d(64, ptr, rows, values)
    if (trim(mode) == 'repeat') then
        call repeated(ptr, rows, values)
        print '(a)', 'complex supernodal LU threads: repeated team factorizations verified'
        stop
    end if
    call compare_teams(ptr, rows, values)
    call public_team(ptr, rows, values)
    if (threaded) then
        print '(a)', 'complex supernodal LU threads: parallel teams verified'
    else
        print '(a)', 'complex supernodal LU threads: library without OpenMP, serial only'
    end if
contains
    real(dp) function uniform()
        seed = int(modulo(1103515245_i8*seed+12345_i8, 2147483648_i8))
        uniform = real(seed, dp)/2147483648.0_dp
    end function uniform

    subroutine grid2d(m, ptr, rows, values)
        ! 9-point convection-diffusion stencil with complex, unsymmetric
        ! off-diagonal values and a weak, column-varying complex diagonal.
        integer, intent(in) :: m
        integer, allocatable, intent(out) :: ptr(:), rows(:)
        complex(dp), allocatable, intent(out) :: values(:)
        integer :: n, ix, iy, dx, dy, j, i, k
        n = m*m
        allocate(ptr(n+1), rows(9*n), values(9*n))
        k = 0
        ptr(1) = 1
        do iy = 1, m
            do ix = 1, m
                j = ix+(iy-1)*m
                do dy = -1, 1
                    do dx = -1, 1
                        if (min(ix+dx, iy+dy) < 1) cycle
                        if (max(ix+dx, iy+dy) > m) cycle
                        i = ix+dx+(iy+dy-1)*m
                        k = k+1
                        rows(k) = i
                        if (i == j) then
                            values(k) = cmplx(2.0_dp+uniform(), 1.0_dp-uniform(), dp) &
                                *(1.0_dp+mod(j, 7))
                        else
                            values(k) = cmplx(-1.0_dp+0.6_dp*dx-0.3_dp*dy+0.5_dp*uniform(), &
                                0.4_dp*dy-0.2_dp*dx+0.3_dp*uniform(), dp)
                        end if
                    end do
                end do
                ptr(j+1) = k+1
            end do
        end do
        rows = rows(:k)
        values = values(:k)
    end subroutine grid2d

    subroutine make_rhs(ptr, rows, values, known, b, op)
        ! b = op(A) known for op 0 (A), 1 (A^T) or 2 (A^H).
        integer, intent(in) :: ptr(:), rows(:), op
        complex(dp), intent(in) :: values(:)
        complex(dp), intent(out) :: known(:), b(:)
        integer :: j, p
        complex(dp) :: v
        known = [(cmplx(1.0_dp+mod(j,5), 0.5_dp*mod(j,3), dp), j=1, size(known))]
        b = (0.0_dp, 0.0_dp)
        do j = 1, size(ptr)-1
            do p = ptr(j), ptr(j+1)-1
                v = values(p)
                if (op == 2) v = conjg(v)
                if (op == 0) then
                    b(rows(p)) = b(rows(p))+v*known(j)
                else
                    b(j) = b(j)+v*known(rows(p))
                end if
            end do
        end do
    end subroutine make_rhs

    real(dp) function backward_error(ptr, rows, values, x, b, op) result(error)
        ! Normwise backward error with max(row, column) sums of |A|, an upper
        ! bound of ||op(A)||_inf for every operation.
        integer, intent(in) :: ptr(:), rows(:), op
        complex(dp), intent(in) :: values(:), x(:), b(:)
        complex(dp) :: ax(size(x))
        real(dp) :: rowsum(size(x)), norm
        integer :: j, p
        call make_product(ptr, rows, values, x, ax, op)
        rowsum = 0.0_dp
        norm = 0.0_dp
        do j = 1, size(ptr)-1
            norm = max(norm, sum(abs(values(ptr(j):ptr(j+1)-1))))
            do p = ptr(j), ptr(j+1)-1
                rowsum(rows(p)) = rowsum(rows(p))+abs(values(p))
            end do
        end do
        norm = max(norm, maxval(rowsum))
        error = maxval(abs(b-ax))/(norm*maxval(abs(x))+maxval(abs(b)))
    end function backward_error

    subroutine make_product(ptr, rows, values, x, y, op)
        integer, intent(in) :: ptr(:), rows(:), op
        complex(dp), intent(in) :: values(:), x(:)
        complex(dp), intent(out) :: y(:)
        integer :: j, p
        complex(dp) :: v
        y = (0.0_dp, 0.0_dp)
        do j = 1, size(ptr)-1
            do p = ptr(j), ptr(j+1)-1
                v = values(p)
                if (op == 2) v = conjg(v)
                if (op == 0) then
                    y(rows(p)) = y(rows(p))+v*x(j)
                else
                    y(j) = y(j)+v*x(rows(p))
                end if
            end do
        end do
    end subroutine make_product

    subroutine compare_teams(ptr, rows, values)
        integer, intent(in) :: ptr(:), rows(:)
        complex(dp), intent(in) :: values(:)
        type(slu_factor_t) :: serial, team
        type(fortsparse_status_t) :: status
        complex(dp), allocatable :: known(:), b(:), x(:)
        real(dp) :: be_serial(0:2), be
        integer :: n, t, op, threads(3)
        n = size(ptr)-1
        allocate(known(n), b(n), x(n))
        call serial%analyze(n, ptr, rows, values, status, nthreads=1)
        call check(status%code == FORTSPARSE_OK, 'serial analysis')
        call serial%factor(ptr, rows, values, status)
        call check(status%code == FORTSPARSE_OK, 'serial factor')
        call check(serial%perturbed_pivots == 0, 'no perturbed pivots')
        call check(serial%block_swaps > 0, 'in-block exchanges exercised')
        do op = 0, 2
            call make_rhs(ptr, rows, values, known, b, op)
            x = b
            call serial%solve(x, status, op)
            be_serial(op) = backward_error(ptr, rows, values, x, b, op)
            call check(be_serial(op) < 1.0e-13_dp, 'serial backward error')
        end do
        threads = [2, 4, 8]
        do t = 1, size(threads)
            call team%analyze(n, ptr, rows, values, status, nthreads=threads(t))
            call check(status%code == FORTSPARSE_OK, 'team analysis')
            if (team%requested_threads > 1) threaded = .true.
            call check(team%task_groups > 0, 'task groups')
            call check(team%plan_bytes == team%memory_bytes(), 'exact plan')
            call team%factor(ptr, rows, values, status)
            call check(status%code == FORTSPARSE_OK, 'team factor')
            call check(team%active_threads >= 1 .and. team%active_threads <= threads(t), &
                'actual team')
            call check(team%perturbed_pivots == serial%perturbed_pivots, &
                'identical perturbation count')
            call check(team%block_swaps == serial%block_swaps, 'identical exchange count')
            call check(team%gate_reason() == serial%gate_reason(), 'identical gate decision')
            call check(abs(team%max_growth-serial%max_growth) &
                <= 1.0e-10_dp*serial%max_growth, 'growth')
            call check(abs(team%min_pivot-serial%min_pivot) &
                <= 1.0e-10_dp*serial%min_pivot, 'minimum pivot')
            do op = 0, 2
                call make_rhs(ptr, rows, values, known, b, op)
                x = b
                call team%solve(x, status, op)
                be = backward_error(ptr, rows, values, x, b, op)
                call check(be < max(4.0_dp*be_serial(op), 1.0e-15_dp), 'team backward error')
            end do
            print '(a,i0,a,i0,a,i0)', 'requested=', threads(t), ' active=', &
                team%active_threads, ' groups=', team%task_groups
            call team%free()
        end do
        call serial%free()
    end subroutine compare_teams

    subroutine public_team(ptr, rows, values)
        integer, intent(in) :: ptr(:), rows(:)
        complex(dp), intent(in) :: values(:)
        type(sparse_solver_t) :: solver
        type(sparse_factor_stats_t) :: stats
        type(fortsparse_status_t) :: status
        type(csc_z_t) :: matrix
        integer, allocatable :: cols(:)
        complex(dp), allocatable :: known(:), b(:), x(:)
        integer :: n, j
        n = size(ptr)-1
        allocate(cols(size(rows)), known(n), b(n), x(n))
        do j = 1, n
            cols(ptr(j):ptr(j+1)-1) = j
        end do
        call csc_from_triplet(n, n, rows, cols, values, matrix, status)
        call check(status_ok(status), 'public: CSC construction')
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%lu_method = FORTSPARSE_LU_SUPERNODAL
        solver%options%ordering = FORTSPARSE_ORDER_AMD
        solver%options%nthreads = 4
        call sparse_factor(solver, matrix, status)
        call check(status_ok(status), 'public: factor '//trim(status%msg))
        call sparse_factor_stats(solver, stats, status)
        call check(stats%native_factor_family == 'LU_SUPERNODAL' .and. &
            .not. stats%lu_fallback, 'public: supernodal factor kept')
        call check(stats%nthreads == 4, 'public: requested team')
        call check(stats%active_threads >= 1 .and. stats%active_threads <= 4, &
            'public: active team')
        if (threaded) call check(stats%active_threads > 1, 'public: parallel team')
        call make_rhs(ptr, rows, values, known, b, 2)
        call sparse_solve_op(solver, b, x, FORTSPARSE_OP_ADJOINT, status)
        call check(status_ok(status), 'public: adjoint solve')
        call check(maxval(abs(x-known)) <= 1.0e-10_dp*maxval(abs(known)), &
            'public: adjoint known solution')
        call sparse_destroy(solver)
    end subroutine public_team

    subroutine repeated(ptr, rows, values)
        ! Repeated factor/refactor with the eight-thread team: repeats of the
        ! same values reproduce the first threaded solution; changed values
        ! keep the backward error (detects task races).
        integer, intent(in) :: ptr(:), rows(:)
        complex(dp), intent(in) :: values(:)
        type(slu_factor_t) :: lu
        type(fortsparse_status_t) :: status
        complex(dp), allocatable :: known(:), b(:), x(:), first(:), scaled(:)
        integer :: n, k
        n = size(ptr)-1
        allocate(known(n), b(n), x(n), first(n))
        call make_rhs(ptr, rows, values, known, b, 0)
        call lu%analyze(n, ptr, rows, values, status, nthreads=8)
        call check(status%code == FORTSPARSE_OK, 'repeat: analysis')
        do k = 1, 30
            if (mod(k, 2) == 0) then
                scaled = values*cmplx(1.0_dp+0.01_dp*k, 0.005_dp*k, dp)
                call lu%factor(ptr, rows, scaled, status)
                call check(status%code == FORTSPARSE_OK, 'repeat: refactor')
                x = b
                call lu%solve(x, status, 0)
                call check(backward_error(ptr, rows, scaled, x, b, 0) < 1.0e-13_dp, &
                    'repeat: refactor backward error')
            else
                call lu%factor(ptr, rows, values, status)
                call check(status%code == FORTSPARSE_OK, 'repeat: factor')
                x = b
                call lu%solve(x, status, 0)
                if (k == 1) first = x
                call check(maxval(abs(x-first)) <= 1.0e-13_dp*maxval(abs(first)), &
                    'repeat: reproduced solution')
            end if
        end do
        call lu%free()
    end subroutine repeated

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (condition) return
        print '(a)', 'FAIL: '//message
        error stop 1
    end subroutine check
end program test_fortsparse_slu_z_threads
