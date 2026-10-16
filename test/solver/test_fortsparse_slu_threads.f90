program test_fortsparse_slu_threads
    ! OpenMP subtree-task supernodal LU against the serial factor of the same
    ! analysis: original-matrix backward errors of known solutions, identical
    ! perturbation and exchange counts and gate decisions, growth and minimum
    ! pivot within rounding, the actual team size, cap-reduced teams,
    ! LOW_MEMORY (one thread) and repeated factor/refactor reuse.
    ! The generated matrices have no near-threshold pivots (pivots are far
    ! above tau*||M|| and in-block maxima are not tied), so the counts must
    ! be identical across team sizes; only update summation order may differ.
    ! Without OpenMP in the library the requested team stays one and the
    ! team-size checks reduce to the serial ones.
    use fortsparse
    use fortsparse_memory, only: memory_budget_t
    use fortsparse_native_slu, only: slu_factor_t
    implicit none
    integer :: seed
    integer, allocatable :: ptr(:), rows(:)
    real(dp), allocatable :: values(:)
    logical :: threaded
    character(16) :: mode

    seed = 4243
    threaded = .false.
    ! "repeat" (a separate ctest entry, so each run stays below one second)
    ! runs only the repeated eight-thread factor/refactor check.
    call get_command_argument(1, mode)
    if (trim(mode) == 'repeat') then
        call grid2d(100, ptr, rows, values)
        call repeated(ptr, rows, values)
        print '(a)', 'supernodal LU threads: repeated team factorizations verified'
        stop
    end if
    call kin_like(9, ptr, rows, values)
    call compare_teams('3D 27-point, 2 coupled dofs', ptr, rows, values)
    call low_memory_team(ptr, rows, values)
    call grid2d(100, ptr, rows, values)
    call compare_teams('2D convection 9-point', ptr, rows, values)
    call capped_team(ptr, rows, values)
    call public_team(ptr, rows, values)
    call disconnected(ptr, rows, values)
    call compare_teams('disconnected 2D blocks', ptr, rows, values)
    if (threaded) then
        print '(a)', 'supernodal LU threads: parallel teams verified'
    else
        print '(a)', 'supernodal LU threads: library without OpenMP, serial checks only'
    end if
contains
    real(dp) function uniform()
        seed = int(modulo(1103515245_i8*seed+12345_i8, 2147483648_i8))
        uniform = real(seed, dp)/2147483648.0_dp
    end function uniform

    subroutine compare_teams(label, ptr, rows, values)
        character(*), intent(in) :: label
        integer, intent(in) :: ptr(:), rows(:)
        real(dp), intent(in) :: values(:)
        type(slu_factor_t) :: serial, team
        type(fortsparse_status_t) :: status
        real(dp), allocatable :: known(:), b(:), x(:)
        real(dp) :: be_serial, be
        integer :: n, t, threads(3)
        n = size(ptr)-1
        allocate(known(n), b(n), x(n))
        call make_rhs(ptr, rows, values, known, b)
        call serial%analyze(n, ptr, rows, values, status, nthreads=1)
        call check(status%code == FORTSPARSE_OK, label//': serial analysis')
        call check(serial%active_threads == 1, label//': serial team')
        call serial%factor(ptr, rows, values, status)
        call check(status%code == FORTSPARSE_OK, label//': serial factor')
        call check(serial%perturbed_pivots == 0, label//': no perturbed pivots')
        call check(serial%block_swaps > 0, label//': in-block exchanges exercised')
        x = b
        call serial%solve(x, status, 0)
        be_serial = backward_error(ptr, rows, values, x, b)
        print '(a,a,es9.2,a,i0,a,es9.2)', label, ': serial backward error=', be_serial, &
            ' swaps=', serial%block_swaps, ' growth=', serial%max_growth
        call check(be_serial < 1.0e-13_dp, label//': serial backward error')
        threads = [2, 4, 8]
        do t = 1, size(threads)
            call team%analyze(n, ptr, rows, values, status, nthreads=threads(t))
            call check(status%code == FORTSPARSE_OK, label//': team analysis')
            if (team%requested_threads > 1) threaded = .true.
            call check(team%requested_threads == threads(t) .or. &
                team%requested_threads == 1, label//': requested team')
            call check(team%active_threads >= 1 .and. &
                team%active_threads <= team%requested_threads, label//': planned team')
            call check(team%task_groups > 0, label//': task groups')
            call check(team%plan_bytes == team%memory_bytes(), label//': exact plan')
            call team%factor(ptr, rows, values, status)
            call check(status%code == FORTSPARSE_OK, label//': team factor')
            call check(team%active_threads <= threads(t), label//': actual team')
            call check(team%perturbed_pivots == serial%perturbed_pivots, &
                label//': identical perturbation count')
            call check(team%block_swaps == serial%block_swaps, &
                label//': identical exchange count')
            call check(team%gate_reason() == serial%gate_reason(), &
                label//': identical gate decision')
            call check(abs(team%max_growth-serial%max_growth) &
                <= 1.0e-10_dp*serial%max_growth, label//': growth')
            call check(abs(team%min_pivot-serial%min_pivot) &
                <= 1.0e-10_dp*serial%min_pivot, label//': minimum pivot')
            x = b
            call team%solve(x, status, 0)
            be = backward_error(ptr, rows, values, x, b)
            call check(be < max(4.0_dp*be_serial, 1.0e-15_dp), label//': backward error')
            x = b
            call team%solve(x, status, 1)
            call check(status%code == FORTSPARSE_OK, label//': transpose solve')
            print '(a,a,i0,a,i0,a,i0,a,es9.2)', label, ': requested=', threads(t), &
                ' active=', team%active_threads, ' groups=', team%task_groups, &
                ' backward error=', be
            call team%free()
        end do
        call serial%free()
    end subroutine compare_teams

    subroutine capped_team(ptr, rows, values)
        ! A cap between the one- and eight-thread layouts halves the team.
        integer, intent(in) :: ptr(:), rows(:)
        real(dp), intent(in) :: values(:)
        type(slu_factor_t) :: one, eight, capped
        type(memory_budget_t) :: budget
        type(fortsparse_status_t) :: status
        real(dp), allocatable :: known(:), b(:), x(:)
        integer(i8) :: per_thread, one_plan, one_peak
        integer :: n
        n = size(ptr)-1
        budget = memory_budget_t()
        budget%limit = huge(1_i8)/4
        call one%analyze(n, ptr, rows, values, status, nthreads=1, budget=budget)
        call check(status%code == FORTSPARSE_OK, 'capped: serial strict analysis')
        one_plan = one%plan_bytes
        one_peak = one%peak_memory_bytes
        call one%free(budget)
        budget = memory_budget_t()
        budget%limit = huge(1_i8)/4
        call eight%analyze(n, ptr, rows, values, status, nthreads=8, budget=budget)
        call check(status%code == FORTSPARSE_OK, 'capped: team strict analysis')
        call check(.not. eight%budget_reduced_threads, 'capped: uncapped team kept')
        if (eight%requested_threads == 1) then
            call eight%free(budget)
            return
        end if
        call check(eight%active_threads >= 4, 'capped: uncapped team of at least four')
        ! Per-thread row map, update scratch and statistics (the schedule
        ! headers also grow with the finer eight-thread grouping, so this
        ! overestimates a thread slightly).
        per_thread = (eight%plan_bytes-one_plan)/(eight%active_threads-1)
        call eight%free(budget)
        budget = memory_budget_t()
        budget%limit = max(one_peak, one_plan+2*per_thread)
        call capped%analyze(n, ptr, rows, values, status, nthreads=8, budget=budget)
        call check(status%code == FORTSPARSE_OK, 'capped: analysis under the cap')
        call check(capped%budget_reduced_threads, 'capped: team reduced by the cap')
        call check(capped%active_threads < 8, 'capped: smaller team')
        call check(budget%peak <= budget%limit, 'capped: ledger within the cap')
        call capped%factor(ptr, rows, values, status)
        call check(status%code == FORTSPARSE_OK, 'capped: factor')
        allocate(known(n), b(n), x(n))
        call make_rhs(ptr, rows, values, known, b)
        x = b
        call capped%solve(x, status, 0)
        call check(backward_error(ptr, rows, values, x, b) < 1.0e-13_dp, &
            'capped: backward error')
        print '(a,i0,a,i0)', 'capped: active=', capped%active_threads, &
            ' cap bytes=', budget%limit
        call capped%free(budget)
        call check(budget%current == 0_i8, 'capped: ledger released')
    end subroutine capped_team

    subroutine low_memory_team(ptr, rows, values)
        integer, intent(in) :: ptr(:), rows(:)
        real(dp), intent(in) :: values(:)
        type(slu_factor_t) :: lu
        type(fortsparse_status_t) :: status
        call lu%analyze(size(ptr)-1, ptr, rows, values, status, low_memory=.true., &
            nthreads=8)
        call check(status%code == FORTSPARSE_OK, 'low memory: analysis')
        call check(lu%requested_threads == 1, 'low memory: one thread requested')
        call check(lu%active_threads == 1, 'low memory: one active thread')
        call lu%factor(ptr, rows, values, status)
        call check(status%code == FORTSPARSE_OK .and. lu%active_threads == 1, &
            'low memory: serial factor')
        call lu%free()
    end subroutine low_memory_team

    subroutine repeated(ptr, rows, values)
        ! Repeated factor/refactor with the eight-thread team: every repeat
        ! of the same values reproduces the first threaded solution, and
        ! changed values keep the backward error (detects task races).
        integer, intent(in) :: ptr(:), rows(:)
        real(dp), intent(in) :: values(:)
        type(slu_factor_t) :: lu
        type(fortsparse_status_t) :: status
        real(dp), allocatable :: known(:), b(:), x(:), first(:), scaled(:)
        integer :: n, k
        n = size(ptr)-1
        allocate(known(n), b(n), x(n), first(n))
        call make_rhs(ptr, rows, values, known, b)
        call lu%analyze(n, ptr, rows, values, status, nthreads=8)
        call check(status%code == FORTSPARSE_OK, 'repeat: analysis')
        do k = 1, 50
            if (mod(k, 2) == 0) then
                scaled = values*(1.0_dp+0.01_dp*k)
                call lu%factor(ptr, rows, scaled, status)
                call check(status%code == FORTSPARSE_OK, 'repeat: refactor')
                x = b
                call lu%solve(x, status, 0)
                call check(backward_error(ptr, rows, scaled, x, b) < 1.0e-13_dp, &
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

    subroutine public_team(ptr, rows, values)
        ! The general adapter passes options%nthreads and reports the team.
        integer, intent(in) :: ptr(:), rows(:)
        real(dp), intent(in) :: values(:)
        type(sparse_solver_t) :: solver
        type(sparse_factor_stats_t) :: stats
        type(fortsparse_status_t) :: status
        type(csc_t) :: matrix
        integer, allocatable :: cols(:)
        real(dp), allocatable :: known(:), b(:), x(:)
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
        call check(stats%lu_method_used == FORTSPARSE_LU_SUPERNODAL .and. &
            .not. stats%lu_fallback, 'public: supernodal factor kept')
        call check(stats%nthreads == 4, 'public: requested team')
        call check(stats%active_threads >= 1 .and. stats%active_threads <= 4, &
            'public: active team')
        if (threaded) call check(stats%active_threads > 1, 'public: parallel team')
        call make_rhs(ptr, rows, values, known, b)
        call sparse_solve(solver, b, x, status)
        call check(status_ok(status), 'public: solve')
        call check(backward_error(ptr, rows, values, x, b) < 1.0e-14_dp, &
            'public: backward error')
        call sparse_free(solver)
    end subroutine public_team

    subroutine make_rhs(ptr, rows, values, known, b)
        integer, intent(in) :: ptr(:), rows(:)
        real(dp), intent(in) :: values(:)
        real(dp), intent(out) :: known(:), b(:)
        integer :: j, p
        do j = 1, size(known)
            known(j) = sin(real(j, dp))+0.25_dp
        end do
        b = 0.0_dp
        do j = 1, size(known)
            do p = ptr(j), ptr(j+1)-1
                b(rows(p)) = b(rows(p))+values(p)*known(j)
            end do
        end do
    end subroutine make_rhs

    real(dp) function backward_error(ptr, rows, values, x, b) result(eta)
        ! ||b - A x||_inf / (||A||_inf ||x||_inf + ||b||_inf), original A.
        integer, intent(in) :: ptr(:), rows(:)
        real(dp), intent(in) :: values(:), x(:), b(:)
        real(dp), allocatable :: r(:), rowsum(:)
        integer :: j, p
        r = b
        allocate(rowsum(size(b)))
        rowsum = 0.0_dp
        do j = 1, size(x)
            do p = ptr(j), ptr(j+1)-1
                r(rows(p)) = r(rows(p))-values(p)*x(j)
                rowsum(rows(p)) = rowsum(rows(p))+abs(values(p))
            end do
        end do
        eta = maxval(abs(r))/(maxval(rowsum)*maxval(abs(x))+maxval(abs(b)))
    end function backward_error

    subroutine grid2d(m, ptr, rows, values)
        ! 9-point convection-diffusion stencil with unsymmetric off-diagonal
        ! values and a weak, column-varying diagonal: the scaled static
        ! pivots need a few in-block exchanges.
        integer, intent(in) :: m
        integer, allocatable, intent(out) :: ptr(:), rows(:)
        real(dp), allocatable, intent(out) :: values(:)
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
                        if (ix+dx < 1 .or. ix+dx > m .or. iy+dy < 1 .or. iy+dy > m) cycle
                        i = ix+dx+(iy+dy-1)*m
                        k = k+1
                        rows(k) = i
                        if (i == j) then
                            values(k) = (2.0_dp+uniform())*(1.0_dp+mod(j, 7))
                        else
                            values(k) = -1.0_dp+0.6_dp*dx-0.3_dp*dy+0.5_dp*uniform()
                        end if
                    end do
                end do
                ptr(j+1) = k+1
            end do
        end do
        rows = rows(:k)
        values = values(:k)
    end subroutine grid2d

    subroutine kin_like(m, ptr, rows, values)
        ! FE-style trilinear (27-point) element coupling on an m^3 grid with
        ! two coupled unknowns per node and an unsymmetric (transport-like)
        ! coupling block, as in the KIN6D consumer pattern.
        integer, intent(in) :: m
        integer, allocatable, intent(out) :: ptr(:), rows(:)
        real(dp), allocatable, intent(out) :: values(:)
        integer :: n, ix, iy, iz, dx, dy, dz, node, other, a, c, k, j
        n = 2*m**3
        allocate(ptr(n+1), rows(54*n), values(54*n))
        k = 0
        ptr(1) = 1
        do iz = 1, m
            do iy = 1, m
                do ix = 1, m
                    node = ix+(iy-1)*m+(iz-1)*m*m
                    do c = 1, 2
                        j = 2*(node-1)+c
                        do dz = -1, 1
                            do dy = -1, 1
                                do dx = -1, 1
                                    if (min(ix+dx, iy+dy, iz+dz) < 1) cycle
                                    if (max(ix+dx, iy+dy, iz+dz) > m) cycle
                                    other = ix+dx+(iy+dy-1)*m+(iz+dz-1)*m*m
                                    do a = 1, 2
                                        k = k+1
                                        rows(k) = 2*(other-1)+a
                                        if (other == node .and. a == c) then
                                            values(k) = 3.0_dp+uniform()
                                        else if (other == node) then
                                            values(k) = 3.0_dp*(a-c)+uniform()
                                        else
                                            values(k) = -0.5_dp+0.4_dp*(dx+2*dy-dz)*(3-2*c) &
                                                +0.2_dp*uniform()
                                        end if
                                    end do
                                end do
                            end do
                        end do
                        ptr(j+1) = k+1
                    end do
                end do
            end do
        end do
        rows = rows(:k)
        values = values(:k)
    end subroutine kin_like

    subroutine disconnected(ptr, rows, values)
        ! Three independent 2D blocks: separate elimination forests.
        integer, allocatable, intent(out) :: ptr(:), rows(:)
        real(dp), allocatable, intent(out) :: values(:)
        integer, allocatable :: bp(:), br(:)
        real(dp), allocatable :: bv(:)
        integer :: b, nb, nnzb
        call grid2d(36, bp, br, bv)
        nb = size(bp)-1
        nnzb = size(br)
        allocate(ptr(3*nb+1), rows(3*nnzb), values(3*nnzb))
        ptr(1) = 1
        do b = 0, 2
            ptr(b*nb+2:(b+1)*nb+1) = bp(2:)+b*nnzb
            rows(b*nnzb+1:(b+1)*nnzb) = br+b*nb
            values(b*nnzb+1:(b+1)*nnzb) = bv*(1.0_dp+b)
        end do
    end subroutine disconnected

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (condition) return
        print '(a)', 'FAIL: '//message
        error stop 1
    end subroutine check
end program test_fortsparse_slu_threads
