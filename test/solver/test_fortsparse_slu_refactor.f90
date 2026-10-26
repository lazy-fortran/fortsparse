program test_fortsparse_slu_refactor
    ! Supernodal LU refactor. A same-pattern refactor assembles through
    ! the analysis' assembly map, any other pattern searches the panel rows:
    ! equal values given once with explicit zeros on the analysed pattern
    ! (map) and once without them (search) give bitwise identical factors.
    ! A pattern outside the analysed structure analyzes again and equals a
    ! fresh factor bitwise. A changed nthreads request is honoured by the
    ! refactor: nthreads = 1 runs one thread and equals a fresh factor, and
    ! the caller's OpenMP thread count survives the one-thread factor.
    use fortsparse
    !$ use omp_lib, only: omp_get_max_threads
    implicit none
    integer, parameter :: g = 12
    integer :: team

    call check_patterns(1, team)
    call check_patterns(2, team)
    !$ call check(team == 2, 'two-thread team planned')
    call check_threads()
    print '(a)', 'supernodal LU refactor passed'
contains
    subroutine grid_matrix(scale, variant, matrix)
        ! Nonsymmetric 7-point operator on a g^3 grid. variant 1 drops a few
        ! off-diagonal pairs, variant 3 stores them as explicit zeros (the
        ! base pattern with variant-1 values), variant 2 adds long-range
        ! couplings.
        real(dp), intent(in) :: scale
        integer, intent(in) :: variant
        type(csc_t), intent(out) :: matrix
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: values(:)
        type(fortsparse_status_t) :: status
        integer :: n, i, j, k, p, q, d, nz
        integer, parameter :: step(3) = [1, g, g*g]
        n = g**3
        allocate(rows(7*n+8), cols(7*n+8), values(7*n+8))
        nz = 0
        do k = 0, g-1
            do j = 0, g-1
                do i = 0, g-1
                    p = 1+i+g*j+g*g*k
                    call add(p, p, scale*(6.5_dp+0.3_dp*sin(real(p,dp))), &
                        rows, cols, values, nz)
                    do d = 1, 3
                        q = p+step(d)
                        if (d == 1 .and. i == g-1) cycle
                        if (d == 2 .and. j == g-1) cycle
                        if (d == 3 .and. k == g-1) cycle
                        if (variant == 1 .and. mod(p, 97) == 5) cycle
                        if (variant == 3 .and. mod(p, 97) == 5) then
                            call add(q, p, 0.0_dp, rows, cols, values, nz)
                            call add(p, q, 0.0_dp, rows, cols, values, nz)
                            cycle
                        end if
                        call add(q, p, scale*(-1.0_dp+0.4_dp*cos(real(p+d,dp))), &
                            rows, cols, values, nz)
                        call add(p, q, scale*(-1.0_dp-0.3_dp*sin(real(q*d,dp))), &
                            rows, cols, values, nz)
                    end do
                end do
            end do
        end do
        if (variant == 2) then
            call add(1, n, 0.25_dp*scale, &
                rows, cols, values, nz)
            call add(n, 1, -0.5_dp*scale, &
                rows, cols, values, nz)
            call add(n/2, 3, 0.125_dp*scale, &
                rows, cols, values, nz)
            call add(7, n-5, -0.25_dp*scale, &
                rows, cols, values, nz)
        end if
        call csc_from_triplet(n, n, rows(1:nz), cols(1:nz), values(1:nz), matrix, &
            status)
        call check(status_ok(status), 'CSC construction')
    end subroutine grid_matrix

    subroutine add(r, c, v, rows, cols, values, nz)
        integer, intent(in) :: r, c
        real(dp), intent(in) :: v
        integer, intent(inout) :: rows(:), cols(:), nz
        real(dp), intent(inout) :: values(:)
        nz = nz+1
        rows(nz) = r
        cols(nz) = c
        values(nz) = v
    end subroutine add

    subroutine supernodal(solver, nthreads)
        type(sparse_solver_t), intent(inout) :: solver
        integer, intent(in) :: nthreads
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry = FORTSPARSE_SYM_GENERAL
        solver%options%lu_method = FORTSPARSE_LU_SUPERNODAL
        solver%options%ordering = FORTSPARSE_ORDER_AMD
        solver%options%nthreads = nthreads
    end subroutine supernodal

    subroutine fresh_solution(matrix, nthreads, x, stats)
        type(csc_t), intent(in) :: matrix
        integer, intent(in) :: nthreads
        real(dp), intent(out) :: x(:)
        type(sparse_factor_stats_t), intent(out) :: stats
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        call supernodal(solver, nthreads)
        call sparse_factor(solver, matrix, status)
        call check(status_ok(status), 'fresh factor: '//trim(status%msg))
        call solve(solver, x, stats)
        call sparse_destroy(solver)
    end subroutine fresh_solution

    subroutine solve(solver, x, stats)
        type(sparse_solver_t), intent(inout) :: solver
        real(dp), intent(out) :: x(:)
        type(sparse_factor_stats_t), intent(out) :: stats
        type(fortsparse_status_t) :: status
        real(dp) :: b(size(x))
        integer :: i
        b = [(cos(0.37_dp*real(i,dp))+0.25_dp, i=1, size(x))]
        call sparse_factor_stats(solver, stats, status)
        call check(status_ok(status), 'factor statistics')
        call sparse_solve(solver, b, x, status)
        call check(status_ok(status), 'solve')
    end subroutine solve

    subroutine same_factor(solver, matrix, nthreads, label, active)
        ! The retained solver's solution and pivot statistics equal those of
        ! a fresh factor of matrix bitwise.
        type(sparse_solver_t), intent(inout) :: solver
        type(csc_t), intent(in) :: matrix
        integer, intent(in) :: nthreads
        character(*), intent(in) :: label
        integer, intent(out) :: active
        type(sparse_factor_stats_t) :: s1, s2
        real(dp) :: x(matrix%ncol), y(matrix%ncol)
        call solve(solver, x, s1)
        call fresh_solution(matrix, nthreads, y, s2)
        call check(s1%native_factor_family == 'LU_SUPERNODAL' .and. &
            s2%native_factor_family == 'LU_SUPERNODAL', label//': supernodal family')
        call check(s1%active_threads == s2%active_threads, label//': thread team')
        call check(all(x == y), label//': bitwise solution')
        call check(s1%min_pivot == s2%min_pivot .and. &
            s1%pivot_growth == s2%pivot_growth .and. &
            s1%perturbed_pivots == s2%perturbed_pivots .and. &
            s1%block_pivot_swaps == s2%block_pivot_swaps, label//': pivot statistics')
        active = s1%active_threads
    end subroutine same_factor

    subroutine check_residual(solver, matrix, label)
        ! A sub-pattern factor lives on the retained (larger) structure, so
        ! a fresh factor is no bitwise oracle: the residual is.
        type(sparse_solver_t), intent(inout) :: solver
        type(csc_t), intent(in) :: matrix
        character(*), intent(in) :: label
        type(sparse_factor_stats_t) :: stats
        real(dp) :: x(matrix%ncol), r(matrix%ncol), b(matrix%ncol), scale
        integer :: i, j, p
        call solve(solver, x, stats)
        b = [(cos(0.37_dp*real(i,dp))+0.25_dp, i=1, size(x))]
        r = b
        scale = maxval(abs(b))
        do j = 1, matrix%ncol
            do p = matrix%col_ptr(j), matrix%col_ptr(j+1)-1
                i = matrix%row_idx(p)
                r(i) = r(i)-matrix%val(p)*x(j)
                scale = max(scale, abs(matrix%val(p)*x(j)))
            end do
        end do
        call check(maxval(abs(r)) <= 64.0_dp*epsilon(1.0_dp)*scale, label//': residual')
    end subroutine check_residual

    subroutine check_patterns(nthreads, active)
        integer, intent(in) :: nthreads
        integer, intent(out) :: active
        type(sparse_solver_t) :: solver, searched
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: s1, s2
        type(csc_t) :: matrix, padded
        real(dp), allocatable :: x(:), y(:)
        logical :: reused
        call supernodal(solver, nthreads)
        call supernodal(searched, nthreads)
        call grid_matrix(1.0_dp, 0, matrix)
        allocate(x(matrix%ncol), y(matrix%ncol))
        call sparse_factor(solver, matrix, status)
        call check(status_ok(status), 'factor: '//trim(status%msg))
        call sparse_factor(searched, matrix, status)
        call check(status_ok(status), 'second factor')
        call same_factor(solver, matrix, nthreads, 'first factor', active)
        ! Same pattern, new values: the analysis and its assembly map.
        call grid_matrix(1.75_dp, 0, matrix)
        matrix%val(1) = matrix%val(1)+0.5_dp
        call sparse_refactor(solver, matrix, status, reused)
        call check(status_ok(status) .and. reused, 'same-pattern refactor reuses')
        call check_residual(solver, matrix, 'same pattern')
        ! Equal values through the map (explicit zeros) and by search.
        call grid_matrix(0.5_dp, 3, padded)
        call grid_matrix(0.5_dp, 1, matrix)
        call check(padded%nnz > matrix%nnz, 'explicit zeros kept')
        call sparse_refactor(solver, padded, status, reused)
        call check(status_ok(status) .and. reused, 'padded refactor reuses')
        call sparse_refactor(searched, matrix, status, reused)
        call check(status_ok(status) .and. reused, 'sub-pattern refactor reuses')
        call solve(solver, x, s1)
        call solve(searched, y, s2)
        call check(all(x == y), 'map and search: bitwise solution')
        call check(s1%min_pivot == s2%min_pivot .and. &
            s1%pivot_growth == s2%pivot_growth .and. &
            s1%perturbed_pivots == s2%perturbed_pivots .and. &
            s1%block_pivot_swaps == s2%block_pivot_swaps, 'map and search: pivots')
        call check_residual(searched, matrix, 'sub-pattern')
        ! Long-range entries outside the structure: a new analysis.
        call grid_matrix(1.25_dp, 2, matrix)
        call sparse_refactor(solver, matrix, status, reused)
        call check(status_ok(status) .and. .not. reused, 'new pattern analyzes again')
        call same_factor(solver, matrix, nthreads, 'new pattern', active)
        ! The base pattern lies inside the new structure: searched assembly.
        call grid_matrix(2.0_dp, 0, matrix)
        call sparse_refactor(solver, matrix, status, reused)
        call check(status_ok(status) .and. reused, 'base pattern after a new analysis')
        call check_residual(solver, matrix, 'base pattern again')
        call sparse_destroy(solver)
        call sparse_destroy(searched)
    end subroutine check_patterns

    subroutine check_threads()
        ! nthreads changed between refactors of the same pattern.
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        type(csc_t) :: matrix
        integer :: active, outer
        logical :: reused
        outer = 1
        !$ outer = omp_get_max_threads()
        call supernodal(solver, 2)
        call grid_matrix(1.0_dp, 0, matrix)
        call sparse_factor(solver, matrix, status)
        call check(status_ok(status), 'two-thread factor')
        call same_factor(solver, matrix, 2, 'two threads', active)
        solver%options%nthreads = 1
        call grid_matrix(1.5_dp, 0, matrix)
        call sparse_refactor(solver, matrix, status, reused)
        call check(status_ok(status), 'one-thread refactor')
        call same_factor(solver, matrix, 1, 'one thread', active)
        call check(active == 1, 'nthreads = 1 runs one thread')
        !$ call check(omp_get_max_threads() == outer, 'OpenMP thread count restored')
        solver%options%nthreads = 2
        call grid_matrix(1.25_dp, 0, matrix)
        call sparse_refactor(solver, matrix, status, reused)
        call check(status_ok(status), 'two-thread refactor')
        call same_factor(solver, matrix, 2, 'two threads again', active)
        call check(active == team, 'nthreads = 2 restores the planned team')
        call sparse_destroy(solver)
    end subroutine check_threads

    subroutine check(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (ok) return
        print '(a)', 'FAIL: '//message
        error stop 1
    end subroutine check
end program test_fortsparse_slu_refactor
