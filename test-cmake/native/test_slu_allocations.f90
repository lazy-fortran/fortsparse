program test_slu_allocations
    ! Linux malloc-wrapper oracle for the supernodal LU: the analysis live
    ! and peak payload equal its ledger, and factor, refactor and all solve
    ! operations call the allocator zero times, also for an OpenMP subtree
    ! task team (FortSparse code only: the shared OpenMP runtime's own task
    ! descriptors are outside this static wrapper).
    use fortsparse_kinds, only: dp, i8
    use fortsparse_memory, only: memory_budget_t
    use fortsparse_status, only: fortsparse_status_t, FORTSPARSE_OK
    use fortsparse_native_slu, only: slu_factor_t
    use, intrinsic :: iso_c_binding, only: c_int64_t
    implicit none
    interface
        subroutine begin(nth) bind(C, name='fsparse_spd_probe_begin')
            import c_int64_t
            integer(c_int64_t), value :: nth
        end subroutine
        subroutine finish(count, current, maximum) bind(C, name='fsparse_spd_probe_end')
            import c_int64_t
            integer(c_int64_t), intent(out) :: count, current, maximum
        end subroutine
    end interface
    integer, parameter :: n = 120
    type(slu_factor_t) :: lu
    type(fortsparse_status_t) :: status
    type(memory_budget_t) :: budget
    integer, allocatable :: ptr(:), rows(:)
    real(dp), allocatable :: values(:), x(:), block(:,:)
    integer(c_int64_t) :: count, current, maximum
    integer :: blas

    call build(ptr, rows, values)
    allocate(x(n), block(n,3))
    do blas = 0, 1
        budget = memory_budget_t()
        call begin(0_c_int64_t)
        call lu%analyze(n, ptr, rows, values, status, budget=budget, force_blas=blas == 1)
        call finish(count, current, maximum)
        call check(status%code == FORTSPARSE_OK, 'analysis')
        call check(current == lu%memory_bytes(), 'live allocations equal the layout')
        call check(maximum == budget%peak, 'allocation peak equals the ledger peak')
        call begin(0_c_int64_t)
        call lu%factor(ptr, rows, values, status)
        call check(status%code == FORTSPARSE_OK, 'factor')
        x = 1.0_dp
        call lu%solve(x, status, 0)
        call lu%solve(x, status, 1)
        block = 1.0_dp
        call lu%solve_block(block, status, 1)
        values = 1.25_dp*values
        call lu%factor(ptr, rows, values, status)
        call lu%solve(x, status, 2)
        call finish(count, current, maximum)
        call check(status%code == FORTSPARSE_OK, 'refactor and solves')
        call check(count == 0, 'factor, refactor and solves call no allocator')
        call lu%free(budget)
        call check(budget%current == 0, 'ledger released')
    end do
    call threaded_team()
    print '(a)', 'supernodal LU allocation oracle passed'
contains
    subroutine threaded_team()
        integer, parameter :: m = 70
        integer, allocatable :: gptr(:), grows(:)
        real(dp), allocatable :: gvalues(:), y(:)
        integer :: ix, iy, dx, dy, i, j, k
        allocate(gptr(m*m+1), grows(9*m*m), gvalues(9*m*m), y(m*m))
        k = 0
        gptr(1) = 1
        do iy = 1, m
            do ix = 1, m
                j = ix+(iy-1)*m
                do dy = -1, 1
                    do dx = -1, 1
                        if (ix+dx < 1 .or. ix+dx > m .or. iy+dy < 1 .or. iy+dy > m) cycle
                        i = ix+dx+(iy+dy-1)*m
                        k = k+1
                        grows(k) = i
                        gvalues(k) = -1.0_dp+0.5_dp*dx-0.25_dp*dy+0.01_dp*mod(i+3*j, 7)
                        if (i == j) gvalues(k) = 10.0_dp
                    end do
                end do
                gptr(j+1) = k+1
            end do
        end do
        budget = memory_budget_t()
        call begin(0_c_int64_t)
        call lu%analyze(m*m, gptr, grows(:k), gvalues(:k), status, budget=budget, nthreads=4)
        call finish(count, current, maximum)
        call check(status%code == FORTSPARSE_OK, 'threaded analysis')
        call check(current == lu%memory_bytes(), 'threaded live allocations equal the layout')
        call check(maximum == budget%peak, 'threaded allocation peak equals the ledger peak')
        call check(lu%requested_threads == 1 .or. lu%active_threads > 1, 'parallel team planned')
        call begin(0_c_int64_t)
        call lu%factor(gptr, grows(:k), gvalues(:k), status)
        call check(status%code == FORTSPARSE_OK, 'threaded factor')
        y = 1.0_dp
        call lu%solve(y, status, 0)
        gvalues = 0.75_dp*gvalues
        call lu%factor(gptr, grows(:k), gvalues(:k), status)
        call lu%solve(y, status, 1)
        call finish(count, current, maximum)
        call check(status%code == FORTSPARSE_OK, 'threaded refactor and solves')
        call check(lu%requested_threads == 1 .or. lu%active_threads > 1, 'parallel team ran')
        call check(count == 0, 'threaded factor, refactor and solves call no allocator')
        call lu%free(budget)
        call check(budget%current == 0, 'threaded ledger released')
    end subroutine threaded_team

    subroutine build(ptr, rows, values)
        ! Near-symmetric 2D-like pattern with a zero diagonal on odd columns.
        integer, allocatable, intent(out) :: ptr(:), rows(:)
        real(dp), allocatable, intent(out) :: values(:)
        integer :: j, k, offsets(5), i, t
        offsets = [-11, -1, 0, 1, 11]
        allocate(ptr(n+1), rows(5*n), values(5*n))
        k = 0
        ptr(1) = 1
        do j = 1, n
            do t = 1, 5
                i = j+offsets(t)
                if (i < 1 .or. i > n) cycle
                if (offsets(t) == 0 .and. mod(j,2) == 1) cycle
                k = k+1
                rows(k) = i
                values(k) = 1.0_dp+real(mod(7*i+3*j, 13), dp)
                if (offsets(t) == 0) values(k) = 20.0_dp
            end do
            ptr(j+1) = k+1
        end do
        rows = rows(1:k)
        values = values(1:k)
    end subroutine build

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (condition) return
        print '(a)', 'FAIL: '//message
        error stop 1
    end subroutine check
end program test_slu_allocations
