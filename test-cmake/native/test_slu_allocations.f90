program test_slu_allocations
    ! Linux malloc-wrapper oracle for the supernodal LU: the analysis live
    ! and peak payload equal its ledger, and factor, refactor and all solve
    ! operations call the allocator zero times.
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
    print '(a)', 'supernodal LU allocation oracle passed'
contains
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
