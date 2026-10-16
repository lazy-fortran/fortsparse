program test_fortsparse_external_copy
    ! External backend context copies keep backend selection and options but
    ! no factors: the destination is unfactored with null handles, the source
    ! keeps solving, and both contexts finalize without double release.
    use fortsparse
    implicit none
    integer :: nfail
    character(4096) :: helper
    integer :: length
    nfail = 0
    call run_backend(FORTSPARSE_BACKEND_SUPERLU, 'SuperLU')
    call get_environment_variable('FORTSPARSE_UMFPACK_HELPER', helper, length)
    if (length > 0) then
        call run_backend(FORTSPARSE_BACKEND_UMFPACK_IPC, 'UMFPACK IPC')
    else
        write (*, '(a)') 'SKIP: FORTSPARSE_UMFPACK_HELPER not set'
    end if
    if (nfail > 0) then
        write (*, '(i0,a)') nfail, ' check(s) failed'
        error stop 1
    end if
    write (*, '(a)') 'PASS'
contains

    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            nfail = nfail + 1
            write (*, '(2a)') 'FAIL: ', label
        end if
    end subroutine check

    subroutine run_backend(id, name)
        integer, intent(in) :: id
        character(*), intent(in) :: name
        integer, parameter :: n = 6
        type(sparse_solver_t) :: s, c
        type(fortsparse_status_t) :: status
        type(csc_t) :: a
        real(dp) :: d(n,n), xe(n), b(n), x(n)
        integer :: rows(3*n-2), cols(3*n-2), i, k, it
        real(dp) :: vals(3*n-2)
        d = 0.0_dp
        k = 0
        do i = 1, n
            d(i,i) = 4.0_dp + i
            if (i > 1) d(i,i-1) = -1.0_dp
            if (i < n) d(i,i+1) = 0.5_dp
        end do
        do i = 1, n
            do it = max(1, i - 1), min(n, i + 1)
                k = k + 1
                rows(k) = it
                cols(k) = i
                vals(k) = d(it,i)
            end do
        end do
        call csc_from_triplet(n, n, rows, cols, vals, a, status)
        if (.not. status_ok(status)) error stop 'external copy matrix'
        xe = [(real(i, dp), i = 1, n)]
        b = matmul(d, xe)
        s%backend_id = id
        call sparse_factor(s, a, status)
        if (status%code == FORTSPARSE_BACKEND_UNAVAILABLE) then
            write (*, '(3a)') 'SKIP: ', name, ' backend is not configured'
            return
        end if
        call check(status_ok(status), name//': source factor')
        do it = 1, 3
            call sparse_copy(c, s, status)
            call check(status%code == FORTSPARSE_NOT_FACTORED, &
                name//': copy reports factors not copied')
            call check(.not. c%factored .and. c%backend_id == id, &
                name//': copy unfactored with selection')
            call sparse_solve(c, b, x, status)
            call check(status%code == FORTSPARSE_NOT_FACTORED, name//': copy solve')
            call sparse_solve(s, b, x, status)
            call check(status_ok(status) .and. maxval(abs(x - xe)) < 1e-10_dp, &
                name//': source still solves')
            c = s
            call check(.not. c%factored, name//': assigned copy unfactored')
            call sparse_factor(c, a, status)
            call check(status_ok(status), name//': copy refactor')
            call sparse_solve(c, b, x, status)
            call check(status_ok(status) .and. maxval(abs(x - xe)) < 1e-10_dp, &
                name//': refactored copy x_true')
            call sparse_destroy(c)
        end do
        c = s
        call sparse_destroy(s)
        call sparse_factor(c, a, status)
        call sparse_solve(c, b, x, status)
        call check(status_ok(status) .and. maxval(abs(x - xe)) < 1e-10_dp, &
            name//': copy outlives source')
        call sparse_free(s)
        call sparse_copy(s, c, status)
        call check(status%code == FORTSPARSE_NOT_FACTORED, name//': copy back')
    end subroutine run_backend
end program test_fortsparse_external_copy
