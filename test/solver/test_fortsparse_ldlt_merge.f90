program test_fortsparse_ldlt_merge
    use fortsparse, only: dp, csc_t, csc_from_triplet, fortsparse_status_t, &
                          FORTSPARSE_OK
    use fortsparse_native_ldlt, only: native_ldlt_t
    implicit none
    integer, parameter :: n = 160
    type(native_ldlt_t) :: factor
    type(csc_t) :: matrix
    type(fortsparse_status_t) :: status
    real(dp) :: a(n, n), values(n*n), x(n), b(n), got(n), reference(n)
    real(dp) :: v(n)
    integer :: rows(n*n), cols(n*n), i, j, entry, mode

    ! A dense rank-one perturbation of a separated signed diagonal has no
    ! structural fill, so LOW and default have identical row capacities and
    ! differ only in row merging, which must not change the result.
    v = [(real(mod(i, 7) + 1, dp)/128, i=1, n)]
    do j = 1, n
        do i = 1, n
            a(i, j) = v(i)*v(j)
        end do
        a(j, j) = a(j, j) + 3*real((-1)**j, dp)
    end do
    entry = 0
    do j = 1, n
        do i = 1, n
            entry = entry + 1
            rows(entry) = i
            cols(entry) = j
            values(entry) = a(i, j)
        end do
    end do
    call csc_from_triplet(n, n, rows, cols, values, matrix, status)
    call checked()
    x = [(sin(real(i, dp)), i=1, n)]
    b = matmul(a, x)
    do mode = 1, 2
        call factor%factor(matrix, status, low_memory=mode == 2, &
                           nthreads=1)
        call checked()
        call factor%solve(b, got, status)
        call checked()
        if (maxval(abs(got - x)) > 2e-13_dp) error stop 'known solution'
        if (maxval(abs(matmul(a, got) - b)) > 2e-13_dp) &
            error stop 'original matrix residual'
        if (mode == 1) reference = got
        if (maxval(abs(got - reference)) > 2e-13_dp) &
            error stop 'merge changed scalar result'
        call factor%free()
    end do

    ! Also exercise independent row ownership with more than one worker.
    call factor%factor(matrix, status, nthreads=2)
    call checked()
!$  if (factor%numeric_threads() /= 2) error stop 'merge team not observed'
    call factor%solve(b, got, status)
    call checked()
    if (maxval(abs(got - reference)) > 2e-13_dp) &
        error stop 'threaded merge changed result'
contains
    subroutine checked()
        if (status%code == FORTSPARSE_OK) return
        print *, status%code, trim(status%msg)
        error stop 'LDLT merge operation'
    end subroutine checked
end program test_fortsparse_ldlt_merge
