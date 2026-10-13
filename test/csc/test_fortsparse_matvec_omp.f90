program test_fortsparse_matvec_omp
    ! Threaded CSC matvec against a serial reference and a dense oracle, real
    ! and complex, for several thread counts, including a skewed matrix with
    ! empty columns and one dense column. Agreement is to summation-order
    ! roundoff (few ulp).
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: dp, csc_t, csc_z_t, csc_from_triplet, csc_matvec, &
        fortsparse_status_t
    !$ use omp_lib, only: omp_set_num_threads
    implicit none

    integer :: nfail, nth
    integer, parameter :: threads(4) = [1, 2, 3, 4]

    nfail = 0
    do nth = 1, size(threads)
        !$ call omp_set_num_threads(threads(nth))
        call run_case(400, 500, 0.2_dp, .false., threads(nth), nfail)
        call run_case(3000, 3000, 0.01_dp, .true., threads(nth), nfail)
        call run_case(50000, 300, 0.001_dp, .false., threads(nth), nfail)
        call run_case(40, 30, 0.3_dp, .false., threads(nth), nfail)
    end do
    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"

contains

    subroutine run_case(nr, nc, dens, skew, nthr, nfail)
        integer, intent(in) :: nr, nc, nthr
        real(dp), intent(in) :: dens
        logical, intent(in) :: skew
        integer, intent(inout) :: nfail
        type(csc_t) :: A
        type(csc_z_t) :: Z
        type(fortsparse_status_t) :: st
        integer, allocatable :: r(:), c(:)
        real(dp), allocatable :: v(:), vi(:), x(:), xi(:), yref(:), y(:), dense(:, :)
        complex(dp), allocatable :: zx(:), zref(:), zy(:), zdense(:, :)
        real(dp), allocatable :: u(:, :)
        integer :: i, j, k, cap
        real(dp) :: tol, e1, e2, e3, e4

        allocate (u(nr, nc))
        call random_number(u)
        if (skew) then
            u(:, 1:nc:3) = 1.0_dp                  ! will be zeroed below
            where (u > 0.9999_dp) u = 0.0_dp
        end if
        cap = nr*nc
        allocate (r(cap), c(cap), v(cap), vi(cap))
        allocate (dense(nr, nc), zdense(nr, nc))
        dense = 0.0_dp; zdense = (0.0_dp, 0.0_dp)
        k = 0
        do j = 1, nc
            if (skew .and. mod(j, 3) == 1 .and. j /= 1) cycle   ! empty columns
            do i = 1, nr
                if (u(i, j) < dens .or. (skew .and. j == 1)) then
                    k = k + 1
                    r(k) = i; c(k) = j
                    call random_number(v(k)); call random_number(vi(k))
                    v(k) = v(k) - 0.5_dp; vi(k) = vi(k) - 0.5_dp
                    dense(i, j) = v(k)
                    zdense(i, j) = cmplx(v(k), vi(k), dp)
                end if
            end do
        end do
        call csc_from_triplet(nr, nc, r(1:k), c(1:k), v(1:k), A, st)
        call csc_from_triplet(nr, nc, r(1:k), c(1:k), cmplx(v(1:k), vi(1:k), dp), Z, st)
        allocate (x(nc), xi(nc), zx(nc), yref(nr), zref(nr))
        call random_number(x); call random_number(xi)
        x = x - 0.5_dp; zx = cmplx(x, xi - 0.5_dp, dp)
        ! serial reference: column-order accumulation from the CSC arrays
        yref = 0.0_dp; zref = (0.0_dp, 0.0_dp)
        do j = 1, nc
            do i = A%col_ptr(j), A%col_ptr(j + 1) - 1
                yref(A%row_idx(i)) = yref(A%row_idx(i)) + A%val(i)*x(j)
                zref(Z%row_idx(i)) = zref(Z%row_idx(i)) + Z%val(i)*zx(j)
            end do
        end do
        y = csc_matvec(A, x)
        zy = csc_matvec(Z, zx)
        tol = 64*epsilon(1.0_dp)*max(1.0_dp, real(k, dp)/nr)
        e1 = maxval(abs(y - yref))/max(1.0_dp, maxval(abs(yref)))
        e2 = maxval(abs(zy - zref))/max(1.0_dp, maxval(abs(zref)))
        e3 = maxval(abs(y - matmul(dense, x)))/max(1.0_dp, maxval(abs(yref)))
        e4 = maxval(abs(zy - matmul(zdense, zx)))/max(1.0_dp, maxval(abs(zref)))
        if (max(e1, e2, e3, e4) > tol .or. size(y) /= nr) then
            nfail = nfail + 1
            write (error_unit, "(a,3i7,a,i0,4es10.2)") "FAIL nr,nc,nnz ", nr, nc, k, &
                " threads ", nthr, e1, e2, e3, e4
        end if
    end subroutine run_case

end program test_fortsparse_matvec_omp
