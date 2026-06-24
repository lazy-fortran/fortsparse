program test_fortsparse_poisson
    ! 1D Poisson with Dirichlet boundaries on the default (SuperLU) backend.
    ! Discretise -u'' = f on (0,1) with the tridiagonal stencil (-1, 2, -1)
    ! scaled by 1/h^2, manufacture u(x) = sin(pi x), and check the relative L2
    ! error and the residual.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: dp, csc_t, csc_from_triplet, sparse_solver_t, &
        sparse_factor, sparse_solve, sparse_free, fortsparse_status_t, status_ok
    implicit none

    integer,  parameter :: n = 64
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: nfail
    nfail = 0

    call run(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    subroutine run(nfail)
        integer, intent(inout) :: nfail

        type(csc_t)               :: A
        type(sparse_solver_t)     :: solver
        type(fortsparse_status_t) :: status
        integer,  allocatable     :: rows(:), cols(:)
        real(dp), allocatable     :: vals(:)
        real(dp)                  :: b(n), x(n), ue(n)
        real(dp)                  :: h, xj, relerr, res
        integer                   :: j, m

        h = 1.0_dp/real(n + 1, dp)
        call build_tridiag(n, h, rows, cols, vals)
        call csc_from_triplet(n, n, rows, cols, vals, A, status)
        call check_true("build_ok", status_ok(status), nfail)

        do j = 1, n
            xj = real(j, dp)*h
            ue(j) = sin(pi*xj)
        end do
        ! Manufacture the RHS from the discrete operator so the discrete
        ! solution equals ue to solver precision, isolating the solver from the
        ! O(h^2) truncation error of the continuous problem.
        do j = 1, n
            b(j) = matrow(n, h, ue, j)
        end do

        call sparse_factor(solver, A, status)
        call check_true("factor_ok", status_ok(status), nfail)
        call sparse_solve(solver, b, x, status)
        call check_true("solve_ok", status_ok(status), nfail)
        call sparse_free(solver)

        relerr = sqrt(sum((x - ue)**2)/sum(ue**2))
        if (relerr > 1.0e-6_dp) then
            nfail = nfail + 1
            write (error_unit, "(a,es13.6)") "FAIL [poisson_relL2] ", relerr
        end if

        res = 0.0_dp
        do m = 1, n
            res = max(res, abs(matrow(n, h, x, m) - b(m)))
        end do
        if (res > 1.0e-9_dp) then
            nfail = nfail + 1
            write (error_unit, "(a,es13.6)") "FAIL [poisson_res] ", res
        end if
    end subroutine run

    ! Triplets for the (-1, 2, -1)/h^2 Dirichlet stencil.
    subroutine build_tridiag(n, h, rows, cols, vals)
        integer,               intent(in)  :: n
        real(dp),              intent(in)  :: h
        integer,  allocatable, intent(out) :: rows(:), cols(:)
        real(dp), allocatable, intent(out) :: vals(:)
        real(dp)                           :: inv
        integer                            :: j, k

        allocate (rows(3*n - 2), cols(3*n - 2), vals(3*n - 2))
        inv = 1.0_dp/(h*h)
        k = 0
        do j = 1, n
            k = k + 1
            rows(k) = j; cols(k) = j; vals(k) = 2.0_dp*inv
            if (j > 1) then
                k = k + 1
                rows(k) = j; cols(k) = j - 1; vals(k) = -inv
            end if
            if (j < n) then
                k = k + 1
                rows(k) = j; cols(k) = j + 1; vals(k) = -inv
            end if
        end do
    end subroutine build_tridiag

    ! Apply row m of the stencil to x.
    function matrow(n, h, x, m) result(y)
        integer,  intent(in) :: n, m
        real(dp), intent(in) :: h, x(:)
        real(dp)             :: y, inv

        inv = 1.0_dp/(h*h)
        y = 2.0_dp*inv*x(m)
        if (m > 1) y = y - inv*x(m - 1)
        if (m < n) y = y - inv*x(m + 1)
    end function matrow

    subroutine check_true(label, cond, nfail)
        character(*), intent(in)    :: label
        logical,      intent(in)    :: cond
        integer,      intent(inout) :: nfail
        if (.not. cond) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a)") "FAIL [", label, "] expected .true."
        end if
    end subroutine check_true

end program test_fortsparse_poisson
