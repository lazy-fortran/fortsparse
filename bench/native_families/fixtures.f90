module native_family_fixtures
    ! Benchmark-only numeric data. No solver, FortSparse or BLAS dependency.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: dp, family_fixture_t, make_family_fixture, fixture_action
    public :: known_rhs, fixture_backward_error, wall_time
    type :: family_fixture_t
        integer :: n = 0, mtype = 11
        logical :: is_complex = .false.
        integer, allocatable :: ptr(:), row(:)
        complex(dp), allocatable :: val(:)
    end type
contains
    subroutine make_family_fixture(name, level, heldout, a, mini5_path)
        character(*), intent(in) :: name
        integer, intent(in) :: level
        logical, intent(in) :: heldout
        type(family_fixture_t), intent(out) :: a
        character(*), intent(in) :: mini5_path
        integer :: n, count, j, i, p, unit, nr, nc, nz, io
        integer, allocatable :: rows(:), cols(:), cursor(:)
        complex(dp), allocatable :: values(:)
        real(dp) :: delta, imag, value
        complex(dp) :: z
        character(256) :: line

        if (level < 1) error stop 'positive fixture level required'
        delta = 0.21_dp
        imag = 0.13_dp
        if (heldout) then
            delta = 0.37_dp
            imag = 0.19_dp
        end if
        select case (trim(name))
        case ('neo2-mini5')
            open (newunit=unit, file=mini5_path, status='old', action='read', iostat=io)
            if (io /= 0) error stop 'NEO-2 mini5 source file unavailable'
            do
                read (unit, '(a)', iostat=io) line
                if (io /= 0) error stop 'mini5 header'
                if (line(1:1) /= '%') exit
            end do
            read (line, *) nr, nc, nz
            if (nr /= 5 .or. nc /= 5 .or. nz /= 22) error stop 'mini5 identity'
            n = nr
            allocate (rows(nz), cols(nz), values(nz))
            do p = 1, nz
                read (unit, *, iostat=io) rows(p), cols(p), value
                if (io /= 0) error stop 'mini5 numeric data'
                values(p) = cmplx(value, 0.0_dp, dp)
            end do
            close (unit)
            count = nz
        case ('drift-real', 'drift-complex')
            if (level > 128) error stop 'drift level outside frozen domain'
            n = level*level
            allocate (rows(5*n), cols(5*n), values(5*n))
            count = 0
            a%is_complex = trim(name) == 'drift-complex'
            if (a%is_complex) a%mtype = 13
            do j = 1, n
                call append(j, j, cmplx(5.5_dp + real(mod(j, 7), dp)/32, 0.0_dp, dp))
               if (a%is_complex) values(count) = values(count) + cmplx(0.0_dp, imag, dp)
                if (mod(j - 1, level) > 0) then
                    z = cmplx(-1.0_dp + delta, 0.0_dp, dp)
                    if (a%is_complex) z = z + cmplx(0.0_dp, imag, dp)
                    call append(j - 1, j, z)
                end if
                if (mod(j, level) /= 0) then
                    z = cmplx(-1.0_dp - delta, 0.0_dp, dp)
                    if (a%is_complex) z = z + cmplx(0.0_dp, -0.7_dp*imag, dp)
                    call append(j + 1, j, z)
                end if
                if (j > level) call append(j - level, j, cmplx(-0.6_dp, 0.0_dp, dp))
               if (j + level <= n) call append(j + level, j, cmplx(-1.2_dp, 0.0_dp, dp))
            end do
        case ('hpd-chain', 'ldlt-chain')
            if (level > 8192) error stop 'chain level outside frozen domain'
            n = 2*level
            allocate (rows(3*n - 2), cols(3*n - 2), values(3*n - 2))
            count = 0
            a%is_complex = trim(name) == 'hpd-chain'
            a%mtype = -2
            if (a%is_complex) a%mtype = 4
            do j = 1, n
                z = cmplx(0.0_dp, 0.0_dp, dp)
                if (a%is_complex) z = cmplx(3.0_dp + delta, 0.0_dp, dp)
                call append(j, j, z)
                if (j == n) cycle
                if (a%is_complex) then
                    z = cmplx(-0.4_dp, imag, dp)
                else if (mod(j, 2) == 1) then
                    z = cmplx(2.0_dp + delta, 0.0_dp, dp)
                else
                    z = cmplx(0.125_dp + delta/8, 0.0_dp, dp)
                end if
                call append(j, j + 1, z)
                call append(j + 1, j, conjg(z))
            end do
        case default
            error stop 'unknown frozen family fixture'
        end select
        a%n = n
        allocate (a%ptr(n + 1), a%row(count), a%val(count), cursor(n))
        a%ptr = 0
        do p = 1, count
            a%ptr(cols(p) + 1) = a%ptr(cols(p) + 1) + 1
        end do
        a%ptr(1) = 1
        do j = 1, n
            a%ptr(j + 1) = a%ptr(j + 1) + a%ptr(j)
        end do
        cursor = a%ptr(:n)
        do p = 1, count
            i = cursor(cols(p))
            a%row(i) = rows(p)
            a%val(i) = values(p)
            cursor(cols(p)) = i + 1
        end do
        do j = 1, n
            do p = a%ptr(j) + 1, a%ptr(j + 1) - 1
                i = p
                do while (i > a%ptr(j))
                    if (a%row(i - 1) <= a%row(i)) exit
                    nr = a%row(i - 1)
                    a%row(i - 1) = a%row(i)
                    a%row(i) = nr
                    z = a%val(i - 1)
                    a%val(i - 1) = a%val(i)
                    a%val(i) = z
                    i = i - 1
                end do
            end do
        end do
    contains
        subroutine append(row, col, z)
            integer, intent(in) :: row, col
            complex(dp), intent(in) :: z
            count = count + 1
            rows(count) = row
            cols(count) = col
            values(count) = z
        end subroutine
    end subroutine make_family_fixture

    subroutine fixture_action(a, x, b, operation)
        type(family_fixture_t), intent(in) :: a
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: b(:)
        integer, intent(in) :: operation
        integer :: j, p, i
        complex(dp) :: z
        b = cmplx(0.0_dp, 0.0_dp, dp)
        do j = 1, a%n
            do p = a%ptr(j), a%ptr(j + 1) - 1
                i = a%row(p)
                z = a%val(p)
                if (operation == 2) z = conjg(z)
                if (operation == 0) then
                    b(i) = b(i) + z*x(j)
                else
                    b(j) = b(j) + z*x(i)
                end if
            end do
        end do
    end subroutine

    subroutine known_rhs(a, name, nrhs, heldout, operation, b, x)
        type(family_fixture_t), intent(in) :: a
        character(*), intent(in) :: name
        integer, intent(in) :: nrhs, operation
        logical, intent(in) :: heldout
        complex(dp), intent(out) :: b(:, :), x(:, :)
        integer :: i, rhs
        real(dp) :: im
        do rhs = 1, nrhs
            do i = 1, a%n
                im = 0.0_dp
                if (a%is_complex) im = real(mod(3*i + rhs, 7) - 3, dp)/8
                x(i, rhs) = cmplx(real(mod(2*i + rhs, 11) - 5, dp)/8, im, dp)
            end do
            if (trim(name) == 'neo2-mini5' .and. .not. heldout) then
                x(:, rhs) = cmplx(0.0_dp, 0.0_dp, dp)
                x(1, rhs) = cmplx(-45.0_dp*rhs/22, 0.0_dp, dp)
                x(5, rhs) = cmplx(real(rhs, dp)/308, 0.0_dp, dp)
            end if
            call fixture_action(a, x(:, rhs), b(:, rhs), operation)
        end do
    end subroutine

    real(dp) function fixture_backward_error(a, x, b, operation) result(error)
        type(family_fixture_t), intent(in) :: a
        complex(dp), intent(in) :: x(:), b(:)
        integer, intent(in) :: operation
        complex(dp) :: residual(a%n)
        real(dp) :: rowsum(a%n), denom
        integer :: j, p, i
        call fixture_action(a, x, residual, operation)
        rowsum = 0.0_dp
        do j = 1, a%n
            do p = a%ptr(j), a%ptr(j + 1) - 1
                i = a%row(p)
                if (operation /= 0) i = j
                rowsum(i) = rowsum(i) + abs(a%val(p))
            end do
        end do
        denom = maxval(rowsum)*maxval(abs(x)) + maxval(abs(b))
        error = maxval(abs(residual - b))/max(denom, tiny(1.0_dp))
    end function

    real(dp) function wall_time() result(seconds)
        integer, parameter :: i8 = selected_int_kind(18)
        integer(i8) :: count, rate
        call system_clock(count, rate)
        seconds = real(count, dp)/real(rate, dp)
    end function
end module native_family_fixtures
