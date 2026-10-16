module spd_fixtures
    ! Deterministic graph Laplacians plus a positive mass shift. Assembly is
    ! deliberately independent of every solver and sparse matrix constructor.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: fixture_t, make_fixture, fixture_product, backward_error, wall_time
    public :: read_symmetric_matrix_market
    public :: fixture_spd_update, fixture_known_rhs

    type :: fixture_t
        integer :: n = 0
        integer, allocatable :: ptr(:), row(:)
        real(dp), allocatable :: val(:)
    end type

contains

    subroutine fixture_spd_update(original, changed)
        ! Positive, nonuniform diagonal addition preserves SPD and the pattern.
        type(fixture_t), intent(in) :: original
        type(fixture_t), intent(out) :: changed
        integer :: j, p

        changed = original
        do j = 1, changed%n
            do p = changed%ptr(j), changed%ptr(j + 1) - 1
                if (changed%row(p) /= j) cycle
                changed%val(p) = changed%val(p) + &
                    0.01_dp*real(1 + mod(j, 7), dp)/7.0_dp
                exit
            end do
        end do
    end subroutine fixture_spd_update

    subroutine fixture_known_rhs(a, nrhs, repetition, rhs, reference)
        type(fixture_t), intent(in) :: a
        integer, intent(in) :: nrhs, repetition
        real(dp), intent(out) :: rhs(:, :), reference(:, :)
        integer :: j

        do j = 1, nrhs
            reference(:, j) = 1.0_dp + real(j, dp)/8.0_dp + &
                real(repetition - 1, dp)/16.0_dp
            reference(1:a%n:2, j) = -0.5_dp - real(repetition - 1, dp)/32.0_dp
            call fixture_product(a, reference(:, j), rhs(:, j))
        end do
    end subroutine fixture_known_rhs

    subroutine make_fixture(family, level, held_out, a, matrix_path)
        character(*), intent(in) :: family
        integer, intent(in) :: level
        logical, intent(in) :: held_out
        type(fixture_t), intent(out) :: a
        character(*), optional, intent(in) :: matrix_path
        integer, parameter :: bands(3) = [128, 512, 2048]
        integer, parameter :: grids2(7) = [16, 24, 32, 48, 64, 96, 128]
        integer, parameter :: grids3(7) = [5, 7, 9, 12, 16, 20, 24]
        integer :: nx, ny, nz, radius, col, r, k, d, axis, stride
        integer :: coordinate, extent, count, position, rows(9), temporary
        real(dp) :: weights(9), diagonal, weight, temp_weight, anisotropy(3)

        if (family == "bcsstk09") then
            if (.not. present(matrix_path)) error stop "matrix path required"
            if (held_out) error stop "archived matrix has no held-out variation"
            call read_symmetric_matrix_market(matrix_path, a, k)
            if (k /= 0) error stop "invalid symmetric MatrixMarket fixture"
            if (a%n /= 1083 .or. size(a%val) /= 18437) &
                error stop "HB/bcsstk09 dimensions differ from frozen metadata"
            return
        end if

        if (level < 1 .or. level > 7) error stop "invalid fixture level"
        nx = 1
        ny = 1
        nz = 1
        radius = 1
        anisotropy = [1.0_dp, 0.01_dp, 0.1_dp]
        select case (family)
        case ("band")
            if (level > 3) error stop "invalid band level"
            nx = bands(level)
            radius = 4
        case ("grid2")
            nx = grids2(level)
            ny = nx
        case ("grid3")
            nx = grids3(level)
            ny = nx
            nz = nx
        case default
            error stop "unknown fixture family"
        end select
        a%n = nx*ny*nz
        allocate (a%ptr(a%n + 1), a%row(9*a%n), a%val(9*a%n))
        position = 1
        do col = 1, a%n
            a%ptr(col) = position
            count = 0
            diagonal = 0.25_dp
            stride = 1
            do axis = 1, 3
                extent = nx
                if (axis == 2) extent = ny
                if (axis == 3) extent = nz
                coordinate = mod((col - 1)/stride, extent)
                do d = -radius, radius
                    if (d == 0) cycle
                    if (coordinate + d < 0) cycle
                    if (coordinate + d >= extent) cycle
                    r = col + d*stride
                    weight = anisotropy(axis)/real(abs(d), dp)
                    if (held_out) weight = weight*(1.0_dp + &
                        0.2_dp*real(mod(min(col, r), 7), dp)/7.0_dp)
                    count = count + 1
                    rows(count) = r
                    weights(count) = -weight
                    diagonal = diagonal + weight
                end do
                stride = stride*extent
            end do
            count = count + 1
            rows(count) = col
            weights(count) = diagonal
            ! At most nine entries: insertion sorting keeps the CSC canonical.
            do k = 2, count
                temporary = rows(k)
                temp_weight = weights(k)
                r = k - 1
                do while (r >= 1)
                    if (rows(r) <= temporary) exit
                    rows(r + 1) = rows(r)
                    weights(r + 1) = weights(r)
                    r = r - 1
                end do
                rows(r + 1) = temporary
                weights(r + 1) = temp_weight
            end do
            a%row(position:position + count - 1) = rows(:count)
            a%val(position:position + count - 1) = weights(:count)
            position = position + count
        end do
        a%ptr(a%n + 1) = position
        a%row = a%row(:position - 1)
        a%val = a%val(:position - 1)
    end subroutine

    subroutine read_symmetric_matrix_market(path, a, status)
        ! Bounded benchmark reader: real symmetric coordinate, lower triangle.
        ! Mirror off-diagonal entries, sort columns and sum duplicate records.
        character(*), intent(in) :: path
        type(fixture_t), intent(out) :: a
        integer, intent(out) :: status
        integer :: unit, io, n, m, stored, count, i, j, k, p, start, finish
        integer :: position, previous, output, temporary
        integer, allocatable :: rows(:), columns(:), counts(:), next(:)
        real(dp), allocatable :: values(:)
        real(dp) :: value, temp_value
        character(512) :: line

        status = 2
        open (newunit=unit, file=path, status="old", action="read", iostat=io)
        if (io /= 0) return
        read (unit, '(a)', iostat=io) line
        if (io /= 0) goto 900
        if (trim(line) /= "%%MatrixMarket matrix coordinate real symmetric") goto 900
        call next_matrix_record(unit, line, io)
        if (io /= 0) goto 900
        read (line, *, iostat=io) n, m, stored
        if (io /= 0) goto 900
        if (n < 1 .or. n == huge(n) .or. n /= m .or. stored < 0) goto 900
        if (stored > (huge(stored) - 1)/2) goto 900
        allocate (rows(2*stored), columns(2*stored), values(2*stored), stat=io)
        if (io /= 0) then
            status = 6
            goto 900
        end if
        count = 0
        do k = 1, stored
            call next_matrix_record(unit, line, io)
            if (io /= 0) goto 900
            read (line, *, iostat=io) i, j, value
            if (io /= 0) goto 900
            if (i < 1 .or. i > n .or. j < 1 .or. j > n) goto 900
            if (i < j .or. .not. ieee_is_finite(value)) goto 900
            count = count + 1
            rows(count) = i
            columns(count) = j
            values(count) = value
            if (i /= j) then
                count = count + 1
                rows(count) = j
                columns(count) = i
                values(count) = value
            end if
        end do
        call next_matrix_record(unit, line, io)
        if (io >= 0) goto 900 ! Reject unannounced extra coordinate records.
        allocate (a%ptr(n + 1), a%row(count), a%val(count), counts(n), &
            next(n), stat=io)
        if (io /= 0) then
            status = 6
            goto 900
        end if
        counts = 0
        do k = 1, count
            counts(columns(k)) = counts(columns(k)) + 1
        end do
        a%ptr(1) = 1
        do j = 1, n
            a%ptr(j + 1) = a%ptr(j) + counts(j)
        end do
        next = a%ptr(:n)
        do k = 1, count
            position = next(columns(k))
            a%row(position) = rows(k)
            a%val(position) = values(k)
            next(columns(k)) = position + 1
        end do
        output = 1
        do j = 1, n
            start = a%ptr(j)
            finish = a%ptr(j + 1) - 1
            do k = start + 1, finish
                temporary = a%row(k)
                temp_value = a%val(k)
                p = k - 1
                do while (p >= start)
                    if (a%row(p) <= temporary) exit
                    a%row(p + 1) = a%row(p)
                    a%val(p + 1) = a%val(p)
                    p = p - 1
                end do
                a%row(p + 1) = temporary
                a%val(p + 1) = temp_value
            end do
            a%ptr(j) = output
            previous = 0
            do k = start, finish
                if (a%row(k) == previous) then
                    a%val(output - 1) = a%val(output - 1) + a%val(k)
                    if (.not. ieee_is_finite(a%val(output - 1))) goto 900
                else
                    previous = a%row(k)
                    a%row(output) = previous
                    a%val(output) = a%val(k)
                    output = output + 1
                end if
            end do
        end do
        a%ptr(n + 1) = output
        a%row = a%row(:output - 1)
        a%val = a%val(:output - 1)
        a%n = n
        status = 0
900     close (unit)
    end subroutine

    subroutine next_matrix_record(unit, line, io)
        integer, intent(in) :: unit
        character(*), intent(out) :: line
        integer, intent(out) :: io
        do
            read (unit, '(a)', iostat=io) line
            if (io /= 0) return
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == "%") cycle
            return
        end do
    end subroutine

    subroutine fixture_product(a, x, y)
        type(fixture_t), intent(in) :: a
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer :: col, p
        y = 0.0_dp
        do col = 1, a%n
            do p = a%ptr(col), a%ptr(col + 1) - 1
                y(a%row(p)) = y(a%row(p)) + a%val(p)*x(col)
            end do
        end do
    end subroutine

    function backward_error(a, x, b) result(error)
        type(fixture_t), intent(in) :: a
        real(dp), intent(in) :: x(:), b(:)
        real(dp) :: error, residual(a%n), row_sum(a%n), denominator
        integer :: col, p
        call fixture_product(a, x, residual)
        row_sum = 0.0_dp
        do col = 1, a%n
            do p = a%ptr(col), a%ptr(col + 1) - 1
                row_sum(a%row(p)) = row_sum(a%row(p)) + abs(a%val(p))
            end do
        end do
        denominator = maxval(row_sum)*maxval(abs(x)) + maxval(abs(b))
        error = maxval(abs(residual - b))/max(denominator, tiny(1.0_dp))
    end function

    function wall_time() result(seconds)
        real(dp) :: seconds
        integer(int64) :: count, rate
        call system_clock(count, rate)
        seconds = real(count, dp)/real(rate, dp)
    end function
end module
