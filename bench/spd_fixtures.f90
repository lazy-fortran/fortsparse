module spd_fixtures
    ! Deterministic graph Laplacians plus a positive mass shift. Assembly is
    ! deliberately independent of every solver and sparse matrix constructor.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    implicit none
    private
    public :: fixture_t, make_fixture, fixture_product, backward_error, wall_time

    type :: fixture_t
        integer :: n = 0
        integer, allocatable :: ptr(:), row(:)
        real(dp), allocatable :: val(:)
    end type

contains

    subroutine make_fixture(family, level, held_out, a)
        character(*), intent(in) :: family
        integer, intent(in) :: level
        logical, intent(in) :: held_out
        type(fixture_t), intent(out) :: a
        integer, parameter :: bands(3) = [128, 512, 2048]
        integer, parameter :: grids2(7) = [16, 24, 32, 48, 64, 96, 128]
        integer, parameter :: grids3(7) = [5, 7, 9, 12, 16, 20, 24]
        integer :: nx, ny, nz, radius, col, r, k, d, axis, stride
        integer :: coordinate, extent, count, position, rows(9), temporary
        real(dp) :: weights(9), diagonal, weight, temp_weight, anisotropy(3)

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
