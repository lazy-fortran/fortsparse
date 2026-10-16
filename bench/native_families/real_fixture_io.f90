module native_family_real_io
    ! Benchmark-only, canonical real MatrixMarket interchange. No solver or BLAS.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: iostat_end
    use native_family_fixtures, only: dp, family_fixture_t, make_family_fixture, &
                                      known_rhs
    implicit none
    private
    public :: read_real_csc, read_real_array, write_real_array
    public :: load_benchmark_fixture, load_benchmark_rhs, export_real_solution
contains
    integer function words(line) result(count)
        character(*), intent(in) :: line
        integer :: i
        logical :: previous, current
        count = 0
        previous = .false.
        do i = 1, len_trim(line)
            current = line(i:i) /= ' ' .and. line(i:i) /= achar(9)
            if (current .and. .not. previous) count = count + 1
            previous = current
        end do
    end function

    subroutine data_line(unit, line, status)
        integer, intent(in) :: unit
        character(*), intent(out) :: line
        integer, intent(out) :: status
        do
            read (unit, '(a)', iostat=status) line
            if (status /= 0) return
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) /= '%') return
        end do
    end subroutine

    subroutine read_real_csc(path, n, ptr, row, val, status)
        character(*), intent(in) :: path
        integer, intent(out) :: n, status
        integer, allocatable, intent(out) :: ptr(:), row(:)
        real(dp), allocatable, intent(out) :: val(:)
        character(1024) :: line
        integer :: unit, io, nr, nc, nz, p, i, j, last_i, last_j, col
        n = 0
        status = 2
        open (newunit=unit, file=path, status='old', action='read', iostat=io)
        if (io /= 0) return
        read (unit, '(a)', iostat=io) line
        if (io /= 0) goto 900
        if (trim(line) /= '%%MatrixMarket matrix coordinate real general') goto 900
        call data_line(unit, line, io)
        if (io /= 0) goto 900
        if (words(line) /= 3) goto 900
        read (line, *, iostat=io) nr, nc, nz
        if (io /= 0) goto 900
        if (nr < 1 .or. nr /= nc .or. nz < 0) goto 900
        if (nr == huge(nr) .or. nz == huge(nz)) goto 900
        allocate (ptr(nr + 1), row(nz), val(nz), stat=io)
        if (io /= 0) goto 900
        ptr(1) = 1
        col = 1
        last_i = 0
        last_j = 0
        do p = 1, nz
            call data_line(unit, line, io)
            if (io /= 0) goto 900
            if (words(line) /= 3) goto 900
            read (line, *, iostat=io) i, j, val(p)
            if (io /= 0) goto 900
            if (i < 1 .or. i > nr .or. j < 1 .or. j > nc) goto 900
            if (.not. ieee_is_finite(val(p))) goto 900
            if (j < last_j) goto 900
            if (j == last_j) then
                if (i <= last_i) goto 900
            end if
            do while (col < j)
                ptr(col + 1) = p
                col = col + 1
            end do
            row(p) = i
            last_i = i
            last_j = j
        end do
        do while (col <= nc)
            ptr(col + 1) = nz + 1
            col = col + 1
        end do
        ! Reject surplus numeric records as well as truncated files.
        call data_line(unit, line, io)
        if (io /= iostat_end) goto 900
        n = nr
        status = 0
900     close (unit)
        if (status /= 0) then
            if (allocated(ptr)) deallocate (ptr)
            if (allocated(row)) deallocate (row)
            if (allocated(val)) deallocate (val)
        end if
    end subroutine

    subroutine read_real_array(path, values, status)
        character(*), intent(in) :: path
        real(dp), allocatable, intent(out) :: values(:, :)
        integer, intent(out) :: status
        character(1024) :: line
        integer :: unit, io, nr, nc, i, j
        status = 2
        open (newunit=unit, file=path, status='old', action='read', iostat=io)
        if (io /= 0) return
        read (unit, '(a)', iostat=io) line
        if (io /= 0) goto 900
        if (trim(line) /= '%%MatrixMarket matrix array real general') goto 900
        call data_line(unit, line, io)
        if (io /= 0) goto 900
        if (words(line) /= 2) goto 900
        read (line, *, iostat=io) nr, nc
        if (io /= 0) goto 900
        if (nr < 1 .or. nc < 1) goto 900
        if (nr > huge(nr)/nc) goto 900
        allocate (values(nr, nc), stat=io)
        if (io /= 0) goto 900
        do j = 1, nc
            do i = 1, nr
                call data_line(unit, line, io)
                if (io /= 0) goto 900
                if (words(line) /= 1) goto 900
                read (line, *, iostat=io) values(i, j)
                if (io /= 0) goto 900
                if (.not. ieee_is_finite(values(i, j))) goto 900
            end do
        end do
        call data_line(unit, line, io)
        if (io /= iostat_end) goto 900
        status = 0
900     close (unit)
        if (status /= 0) then
            if (allocated(values)) deallocate (values)
        end if
    end subroutine

    subroutine write_real_array(path, values, status)
        character(*), intent(in) :: path
        real(dp), intent(in) :: values(:, :)
        integer, intent(out) :: status
        integer :: unit, io, i, j
        status = 2
        if (.not. all(ieee_is_finite(values))) return
        open (newunit=unit, file=path, status='new', action='write', iostat=io)
        if (io /= 0) return
        write (unit, '(a)', iostat=io) '%%MatrixMarket matrix array real general'
        if (io == 0) write (unit, *, iostat=io) size(values, 1), size(values, 2)
        do j = 1, size(values, 2)
            do i = 1, size(values, 1)
                if (io /= 0) exit
                write (unit, '(es26.17e3)', iostat=io) values(i, j)
            end do
        end do
        if (io == 0) status = 0
        close (unit, iostat=io)
        if (io /= 0) status = 2
    end subroutine

    subroutine load_benchmark_fixture(name, level, heldout, a, source_path)
        character(*), intent(in) :: name, source_path
        integer, intent(in) :: level
        logical, intent(in) :: heldout
        type(family_fixture_t), intent(out) :: a
        real(dp), allocatable :: values(:)
        integer :: status
        if (name /= 'real-mm') then
            call make_family_fixture(name, level, heldout, a, source_path)
            return
        end if
        call read_real_csc(source_path//'/matrix.mtx', a%n, a%ptr, a%row, &
                           values, status)
        if (status /= 0) error stop 'invalid canonical real MatrixMarket CSC'
        a%val = cmplx(values, 0.0_dp, dp)
    end subroutine

    subroutine load_benchmark_rhs(a, name, nrhs, heldout, operation, b, x, source_path)
        type(family_fixture_t), intent(in) :: a
        character(*), intent(in) :: name, source_path
        integer, intent(in) :: nrhs, operation
        logical, intent(in) :: heldout
        complex(dp), intent(out) :: b(:, :), x(:, :)
        real(dp), allocatable :: values(:, :)
        integer :: status
        if (name /= 'real-mm') then
            call known_rhs(a, name, nrhs, heldout, operation, b, x)
            return
        end if
        if (operation /= 0 .or. nrhs /= 4) &
            error stop 'real-mm requires normal rhs4 workload'
        call read_real_array(source_path//'/physical-rhs.mtx', values, status)
        call require_shape(3)
        b(:, :3) = cmplx(values, 0.0_dp, dp)
        call read_real_array(source_path//'/known-rhs.mtx', values, status)
        call require_shape(1)
        b(:, 4:) = cmplx(values, 0.0_dp, dp)
        call read_real_array(source_path//'/physical-reference.mtx', values, status)
        call require_shape(3)
        x(:, :3) = cmplx(values, 0.0_dp, dp)
        call read_real_array(source_path//'/known-x.mtx', values, status)
        call require_shape(1)
        x(:, 4:) = cmplx(values, 0.0_dp, dp)
    contains
        subroutine require_shape(columns)
            integer, intent(in) :: columns
            if (status /= 0) error stop 'invalid MatrixMarket RHS/reference'
            if (size(values, 1) /= a%n .or. size(values, 2) /= columns) &
                error stop 'MatrixMarket RHS/reference shape mismatch'
        end subroutine
    end subroutine

    subroutine export_real_solution(name, x)
        character(*), intent(in) :: name
        complex(dp), intent(in) :: x(:, :)
        character(4096) :: path
        integer :: status
        if (name /= 'real-mm') return
        call get_environment_variable('FORTSPARSE_BENCH_SOLUTION', path, &
                                      status=status)
        if (status /= 0 .or. len_trim(path) == 0) &
            error stop 'real-mm requires independent-audit solution export path'
        if (any(aimag(x) /= 0.0_dp)) error stop 'real-mm complex solution'
        call write_real_array(trim(path), real(x, dp), status)
        if (status /= 0) error stop 'MatrixMarket solution export failed'
    end subroutine
end module native_family_real_io
