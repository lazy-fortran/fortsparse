program test_complex_fixture_io
    use, intrinsic :: iso_fortran_env, only: int64
    use native_family_fixtures, only: dp
    use native_family_real_io, only: read_complex_csc, read_complex_array, &
                                    write_complex_array
    implicit none
    integer, allocatable :: ptr(:), row(:)
    complex(dp), allocatable :: val(:), values(:, :)
    complex(dp) :: original(3, 2), expected(3)
    integer :: unit, status, n
    character(*), parameter :: path = 'test-complex-general.mtx'

    call matrix(['3 3 3       ', '1 1 2 -1    ', '1 3 1 0.5   ', '3 3 4 2     '])
    call read_complex_csc(path, n, ptr, row, val, status)
    if (status /= 0 .or. n /= 3) error stop 'valid complex CSC rejected'
    if (any(ptr /= [1, 2, 2, 4])) error stop 'complex empty column offsets'
    if (any(row /= [1, 1, 3])) error stop 'complex CSC rows'
    expected = [cmplx(2.0_dp, -1.0_dp, dp), &
                cmplx(1.0_dp, 0.5_dp, dp), cmplx(4.0_dp, 2.0_dp, dp)]
    if (any(transfer(val, [0_int64], 6) /= &
            transfer(expected, [0_int64], 6))) error stop 'complex CSC components'
    call matrix(['3 3 1       ', '4 1 2 1     '])
    call rejected_matrix()
    call matrix(['3 3 2       ', '1 1 2 1     '])
    call rejected_matrix()
    call matrix(['3 3 1       ', '1 1 2 1     ', '2 2 4 1     '])
    call rejected_matrix()
    call matrix(['3 3 2       ', '1 1 2 1     ', '1 1 4 1     '])
    call rejected_matrix()
    call matrix(['3 3 2       ', '1 2 2 1     ', '1 1 4 1     '])
    call rejected_matrix()
    call matrix(['3 3 1       ', '1 0 2 1     '])
    call rejected_matrix()
    call matrix(['3 3 1       ', '1 1 2 NaN   '])
    call rejected_matrix()
    call matrix(['3 2 1       ', '1 1 2 1     '])
    call rejected_matrix()
    call matrix(['3 3 1       ', '1 1 2       '])
    call rejected_matrix()
    call matrix(['3 3 1       ', '1 1 2 1 x   '])
    call rejected_matrix()
    call remove_file()

    original = reshape([cmplx(0.0_dp, -0.5_dp, dp), &
        cmplx(tiny(1.0_dp), huge(1.0_dp), dp), cmplx(1.0_dp/3, -7.0_dp, dp), &
        cmplx(-2.0_dp, 0.25_dp, dp), cmplx(3.0_dp, -1.0_dp, dp), &
        cmplx(2.0_dp, 4.0_dp, dp)], [3, 2])
    call write_complex_array(path, original, status)
    if (status /= 0) error stop 'complex array write'
    call write_complex_array(path, original, status)
    if (status == 0) error stop 'existing complex export overwritten'
    call read_complex_array(path, values, status)
    if (status /= 0) error stop 'complex array read'
    if (any(transfer(values, [0_int64], 12) /= &
            transfer(original, [0_int64], 12))) error stop 'complex binary64 roundtrip'
    call remove_file()
    open (newunit=unit, file=path, status='new')
    write (unit, '(a)') '%%MatrixMarket matrix array complex general', '3 2', '1 2'
    close (unit)
    call read_complex_array(path, values, status)
    if (status == 0 .or. allocated(values)) error stop 'truncated complex array'
    call remove_file()
    print '(a)', 'PASS canonical complex CSC, malformed input and array roundtrip'
contains
    subroutine matrix(lines)
        character(*), intent(in) :: lines(:)
        integer :: i
        open (newunit=unit, file=path, status='replace')
        write (unit, '(a)') '%%MatrixMarket matrix coordinate complex general'
        do i = 1, size(lines)
            write (unit, '(a)') trim(lines(i))
        end do
        close (unit)
    end subroutine

    subroutine rejected_matrix()
        call read_complex_csc(path, n, ptr, row, val, status)
        if (status == 0) error stop 'malformed complex CSC accepted'
        if (allocated(ptr) .or. allocated(row) .or. allocated(val)) &
            error stop 'rejected complex CSC retained arrays'
    end subroutine

    subroutine remove_file()
        open (newunit=unit, file=path, status='old')
        close (unit, status='delete')
    end subroutine
end program
