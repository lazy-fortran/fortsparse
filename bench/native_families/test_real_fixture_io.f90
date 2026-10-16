program test_real_fixture_io
    use native_family_fixtures, only: dp
    use native_family_real_io, only: read_real_csc, read_real_array, write_real_array
    implicit none
    integer, allocatable :: ptr(:), row(:)
    real(dp), allocatable :: val(:), values(:, :)
    real(dp) :: original(3, 2)
    integer :: unit, status, n
    character(*), parameter :: path = 'test-real-general.mtx'

    call matrix(['3 3 3   ', '1 1 2   ', '1 3 1   ', '3 3 4   '])
    call read_real_csc(path, n, ptr, row, val, status)
    if (status /= 0 .or. n /= 3) error stop 'valid general CSC rejected'
    if (any(ptr /= [1, 2, 2, 4])) error stop 'empty column offsets'
    if (any(row /= [1, 1, 3])) error stop 'general CSC rows'
    if (any(val /= [2.0_dp, 1.0_dp, 4.0_dp])) error stop 'general CSC values'
    call matrix(['3 3 1   ', '4 1 2   '])
    call rejected_matrix()
    call matrix(['3 3 2   ', '1 1 2   '])
    call rejected_matrix()
    call matrix(['3 3 1   ', '1 1 2   ', '2 2 4   '])
    call rejected_matrix()
    call matrix(['3 3 2   ', '1 1 2   ', '1 1 4   '])
    call rejected_matrix()
    call matrix(['3 3 2   ', '1 2 2   ', '1 1 4   '])
    call rejected_matrix()
    call matrix(['3 3 1   ', '1 0 2   '])
    call rejected_matrix()
    call matrix(['3 3 1   ', '1 1 NaN '])
    call rejected_matrix()
    call matrix(['3 2 1   ', '1 1 2   '])
    call rejected_matrix()
    call matrix(['3 3 1   ', '1 1 2 x '])
    call rejected_matrix()
    call remove_file()

    original = reshape([0.0_dp, -0.5_dp, tiny(1.0_dp), &
                        huge(1.0_dp), 1.0_dp/3, -7.0_dp], [3, 2])
    call write_real_array(path, original, status)
    if (status /= 0) error stop 'array write'
    call read_real_array(path, values, status)
    if (status /= 0) error stop 'array read'
    if (any(values /= original)) error stop 'binary64 decimal roundtrip'
    call remove_file()
    open (newunit=unit, file=path, status='new')
    write (unit, '(a)') '%%MatrixMarket matrix array real general', '3 2', '1'
    close (unit)
    call read_real_array(path, values, status)
    if (status == 0) error stop 'truncated array accepted'
    call remove_file()
    print '(a)', 'PASS canonical real CSC, malformed input, array roundtrip'
contains
    subroutine matrix(lines)
        character(*), intent(in) :: lines(:)
        integer :: i
        open (newunit=unit, file=path, status='replace')
        write (unit, '(a)') '%%MatrixMarket matrix coordinate real general'
        do i = 1, size(lines)
            write (unit, '(a)') trim(lines(i))
        end do
        close (unit)
    end subroutine

    subroutine rejected_matrix()
        call read_real_csc(path, n, ptr, row, val, status)
        if (status == 0) error stop 'malformed CSC accepted'
        if (allocated(ptr) .or. allocated(row) .or. allocated(val)) &
            error stop 'rejected CSC retained arrays'
    end subroutine

    subroutine remove_file()
        open (newunit=unit, file=path, status='old')
        close (unit, status='delete')
    end subroutine
end program
