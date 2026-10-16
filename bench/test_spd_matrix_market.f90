program test_spd_matrix_market
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spd_fixtures, only: fixture_t, read_symmetric_matrix_market, &
        fixture_product, make_fixture
    implicit none
    type(fixture_t) :: a
    integer :: unit, status
    real(dp) :: result(3)
    character(*), parameter :: path = "test-mm-fixture.mtx"

    open (newunit=unit, file=path, status="replace")
    write (unit, '(a)') "%%MatrixMarket matrix coordinate real symmetric"
    write (unit, '(a)') "% independent, unsorted lower triangle with duplicate"
    write (unit, '(a)') "3 3 6", "3 2 -0.4", "2 1 1", "2 2 3", &
        "1 1 4", "3 3 2", "3 2 -0.6"
    close (unit)
    call read_symmetric_matrix_market(path, a, status)
    if (status /= 0) error stop "valid real symmetric input rejected"
    if (a%n /= 3 .or. size(a%val) /= 7) error stop "mirrored size mismatch"
    if (any(a%ptr /= [1, 3, 6, 8])) error stop "wrong canonical column offsets"
    if (any(a%row /= [1, 2, 1, 2, 3, 2, 3])) error stop "wrong sorted rows"
    call fixture_product(a, [1.0_dp, 2.0_dp, 3.0_dp], result)
    if (maxval(abs(result - [6.0_dp, 4.0_dp, 4.0_dp])) > 1.0e-14_dp) &
        error stop "independent dense product mismatch"
    call make_fixture('matrix-market', 1, .false., a, path)
    call fixture_product(a, [1.0_dp, 2.0_dp, 3.0_dp], result)
    if (maxval(abs(result - [6.0_dp, 4.0_dp, 4.0_dp])) > 1.0e-14_dp) &
        error stop "generic MatrixMarket fixture product mismatch"

    open (newunit=unit, file=path, status="replace")
    write (unit, '(a)') "%%MatrixMarket matrix coordinate complex symmetric"
    close (unit)
    call read_symmetric_matrix_market(path, a, status)
    if (status /= 2) error stop "unsupported complex header accepted"
    open (newunit=unit, file=path, status="replace")
    write (unit, '(a)') "%%MatrixMarket matrix coordinate real symmetric"
    write (unit, '(a)') "3 3 1", "1 2 1"
    close (unit)
    call read_symmetric_matrix_market(path, a, status)
    if (status /= 2) error stop "unsupported upper triangle accepted"
    open (newunit=unit, file=path, status="replace")
    write (unit, '(a)') "%%MatrixMarket matrix coordinate real symmetric"
    write (unit, '(a)') "3 3 1", "4 2 1"
    close (unit)
    call read_symmetric_matrix_market(path, a, status)
    if (status /= 2) error stop "out of bounds row accepted"
    open (newunit=unit, file=path, status="old")
    close (unit, status="delete")
    print '(a)', "PASS MatrixMarket mirroring, sorting, duplicates and rejection"
end program
