program test_fortsparse_order_wide
    ! AMD with 64-bit counts (taken by the LU when METIS or 32-bit AMD would
    ! overflow) must order a small unsymmetric pattern exactly as amd_order.
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    interface
        integer(c_int) function order_amd_general(n, nnz, ptr, row, perm, inverse) &
                bind(C, name='fsparse_native_order_amd_general')
            import :: c_int
            integer(c_int), value :: n, nnz
            integer(c_int), intent(in) :: ptr(*), row(*)
            integer(c_int), intent(out) :: perm(*), inverse(*)
        end function
        integer(c_int) function order_amd_wide(n, nnz, ptr, row, perm, inverse) &
                bind(C, name='fsparse_native_order_amd_general_wide')
            import :: c_int
            integer(c_int), value :: n, nnz
            integer(c_int), intent(in) :: ptr(*), row(*)
            integer(c_int), intent(out) :: perm(*), inverse(*)
        end function
    end interface
    integer, parameter :: n = 40
    integer(c_int) :: ptr(n + 1), row(6*n), p32(n), i32(n), p64(n), i64(n), code
    integer :: j, k, nnz
    ! Unsymmetric 1D-like band with a few long-range entries.
    nnz = 0
    ptr(1) = 1
    do j = 1, n
        do k = 1, n
            if (k == j .or. k == j + 1 .or. k == j - 2 .or. k == mod(7*j, n) + 1) then
                if (count(row(ptr(j):nnz) == k) == 0) then
                    nnz = nnz + 1
                    row(nnz) = k
                end if
            end if
        end do
        ptr(j + 1) = nnz + 1
    end do
    code = order_amd_wide(n, nnz, ptr, row, p64, i64)
    if (code == 3) then
        print '(a)', 'SKIP: external AMD not built in'
        stop
    end if
    if (code /= 0) error stop 'wide AMD failed'
    code = order_amd_general(n, nnz, ptr, row, p32, i32)
    if (code /= 0) error stop 'AMD failed'
    if (any(p64 /= p32) .or. any(i64 /= i32)) error stop 'wide AMD order differs'
    do j = 1, n
        if (i64(p64(j)) /= j) error stop 'inverse of the wide order'
    end do
    print '(a)', 'PASS: 64-bit AMD orders like amd_order'
end program
