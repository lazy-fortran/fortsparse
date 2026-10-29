module fortsparse_native_amd
    ! Shared thin bindings to the audited, allocation-free caller-workspace ABI.
    use, intrinsic :: iso_c_binding, only: c_int, c_int64_t
    implicit none
    private
    public :: bounded_amd_available, bounded_amd_size, bounded_amd_order
    public :: default_prefers_amd
    interface
        integer(c_int) function bounded_amd_available() &
            bind(C, name='fsparse_native_amd_bounded_available')
            import :: c_int
        end function bounded_amd_available
        integer(c_int) function bounded_amd_size(n, nnz, ptr, row, low, &
                                                 arena_bytes, output_bytes) &
            bind(C, name='fsparse_native_amd_bounded_size')
            import :: c_int, c_int64_t
            integer(c_int), value :: n, nnz, low
            integer(c_int), intent(in) :: ptr(*), row(*)
            integer(c_int64_t), intent(out) :: arena_bytes, output_bytes
        end function bounded_amd_size
        integer(c_int) function bounded_amd_order(n, nnz, ptr, row, low, &
                                                  arena, &
                                                  arena_bytes, perm, inverse) &
            bind(C, name='fsparse_native_amd_bounded')
            import :: c_int, c_int64_t
            integer(c_int), value :: n, nnz, low
            integer(c_int), intent(in) :: ptr(*), row(*)
            integer(c_int), intent(inout) :: arena(*)
            integer(c_int64_t), value :: arena_bytes
            integer(c_int), intent(out) :: perm(*), inverse(*)
        end function bounded_amd_order
    end interface
contains
    pure logical function default_prefers_amd(ptr, rows) result(prefers)
        ! DEFAULT ordering policy shared by all native CPU factorizations.
        ! Natural order is already near-optimal for small or narrow-band
        ! patterns, so AMD is skipped there. Tolerates malformed CSC; the
        ! factorization validates it.
        integer, intent(in) :: ptr(:), rows(:)
        integer :: j, p
        prefers = .false.
        if (size(ptr) - 1 < 4096) return
        do j = 1, size(ptr) - 1
            do p = max(ptr(j), 1), min(ptr(j + 1) - 1, size(rows))
                if (abs(rows(p) - j) > 12) then
                    prefers = .true.
                    return
                end if
            end do
        end do
    end function default_prefers_amd
end module fortsparse_native_amd
