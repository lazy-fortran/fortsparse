module fortsparse_csc
    use fortsparse_kinds, only: dp
    use fortsparse_status, only: fortsparse_status_t, status_set, &
        FORTSPARSE_OK, FORTSPARSE_INVALID_MATRIX
    implicit none
    private

    ! Compressed-sparse-column matrix, real double values.
    ! Indices are 1-based Fortran indices: col_ptr has size ncol+1,
    ! row_idx and val have size nnz, row_idx values lie in [1, nrow].
    type, public :: csc_t
        integer               :: nrow = 0
        integer               :: ncol = 0
        integer               :: nnz  = 0
        integer,  allocatable :: col_ptr(:)
        integer,  allocatable :: row_idx(:)
        real(dp), allocatable :: val(:)
    end type csc_t

    ! Compressed-sparse-column matrix, complex double values.
    type, public :: csc_z_t
        integer                  :: nrow = 0
        integer                  :: ncol = 0
        integer                  :: nnz  = 0
        integer,     allocatable :: col_ptr(:)
        integer,     allocatable :: row_idx(:)
        complex(dp), allocatable :: val(:)
    end type csc_z_t

    public :: csc_from_triplet
    public :: csc_is_valid
    public :: csc_matvec

    ! Sparse matrix-vector product y = A x for real and complex matrices.
    interface csc_matvec
        module procedure csc_matvec_real
        module procedure csc_matvec_complex
    end interface csc_matvec

    ! COO -> CSC. Sums duplicate (row, col) entries and sorts row indices
    ! ascending within each column. Real and complex value arrays dispatch
    ! through this generic name.
    interface csc_from_triplet
        module procedure csc_from_triplet_real
        module procedure csc_from_triplet_complex
    end interface csc_from_triplet

    ! Structural validity check covering both real and complex matrices.
    interface csc_is_valid
        module procedure csc_is_valid_real
        module procedure csc_is_valid_complex
    end interface csc_is_valid

contains

    ! Build a real CSC matrix from coordinate triplets.
    subroutine csc_from_triplet_real(nrow, ncol, rows, cols, vals, A, status)
        integer,                   intent(in)  :: nrow, ncol
        integer,                   intent(in)  :: rows(:), cols(:)
        real(dp),                  intent(in)  :: vals(:)
        type(csc_t),               intent(out) :: A
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: perm(:)
        integer              :: ntrip, nnz

        call validate_triplet(nrow, ncol, rows, cols, size(vals), status)
        if (status%code /= FORTSPARSE_OK) return

        ntrip = size(rows)
        call column_order(ncol, rows, cols, ntrip, perm)
        A%nrow = nrow
        A%ncol = ncol
        allocate (A%col_ptr(ncol + 1))
        allocate (A%row_idx(ntrip))
        allocate (A%val(ntrip))
        call compress_real(rows, cols, vals, perm, ncol, A%col_ptr, &
            A%row_idx, A%val, nnz)
        A%nnz = nnz
        call shrink_real(A, nnz)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine csc_from_triplet_real

    ! Build a complex CSC matrix from coordinate triplets.
    subroutine csc_from_triplet_complex(nrow, ncol, rows, cols, vals, A, status)
        integer,                   intent(in)  :: nrow, ncol
        integer,                   intent(in)  :: rows(:), cols(:)
        complex(dp),               intent(in)  :: vals(:)
        type(csc_z_t),             intent(out) :: A
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: perm(:)
        integer              :: ntrip, nnz

        call validate_triplet(nrow, ncol, rows, cols, size(vals), status)
        if (status%code /= FORTSPARSE_OK) return

        ntrip = size(rows)
        call column_order(ncol, rows, cols, ntrip, perm)
        A%nrow = nrow
        A%ncol = ncol
        allocate (A%col_ptr(ncol + 1))
        allocate (A%row_idx(ntrip))
        allocate (A%val(ntrip))
        call compress_complex(rows, cols, vals, perm, ncol, A%col_ptr, &
            A%row_idx, A%val, nnz)
        A%nnz = nnz
        call shrink_complex(A, nnz)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine csc_from_triplet_complex

    ! Reject malformed triplet input before any allocation.
    subroutine validate_triplet(nrow, ncol, rows, cols, nval, status)
        integer,                   intent(in)  :: nrow, ncol, nval
        integer,                   intent(in)  :: rows(:), cols(:)
        type(fortsparse_status_t), intent(out) :: status
        integer                                :: k

        if (nrow < 1 .or. ncol < 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "csc_from_triplet: nrow and ncol must be positive")
            return
        end if
        if (size(rows) /= size(cols) .or. size(rows) /= nval) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "csc_from_triplet: rows, cols, vals length mismatch")
            return
        end if
        do k = 1, size(rows)
            if (rows(k) < 1 .or. rows(k) > nrow) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "csc_from_triplet: row index out of range")
                return
            end if
            if (cols(k) < 1 .or. cols(k) > ncol) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "csc_from_triplet: col index out of range")
                return
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_triplet

    ! Permutation ordering triplets by (col, row) ascending. Two stable
    ! counting passes (minor key rows first, major key cols second) leave the
    ! entries ordered by (col, row).
    subroutine column_order(ncol, rows, cols, ntrip, perm)
        integer,              intent(in)  :: ncol, ntrip
        integer,              intent(in)  :: rows(:), cols(:)
        integer, allocatable, intent(out) :: perm(:)

        integer, allocatable :: tmp(:)
        integer              :: maxrow

        allocate (perm(ntrip))
        allocate (tmp(ntrip))
        call iota(perm)
        maxrow = 0
        if (ntrip > 0) maxrow = maxval(rows)
        call counting_sort(perm, tmp, rows, maxrow)
        call counting_sort(perm, tmp, cols, ncol)
    end subroutine column_order

    ! Stable counting sort of perm by key(perm(:)); keys lie in [1, nkey].
    subroutine counting_sort(perm, tmp, key, nkey)
        integer, intent(inout) :: perm(:)
        integer, intent(inout) :: tmp(:)
        integer, intent(in)    :: key(:)
        integer, intent(in)    :: nkey

        integer, allocatable :: count(:)
        integer              :: i, k, pos, n

        n = size(perm)
        if (n == 0) return
        allocate (count(nkey + 1))
        count = 0
        do i = 1, n
            k = key(perm(i))
            count(k + 1) = count(k + 1) + 1
        end do
        do k = 2, nkey + 1
            count(k) = count(k) + count(k - 1)
        end do
        do i = 1, n
            k = key(perm(i))
            pos = count(k) + 1
            tmp(pos) = perm(i)
            count(k) = pos
        end do
        perm = tmp
    end subroutine counting_sort

    ! Compress sorted real triplets, summing duplicate (row, col) pairs.
    subroutine compress_real(rows, cols, vals, perm, ncol, col_ptr, &
            row_idx, val, nnz)
        integer,  intent(in)  :: rows(:), cols(:), perm(:), ncol
        real(dp), intent(in)  :: vals(:)
        integer,  intent(out) :: col_ptr(:), row_idx(:)
        real(dp), intent(out) :: val(:)
        integer,  intent(out) :: nnz

        integer :: i, p, r, c, prev_r, prev_c

        nnz = 0
        col_ptr(1) = 1
        prev_r = 0
        prev_c = 0
        do i = 1, size(perm)
            p = perm(i)
            r = rows(p)
            c = cols(p)
            if (nnz > 0 .and. c == prev_c) then
                if (r == prev_r) then
                    val(nnz) = val(nnz) + vals(p)
                    cycle
                end if
            end if
            call close_columns(col_ptr, prev_c, c, nnz)
            nnz = nnz + 1
            row_idx(nnz) = r
            val(nnz) = vals(p)
            prev_r = r
            prev_c = c
        end do
        call finish_columns(col_ptr, prev_c, ncol, nnz)
    end subroutine compress_real

    ! Compress sorted complex triplets, summing duplicate (row, col) pairs.
    subroutine compress_complex(rows, cols, vals, perm, ncol, col_ptr, &
            row_idx, val, nnz)
        integer,     intent(in)  :: rows(:), cols(:), perm(:), ncol
        complex(dp), intent(in)  :: vals(:)
        integer,     intent(out) :: col_ptr(:), row_idx(:)
        complex(dp), intent(out) :: val(:)
        integer,     intent(out) :: nnz

        integer :: i, p, r, c, prev_r, prev_c

        nnz = 0
        col_ptr(1) = 1
        prev_r = 0
        prev_c = 0
        do i = 1, size(perm)
            p = perm(i)
            r = rows(p)
            c = cols(p)
            if (nnz > 0 .and. c == prev_c) then
                if (r == prev_r) then
                    val(nnz) = val(nnz) + vals(p)
                    cycle
                end if
            end if
            call close_columns(col_ptr, prev_c, c, nnz)
            nnz = nnz + 1
            row_idx(nnz) = r
            val(nnz) = vals(p)
            prev_r = r
            prev_c = c
        end do
        call finish_columns(col_ptr, prev_c, ncol, nnz)
    end subroutine compress_complex

    ! Set col_ptr starts for every column opened between prev_c and c.
    subroutine close_columns(col_ptr, prev_c, c, nnz)
        integer, intent(inout) :: col_ptr(:)
        integer, intent(in)    :: prev_c, c, nnz
        integer                :: j

        do j = prev_c + 1, c
            col_ptr(j) = nnz + 1
        end do
    end subroutine close_columns

    ! Close every remaining column up to ncol with the final pointer.
    subroutine finish_columns(col_ptr, prev_c, ncol, nnz)
        integer, intent(inout) :: col_ptr(:)
        integer, intent(in)    :: prev_c, ncol, nnz
        integer                :: j

        do j = prev_c + 1, ncol + 1
            col_ptr(j) = nnz + 1
        end do
    end subroutine finish_columns

    ! Trim a real matrix's storage arrays to the deduplicated length.
    subroutine shrink_real(A, nnz)
        type(csc_t), intent(inout) :: A
        integer,     intent(in)    :: nnz

        integer,  allocatable :: ri(:)
        real(dp), allocatable :: vv(:)

        if (nnz == size(A%row_idx)) return
        allocate (ri(nnz))
        allocate (vv(nnz))
        ri = A%row_idx(1:nnz)
        vv = A%val(1:nnz)
        call move_alloc(ri, A%row_idx)
        call move_alloc(vv, A%val)
    end subroutine shrink_real

    ! Trim a complex matrix's storage arrays to the deduplicated length.
    subroutine shrink_complex(A, nnz)
        type(csc_z_t), intent(inout) :: A
        integer,       intent(in)    :: nnz

        integer,     allocatable :: ri(:)
        complex(dp), allocatable :: vv(:)

        if (nnz == size(A%row_idx)) return
        allocate (ri(nnz))
        allocate (vv(nnz))
        ri = A%row_idx(1:nnz)
        vv = A%val(1:nnz)
        call move_alloc(ri, A%row_idx)
        call move_alloc(vv, A%val)
    end subroutine shrink_complex

    ! Fill perm with 1, 2, ..., n.
    pure subroutine iota(perm)
        integer, intent(out) :: perm(:)
        integer              :: i

        do i = 1, size(perm)
            perm(i) = i
        end do
    end subroutine iota

    ! Sparse matvec y = A x for a real CSC matrix.
    pure function csc_matvec_real(A, x) result(y)
        type(csc_t), intent(in) :: A
        real(dp),    intent(in) :: x(:)
        real(dp)                :: y(A%nrow)
        integer                 :: j, p

        y = 0.0_dp
        do j = 1, A%ncol
            do p = A%col_ptr(j), A%col_ptr(j + 1) - 1
                y(A%row_idx(p)) = y(A%row_idx(p)) + A%val(p)*x(j)
            end do
        end do
    end function csc_matvec_real

    ! Sparse matvec y = A x for a complex CSC matrix.
    pure function csc_matvec_complex(A, x) result(y)
        type(csc_z_t), intent(in) :: A
        complex(dp),   intent(in) :: x(:)
        complex(dp)               :: y(A%nrow)
        integer                   :: j, p

        y = (0.0_dp, 0.0_dp)
        do j = 1, A%ncol
            do p = A%col_ptr(j), A%col_ptr(j + 1) - 1
                y(A%row_idx(p)) = y(A%row_idx(p)) + A%val(p)*x(j)
            end do
        end do
    end function csc_matvec_complex

    ! Structural validity for a real CSC matrix.
    pure logical function csc_is_valid_real(A) result(ok)
        type(csc_t), intent(in) :: A
        ok = csc_structure_ok(A%nrow, A%ncol, A%nnz, A%col_ptr, A%row_idx)
    end function csc_is_valid_real

    ! Structural validity for a complex CSC matrix.
    pure logical function csc_is_valid_complex(A) result(ok)
        type(csc_z_t), intent(in) :: A
        ok = csc_structure_ok(A%nrow, A%ncol, A%nnz, A%col_ptr, A%row_idx)
    end function csc_is_valid_complex

    ! Shared structure check: col_ptr monotone with the canonical endpoints
    ! and every row index inside [1, nrow].
    pure logical function csc_structure_ok(nrow, ncol, nnz, col_ptr, row_idx) &
            result(ok)
        integer, intent(in) :: nrow, ncol, nnz
        integer, intent(in) :: col_ptr(:), row_idx(:)
        integer             :: j

        ok = .false.
        if (nrow < 1 .or. ncol < 1 .or. nnz < 0) return
        if (size(col_ptr) /= ncol + 1) return
        if (size(row_idx) /= nnz) return
        if (col_ptr(1) /= 1) return
        if (col_ptr(ncol + 1) /= nnz + 1) return
        do j = 1, ncol
            if (col_ptr(j + 1) < col_ptr(j)) return
        end do
        do j = 1, nnz
            if (row_idx(j) < 1 .or. row_idx(j) > nrow) return
        end do
        ok = .true.
    end function csc_structure_ok

end module fortsparse_csc
