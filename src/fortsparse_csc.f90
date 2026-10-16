module fortsparse_csc
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_positive_inf
    use fortsparse_kinds, only: dp
    use fortsparse_status, only: fortsparse_status_t, status_set, &
        FORTSPARSE_OK, FORTSPARSE_INVALID_MATRIX, FORTSPARSE_INDEX_OVERFLOW
    !$ use omp_lib, only: omp_get_thread_num, omp_get_max_threads, omp_get_num_threads
    implicit none
    private

    ! Matrices with fewer stored entries than this are multiplied serially.
    integer, parameter :: MATVEC_MIN_NNZ = 20000
    ! Upper bound on the team size used by the matvec.
    integer, parameter :: MAXT = 512

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

    public :: csc_normalize_symmetric_raw
    public :: csc_from_triplet
    public :: csc_is_valid
    public :: csc_matvec
    public :: csc_matmul
    public :: csc_transpose
    public :: csc_conjugate_transpose

    interface csc_matmul
        module procedure csc_matmul_real
    end interface csc_matmul

    interface csc_transpose
        module procedure csc_transpose_real
    end interface csc_transpose

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

    subroutine csc_normalize_symmetric_raw(n, col_ptr, row_idx, values, rtol, &
                                           max_entry_correction, status)
        ! Explicit caller-owned projection, never an implicit factor tolerance.
        ! Full, sorted, unique square CSC is checked before the first store.
        integer, intent(in) :: n, col_ptr(:), row_idx(:)
        real(dp), intent(inout) :: values(:)
        real(dp), intent(in) :: rtol
        real(dp), intent(out) :: max_entry_correction
        type(fortsparse_status_t), intent(out) :: status
        integer :: nz, column, entry, row, previous, reflected
        real(dp) :: scale, tolerance, correction
        max_entry_correction = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        'symmetric normalization: invalid tolerance or CSC')
        if (.not. ieee_is_finite(rtol)) return
        if (rtol < 0.0_dp .or. rtol > 1.0_dp .or. n < 0) return
        if (n == huge(n) .or. size(values) == huge(n)) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                            'symmetric normalization: pointer extent overflow')
            return
        end if
        nz = size(values)
        if (size(col_ptr) /= n + 1 .or. size(row_idx) /= nz) return
        if (col_ptr(1) /= 1 .or. col_ptr(n + 1) /= nz + 1) return
        ! Validate every pointer before using any pointer as an array bound.
        do column = 1, n
            if (col_ptr(column) < 1 .or. col_ptr(column) > nz + 1) return
            if (col_ptr(column + 1) < col_ptr(column)) return
        end do
        scale = 0.0_dp
        do column = 1, n
            previous = 0
            do entry = col_ptr(column), col_ptr(column + 1) - 1
                row = row_idx(entry)
                if (row <= previous .or. row > n) return
                if (.not. ieee_is_finite(values(entry))) return
                previous = row
                scale = max(scale, abs(values(entry)))
            end do
        end do
        tolerance = rtol*scale
        do column = 1, n
            do entry = col_ptr(column), col_ptr(column + 1) - 1
                row = row_idx(entry)
                if (row == column) cycle
                reflected = find_entry(column, row)
                if (reflected == 0) return
                correction = difference_bound(values(entry), values(reflected))
                max_entry_correction = max(max_entry_correction, correction)
            end do
        end do
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                        'symmetric normalization: asymmetry exceeds tolerance')
        if (max_entry_correction > tolerance) return
        ! Only lower entries are sources, so a store never changes a later source.
        do column = 1, n
            do entry = col_ptr(column), col_ptr(column + 1) - 1
                row = row_idx(entry)
                if (row <= column) cycle
                reflected = find_entry(column, row)
                values(reflected) = values(entry)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, '')
    contains
        integer function find_entry(wanted_row, wanted_column) result(found)
            integer, intent(in) :: wanted_row, wanted_column
            integer :: lo, hi, mid
            found = 0
            lo = col_ptr(wanted_column)
            hi = col_ptr(wanted_column + 1) - 1
            do while (lo <= hi)
                mid = lo + (hi - lo)/2
                if (row_idx(mid) == wanted_row) then
                    found = mid
                    return
                else if (row_idx(mid) < wanted_row) then
                    lo = mid + 1
                else
                    hi = mid - 1
                end if
            end do
        end function find_entry

        real(dp) function difference_bound(a, b) result(bound)
            ! Outward one-ulp bound on subtraction of finite binary64 inputs.
            ! Opposite signs may have an unrepresentable difference; create its
            ! diagnostic infinity without performing an overflowing operation.
            real(dp), intent(in) :: a, b
            real(dp) :: aa, bb, delta
            logical :: opposite
            if (a == b) then
                bound = 0.0_dp
                return
            end if
            aa = abs(a)
            bb = abs(b)
            opposite = (a > 0.0_dp .and. b < 0.0_dp) .or. &
                       (a < 0.0_dp .and. b > 0.0_dp)
            if (opposite) then
                if (aa > huge(a) - bb) then
                    bound = ieee_value(a, ieee_positive_inf)
                    return
                end if
            end if
            delta = abs(a - b)
            if (delta == huge(a)) then
                bound = delta
                if (opposite) bound = ieee_value(a, ieee_positive_inf)
            else
                bound = nearest(delta, 1.0_dp)
            end if
        end function difference_bound
    end subroutine csc_normalize_symmetric_raw

    subroutine csc_conjugate_transpose(A, adjoint_A, status)
        type(csc_z_t), intent(in) :: A
        type(csc_z_t), intent(out) :: adjoint_A
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: counts(:), next(:)
        integer :: column, entry, position, row

        if (.not. csc_is_valid_complex(A)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "csc_conjugate_transpose: invalid input matrix")
            return
        end if
        adjoint_A%nrow = A%ncol
        adjoint_A%ncol = A%nrow
        adjoint_A%nnz = A%nnz
        allocate ( &
            adjoint_A%col_ptr(adjoint_A%ncol + 1), &
            adjoint_A%row_idx(A%nnz), adjoint_A%val(A%nnz), &
            counts(A%nrow), next(A%nrow))
        counts = 0
        do entry = 1, A%nnz
            counts(A%row_idx(entry)) = counts(A%row_idx(entry)) + 1
        end do
        adjoint_A%col_ptr(1) = 1
        do column = 1, adjoint_A%ncol
            adjoint_A%col_ptr(column + 1) = &
                adjoint_A%col_ptr(column) + counts(column)
        end do
        next = adjoint_A%col_ptr(:adjoint_A%ncol)
        do column = 1, A%ncol
            do entry = A%col_ptr(column), A%col_ptr(column + 1) - 1
                row = A%row_idx(entry)
                position = next(row)
                adjoint_A%row_idx(position) = column
                adjoint_A%val(position) = conjg(A%val(entry))
                next(row) = position + 1
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine csc_conjugate_transpose

    subroutine csc_transpose_real(A, transpose_A, status)
        type(csc_t), intent(in) :: A
        type(csc_t), intent(out) :: transpose_A
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: counts(:), next(:)
        integer :: column, entry, position, row

        if (.not. csc_is_valid_real(A)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "csc_transpose: invalid input matrix")
            return
        end if
        transpose_A%nrow = A%ncol
        transpose_A%ncol = A%nrow
        transpose_A%nnz = A%nnz
        allocate( &
            transpose_A%col_ptr(transpose_A%ncol + 1), &
            transpose_A%row_idx(A%nnz), transpose_A%val(A%nnz), &
            counts(A%nrow), next(A%nrow))
        counts = 0
        do entry = 1, A%nnz
            counts(A%row_idx(entry)) = counts(A%row_idx(entry)) + 1
        end do
        transpose_A%col_ptr(1) = 1
        do column = 1, transpose_A%ncol
            transpose_A%col_ptr(column + 1) = &
                transpose_A%col_ptr(column) + counts(column)
        end do
        next = transpose_A%col_ptr(:transpose_A%ncol)
        do column = 1, A%ncol
            do entry = A%col_ptr(column), A%col_ptr(column + 1) - 1
                row = A%row_idx(entry)
                position = next(row)
                transpose_A%row_idx(position) = column
                transpose_A%val(position) = A%val(entry)
                next(row) = position + 1
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine csc_transpose_real

    subroutine csc_matmul_real(A, B, C, status)
        type(csc_t), intent(in) :: A, B
        type(csc_t), intent(out) :: C
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: marker(:), row_list(:)
        real(dp), allocatable :: workspace(:)
        integer :: b_entry, column, entry, list_count, position, row

        if (.not. csc_is_valid_real(A) .or. &
            .not. csc_is_valid_real(B) .or. A%ncol /= B%nrow) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "csc_matmul: invalid or incompatible input matrices")
            return
        end if
        C%nrow = A%nrow
        C%ncol = B%ncol
        C%nnz = 0
        allocate(C%col_ptr(C%ncol + 1))
        allocate(C%row_idx(max(1, A%nnz + B%nnz)))
        allocate(C%val(size(C%row_idx)))
        allocate(marker(A%nrow), row_list(A%nrow), workspace(A%nrow))
        marker = 0
        workspace = 0.0_dp
        C%col_ptr(1) = 1
        do column = 1, B%ncol
            list_count = 0
            do b_entry = B%col_ptr(column), B%col_ptr(column + 1) - 1
                do entry = A%col_ptr(B%row_idx(b_entry)), &
                        A%col_ptr(B%row_idx(b_entry) + 1) - 1
                    row = A%row_idx(entry)
                    if (marker(row) /= column) then
                        marker(row) = column
                        list_count = list_count + 1
                        row_list(list_count) = row
                        workspace(row) = 0.0_dp
                    end if
                    workspace(row) = workspace(row) + &
                        A%val(entry)*B%val(b_entry)
                end do
            end do
            call sort_integer_prefix(row_list, list_count)
            do position = 1, list_count
                row = row_list(position)
                if (workspace(row) == 0.0_dp) cycle
                call ensure_real_capacity(C, C%nnz + 1)
                C%nnz = C%nnz + 1
                C%row_idx(C%nnz) = row
                C%val(C%nnz) = workspace(row)
            end do
            C%col_ptr(column + 1) = C%nnz + 1
        end do
        call shrink_real(C, C%nnz)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine csc_matmul_real

    subroutine ensure_real_capacity(A, required)
        type(csc_t), intent(inout) :: A
        integer, intent(in) :: required

        integer, allocatable :: rows(:)
        real(dp), allocatable :: values(:)
        integer :: capacity

        if (required <= size(A%row_idx)) return
        capacity = max(required, 2*size(A%row_idx))
        allocate(rows(capacity), values(capacity))
        if (A%nnz > 0) then
            rows(:A%nnz) = A%row_idx(:A%nnz)
            values(:A%nnz) = A%val(:A%nnz)
        end if
        call move_alloc(rows, A%row_idx)
        call move_alloc(values, A%val)
    end subroutine ensure_real_capacity

    pure subroutine sort_integer_prefix(values, count)
        integer, intent(inout) :: values(:)
        integer, intent(in) :: count

        integer :: current, index, position

        do index = 2, count
            current = values(index)
            position = index - 1
            do while (position >= 1)
                if (values(position) <= current) exit
                values(position + 1) = values(position)
                position = position - 1
            end do
            values(position + 1) = current
        end do
    end subroutine sort_integer_prefix

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

    ! Sparse matvec y = A x for a real CSC matrix. Columns are split into
    ! nnz-balanced chunks, one per thread; each thread accumulates into a
    ! private copy of y and the copies are summed row-parallel. Results differ
    ! from the serial sum by summation-order roundoff only. Small matrices run
    ! serially (threading cost exceeds the work).
    function csc_matvec_real(A, x) result(y)
        type(csc_t), intent(in) :: A
        real(dp), intent(in) :: x(:)
        real(dp) :: y(A%nrow)
        real(dp), allocatable :: yp(:, :)
        integer :: j, p, i, t, nt, nta, jlo, jhi, nrow, rlo, rhi
        integer :: rmin(MAXT), rmax(MAXT)

        nrow = A%nrow
        y = 0.0_dp
        nt = min(matvec_thread_count(A%nnz), MAXT)
        if (nt <= 1) then
            do j = 1, A%ncol
                do p = A%col_ptr(j), A%col_ptr(j + 1) - 1
                    y(A%row_idx(p)) = y(A%row_idx(p)) + A%val(p)*x(j)
                end do
            end do
            return
        end if
        allocate (yp(nrow, nt))
        !$omp parallel num_threads(nt) default(shared) &
        !$omp private(t, nta, j, p, i, jlo, jhi, rlo, rhi)
        t = 1
        nta = 1
        !$ t = omp_get_thread_num() + 1
        !$ nta = omp_get_num_threads()
        call column_chunk(A%col_ptr, A%ncol, A%nnz, t, nta, jlo, jhi)
        rlo = nrow + 1
        rhi = 0
        do p = A%col_ptr(jlo), A%col_ptr(max(jlo, jhi + 1)) - 1
            rlo = min(rlo, A%row_idx(p))
            rhi = max(rhi, A%row_idx(p))
        end do
        rmin(t) = rlo
        rmax(t) = rhi
        yp(rlo:rhi, t) = 0.0_dp
        do j = jlo, jhi
            do p = A%col_ptr(j), A%col_ptr(j + 1) - 1
                yp(A%row_idx(p), t) = yp(A%row_idx(p), t) + A%val(p)*x(j)
            end do
        end do
        !$omp barrier
        !$omp do schedule(static)
        do i = 1, nrow
            y(i) = 0.0_dp
            do j = 1, nta
                if (i >= rmin(j) .and. i <= rmax(j)) y(i) = y(i) + yp(i, j)
            end do
        end do
        !$omp end do
        !$omp end parallel
    end function csc_matvec_real

    ! Sparse matvec y = A x for a complex CSC matrix; see csc_matvec_real.
    function csc_matvec_complex(A, x) result(y)
        type(csc_z_t), intent(in) :: A
        complex(dp), intent(in) :: x(:)
        complex(dp) :: y(A%nrow)
        complex(dp), allocatable :: yp(:, :)
        integer :: j, p, i, t, nt, nta, jlo, jhi, nrow, rlo, rhi
        integer :: rmin(MAXT), rmax(MAXT)

        nrow = A%nrow
        y = (0.0_dp, 0.0_dp)
        nt = min(matvec_thread_count(A%nnz), MAXT)
        if (nt <= 1) then
            do j = 1, A%ncol
                do p = A%col_ptr(j), A%col_ptr(j + 1) - 1
                    y(A%row_idx(p)) = y(A%row_idx(p)) + A%val(p)*x(j)
                end do
            end do
            return
        end if
        allocate (yp(nrow, nt))
        !$omp parallel num_threads(nt) default(shared) &
        !$omp private(t, nta, j, p, i, jlo, jhi, rlo, rhi)
        t = 1
        nta = 1
        !$ t = omp_get_thread_num() + 1
        !$ nta = omp_get_num_threads()
        call column_chunk(A%col_ptr, A%ncol, A%nnz, t, nta, jlo, jhi)
        rlo = nrow + 1
        rhi = 0
        do p = A%col_ptr(jlo), A%col_ptr(max(jlo, jhi + 1)) - 1
            rlo = min(rlo, A%row_idx(p))
            rhi = max(rhi, A%row_idx(p))
        end do
        rmin(t) = rlo
        rmax(t) = rhi
        yp(rlo:rhi, t) = (0.0_dp, 0.0_dp)
        do j = jlo, jhi
            do p = A%col_ptr(j), A%col_ptr(j + 1) - 1
                yp(A%row_idx(p), t) = yp(A%row_idx(p), t) + A%val(p)*x(j)
            end do
        end do
        !$omp barrier
        !$omp do schedule(static)
        do i = 1, nrow
            y(i) = 0.0_dp
            do j = 1, nta
                if (i >= rmin(j) .and. i <= rmax(j)) y(i) = y(i) + yp(i, j)
            end do
        end do
        !$omp end do
        !$omp end parallel
    end function csc_matvec_complex

    ! Threads used for a matvec with nnz stored entries: serial below
    ! MATVEC_MIN_NNZ, else the OpenMP team size (one thread per
    ! MATVEC_MIN_NNZ/4 entries at most).
    function matvec_thread_count(nnz) result(nt)
        integer, intent(in) :: nnz
        integer             :: nt

        nt = 1
        !$ nt = min(omp_get_max_threads(), max(1, nnz/(MATVEC_MIN_NNZ/4)))
        if (nnz < MATVEC_MIN_NNZ) nt = 1
    end function matvec_thread_count

    ! Columns jlo:jhi of chunk t of nt, balanced by stored entries (bisection
    ! on col_ptr). Empty chunks give jlo > jhi.
    subroutine column_chunk(col_ptr, ncol, nnz, t, nt, jlo, jhi)
        integer, intent(in)  :: col_ptr(:), ncol, nnz, t, nt
        integer, intent(out) :: jlo, jhi

        jlo = first_column_at(col_ptr, ncol, int((int(t - 1, 8)*nnz)/nt) + 1)
        jhi = first_column_at(col_ptr, ncol, int((int(t, 8)*nnz)/nt) + 1) - 1
        if (t == nt) jhi = ncol
    end subroutine column_chunk

    ! Smallest column j in [1, ncol+1] whose entries reach position target
    ! (ncol+1 if none).
    pure function first_column_at(col_ptr, ncol, target) result(j)
        integer, intent(in) :: col_ptr(:), ncol, target
        integer             :: j, lo, hi, mid

        lo = 1
        hi = ncol + 1
        do while (lo < hi)
            mid = (lo + hi)/2
            if (col_ptr(mid + 1) > col_ptr(1) + target - 1) then
                hi = mid
            else
                lo = mid + 1
            end if
        end do
        j = lo
    end function first_column_at

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
