module fortsparse_triangle
    ! One-triangle symmetric input. A caller that stores only the lower or the
    ! upper triangle of a symmetric matrix (the uplo convention of LAPACK; the
    ! lower triangle of MUMPS and HSL) is expanded here into the full CSC
    ! matrix that the backends factor. The unreferenced triangle is never
    ! read, so the expanded matrix is exactly symmetric by construction, with
    ! ascending rows in every column.
    !
    ! The expansion remembers the pattern it was built for and the source
    ! position of every expanded entry. A later matrix with the same pattern
    ! (a refactor) only gathers its values, so the pattern of the expanded
    ! matrix, which the retained symbolic analysis is compared with, does not
    ! change and the analysis is reused.
    use fortsparse_kinds, only: dp
    use fortsparse_csc, only: csc_t
    use fortsparse_status, only: fortsparse_status_t, status_set, FORTSPARSE_OK, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_INDEX_OVERFLOW
    implicit none
    private

    integer, parameter, public :: TRIANGLE_LOWER = 1
    integer, parameter, public :: TRIANGLE_UPPER = 2

    type, public :: triangle_cache_t
        integer :: n = -1, nz = -1, uplo = 0
        ! pattern of the stored triangle the expansion was built for
        integer, allocatable :: src_ptr(:), src_row(:)
        ! position in the stored triangle of every expanded entry
        integer, allocatable :: source(:)
        ! the expanded full matrix; val is refreshed by every expansion
        type(csc_t) :: full
    end type triangle_cache_t

    public :: triangle_expand
contains

    ! Expand the stored triangle (col_ptr, row_idx, val) of an nrow x ncol
    ! symmetric matrix into cache%full. Entries in the other triangle are
    ! ignored. reused is .true. when the cached pattern was valid.
    subroutine triangle_expand(cache, nrow, ncol, nz, col_ptr, row_idx, val, uplo, &
            reused, status)
        type(triangle_cache_t), intent(inout) :: cache
        integer, intent(in) :: nrow, ncol, nz, col_ptr(:), row_idx(:), uplo
        real(dp), intent(in) :: val(:)
        logical, intent(out) :: reused
        type(fortsparse_status_t), intent(out) :: status
        integer :: j, p

        reused = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            'one-triangle input: square CSC with valid 1-based indices required')
        if (nrow /= ncol .or. ncol < 1 .or. nz < 0) return
        if (size(col_ptr) < ncol + 1 .or. size(row_idx) < nz .or. size(val) < nz) return
        if (col_ptr(1) /= 1 .or. col_ptr(ncol + 1) /= nz + 1) return
        do j = 1, ncol
            if (col_ptr(j + 1) < col_ptr(j)) return
        end do
        do p = 1, nz
            if (row_idx(p) < 1 .or. row_idx(p) > nrow) return
        end do
        call status_set(status, FORTSPARSE_OK, '')
        if (cache%uplo == uplo .and. cache%n == ncol .and. cache%nz == nz) then
            if (all(cache%src_ptr == col_ptr(1:ncol + 1))) then
                reused = all(cache%src_row == row_idx(1:nz))
            end if
        end if
        if (.not. reused) then
            call build(cache, ncol, nz, col_ptr, row_idx, uplo, status)
            if (status%code /= FORTSPARSE_OK) return
        end if
        !$omp parallel do schedule(static) if (cache%full%nnz > 100000)
        do p = 1, cache%full%nnz
            cache%full%val(p) = val(cache%source(p))
        end do
        !$omp end parallel do
    end subroutine triangle_expand

    subroutine build(cache, n, nz, col_ptr, row_idx, uplo, status)
        type(triangle_cache_t), intent(inout) :: cache
        integer, intent(in) :: n, nz, col_ptr(:), row_idx(:), uplo
        type(fortsparse_status_t), intent(out) :: status
        integer, allocatable :: cnt(:), fill(:)
        integer(kind=8) :: total
        integer :: j, p, i, q, nfull
        logical :: keep

        cache%n = -1
        allocate (cnt(n + 1), fill(n))
        cnt = 0
        total = 0
        do j = 1, n
            do p = col_ptr(j), col_ptr(j + 1) - 1
                i = row_idx(p)
                if (.not. kept(i, j, uplo)) cycle
                cnt(j) = cnt(j) + 1
                total = total + 1
                if (i /= j) then
                    cnt(i) = cnt(i) + 1
                    total = total + 1
                end if
            end do
        end do
        if (total > huge(0) - 1) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                'one-triangle input: expanded matrix exceeds 32-bit indices')
            return
        end if
        nfull = int(total)
        if (allocated(cache%full%col_ptr)) deallocate (cache%full%col_ptr)
        if (allocated(cache%full%row_idx)) deallocate (cache%full%row_idx)
        if (allocated(cache%full%val)) deallocate (cache%full%val)
        if (allocated(cache%source)) deallocate (cache%source)
        allocate (cache%full%col_ptr(n + 1), cache%full%row_idx(nfull), &
            cache%full%val(nfull), cache%source(nfull))
        cache%full%col_ptr(1) = 1
        do j = 1, n
            cache%full%col_ptr(j + 1) = cache%full%col_ptr(j) + cnt(j)
        end do
        fill = cache%full%col_ptr(1:n)
        ! Columns in ascending order: the mirrored entries that a column
        ! receives from earlier columns, then its own, give ascending rows
        ! for the lower triangle; own entries first, then mirrors from later
        ! columns, for the upper one.
        do j = 1, n
            do p = col_ptr(j), col_ptr(j + 1) - 1
                i = row_idx(p)
                if (.not. kept(i, j, uplo)) cycle
                q = fill(j)
                cache%full%row_idx(q) = i
                cache%source(q) = p
                fill(j) = q + 1
                if (i /= j) then
                    q = fill(i)
                    cache%full%row_idx(q) = j
                    cache%source(q) = p
                    fill(i) = q + 1
                end if
            end do
        end do
        cache%full%nrow = n
        cache%full%ncol = n
        cache%full%nnz = nfull
        call sort_columns(cache%full, cache%source, keep)
        if (.not. keep) then
            cache%full%nnz = 0
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'one-triangle input: duplicate entry in the stored triangle')
            return
        end if
        cache%src_ptr = col_ptr(1:n + 1)
        cache%src_row = row_idx(1:nz)
        cache%n = n
        cache%nz = nz
        cache%uplo = uplo
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine build

    pure logical function kept(i, j, uplo)
        integer, intent(in) :: i, j, uplo
        if (uplo == TRIANGLE_LOWER) then
            kept = i >= j
        else
            kept = i <= j
        end if
    end function kept

    ! Ascending rows in every column (insertion sort of unsorted input);
    ! ok is .false. on a repeated row.
    subroutine sort_columns(m, source, ok)
        type(csc_t), intent(inout) :: m
        integer, intent(inout) :: source(:)
        logical, intent(out) :: ok
        integer :: j, p, q, r, s
        ok = .true.
        do j = 1, m%ncol
            do p = m%col_ptr(j) + 1, m%col_ptr(j + 1) - 1
                if (m%row_idx(p) > m%row_idx(p - 1)) cycle
                r = m%row_idx(p)
                s = source(p)
                q = p - 1
                do while (q >= m%col_ptr(j))
                    if (m%row_idx(q) <= r) exit
                    m%row_idx(q + 1) = m%row_idx(q)
                    source(q + 1) = source(q)
                    q = q - 1
                end do
                m%row_idx(q + 1) = r
                source(q + 1) = s
            end do
            do p = m%col_ptr(j) + 1, m%col_ptr(j + 1) - 1
                if (m%row_idx(p) == m%row_idx(p - 1)) ok = .false.
            end do
        end do
    end subroutine sort_columns
end module fortsparse_triangle
