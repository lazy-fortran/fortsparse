module fortsparse_umfpack_ipc
    ! Out-of-process UMFPACK backend (MIT). The library never links UMFPACK;
    ! it drives a separate GPL helper binary over shared memory and named
    ! semaphores through the MIT IPC shim (src/ipc/fsparse_ipc.c). This module
    ! holds no UMFPACK symbol. It discovers the helper from the
    ! FORTSPARSE_UMFPACK_HELPER environment variable or from PATH; when no
    ! helper is present it reports FORTSPARSE_BACKEND_UNAVAILABLE, which is the
    ! expected, testable state under the MIT-only fpm/fo build.
    use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, &
        c_char, c_int, c_int32_t, c_int64_t, c_double, c_f_pointer, c_loc, &
        c_null_char, c_size_t
    use fortsparse_kinds, only: dp, i8
    use fortsparse_backend, only: sparse_backend_t
    use fortsparse_csc, only: csc_t, csc_z_t
    use fortsparse_status, only: fortsparse_status_t, status_set, &
        FORTSPARSE_OK, FORTSPARSE_SINGULAR, FORTSPARSE_INTERNAL_ERROR, &
        FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_INVALID_MATRIX, &
        FORTSPARSE_NOT_FACTORED
    use fortsparse_ipc_proto, only: shm_header_t, &
        OP_FACTOR_REAL, OP_FACTOR_COMPLEX, OP_SOLVE_REAL, OP_SOLVE_COMPLEX, &
        ST_OK, ST_SINGULAR, find_helper
    implicit none
    private

    public :: umfpack_ipc_backend_t

    ! Out-of-process UMFPACK backend. Owns one helper session and the shared
    ! mapping while a factorization is live; the resident factorization lives in
    ! the helper process. sparse_free releases the whole session (helper and
    ! mapping): the mapping is matrix-sized, so an idle retained session would
    ! raise the caller's peak footprint for the price of one saved process
    ! spawn. The FINAL shuts a still-live helper down when the backend is
    ! deallocated, so a solver needs no explicit teardown.
    type, extends(sparse_backend_t) :: umfpack_ipc_backend_t
        type(c_ptr) :: sess = c_null_ptr
        integer(i8) :: mapped = 0_i8
        integer(i8) :: o_pool = 0_i8 ! byte offset of the vector pool
        integer     :: pool_next = 0 ! next free pool slot
        integer     :: n = 0
        logical     :: is_complex = .false.
    contains
        procedure :: factor_real => umf_factor_real
        procedure :: factor_complex => umf_factor_complex
        procedure :: factor_real_raw => umf_factor_real_raw
        procedure :: factor_complex_raw => umf_factor_complex_raw
        procedure :: solve_real => umf_solve_real
        procedure :: solve_complex => umf_solve_complex
        procedure :: solve_real_inplace => umf_solve_real_inplace
        procedure :: solve_complex_inplace => umf_solve_complex_inplace
        procedure :: vector => umf_vector
        procedure :: free => umf_free
        final :: umf_final
    end type umfpack_ipc_backend_t

    ! Slots reserved in the shared mapping for zero-copy solve vectors. One per
    ! concurrently-live solve vector a client holds; eight covers typical use.
    integer, parameter :: POOL_SLOTS = 8

    interface

        function fsparse_ipc_header_bytes() &
                bind(c, name="fsparse_ipc_header_bytes") result(b)
            import :: c_int64_t
            integer(c_int64_t) :: b
        end function fsparse_ipc_header_bytes

        function fsparse_ipc_start(helper, bytes, err) &
                bind(c, name="fsparse_ipc_start") result(sess)
            import :: c_ptr, c_char, c_int64_t, c_int
            character(kind=c_char), intent(in) :: helper(*)
            integer(c_int64_t),     value      :: bytes
            integer(c_int),         intent(out):: err
            type(c_ptr)                         :: sess
        end function fsparse_ipc_start

        function fsparse_ipc_data(sess) bind(c, name="fsparse_ipc_data") &
                result(p)
            import :: c_ptr
            type(c_ptr), value :: sess
            type(c_ptr)        :: p
        end function fsparse_ipc_data

        function fsparse_ipc_offset(sess, ptr) &
                bind(c, name="fsparse_ipc_offset") result(off)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: sess
            type(c_ptr), value :: ptr
            integer(c_int64_t) :: off
        end function fsparse_ipc_offset

        integer(c_int) function fsparse_ipc_call(sess) &
                bind(c, name="fsparse_ipc_call")
            import :: c_ptr, c_int
            type(c_ptr), value :: sess
        end function fsparse_ipc_call

        subroutine fsparse_ipc_stop(sess) bind(c, name="fsparse_ipc_stop")
            import :: c_ptr
            type(c_ptr), value :: sess
        end subroutine fsparse_ipc_stop

    end interface

contains

    ! Factor a real matrix in the helper.
    subroutine umf_factor_real(self, A, refine, status)
        class(umfpack_ipc_backend_t), intent(inout) :: self
        type(csc_t),                  intent(in)    :: A
        logical,                      intent(in)    :: refine
        type(fortsparse_status_t),    intent(out)   :: status

        call umf_factor_real_raw(self, A%nrow, A%ncol, A%nnz, A%col_ptr, &
            A%row_idx, A%val, refine, status)
    end subroutine umf_factor_real

    ! Factor a complex matrix in the helper.
    subroutine umf_factor_complex(self, A, refine, status)
        class(umfpack_ipc_backend_t), intent(inout) :: self
        type(csc_z_t),                intent(in)    :: A
        logical,                      intent(in)    :: refine
        type(fortsparse_status_t),    intent(out)   :: status

        call umf_factor_complex_raw(self, A%nrow, A%ncol, A%nnz, A%col_ptr, &
            A%row_idx, A%val, refine, status)
    end subroutine umf_factor_complex

    ! Factor a real matrix straight from caller-owned CSC arrays. Lays out the
    ! shared region, streams the arrays into it (the only copy the parent
    ! makes), then issues FACTOR_REAL. UMFPACK factors square systems only.
    subroutine umf_factor_real_raw(self, nrow, ncol, nz, col_ptr, row_idx, &
            val, refine, status)
        class(umfpack_ipc_backend_t), intent(inout) :: self
        integer,                      intent(in)    :: nrow, ncol, nz
        integer,                      intent(in)    :: col_ptr(:), row_idx(:)
        real(dp),                     intent(in)    :: val(:)
        logical,                      intent(in)    :: refine
        type(fortsparse_status_t),    intent(out)   :: status

        type(shm_header_t), pointer :: h
        integer(i8) :: o_cp, o_ri, o_ax, o_b, o_x, o_pool, total

        if (nrow /= ncol) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "umfpack ipc: matrix must be square")
            return
        end if
        call layout_real(ncol, nz, o_b, o_x, o_pool, o_cp, o_ri, o_ax, total)
        call ensure_session(self, total, status)
        if (status%code /= FORTSPARSE_OK) return
        call header_of(self%sess, h)
        h%n = int(ncol, c_int64_t)
        h%nnz = int(nz, c_int64_t)
        h%refine = merge(1_c_int32_t, 0_c_int32_t, refine)
        h%is_complex = 0_c_int32_t
        call set_offsets(h, o_cp, o_ri, o_ax, 0_i8, o_b, 0_i8, o_x, 0_i8, total)
        call write_index(self%sess, o_cp, col_ptr, ncol + 1)
        call write_index(self%sess, o_ri, row_idx, nz)
        call write_real(self%sess, o_ax, val(1:nz))
        self%n = ncol
        self%o_pool = o_pool
        self%pool_next = 0
        self%is_complex = .false.
        call run(self, OP_FACTOR_REAL, status)
    end subroutine umf_factor_real_raw

    ! Factor a complex matrix straight from caller-owned CSC arrays, using
    ! split real/imag value buffers in the shared region.
    subroutine umf_factor_complex_raw(self, nrow, ncol, nz, col_ptr, row_idx, &
            val, refine, status)
        class(umfpack_ipc_backend_t), intent(inout) :: self
        integer,                      intent(in)    :: nrow, ncol, nz
        integer,                      intent(in)    :: col_ptr(:), row_idx(:)
        complex(dp),                  intent(in)    :: val(:)
        logical,                      intent(in)    :: refine
        type(fortsparse_status_t),    intent(out)   :: status

        type(shm_header_t), pointer :: h
        integer(i8) :: o_cp, o_ri, o_ax, o_az, o_b, o_bz, o_x, o_xz, total

        if (nrow /= ncol) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "umfpack ipc: matrix must be square")
            return
        end if
        call layout_complex(ncol, nz, o_cp, o_ri, o_ax, o_az, o_b, o_bz, &
            o_x, o_xz, total)
        call ensure_session(self, total, status)
        if (status%code /= FORTSPARSE_OK) return
        call header_of(self%sess, h)
        h%n = int(ncol, c_int64_t)
        h%nnz = int(nz, c_int64_t)
        h%refine = merge(1_c_int32_t, 0_c_int32_t, refine)
        h%is_complex = 1_c_int32_t
        call set_offsets(h, o_cp, o_ri, o_ax, o_az, o_b, o_bz, o_x, o_xz, total)
        call write_index(self%sess, o_cp, col_ptr, ncol + 1)
        call write_index(self%sess, o_ri, row_idx, nz)
        call write_split(self%sess, o_ax, o_az, val(1:nz))
        self%n = ncol
        self%is_complex = .true.
        call run(self, OP_FACTOR_COMPLEX, status)
    end subroutine umf_factor_complex_raw

    ! Solve a real RHS through the resident factorization.
    subroutine umf_solve_real(self, b, x, status)
        class(umfpack_ipc_backend_t), intent(inout) :: self
        real(dp), target, contiguous, intent(in)    :: b(:)
        real(dp), target, contiguous, intent(out)   :: x(:)
        type(fortsparse_status_t),    intent(out)   :: status

        type(shm_header_t), pointer :: h
        integer(i8)                 :: o_b_fixed, o_x_fixed
        integer(c_int64_t)          :: ob, ox

        if (no_session(self, status)) return
        call header_of(self%sess, h)
        o_b_fixed = fsparse_ipc_header_bytes()
        o_x_fixed = o_b_fixed + int(self%n, i8)*8_i8
        ! RHS: when b is a fortsparse vector it already lives in the mapping, so
        ! point the solve at it; a plain array is copied into the fixed RHS slot.
        ob = fsparse_ipc_offset(self%sess, c_loc(b(1)))
        if (ob >= 0_c_int64_t) then
            h%off_b = ob
        else
            h%off_b = int(o_b_fixed, c_int64_t)
            call write_real(self%sess, o_b_fixed, b)
        end if
        ! Solution: write straight into x's slot when x is a fortsparse vector;
        ! otherwise the fixed slot, copied out after the solve.
        ox = fsparse_ipc_offset(self%sess, c_loc(x(1)))
        if (ox >= 0_c_int64_t) then
            h%off_x = ox
        else
            h%off_x = int(o_x_fixed, c_int64_t)
        end if
        call run(self, OP_SOLVE_REAL, status)
        if (status%code /= FORTSPARSE_OK) return
        if (ox < 0_c_int64_t) call read_real(self%sess, o_x_fixed, x)
    end subroutine umf_solve_real

    ! Solve a complex RHS through the resident factorization.
    subroutine umf_solve_complex(self, b, x, status)
        class(umfpack_ipc_backend_t), intent(inout) :: self
        complex(dp),                  intent(in)    :: b(:)
        complex(dp),                  intent(out)   :: x(:)
        type(fortsparse_status_t),    intent(out)   :: status

        type(shm_header_t), pointer :: h

        if (no_session(self, status)) return
        call header_of(self%sess, h)
        call write_split(self%sess, int(h%off_b, i8), int(h%off_bz, i8), b)
        call run(self, OP_SOLVE_COMPLEX, status)
        if (status%code /= FORTSPARSE_OK) return
        call read_split(self%sess, int(h%off_x, i8), int(h%off_xz, i8), x)
    end subroutine umf_solve_complex

    ! In-place real solve: write the RHS from b, then read the solution back
    ! into b. Saves the caller a temporary and the final copy a separate-output
    ! solve forces; b crosses the boundary once each way, the irreducible cost.
    subroutine umf_solve_real_inplace(self, b, status)
        class(umfpack_ipc_backend_t), intent(inout) :: self
        real(dp),                     intent(inout) :: b(:)
        type(fortsparse_status_t),    intent(out)   :: status

        type(shm_header_t), pointer :: h
        integer(i8)                 :: o_b_fixed, o_x_fixed

        if (no_session(self, status)) return
        call header_of(self%sess, h)
        o_b_fixed = fsparse_ipc_header_bytes()
        o_x_fixed = o_b_fixed + int(self%n, i8)*8_i8
        h%off_b = int(o_b_fixed, c_int64_t)
        h%off_x = int(o_x_fixed, c_int64_t)
        call write_real(self%sess, o_b_fixed, b)
        call run(self, OP_SOLVE_REAL, status)
        if (status%code /= FORTSPARSE_OK) return
        call read_real(self%sess, o_x_fixed, b)
    end subroutine umf_solve_real_inplace

    ! In-place complex solve; b is the RHS on entry, the solution on return.
    subroutine umf_solve_complex_inplace(self, b, status)
        class(umfpack_ipc_backend_t), intent(inout) :: self
        complex(dp),                  intent(inout) :: b(:)
        type(fortsparse_status_t),    intent(out)   :: status

        type(shm_header_t), pointer :: h

        if (no_session(self, status)) return
        call header_of(self%sess, h)
        call write_split(self%sess, int(h%off_b, i8), int(h%off_bz, i8), b)
        call run(self, OP_SOLVE_COMPLEX, status)
        if (status%code /= FORTSPARSE_OK) return
        call read_split(self%sess, int(h%off_x, i8), int(h%off_xz, i8), b)
    end subroutine umf_solve_complex_inplace

    ! Hand out the next free pool slot as a length-n real array aliasing the
    ! shared mapping. A solve whose RHS and solution are such vectors crosses the
    ! process boundary with no copy. Returns null if there is no factorization
    ! yet, the size does not match it, or the pool is exhausted; valid until the
    ! next factor (which resets the pool) or the backend's teardown.
    function umf_vector(self, n) result(p)
        class(umfpack_ipc_backend_t), intent(inout) :: self
        integer,                      intent(in)    :: n
        real(dp), pointer                           :: p(:)

        integer(i8) :: off

        p => null()
        if (.not. c_associated(self%sess)) return
        if (n /= self%n .or. self%pool_next >= POOL_SLOTS) return
        off = self%o_pool + int(self%pool_next, i8)*int(n, i8)*8_i8
        call c_f_pointer(region_at(self%sess, off), p, [n])
        self%pool_next = self%pool_next + 1
    end function umf_vector

    ! Release the resident factorization together with its session: the helper
    ! exits and the matrix-sized shared mapping is unmapped, so a freed solver
    ! holds no memory. The next factor spawns a fresh helper; that spawn is
    ! milliseconds against factorizations worth keeping out-of-process.
    subroutine umf_free(self)
        class(umfpack_ipc_backend_t), intent(inout) :: self

        if (c_associated(self%sess)) call fsparse_ipc_stop(self%sess)
        self%sess = c_null_ptr
        self%mapped = 0_i8
        self%o_pool = 0_i8
        self%pool_next = 0
        self%n = 0
        self%is_complex = .false.
    end subroutine umf_free

    ! Shut the helper down when the backend is destroyed. Runs automatically when
    ! a solver (and its allocatable backend) goes out of scope, so client code
    ! needs no explicit teardown.
    subroutine umf_final(self)
        type(umfpack_ipc_backend_t), intent(inout) :: self

        if (c_associated(self%sess)) call fsparse_ipc_stop(self%sess)
        self%sess = c_null_ptr
        self%mapped = 0_i8
    end subroutine umf_final

    ! Ensure a helper session whose mapping holds at least `total` bytes. Reuses
    ! the resident session when it already fits, so a steady problem size spawns
    ! the helper once; a larger matrix respawns it bigger. A quarter of headroom
    ! absorbs the small per-factorization variation in nonzero count. On a missing
    ! helper this reports FORTSPARSE_BACKEND_UNAVAILABLE, the expected state.
    subroutine ensure_session(self, total, status)
        class(umfpack_ipc_backend_t), intent(inout) :: self
        integer(i8),                  intent(in)    :: total
        type(fortsparse_status_t),    intent(out)   :: status

        character(:), allocatable :: helper
        integer(i8)               :: want
        integer                   :: err

        if (c_associated(self%sess) .and. total <= self%mapped) then
            call status_set(status, FORTSPARSE_OK, "")
            return
        end if
        if (c_associated(self%sess)) then
            call fsparse_ipc_stop(self%sess)
            self%sess = c_null_ptr
            self%mapped = 0_i8
        end if
        helper = find_helper()
        if (len(helper) == 0) then
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                "umfpack ipc: helper not found via "// &
                "FORTSPARSE_UMFPACK_HELPER or PATH")
            return
        end if
        want = total + total/4_i8
        self%sess = fsparse_ipc_start(helper//c_null_char, &
            int(want, c_int64_t), err)
        if (err /= 0 .or. .not. c_associated(self%sess)) then
            self%sess = c_null_ptr
            self%mapped = 0_i8
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                "umfpack ipc: failed to start helper "//helper)
            return
        end if
        self%mapped = want
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine ensure_session

    ! True (with status set) when no helper session is live, i.e. before any
    ! factor or after free. The solver's factored flag already rejects these
    ! calls at the public API; this guard keeps a direct backend call from
    ! dereferencing a torn-down session.
    logical function no_session(self, status)
        class(umfpack_ipc_backend_t), intent(in)  :: self
        type(fortsparse_status_t),    intent(out) :: status

        no_session = .not. c_associated(self%sess)
        if (no_session) then
            call status_set(status, FORTSPARSE_NOT_FACTORED, &
                "umfpack ipc: no live factorization")
        else
            call status_set(status, FORTSPARSE_OK, "")
        end if
    end function no_session

    ! Post an opcode and map the helper's normalized status onto fortsparse.
    subroutine run(self, opcode, status)
        class(umfpack_ipc_backend_t), intent(inout) :: self
        integer(c_int32_t),           intent(in)    :: opcode
        type(fortsparse_status_t),    intent(out)   :: status

        type(shm_header_t), pointer :: h
        integer(c_int)              :: st

        call header_of(self%sess, h)
        h%opcode = opcode
        st = fsparse_ipc_call(self%sess)
        if (st == ST_OK) then
            call status_set(status, FORTSPARSE_OK, "")
        else if (st == ST_SINGULAR) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "umfpack ipc: matrix is singular")
        else
            call status_set(status, FORTSPARSE_INTERNAL_ERROR, &
                "umfpack ipc: helper reported an error")
        end if
    end subroutine run

    ! Bind the protocol header at offset 0 of the mapped region.
    subroutine header_of(sess, h)
        type(c_ptr),                 intent(in)  :: sess
        type(shm_header_t), pointer, intent(out) :: h

        call c_f_pointer(fsparse_ipc_data(sess), h)
    end subroutine header_of

    ! Real layout: header, b(n), x(n), pool(POOL_SLOTS*n), then the matrix
    ! colptr(n+1), rowidx(nnz), ax(nnz). The fixed b/x slots and the vector pool
    ! sit before the matrix so their offsets do not shift when the nonzero count
    ! varies, keeping handed-out vector pointers valid across factorizations of
    ! the same size.
    subroutine layout_real(n, nnz, o_b, o_x, o_pool, o_cp, o_ri, o_ax, total)
        integer,     intent(in)  :: n, nnz
        integer(i8), intent(out) :: o_b, o_x, o_pool, o_cp, o_ri, o_ax, total

        o_b = fsparse_ipc_header_bytes()
        o_x = o_b + int(n, i8)*8_i8
        o_pool = o_x + int(n, i8)*8_i8
        o_cp = o_pool + int(POOL_SLOTS, i8)*int(n, i8)*8_i8
        o_ri = o_cp + int(n + 1, i8)*8_i8
        o_ax = o_ri + int(nnz, i8)*8_i8
        total = o_ax + int(nnz, i8)*8_i8
    end subroutine layout_real

    ! Complex layout adds the imaginary arrays az, bz, xz after their reals.
    subroutine layout_complex(n, nnz, o_cp, o_ri, o_ax, o_az, o_b, o_bz, &
            o_x, o_xz, total)
        integer,     intent(in)  :: n, nnz
        integer(i8), intent(out) :: o_cp, o_ri, o_ax, o_az, o_b, o_bz
        integer(i8), intent(out) :: o_x, o_xz, total

        o_cp = fsparse_ipc_header_bytes()
        o_ri = o_cp + int(n + 1, i8)*8_i8
        o_ax = o_ri + int(nnz, i8)*8_i8
        o_az = o_ax + int(nnz, i8)*8_i8
        o_b = o_az + int(nnz, i8)*8_i8
        o_bz = o_b + int(n, i8)*8_i8
        o_x = o_bz + int(n, i8)*8_i8
        o_xz = o_x + int(n, i8)*8_i8
        total = o_xz + int(n, i8)*8_i8
    end subroutine layout_complex

    ! Record offsets and the region size in the header.
    subroutine set_offsets(h, o_cp, o_ri, o_ax, o_az, o_b, o_bz, o_x, o_xz, &
            total)
        type(shm_header_t), intent(inout) :: h
        integer(i8),        intent(in)    :: o_cp, o_ri, o_ax, o_az
        integer(i8),        intent(in)    :: o_b, o_bz, o_x, o_xz, total

        h%off_colptr = int(o_cp, c_int64_t)
        h%off_rowidx = int(o_ri, c_int64_t)
        h%off_ax = int(o_ax, c_int64_t)
        h%off_az = int(o_az, c_int64_t)
        h%off_b = int(o_b, c_int64_t)
        h%off_bz = int(o_bz, c_int64_t)
        h%off_x = int(o_x, c_int64_t)
        h%off_xz = int(o_xz, c_int64_t)
        h%data_bytes = int(total, c_int64_t)
    end subroutine set_offsets

    ! Write the first n entries of a 1-based integer index array into the
    ! region as 0-based int64. The explicit loop converts element by element;
    ! an array expression here would give the compiler license to build an
    ! n-element temporary, which for a large matrix is hundreds of megabytes.
    subroutine write_index(sess, off, idx, n)
        type(c_ptr), intent(in) :: sess
        integer(i8), intent(in) :: off
        integer,     intent(in) :: idx(:)
        integer,     intent(in) :: n

        integer(c_int64_t), pointer :: dst(:)
        integer                     :: i

        call c_f_pointer(region_at(sess, off), dst, [n])
        do i = 1, n
            dst(i) = int(idx(i), c_int64_t) - 1_c_int64_t
        end do
    end subroutine write_index

    ! Copy a real array into the region.
    subroutine write_real(sess, off, v)
        type(c_ptr), intent(in) :: sess
        integer(i8), intent(in) :: off
        real(dp),    intent(in) :: v(:)

        real(c_double), pointer :: dst(:)
        integer                 :: i

        call c_f_pointer(region_at(sess, off), dst, [size(v)])
        do i = 1, size(v)
            dst(i) = v(i)
        end do
    end subroutine write_real

    ! Read a real array out of the region.
    subroutine read_real(sess, off, v)
        type(c_ptr), intent(in)  :: sess
        integer(i8), intent(in)  :: off
        real(dp),    intent(out) :: v(:)

        real(c_double), pointer :: src(:)
        integer                  :: i

        call c_f_pointer(region_at(sess, off), src, [size(v)])
        do i = 1, size(v)
            v(i) = src(i)
        end do
    end subroutine read_real

    ! Split a complex array into separate real and imaginary buffers, element
    ! by element so no value-sized temporary is materialized.
    subroutine write_split(sess, off_re, off_im, z)
        type(c_ptr), intent(in) :: sess
        integer(i8), intent(in) :: off_re, off_im
        complex(dp), intent(in) :: z(:)

        real(c_double), pointer :: re(:), im(:)
        integer                 :: i

        call c_f_pointer(region_at(sess, off_re), re, [size(z)])
        call c_f_pointer(region_at(sess, off_im), im, [size(z)])
        do i = 1, size(z)
            re(i) = real(z(i), dp)
            im(i) = aimag(z(i))
        end do
    end subroutine write_split

    ! Reassemble a complex array from separate real and imaginary buffers.
    subroutine read_split(sess, off_re, off_im, z)
        type(c_ptr), intent(in)  :: sess
        integer(i8), intent(in)  :: off_re, off_im
        complex(dp), intent(out) :: z(:)

        real(c_double), pointer :: re(:), im(:)
        integer                  :: i

        call c_f_pointer(region_at(sess, off_re), re, [size(z)])
        call c_f_pointer(region_at(sess, off_im), im, [size(z)])
        do i = 1, size(z)
            z(i) = cmplx(re(i), im(i), dp)
        end do
    end subroutine read_split

    ! C pointer to the mapped region advanced by off bytes.
    function region_at(sess, off) result(p)
        use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_loc
        type(c_ptr), intent(in) :: sess
        integer(i8), intent(in) :: off
        type(c_ptr)             :: p

        character(kind=c_char), pointer :: base(:)

        call c_f_pointer(fsparse_ipc_data(sess), base, [off + 1_i8])
        p = c_loc(base(off + 1_i8))
    end function region_at

end module fortsparse_umfpack_ipc
