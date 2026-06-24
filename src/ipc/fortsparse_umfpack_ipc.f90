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
        FORTSPARSE_BACKEND_UNAVAILABLE
    use fortsparse_ipc_proto, only: shm_header_t, &
        OP_FACTOR_REAL, OP_FACTOR_COMPLEX, OP_SOLVE_REAL, OP_SOLVE_COMPLEX, &
        ST_OK, ST_SINGULAR, find_helper
    implicit none
    private

    public :: umfpack_ipc_backend_t

    ! Out-of-process UMFPACK backend. Owns one helper session and the shared
    ! mapping; the resident factorization lives in the helper process.
    type, extends(sparse_backend_t) :: umfpack_ipc_backend_t
        type(c_ptr) :: sess = c_null_ptr
        integer     :: n = 0
        logical     :: is_complex = .false.
    contains
        procedure :: factor_real => umf_factor_real
        procedure :: factor_complex => umf_factor_complex
        procedure :: solve_real => umf_solve_real
        procedure :: solve_complex => umf_solve_complex
        procedure :: free => umf_free
    end type umfpack_ipc_backend_t

    interface

        function fsparse_ipc_header_bytes() &
                bind(c, name="fsparse_ipc_header_bytes") result(b)
            import :: c_int64_t
            integer(c_int64_t) :: b
        end function fsparse_ipc_header_bytes

        function fsparse_ipc_acquire(helper, bytes, err) &
                bind(c, name="fsparse_ipc_acquire") result(sess)
            import :: c_ptr, c_char, c_int64_t, c_int
            character(kind=c_char), intent(in) :: helper(*)
            integer(c_int64_t),     value      :: bytes
            integer(c_int),         intent(out):: err
            type(c_ptr)                         :: sess
        end function fsparse_ipc_acquire

        function fsparse_ipc_data(sess) bind(c, name="fsparse_ipc_data") &
                result(p)
            import :: c_ptr
            type(c_ptr), value :: sess
            type(c_ptr)        :: p
        end function fsparse_ipc_data

        integer(c_int) function fsparse_ipc_call(sess) &
                bind(c, name="fsparse_ipc_call")
            import :: c_ptr, c_int
            type(c_ptr), value :: sess
        end function fsparse_ipc_call

        subroutine fsparse_ipc_release(sess) &
                bind(c, name="fsparse_ipc_release")
            import :: c_ptr
            type(c_ptr), value :: sess
        end subroutine fsparse_ipc_release

    end interface

contains

    ! Factor a real matrix in the helper. Lays out the shared region, writes
    ! the 0-based CSC arrays into it, then issues FACTOR_REAL.
    subroutine umf_factor_real(self, A, refine, status)
        class(umfpack_ipc_backend_t), intent(inout) :: self
        type(csc_t),                  intent(in)    :: A
        logical,                      intent(in)    :: refine
        type(fortsparse_status_t),    intent(out)   :: status

        type(shm_header_t), pointer :: h
        integer(i8) :: o_cp, o_ri, o_ax, o_b, o_x, total
        integer                     :: err

        call self%free()
        call layout_real(A%ncol, A%nnz, o_cp, o_ri, o_ax, o_b, o_x, total)
        call start_session(self, total, status)
        if (status%code /= FORTSPARSE_OK) return
        call header_of(self%sess, h)
        h%n = int(A%ncol, c_int64_t)
        h%nnz = int(A%nnz, c_int64_t)
        h%refine = merge(1_c_int32_t, 0_c_int32_t, refine)
        h%is_complex = 0_c_int32_t
        call set_offsets(h, o_cp, o_ri, o_ax, 0_i8, o_b, 0_i8, o_x, 0_i8, total)
        call write_index(self%sess, o_cp, A%col_ptr)
        call write_index(self%sess, o_ri, A%row_idx)
        call write_real(self%sess, o_ax, A%val)
        self%n = A%ncol
        self%is_complex = .false.
        call run(self, OP_FACTOR_REAL, status)
    end subroutine umf_factor_real

    ! Factor a complex matrix in the helper using split real/imag arrays.
    subroutine umf_factor_complex(self, A, refine, status)
        class(umfpack_ipc_backend_t), intent(inout) :: self
        type(csc_z_t),                intent(in)    :: A
        logical,                      intent(in)    :: refine
        type(fortsparse_status_t),    intent(out)   :: status

        type(shm_header_t), pointer :: h
        integer(i8) :: o_cp, o_ri, o_ax, o_az, o_b, o_bz, o_x, o_xz, total
        integer                     :: err

        call self%free()
        call layout_complex(A%ncol, A%nnz, o_cp, o_ri, o_ax, o_az, o_b, o_bz, &
            o_x, o_xz, total)
        call start_session(self, total, status)
        if (status%code /= FORTSPARSE_OK) return
        call header_of(self%sess, h)
        h%n = int(A%ncol, c_int64_t)
        h%nnz = int(A%nnz, c_int64_t)
        h%refine = merge(1_c_int32_t, 0_c_int32_t, refine)
        h%is_complex = 1_c_int32_t
        call set_offsets(h, o_cp, o_ri, o_ax, o_az, o_b, o_bz, o_x, o_xz, total)
        call write_index(self%sess, o_cp, A%col_ptr)
        call write_index(self%sess, o_ri, A%row_idx)
        call write_split(self%sess, o_ax, o_az, A%val)
        self%n = A%ncol
        self%is_complex = .true.
        call run(self, OP_FACTOR_COMPLEX, status)
    end subroutine umf_factor_complex

    ! Solve a real RHS through the resident factorization.
    subroutine umf_solve_real(self, b, x, status)
        class(umfpack_ipc_backend_t), intent(inout) :: self
        real(dp),                     intent(in)    :: b(:)
        real(dp),                     intent(out)   :: x(:)
        type(fortsparse_status_t),    intent(out)   :: status

        type(shm_header_t), pointer :: h

        call header_of(self%sess, h)
        call write_real(self%sess, int(h%off_b, i8), b)
        call run(self, OP_SOLVE_REAL, status)
        if (status%code /= FORTSPARSE_OK) return
        call read_real(self%sess, int(h%off_x, i8), x)
    end subroutine umf_solve_real

    ! Solve a complex RHS through the resident factorization.
    subroutine umf_solve_complex(self, b, x, status)
        class(umfpack_ipc_backend_t), intent(inout) :: self
        complex(dp),                  intent(in)    :: b(:)
        complex(dp),                  intent(out)   :: x(:)
        type(fortsparse_status_t),    intent(out)   :: status

        type(shm_header_t), pointer :: h

        call header_of(self%sess, h)
        call write_split(self%sess, int(h%off_b, i8), int(h%off_bz, i8), b)
        call run(self, OP_SOLVE_COMPLEX, status)
        if (status%code /= FORTSPARSE_OK) return
        call read_split(self%sess, int(h%off_x, i8), int(h%off_xz, i8), x)
    end subroutine umf_solve_complex

    ! Detach from the persistent helper session. The helper stays resident for
    ! the next factorization; fsparse_ipc_release does not tear it down.
    subroutine umf_free(self)
        class(umfpack_ipc_backend_t), intent(inout) :: self

        if (c_associated(self%sess)) call fsparse_ipc_release(self%sess)
        self%sess = c_null_ptr
        self%n = 0
        self%is_complex = .false.
    end subroutine umf_free

    ! Discover the helper, start the session, size the mapping. On a missing
    ! helper this reports FORTSPARSE_BACKEND_UNAVAILABLE, the expected state.
    subroutine start_session(self, total, status)
        class(umfpack_ipc_backend_t), intent(inout) :: self
        integer(i8),                  intent(in)    :: total
        type(fortsparse_status_t),    intent(out)   :: status

        character(:), allocatable :: helper
        integer                   :: err

        helper = find_helper()
        if (len(helper) == 0) then
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                "umfpack ipc: helper not found via "// &
                "FORTSPARSE_UMFPACK_HELPER or PATH")
            return
        end if
        self%sess = fsparse_ipc_acquire(helper//c_null_char, &
            int(total, c_int64_t), err)
        if (err /= 0 .or. .not. c_associated(self%sess)) then
            self%sess = c_null_ptr
            call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
                "umfpack ipc: failed to start helper "//helper)
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine start_session

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

    ! Real layout: header, colptr(n+1), rowidx(nnz), ax(nnz), b(n), x(n).
    subroutine layout_real(n, nnz, o_cp, o_ri, o_ax, o_b, o_x, total)
        integer,     intent(in)  :: n, nnz
        integer(i8), intent(out) :: o_cp, o_ri, o_ax, o_b, o_x, total

        o_cp = fsparse_ipc_header_bytes()
        o_ri = o_cp + int(n + 1, i8)*8_i8
        o_ax = o_ri + int(nnz, i8)*8_i8
        o_b = o_ax + int(nnz, i8)*8_i8
        o_x = o_b + int(n, i8)*8_i8
        total = o_x + int(n, i8)*8_i8
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

    ! Write a 1-based integer index array into the region as 0-based int64.
    subroutine write_index(sess, off, idx)
        type(c_ptr), intent(in) :: sess
        integer(i8), intent(in) :: off
        integer,     intent(in) :: idx(:)

        integer(c_int64_t), pointer :: dst(:)

        call c_f_pointer(region_at(sess, off), dst, [size(idx)])
        dst = int(idx, c_int64_t) - 1_c_int64_t
    end subroutine write_index

    ! Copy a real array into the region.
    subroutine write_real(sess, off, v)
        type(c_ptr), intent(in) :: sess
        integer(i8), intent(in) :: off
        real(dp),    intent(in) :: v(:)

        real(c_double), pointer :: dst(:)

        call c_f_pointer(region_at(sess, off), dst, [size(v)])
        dst = v
    end subroutine write_real

    ! Read a real array out of the region.
    subroutine read_real(sess, off, v)
        type(c_ptr), intent(in)  :: sess
        integer(i8), intent(in)  :: off
        real(dp),    intent(out) :: v(:)

        real(c_double), pointer :: src(:)

        call c_f_pointer(region_at(sess, off), src, [size(v)])
        v = src
    end subroutine read_real

    ! Split a complex array into separate real and imaginary buffers.
    subroutine write_split(sess, off_re, off_im, z)
        type(c_ptr), intent(in) :: sess
        integer(i8), intent(in) :: off_re, off_im
        complex(dp), intent(in) :: z(:)

        real(c_double), pointer :: re(:), im(:)

        call c_f_pointer(region_at(sess, off_re), re, [size(z)])
        call c_f_pointer(region_at(sess, off_im), im, [size(z)])
        re = real(z, dp)
        im = aimag(z)
    end subroutine write_split

    ! Reassemble a complex array from separate real and imaginary buffers.
    subroutine read_split(sess, off_re, off_im, z)
        type(c_ptr), intent(in)  :: sess
        integer(i8), intent(in)  :: off_re, off_im
        complex(dp), intent(out) :: z(:)

        real(c_double), pointer :: re(:), im(:)

        call c_f_pointer(region_at(sess, off_re), re, [size(z)])
        call c_f_pointer(region_at(sess, off_im), im, [size(z)])
        z = cmplx(re, im, dp)
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
