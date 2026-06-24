module fortsparse_superlu
    ! SuperLU (BSD) in-process backend. Binds the MIT C shim
    ! (src/superlu/fsparse_superlu.c) through iso_c_binding. CSC index arrays
    ! are passed 1-based as 64-bit integers; the shim converts to SuperLU's
    ! 0-based int_t internally. The opaque factorization lives behind a c_ptr
    ! held in the backend object and is reused for every solve.
    use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, &
        c_int, c_int64_t, c_double
    use fortsparse_kinds, only: dp, i8
    use fortsparse_backend, only: sparse_backend_t
    use fortsparse_csc, only: csc_t, csc_z_t, csc_matvec
    use fortsparse_status, only: fortsparse_status_t, status_set, &
        FORTSPARSE_OK, FORTSPARSE_SINGULAR, FORTSPARSE_INTERNAL_ERROR
    implicit none
    private

    public :: superlu_backend_t

    ! SuperLU backend handle. Exactly one of h_real / h_cplx is associated after
    ! a successful factorization; is_complex records which.
    ! One iterative-refinement sweep when refinement is enabled. SuperLU's
    ! simple driver does no refinement, so the backend retains the matrix and
    ! sharpens the solution from the residual r = b - A x.
    integer, parameter :: REFINE_STEPS = 2

    type, extends(sparse_backend_t) :: superlu_backend_t
        type(c_ptr) :: h_real = c_null_ptr
        type(c_ptr) :: h_cplx = c_null_ptr
        logical     :: is_complex = .false.
        logical     :: refine = .false.
        type(csc_t)   :: a_real
        type(csc_z_t) :: a_cplx
    contains
        procedure :: factor_real => slu_factor_real
        procedure :: factor_complex => slu_factor_complex
        procedure :: solve_real => slu_solve_real
        procedure :: solve_complex => slu_solve_complex
        procedure :: free => slu_free
    end type superlu_backend_t

    interface

        function fsparse_slu_factor_d(n, nnz, colptr1, rowidx1, val, info) &
                bind(c, name="fsparse_slu_factor_d") result(h)
            import :: c_ptr, c_int, c_int64_t, c_double
            integer(c_int64_t), value       :: n, nnz
            integer(c_int64_t), intent(in)  :: colptr1(*), rowidx1(*)
            real(c_double),     intent(in)  :: val(*)
            integer(c_int),     intent(out) :: info
            type(c_ptr)                     :: h
        end function fsparse_slu_factor_d

        integer(c_int) function fsparse_slu_solve_d(h, b, x, n) &
                bind(c, name="fsparse_slu_solve_d")
            import :: c_ptr, c_int, c_int64_t, c_double
            type(c_ptr),        value       :: h
            real(c_double),     intent(in)  :: b(*)
            real(c_double),     intent(out) :: x(*)
            integer(c_int64_t), value       :: n
        end function fsparse_slu_solve_d

        subroutine fsparse_slu_free_d(h) bind(c, name="fsparse_slu_free_d")
            import :: c_ptr
            type(c_ptr), value :: h
        end subroutine fsparse_slu_free_d

        function fsparse_slu_factor_z(n, nnz, colptr1, rowidx1, val, info) &
                bind(c, name="fsparse_slu_factor_z") result(h)
            import :: c_ptr, c_int, c_int64_t, c_double
            integer(c_int64_t), value       :: n, nnz
            integer(c_int64_t), intent(in)  :: colptr1(*), rowidx1(*)
            real(c_double),     intent(in)  :: val(*)
            integer(c_int),     intent(out) :: info
            type(c_ptr)                     :: h
        end function fsparse_slu_factor_z

        integer(c_int) function fsparse_slu_solve_z(h, b, x, n) &
                bind(c, name="fsparse_slu_solve_z")
            import :: c_ptr, c_int, c_int64_t, c_double
            type(c_ptr),        value       :: h
            real(c_double),     intent(in)  :: b(*)
            real(c_double),     intent(out) :: x(*)
            integer(c_int64_t), value       :: n
        end function fsparse_slu_solve_z

        subroutine fsparse_slu_free_z(h) bind(c, name="fsparse_slu_free_z")
            import :: c_ptr
            type(c_ptr), value :: h
        end subroutine fsparse_slu_free_z

    end interface

contains

    ! Factor a real matrix. SuperLU info: 0 ok, 1 singular, <0 internal error.
    subroutine slu_factor_real(self, A, refine, status)
        class(superlu_backend_t),  intent(inout) :: self
        type(csc_t),               intent(in)    :: A
        logical,                   intent(in)    :: refine
        type(fortsparse_status_t), intent(out)   :: status

        integer(i8), allocatable :: colptr1(:), rowidx1(:)
        integer(c_int)           :: info

        call self%free()
        colptr1 = int(A%col_ptr, i8)
        rowidx1 = int(A%row_idx, i8)
        self%h_real = fsparse_slu_factor_d(int(A%ncol, i8), int(A%nnz, i8), &
            colptr1, rowidx1, A%val, info)
        call map_factor_status(info, status)
        if (info == 0) then
            self%is_complex = .false.
            self%refine = refine
            self%a_real = A
        end if
    end subroutine slu_factor_real

    ! Factor a complex matrix. complex(dp) maps onto SuperLU's interleaved
    ! doublecomplex directly, so values pass through as a real(dp) buffer.
    subroutine slu_factor_complex(self, A, refine, status)
        class(superlu_backend_t),  intent(inout) :: self
        type(csc_z_t),             intent(in)    :: A
        logical,                   intent(in)    :: refine
        type(fortsparse_status_t), intent(out)   :: status

        integer(i8), allocatable :: colptr1(:), rowidx1(:)
        real(dp),    allocatable :: val(:)
        integer(c_int)           :: info

        call self%free()
        colptr1 = int(A%col_ptr, i8)
        rowidx1 = int(A%row_idx, i8)
        val = interleave_complex(A%val)
        self%h_cplx = fsparse_slu_factor_z(int(A%ncol, i8), int(A%nnz, i8), &
            colptr1, rowidx1, val, info)
        call map_factor_status(info, status)
        if (info == 0) then
            self%is_complex = .true.
            self%refine = refine
            self%a_cplx = A
        end if
    end subroutine slu_factor_complex

    ! Solve A x = b for a real RHS using the retained factorization.
    subroutine slu_solve_real(self, b, x, status)
        class(superlu_backend_t),  intent(inout) :: self
        real(dp),                  intent(in)    :: b(:)
        real(dp),                  intent(out)   :: x(:)
        type(fortsparse_status_t), intent(out)   :: status

        real(dp), allocatable :: r(:), dx(:)
        integer(c_int)        :: info
        integer               :: step

        info = fsparse_slu_solve_d(self%h_real, b, x, int(size(b), i8))
        call map_solve_status(info, status)
        if (info /= 0 .or. .not. self%refine) return
        allocate (dx(size(b)))
        do step = 1, REFINE_STEPS
            r = b - csc_matvec(self%a_real, x)
            info = fsparse_slu_solve_d(self%h_real, r, dx, int(size(b), i8))
            if (info /= 0) exit
            x = x + dx
        end do
        call map_solve_status(info, status)
    end subroutine slu_solve_real

    ! Solve A x = b for a complex RHS using the retained factorization.
    subroutine slu_solve_complex(self, b, x, status)
        class(superlu_backend_t),  intent(inout) :: self
        complex(dp),               intent(in)    :: b(:)
        complex(dp),               intent(out)   :: x(:)
        type(fortsparse_status_t), intent(out)   :: status

        real(dp), allocatable    :: bi(:), xi(:), ri(:), dxi(:)
        complex(dp), allocatable :: r(:), dx(:)
        integer(c_int)           :: info
        integer                  :: n, step

        n = size(b)
        bi = interleave_complex(b)
        allocate (xi(2*n))
        info = fsparse_slu_solve_z(self%h_cplx, bi, xi, int(n, i8))
        call map_solve_status(info, status)
        if (info /= 0) return
        x = deinterleave_complex(xi)
        if (.not. self%refine) return
        allocate (dxi(2*n))
        do step = 1, REFINE_STEPS
            r = b - csc_matvec(self%a_cplx, x)
            ri = interleave_complex(r)
            info = fsparse_slu_solve_z(self%h_cplx, ri, dxi, int(n, i8))
            if (info /= 0) exit
            dx = deinterleave_complex(dxi)
            x = x + dx
        end do
        call map_solve_status(info, status)
    end subroutine slu_solve_complex

    ! Release the retained factorization, if any.
    subroutine slu_free(self)
        class(superlu_backend_t), intent(inout) :: self

        if (c_associated(self%h_real)) call fsparse_slu_free_d(self%h_real)
        if (c_associated(self%h_cplx)) call fsparse_slu_free_z(self%h_cplx)
        self%h_real = c_null_ptr
        self%h_cplx = c_null_ptr
        self%is_complex = .false.
        self%refine = .false.
    end subroutine slu_free

    ! Pack a complex array into interleaved (real, imag) real(dp) pairs.
    pure function interleave_complex(z) result(r)
        complex(dp), intent(in) :: z(:)
        real(dp), allocatable   :: r(:)
        integer                 :: k

        allocate (r(2*size(z)))
        do k = 1, size(z)
            r(2*k - 1) = real(z(k), dp)
            r(2*k) = aimag(z(k))
        end do
    end function interleave_complex

    ! Inverse of interleave_complex.
    pure function deinterleave_complex(r) result(z)
        real(dp), intent(in)     :: r(:)
        complex(dp), allocatable :: z(:)
        integer                  :: k, n

        n = size(r)/2
        allocate (z(n))
        do k = 1, n
            z(k) = cmplx(r(2*k - 1), r(2*k), dp)
        end do
    end function deinterleave_complex

    ! Map the shim factor info code onto a fortsparse status.
    subroutine map_factor_status(info, status)
        integer(c_int),            intent(in)  :: info
        type(fortsparse_status_t), intent(out) :: status

        if (info == 0) then
            call status_set(status, FORTSPARSE_OK, "")
        else if (info > 0) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "superlu factor: matrix is singular")
        else
            call status_set(status, FORTSPARSE_INTERNAL_ERROR, &
                "superlu factor: internal error")
        end if
    end subroutine map_factor_status

    ! Map the shim solve info code onto a fortsparse status.
    subroutine map_solve_status(info, status)
        integer(c_int),            intent(in)  :: info
        type(fortsparse_status_t), intent(out) :: status

        if (info == 0) then
            call status_set(status, FORTSPARSE_OK, "")
        else
            call status_set(status, FORTSPARSE_INTERNAL_ERROR, &
                "superlu solve: internal error")
        end if
    end subroutine map_solve_status

end module fortsparse_superlu
