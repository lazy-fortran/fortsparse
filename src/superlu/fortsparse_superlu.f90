module fortsparse_superlu
    ! SuperLU (BSD) in-process backend. Binds the MIT C shim
    ! (src/superlu/fsparse_superlu.c) through iso_c_binding. CSC index arrays
    ! are passed 1-based; the shim converts to SuperLU's
    ! 0-based int_t internally. The opaque factorization lives behind a c_ptr
    ! held in the backend object and is reused for every solve.
    use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, &
                                            c_int, c_int64_t, c_double, c_double_complex
    use fortsparse_kinds, only: dp, i8
    use fortsparse_backend, only: sparse_backend_t, FORTSPARSE_SYM_GENERAL, &
                                  FORTSPARSE_OP_NORMAL
    use fortsparse_csc, only: csc_t, csc_z_t, csc_matvec
    use fortsparse_status, only: fortsparse_status_t, status_set, &
        FORTSPARSE_OK, FORTSPARSE_SINGULAR, FORTSPARSE_INTERNAL_ERROR, &
        FORTSPARSE_OUT_OF_MEMORY, FORTSPARSE_INDEX_OVERFLOW
    implicit none
    private

    public :: superlu_backend_t
    public :: slu_map_factor_status

    ! SuperLU backend handle. Exactly one of h_real / h_cplx is associated after
    ! a successful factorization; is_complex records which.
    ! One iterative-refinement sweep when refinement is enabled. SuperLU's
    ! simple driver does no refinement, so the backend retains the matrix and
    ! sharpens the solution from the residual r = b - A x.
    integer, parameter :: REFINE_STEPS = 2

    ! Owns one heap-allocated solve vector handed out by the backend.
    type :: rvec_t
        real(dp), pointer :: p(:) => null()
    end type rvec_t

    type, extends(sparse_backend_t) :: superlu_backend_t
        type(c_ptr) :: h_real = c_null_ptr
        type(c_ptr) :: h_cplx = c_null_ptr
        logical     :: is_complex = .false.
        logical     :: refine = .false.
        type(csc_t)   :: a_real
        type(csc_z_t) :: a_cplx
        real(dp),    allocatable :: xbuf(:)
        complex(dp), allocatable :: zbuf(:)
        real(dp), allocatable :: residual_real(:), correction_real(:)
        complex(dp), allocatable :: residual_complex(:), correction_complex(:)
        type(rvec_t), allocatable :: vlist(:) ! vectors handed out by vector()
        integer                   :: nvec = 0
    contains
        procedure :: factor_real => slu_factor_real
        procedure :: factor_real_raw => slu_factor_real_raw
        procedure :: refactor_real => slu_refactor_real
        procedure :: refactor_real_raw => slu_refactor_real_raw
        procedure :: factor_complex => slu_factor_complex
        procedure :: solve_real => slu_solve_real
        procedure :: solve_complex => slu_solve_complex
        procedure :: solve_real_inplace => slu_solve_real_inplace
        procedure :: solve_complex_inplace => slu_solve_complex_inplace
        procedure :: solve_real_op => slu_solve_real_op
        procedure :: solve_complex_op => slu_solve_complex_op
        procedure :: vector => slu_vector
        procedure :: free => slu_free
        final :: slu_final
    end type superlu_backend_t

    interface

        function fsparse_slu_factor_d(n, nnz, colptr1, rowidx1, val, symmetric, &
                info) &
                bind(c, name="fsparse_slu_factor_d") result(h)
            import :: c_ptr, c_int, c_double
            integer(c_int), value          :: n, nnz
            integer(c_int), intent(in)     :: colptr1(*), rowidx1(*)
            real(c_double),     intent(in)  :: val(*)
            integer(c_int),     value       :: symmetric
            integer(c_int),     intent(out) :: info
            type(c_ptr)                     :: h
        end function fsparse_slu_factor_d

        integer(c_int) function fsparse_slu_refactor_d(h, n, nnz, colptr1, &
                rowidx1, val, symmetric) bind(c, name="fsparse_slu_refactor_d")
            import :: c_ptr, c_int, c_double
            type(c_ptr),    value      :: h
            integer(c_int), value      :: n, nnz
            integer(c_int), intent(in) :: colptr1(*), rowidx1(*)
            real(c_double), intent(in) :: val(*)
            integer(c_int), value      :: symmetric
        end function fsparse_slu_refactor_d

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

        function fsparse_slu_factor_z(n, nnz, colptr1, rowidx1, val, symmetric, &
                info) &
                bind(c, name="fsparse_slu_factor_z") result(h)
            import :: c_ptr, c_int, c_int64_t, c_double
            integer(c_int64_t), value       :: n, nnz
            integer(c_int64_t), intent(in)  :: colptr1(*), rowidx1(*)
            real(c_double),     intent(in)  :: val(*)
            integer(c_int),     value       :: symmetric
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

        integer(c_int) function fsparse_slu_solve_op_d(h, b, x, n, op) &
            bind(c, name="fsparse_slu_solve_op_d")
            import :: c_ptr, c_int, c_int64_t, c_double, c_double_complex
            type(c_ptr), value :: h
            real(c_double), intent(in) :: b(*)
            real(c_double), intent(out) :: x(*)
            integer(c_int64_t), value :: n
            integer(c_int), value :: op
        end function fsparse_slu_solve_op_d

        integer(c_int) function fsparse_slu_solve_op_z(h, b, x, n, op) &
            bind(c, name="fsparse_slu_solve_op_z")
            import :: c_ptr, c_int, c_int64_t, c_double, c_double_complex
            type(c_ptr), value :: h
            complex(c_double_complex), intent(in) :: b(*)
            complex(c_double_complex), intent(out) :: x(*)
            integer(c_int64_t), value :: n
            integer(c_int), value :: op
        end function fsparse_slu_solve_op_z

    end interface

contains

    ! Factor a real matrix. Shim info: 0 ok, 1 singular,
    ! 2 allocation failure, 3 index overflow, <0 internal error.
    subroutine slu_factor_real(self, A, refine, status)
        class(superlu_backend_t),  intent(inout) :: self
        type(csc_t),               intent(in)    :: A
        logical,                   intent(in)    :: refine
        type(fortsparse_status_t), intent(out)   :: status

        call self%factor_real_raw(A%nrow, A%ncol, A%nnz, A%col_ptr, &
            A%row_idx, A%val, refine, status)
    end subroutine slu_factor_real

    ! Convert caller-owned CSC directly into the C factor storage. Retain a
    ! Fortran copy only when it is needed for refinement residuals.
    subroutine slu_factor_real_raw(self, nrow, ncol, nz, col_ptr, &
            row_idx, val, refine, status)
        class(superlu_backend_t),  intent(inout) :: self
        integer,                   intent(in)    :: nrow, ncol, nz
        integer,                   intent(in)    :: col_ptr(:), row_idx(:)
        real(dp),                  intent(in)    :: val(:)
        logical,                   intent(in)    :: refine
        type(fortsparse_status_t), intent(out)   :: status

        integer(c_int) :: info

        call self%free()
        self%h_real = fsparse_slu_factor_d(ncol, nz, col_ptr, row_idx, val, &
            symmetric_hint(self), info)
        call slu_map_factor_status(info, status)
        if (info /= 0) return
        self%refine = refine
        if (refine) then
            self%a_real%nrow = nrow
            self%a_real%ncol = ncol
            self%a_real%nnz = nz
            self%a_real%col_ptr = col_ptr(1:ncol + 1)
            self%a_real%row_idx = row_idx(1:nz)
            self%a_real%val = val(1:nz)
        end if
    end subroutine slu_factor_real_raw

    ! Same pattern and symmetry hint as the retained real factorization: keep
    ! the column ordering (pattern-only, the minimum-degree cost) and redo the
    ! numerical factorization; the factor equals a fresh one on the same input.
    ! Otherwise a full factor.
    subroutine slu_refactor_real(self, A, refine, status, reused)
        class(superlu_backend_t),  intent(inout) :: self
        type(csc_t),               intent(in)    :: A
        logical,                   intent(in)    :: refine
        type(fortsparse_status_t), intent(out)   :: status
        logical,                   intent(out)   :: reused

        call self%refactor_real_raw(A%nrow, A%ncol, A%nnz, A%col_ptr, &
            A%row_idx, A%val, refine, status, reused)
    end subroutine slu_refactor_real

    subroutine slu_refactor_real_raw(self, nrow, ncol, nz, col_ptr, row_idx, &
            val, refine, status, reused)
        class(superlu_backend_t),  intent(inout) :: self
        integer,                   intent(in)    :: nrow, ncol, nz
        integer,                   intent(in)    :: col_ptr(:), row_idx(:)
        real(dp),                  intent(in)    :: val(:)
        logical,                   intent(in)    :: refine
        type(fortsparse_status_t), intent(out)   :: status
        logical,                   intent(out)   :: reused

        integer(c_int) :: info

        reused = .false.
        info = 4
        if (c_associated(self%h_real)) info = fsparse_slu_refactor_d( &
            self%h_real, ncol, nz, col_ptr, row_idx, val, symmetric_hint(self))
        if (info == 4) then
            call self%factor_real_raw(nrow, ncol, nz, col_ptr, row_idx, val, &
                refine, status)
            return
        end if
        call slu_map_factor_status(info, status)
        if (info /= 0) then
            call self%free()
            return
        end if
        reused = .true.
        self%refine = refine
        if (refine) then
            self%a_real%nrow = nrow
            self%a_real%ncol = ncol
            self%a_real%nnz = nz
            self%a_real%col_ptr = col_ptr(1:ncol + 1)
            self%a_real%row_idx = row_idx(1:nz)
            self%a_real%val = val(1:nz)
        else
            self%a_real = csc_t()
        end if
    end subroutine slu_refactor_real_raw

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
            colptr1, rowidx1, val, symmetric_hint(self), info)
        call slu_map_factor_status(info, status)
        if (info == 0) then
            self%is_complex = .true.
            self%refine = refine
            if (refine) self%a_cplx = A
        end if
    end subroutine slu_factor_complex

    ! Solve A x = b for a real RHS using the retained factorization.
    subroutine slu_solve_real(self, b, x, status)
        class(superlu_backend_t),     intent(inout) :: self
        real(dp), target, contiguous, intent(in)    :: b(:)
        real(dp), target, contiguous, intent(out)   :: x(:)
        type(fortsparse_status_t),    intent(out)   :: status
        call slu_solve_real_op(self, b, x, FORTSPARSE_OP_NORMAL, status)
    end subroutine slu_solve_real

    subroutine slu_solve_real_op(self, b, x, operation, status)
        class(superlu_backend_t), intent(inout) :: self
        real(dp), target, contiguous, intent(in) :: b(:)
        real(dp), target, contiguous, intent(out) :: x(:)
        integer, intent(in) :: operation
        type(fortsparse_status_t), intent(out) :: status
        integer(c_int)        :: info
        integer :: step, j, k
        real(dp) :: coefficient
        info = fsparse_slu_solve_op_d(self%h_real, b, x, &
                                      int(size(b), i8), int(operation, c_int))
        call map_solve_status(info, status)
        if (info /= 0 .or. .not. self%refine) return
        if (.not. allocated(self%residual_real)) then
            allocate (self%residual_real(size(b)), self%correction_real(size(b)))
        end if
        do step = 1, REFINE_STEPS
            self%residual_real = b
            do j = 1, self%a_real%ncol
                do k = self%a_real%col_ptr(j), &
                    self%a_real%col_ptr(j + 1) - 1
                    coefficient = self%a_real%val(k)
                    if (operation == 0) then
                        self%residual_real(self%a_real%row_idx(k)) = &
                            self%residual_real(self%a_real%row_idx(k)) - &
                            coefficient*x(j)
                    else
                        self%residual_real(j) = self%residual_real(j) - &
                                                coefficient*x(self%a_real%row_idx(k))
                    end if
                end do
            end do
            info = fsparse_slu_solve_op_d(self%h_real, &
                           self%residual_real, self%correction_real, int(size(b), i8), &
                                          int(operation, c_int))
            if (info /= 0) exit
            x = x + self%correction_real
        end do
        call map_solve_status(info, status)
    end subroutine slu_solve_real_op

    subroutine slu_solve_complex(self, b, x, status)
        class(superlu_backend_t),  intent(inout) :: self
        complex(dp),               intent(in)    :: b(:)
        complex(dp),               intent(out)   :: x(:)
        type(fortsparse_status_t), intent(out)   :: status
        call slu_solve_complex_op(self, b, x, FORTSPARSE_OP_NORMAL, status)
    end subroutine slu_solve_complex

    subroutine slu_solve_complex_op(self, b, x, operation, status)
        class(superlu_backend_t), intent(inout) :: self
        complex(dp), intent(in) :: b(:)
        complex(dp), intent(out) :: x(:)
        integer, intent(in) :: operation
        type(fortsparse_status_t), intent(out) :: status
        integer(c_int)           :: info
        integer :: step, j, k
        complex(dp) :: coefficient
        info = fsparse_slu_solve_op_z(self%h_cplx, b, x, &
                                      int(size(b), i8), int(operation, c_int))
        call map_solve_status(info, status)
        if (info /= 0 .or. .not. self%refine) return
        if (.not. allocated(self%residual_complex)) then
            allocate (self%residual_complex(size(b)), self%correction_complex(size(b)))
        end if
        do step = 1, REFINE_STEPS
            self%residual_complex = b
            do j = 1, self%a_cplx%ncol
                do k = self%a_cplx%col_ptr(j), &
                    self%a_cplx%col_ptr(j + 1) - 1
                    coefficient = self%a_cplx%val(k)
                    if (operation == 2) coefficient = conjg(coefficient)
                    if (operation == 0) then
                        self%residual_complex(self%a_cplx%row_idx(k)) = &
                            self%residual_complex(self%a_cplx%row_idx(k)) - &
                            coefficient*x(j)
                    else
                        self%residual_complex(j) = self%residual_complex(j) - &
                                                   coefficient*x(self%a_cplx%row_idx(k))
                    end if
                end do
            end do
            info = fsparse_slu_solve_op_z(self%h_cplx, &
                     self%residual_complex, self%correction_complex, int(size(b), i8), &
                                          int(operation, c_int))
            if (info /= 0) exit
            x = x + self%correction_complex
        end do
        call map_solve_status(info, status)
    end subroutine slu_solve_complex_op

    ! In-place real solve: b is the RHS on entry, the solution on return. Uses a
    ! backend-held scratch (grown once, not per solve) so the in-process solve
    ! needs no caller temporary.
    subroutine slu_solve_real_inplace(self, b, status)
        class(superlu_backend_t),  intent(inout) :: self
        real(dp),                  intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out)   :: status

        if (.not. allocated(self%xbuf)) then
            allocate (self%xbuf(size(b)))
        else if (size(self%xbuf) < size(b)) then
            deallocate (self%xbuf)
            allocate (self%xbuf(size(b)))
        end if
        call slu_solve_real(self, b, self%xbuf(1:size(b)), status)
        if (status%code == FORTSPARSE_OK) b = self%xbuf(1:size(b))
    end subroutine slu_solve_real_inplace

    ! In-place complex solve; b is the RHS on entry, the solution on return.
    subroutine slu_solve_complex_inplace(self, b, status)
        class(superlu_backend_t),  intent(inout) :: self
        complex(dp),               intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out)   :: status

        if (.not. allocated(self%zbuf)) then
            allocate (self%zbuf(size(b)))
        else if (size(self%zbuf) < size(b)) then
            deallocate (self%zbuf)
            allocate (self%zbuf(size(b)))
        end if
        call slu_solve_complex(self, b, self%zbuf(1:size(b)), status)
        if (status%code == FORTSPARSE_OK) b = self%zbuf(1:size(b))
    end subroutine slu_solve_complex_inplace

    ! In-process backend: a solve vector is a plain heap array the backend owns
    ! and frees on teardown (the in-process solve already copies nothing). Kept
    ! uniform with the out-of-process backend so client code is identical.
    function slu_vector(self, n) result(p)
        class(superlu_backend_t), intent(inout) :: self
        integer,                  intent(in)    :: n
        real(dp), pointer                       :: p(:)

        type(rvec_t), allocatable :: tmp(:)

        allocate (p(n))
        if (.not. allocated(self%vlist)) allocate (self%vlist(8))
        if (self%nvec >= size(self%vlist)) then
            allocate (tmp(2*size(self%vlist)))
            tmp(1:self%nvec) = self%vlist(1:self%nvec)
            call move_alloc(tmp, self%vlist)
        end if
        self%nvec = self%nvec + 1
        self%vlist(self%nvec)%p => p
    end function slu_vector

    ! Release the retained factorization, if any.
    subroutine slu_free(self)
        class(superlu_backend_t), intent(inout) :: self

        if (c_associated(self%h_real)) call fsparse_slu_free_d(self%h_real)
        if (c_associated(self%h_cplx)) call fsparse_slu_free_z(self%h_cplx)
        self%h_real = c_null_ptr
        self%h_cplx = c_null_ptr
        self%is_complex = .false.
        self%refine = .false.
        if (allocated(self%residual_real)) then
            deallocate (self%residual_real, self%correction_real)
        end if
        if (allocated(self%residual_complex)) then
            deallocate (self%residual_complex, self%correction_complex)
        end if
        self%a_real = csc_t()
        self%a_cplx = csc_z_t()
    end subroutine slu_free

    ! Free the retained factorization when the backend is destroyed. Runs
    ! automatically when a solver (and its allocatable backend) goes out of
    ! scope, so client code needs no explicit teardown.
    subroutine slu_final(self)
        type(superlu_backend_t), intent(inout) :: self

        integer :: i

        if (c_associated(self%h_real)) call fsparse_slu_free_d(self%h_real)
        if (c_associated(self%h_cplx)) call fsparse_slu_free_z(self%h_cplx)
        self%h_real = c_null_ptr
        self%h_cplx = c_null_ptr
        if (allocated(self%vlist)) then
            do i = 1, self%nvec
                if (associated(self%vlist(i)%p)) deallocate (self%vlist(i)%p)
            end do
        end if
    end subroutine slu_final

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

    ! 1 when the caller declared a symmetric matrix (SPD, INDEFINITE or
    ! HERMITIAN): SuperLU then uses its symmetric-mode ordering and pivoting
    ! on the full matrix as given; 0 for GENERAL.
    integer(c_int) function symmetric_hint(self)
        class(superlu_backend_t), intent(in) :: self
        symmetric_hint = merge(1_c_int, 0_c_int, &
            self%options%symmetry /= FORTSPARSE_SYM_GENERAL)
    end function symmetric_hint

    ! Map the shim factor info code onto a fortsparse status.
    subroutine slu_map_factor_status(info, status)
        integer(c_int),            intent(in)  :: info
        type(fortsparse_status_t), intent(out) :: status

        if (info == 0) then
            call status_set(status, FORTSPARSE_OK, "")
        else if (info == 1) then
            call status_set(status, FORTSPARSE_SINGULAR, &
                "superlu factor: matrix is singular (exact zero pivot)")
        else if (info == 2) then
            call status_set(status, FORTSPARSE_OUT_OF_MEMORY, &
                "superlu factor: memory allocation failed (LU fill too large)")
        else if (info == 3) then
            call status_set(status, FORTSPARSE_INDEX_OVERFLOW, &
                "superlu factor: matrix size exceeds SuperLU index range")
        else if (info > 3) then
            call status_set(status, FORTSPARSE_INTERNAL_ERROR, &
                "superlu factor: unknown positive info code")
        else
            call status_set(status, FORTSPARSE_INTERNAL_ERROR, &
                "superlu factor: internal error")
        end if
    end subroutine slu_map_factor_status

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
