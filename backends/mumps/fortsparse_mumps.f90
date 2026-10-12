module fortsparse_mumps
    ! MUMPS (CeCILL-C) in-process backend: supernodal/multifrontal Cholesky,
    ! LDL^T and LU with METIS nested-dissection ordering, an OpenMP-threaded
    ! factorization and solve, one factorization reused for any number of
    ! right-hand sides, a native multi-RHS block solve and optional iterative
    ! refinement against the original matrix. Binds the MIT C shim
    ! (backends/mumps/fsparse_mumps.c) through iso_c_binding; real matrices only.
    use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, &
        c_int, c_long, c_double
    use fortsparse_kinds, only: dp, i8
    use fortsparse_backend, only: sparse_backend_t, sparse_factor_stats_t, &
        FORTSPARSE_SYM_GENERAL, FORTSPARSE_SYM_SPD, FORTSPARSE_SYM_INDEFINITE, &
        FORTSPARSE_ORDER_AMD
    use fortsparse_csc, only: csc_t, csc_z_t, csc_matvec
    use fortsparse_status, only: fortsparse_status_t, status_set, &
        FORTSPARSE_OK, FORTSPARSE_SINGULAR, FORTSPARSE_INTERNAL_ERROR, &
        FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_INVALID_MATRIX
    implicit none
    private

    public :: mumps_backend_t

    ! Refinement sweeps at most; stops earlier at the componentwise-free
    ! backward error below REFINE_TOL.
    integer,  parameter :: REFINE_STEPS = 4
    real(dp), parameter :: REFINE_TOL = 4.0_dp*epsilon(1.0_dp)

    ! MUMPS ICNTL(7) values.
    integer, parameter :: MUMPS_ORDER_AMD = 0
    integer, parameter :: MUMPS_ORDER_METIS = 5

    ! Owns one heap-allocated solve vector handed out by vector().
    type :: rvec_t
        real(dp), pointer :: p(:) => null()
    end type rvec_t

    type, extends(sparse_backend_t) :: mumps_backend_t
        type(rvec_t), allocatable :: vlist(:)
        integer :: nvec = 0
        type(c_ptr) :: h = c_null_ptr
        integer     :: n = 0
        logical     :: refine = .false.
        logical     :: symmetric = .false.
        real(dp)    :: anorm = 0.0_dp
        type(csc_t) :: a_real ! kept only when refinement is on
    contains
        procedure :: factor_real => mumps_factor_real
        procedure :: factor_complex => mumps_factor_complex
        procedure :: solve_real => mumps_solve_real
        procedure :: solve_complex => mumps_solve_complex
        procedure :: solve_real_inplace => mumps_solve_real_inplace
        procedure :: solve_complex_inplace => mumps_solve_complex_inplace
        procedure :: solve_real_block => mumps_solve_real_block
        procedure :: factor_stats => mumps_factor_stats
        procedure :: vector => mumps_vector
        procedure :: free => mumps_free
        final :: mumps_final
    end type mumps_backend_t

    interface

        function fsparse_mumps_factor_d(n, nnz, colptr1, rowidx1, val, sym, &
                ordering, nthreads, verbose, info1, info2) &
                bind(c, name="fsparse_mumps_factor_d") result(h)
            import :: c_ptr, c_int, c_long, c_double
            integer(c_long), value       :: n, nnz
            integer(c_long), intent(in)  :: colptr1(*), rowidx1(*)
            real(c_double),  intent(in)  :: val(*)
            integer(c_int),  value       :: sym, ordering, nthreads, verbose
            integer(c_int),  intent(out) :: info1, info2
            type(c_ptr)                  :: h
        end function fsparse_mumps_factor_d

        integer(c_int) function fsparse_mumps_solve_d(h, nrhs, rhs, ld, info2) &
                bind(c, name="fsparse_mumps_solve_d")
            import :: c_ptr, c_int, c_long, c_double
            type(c_ptr),     value        :: h
            integer(c_int),  value        :: nrhs
            real(c_double),  intent(inout) :: rhs(*)
            integer(c_long), value        :: ld
            integer(c_int),  intent(out)  :: info2
        end function fsparse_mumps_solve_d

        subroutine fsparse_mumps_stats_d(h, istat, rstat) &
                bind(c, name="fsparse_mumps_stats_d")
            import :: c_ptr, c_long, c_double
            type(c_ptr),     value       :: h
            integer(c_long), intent(out) :: istat(6)
            real(c_double),  intent(out) :: rstat(5)
        end subroutine fsparse_mumps_stats_d

        subroutine fsparse_mumps_free_d(h) bind(c, name="fsparse_mumps_free_d")
            import :: c_ptr
            type(c_ptr), value :: h
        end subroutine fsparse_mumps_free_d

    end interface

contains

    subroutine mumps_factor_real(self, A, refine, status)
        class(mumps_backend_t),    intent(inout) :: self
        type(csc_t),               intent(in)    :: A
        logical,                   intent(in)    :: refine
        type(fortsparse_status_t), intent(out)   :: status

        integer(c_long), allocatable :: colptr1(:), rowidx1(:)
        integer(c_int) :: info1, info2, sym, ordering
        character(80)  :: msg

        call self%free()
        if (A%nrow /= A%ncol) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "mumps factor: matrix is not square")
            return
        end if
        sym = int(self%options%symmetry, c_int)
        select case (self%options%symmetry)
        case (FORTSPARSE_SYM_GENERAL, FORTSPARSE_SYM_SPD, &
              FORTSPARSE_SYM_INDEFINITE)
        case default
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "mumps factor: unknown symmetry option")
            return
        end select
        if (self%options%ordering == FORTSPARSE_ORDER_AMD) then
            ordering = MUMPS_ORDER_AMD
        else
            ordering = MUMPS_ORDER_METIS
        end if
        colptr1 = int(A%col_ptr, c_long)
        rowidx1 = int(A%row_idx, c_long)
        self%h = fsparse_mumps_factor_d(int(A%ncol, c_long), &
            int(A%nnz, c_long), colptr1, rowidx1, A%val, sym, ordering, &
            int(self%options%nthreads, c_int), &
            merge(1_c_int, 0_c_int, self%options%verbose), info1, info2)
        if (.not. c_associated(self%h)) then
            write (msg, '(a,i0,a,i0)') "mumps factor failed: INFOG(1)=", &
                info1, " INFOG(2)=", info2
            select case (info1)
            case (-10, -40)
                call status_set(status, FORTSPARSE_SINGULAR, trim(msg))
            case default
                call status_set(status, FORTSPARSE_INTERNAL_ERROR, trim(msg))
            end select
            return
        end if
        self%n = A%ncol
        self%refine = refine
        self%symmetric = (sym /= 0)
        if (refine) then
            self%a_real = A
            self%anorm = inf_norm(A)
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine mumps_factor_real

    subroutine mumps_factor_complex(self, A, refine, status)
        class(mumps_backend_t),    intent(inout) :: self
        type(csc_z_t),             intent(in)    :: A
        logical,                   intent(in)    :: refine
        type(fortsparse_status_t), intent(out)   :: status

        call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
            "mumps backend supports real matrices only")
    end subroutine mumps_factor_complex

    ! Solve for nrhs columns in place; refinement acts on every column.
    subroutine solve_block_inplace(self, X, B, status)
        class(mumps_backend_t),    intent(inout) :: self
        real(dp), contiguous,      intent(inout) :: X(:, :)
        real(dp), contiguous,      intent(in)    :: B(:, :)
        type(fortsparse_status_t), intent(out)   :: status

        real(dp), allocatable :: R(:, :)
        real(dp) :: berr
        integer(c_int) :: info1, info2
        integer :: n, nrhs, j, step
        character(80) :: msg

        n = size(X, 1)
        nrhs = size(X, 2)
        if (n /= self%n) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "mumps solve: right-hand side length differs from the matrix")
            return
        end if
        if (nrhs == 0) then
            call status_set(status, FORTSPARSE_OK, "")
            return
        end if
        X = B
        info1 = fsparse_mumps_solve_d(self%h, int(nrhs, c_int), X, &
            int(n, c_long), info2)
        if (info1 /= 0) go to 900
        if (self%refine) then
            allocate (R(n, nrhs))
            do step = 1, REFINE_STEPS
                berr = 0.0_dp
                do j = 1, nrhs
                    call residual(self, B(:, j), X(:, j), R(:, j), berr)
                end do
                if (berr <= REFINE_TOL) exit
                info1 = fsparse_mumps_solve_d(self%h, int(nrhs, c_int), R, &
                    int(n, c_long), info2)
                if (info1 /= 0) go to 900
                X = X + R
            end do
        end if
        call status_set(status, FORTSPARSE_OK, "")
        return
900     write (msg, '(a,i0,a,i0)') "mumps solve failed: INFOG(1)=", info1, &
            " INFOG(2)=", info2
        call status_set(status, FORTSPARSE_INTERNAL_ERROR, trim(msg))
    end subroutine solve_block_inplace

    ! r = b - A x and the scaled backward error max(berr, |r|/(|A||x|+|b|)).
    subroutine residual(self, b, x, r, berr)
        class(mumps_backend_t), intent(in)    :: self
        real(dp),               intent(in)    :: b(:), x(:)
        real(dp),               intent(out)   :: r(:)
        real(dp),               intent(inout) :: berr

        real(dp) :: scale

        call apply(self, x, r)
        r = b - r
        scale = self%anorm*maxval(abs(x)) + maxval(abs(b))
        if (scale > 0.0_dp) berr = max(berr, maxval(abs(r))/scale)
    end subroutine residual

    ! y = A x; a symmetric matrix is applied as A^T x, whose per-column dot
    ! products are independent and thread-invariant.
    subroutine apply(self, x, y)
        class(mumps_backend_t), intent(in)  :: self
        real(dp),               intent(in)  :: x(:)
        real(dp),               intent(out) :: y(:)

        integer  :: j, e
        real(dp) :: s

        if (self%symmetric) then
            !$omp parallel do private(e, s) schedule(static)
            do j = 1, self%a_real%ncol
                s = 0.0_dp
                do e = self%a_real%col_ptr(j), self%a_real%col_ptr(j + 1) - 1
                    s = s + self%a_real%val(e)*x(self%a_real%row_idx(e))
                end do
                y(j) = s
            end do
            !$omp end parallel do
        else
            y = csc_matvec(self%a_real, x)
        end if
    end subroutine apply

    pure function inf_norm(A) result(nrm)
        type(csc_t), intent(in) :: A
        real(dp) :: nrm
        real(dp), allocatable :: rs(:)
        integer :: e

        allocate (rs(A%nrow), source=0.0_dp)
        do e = 1, A%nnz
            rs(A%row_idx(e)) = rs(A%row_idx(e)) + abs(A%val(e))
        end do
        nrm = maxval(rs)
    end function inf_norm

    subroutine mumps_solve_real(self, b, x, status)
        class(mumps_backend_t),       intent(inout) :: self
        real(dp), target, contiguous, intent(in)    :: b(:)
        real(dp), target, contiguous, intent(out)   :: x(:)
        type(fortsparse_status_t),    intent(out)   :: status

        real(dp), pointer :: b2(:, :), x2(:, :)

        b2(1:size(b), 1:1) => b
        x2(1:size(x), 1:1) => x
        call solve_block_inplace(self, x2, b2, status)
    end subroutine mumps_solve_real

    subroutine mumps_solve_real_inplace(self, b, status)
        class(mumps_backend_t),    intent(inout) :: self
        real(dp),                  intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out)   :: status

        real(dp), allocatable :: rhs(:), x(:)

        rhs = b
        allocate (x(size(b)))
        call mumps_solve_real(self, rhs, x, status)
        if (status%code == FORTSPARSE_OK) b = x
    end subroutine mumps_solve_real_inplace

    subroutine mumps_solve_real_block(self, B, X, status)
        class(mumps_backend_t),    intent(inout) :: self
        real(dp), contiguous,      intent(in)    :: B(:, :)
        real(dp), contiguous,      intent(out)   :: X(:, :)
        type(fortsparse_status_t), intent(out)   :: status

        call solve_block_inplace(self, X, B, status)
    end subroutine mumps_solve_real_block

    subroutine mumps_solve_complex(self, b, x, status)
        class(mumps_backend_t),    intent(inout) :: self
        complex(dp),               intent(in)    :: b(:)
        complex(dp),               intent(out)   :: x(:)
        type(fortsparse_status_t), intent(out)   :: status

        call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
            "mumps backend supports real matrices only")
    end subroutine mumps_solve_complex

    subroutine mumps_solve_complex_inplace(self, b, status)
        class(mumps_backend_t),    intent(inout) :: self
        complex(dp),               intent(inout) :: b(:)
        type(fortsparse_status_t), intent(out)   :: status

        call status_set(status, FORTSPARSE_BACKEND_UNAVAILABLE, &
            "mumps backend supports real matrices only")
    end subroutine mumps_solve_complex_inplace

    subroutine mumps_factor_stats(self, stats, status)
        class(mumps_backend_t),      intent(inout) :: self
        type(sparse_factor_stats_t), intent(out)   :: stats
        type(fortsparse_status_t),   intent(out)   :: status

        integer(c_long) :: istat(6)
        real(c_double)  :: rstat(5)

        if (.not. c_associated(self%h)) then
            call status_set(status, FORTSPARSE_INTERNAL_ERROR, &
                "mumps stats: no factorization")
            return
        end if
        call fsparse_mumps_stats_d(self%h, istat, rstat)
        stats%factor_nnz = int(istat(1), i8)
        stats%factor_mem_mb = real(istat(2), dp)
        stats%flops = rstat(3)
        stats%analysis_s = rstat(4)
        stats%factor_s = rstat(5)
        stats%nthreads = int(istat(5))
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine mumps_factor_stats

    ! A plain heap array owned by the backend and freed at teardown, uniform
    ! with the other backends.
    function mumps_vector(self, n) result(p)
        class(mumps_backend_t), intent(inout) :: self
        integer,                intent(in)    :: n
        real(dp), pointer                     :: p(:)

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
    end function mumps_vector

    subroutine mumps_free(self)
        class(mumps_backend_t), intent(inout) :: self

        if (c_associated(self%h)) call fsparse_mumps_free_d(self%h)
        self%h = c_null_ptr
        self%n = 0
        self%refine = .false.
        self%symmetric = .false.
        if (allocated(self%a_real%col_ptr)) deallocate (self%a_real%col_ptr)
        if (allocated(self%a_real%row_idx)) deallocate (self%a_real%row_idx)
        if (allocated(self%a_real%val)) deallocate (self%a_real%val)
    end subroutine mumps_free

    subroutine mumps_final(self)
        type(mumps_backend_t), intent(inout) :: self

        integer :: i

        if (c_associated(self%h)) call fsparse_mumps_free_d(self%h)
        self%h = c_null_ptr
        if (allocated(self%vlist)) then
            do i = 1, self%nvec
                if (associated(self%vlist(i)%p)) deallocate (self%vlist(i)%p)
            end do
        end if
    end subroutine mumps_final

end module fortsparse_mumps
