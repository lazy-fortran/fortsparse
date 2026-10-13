program test_fortsparse_mumps_complex
    ! Complex (ZMUMPS) path of the MUMPS backend against dense LAPACK zgesv:
    ! general non-Hermitian, Hermitian indefinite and positive definite,
    ! complex symmetric, factor-once/solve-many, block solve, in-place solve,
    ! concurrent solver objects under OpenMP, and singular detection.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: dp, csc_z_t, csc_from_triplet, sparse_solver_t, &
        sparse_factor, sparse_refactor, sparse_solve, sparse_solve_block, sparse_destroy, &
        fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_MUMPS, &
        FORTSPARSE_SYM_GENERAL, FORTSPARSE_SYM_SPD, FORTSPARSE_SYM_INDEFINITE, &
        FORTSPARSE_SYM_HERMITIAN, FORTSPARSE_SINGULAR
    use fs_test_mpi, only: test_mpi_init, test_mpi_finalize
    implicit none

    integer :: nfail
    integer :: seed_state = 12345

    call test_mpi_init()
    nfail = 0
    call test_kind("general", FORTSPARSE_SYM_GENERAL, 0, 300, nfail)
    call test_kind("hermitian-indefinite", FORTSPARSE_SYM_HERMITIAN, 1, 300, nfail)
    call test_kind("hermitian-posdef", FORTSPARSE_SYM_SPD, 2, 300, nfail)
    call test_kind("complex-symmetric", FORTSPARSE_SYM_INDEFINITE, 3, 300, nfail)
    call test_refactor("general", FORTSPARSE_SYM_GENERAL, 0, 300, nfail)
    call test_refactor("hermitian-indefinite", FORTSPARSE_SYM_HERMITIAN, 1, 300, nfail)
    call test_refactor("hermitian-posdef", FORTSPARSE_SYM_SPD, 2, 300, nfail)
    call test_refactor("complex-symmetric", FORTSPARSE_SYM_INDEFINITE, 3, 300, nfail)
    call test_block_and_inplace(nfail)
    call test_concurrent_objects(nfail)
    call test_singular(nfail)

    call test_mpi_finalize()
    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"

contains

    subroutine check_true(name, ok, nfail)
        character(*), intent(in)    :: name
        logical,      intent(in)    :: ok
        integer,      intent(inout) :: nfail

        if (.not. ok) then
            write (error_unit, "(a,a)") "FAIL: ", name
            nfail = nfail + 1
        end if
    end subroutine check_true

    ! Deterministic uniform number in (-1, 1); no global RNG state.
    function urand() result(u)
        real(dp) :: u

        seed_state = mod(seed_state*1103515245_8 + 12345_8, 2147483648_8)
        u = 2.0_dp*real(seed_state, dp)/2147483648.0_dp - 1.0_dp
    end function urand

    ! Sparse test matrix. kind 0: general complex; 1: Hermitian with diagonal
    ! of both signs (indefinite); 2: Hermitian positive definite; 3: complex
    ! symmetric (A = A^T) with a complex diagonal. Returns the CSC matrix and
    ! the dense copy used as the zgesv reference.
    subroutine build(n, kind, A, D)
        integer,       intent(in)  :: n, kind
        type(csc_z_t), intent(out) :: A
        complex(dp), allocatable, intent(out) :: D(:, :)

        integer, allocatable     :: r(:), c(:)
        complex(dp), allocatable :: v(:)
        type(fortsparse_status_t) :: st
        complex(dp) :: z
        integer :: i, j, k, m, nt

        allocate (D(n, n), source=(0.0_dp, 0.0_dp))
        do i = 1, n
            do m = 1, 4
                j = 1 + int((urand() + 1.0_dp)*0.5_dp*real(n, dp)*0.999_dp)
                if (kind == 0) then
                    D(i, j) = D(i, j) + cmplx(urand(), urand(), dp)
                else if (j < i) then
                    z = cmplx(urand(), urand(), dp)
                    D(i, j) = D(i, j) + z
                    if (kind == 3) then
                        D(j, i) = D(j, i) + z
                    else
                        D(j, i) = D(j, i) + conjg(z)
                    end if
                end if
            end do
        end do
        do i = 1, n
            select case (kind)
            case (0)
                D(i, i) = D(i, i) + cmplx(3.0_dp + urand(), urand(), dp)
            case (1)
                D(i, i) = cmplx((-1.0_dp)**i*(2.0_dp + urand()), 0.0_dp, dp)
            case (2)
                D(i, i) = cmplx(8.0_dp + urand(), 0.0_dp, dp)
            case (3)
                D(i, i) = cmplx(3.0_dp + urand(), 1.0_dp + urand(), dp)
            end select
        end do
        ! Hermitian kinds with a dominant diagonal stay nonsingular; make the
        ! indefinite one genuinely dominant so the dense reference is accurate.
        if (kind == 1) then
            do i = 1, n
                D(i, i) = D(i, i)*3.0_dp
            end do
        end if
        nt = count(D /= (0.0_dp, 0.0_dp))
        allocate (r(nt), c(nt), v(nt))
        k = 0
        do j = 1, n
            do i = 1, n
                if (D(i, j) /= (0.0_dp, 0.0_dp)) then
                    k = k + 1
                    r(k) = i
                    c(k) = j
                    v(k) = D(i, j)
                end if
            end do
        end do
        call csc_from_triplet(n, n, r, c, v, A, st)
    end subroutine build

    subroutine zgesv_ref(D, B, X)
        complex(dp), intent(in)  :: D(:, :), B(:, :)
        complex(dp), allocatable, intent(out) :: X(:, :)

        complex(dp), allocatable :: M(:, :)
        integer, allocatable :: ipiv(:)
        integer :: info, n

        interface
            subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
                import :: dp
                integer :: n, nrhs, lda, ldb, info, ipiv(*)
                complex(dp) :: a(lda, *), b(ldb, *)
            end subroutine zgesv
        end interface

        n = size(D, 1)
        M = D
        X = B
        allocate (ipiv(n))
        call zgesv(n, size(B, 2), M, n, ipiv, X, n, info)
        if (info /= 0) error stop "zgesv failed"
    end subroutine zgesv_ref

    function rhs(n, nrhs) result(B)
        integer, intent(in) :: n, nrhs
        complex(dp) :: B(n, nrhs)
        integer :: i, j

        do j = 1, nrhs
            do i = 1, n
                B(i, j) = cmplx(urand(), urand(), dp)
            end do
        end do
    end function rhs

    subroutine make_solver(solver, sym, refine, nthreads)
        type(sparse_solver_t), intent(out) :: solver
        integer, intent(in) :: sym, nthreads
        logical, intent(in) :: refine

        solver%backend_id = FORTSPARSE_BACKEND_MUMPS
        solver%refine = refine
        solver%options%symmetry = sym
        solver%options%nthreads = nthreads
    end subroutine make_solver

    ! Factor once, solve several right-hand sides, compare with zgesv; both
    ! with and without refinement, and for two thread counts.
    subroutine test_kind(name, sym, kind, n, nfail)
        character(*), intent(in)    :: name
        integer,      intent(in)    :: sym, kind, n
        integer,      intent(inout) :: nfail

        type(csc_z_t) :: A
        complex(dp), allocatable :: D(:, :), B(:, :), Xref(:, :), x(:)
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        real(dp) :: err
        integer :: j, ithr, iref
        logical :: refine

        call build(n, kind, A, D)
        B = rhs(n, 5)
        call zgesv_ref(D, B, Xref)
        allocate (x(n))
        do iref = 0, 1
            refine = (iref == 1)
            do ithr = 1, 2
                call make_solver(solver, sym, refine, ithr)
                call sparse_factor(solver, A, status)
                call check_true(name//" factor", status_ok(status), nfail)
                if (.not. status_ok(status)) return
                err = 0.0_dp
                do j = 1, 5
                    call sparse_solve(solver, B(:, j), x, status)
                    call check_true(name//" solve", status_ok(status), nfail)
                    err = max(err, maxval(abs(x - Xref(:, j))) &
                        /maxval(abs(Xref(:, j))))
                end do
                write (*, "(a,a,l2,i2,es10.2)") "  ", name, refine, ithr, err
                call check_true(name//" matches zgesv", err < 1.0e-10_dp, nfail)
                if (refine) call check_true(name//" refined to roundoff", &
                    err < 1.0e-12_dp, nfail)
                call sparse_destroy(solver)
            end do
        end do
    end subroutine test_kind

    ! New values on the same pattern (off-diagonals scaled, real diagonal
    ! shifted, so Hermitian and complex symmetry are kept) reuse the analysis
    ! and match zgesv; an independent matrix of the same size has another
    ! pattern, which is detected and handled by a full analysis.
    subroutine test_refactor(name, sym, kind, n, nfail)
        character(*), intent(in)    :: name
        integer,      intent(in)    :: sym, kind, n
        integer,      intent(inout) :: nfail

        type(csc_z_t) :: A, A2, A3
        complex(dp), allocatable :: D(:, :), D2(:, :), D3(:, :), B(:, :), &
            Xref(:, :), x(:)
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        logical :: reused
        integer :: i, j, e
        real(dp) :: err

        call build(n, kind, A, D)
        A2 = A
        D2 = 0.7_dp*D
        do j = 1, n
            do e = A%col_ptr(j), A%col_ptr(j + 1) - 1
                if (A%row_idx(e) == j) then
                    A2%val(e) = 0.7_dp*A%val(e) + 2.5_dp
                else
                    A2%val(e) = 0.7_dp*A%val(e)
                end if
            end do
            D2(j, j) = D2(j, j) + 2.5_dp
        end do
        B = rhs(n, 3)
        allocate (x(n))
        call make_solver(solver, sym, .true., 2)
        call sparse_factor(solver, A, status)
        call check_true(name//" refactor: initial factor", status_ok(status), nfail)
        if (.not. status_ok(status)) return
        call sparse_refactor(solver, A2, status, reused)
        call check_true(name//" refactor status", status_ok(status), nfail)
        call check_true(name//" refactor reuses analysis", reused, nfail)
        call zgesv_ref(D2, B, Xref)
        err = 0.0_dp
        do i = 1, 3
            call sparse_solve(solver, B(:, i), x, status)
            err = max(err, maxval(abs(x - Xref(:, i)))/maxval(abs(Xref(:, i))))
        end do
        call check_true(name//" refactor matches zgesv", err < 1.0e-10_dp, nfail)
        call build(n, kind, A3, D3)
        call sparse_refactor(solver, A3, status, reused)
        call check_true(name//" pattern change status", status_ok(status), nfail)
        call check_true(name//" pattern change detected", .not. reused, nfail)
        call zgesv_ref(D3, B, Xref)
        call sparse_solve(solver, B(:, 1), x, status)
        call check_true(name//" pattern change matches zgesv", &
            maxval(abs(x - Xref(:, 1))) < 1.0e-10_dp*maxval(abs(Xref(:, 1))), nfail)
        call sparse_destroy(solver)
        call test_raw_ownership(name, sym, A, A2, D, D2, B, nfail)
    end subroutine test_refactor

    ! Exercise the raw entry points against the same independent dense oracle.
    subroutine test_raw_ownership(name, sym, A, A2, D, D2, B, nfail)
        character(*), intent(in) :: name
        integer, intent(in) :: sym
        type(csc_z_t), intent(in) :: A, A2
        complex(dp), intent(in) :: D(:, :), D2(:, :), B(:, :)
        integer, intent(inout) :: nfail

        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        integer, allocatable :: ptr(:), idx(:)
        complex(dp), allocatable :: val(:), x(:), Xref(:, :)
        integer :: n, nz, iref, stride, phase
        logical :: reused

        n = A%ncol
        nz = A%nnz
        allocate (x(n))
        do iref = 0, 1
            do stride = 1, 2
                call make_solver(solver, sym, iref == 1, 2)
                ! Exercise refinement off -> on -> off on the retained analysis.
                do phase = 1, 3
                    allocate (ptr(stride*(n + 2)), idx(stride*(nz + 1)), &
                        val(stride*(nz + 1)))
                    ptr = -1; idx = -1; val = (-999.0_dp, 0.0_dp)
                    ptr(1:stride*n + 1:stride) = A%col_ptr
                    idx(1:stride*(nz - 1) + 1:stride) = A%row_idx
                    if (phase == 1) then
                        val(1:stride*(nz - 1) + 1:stride) = A%val
                        call zgesv_ref(D, B(:, 1:1), Xref)
                        call sparse_factor(solver, n, n, nz, ptr(1::stride), &
                            idx(1::stride), val(1::stride), status)
                    else
                        solver%refine = (phase == 2)
                        val(1:stride*(nz - 1) + 1:stride) = A2%val
                        call zgesv_ref(D2, B(:, 1:1), Xref)
                        call sparse_refactor(solver, n, n, nz, ptr(1::stride), &
                            idx(1::stride), val(1::stride), status, reused)
                        call check_true(name//" raw refactor reuses saved pattern", &
                            reused, nfail)
                    end if
                    call check_true(name//" raw factor/refactor", &
                        status_ok(status), nfail)
                    ptr = 0; idx = 0; val = (0.0_dp, 0.0_dp)
                    deallocate (ptr, idx, val)
                    if (.not. status_ok(status)) exit
                    call sparse_solve(solver, B(:, 1), x, status)
                    call check_true(name//" raw solve after input release", &
                        status_ok(status), nfail)
                    if (status_ok(status)) call check_true(name//" raw matches zgesv", &
                        maxval(abs(x - Xref(:, 1))) &
                        < 1.0e-10_dp*maxval(abs(Xref(:, 1))), nfail)
                end do
                call sparse_destroy(solver)
            end do
        end do
    end subroutine test_raw_ownership

    subroutine test_block_and_inplace(nfail)
        integer, intent(inout) :: nfail

        type(csc_z_t) :: A
        complex(dp), allocatable :: D(:, :), B(:, :), Xref(:, :), X(:, :), &
            x1(:), b1(:)
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        integer, parameter :: n = 200

        call build(n, 0, A, D)
        B = rhs(n, 7)
        call zgesv_ref(D, B, Xref)
        call make_solver(solver, FORTSPARSE_SYM_GENERAL, .true., 2)
        call sparse_factor(solver, A, status)
        call check_true("block factor", status_ok(status), nfail)
        allocate (X(n, 7), x1(n))
        call sparse_solve_block(solver, B, X, status)
        call check_true("block solve", status_ok(status), nfail)
        call check_true("block vs zgesv", &
            maxval(abs(X - Xref)) < 1.0e-12_dp*maxval(abs(Xref)), nfail)
        b1 = B(:, 3)
        call sparse_solve(solver, b1, status) ! in place
        call check_true("in-place solve", status_ok(status), nfail)
        call check_true("in-place vs block", &
            maxval(abs(b1 - X(:, 3))) < 1.0e-12_dp*maxval(abs(Xref)), nfail)
        call sparse_destroy(solver)
        ! A freed solver refuses further solves instead of crashing.
        call sparse_solve(solver, B(:, 1), x1, status)
        call check_true("solve after destroy fails", .not. status_ok(status), nfail)
    end subroutine test_block_and_inplace

    ! Independent solver objects used from OpenMP threads; the shim serializes
    ! MUMPS calls, so the results stay exact.
    subroutine test_concurrent_objects(nfail)
        integer, intent(inout) :: nfail

        integer, parameter :: n = 120, nobj = 4
        type(csc_z_t) :: A(nobj)
        complex(dp), allocatable :: D(:, :, :), B(:, :), Xref(:, :, :)
        complex(dp), allocatable :: x(:, :)
        complex(dp), allocatable :: Dk(:, :), Xk(:, :)
        type(sparse_solver_t) :: solver(nobj)
        type(fortsparse_status_t) :: status(nobj)
        integer :: k
        logical :: ok(nobj)

        allocate (D(n, n, nobj), Xref(n, 1, nobj), x(n, nobj))
        B = rhs(n, 1)
        do k = 1, nobj
            call build(n, 0, A(k), Dk)
            D(:, :, k) = Dk
            call zgesv_ref(Dk, B, Xk)
            Xref(:, :, k) = Xk
        end do
        ok = .false.
        !$omp parallel do private(k) schedule(static, 1)
        do k = 1, nobj
            call make_solver(solver(k), FORTSPARSE_SYM_GENERAL, .true., 1)
            call sparse_factor(solver(k), A(k), status(k))
            if (status_ok(status(k))) &
                call sparse_solve(solver(k), B(:, 1), x(:, k), status(k))
            ok(k) = status_ok(status(k))
        end do
        !$omp end parallel do
        do k = 1, nobj
            call check_true("concurrent object status", ok(k), nfail)
            if (ok(k)) call check_true("concurrent object result", &
                maxval(abs(x(:, k) - Xref(:, 1, k))) &
                < 1.0e-12_dp*maxval(abs(Xref(:, 1, k))), nfail)
            call sparse_destroy(solver(k))
        end do
    end subroutine test_concurrent_objects

    ! A structurally singular matrix (empty last column) is reported as
    ! FORTSPARSE_SINGULAR, for the general and the Hermitian path.
    subroutine test_singular(nfail)
        integer, intent(inout) :: nfail

        type(csc_z_t) :: A
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        integer, parameter :: r(4) = [1, 1, 2, 2], c(4) = [1, 2, 1, 2]
        complex(dp), parameter :: v(4) = [(1.0_dp, 1.0_dp), (2.0_dp, 2.0_dp), &
            (2.0_dp, -2.0_dp), (4.0_dp, 4.0_dp)]

        call csc_from_triplet(3, 3, r, c, v, A, status)
        call make_solver(solver, FORTSPARSE_SYM_GENERAL, .false., 1)
        call sparse_factor(solver, A, status)
        call check_true("singular complex detected (general)", &
            status%code == FORTSPARSE_SINGULAR, nfail)
        call make_solver(solver, FORTSPARSE_SYM_HERMITIAN, .false., 1)
        call sparse_factor(solver, A, status)
        call check_true("singular complex detected (Hermitian)", &
            status%code == FORTSPARSE_SINGULAR, nfail)
    end subroutine test_singular

end program test_fortsparse_mumps_complex
