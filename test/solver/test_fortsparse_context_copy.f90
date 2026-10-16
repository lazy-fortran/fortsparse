program test_fortsparse_context_copy
    ! Public context copy contract for native CPU backends. Every copy must
    ! solve manufactured systems against x_true after its source is destroyed,
    ! own its storage (no aliasing, no double free) and respect its host cap.
    use fortsparse
    implicit none
    integer :: nfail, amd_order
    nfail = 0
    ! fpm builds without the AMD library; pivoting is still exercised.
    amd_order = FORTSPARSE_ORDER_DEFAULT
    if (native_ordering_available(FORTSPARSE_ORDER_AMD)) &
        amd_order = FORTSPARSE_ORDER_AMD
    call spd_cases()
    call lu_cases()
    call complex_cases()
    call unfactored_cases()
    call cap_cases()
    call loop_cases()
    if (nfail > 0) then
        write (*, '(i0,a)') nfail, ' check(s) failed'
        error stop 1
    end if
    write (*, '(a)') 'PASS'
contains

    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            nfail = nfail + 1
            write (*, '(2a)') 'FAIL: ', label
        end if
    end subroutine check

    subroutine real_csc(d, a)
        real(dp), intent(in) :: d(:,:)
        type(csc_t), intent(out) :: a
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        type(fortsparse_status_t) :: status
        integer :: i, j, k
        k = count(d /= 0.0_dp)
        allocate (rows(k), cols(k), vals(k))
        k = 0
        do j = 1, size(d, 2)
            do i = 1, size(d, 1)
                if (d(i,j) == 0.0_dp) cycle
                k = k + 1
                rows(k) = i
                cols(k) = j
                vals(k) = d(i,j)
            end do
        end do
        call csc_from_triplet(size(d,1), size(d,2), rows, cols, vals, a, status)
        if (.not. status_ok(status)) error stop 'real CSC construction'
    end subroutine real_csc

    subroutine complex_csc(d, a)
        complex(dp), intent(in) :: d(:,:)
        type(csc_z_t), intent(out) :: a
        integer, allocatable :: rows(:), cols(:)
        complex(dp), allocatable :: vals(:)
        type(fortsparse_status_t) :: status
        integer :: i, j, k
        k = count(d /= (0.0_dp, 0.0_dp))
        allocate (rows(k), cols(k), vals(k))
        k = 0
        do j = 1, size(d, 2)
            do i = 1, size(d, 1)
                if (d(i,j) == (0.0_dp, 0.0_dp)) cycle
                k = k + 1
                rows(k) = i
                cols(k) = j
                vals(k) = d(i,j)
            end do
        end do
        call csc_from_triplet(size(d,1), size(d,2), rows, cols, vals, a, status)
        if (.not. status_ok(status)) error stop 'complex CSC construction'
    end subroutine complex_csc

    ! Solve op(D) x = b for x_true with vectors and a 3-column block.
    subroutine check_real(solver, d, label)
        type(sparse_solver_t), intent(inout) :: solver
        real(dp), intent(in) :: d(:,:)
        character(*), intent(in) :: label
        real(dp) :: xe(size(d,1)), b(size(d,1)), x(size(d,1))
        real(dp) :: bb(size(d,1),3), xb(size(d,1),3)
        type(fortsparse_status_t) :: status
        integer :: i, op
        call check(solver%factored, label//': factored')
        xe = [(sin(real(i, dp)) + 0.5_dp, i = 1, size(d,1))]
        do op = FORTSPARSE_OP_NORMAL, FORTSPARSE_OP_TRANSPOSE
            if (op == FORTSPARSE_OP_NORMAL) then
                b = matmul(d, xe)
            else
                b = matmul(transpose(d), xe)
            end if
            call sparse_solve_op(solver, b, x, op, status)
            call check(status_ok(status), label//': solve status')
            call check(maxval(abs(x - xe)) < 1e-10_dp, label//': x_true')
            do i = 1, 3
                bb(:,i) = real(i, dp)*b
            end do
            call sparse_solve_block_op(solver, bb, xb, op, status)
            call check(status_ok(status), label//': block status')
            call check(maxval(abs(xb(:,3) - 3.0_dp*xe)) < 1e-9_dp, &
                label//': block x_true')
        end do
    end subroutine check_real

    subroutine check_complex(solver, d, label)
        type(sparse_solver_t), intent(inout) :: solver
        complex(dp), intent(in) :: d(:,:)
        character(*), intent(in) :: label
        complex(dp) :: xe(size(d,1)), b(size(d,1)), x(size(d,1))
        complex(dp) :: bb(size(d,1),2), xb(size(d,1),2)
        type(fortsparse_status_t) :: status
        integer :: i, op
        call check(solver%factored, label//': factored')
        xe = [(cmplx(cos(real(i, dp)), 0.3_dp*i, dp), i = 1, size(d,1))]
        do op = FORTSPARSE_OP_NORMAL, FORTSPARSE_OP_ADJOINT
            select case (op)
            case (FORTSPARSE_OP_NORMAL)
                b = matmul(d, xe)
            case (FORTSPARSE_OP_TRANSPOSE)
                b = matmul(transpose(d), xe)
            case default
                b = matmul(conjg(transpose(d)), xe)
            end select
            call sparse_solve_op(solver, b, x, op, status)
            call check(status_ok(status), label//': solve status')
            call check(maxval(abs(x - xe)) < 1e-10_dp, label//': x_true')
            bb(:,1) = b
            bb(:,2) = (0.0_dp, 2.0_dp)*b
            call sparse_solve_block_op(solver, bb, xb, op, status)
            call check(status_ok(status), label//': block status')
            call check(maxval(abs(xb(:,2) - (0.0_dp, 2.0_dp)*xe)) < 1e-9_dp, &
                label//': block x_true')
        end do
    end subroutine check_complex

    ! Copy src twice (assignment and sparse_copy), destroy src, then check
    ! both copies, destroy one and check the other again.
    subroutine copy_real(src, d, family, label)
        type(sparse_solver_t), intent(inout) :: src
        real(dp), intent(in) :: d(:,:)
        character(*), intent(in) :: family, label
        type(sparse_solver_t) :: c1, c2
        type(sparse_factor_stats_t) :: stats
        type(fortsparse_status_t) :: status
        call sparse_factor_stats(src, stats, status)
        call check(stats%native_factor_family == family, label//': source family')
        c1 = src
        call sparse_copy(c2, src, status)
        call check(status_ok(status), label//': sparse_copy status')
        call sparse_destroy(src)
        call check_real(c1, d, label//' assigned')
        call check_real(c2, d, label//' copied')
        call sparse_factor_stats(c1, stats, status)
        call check(stats%native_factor_family == family, label//': copy family')
        call sparse_destroy(c1)
        call check_real(c2, d, label//' survivor')
    end subroutine copy_real

    subroutine copy_complex(src, d, family, label)
        type(sparse_solver_t), intent(inout) :: src
        complex(dp), intent(in) :: d(:,:)
        character(*), intent(in) :: family, label
        type(sparse_solver_t) :: c1, c2
        type(sparse_factor_stats_t) :: stats
        type(fortsparse_status_t) :: status
        c1 = src
        call sparse_copy(c2, src, status)
        call check(status_ok(status), label//': sparse_copy status')
        call sparse_destroy(src)
        call check_complex(c1, d, label//' assigned')
        call sparse_factor_stats(c1, stats, status)
        call check(stats%native_factor_family == family, label//': copy family')
        call sparse_destroy(c1)
        call check_complex(c2, d, label//' survivor')
    end subroutine copy_complex

    subroutine spd_matrix(n, band, d)
        integer, intent(in) :: n, band
        real(dp), intent(out) :: d(n,n)
        integer :: i, j
        d = 0.0_dp
        do j = 1, n
            do i = max(1, j - band), min(n, j + band)
                if (i /= j) d(i,j) = -1.0_dp/real(1 + abs(i - j), dp)
            end do
            d(j,j) = 2.0_dp*band + 1.0_dp
        end do
        d(1,n) = -0.5_dp
        d(n,1) = -0.5_dp
    end subroutine spd_matrix

    subroutine spd_cases()
        type(sparse_solver_t) :: s, c
        type(sparse_factor_stats_t) :: st_src, st_copy
        type(fortsparse_status_t) :: status
        type(csc_t) :: a
        real(dp) :: d20(20,20), d80(80,80), b(20), x(20)
        real(dp), pointer :: p1(:), p2(:)
        integer :: i

        call spd_matrix(20, 1, d20)
        call real_csc(d20, a)
        s%backend_id = FORTSPARSE_BACKEND_NATIVE
        s%options%symmetry = FORTSPARSE_SYM_SPD
        call sparse_factor(s, a, status)
        call check(status_ok(status), 'SPD simplicial factor')
        call copy_real(s, d20, 'CHOLESKY_SIMPLE', 'SPD simplicial')

        ! Bounded LOW_MEMORY with refinement: budget and checker state copy.
        s%backend_id = FORTSPARSE_BACKEND_NATIVE
        s%options%symmetry = FORTSPARSE_SYM_SPD
        s%options%memory_policy = FORTSPARSE_MEMORY_LOW
        s%refine = .true.
        call sparse_factor(s, a, status)
        call check(status_ok(status), 'SPD LOW refined factor')
        p1 => sparse_vector(s, 20)
        c = s
        p2 => sparse_vector(c, 20)
        call check(associated(p1) .and. associated(p2), 'SPD vectors')
        if (associated(p1) .and. associated(p2)) &
            call check(.not. associated(p1, p2), 'SPD vectors do not alias')
        call sparse_factor_stats(s, st_src, status)
        call sparse_factor_stats(c, st_copy, status)
        call check(st_copy%strict_memory_supported, 'SPD LOW copy strict')
        call check(st_copy%host_current_bytes == st_src%host_current_bytes, &
            'SPD LOW copy charges retained payload')
        call sparse_destroy(s)
        p2 = matmul(d20, [(real(i, dp), i = 1, 20)])
        call sparse_solve(c, p2, x, status)
        call check(status_ok(status), 'SPD LOW refined copy solve')
        call check(maxval(abs(x - [(real(i, dp), i = 1, 20)])) < 1e-10_dp, &
            'SPD LOW refined copy x_true')
        call check_real(c, d20, 'SPD LOW refined copy')
        call sparse_destroy(c)
        s = sparse_solver_t()

        call spd_matrix(80, 15, d80)
        call real_csc(d80, a)
        s%backend_id = FORTSPARSE_BACKEND_NATIVE
        s%options%symmetry = FORTSPARSE_SYM_SPD
        call sparse_factor(s, a, status)
        call check(status_ok(status), 'SPD supernodal factor')
        call copy_real(s, d80, 'CHOLESKY_BLOCK', 'SPD supernodal')

        ! Self-assignment and self-copy keep the factors.
        s%backend_id = FORTSPARSE_BACKEND_NATIVE
        s%options%symmetry = FORTSPARSE_SYM_SPD
        call sparse_factor(s, a, status)
        s = s
        call check_real(s, d80, 'SPD self-assignment')
        call sparse_copy(s, s, status)
        call check(status_ok(status), 'SPD self-copy status')
        call check_real(s, d80, 'SPD self-copy')
        call sparse_destroy(s)
        b = 0.0_dp
        call sparse_solve(s, b, x, status)
        call check(status%code == FORTSPARSE_NOT_FACTORED, 'destroyed source')
    end subroutine spd_cases

    subroutine lu_cases()
        type(sparse_solver_t) :: s, c
        type(fortsparse_status_t) :: status
        type(csc_t) :: a, adot
        real(dp) :: d(12,12), xe(12), x(12), dx(12)
        integer :: i
        d = 0.0_dp
        do i = 1, 12
            d(i,i) = 4.0_dp + real(i, dp)
            if (i > 1) d(i,i-1) = 0.5_dp
            if (i < 11) d(i,i+2) = -0.2_dp*real(i, dp)
        end do
        d(1,1) = 0.0_dp
        d(1,2) = 1.0_dp
        d(2,1) = 3.0_dp
        d(12,1) = 0.7_dp
        call real_csc(d, a)
        s%backend_id = FORTSPARSE_BACKEND_NATIVE
        s%options%ordering = amd_order
        call sparse_factor(s, a, status)
        call check(status_ok(status), 'real LU AMD pivoting factor')
        c = s
        call copy_real(s, d, 'LU_REAL', 'real LU AMD')
        ! Derivative on a copy whose source is gone: A_dot = A/8, dx = -x/8.
        xe = [(1.0_dp + 0.1_dp*i, i = 1, 12)]
        call sparse_solve(c, matmul(d, xe), x, status)
        adot = a
        adot%val = 0.125_dp*a%val
        call sparse_solve_jvp(c, adot, x, [(0.0_dp, i = 1, 12)], dx, status)
        call check(status_ok(status), 'LU copy JVP status')
        call check(maxval(abs(dx + 0.125_dp*xe)) < 1e-10_dp, 'LU copy JVP oracle')

        ! Real symmetric indefinite LDLT with a zero-diagonal 2x2 pivot.
        d = 0.0_dp
        do i = 1, 12
            d(i,i) = 3.0_dp*(-1.0_dp)**i
            if (i > 1) then
                d(i,i-1) = 1.0_dp
                d(i-1,i) = 1.0_dp
            end if
        end do
        d(1,1) = 0.0_dp
        d(1,2) = 2.0_dp
        d(2,1) = 2.0_dp
        call real_csc(d, a)
        s%backend_id = FORTSPARSE_BACKEND_NATIVE
        s%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        call sparse_factor(s, a, status)
        call check(status_ok(status), 'real LDLT factor')
        call copy_real(s, d, 'LDLT_REAL', 'real LDLT')
    end subroutine lu_cases

    subroutine complex_cases()
        type(sparse_solver_t) :: s
        type(fortsparse_status_t) :: status
        type(csc_z_t) :: a
        complex(dp) :: d(12,12)
        integer :: i
        ! Complex general LU.
        d = (0.0_dp, 0.0_dp)
        do i = 1, 12
            d(i,i) = cmplx(4.0_dp + i, 1.0_dp, dp)
            if (i > 1) d(i,i-1) = (0.3_dp, -0.2_dp)
            if (i < 10) d(i,i+3) = (-0.1_dp, 0.4_dp)
        end do
        d(1,1) = (0.0_dp, 0.0_dp)
        d(2,1) = (2.0_dp, 1.0_dp)
        d(1,2) = (1.0_dp, -1.0_dp)
        call complex_csc(d, a)
        s%backend_id = FORTSPARSE_BACKEND_NATIVE
        call sparse_factor(s, a, status)
        call check(status_ok(status), 'complex LU factor')
        call copy_complex(s, d, 'LU_COMPLEX', 'complex LU')

        ! Hermitian positive definite.
        d = (0.0_dp, 0.0_dp)
        do i = 1, 12
            d(i,i) = (6.0_dp, 0.0_dp)
            if (i > 1) then
                d(i,i-1) = (-1.0_dp, 0.5_dp)
                d(i-1,i) = conjg(d(i,i-1))
            end if
        end do
        call complex_csc(d, a)
        s%backend_id = FORTSPARSE_BACKEND_NATIVE
        s%options%symmetry = FORTSPARSE_SYM_SPD
        call sparse_factor(s, a, status)
        call check(status_ok(status), 'HPD factor')
        call copy_complex(s, d, 'CHOLESKY_HPD', 'HPD')

        ! Complex symmetric (A = A^T) indefinite LDLT.
        d = (0.0_dp, 0.0_dp)
        do i = 1, 12
            d(i,i) = cmplx(3.0_dp*(-1.0_dp)**i, 1.0_dp, dp)
            if (i > 1) then
                d(i,i-1) = (1.0_dp, 0.5_dp)
                d(i-1,i) = d(i,i-1)
            end if
        end do
        call complex_csc(d, a)
        s%backend_id = FORTSPARSE_BACKEND_NATIVE
        s%options%symmetry = FORTSPARSE_SYM_INDEFINITE
        call sparse_factor(s, a, status)
        call check(status_ok(status), 'complex LDLT factor')
        call copy_complex(s, d, 'LDLT_COMPLEX', 'complex LDLT')
    end subroutine complex_cases

    subroutine unfactored_cases()
        type(sparse_solver_t) :: s, c
        type(fortsparse_status_t) :: status
        type(csc_t) :: a
        real(dp) :: d(20,20)
        call spd_matrix(20, 2, d)
        call real_csc(d, a)
        s%backend_id = FORTSPARSE_BACKEND_NATIVE
        s%refine = .true.
        s%options%symmetry = FORTSPARSE_SYM_SPD
        s%options%ordering = amd_order
        call sparse_copy(c, s, status)
        call check(status_ok(status), 'unfactored copy status')
        call check(.not. c%factored .and. .not. allocated(c%backend), &
            'unfactored copy has no backend')
        call check(c%backend_id == FORTSPARSE_BACKEND_NATIVE .and. c%refine .and. &
            c%options%ordering == amd_order, 'unfactored copy options')
        call sparse_factor(c, a, status)
        call check(status_ok(status), 'unfactored copy factor')
        call check_real(c, d, 'unfactored copy refactored')
        ! A freed (allocated, unfactored) backend copies as unfactored too.
        call sparse_free(c)
        s = c
        call check(.not. s%factored .and. allocated(s%backend), 'freed backend copy')
        call sparse_factor(s, a, status)
        call check_real(s, d, 'freed copy factor')
    end subroutine unfactored_cases

    subroutine cap_cases()
        type(sparse_solver_t) :: s, c
        type(sparse_factor_stats_t) :: stats
        type(fortsparse_status_t) :: status
        type(csc_t) :: a
        real(dp) :: d(20,20), b(20), x(20)
        integer(i8) :: live
        call spd_matrix(20, 3, d)
        d(1,1) = 0.0_dp
        d(5,2) = 1.5_dp
        call real_csc(d, a)
        s%backend_id = FORTSPARSE_BACKEND_NATIVE
        call sparse_factor(s, a, status)
        call check(status_ok(status), 'cap source factor')
        call sparse_factor_stats(s, stats, status)
        live = stats%host_current_bytes
        call check(live > 64_i8, 'cap source payload')
        ! Factor c first so the failed copy must also drop old factors.
        c%backend_id = FORTSPARSE_BACKEND_NATIVE
        call sparse_factor(c, a, status)
        c%options%host_memory_limit_bytes = 64_i8
        call sparse_copy(c, s, status)
        call check(status%code == FORTSPARSE_RESOURCE_LIMIT, 'too-small cap refused')
        call check(.not. c%factored, 'refused copy unfactored')
        call check(c%options%host_memory_limit_bytes == 64_i8, 'destination cap kept')
        call sparse_factor_stats(c, stats, status)
        call check(stats%budget_required_bytes == live .and. &
            stats%budget_available_bytes == 64_i8, 'refused copy diagnostics')
        b = 1.0_dp
        call sparse_solve(c, b, x, status)
        call check(status%code == FORTSPARSE_NOT_FACTORED, 'refused copy solve')
        c%options%host_memory_limit_bytes = live
        call sparse_copy(c, s, status)
        call check(status_ok(status), 'exact cap copy')
        call check_real(c, d, 'exact cap copy')
        c%options%host_memory_limit_bytes = 0_i8
        call sparse_factor(c, a, status)
        call check(status_ok(status), 'destination refactor after refusal')
        call check_real(c, d, 'destination refactor')
    end subroutine cap_cases

    subroutine loop_cases()
        type(sparse_solver_t) :: s, c
        type(sparse_solver_t) :: many(3)
        type(fortsparse_status_t) :: status
        type(csc_t) :: a
        type(csc_z_t) :: z
        real(dp) :: d(80,80)
        complex(dp) :: dz(80,80)
        integer :: it, k
        call spd_matrix(80, 15, d)
        call real_csc(d, a)
        dz = cmplx(d, 0.0_dp, dp)
        dz(3,1) = (0.0_dp, 0.25_dp)
        call complex_csc(dz, z)
        s%backend_id = FORTSPARSE_BACKEND_NATIVE
        s%options%symmetry = FORTSPARSE_SYM_SPD
        call sparse_factor(s, a, status)
        do it = 1, 40
            c = s
            call sparse_copy(c, s, status)
            call check(status_ok(status), 'loop copy')
            s = c
            if (mod(it, 2) == 0) call sparse_destroy(c)
        end do
        many = s
        do k = 1, size(many)
            call check_real(many(k), d, 'elemental assignment')
        end do
        call sparse_destroy(s)
        call sparse_destroy(many(1))
        call check_real(many(3), d, 'elemental copy survivor')
        s%backend_id = FORTSPARSE_BACKEND_NATIVE
        s%options%symmetry = FORTSPARSE_SYM_GENERAL
        call sparse_factor(s, z, status)
        call check(status_ok(status), 'loop complex factor')
        do it = 1, 40
            call sparse_copy(c, s, status)
            s = c
            call sparse_free(c)
        end do
        call check_complex(s, dz, 'complex loop copy')
    end subroutine loop_cases
end program test_fortsparse_context_copy
