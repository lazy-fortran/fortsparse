program bench_block_tridiagonal
    ! Symmetric indefinite block-tridiagonal pencils with dense blocks, the
    ! shape of radial finite-element stability matrices. A = S L D L^T S with
    ! unit lower block-bidiagonal L, diagonal D of prescribed signs and a
    ! diagonal scaling S, so Sylvester's law gives the exact inertia. Times
    ! the native factorization against a dense block elimination with LAPACK
    ! DSYTRF/DSYTRS/DGEMM on the same blocks.
    use, intrinsic :: iso_fortran_env, only: int64
    use fortsparse
    !$ use omp_lib, only: omp_get_max_threads
    implicit none
    integer :: nb, w, repeats, n, nnz, negatives, threads, r, info
    real(dp), allocatable :: diag(:,:,:), lower(:,:,:), b(:), x(:)
    type(csc_t) :: a
    real(dp) :: decades, t_factor, t_refactor, t_solve, t_block, residual
    character(len=32) :: argument, label
    type(sparse_solver_t) :: solver
    type(fortsparse_status_t) :: status
    type(sparse_factor_stats_t) :: stats
    integer(i8) :: nneg, nzero, npos
    integer :: reliability, block_negative
    logical :: reused

    if (command_argument_count() /= 5) error stop &
        'usage: bench_block_tridiagonal BLOCKS WIDTH DECADES MODE REPEATS'
    call get_command_argument(1, argument); read (argument, *) nb
    call get_command_argument(2, argument); read (argument, *) w
    call get_command_argument(3, argument); read (argument, *) decades
    call get_command_argument(4, label)
    call get_command_argument(5, argument); read (argument, *) repeats
    n = nb*w
    call build(nb, w, decades, diag, lower, negatives)
    call to_csc(nb, w, diag, lower, a)
    nnz = a%nnz
    threads = 1
    !$ threads = omp_get_max_threads()
    allocate (b(n), x(n))
    call random_number(b)

    solver%backend_id = FORTSPARSE_BACKEND_NATIVE
    solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
    select case (trim(label))
    case ('scalar')
    case ('blocked')
        solver%options%ldlt_blocked = .true.
    case ('supernodal')
        solver%options%ldlt_supernodal = .true.
    case ('supernodal-natural')
        solver%options%ldlt_supernodal = .true.
        solver%options%ordering = FORTSPARSE_ORDER_NATURAL
    case ('supernodal-nd')
        solver%options%ldlt_supernodal = .true.
        solver%options%ordering = FORTSPARSE_ORDER_ND
    case default
        error stop 'MODE: scalar|blocked|supernodal|supernodal-natural|supernodal-nd'
    end select

    t_factor = huge(1.0_dp)
    t_refactor = huge(1.0_dp)
    t_solve = huge(1.0_dp)
    t_block = huge(1.0_dp)
    do r = 1, repeats
        call sparse_destroy(solver)
        t_factor = min(t_factor, timed_factor(.false.))
        if (.not. status_ok(status)) then
            print '(a)', trim(status%msg)
            error stop 'factor failed'
        end if
        t_refactor = min(t_refactor, timed_factor(.true.))
        if (.not. status_ok(status)) error stop 'refactor failed'
        t_solve = min(t_solve, timed_solve())
        t_block = min(t_block, timed_block())
    end do
    call sparse_inertia(solver, nneg, nzero, npos, reliability, status)
    call sparse_factor_stats(solver, stats, status)
    residual = backward_error()
    print '(a)', 'mode,threads,blocks,width,n,nnz,decades,factor_s,refactor_s,' &
        //'analysis_s,numeric_s,solve_s,block_lapack_s,negative,expected,' &
        //'block_negative,reliability,uncertain,ordering,backward_error'
    print '(a,",",5(i0,","),f0.1,",",6(es11.4,","),6(i0,","),es10.3)', &
        trim(label), threads, nb, w, n, nnz, decades, t_factor, t_refactor, &
        stats%analysis_s, stats%factor_s, t_solve, t_block, nneg, negatives, &
        block_negative, reliability, stats%inertia_uncertain, &
        stats%ordering_used, residual
    if (nneg /= negatives .or. nzero /= 0) error stop 'inertia mismatch'
    if (.not. residual < 1.0e-9_dp) error stop 'backward error'

contains

    real(dp) function seconds()
        integer(int64) :: count, rate
        call system_clock(count, rate)
        seconds = real(count, dp)/real(rate, dp)
    end function seconds

    real(dp) function timed_factor(again) result(elapsed)
        logical, intent(in) :: again
        real(dp) :: start
        start = seconds()
        if (again) then
            call sparse_refactor(solver, a, status, reused)
        else
            call sparse_factor(solver, a, status)
        end if
        elapsed = seconds() - start
    end function timed_factor

    real(dp) function timed_solve() result(elapsed)
        real(dp) :: start
        x = b
        start = seconds()
        call sparse_solve(solver, x, status)
        elapsed = seconds() - start
        if (.not. status_ok(status)) error stop 'solve failed'
    end function timed_solve

    ! max_i |A x - b|_i / (|A| |x| + |b|)_i against the input matrix.
    real(dp) function backward_error() result(omega)
        real(dp), allocatable :: residual_vector(:), bound(:)
        integer :: j, p
        allocate (residual_vector(n), bound(n))
        residual_vector = -b
        bound = abs(b)
        do j = 1, n
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                residual_vector(a%row_idx(p)) = residual_vector(a%row_idx(p)) &
                    + a%val(p)*x(j)
                bound(a%row_idx(p)) = bound(a%row_idx(p)) + abs(a%val(p)*x(j))
            end do
        end do
        omega = maxval(abs(residual_vector)/bound)
    end function backward_error

    ! Reference: S_i = A_ii - A_{i,i-1} S_{i-1}^{-1} A_{i,i-1}^T with pivots
    ! restricted to each block (DSYTRF), as in a block-tridiagonal solver.
    real(dp) function timed_block() result(elapsed)
        real(dp), allocatable :: schur(:,:,:), coupled(:,:), work(:)
        integer, allocatable :: pivots(:,:)
        real(dp) :: start
        integer :: i, j
        external :: dsytrf, dsytrs, dgemm
        allocate (schur(w, w, nb), coupled(w, w), work(64*w), pivots(w, nb))
        start = seconds()
        schur = diag
        block_negative = 0
        do i = 1, nb
            if (i > 1) then
                coupled = transpose(lower(:, :, i - 1))
                call dsytrs('L', w, w, schur(1, 1, i - 1), w, pivots(1, i - 1), &
                    coupled, w, info)
                call dgemm('N', 'N', w, w, w, -1.0_dp, lower(1, 1, i - 1), w, &
                    coupled, w, 1.0_dp, schur(1, 1, i), w)
            end if
            call dsytrf('L', w, schur(1, 1, i), w, pivots(1, i), work, size(work), info)
            if (info /= 0) error stop 'block reference failed'
            j = 1
            do while (j <= w)
                if (pivots(j, i) > 0) then
                    if (schur(j, j, i) < 0.0_dp) block_negative = block_negative + 1
                    j = j + 1
                else
                    ! A 2x2 pivot has one negative eigenvalue when its
                    ! determinant is negative, else two of the trace's sign.
                    if (schur(j, j, i)*schur(j + 1, j + 1, i) &
                        - schur(j + 1, j, i)**2 < 0.0_dp) then
                        block_negative = block_negative + 1
                    else if (schur(j, j, i) < 0.0_dp) then
                        block_negative = block_negative + 2
                    end if
                    j = j + 2
                end if
            end do
        end do
        elapsed = seconds() - start
    end function timed_block

    subroutine build(nb, w, decades, diag, lower, negatives)
        integer, intent(in) :: nb, w
        real(dp), intent(in) :: decades
        real(dp), allocatable, intent(out) :: diag(:,:,:), lower(:,:,:)
        integer, intent(out) :: negatives
        real(dp), allocatable :: l(:,:,:), c(:,:,:), d(:,:), s(:,:), t(:,:)
        integer, allocatable :: seed(:)
        integer :: i, j, k, m
        call random_seed(size=m)
        allocate (seed(m))
        seed = [(7919*k + 104729, k=1, m)]
        call random_seed(put=seed)
        allocate (diag(w, w, nb), lower(w, w, max(1, nb - 1)), l(w, w, nb), &
            c(w, w, nb), d(w, nb), s(w, nb), t(w, w))
        call random_number(l)
        call random_number(c)
        call random_number(d)
        call random_number(s)
        l = (l - 0.5_dp)/real(w, dp)
        c = (c - 0.5_dp)/real(w, dp)
        d = 0.5_dp + d
        negatives = 0
        do i = 1, nb
            do j = 1, w
                l(j, j, i) = 1.0_dp
                l(1:j - 1, j, i) = 0.0_dp
                ! One negative direction in every seventh block.
                if (mod(i, 7) == 3 .and. j == (w + 1)/2) then
                    d(j, i) = -d(j, i)
                    negatives = negatives + 1
                end if
            end do
        end do
        s = 10.0_dp**(decades*(s - 0.5_dp))
        do i = 1, nb
            do j = 1, w
                t(:, j) = l(:, j, i)*d(j, i)
            end do
            diag(:, :, i) = matmul(t, transpose(l(:, :, i)))
            if (i > 1) then
                do j = 1, w
                    t(:, j) = c(:, j, i - 1)*d(j, i - 1)
                end do
                diag(:, :, i) = diag(:, :, i) + matmul(t, transpose(c(:, :, i - 1)))
                lower(:, :, i - 1) = matmul(t, transpose(l(:, :, i - 1)))
            end if
        end do
        do i = 1, nb
            do j = 1, w
                do k = j, w
                    diag(k, j, i) = diag(k, j, i)*s(k, i)*s(j, i)
                    diag(j, k, i) = diag(k, j, i)
                end do
                if (i < nb) lower(:, j, i) = lower(:, j, i)*s(:, i + 1)*s(j, i)
            end do
        end do
    end subroutine build

    subroutine to_csc(nb, w, diag, lower, matrix)
        integer, intent(in) :: nb, w
        real(dp), intent(in) :: diag(:,:,:), lower(:,:,:)
        type(csc_t), intent(out) :: matrix
        integer, allocatable :: ptr(:), rows(:)
        real(dp), allocatable :: values(:)
        integer :: i, j, k, p
        allocate (ptr(nb*w + 1), rows(w*w*(3*nb - 2)), values(w*w*(3*nb - 2)))
        p = 1
        do i = 1, nb
            do j = 1, w
                ptr((i - 1)*w + j) = p
                if (i > 1) then
                    do k = 1, w
                        rows(p) = (i - 2)*w + k
                        values(p) = lower(j, k, i - 1)
                        p = p + 1
                    end do
                end if
                do k = 1, w
                    rows(p) = (i - 1)*w + k
                    values(p) = diag(k, j, i)
                    p = p + 1
                end do
                if (i < nb) then
                    do k = 1, w
                        rows(p) = i*w + k
                        values(p) = lower(k, j, i)
                        p = p + 1
                    end do
                end if
            end do
        end do
        ptr(nb*w + 1) = p
        matrix%nrow = nb*w
        matrix%ncol = nb*w
        matrix%nnz = p - 1
        call move_alloc(ptr, matrix%col_ptr)
        call move_alloc(rows, matrix%row_idx)
        call move_alloc(values, matrix%val)
    end subroutine to_csc
end program bench_block_tridiagonal
