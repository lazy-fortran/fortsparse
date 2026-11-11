program bench_block_tridiagonal
    ! Symmetric indefinite block-tridiagonal pencils with dense blocks, the
    ! shape of radial finite-element stability matrices. A = S L D L^T S with
    ! unit lower block-bidiagonal L, diagonal D of prescribed signs and a
    ! diagonal scaling S, so Sylvester's law gives the exact inertia. One
    ! process times one MODE: a FortSparse backend, or `lapack`, a dense block
    ! elimination with DSYTRF/DSYTRS/DGEMM on the same blocks.
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: iso_c_binding, only: c_double
    use fortsparse
    !$ use omp_lib, only: omp_get_max_threads
    implicit none
    interface
        real(c_double) function process_cpu() bind(C, name='fsbench_parent_cpu_seconds')
            import :: c_double
        end function process_cpu
        subroutine process_metrics(rss_kib, child_rss_kib, child_cpu_s) &
                bind(C, name='fsbench_process_metrics')
            import :: c_double
            real(c_double), intent(out) :: rss_kib, child_rss_kib, child_cpu_s
        end subroutine process_metrics
    end interface
    integer :: nb, w, repeats, n, nnz, negatives, threads, r, info
    real(dp), allocatable :: diag(:,:,:), lower(:,:,:), b(:), x(:)
    type(csc_t) :: a
    real(dp) :: decades, wall, cpu, residual
    ! Wall and process CPU seconds: the first factorization of the process,
    ! then the fastest of REPEATS refactorizations and solves.
    real(dp) :: t_factor, c_factor, t_refactor, c_refactor, t_numeric, t_solve
    real(dp) :: t_analysis
    real(c_double) :: rss_kib, child_rss_kib, child_cpu_s
    character(len=32) :: argument, label
    type(sparse_solver_t) :: solver
    type(fortsparse_status_t) :: status
    type(sparse_factor_stats_t) :: stats
    integer(i8) :: nneg, nzero, npos
    integer :: reliability, block_negative, active, ordering
    logical :: reused, reference

    if (command_argument_count() /= 5) error stop &
        'usage: bench_block_tridiagonal BLOCKS WIDTH DECADES MODE REPEATS'
    call get_command_argument(1, argument); read (argument, *) nb
    call get_command_argument(2, argument); read (argument, *) w
    call get_command_argument(3, argument); read (argument, *) decades
    call get_command_argument(4, label)
    call get_command_argument(5, argument); read (argument, *) repeats
    n = nb*w
    call build(nb, w, decades, diag, lower, negatives)
    threads = 1
    !$ threads = omp_get_max_threads()

    solver%backend_id = FORTSPARSE_BACKEND_NATIVE
    solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
    reference = .false.
    select case (trim(label))
    case ('lapack')
        reference = .true.
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
    case ('mumps')
        solver%backend_id = FORTSPARSE_BACKEND_MUMPS
    case ('superlu')
        solver%backend_id = FORTSPARSE_BACKEND_SUPERLU
    case ('umfpack')
        solver%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC
    case default
        error stop 'MODE: lapack|scalar|blocked|supernodal|supernodal-natural|' &
            //'supernodal-nd|mumps|superlu|umfpack'
    end select

    t_refactor = huge(1.0_dp)
    c_refactor = -1.0_dp
    t_numeric = -1.0_dp
    t_analysis = -1.0_dp
    t_solve = -1.0_dp
    nneg = -1
    reliability = -1
    active = -1
    ordering = -1
    residual = -1.0_dp
    nnz = w*w*(3*nb - 2)
    if (reference) then
        call timed_block(t_factor, c_factor)
        do r = 1, repeats
            call timed_block(wall, cpu)
            if (wall < t_refactor) then
                t_refactor = wall
                c_refactor = cpu
            end if
        end do
        nneg = block_negative
    else
        call to_csc(nb, w, diag, lower, a)
        allocate (b(n), x(n))
        call random_number(b)
        call timed_factor(.false., t_factor, c_factor)
        if (.not. status_ok(status)) then
            print '(a)', trim(status%msg)
            error stop 'factor failed'
        end if
        call sparse_factor_stats(solver, stats, status)
        t_analysis = stats%analysis_s
        t_solve = huge(1.0_dp)
        t_numeric = huge(1.0_dp)
        do r = 1, repeats
            call timed_factor(.true., wall, cpu)
            if (.not. status_ok(status)) error stop 'refactor failed'
            call sparse_factor_stats(solver, stats, status)
            if (wall < t_refactor) then
                t_refactor = wall
                c_refactor = cpu
            end if
            t_numeric = min(t_numeric, stats%factor_s)
            t_solve = min(t_solve, timed_solve())
        end do
        active = stats%active_threads
        ordering = stats%ordering_used
        residual = backward_error()
        ! Backends without inertia keep the fields at -1.
        call sparse_inertia(solver, nneg, nzero, npos, reliability, status)
        if (.not. status_ok(status)) then
            nneg = -1
            reliability = -1
        end if
    end if
    call process_metrics(rss_kib, child_rss_kib, child_cpu_s)
    print '(a)', 'mode,threads,active_threads,blocks,width,n,nnz,decades,repeats,' &
        //'factor_s,factor_cpu_s,refactor_s,refactor_cpu_s,analysis_s,numeric_s,' &
        //'solve_s,negative,expected,reliability,ordering,backward_error,' &
        //'process_cpu_s,peak_rss_kib,child_cpu_s,child_peak_rss_kib'
    print '(a,",",6(i0,","),f0.1,",",i0,",",7(es11.4,","),4(i0,","),' &
        //'es10.3,",",es11.4,",",f0.0,",",es11.4,",",f0.0)', &
        trim(label), threads, active, nb, w, n, nnz, decades, repeats, &
        t_factor, c_factor, t_refactor, c_refactor, t_analysis, t_numeric, &
        t_solve, nneg, negatives, reliability, ordering, residual, &
        process_cpu(), rss_kib, child_cpu_s, child_rss_kib
    if (nneg >= 0 .and. nneg /= negatives) error stop 'inertia mismatch'
    ! External backends keep their row: the backward error column tells
    ! whether the comparison is at matched accuracy.
    if (solver%backend_id == FORTSPARSE_BACKEND_NATIVE .and. .not. reference &
        .and. .not. residual < 1.0e-9_dp) error stop 'backward error'

contains

    real(dp) function seconds()
        integer(int64) :: count, rate
        call system_clock(count, rate)
        seconds = real(count, dp)/real(rate, dp)
    end function seconds

    subroutine timed_factor(again, elapsed, cpu_used)
        logical, intent(in) :: again
        real(dp), intent(out) :: elapsed, cpu_used
        real(dp) :: start, cpu_start
        cpu_start = process_cpu()
        start = seconds()
        if (again) then
            call sparse_refactor(solver, a, status, reused)
        else
            call sparse_factor(solver, a, status)
        end if
        elapsed = seconds() - start
        cpu_used = process_cpu() - cpu_start
    end subroutine timed_factor

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
    subroutine timed_block(elapsed, cpu_used)
        real(dp), intent(out) :: elapsed, cpu_used
        real(dp), allocatable :: schur(:,:,:), coupled(:,:), work(:)
        integer, allocatable :: pivots(:,:)
        real(dp) :: start, cpu_start
        integer :: i, j
        external :: dsytrf, dsytrs, dgemm
        allocate (schur(w, w, nb), coupled(w, w), work(64*w), pivots(w, nb))
        cpu_start = process_cpu()
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
        cpu_used = process_cpu() - cpu_start
    end subroutine timed_block

    subroutine build(nb, w, decades, diag, lower, negatives)
        integer, intent(in) :: nb, w
        real(dp), intent(in) :: decades
        real(dp), allocatable, intent(out) :: diag(:,:,:), lower(:,:,:)
        integer, intent(out) :: negatives
        real(dp), allocatable :: l(:,:,:), c(:,:,:), d(:,:), s(:,:), t(:,:)
        integer, allocatable :: seed(:)
        integer :: i, j, k, m
        external :: dgemm
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
            call dgemm('N', 'T', w, w, w, 1.0_dp, t, w, l(1, 1, i), w, 0.0_dp, &
                diag(1, 1, i), w)
            if (i > 1) then
                do j = 1, w
                    t(:, j) = c(:, j, i - 1)*d(j, i - 1)
                end do
                call dgemm('N', 'T', w, w, w, 1.0_dp, t, w, c(1, 1, i - 1), w, 1.0_dp, &
                    diag(1, 1, i), w)
                call dgemm('N', 'T', w, w, w, 1.0_dp, t, w, l(1, 1, i - 1), w, 0.0_dp, &
                    lower(1, 1, i - 1), w)
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
