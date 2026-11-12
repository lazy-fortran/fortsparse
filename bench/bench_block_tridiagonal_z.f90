program bench_block_tridiagonal_z
    ! Complex Hermitian indefinite block-tridiagonal pencils with dense blocks,
    ! the shape of Fourier-coupled stability matrices. Usage: BLOCKS WIDTH
    ! MODE with MODE 0 the scalar LDL^H and 1 the supernodal real embedding;
    ! prints wall and CPU seconds of the factorization and the inertia.
    use fortsparse
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    integer :: nb, w, n, i, j, k, b, mode, rel
    integer, allocatable :: rows(:), cols(:)
    complex(dp), allocatable :: vals(:), blk(:, :, :), off(:, :, :)
    real(dp) :: t0, t1
    integer(int64) :: c0, c1, rate, neg, zer, pos
    type(csc_z_t) :: a
    type(sparse_solver_t) :: s
    type(fortsparse_status_t) :: st
    character(16) :: arg
    call get_command_argument(1, arg); read(arg, *) nb
    call get_command_argument(2, arg); read(arg, *) w
    call get_command_argument(3, arg); read(arg, *) mode
    n = nb*w
    ! A = diagonally shifted Hermitian blocks with off-diagonal coupling;
    ! indefinite through negative shifts on some blocks.
    allocate(blk(w, w, nb), off(w, w, nb))
    do b = 1, nb
        do j = 1, w
            do i = 1, w
                blk(i, j, b) = cmplx(sin(0.1_dp*i*j + b), cos(0.3_dp*i - 0.2_dp*j + b), dp)/w
                off(i, j, b) = cmplx(cos(0.7_dp*i + 0.1_dp*j*b), sin(0.4_dp*i*j), dp)/w
            end do
        end do
        blk(:, :, b) = (blk(:, :, b) + conjg(transpose(blk(:, :, b))))/2
        do i = 1, w
            blk(i, i, b) = blk(i, i, b) + merge(-3.0_dp, 3.0_dp, mod(i + b, 11) == 0)
        end do
    end do
    allocate(rows(w*w*(3*nb)), cols(w*w*(3*nb)), vals(w*w*(3*nb)))
    k = 0
    do b = 1, nb
        do j = 1, w
            do i = 1, w
                k = k + 1; rows(k) = (b-1)*w + i; cols(k) = (b-1)*w + j; vals(k) = blk(i, j, b)
                if (b < nb) then
                    k = k + 1; rows(k) = b*w + i; cols(k) = (b-1)*w + j; vals(k) = off(i, j, b)
                    k = k + 1; rows(k) = (b-1)*w + j; cols(k) = b*w + i; vals(k) = conjg(off(i, j, b))
                end if
            end do
        end do
    end do
    call csc_from_triplet(n, n, rows(:k), cols(:k), vals(:k), a, st)
    s%backend_id = FORTSPARSE_BACKEND_NATIVE
    s%options%symmetry = FORTSPARSE_SYM_HERMITIAN
    s%options%ldlt_supernodal = mode == 1
    call system_clock(c0, rate); call cpu_time(t0)
    call sparse_factor(s, a, st)
    call system_clock(c1); call cpu_time(t1)
    if (.not. status_ok(st)) error stop trim(st%msg)
    call sparse_inertia(s, neg, zer, pos, rel, st)
    print '(a,i2,a,i7,a,i9,a,f9.3,a,f9.2,a,3i7)', 'mode', mode, ' n', n, ' nnz', a%nnz, &
        ' wall', real(c1 - c0)/rate, ' cpu', t1 - t0, ' inertia', neg, zer, pos
end program bench_block_tridiagonal_z
