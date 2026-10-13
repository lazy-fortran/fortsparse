program bench_fortsparse_mumps
    ! Wall-clock factor and solve timings for the MUMPS backend.
    !   bench_fortsparse_mumps <p3|lap> <size> <threads> [nrhs] [both]
    ! FS_BENCH_BACKEND=<any> selects SuperLU for a baseline; FS_BENCH_VERBOSE=<any>
    ! turns on MUMPS output.
    ! With the optional fifth argument `both` the run is repeated without
    ! iterative refinement first.
    ! p3: size = cells per side (n ~ 9 size^2); lap: size = grid side.
    use, intrinsic :: iso_fortran_env, only: int64
    use fortsparse, only: dp, i8, csc_t, csc_matvec, sparse_solver_t, &
        sparse_factor, sparse_solve, sparse_solve_block, sparse_destroy, &
        sparse_factor_stats, sparse_factor_stats_t, fortsparse_status_t, &
        status_ok, FORTSPARSE_BACKEND_MUMPS, FORTSPARSE_BACKEND_SUPERLU, FORTSPARSE_SYM_SPD
    use fs_test_matrices, only: laplace2d, p3_poisson, inf_residual
    use fs_test_mpi, only: test_mpi_init, test_mpi_finalize
    implicit none

    character(16) :: kind, arg
    integer :: size_, nthreads, nrhs, i, k, nrep, first
    type(csc_t) :: A
    type(sparse_solver_t) :: solver
    type(sparse_factor_stats_t) :: st
    type(fortsparse_status_t) :: status
    real(dp), allocatable :: b(:), x(:), B2(:, :), X2(:, :)
    real(dp) :: t0, t_asm, t_fac, t_solve, t_solve_raw, t_block, res
    logical :: refine
    integer :: nv, j

    call test_mpi_init()
    call get_command_argument(1, kind)
    call get_command_argument(2, arg); read (arg, *) size_
    call get_command_argument(3, arg); read (arg, *) nthreads
    nrhs = 16
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg); read (arg, *) nrhs
    end if
    first = 2
    if (command_argument_count() >= 5) then
        call get_command_argument(5, arg)
        if (trim(arg) == "both") first = 1
    end if
    nrep = 3

    t0 = wtime()
    if (trim(kind) == "p3") then
        call p3_poisson(size_, A)
    else
        call laplace2d(size_, A)
    end if
    t_asm = wtime() - t0

    allocate (b(A%ncol), x(A%ncol), B2(A%ncol, nrhs), X2(A%ncol, nrhs))
    do i = 1, A%ncol
        b(i) = sin(1.234_dp*i) + 0.5_dp*cos(0.37_dp*i)
    end do
    do k = 1, nrhs
        do i = 1, A%ncol
            B2(i, k) = sin(1.234_dp*i + 0.7_dp*k) + 0.5_dp*cos(0.37_dp*i*k)
        end do
    end do

    do i = first, 2
        refine = (i == 2)
        solver%backend_id = FORTSPARSE_BACKEND_MUMPS
        call get_environment_variable("FS_BENCH_BACKEND", length=nv)
        if (nv > 0) solver%backend_id = FORTSPARSE_BACKEND_SUPERLU
        solver%refine = refine
        solver%options%symmetry = FORTSPARSE_SYM_SPD
        solver%options%nthreads = nthreads
        call get_environment_variable("FS_BENCH_VERBOSE", length=nv)
        solver%options%verbose = (nv > 0)
        t0 = wtime()
        call sparse_factor(solver, A, status)
        t_fac = wtime() - t0
        if (.not. status_ok(status)) error stop trim(status%msg)
        call sparse_factor_stats(solver, st, status) ! unavailable: SuperLU
        t0 = wtime()
        do k = 1, nrep
            call sparse_solve(solver, b, x, status)
        end do
        t_solve = (wtime() - t0)/nrep
        res = inf_residual(A, x, b)
        t0 = wtime()
        call sparse_solve_block(solver, B2, X2, status)
        t_block = (wtime() - t0)
        t_solve_raw = t_solve
        write (*, "(a,a,a,i0,a,i0,a,l1)") "kind=", trim(kind), " n=", A%ncol, &
            " nnz=", A%nnz, " refine=", refine
        write (*, "(a,i0,a,f9.3,a,f9.3,a,f9.3,a,f9.3)") " threads=", nthreads, &
            " assemble_s=", t_asm, " factor_s=", t_fac, " solve1_s=", t_solve, &
            " block_s(" // trim(itoa(nrhs)) // ")=", t_block
        write (*, "(a,es10.2,a,i0,a,f8.2,a,f9.1,a)") " residual=", res, &
            " factor_nnz=", st%factor_nnz, " fill=", &
            real(st%factor_nnz, dp)/real(A%nnz, dp), " flops(G)=", st%flops/1.0e9_dp, ""
        write (*, "(a,es25.17,a,es25.17)") " norm2(x)=", norm2(x), &
            " checksum=", sum(x*[(1.0_dp + mod(j, 7), j=1, size(x))])
        write (*, "(a,f9.3,a,f9.3)") " analysis_s=", st%analysis_s, &
            " numeric_factor_s=", st%factor_s
        write (*, "(a,f9.1,a,f9.1,a)") " factor_mem_MB=", st%factor_mem_mb, &
            " factor_storage_MB=", real(st%factor_nnz, dp)*8.0_dp/1.0e6_dp, ""
        call sparse_destroy(solver)
    end do

    call test_mpi_finalize()

contains

    function itoa(n) result(s)
        integer, intent(in) :: n
        character(16) :: s
        write (s, "(i0)") n
    end function itoa

    function wtime() result(t)
        real(dp) :: t
        integer(int64) :: c, r
        call system_clock(c, r)
        t = real(c, dp)/real(r, dp)
    end function wtime

end program bench_fortsparse_mumps
