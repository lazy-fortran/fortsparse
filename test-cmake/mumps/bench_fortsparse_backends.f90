program bench_fortsparse_backends
    ! Factor / refactor / solve wall times of one backend on one FE-shaped
    ! matrix (see fs_bench_matrices for the kinds):
    !   bench_fortsparse_backends <kind> <size> <threads> <superlu|mumps|auto>
    !                             [default|amd|nd]
    ! threads sets options%nthreads (MUMPS) and is reported; run with the same
    ! OMP_NUM_THREADS. One output line per run:
    !   kind n nnz backend used ord thr once_us factor_s refactor_s solve_s resid
    !   flops analysis_s
    ! once_us is the mean wall time of factor + one solve + destroy, the
    ! per-call cost of sparse_solve_once; factor_s includes the analysis,
    ! refactor_s keeps the pattern (MUMPS reuses its analysis).
    use, intrinsic :: iso_fortran_env, only: int64
    use fortsparse, only: dp, csc_t, csc_z_t, sparse_solver_t, &
        sparse_factor, sparse_refactor, sparse_solve, sparse_destroy, &
        sparse_factor_stats, sparse_factor_stats_t, &
        fortsparse_status_t, status_ok, FORTSPARSE_BACKEND_MUMPS, &
        FORTSPARSE_BACKEND_SUPERLU, FORTSPARSE_BACKEND_AUTO, &
        FORTSPARSE_SYM_SPD, FORTSPARSE_SYM_GENERAL, FORTSPARSE_SYM_INDEFINITE, &
        FORTSPARSE_ORDER_DEFAULT, FORTSPARSE_ORDER_AMD, FORTSPARSE_ORDER_ND
    use fortsparse_mumps, only: mumps_backend_t
    use fs_bench_matrices, only: bench_matrix_real, bench_matrix_complex, &
        bench_kind_is_complex
    implicit none

    character(16) :: kind, arg, bname, oname
    integer :: size_, nthreads, bid, ord, sym, nrep, k
    type(sparse_solver_t) :: solver
    type(fortsparse_status_t) :: status
    type(sparse_factor_stats_t) :: st
    type(csc_t) :: A, A2
    type(csc_z_t) :: Z, Z2
    real(dp), allocatable :: b(:), x(:)
    complex(dp), allocatable :: bz(:), xz(:)
    real(dp) :: t0, t_once, t_fac, t_refac, t_solve, res
    integer :: n, nnz
    logical :: cplx, reused

    call get_command_argument(1, kind)
    call get_command_argument(2, arg); read (arg, *) size_
    call get_command_argument(3, arg); read (arg, *) nthreads
    call get_command_argument(4, bname)
    oname = "default"
    if (command_argument_count() >= 5) call get_command_argument(5, oname)
    select case (trim(bname))
    case ("superlu"); bid = FORTSPARSE_BACKEND_SUPERLU
    case ("mumps"); bid = FORTSPARSE_BACKEND_MUMPS
    case ("auto"); bid = FORTSPARSE_BACKEND_AUTO
    case default; error stop "backend: superlu|mumps|auto"
    end select
    select case (trim(oname))
    case ("amd"); ord = FORTSPARSE_ORDER_AMD
    case ("nd"); ord = FORTSPARSE_ORDER_ND
    case default; ord = FORTSPARSE_ORDER_DEFAULT
    end select
    select case (trim(kind))
    case ("p1", "p2", "p3", "tet1"); sym = FORTSPARSE_SYM_SPD
    case ("ideal"); sym = FORTSPARSE_SYM_INDEFINITE
    case default; sym = FORTSPARSE_SYM_GENERAL
    end select

    cplx = bench_kind_is_complex(kind)
    if (cplx) then
        call bench_matrix_complex(kind, size_, Z)
        call bench_matrix_complex(kind, size_, Z2, shift=-0.01_dp)
        n = Z%ncol; nnz = Z%nnz
        allocate (bz(n), xz(n))
        bz = cmplx([(sin(1.3_dp*k), k=1, n)], [(cos(0.7_dp*k), k=1, n)], dp)
    else
        call bench_matrix_real(kind, size_, A)
        call bench_matrix_real(kind, size_, A2, shift=-0.01_dp)
        n = A%ncol; nnz = A%nnz
        allocate (b(n), x(n))
        b = [(sin(1.3_dp*k), k=1, n)]
    end if
    nrep = max(1, min(2000, int(2.0e5_dp/real(nnz, dp))))

    ! Per-call cost: factor + one solve + destroy, averaged.
    t0 = wtime()
    do k = 1, nrep
        call configure()
        call factor_once()
        call solve_once()
        call sparse_destroy(solver)
    end do
    t_once = (wtime() - t0)/nrep

    ! Separate phases on a retained solver.
    call configure()
    t0 = wtime()
    call factor_once()
    t_fac = wtime() - t0
    call sparse_factor_stats(solver, st, status) ! -1 when not reported
    t0 = wtime()
    if (cplx) then
        call sparse_refactor(solver, Z2, status, reused)
    else
        call sparse_refactor(solver, A2, status, reused)
    end if
    t_refac = wtime() - t0
    if (.not. status_ok(status)) error stop trim(status%msg)
    t0 = wtime()
    do k = 1, nrep
        call solve_once()
    end do
    t_solve = (wtime() - t0)/nrep
    if (cplx) then
        res = maxval(abs(matvec_z(Z2, xz) - bz))/maxval(abs(bz))
    else
        res = maxval(abs(matvec_r(A2, x) - b))/maxval(abs(b))
    end if
    write (*, "(a,1x,i0,1x,i0,1x,a,1x,a,1x,a,1x,i0,1x,f12.1,3(1x,es10.3)," &
        //"1x,es9.2,2(1x,es10.3))") &
        trim(kind), n, nnz, trim(bname), used_name(), trim(oname), nthreads, &
        1.0e6_dp*t_once, t_fac, t_refac, t_solve, res, st%flops, st%analysis_s
    call sparse_destroy(solver)

contains

    subroutine configure()
        solver%backend_id = bid
        solver%options%symmetry = sym
        solver%options%ordering = ord
        solver%options%nthreads = nthreads
    end subroutine configure

    subroutine factor_once()
        if (cplx) then
            call sparse_factor(solver, Z, status)
        else
            call sparse_factor(solver, A, status)
        end if
        if (.not. status_ok(status)) error stop trim(status%msg)
    end subroutine factor_once

    subroutine solve_once()
        if (cplx) then
            call sparse_solve(solver, bz, xz, status)
        else
            call sparse_solve(solver, b, x, status)
        end if
        if (.not. status_ok(status)) error stop trim(status%msg)
    end subroutine solve_once

    function used_name() result(s)
        character(8) :: s

        s = "superlu"
        if (.not. allocated(solver%backend)) then
            s = "none"
            return
        end if
        select type (be => solver%backend)
        type is (mumps_backend_t)
            s = "mumps"
        end select
    end function used_name

    function matvec_r(M, v) result(y)
        type(csc_t), intent(in) :: M
        real(dp), intent(in) :: v(:)
        real(dp) :: y(M%nrow)
        integer :: j, p

        y = 0.0_dp
        do j = 1, M%ncol
            do p = M%col_ptr(j), M%col_ptr(j + 1) - 1
                y(M%row_idx(p)) = y(M%row_idx(p)) + M%val(p)*v(j)
            end do
        end do
    end function matvec_r

    function matvec_z(M, v) result(y)
        type(csc_z_t), intent(in) :: M
        complex(dp), intent(in) :: v(:)
        complex(dp) :: y(M%nrow)
        integer :: j, p

        y = (0.0_dp, 0.0_dp)
        do j = 1, M%ncol
            do p = M%col_ptr(j), M%col_ptr(j + 1) - 1
                y(M%row_idx(p)) = y(M%row_idx(p)) + M%val(p)*v(j)
            end do
        end do
    end function matvec_z

    function wtime() result(t)
        real(dp) :: t
        integer(int64) :: c, r
        call system_clock(c, r)
        t = real(c, dp)/real(r, dp)
    end function wtime

end program bench_fortsparse_backends
