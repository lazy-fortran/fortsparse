program test_fortsparse_umfpack_ipc
    ! Round-trip the out-of-process UMFPACK backend through the GPL helper. The
    ! ctest harness sets FORTSPARSE_UMFPACK_HELPER to the built helper. Selects
    ! FORTSPARSE_BACKEND_UMFPACK_IPC and solves a known real and complex system
    ! over shared memory, checking the solution and residual.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortsparse, only: dp, csc_t, csc_z_t, csc_from_triplet, csc_matvec, &
        sparse_solver_t, sparse_factor, sparse_solve, sparse_free, &
        sparse_vector, fortsparse_status_t, status_ok, &
        FORTSPARSE_BACKEND_UMFPACK_IPC
    implicit none

    integer :: nfail
    nfail = 0

    if (.not. helper_present()) then
        ! No GPL helper in this build (the MIT-only fpm/fo path). The
        ! out-of-process backend is unavailable by design, so there is nothing
        ! to round-trip. ctest builds the helper and sets the variable.
        write (*, "(a)") "SKIP: FORTSPARSE_UMFPACK_HELPER not set"
        write (*, "(a)") "PASS"
        stop 0
    end if

    call run_real(nfail)
    call run_complex(nfail)
    call run_concurrent(nfail)
    call run_reuse(nfail)
    call run_inplace(nfail)
    call run_zerocopy(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    ! True when the ctest harness points at the built GPL helper binary.
    logical function helper_present()
        character(4096) :: buf
        integer         :: n

        call get_environment_variable("FORTSPARSE_UMFPACK_HELPER", buf, n)
        helper_present = (n > 0)
    end function helper_present

    ! Real system A x = b solved through the helper, second RHS reuses factors.
    subroutine run_real(nfail)
        integer, intent(inout) :: nfail

        type(csc_t)               :: A
        type(sparse_solver_t)     :: solver
        type(fortsparse_status_t) :: status
        integer                   :: rows(7), cols(7)
        real(dp)                  :: vals(7)
        real(dp)                  :: b(3), x(3), x2(3), xe(3), xe2(3)

        rows = [1, 1, 2, 2, 2, 3, 3]
        cols = [1, 2, 1, 2, 3, 2, 3]
        vals = [4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]
        call csc_from_triplet(3, 3, rows, cols, vals, A, status)
        call check_true("real_build", status_ok(status), nfail)

        solver%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC
        xe = [1.0_dp, 2.0_dp, 3.0_dp]
        b = [6.0_dp, 10.0_dp, 8.0_dp]
        call sparse_factor(solver, A, status)
        call check_true("real_factor", status_ok(status), nfail)
        ! The shm/semaphore names must already be unlinked once the helper has
        ! attached, even though this session is still live: a hard kill now
        ! leaves no /dev/shm residue. Names gone == leak-free under SIGKILL.
        call check_no_shm_residue("real_active_session", nfail)

        call sparse_solve(solver, b, x, status)
        call check_true("real_solve", status_ok(status), nfail)
        call check_err("real_x", x, xe, nfail)

        xe2 = [0.0_dp, 1.0_dp, 0.0_dp]
        call sparse_solve(solver, csc_matvec(A, xe2), x2, status)
        call check_true("real_solve2", status_ok(status), nfail)
        call check_err("real_x2", x2, xe2, nfail)
        call sparse_free(solver)
        call check_no_shm_residue("real_after_free", nfail)
    end subroutine run_real

    ! Complex system A x = b solved through the helper.
    subroutine run_complex(nfail)
        integer, intent(inout) :: nfail

        type(csc_z_t)             :: A
        type(sparse_solver_t)     :: solver
        type(fortsparse_status_t) :: status
        integer                   :: rows(5), cols(5)
        complex(dp)               :: vals(5), b(3), x(3), xe(3)

        rows = [1, 2, 2, 3, 3]
        cols = [1, 1, 2, 2, 3]
        vals = [cmplx(2.0_dp, 1.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
            cmplx(3.0_dp, -1.0_dp, dp), cmplx(1.0_dp, 1.0_dp, dp), &
            cmplx(2.0_dp, 0.0_dp, dp)]
        call csc_from_triplet(3, 3, rows, cols, vals, A, status)
        call check_true("cplx_build", status_ok(status), nfail)

        solver%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC
        xe = [cmplx(1.0_dp, -1.0_dp, dp), cmplx(0.0_dp, 2.0_dp, dp), &
            cmplx(3.0_dp, 0.0_dp, dp)]
        b = csc_matvec(A, xe)
        call sparse_factor(solver, A, status)
        call check_true("cplx_factor", status_ok(status), nfail)

        call sparse_solve(solver, b, x, status)
        call check_true("cplx_solve", status_ok(status), nfail)
        call check_zerr("cplx_x", x, xe, nfail)
        call sparse_free(solver)
    end subroutine run_complex

    ! Two factorizations held live at once, a real and a complex solver, must
    ! not share resident factors. Both are factored before either is solved, so
    ! a single shared helper session would let the complex factor overwrite the
    ! real one and corrupt the real solve. The persistent-session pool hands each
    ! live factorization its own session; solving both, then the real one again
    ! after the complex solve, proves the isolation.
    subroutine run_concurrent(nfail)
        integer, intent(inout) :: nfail

        type(csc_t)               :: ar
        type(csc_z_t)             :: az
        type(sparse_solver_t)     :: sr, sz
        type(fortsparse_status_t) :: status
        integer                   :: rrows(7), rcols(7), zrows(5), zcols(5)
        real(dp)                  :: rvals(7), rb(3), rx(3), rxe(3)
        complex(dp)               :: zvals(5), zb(3), zx(3), zxe(3)

        rrows = [1, 1, 2, 2, 2, 3, 3]
        rcols = [1, 2, 1, 2, 3, 2, 3]
        rvals = [4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]
        call csc_from_triplet(3, 3, rrows, rcols, rvals, ar, status)
        rxe = [1.0_dp, 2.0_dp, 3.0_dp]
        rb = csc_matvec(ar, rxe)

        zrows = [1, 2, 2, 3, 3]
        zcols = [1, 1, 2, 2, 3]
        zvals = [cmplx(2.0_dp, 1.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
            cmplx(3.0_dp, -1.0_dp, dp), cmplx(1.0_dp, 1.0_dp, dp), &
            cmplx(2.0_dp, 0.0_dp, dp)]
        call csc_from_triplet(3, 3, zrows, zcols, zvals, az, status)
        zxe = [cmplx(1.0_dp, -1.0_dp, dp), cmplx(0.0_dp, 2.0_dp, dp), &
            cmplx(3.0_dp, 0.0_dp, dp)]
        zb = csc_matvec(az, zxe)

        sr%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC
        sz%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC

        call sparse_factor(sr, ar, status)
        call check_true("conc_real_factor", status_ok(status), nfail)
        call sparse_factor(sz, az, status)
        call check_true("conc_cplx_factor", status_ok(status), nfail)

        call sparse_solve(sr, rb, rx, status)
        call check_true("conc_real_solve", status_ok(status), nfail)
        call check_err("conc_real_x", rx, rxe, nfail)
        call sparse_solve(sz, zb, zx, status)
        call check_true("conc_cplx_solve", status_ok(status), nfail)
        call check_zerr("conc_cplx_x", zx, zxe, nfail)
        call sparse_solve(sr, rb, rx, status)
        call check_err("conc_real_x_again", rx, rxe, nfail)

        call sparse_free(sr)
        call sparse_free(sz)
    end subroutine run_concurrent

    ! Many factor/solve/free cycles reuse pooled helpers without respawning one
    ! per factorization. Every cycle must still return the right answer: a stale
    ! reused mapping or leftover resident factors would corrupt a later solve.
    subroutine run_reuse(nfail)
        integer, intent(inout) :: nfail

        type(csc_t)               :: A
        type(sparse_solver_t)     :: solver
        type(fortsparse_status_t) :: status
        integer                   :: rows(7), cols(7), i
        real(dp)                  :: vals(7), b(3), x(3), xe(3)

        rows = [1, 1, 2, 2, 2, 3, 3]
        cols = [1, 2, 1, 2, 3, 2, 3]
        vals = [4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]
        call csc_from_triplet(3, 3, rows, cols, vals, A, status)
        solver%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC

        do i = 1, 50
            xe = [real(i, dp), 2.0_dp, 3.0_dp]
            b = csc_matvec(A, xe)
            call sparse_factor(solver, A, status)
            call check_true("reuse_factor", status_ok(status), nfail)
            call sparse_solve(solver, b, x, status)
            call check_err("reuse_x", x, xe, nfail)
            call sparse_free(solver)
        end do
    end subroutine run_reuse

    ! In-place real solve: b carries the RHS in and the solution out, with no
    ! caller temporary. Must match the two-vector solve.
    subroutine run_inplace(nfail)
        integer, intent(inout) :: nfail

        type(csc_t)               :: A
        type(sparse_solver_t)     :: solver
        type(fortsparse_status_t) :: status
        integer                   :: rows(7), cols(7)
        real(dp)                  :: vals(7), b(3), xe(3)

        rows = [1, 1, 2, 2, 2, 3, 3]
        cols = [1, 2, 1, 2, 3, 2, 3]
        vals = [4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]
        call csc_from_triplet(3, 3, rows, cols, vals, A, status)
        solver%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC
        xe = [1.0_dp, 2.0_dp, 3.0_dp]
        b = [6.0_dp, 10.0_dp, 8.0_dp]
        call sparse_factor(solver, A, status)
        call check_true("inplace_factor", status_ok(status), nfail)
        call sparse_solve(solver, b, status)
        call check_true("inplace_solve", status_ok(status), nfail)
        call check_err("inplace_x", b, xe, nfail)
        ! A second in-place solve reuses the factorization.
        b = csc_matvec(A, [0.0_dp, 1.0_dp, 0.0_dp])
        call sparse_solve(solver, b, status)
        call check_err("inplace_x2", b, [0.0_dp, 1.0_dp, 0.0_dp], nfail)
        call sparse_free(solver)
    end subroutine run_inplace

    ! Zero-copy solve: the RHS and the solution are fortsparse vectors, i.e.
    ! slots in the shared mapping, so the helper reads and writes them directly
    ! with no marshalling. Exercises reuse and a mixed pool/plain-array call; a
    ! wrong slot offset would corrupt the result, so correctness here certifies
    ! the offset and pool logic.
    subroutine run_zerocopy(nfail)
        integer, intent(inout) :: nfail

        type(csc_t)               :: A
        type(sparse_solver_t)     :: solver
        type(fortsparse_status_t) :: status
        integer                   :: rows(7), cols(7)
        real(dp)                  :: vals(7), xe(3), xreg(3)
        real(dp), pointer         :: b(:), x(:)

        rows = [1, 1, 2, 2, 2, 3, 3]
        cols = [1, 2, 1, 2, 3, 2, 3]
        vals = [4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]
        call csc_from_triplet(3, 3, rows, cols, vals, A, status)
        solver%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC
        call sparse_factor(solver, A, status)
        call check_true("zc_factor", status_ok(status), nfail)

        b => sparse_vector(solver, 3)
        x => sparse_vector(solver, 3)
        call check_true("zc_vec_b", associated(b), nfail)
        call check_true("zc_vec_x", associated(x), nfail)
        if (.not. (associated(b) .and. associated(x))) return

        xe = [1.0_dp, 2.0_dp, 3.0_dp]
        b = csc_matvec(A, xe)
        call sparse_solve(solver, b, x, status)
        call check_true("zc_solve", status_ok(status), nfail)
        call check_err("zc_x", x, xe, nfail)

        ! Reuse the factorization, still zero-copy.
        b = csc_matvec(A, [0.0_dp, 1.0_dp, 0.0_dp])
        call sparse_solve(solver, b, x, status)
        call check_err("zc_x2", x, [0.0_dp, 1.0_dp, 0.0_dp], nfail)

        ! Mixed: shared-memory RHS, plain-array solution (marshalled out).
        b = csc_matvec(A, xe)
        call sparse_solve(solver, b, xreg, status)
        call check_err("zc_mixed", xreg, xe, nfail)

        call sparse_free(solver)
    end subroutine run_zerocopy

    subroutine check_err(label, got, want, nfail)
        character(*), intent(in)    :: label
        real(dp),     intent(in)    :: got(:), want(:)
        integer,      intent(inout) :: nfail
        if (maxval(abs(got - want)) > 1.0e-9_dp) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,es13.6)") "FAIL [", label, &
                "] max|err| ", maxval(abs(got - want))
        end if
    end subroutine check_err

    subroutine check_zerr(label, got, want, nfail)
        character(*), intent(in)    :: label
        complex(dp),  intent(in)    :: got(:), want(:)
        integer,      intent(inout) :: nfail
        if (maxval(abs(got - want)) > 1.0e-9_dp) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,es13.6)") "FAIL [", label, &
                "] max|err| ", maxval(abs(got - want))
        end if
    end subroutine check_zerr

    subroutine check_true(label, cond, nfail)
        character(*), intent(in)    :: label
        logical,      intent(in)    :: cond
        integer,      intent(inout) :: nfail
        if (.not. cond) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a)") "FAIL [", label, "] expected .true."
        end if
    end subroutine check_true

    ! Fail if any /dev/shm/fsparse_* object exists. The early-unlink handshake
    ! drops the names as soon as the helper attaches, so a residue here means a
    ! hard kill would leak shared memory on the node. Skipped where /dev/shm is
    ! absent (non-Linux), since the check is platform-specific.
    subroutine check_no_shm_residue(label, nfail)
        character(*), intent(in)    :: label
        integer,      intent(inout) :: nfail
        character(*), parameter     :: out = "fsparse_shm_residue.txt"
        integer                     :: u, n, ios, cmdstat
        logical                     :: have_shm

        inquire (file="/dev/shm/.", exist=have_shm)
        if (.not. have_shm) return

        call execute_command_line( &
            "ls /dev/shm 2>/dev/null | grep -c fsparse > " // out, &
            exitstat=ios, cmdstat=cmdstat)
        if (cmdstat /= 0) return

        n = -1
        open (newunit=u, file=out, status="old", action="read", iostat=ios)
        if (ios == 0) then
            read (u, *, iostat=ios) n
            close (u, status="delete")
        end if
        if (n /= 0) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,i0,a)") "FAIL [", label, &
                "] /dev/shm fsparse residue count ", n, " (want 0)"
        end if
    end subroutine check_no_shm_residue

end program test_fortsparse_umfpack_ipc
