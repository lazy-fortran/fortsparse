program test_fortsparse_superlu_status
    ! SuperLU factor status truthfulness. The mapping of shim codes is checked
    ! with injected values; a real allocation failure is forced by factoring a
    ! 2D Laplacian (large COLAMD fill) under an address-space limit (the
    ! FORTSPARSE_TEST_ULIMIT_V environment is set by the ctest oom case; if the factor
    ! unexpectedly succeeds the limit was not effective and the case is
    ! reported as skipped). Allocation failure must never be SINGULAR.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use, intrinsic :: iso_c_binding, only: c_int
    use fortsparse, only: dp, csc_t, csc_from_triplet, sparse_solver_t, &
        sparse_factor, sparse_free, fortsparse_status_t, &
        FORTSPARSE_OK, FORTSPARSE_SINGULAR, FORTSPARSE_OUT_OF_MEMORY, &
        FORTSPARSE_INDEX_OVERFLOW, FORTSPARSE_INTERNAL_ERROR, &
        FORTSPARSE_BACKEND_SUPERLU, FORTSPARSE_BACKEND_UNAVAILABLE
#ifdef FORTSPARSE_HAVE_SUPERLU
    use fortsparse_superlu, only: slu_map_factor_status
#endif
    implicit none

    integer :: nfail
    block
        type(csc_t) :: probe_matrix
        type(sparse_solver_t) :: probe_solver
        type(fortsparse_status_t) :: probe_status
        call csc_from_triplet(1, 1, [1], [1], [1.0_dp], &
            probe_matrix, probe_status)
        if (probe_status%code /= 0) error stop "backend probe matrix failed"
        probe_solver%backend_id = FORTSPARSE_BACKEND_SUPERLU
        call sparse_factor(probe_solver, probe_matrix, probe_status)
        call sparse_free(probe_solver)
        if (probe_status%code == FORTSPARSE_BACKEND_UNAVAILABLE) then
            write (*, '(a)') 'SKIP: SuperLU backend is not configured'
            stop 0
        end if
        if (probe_status%code /= 0) error stop "configured SuperLU probe failed"
    end block
    nfail = 0
#ifdef FORTSPARSE_HAVE_SUPERLU
    call test_mapping(nfail)
#endif
    call test_real_allocation_failure(nfail)
    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"

contains

#ifdef FORTSPARSE_HAVE_SUPERLU
    subroutine expect(label, info, code, nfail)
        character(*), intent(in)    :: label
        integer(c_int), intent(in)  :: info
        integer, intent(in)         :: code
        integer, intent(inout)      :: nfail
        type(fortsparse_status_t) :: st
        call slu_map_factor_status(info, st)
        if (st%code /= code) then
            nfail = nfail + 1
            write (error_unit, "(a,a,i0,a,i0)") label, ": got ", st%code, &
                " expected ", code
        end if
        if (code /= FORTSPARSE_OK .and. len_trim(st%msg) == 0) then
            nfail = nfail + 1
            write (error_unit, "(a,a)") label, ": empty message"
        end if
    end subroutine expect

    subroutine test_mapping(nfail)
        integer, intent(inout) :: nfail
        call expect("ok", 0_c_int, FORTSPARSE_OK, nfail)
        call expect("singular", 1_c_int, FORTSPARSE_SINGULAR, nfail)
        call expect("alloc", 2_c_int, FORTSPARSE_OUT_OF_MEMORY, nfail)
        call expect("overflow", 3_c_int, FORTSPARSE_INDEX_OVERFLOW, nfail)
        call expect("internal", -1_c_int, FORTSPARSE_INTERNAL_ERROR, nfail)
        call expect("unknown_pos", 77_c_int, FORTSPARSE_INTERNAL_ERROR, nfail)
    end subroutine test_mapping

#endif

    subroutine test_real_allocation_failure(nfail)
        integer, intent(inout) :: nfail
        integer, parameter :: k = 400, n = k*k
        type(csc_t) :: A
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: status
        integer, allocatable :: r(:), c(:)
        real(dp), allocatable :: v(:)
        integer :: i, j, m, id
        character(len=16) :: env
        integer :: elen

        call get_environment_variable("FORTSPARSE_TEST_ULIMIT_V", env, elen)
        if (elen == 0) then
            write (*, "(a)") "SKIP real allocation failure: no address limit"
            return
        end if
        allocate (r(5*n), c(5*n), v(5*n))
        m = 0
        do j = 1, k
            do i = 1, k
                id = (j - 1)*k + i
                m = m + 1; r(m) = id; c(m) = id; v(m) = 4.0_dp
                if (i > 1) then
                    m = m + 1; r(m) = id; c(m) = id - 1; v(m) = -1.0_dp
                end if
                if (i < k) then
                    m = m + 1; r(m) = id; c(m) = id + 1; v(m) = -1.0_dp
                end if
                if (j > 1) then
                    m = m + 1; r(m) = id; c(m) = id - k; v(m) = -1.0_dp
                end if
                if (j < k) then
                    m = m + 1; r(m) = id; c(m) = id + k; v(m) = -1.0_dp
                end if
            end do
        end do
        call csc_from_triplet(n, n, r(1:m), c(1:m), v(1:m), A, status)
        solver%backend_id = FORTSPARSE_BACKEND_SUPERLU
        call sparse_factor(solver, A, status)
        if (status%code == FORTSPARSE_OK) then
            write (*, "(a)") "SKIP real allocation failure: limit not effective"
            call sparse_free(solver)
        else if (status%code /= FORTSPARSE_OUT_OF_MEMORY) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,2a)") "real_alloc: code ", status%code, &
                " msg ", trim(status%msg)
        end if
    end subroutine test_real_allocation_failure

end program test_fortsparse_superlu_status
