program test_fortsparse_native_hpd
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csc, only: csc_z_t, csc_from_triplet
    use fortsparse_status, only: fortsparse_status_t, FORTSPARSE_OK, &
                          FORTSPARSE_INVALID_MATRIX, FORTSPARSE_NOT_POSITIVE_DEFINITE, &
                                 FORTSPARSE_NOT_FACTORED, FORTSPARSE_RESOURCE_LIMIT, &
        FORTSPARSE_INTERNAL_ERROR
    use fortsparse_native_hpd, only: native_hpd_t
    implicit none
    type(native_hpd_t) :: solver, capped
    type(csc_z_t) :: A, bad
    type(fortsparse_status_t) :: status
    integer, parameter :: n = 7
    complex(dp) :: dense(n, n), truth(n), rhs(n), solution(n), u(n), v(n)
    complex(dp) :: blocks(n, 10), answers(n, 10), reference(n, 10), adjoint(n)
    real(dp) :: lhs, right
    integer :: i, j, k, nz, rows(n*n), cols(n*n)
    complex(dp) :: values(n*n)
    integer(i8) :: bytes, peak, fb, pb, wb
    logical :: reused

    dense = cmplx(0.0_dp, 0.0_dp, dp)
    do i = 1, n
        dense(i, i) = cmplx(6.0_dp + real(i, dp)/10, 0.0_dp, dp)
        truth(i) = cmplx(real(i, dp)/4, 0.3_dp - real(i, dp)/7, dp)
        u(i) = cmplx(sin(real(i, dp)), cos(real(i, dp)), dp)
        v(i) = cmplx(0.2_dp*real(i, dp), sin(real(i, dp)/2), dp)
    end do
    do i = 1, 4
        dense(i, i + 1) = cmplx(-0.3_dp, real(i, dp)/6, dp)
    end do
    dense(1, 5) = cmplx(0.7_dp, -0.5_dp, dp)
    dense(6, 7) = cmplx(-0.2_dp, 0.3_dp, dp)
    do j = 1, n
        do i = j + 1, n
            dense(i, j) = conjg(dense(j, i))
        end do
    end do
    nz = 0
    do j = 1, n
        do i = 1, n
            if (abs(dense(i, j)) == 0.0_dp) cycle
            nz = nz + 1
            rows(nz) = i
            cols(nz) = j
            values(nz) = dense(i, j)
        end do
    end do
    call csc_from_triplet(n, n, rows(:nz), cols(:nz), values(:nz), A, status)
    call require(status%code == FORTSPARSE_OK, 'CSC fixture')
    call solver%factor(A, status)
    call require(status%code == FORTSPARSE_OK, 'HPD factor')
    call require(solver%is_factored(), 'ready')
    call require(solver%dimension() == n, 'dimension')
    ! The cycle requires genuine fill; disconnected block is retained.
    call require(solver%factor_nnz() > int((nz + n)/2, i8), 'symbolic fill')
    rhs = matmul(dense, truth)
    call solver%solve(rhs, solution, status)
    call require(status%code == FORTSPARSE_OK, 'normal solve')
    call require(maxval(abs(solution - truth)) < 2.e-14_dp, 'known normal solution')
    solution = rhs
    call solver%solve_inplace(solution, status)
    call require(status%code == FORTSPARSE_OK, 'in-place normal status')
    call require(maxval(abs(solution - truth)) < 2.e-14_dp, 'in-place normal solution')
    call solver%solve(rhs, solution, status, operation=2)
    call require(maxval(abs(solution - truth)) < 2.e-14_dp, 'Hermitian solve')
    rhs = matmul(transpose(dense), truth)
    call solver%solve(rhs, solution, status, operation=1)
    call require(maxval(abs(solution - truth)) < 2.e-14_dp, 'transpose conjugation')
    solution = rhs
    call solver%solve_inplace(solution, status, operation=1)
    call require(status%code == FORTSPARSE_OK, 'in-place transpose status')
    call require(maxval(abs(solution - truth)) < 2.e-14_dp, &
                 'in-place transpose solution')
    solution = matmul(dense, truth)
    call solver%solve_inplace(solution, status, operation=2)
    call require(status%code == FORTSPARSE_OK, 'in-place adjoint status')
    call require(maxval(abs(solution - truth)) < 2.e-14_dp, 'in-place adjoint solution')
    call solver%solve_inplace(solution(:n - 1), status)
    call require(status%code == FORTSPARSE_INVALID_MATRIX, 'in-place shape guard')
    call solver%solve_inplace(solution, status, operation=3)
    call require(status%code == FORTSPARSE_INVALID_MATRIX, 'in-place operation guard')
    solution(1) = cmplx(ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
    call solver%solve_inplace(solution, status)
    call require(status%code == FORTSPARSE_INVALID_MATRIX, 'in-place finite RHS guard')
    do k = 1, 10
        reference(:, k) = truth*cmplx(real(k, dp), -0.2_dp, dp)
        blocks(:, k) = matmul(dense, reference(:, k))
    end do
    call solver%solve(blocks, answers, status)
    call require(status%code == FORTSPARSE_OK, 'block solve')
    call require(maxval(abs(answers - reference)) < 2.e-13_dp, 'known block solution')
    call solver%solve(v, solution, status)
    call solver%solve(u, adjoint, status, operation=2)
    lhs = real(sum(conjg(u)*solution), dp)
    right = real(sum(conjg(adjoint)*v), dp)
    call require(abs(lhs - right) < 4.e-15_dp, 'real-pairing adjoint identity')
    call solver%memory_payloads(fb, pb, wb)
    bytes = solver%memory_current_bytes()
    peak = solver%memory_peak_bytes()
    call require(bytes == fb + pb + wb, 'disjoint exact payload accounting')
    call require(bytes == 20_i8*solver%factor_nnz() + 44_i8*n + 4_i8*nz + 8_i8, &
                 'independent exact layout formula')
    call require(peak == bytes, 'no transient overlap')
    call capped%factor(A, status, memory_limit_bytes=bytes - 1_i8, low_memory=.true.)
    call require(status%code == FORTSPARSE_RESOURCE_LIMIT, 'cap rejects one byte short')
    call require(capped%memory_current_bytes() == 0_i8, 'failed allocation releases')
    call require(capped%memory_peak_bytes() < bytes, 'rejected allocation not charged')
    call capped%factor(A, status, memory_limit_bytes=bytes, low_memory=.true.)
    call require(status%code == FORTSPARSE_OK, 'exact boundary cap fits')
    A%val = 2.0_dp*A%val
    call solver%refactor(A, status, reused)
    call require(reused .and. status%code == FORTSPARSE_OK, 'same-pattern reuse')
    call require(solver%memory_current_bytes() == bytes, 'refactor no allocations')
    call require(solver%memory_peak_bytes() == peak, 'refactor no hidden overlap')
    rhs = 2.0_dp*matmul(dense, truth)
    call solver%solve(rhs, solution, status)
    call require(maxval(abs(solution - truth)) < 3.e-14_dp, 'refactor numeric values')
    call solver%refactor(A, status, reused, memory_limit_bytes=bytes - 1_i8)
    call require(status%code == FORTSPARSE_RESOURCE_LIMIT, 'tightened cap on reuse')
    call require(.not. solver%is_factored(), 'capped refactor invalidates readiness')
    call solver%refactor(A, status, reused, memory_limit_bytes=0_i8)
    call require(status%code == FORTSPARSE_OK .and. reused, 'recover cap failure')
    bad = A
    bad%val(1) = cmplx(-1.0_dp, 0.0_dp, dp)
    call solver%refactor(bad, status, reused)
    call require(status%code == FORTSPARSE_NOT_POSITIVE_DEFINITE, 'negative pivot')
    call solver%solve(rhs, solution, status)
    call require(status%code == FORTSPARSE_NOT_FACTORED, 'failed numeric unusable')
    bad = A
    bad%val(1) = cmplx(0.0_dp, 0.0_dp, dp)
    call solver%refactor(bad, status, reused)
    call require(status%code == FORTSPARSE_NOT_POSITIVE_DEFINITE, 'zero pivot')
    call solver%refactor(A, status, reused)
    call require(status%code == FORTSPARSE_OK, 'recover failed pivot')
    bad = A
    bad%val(1) = cmplx(real(A%val(1), dp), 0.1_dp, dp)
    call solver%factor(bad, status)
    call require(status%code == FORTSPARSE_INVALID_MATRIX, 'complex diagonal rejected')
    bad = A
    bad%val(2) = bad%val(2) + cmplx(0.0_dp, 0.1_dp, dp)
    call solver%factor(bad, status)
    call require(status%code == FORTSPARSE_INVALID_MATRIX, 'non-Hermitian rejected')
    bad = A
    bad%row_idx(2) = bad%row_idx(1)
    call solver%factor(bad, status)
    call require(status%code == FORTSPARSE_INVALID_MATRIX, 'duplicate rejected')
    bad = A
    bad%val(2) = cmplx(ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
    call solver%factor(bad, status)
    call require(status%code == FORTSPARSE_INVALID_MATRIX, 'nonfinite rejected')
    call solver%factor(A, status, memory_limit_bytes=-1_i8)
    call require(status%code == FORTSPARSE_INVALID_MATRIX, 'negative cap rejected')
    call solver%free()
    call capped%free()
    call require(solver%memory_current_bytes() == 0_i8, 'free exact payload zero')
    call changed_pattern_cap_fixture()
    call sparse_storage_fixture()
    print '(a)', 'native HPD known-solution, transpose, reuse, cap and failures: PASS'
contains
    subroutine changed_pattern_cap_fixture()
        type(native_hpd_t) :: factor
        type(csc_z_t) :: small, larger
        type(fortsparse_status_t) :: st
        complex(dp) :: rhs(3)
        integer(i8) :: required, available
        logical :: did_reuse
        integer(i8), parameter :: small_bytes = 144_i8, large_bytes = 212_i8

        call csc_from_triplet(2, 2, [1, 2], [1, 2], &
                              cmplx([2.0_dp, 3.0_dp], 0.0_dp, dp), small, st)
        call require(st%code == FORTSPARSE_OK, 'small diagonal fixture')
        call csc_from_triplet(3, 3, [1, 2, 3], [1, 2, 3], &
                              cmplx([2.0_dp, 3.0_dp, 4.0_dp], 0.0_dp, dp), larger, st)
        call require(st%code == FORTSPARSE_OK, 'larger diagonal fixture')
        call factor%factor(small, st, memory_limit_bytes=small_bytes, low_memory=.true.)
        call require(st%code == FORTSPARSE_OK, 'changed-pattern initial cap')
        call factor%refactor(larger, st, did_reuse)
       call require(st%code == FORTSPARSE_RESOURCE_LIMIT, 'changed pattern retains cap')
        call require(.not. did_reuse, 'changed pattern rebuilds analysis')
     call require(factor%memory_current_bytes() == 0_i8, 'changed cap releases payload')
     call require(factor%memory_peak_bytes() <= small_bytes, 'changed cap peak bounded')
        call factor%memory_failure_bytes(required, available)
        call require(required > available, 'denied allocation diagnostics')
        call factor%factor(small, st, memory_limit_bytes=0_i8)
        call require(st%code == FORTSPARSE_OK, 'changed-pattern reset fixture')
        call factor%refactor(larger, st, did_reuse, &
                             memory_limit_bytes=large_bytes - 1_i8, low_memory=.true.)
      call require(st%code == FORTSPARSE_RESOURCE_LIMIT, 'changed pattern supplied cap')
        call factor%refactor(larger, st, did_reuse, memory_limit_bytes=large_bytes)
        call require(st%code == FORTSPARSE_OK, 'changed pattern exact cap fits')
        call require(factor%memory_current_bytes() == large_bytes, 'new exact payload')
        rhs = cmplx([2.0_dp, 6.0_dp, 12.0_dp], 0.0_dp, dp)
        call factor%solve_inplace(rhs, st)
        call require(st%code == FORTSPARSE_OK, 'changed pattern in-place solve')
        call require(maxval(abs(rhs - cmplx([1.0_dp, 2.0_dp, 3.0_dp], 0.0_dp, dp))) &
                     < 1.e-15_dp, 'changed pattern known solution')
    call require(factor%memory_current_bytes() == large_bytes, 'in-place no allocation')
      call require(factor%memory_peak_bytes() == large_bytes, 'in-place no peak growth')
        call factor%free()
        call factor%solve_inplace(rhs, st)
        call require(st%code == FORTSPARSE_NOT_FACTORED, 'in-place readiness guard')
        larger%val(1) = cmplx(1.e-300_dp,0.0_dp,dp)
        call factor%factor(larger,st)
        call require(st%code == FORTSPARSE_OK, 'small positive diagonal fixture')
        rhs = cmplx(1.0_dp,0.0_dp,dp)
        rhs(1) = cmplx(huge(1.0_dp),0.0_dp,dp)
        call factor%solve_inplace(rhs,st)
        call require(st%code == FORTSPARSE_INTERNAL_ERROR, &
                     'in-place finite result guard')
    end subroutine changed_pattern_cap_fixture

    subroutine sparse_storage_fixture()
        integer, parameter :: large_n = 2000, large_nz = 3*large_n - 2
        type(native_hpd_t) :: sparse
        type(csc_z_t) :: matrix, empty
        type(fortsparse_status_t) :: st
        integer :: r(large_nz), c(large_nz), k, p
        complex(dp) :: a(large_nz), x(large_n), b(large_n), y(large_n)
        complex(dp) :: nothing(0), nothing_out(0)
        complex(dp), parameter :: offdiag = cmplx(-0.2_dp, 0.3_dp, dp)
        integer(i8) :: payload

        p = 0
        do k = 1, large_n
            x(k) = cmplx(sin(real(k, dp)), cos(real(k, dp)), dp)
            p = p + 1
            r(p) = k
            c(p) = k
            a(p) = cmplx(2.0_dp, 0.0_dp, dp)
            if (k == large_n) cycle
            p = p + 1
            r(p) = k
            c(p) = k + 1
            a(p) = offdiag
            p = p + 1
            r(p) = k + 1
            c(p) = k
            a(p) = conjg(offdiag)
        end do
        ! Independent tridiagonal action; no dense oracle or global Jacobian.
        b = 2.0_dp*x
        do k = 1, large_n - 1
            b(k) = b(k) + offdiag*x(k + 1)
            b(k + 1) = b(k + 1) + conjg(offdiag)*x(k)
        end do
        call csc_from_triplet(large_n, large_n, r, c, a, matrix, st)
        call require(st%code == FORTSPARSE_OK, 'large sparse fixture')
        payload = 96_i8*large_n - 20_i8
        call sparse%factor(matrix, st, memory_limit_bytes=payload, low_memory=.true.)
        call require(st%code == FORTSPARSE_OK, 'sparse linear-memory cap')
        call require(sparse%factor_nnz() == 2_i8*large_n - 1_i8, 'tridiagonal fill')
        call require(sparse%memory_current_bytes() == payload, 'linear memory exact')
        call sparse%solve(b, y, st)
        call require(maxval(abs(x - y)) < 1.e-15_dp, 'large sparse known solution')
        allocate (empty%col_ptr(1), empty%row_idx(0), empty%val(0))
        empty%col_ptr = 1
        call sparse%factor(empty, st)
        call require(st%code == FORTSPARSE_OK, 'empty HPD factor')
        call sparse%solve(nothing, nothing_out, st)
        call require(st%code == FORTSPARSE_OK, 'empty solve')
        call require(sparse%memory_current_bytes() == 8_i8, 'empty metadata bytes')
        call sparse%free()
    end subroutine sparse_storage_fixture

    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (.not. condition) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
    end subroutine require
end program test_fortsparse_native_hpd
