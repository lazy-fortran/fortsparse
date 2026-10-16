program test_fortsparse_slu_fallback
    ! Gates, mandatory validation and scalar fallback of the supernodal LU.
    ! Block-diagonal 2x2 matrices (n = 1200, perturbation limit 1): a special
    ! block [1 .9; .9 1] has the static second pivot 0.19 and ||A||_inf = 1.9
    ! (the identity matching and unit scaling are optimal). With
    ! pivot_perturbation just above 0.1 that pivot is perturbed by a relative
    ! 1e-6, so refinement converges; with twice the threshold the
    ! perturbation halves the convergence rate and the validation fails.
    use fortsparse
    implicit none
    integer, parameter :: n = 1200
    real(dp), parameter :: tight = 0.1_dp*(1.0_dp+1.0e-6_dp)
    type(sparse_factor_stats_t) :: stats
    type(fortsparse_status_t) :: status
    type(csc_t) :: a
    integer(i8) :: supernodal_peak, scalar_peak
    integer :: seed

    seed = 11
    ! A non-binding cap selects the strict (capped) layout whose peak is the
    ! exact admission threshold of later capped runs.
    call blocks(1, a)
    call run(a, FORTSPARSE_LU_SUPERNODAL, tight, huge(1_i8)/4, status, stats)
    call check(status_ok(status), 'one perturbation within the limit: '//trim(status%msg))
    call check(stats%native_factor_family == 'LU_SUPERNODAL', 'supernodal kept')
    call check(stats%perturbed_pivots == 1, 'exactly one perturbed pivot')
    call check(.not. stats%lu_fallback, 'validation passes after refinement')
    supernodal_peak = stats%host_peak_bytes
    call check(supernodal_peak > 0, 'exact supernodal peak')
    ! Without refine, a perturbed factor retains the original A for its
    ! refinement; when that no longer fits the cap the scalar LU is used.
    call run(a, FORTSPARSE_LU_SUPERNODAL, tight, supernodal_peak-1, status, stats)
    call check(status_ok(status) .and. stats%lu_fallback_reason == 'PERTURBATION', &
        'refinement retention beyond the cap falls back')
    call run(a, FORTSPARSE_LU_SUPERNODAL, tight, supernodal_peak, status, stats)
    call check(status_ok(status) .and. .not. stats%lu_fallback, &
        'supernodal fits its own measured peak')
    call check(stats%host_peak_bytes <= supernodal_peak, 'capped peak within the cap')

    call run(a, FORTSPARSE_LU_SUPERNODAL, 2.0_dp*tight, 0_i8, status, stats)
    call check(status_ok(status), 'refinement fallback succeeds')
    call check(stats%lu_fallback .and. stats%lu_fallback_reason == 'REFINEMENT', &
        'fallback reason REFINEMENT')

    call blocks(2, a)
    call run(a, FORTSPARSE_LU_SUPERNODAL, tight, 0_i8, status, stats)
    call check(status_ok(status), 'perturbation fallback succeeds')
    call check(stats%lu_fallback .and. stats%lu_fallback_reason == 'PERTURBATION', &
        'fallback reason PERTURBATION')
    call check(stats%native_factor_family == 'LU_REAL', 'fallback family')
    call check(stats%lu_method_used == FORTSPARSE_LU_SCALAR, 'fallback method')
    call check(stats%matching_used == FORTSPARSE_MATCHING_STRUCTURAL, 'fallback matching')
    call check(stats%perturbed_pivots == 2, 'attempt counters kept after fallback')

    ! Banded, badly scaled: threshold pivoting needs more than the static
    ! structure, so a cap admitting the supernodal attempt refuses the
    ! fallback. A large pivot_perturbation forces the fallback.
    call banded(300, a)
    call run(a, FORTSPARSE_LU_SUPERNODAL, 1.0e-13_dp, huge(1_i8)/4, status, stats)
    call check(status_ok(status) .and. .not. stats%lu_fallback, 'banded supernodal factor')
    supernodal_peak = stats%host_peak_bytes
    call run(a, FORTSPARSE_LU_SUPERNODAL, 1.0e-13_dp, supernodal_peak-1, status, stats)
    call check(status%code == FORTSPARSE_RESOURCE_LIMIT, 'supernodal refused under a cap')
    call check(.not. stats%lu_fallback, 'cap refusal is not a numerical fallback')
    call check(stats%budget_required_bytes > stats%budget_available_bytes, &
        'supernodal refusal reports required/available')
    call run(a, FORTSPARSE_LU_SUPERNODAL, 1.0e-13_dp, supernodal_peak, status, stats)
    call check(status_ok(status) .and. .not. stats%lu_fallback, &
        'banded supernodal fits its measured peak')
    call run(a, FORTSPARSE_LU_SCALAR, 0.0_dp, -1_i8, status, stats)
    scalar_peak = stats%host_peak_bytes
    call check(scalar_peak > supernodal_peak, 'scalar fallback needs more memory here')
    call run(a, FORTSPARSE_LU_SUPERNODAL, 0.5_dp, 0_i8, status, stats)
    call check(status_ok(status) .and. stats%lu_fallback_reason == 'PERTURBATION', &
        'forced fallback without a cap')
    call run(a, FORTSPARSE_LU_SUPERNODAL, 0.5_dp, supernodal_peak, status, stats)
    call check(status%code == FORTSPARSE_RESOURCE_LIMIT, 'fallback refused under the cap')
    call check(stats%lu_fallback .and. stats%lu_fallback_reason == 'PERTURBATION', &
        'refused fallback keeps its reason')
    call check(stats%budget_required_bytes > stats%budget_available_bytes, &
        'fallback refusal reports required/available')
    print '(a)', 'supernodal LU gates and fallback passed'
contains
    real(dp) function uniform()
        seed = modulo(1103515245_i8*seed+12345_i8, 2147483648_i8)
        uniform = real(seed, dp)/2147483648.0_dp
    end function uniform

    subroutine blocks(special, a)
        integer, intent(in) :: special
        type(csc_t), intent(out) :: a
        integer :: rows(2*n), cols(2*n), k, j
        real(dp) :: values(2*n), off
        type(fortsparse_status_t) :: status
        k = 0
        do j = 1, n, 2
            off = 0.1_dp
            if (j/2 < special) off = 0.9_dp
            rows(k+1:k+4) = [j, j+1, j, j+1]
            cols(k+1:k+4) = [j, j, j+1, j+1]
            values(k+1:k+4) = [1.0_dp, off, off, 1.0_dp]
            k = k+4
        end do
        call csc_from_triplet(n, n, rows, cols, values, a, status)
        call check(status_ok(status), 'block matrix construction')
    end subroutine blocks

    subroutine banded(m, a)
        ! Near-symmetric band of half width 6, magnitudes over three decades,
        ! weak diagonal where no entry was drawn.
        integer, intent(in) :: m
        type(csc_t), intent(out) :: a
        integer, allocatable :: r(:), c(:)
        real(dp), allocatable :: v(:)
        logical, allocatable :: used(:,:)
        integer :: k, i, j
        type(fortsparse_status_t) :: status
        allocate(r(6*m), c(6*m), v(6*m), used(m,m))
        used = .false.
        k = 0
        do while (k < 5*m-1)
            i = 1+int(uniform()*m)
            j = min(m, max(1, i+int((uniform()-0.5_dp)*12.0_dp)))
            if (used(i,j) .or. used(j,i)) cycle
            used(i,j) = .true.
            k = k+1
            r(k) = i
            c(k) = j
            v(k) = (2.0_dp*uniform()-1.0_dp)*10.0_dp**(3.0_dp*uniform())
            if (uniform() < 0.9_dp .and. i /= j) then
                used(j,i) = .true.
                k = k+1
                r(k) = j
                c(k) = i
                v(k) = (2.0_dp*uniform()-1.0_dp)*10.0_dp**(3.0_dp*uniform())
            end if
        end do
        do j = 1, m
            if (used(j,j)) cycle
            k = k+1
            r(k) = j
            c(k) = j
            v(k) = 1.0e-3_dp
        end do
        call csc_from_triplet(m, m, r(:k), c(:k), v(:k), a, status)
        call check(status_ok(status), 'banded matrix construction')
    end subroutine banded

    subroutine run(a, method, tau, limit, status, stats)
        ! limit < 0: scalar fallback settings under a non-binding cap, which
        ! selects the bounded AMD so that the reported peak is exact.
        type(csc_t), intent(in) :: a
        integer, intent(in) :: method
        real(dp), intent(in) :: tau
        integer(i8), intent(in) :: limit
        type(fortsparse_status_t), intent(out) :: status
        type(sparse_factor_stats_t), intent(out) :: stats
        type(sparse_solver_t) :: solver
        type(fortsparse_status_t) :: factor_status
        real(dp), allocatable :: exact(:), b(:), x(:), r(:)
        integer :: j, k
        solver%backend_id = FORTSPARSE_BACKEND_NATIVE
        solver%options%lu_method = method
        if (method == FORTSPARSE_LU_SUPERNODAL) solver%options%pivot_perturbation = tau
        solver%options%host_memory_limit_bytes = limit
        if (limit < 0_i8) then
            solver%options%structural_matching = .true.
            solver%options%pivot_tolerance = 0.1_dp
            solver%options%ordering = FORTSPARSE_ORDER_AMD
            solver%options%host_memory_limit_bytes = huge(1_i8)/4
        end if
        call sparse_factor(solver, a, factor_status)
        call sparse_factor_stats(solver, stats, status)
        if (status_ok(factor_status)) then
            allocate(exact(a%ncol), b(a%ncol), x(a%ncol), r(a%ncol))
            do j = 1, a%ncol
                exact(j) = 1.0_dp+real(mod(3*j, 11), dp)
            end do
            b = 0.0_dp
            do j = 1, a%ncol
                do k = a%col_ptr(j), a%col_ptr(j+1)-1
                    b(a%row_idx(k)) = b(a%row_idx(k))+a%val(k)*exact(j)
                end do
            end do
            call sparse_solve(solver, b, x, status)
            call check(status_ok(status), 'solve')
            r = b
            do j = 1, a%ncol
                do k = a%col_ptr(j), a%col_ptr(j+1)-1
                    r(a%row_idx(k)) = r(a%row_idx(k))-a%val(k)*x(j)
                end do
            end do
            ! Independent normwise backward error of the returned solution.
            call check(maxval(abs(r)) <= 64.0_dp*epsilon(1.0_dp)*(maxval(abs(b)) &
                +norm_inf(a)*maxval(abs(x))), 'independent original-matrix residual')
        end if
        status = factor_status
        call sparse_destroy(solver)
    end subroutine run

    real(dp) function norm_inf(a)
        type(csc_t), intent(in) :: a
        real(dp), allocatable :: sums(:)
        integer :: j, k
        allocate(sums(a%nrow))
        sums = 0.0_dp
        do j = 1, a%ncol
            do k = a%col_ptr(j), a%col_ptr(j+1)-1
                sums(a%row_idx(k)) = sums(a%row_idx(k))+abs(a%val(k))
            end do
        end do
        norm_inf = maxval(sums)
    end function norm_inf

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (condition) return
        print '(a)', 'FAIL: '//message
        error stop 1
    end subroutine check
end program test_fortsparse_slu_fallback
