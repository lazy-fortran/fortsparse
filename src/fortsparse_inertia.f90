module fortsparse_inertia
    ! Sylvester inertia of a symmetric/Hermitian block-diagonal pivot matrix D
    ! from P^T A P = L D L^T (or L D L^H) with 1x1 and 2x2 pivots. The counts
    ! are the signs of the actual computed pivots; reliability compares every
    ! pivot eigenvalue magnitude with a rounding threshold (docs/native-cpu.md).
    use fortsparse_kinds, only: dp, i8
    implicit none
    private

    ! sparse_factor_stats_t%inertia_status and sparse_inertia reliability.
    ! UNAVAILABLE: the factor family does not provide inertia (LU, perturbed
    ! pivots, complex symmetric, external backends) or no factorization ran.
    integer, parameter, public :: FORTSPARSE_INERTIA_UNAVAILABLE = -1
    ! Every pivot eigenvalue magnitude exceeds the rounding threshold.
    integer, parameter, public :: FORTSPARSE_INERTIA_RELIABLE = 0
    ! Some nonzero pivot eigenvalue lies at or below the threshold; its sign,
    ! and so the counts, may differ from the exact inertia of A.
    integer, parameter, public :: FORTSPARSE_INERTIA_UNCERTAIN = 1
    ! An exactly zero pivot column occurred: n_zero > 0, solves unavailable.
    integer, parameter, public :: FORTSPARSE_INERTIA_SINGULAR = 2

    ! Threshold multiplier: tolerance = INERTIA_SAFETY * n * eps * growth,
    ! growth = max |entry| of A and of every formed Schur complement entry.
    real(dp), parameter, public :: INERTIA_SAFETY = 1.0_dp

    type, public :: inertia_count_t
        integer(i8) :: negative = 0_i8, zero = 0_i8, positive = 0_i8
        integer(i8) :: uncertain = 0_i8 ! nonzero pivot eigenvalues <= tolerance
        real(dp) :: min_magnitude = huge(1.0_dp) ! smallest pivot eigenvalue
    end type inertia_count_t

    public :: inertia_add_pivot, inertia_status_of, inertia_tolerance
contains
    ! Add one 1x1 pivot a (step 1) or 2x2 pivot [[a, b], [conj(b), c]] with
    ! real a, c and |b| = babs (step 2). Eigenvalues are formed after scaling
    ! by the largest entry: the larger-magnitude one is m + sign(m) r, the
    ! smaller one det/larger, so neither suffers cancellation.
    pure subroutine inertia_add_pivot(count, step, a, babs, c, tolerance)
        type(inertia_count_t), intent(inout) :: count
        integer, intent(in) :: step
        real(dp), intent(in) :: a, babs, c, tolerance
        real(dp) :: s, as, bs, cs, m, r, large, small
        if (step == 1) then
            call add_eigenvalue(count, a, tolerance)
            return
        end if
        s = max(abs(a), babs, abs(c))
        if (s <= 0.0_dp) then
            call add_eigenvalue(count, 0.0_dp, tolerance)
            call add_eigenvalue(count, 0.0_dp, tolerance)
            return
        end if
        as = a/s
        bs = babs/s
        cs = c/s
        m = 0.5_dp*(as + cs)
        r = hypot(0.5_dp*(as - cs), bs)
        large = m + sign(r, m)
        small = (as*cs - bs*bs)/large
        call add_eigenvalue(count, s*large, tolerance)
        call add_eigenvalue(count, s*small, tolerance)
    end subroutine inertia_add_pivot

    pure subroutine add_eigenvalue(count, lambda, tolerance)
        type(inertia_count_t), intent(inout) :: count
        real(dp), intent(in) :: lambda, tolerance
        count%min_magnitude = min(count%min_magnitude, abs(lambda))
        if (lambda > 0.0_dp) then
            count%positive = count%positive + 1_i8
        else if (lambda < 0.0_dp) then
            count%negative = count%negative + 1_i8
        else
            count%zero = count%zero + 1_i8
            return
        end if
        if (abs(lambda) <= tolerance) count%uncertain = count%uncertain + 1_i8
    end subroutine add_eigenvalue

    pure real(dp) function inertia_tolerance(n, growth) result(tolerance)
        integer, intent(in) :: n
        real(dp), intent(in) :: growth
        tolerance = INERTIA_SAFETY*real(n, dp)*epsilon(1.0_dp)*growth
    end function inertia_tolerance

    pure integer function inertia_status_of(count) result(state)
        type(inertia_count_t), intent(in) :: count
        state = FORTSPARSE_INERTIA_RELIABLE
        if (count%uncertain > 0_i8) state = FORTSPARSE_INERTIA_UNCERTAIN
        if (count%zero > 0_i8) state = FORTSPARSE_INERTIA_SINGULAR
    end function inertia_status_of
end module fortsparse_inertia
