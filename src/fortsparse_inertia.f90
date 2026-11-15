module fortsparse_inertia
    ! Sylvester inertia of a symmetric/Hermitian block-diagonal pivot matrix D
    ! from P^T A P = L D L^T (or L D L^H) with 1x1 and 2x2 pivots. The counts
    ! are the signs of the actual computed pivots; reliability compares every
    ! pivot eigenvalue magnitude with a local rounding threshold
    ! (docs/native-cpu.md).
    use fortsparse_kinds, only: dp, i8
    implicit none
    private

    ! sparse_factor_stats_t%inertia_status and sparse_inertia reliability.
    ! UNAVAILABLE: the factor family does not provide inertia (LU, perturbed
    ! pivots, complex symmetric, external backends) or no factorization ran.
    integer, parameter, public :: FORTSPARSE_INERTIA_UNAVAILABLE = -1
    ! Every pivot eigenvalue magnitude exceeds its rounding threshold.
    integer, parameter, public :: FORTSPARSE_INERTIA_RELIABLE = 0
    ! Some nonzero pivot eigenvalue lies at or below its threshold; its sign,
    ! and so the counts, may differ from the exact inertia of A.
    integer, parameter, public :: FORTSPARSE_INERTIA_UNCERTAIN = 1
    ! An exactly zero pivot column occurred: n_zero > 0, solves unavailable.
    integer, parameter, public :: FORTSPARSE_INERTIA_SINGULAR = 2

    ! Threshold multiplier: tolerance = INERTIA_SAFETY * n * eps * row_scale,
    ! row_scale = max |entry| of A and of every Schur entry formed in the
    ! pivot rows.
    real(dp), parameter, public :: INERTIA_SAFETY = 1.0_dp

    type, public :: inertia_count_t
        integer(i8) :: negative = 0_i8, zero = 0_i8, positive = 0_i8
        integer(i8) :: uncertain = 0_i8 ! nonzero pivot eigenvalues <= threshold
        real(dp) :: min_magnitude = huge(1.0_dp) ! smallest pivot eigenvalue
        real(dp) :: tolerance = -1.0_dp ! threshold of that smallest pivot
    end type inertia_count_t

    public :: inertia_add_pivot, inertia_add_scaled_block, inertia_status_of, inertia_tolerance
contains
    ! Add one 1x1 pivot a (step 1) or 2x2 pivot [[a, b], [conj(b), c]] with
    ! real a, c and |b| = babs (step 2), both eigenvalues against tolerance.
    pure subroutine inertia_add_pivot(count, step, a, babs, c, tolerance)
        type(inertia_count_t), intent(inout) :: count
        integer, intent(in) :: step
        real(dp), intent(in) :: a, babs, c, tolerance
        real(dp) :: lambda(2)
        if (step == 1) then
            call add_eigenvalue(count, a, tolerance)
            return
        end if
        call block_eigenvalues(a, babs, c, lambda)
        call add_eigenvalue(count, lambda(1), tolerance)
        call add_eigenvalue(count, lambda(2), tolerance)
    end subroutine inertia_add_pivot

    ! Add a 2x2 pivot B = [[a, b], [conj(b), c]] of the rows i, j whose
    ! computed entries carry rounding errors |F_ii| <= gamma g_i,
    ! |F_jj| <= gamma g_j and |F_ij| <= gamma sqrt(g_i g_j), g = scale_i,
    ! scale_j, gamma ~ n u = t/2 with t = n eps (inertia_tolerance per unit
    ! scale). With G = diag(g_i, g_j), the congruence G^-1/2 B G^-1/2 has the
    ! inertia of B (Sylvester), and its error G^-1/2 F G^-1/2 has entries at
    ! most t/2, so 2-norm at most t. By Weyl, an eigenvalue mu of the scaled
    ! block with |mu| > t keeps its sign: the margin of inertia_add_pivot
    ! with tolerance t max(g_i, g_j), whose bound |F| <= gamma max(g_i, g_j)
    ! the scaled one implies. Each eigenvalue lambda of B (ascending, paired
    ! with the ascending mu by Ostrowski's theorem) is counted against
    ! t |lambda/mu|, so that it is uncertain exactly when |mu| <= t. Rows of
    ! equal scale give the absolute rule. Nonpositive scales fall back to
    ! it; a zero mu or a sign that the scaling changed counts uncertain.
    pure subroutine inertia_add_scaled_block(count, n, a, babs, c, scale_i, scale_j)
        type(inertia_count_t), intent(inout) :: count
        integer, intent(in) :: n
        real(dp), intent(in) :: a, babs, c, scale_i, scale_j
        real(dp) :: lambda(2), mu(2), ri, rj, t, tolerance
        integer :: k
        if (.not. (scale_i > 0.0_dp .and. scale_j > 0.0_dp)) then
            call inertia_add_pivot(count, 2, a, babs, c, &
                inertia_tolerance(n, max(scale_i, scale_j)))
            return
        end if
        t = inertia_tolerance(n, 1.0_dp)
        ri = sqrt(scale_i)
        rj = sqrt(scale_j)
        call block_eigenvalues(a, babs, c, lambda)
        call block_eigenvalues((a/ri)/ri, (babs/ri)/rj, (c/rj)/rj, mu)
        do k = 1, 2
            tolerance = huge(1.0_dp)
            if (lambda(k) > 0.0_dp .eqv. mu(k) > 0.0_dp) then
                if (mu(k) /= 0.0_dp) tolerance = min(huge(1.0_dp), t*abs(lambda(k)/mu(k)))
            end if
            call add_eigenvalue(count, lambda(k), tolerance)
        end do
    end subroutine inertia_add_scaled_block

    ! Eigenvalues of [[a, b], [conj(b), c]], |b| = babs, in ascending order.
    ! They are formed after scaling by the largest entry: the larger-magnitude
    ! one is m + sign(m) r, the smaller one det/larger, so neither suffers
    ! cancellation.
    pure subroutine block_eigenvalues(a, babs, c, lambda)
        real(dp), intent(in) :: a, babs, c
        real(dp), intent(out) :: lambda(2)
        real(dp) :: s, as, bs, cs, m, r, large, small
        lambda = 0.0_dp
        s = max(abs(a), babs, abs(c))
        if (s <= 0.0_dp) return
        as = a/s
        bs = babs/s
        cs = c/s
        m = 0.5_dp*(as + cs)
        r = hypot(0.5_dp*(as - cs), bs)
        large = m + sign(r, m)
        small = (as*cs - bs*bs)/large
        lambda = [min(s*large, s*small), max(s*large, s*small)]
    end subroutine block_eigenvalues

    pure subroutine add_eigenvalue(count, lambda, tolerance)
        type(inertia_count_t), intent(inout) :: count
        real(dp), intent(in) :: lambda, tolerance
        if (abs(lambda) < count%min_magnitude) then
            count%min_magnitude = abs(lambda)
            count%tolerance = tolerance
        end if
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

    pure real(dp) function inertia_tolerance(n, row_scale) result(tolerance)
        integer, intent(in) :: n
        real(dp), intent(in) :: row_scale
        tolerance = INERTIA_SAFETY*real(n, dp)*epsilon(1.0_dp)*row_scale
    end function inertia_tolerance

    pure integer function inertia_status_of(count) result(state)
        type(inertia_count_t), intent(in) :: count
        state = FORTSPARSE_INERTIA_RELIABLE
        if (count%uncertain > 0_i8) state = FORTSPARSE_INERTIA_UNCERTAIN
        if (count%zero > 0_i8) state = FORTSPARSE_INERTIA_SINGULAR
    end function inertia_status_of
end module fortsparse_inertia
