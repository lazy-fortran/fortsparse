module fortsparse_ams
    !! Auxiliary-space Maxwell solver (AMS) preconditioner for H(curl)
    !! matrices A = C^T M C (+ beta M1): smoothing on A and corrections in the
    !! gradient space G (n1 x n0) and the nodal vector space Pi (n1 x 3 n0,
    !! x, y, z column blocks), each solved by smoothed aggregation AMG.
    !! Hiptmair and Xu, SIAM J. Numer. Anal. 45, 2483-2509 (2007); Kolev and
    !! Vassilevski, J. Comput. Math. 27, 604-623 (2009), cycle "01210".
    use fortsparse_kinds, only: dp
    use fortsparse_csr, only: csr_t, csr_budget_t, csr_matvec, csr_transpose, &
        csr_multiply, csr_valid
    use fortsparse_status, only: fortsparse_status_t, status_set, FORTSPARSE_OK, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_NOT_FACTORED
    use fortsparse_amg, only: amg_params_t, amg_handle_t, amg_setup, amg_apply, &
        amg_destroy
    use fortsparse_smoothers, only: csr_smoother_t, smoother_setup, smoother_apply, &
        smoother_destroy, SMOOTHER_L1_JACOBI
    implicit none
    private
    public :: ams_params_t, ams_handle_t, ams_setup_csr, ams_setup, ams_apply, &
        ams_destroy, ams_uses_gradient
    integer, parameter, public :: AMS_MULTIPLICATIVE = 1, AMS_ADDITIVE = 2

    type :: ams_params_t
        integer :: cycle = AMS_MULTIPLICATIVE
        integer :: smoothing_steps = 1
        !! The auxiliary operators G^T A G and Pi^T A Pi get this multiple of
        !! their diagonal added, so that a semidefinite A (curl-curl without a
        !! mass term) still gives definite AMG problems.
        real(dp) :: aux_shift = 1.0e-6_dp
        !! True: three scalar AMG solves on Pi_d^T A Pi_d, d = x, y, z (block
        !! diagonal); false: one AMG on the vector problem Pi^T A Pi.
        logical :: nodal_components = .true.
        !! With the 1-form mass matrix M given to setup: A = C^T M2 C + beta M,
        !! and the nodal components are solved with G^T M G + beta Pi_d^T M Pi_d.
        real(dp) :: beta = 0
        type(amg_params_t) :: amg_gradient, amg_nodal
    end type

    type :: ams_handle_t
        private
        type(csr_t) :: a, g, gt, pi, pit, ag, api
        ! Allocatable for the same finalization reason as in fortsparse_pcg.
        type(amg_handle_t), allocatable :: bg, bpi, bcomp(:)
        type(csr_t) :: acomp(3)
        logical :: components = .false.
        integer :: n0 = 0
        type(csr_smoother_t) :: smoother
        integer :: cycle = AMS_MULTIPLICATIVE, steps = 1
        logical :: use_gradient = .false., ready = .false.
        real(dp), allocatable :: res(:), z0(:), r0(:), z0n(:), rpi(:), zpi(:)
    end type

contains

    subroutine ams_setup_csr(n1, a_ptr, a_col, a_val, n0, g_ptr, g_col, g_val, &
            pi_ptr, pi_col, pi_val, params, handle, info)
        !! One-based CSR arrays: A (n1 x n1, both triangles), G (n1 x n0) and
        !! Pi (n1 x 3 n0), columns sorted in each row.
        integer, intent(in) :: n1, n0
        integer, intent(in) :: a_ptr(:), a_col(:), g_ptr(:), g_col(:), pi_ptr(:), pi_col(:)
        real(dp), intent(in) :: a_val(:), g_val(:), pi_val(:)
        type(ams_params_t), intent(in) :: params
        type(ams_handle_t), intent(inout) :: handle
        type(fortsparse_status_t), intent(out) :: info
        type(csr_t) :: a, g, pi

        a = csr_t(n1, n1, a_ptr, a_col, a_val)
        g = csr_t(n1, n0, g_ptr, g_col, g_val)
        pi = csr_t(n1, 3*n0, pi_ptr, pi_col, pi_val)
        call ams_setup(a, g, pi, params, handle, info)
    end subroutine ams_setup_csr

    subroutine ams_setup(a, g, pi, params, handle, info, mass)
        !! Each Cartesian component of the nodal space is solved by AMG: with
        !! the 1-form mass matrix `mass` (n1 x n1) on G^T mass G + beta
        !! Pi_d^T mass Pi_d (params%beta), the vector Laplacian auxiliary
        !! operator of Hiptmair and Xu; without it on the block Pi_d^T A Pi_d,
        !! which for a curl-curl A lacks a derivative and suits AMG poorly.
        type(csr_t), intent(in) :: a, g, pi
        type(ams_params_t), intent(in) :: params
        type(ams_handle_t), intent(inout) :: handle
        type(fortsparse_status_t), intent(out) :: info
        type(csr_t), intent(in), optional :: mass
        type(csr_t) :: laplacian, nodal_mass, block
        type(csr_budget_t) :: budget
        type(csr_t) :: work
        type(amg_params_t) :: nodal
        integer :: n0, i

        call ams_destroy(handle)
        call status_set(info, FORTSPARSE_INVALID_MATRIX, 'AMS: invalid A, G or Pi')
        if (.not. (csr_valid(a) .and. csr_valid(g) .and. csr_valid(pi))) return
        if (a%ncol /= a%nrow .or. g%nrow /= a%nrow .or. pi%nrow /= a%nrow) return
        n0 = g%ncol
        if (pi%ncol /= 3*n0 .or. params%smoothing_steps < 1) return
        if (params%cycle /= AMS_MULTIPLICATIVE .and. params%cycle /= AMS_ADDITIVE) return
        handle%a = a
        handle%g = g
        handle%pi = pi
        handle%cycle = params%cycle
        handle%steps = params%smoothing_steps

        call smoother_setup(a, SMOOTHER_L1_JACOBI, handle%smoother, info)
        if (info%code /= FORTSPARSE_OK) return

        ! Galerkin operators of the auxiliary spaces.
        call csr_transpose(g, handle%gt, budget, info)
        if (info%code /= FORTSPARSE_OK) return
        call csr_transpose(pi, handle%pit, budget, info)
        if (info%code /= FORTSPARSE_OK) return
        call csr_multiply(a, g, work, budget, info)
        if (info%code /= FORTSPARSE_OK) return
        call csr_multiply(handle%gt, work, handle%ag, budget, info)
        if (info%code /= FORTSPARSE_OK) return
        call csr_multiply(a, pi, work, budget, info)
        if (info%code /= FORTSPARSE_OK) return
        call csr_multiply(handle%pit, work, handle%api, budget, info)
        if (info%code /= FORTSPARSE_OK) return

        ! A G = 0 (curl-curl without mass): the gradient correction of a
        ! consistent residual vanishes, so it is left out.
        handle%use_gradient = maxval(abs(handle%ag%val)) > &
            1.0e-12_dp*maxval(abs(a%val))*maxval(abs(g%val))**2
        if (handle%use_gradient) then
            call regularize(handle%ag, params%aux_shift)
            allocate (handle%bg)
            call amg_setup(handle%ag, params%amg_gradient, handle%bg, info)
            if (info%code /= FORTSPARSE_OK) then
                info%msg = 'AMS: gradient AMG: '//trim(info%msg)
                return
            end if
        end if
        handle%components = params%nodal_components
        handle%n0 = n0
        if (handle%components) then
            allocate (handle%bcomp(3))
            if (present(mass)) then
                call csr_multiply(mass, g, work, budget, info)
                if (info%code /= FORTSPARSE_OK) return
                call csr_multiply(handle%gt, work, laplacian, budget, info)
                if (info%code /= FORTSPARSE_OK) return
            end if
            if (present(mass) .and. params%beta /= 0) then
                call csr_multiply(mass, pi, work, budget, info)
                if (info%code /= FORTSPARSE_OK) return
                call csr_multiply(handle%pit, work, nodal_mass, budget, info)
                if (info%code /= FORTSPARSE_OK) return
                nodal_mass%val = params%beta*nodal_mass%val
            end if
            do i = 1, 3
                if (present(mass)) then
                    ! The vector Laplacian of Hiptmair and Xu, component by
                    ! component: the blocks of Pi^T A Pi lack a derivative.
                    handle%acomp(i) = laplacian
                    if (params%beta /= 0) then
                        call diagonal_block(nodal_mass, (i - 1)*n0, n0, block)
                        call add_into(handle%acomp(i), block)
                    end if
                else
                    call diagonal_block(handle%api, (i - 1)*n0, n0, handle%acomp(i))
                end if
                call regularize(handle%acomp(i), params%aux_shift)
                call amg_setup(handle%acomp(i), params%amg_nodal, handle%bcomp(i), info)
                if (info%code /= FORTSPARSE_OK) then
                    info%msg = 'AMS: nodal component AMG: '//trim(info%msg)
                    return
                end if
            end do
            allocate (handle%res(a%nrow), handle%z0(a%nrow), handle%r0(n0), &
                handle%z0n(n0), handle%rpi(3*n0), handle%zpi(3*n0))
            handle%ready = .true.
            call status_set(info, FORTSPARSE_OK, '')
            return
        end if
        call regularize(handle%api, params%aux_shift)
        nodal = params%amg_nodal
        if (.not. allocated(nodal%near_nullspace)) then
            ! The Cartesian constant vectors, one per component block.
            allocate (nodal%near_nullspace(3*n0, 3))
            nodal%near_nullspace = 0
            do i = 1, 3
                nodal%near_nullspace((i - 1)*n0 + 1:i*n0, i) = 1
            end do
        end if
        allocate (handle%bpi)
        call amg_setup(handle%api, nodal, handle%bpi, info)
        if (info%code /= FORTSPARSE_OK) then
            info%msg = 'AMS: nodal AMG: '//trim(info%msg)
            return
        end if
        allocate (handle%res(a%nrow), handle%z0(a%nrow), handle%r0(n0), handle%z0n(n0), &
            handle%rpi(3*n0), handle%zpi(3*n0))
        handle%ready = .true.
        call status_set(info, FORTSPARSE_OK, '')
    end subroutine ams_setup

    logical function ams_uses_gradient(handle)
        !! False when A G = 0 and the gradient correction is left out.
        type(ams_handle_t), intent(in) :: handle
        ams_uses_gradient = handle%use_gradient
    end function ams_uses_gradient

    subroutine ams_apply(handle, r, z, info)
        !! z = B r, a symmetric positive definite approximation of A^{-1} on
        !! the range of A. One handle owns its workspace.
        type(ams_handle_t), intent(inout) :: handle
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        type(fortsparse_status_t), intent(out) :: info
        integer :: step

        z = 0
        call status_set(info, FORTSPARSE_NOT_FACTORED, 'AMS: setup required')
        if (.not. handle%ready) return
        call status_set(info, FORTSPARSE_INVALID_MATRIX, 'AMS: vector shape')
        if (size(r) /= handle%a%nrow .or. size(z) /= size(r)) return
        call status_set(info, FORTSPARSE_OK, '')
        if (handle%cycle == AMS_ADDITIVE) then
            do step = 1, handle%steps
                call smoother_apply(handle%a, handle%smoother, r, z, info)
                if (info%code /= FORTSPARSE_OK) return
            end do
            handle%res = r
            if (handle%use_gradient) call correct_gradient(handle, z, info)
            if (info%code /= FORTSPARSE_OK) return
            call correct_nodal(handle, z, info)
            return
        end if
        do step = 1, handle%steps
            call smoother_apply(handle%a, handle%smoother, r, z, info)
            if (info%code /= FORTSPARSE_OK) return
        end do
        if (handle%use_gradient) then
            call residual(handle, r, z)
            call correct_gradient(handle, z, info)
            if (info%code /= FORTSPARSE_OK) return
        end if
        call residual(handle, r, z)
        call correct_nodal(handle, z, info)
        if (info%code /= FORTSPARSE_OK) return
        if (handle%use_gradient) then
            call residual(handle, r, z)
            call correct_gradient(handle, z, info)
            if (info%code /= FORTSPARSE_OK) return
        end if
        do step = 1, handle%steps
            call smoother_apply(handle%a, handle%smoother, r, z, info)
            if (info%code /= FORTSPARSE_OK) return
        end do
    end subroutine ams_apply

    subroutine ams_destroy(handle)
        type(ams_handle_t), intent(inout) :: handle

        integer :: d

        if (allocated(handle%bg)) then
            call amg_destroy(handle%bg)
            deallocate (handle%bg)
        end if
        if (allocated(handle%bcomp)) then
            do d = 1, size(handle%bcomp)
                call amg_destroy(handle%bcomp(d))
            end do
            deallocate (handle%bcomp)
        end if
        if (allocated(handle%bpi)) then
            call amg_destroy(handle%bpi)
            deallocate (handle%bpi)
        end if
        call smoother_destroy(handle%smoother)
        ! Released component by component, as pcg_destroy does.
        call release(handle%a)
        call release(handle%g)
        call release(handle%gt)
        call release(handle%pi)
        call release(handle%pit)
        call release(handle%ag)
        call release(handle%api)
        do d = 1, 3
            call release(handle%acomp(d))
        end do
        if (allocated(handle%res)) deallocate (handle%res, handle%z0, handle%r0, &
            handle%z0n, handle%rpi, handle%zpi)
        handle%use_gradient = .false.
        handle%components = .false.
        handle%ready = .false.
    end subroutine ams_destroy

    subroutine release(a)
        type(csr_t), intent(inout) :: a
        if (allocated(a%row_ptr)) deallocate (a%row_ptr)
        if (allocated(a%col)) deallocate (a%col)
        if (allocated(a%val)) deallocate (a%val)
        a%nrow = 0
        a%ncol = 0
    end subroutine release

    subroutine residual(handle, r, z)
        type(ams_handle_t), intent(inout) :: handle
        real(dp), intent(in) :: r(:), z(:)

        call csr_matvec(handle%a, z, handle%res)
        handle%res = r - handle%res
    end subroutine residual

    subroutine correct_gradient(handle, z, info)
        !! z += G B_G G^T res.
        type(ams_handle_t), intent(inout) :: handle
        real(dp), intent(inout) :: z(:)
        type(fortsparse_status_t), intent(out) :: info

        call csr_matvec(handle%gt, handle%res, handle%r0)
        call amg_apply(handle%bg, handle%r0, handle%z0n, info)
        if (info%code /= FORTSPARSE_OK) return
        call csr_matvec(handle%g, handle%z0n, handle%z0)
        z = z + handle%z0
    end subroutine correct_gradient

    subroutine correct_nodal(handle, z, info)
        !! z += Pi B_Pi Pi^T res.
        type(ams_handle_t), intent(inout) :: handle
        real(dp), intent(inout) :: z(:)
        type(fortsparse_status_t), intent(out) :: info

        integer :: d, n0

        call csr_matvec(handle%pit, handle%res, handle%rpi)
        if (handle%components) then
            n0 = handle%n0
            do d = 1, 3
                call amg_apply(handle%bcomp(d), handle%rpi((d - 1)*n0 + 1:d*n0), &
                    handle%zpi((d - 1)*n0 + 1:d*n0), info)
                if (info%code /= FORTSPARSE_OK) return
            end do
        else
            call amg_apply(handle%bpi, handle%rpi, handle%zpi, info)
            if (info%code /= FORTSPARSE_OK) return
        end if
        call csr_matvec(handle%pi, handle%zpi, handle%z0)
        z = z + handle%z0
    end subroutine correct_nodal

    subroutine add_into(a, b)
        !! a = a + b on the union of both patterns, rows sorted.
        type(csr_t), intent(inout) :: a
        type(csr_t), intent(in) :: b
        type(csr_t) :: c
        integer :: i, p, q, k

        c%nrow = a%nrow
        c%ncol = a%ncol
        allocate (c%row_ptr(a%nrow + 1), c%col(size(a%col) + size(b%col)), &
            c%val(size(a%col) + size(b%col)))
        k = 0
        c%row_ptr(1) = 1
        do i = 1, a%nrow
            p = a%row_ptr(i)
            q = b%row_ptr(i)
            do while (p < a%row_ptr(i + 1) .or. q < b%row_ptr(i + 1))
                k = k + 1
                if (q >= b%row_ptr(i + 1)) then
                    c%col(k) = a%col(p)
                    c%val(k) = a%val(p)
                    p = p + 1
                else if (p >= a%row_ptr(i + 1)) then
                    c%col(k) = b%col(q)
                    c%val(k) = b%val(q)
                    q = q + 1
                else if (a%col(p) == b%col(q)) then
                    c%col(k) = a%col(p)
                    c%val(k) = a%val(p) + b%val(q)
                    p = p + 1
                    q = q + 1
                else if (a%col(p) < b%col(q)) then
                    c%col(k) = a%col(p)
                    c%val(k) = a%val(p)
                    p = p + 1
                else
                    c%col(k) = b%col(q)
                    c%val(k) = b%val(q)
                    q = q + 1
                end if
            end do
            c%row_ptr(i + 1) = k + 1
        end do
        c%col = c%col(:k)
        c%val = c%val(:k)
        call move_alloc(c%row_ptr, a%row_ptr)
        call move_alloc(c%col, a%col)
        call move_alloc(c%val, a%val)
    end subroutine add_into

    subroutine diagonal_block(b, offset, n, c)
        !! c = b(offset+1:offset+n, offset+1:offset+n), rows kept sorted.
        type(csr_t), intent(in) :: b
        integer, intent(in) :: offset, n
        type(csr_t), intent(out) :: c
        integer :: i, k, m

        c%nrow = n
        c%ncol = n
        allocate (c%row_ptr(n + 1))
        c%row_ptr(1) = 1
        do i = 1, n
            m = 0
            do k = b%row_ptr(offset + i), b%row_ptr(offset + i + 1) - 1
                if (b%col(k) > offset .and. b%col(k) <= offset + n) m = m + 1
            end do
            c%row_ptr(i + 1) = c%row_ptr(i) + m
        end do
        allocate (c%col(c%row_ptr(n + 1) - 1), c%val(c%row_ptr(n + 1) - 1))
        m = 0
        do i = 1, n
            do k = b%row_ptr(offset + i), b%row_ptr(offset + i + 1) - 1
                if (b%col(k) > offset .and. b%col(k) <= offset + n) then
                    m = m + 1
                    c%col(m) = b%col(k) - offset
                    c%val(m) = b%val(k)
                end if
            end do
        end do
    end subroutine diagonal_block

    subroutine regularize(b, shift)
        !! Symmetrizes the Galerkin product's rounding and adds shift times a
        !! diagonal floored at 1e-3 of its mean, so empty rows stay definite.
        type(csr_t), intent(inout) :: b
        real(dp), intent(in) :: shift
        integer :: i, k, j, m
        real(dp) :: floor
        real(dp), allocatable :: d(:)

        allocate (d(b%nrow))
        d = 0
        do i = 1, b%nrow
            do k = b%row_ptr(i), b%row_ptr(i + 1) - 1
                if (b%col(k) == i) d(i) = b%val(k)
            end do
        end do
        floor = 1.0e-3_dp*sum(abs(d))/b%nrow
        do i = 1, b%nrow
            do k = b%row_ptr(i), b%row_ptr(i + 1) - 1
                j = b%col(k)
                if (j > i) then
                    ! The mirrored entry; patterns of Galerkin products are symmetric.
                    do m = b%row_ptr(j), b%row_ptr(j + 1) - 1
                        if (b%col(m) == i) then
                            b%val(k) = 0.5_dp*(b%val(k) + b%val(m))
                            b%val(m) = b%val(k)
                            exit
                        end if
                    end do
                else if (j == i) then
                    b%val(k) = b%val(k) + shift*max(d(i), floor)
                end if
            end do
        end do
    end subroutine regularize
end module fortsparse_ams
