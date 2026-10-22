module fortsparse_amg
    !! Smoothed aggregation for caller-declared real symmetric positive definite CSR.
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csr, only: csr_t, csr_budget_t, csr_valid, csr_find, &
            csr_allocate, csr_copy, csr_release, csr_matvec, csr_transpose, csr_multiply
    use fortsparse_csc, only: csc_t
    use fortsparse_status, only: fortsparse_status_t, status_set, &
                    FORTSPARSE_OK, FORTSPARSE_INVALID_MATRIX, FORTSPARSE_NOT_FACTORED, &
                             FORTSPARSE_RESOURCE_LIMIT, FORTSPARSE_NOT_POSITIVE_DEFINITE
    use fortsparse_inertia, only: inertia_count_t
    use fortsparse_native_ldlt, only: native_ldlt_t
    use fortsparse_smoothers, only: csr_smoother_t, smoother_setup, &
                                  smoother_apply, SMOOTHER_CHEBYSHEV, SMOOTHER_L1_JACOBI
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: amg_params_t, amg_handle_t, amg_setup, amg_apply, amg_destroy
    public :: amg_setup_csr, amg_tentative_prolongator
    type :: amg_params_t
        real(dp) :: strength_threshold = 0.02_dp
        real(dp) :: prolongator_damping = 4.0_dp/3.0_dp
        integer :: coarse_size = 64, max_levels = 20
        integer :: smoother = SMOOTHER_CHEBYSHEV, smoothing_steps = 2
        real(dp) :: chebyshev_ratio = 20
        integer :: cycle = 1 ! 1: V, 2: W
        ! Strict whole-hierarchy cap is not yet supported, positive requests fail.
        integer(i8) :: memory_limit_bytes = 0
        real(dp), allocatable :: near_nullspace(:, :)
    end type
    type :: level_t
        type(csr_t) :: a, p, r
        type(csr_smoother_t) :: smoother
        real(dp), allocatable :: b(:), x(:), residual(:)
    end type
    type :: amg_handle_t
        private
        type(level_t), allocatable :: level(:)
        type(native_ldlt_t) :: coarse
        integer :: count = 0, visits = 1
        logical :: ready = .false.
    contains
        procedure :: levels => level_count
        procedure :: rows => level_rows
        procedure :: operator_complexity
    end type
contains
    subroutine amg_setup_csr(n, row_ptr, columns, values, params, handle, info)
        !! Raw-array bridge for consumers without a FortSparse matrix object.
        integer, intent(in) :: n, row_ptr(:), columns(:)
        real(dp), intent(in) :: values(:)
        type(amg_params_t), intent(in) :: params
        type(amg_handle_t), intent(inout) :: handle
        type(fortsparse_status_t), intent(out) :: info
        type(csr_t) :: a
        if (params%memory_limit_bytes /= 0) then
            call amg_destroy(handle)
            call status_set(info, FORTSPARSE_RESOURCE_LIMIT, &
                            'AMG: strict hierarchy memory cap unsupported')
            return
        end if
        a%nrow = n
        a%ncol = n
        a%row_ptr = row_ptr
        a%col = columns
        a%val = values
        call amg_setup(a, params, handle, info)
    end subroutine

    subroutine amg_setup(a, params, handle, info)
        type(csr_t), intent(in) :: a
        type(amg_params_t), intent(in) :: params
        type(amg_handle_t), intent(inout) :: handle
        type(fortsparse_status_t), intent(out) :: info
        type(csr_budget_t) :: budget
        type(csr_t) :: tentative, product, temporary
        type(csc_t) :: csc
        type(inertia_count_t) :: inertia
        real(dp), allocatable :: basis(:, :), next_basis(:, :), diagonal(:)
        real(dp) :: bound, scale, asymmetry
        integer :: l, n, i, j, k, reverse, invalid, inertia_state
        call amg_destroy(handle)
        call status_set(info, FORTSPARSE_INVALID_MATRIX, 'AMG: CSR/parameters')
        if (.not. csr_valid(a)) return
        if (a%nrow /= a%ncol) return
        if (.not. ieee_is_finite(params%strength_threshold)) return
        if (.not. ieee_is_finite(params%prolongator_damping)) return
        if (.not. ieee_is_finite(params%chebyshev_ratio)) return
        if (params%strength_threshold < 0 .or. params%strength_threshold > 1) return
        if (params%prolongator_damping <= 0 .or. &
            params%prolongator_damping >= 2) return
        if (params%coarse_size < 1 .or. params%max_levels < 1) return
        if (params%smoothing_steps < 1 .or. params%chebyshev_ratio <= 1) return
        if (params%smoother /= SMOOTHER_L1_JACOBI .and. &
            params%smoother /= SMOOTHER_CHEBYSHEV) return
        if (params%cycle < 1 .or. params%cycle > 2) return
        if (params%memory_limit_bytes /= 0) then
            call status_set(info, FORTSPARSE_RESOURCE_LIMIT, &
                            'AMG: strict hierarchy memory cap unsupported')
            return
        end if
        n = a%nrow
        invalid = 0
        !$omp parallel do private(k,reverse,asymmetry,scale) reduction(max:invalid)
        do i = 1, n
            k = csr_find(a, i, i)
            if (k == 0) then
                invalid = 1
            else if (a%val(k) <= 0) then
                invalid = 1
            end if
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                reverse = csr_find(a, a%col(k), i)
                if (reverse == 0) then
                    invalid = 1
                else
                    scale = max(abs(a%val(k)), abs(a%val(reverse)), tiny(1.0_dp))
                    asymmetry = abs(a%val(k) - a%val(reverse))/scale
                    if (asymmetry > 64*epsilon(1.0_dp)) invalid = 1
                end if
            end do
        end do
        if (invalid /= 0) then
            call status_set(info, FORTSPARSE_INVALID_MATRIX, &
                            'AMG: symmetric storage and positive diagonal required')
            return
        end if
        if (allocated(params%near_nullspace)) then
            if (size(params%near_nullspace, 1) /= n) return
            if (size(params%near_nullspace, 2) < 1) return
            if (.not. all(ieee_is_finite(params%near_nullspace))) return
            basis = params%near_nullspace
        else
            allocate (basis(n, 1))
            basis = 1
        end if
        allocate (handle%level(params%max_levels))
        handle%visits = params%cycle
        call csr_copy(a, handle%level(1)%a, budget, info)
        if (info%code /= 0) goto 900
        do l = 1, params%max_levels
            handle%count = l
            n = handle%level(l)%a%nrow
            allocate (handle%level(l)%b(n), handle%level(l)%x(n), &
                      handle%level(l)%residual(n))
            if (n <= params%coarse_size .or. l == params%max_levels) exit
            call amg_tentative_prolongator(handle%level(l)%a, &
                          params%strength_threshold, basis, tentative, next_basis, info)
            if (info%code /= 0) goto 900
            if (tentative%ncol >= n) exit
            ! Jacobi smoothing with a safe bound on rho(D^-1 A).
            allocate (diagonal(n))
            bound = 0
            !$omp parallel do private(k,j,scale) reduction(max:bound)
            do i = 1, n
                k = csr_find(handle%level(l)%a, i, i)
                diagonal(i) = handle%level(l)%a%val(k)
                scale = 0
                do j = handle%level(l)%a%row_ptr(i), &
                    handle%level(l)%a%row_ptr(i + 1) - 1
                    scale = scale + abs(handle%level(l)%a%val(j))
                end do
                bound = max(bound, scale/diagonal(i))
            end do
            call csr_multiply(handle%level(l)%a, tentative, product, budget, info)
            if (info%code /= 0) goto 900
            !$omp parallel do private(k,j)
            do i = 1, n
                do k = product%row_ptr(i), product%row_ptr(i + 1) - 1
                    product%val(k) = -params%prolongator_damping/bound &
                                     *product%val(k)/diagonal(i)
                    j = csr_find(tentative, i, product%col(k))
                    if (j /= 0) product%val(k) = product%val(k) + tentative%val(j)
                end do
            end do
            deallocate (diagonal)
            call csr_release(tentative)
            call move_csr(product, handle%level(l)%p)
            call csr_transpose(handle%level(l)%p, handle%level(l)%r, budget, info)
            if (info%code /= 0) goto 900
            call csr_multiply(handle%level(l)%a, handle%level(l)%p, &
                              temporary, budget, info)
            if (info%code /= 0) goto 900
            call csr_multiply(handle%level(l)%r, temporary, &
                              handle%level(l + 1)%a, budget, info)
            if (info%code /= 0) goto 900
            call csr_release(temporary, budget)
            ! Average roundoff in RAP with its transpose, keeping a symmetric PC.
            call symmetrize(handle%level(l + 1)%a)
            call smoother_setup(handle%level(l)%a, params%smoother, &
                               handle%level(l)%smoother, info, params%smoothing_steps, &
                                params%chebyshev_ratio)
            if (info%code /= 0) goto 900
            call move_alloc(next_basis, basis)
        end do
        ! Symmetric CSR is also its CSC transpose, so no numeric conversion.
        csc%nrow = n
        csc%ncol = n
        csc%nnz = size(handle%level(l)%a%val)
        csc%col_ptr = handle%level(l)%a%row_ptr
        csc%row_idx = handle%level(l)%a%col
        csc%val = handle%level(l)%a%val
        call handle%coarse%factor(csc, info)
        if (info%code /= 0) goto 900
        call handle%coarse%inertia(inertia, inertia_state)
        if (inertia%positive /= n) then
            call status_set(info, FORTSPARSE_NOT_POSITIVE_DEFINITE, &
                            'AMG: coarse operator is not positive definite')
            goto 900
        end if
        handle%ready = .true.
        return
900     call amg_destroy(handle)
    end subroutine

    subroutine move_csr(a, b)
        type(csr_t), intent(inout) :: a
        type(csr_t), intent(out) :: b
        b%nrow = a%nrow
        b%ncol = a%ncol
        call move_alloc(a%row_ptr, b%row_ptr)
        call move_alloc(a%col, b%col)
        call move_alloc(a%val, b%val)
        a%nrow = 0
        a%ncol = 0
    end subroutine

    subroutine symmetrize(a)
        type(csr_t), intent(inout) :: a
        integer :: i, k, q
        !$omp parallel do private(k,q)
        do i = 1, a%nrow
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                if (a%col(k) <= i) cycle
                q = csr_find(a, a%col(k), i)
                a%val(k) = 0.5_dp*(a%val(k) + a%val(q))
                a%val(q) = a%val(k)
            end do
        end do
    end subroutine

    subroutine amg_tentative_prolongator(a, threshold, basis, t, coarse_basis, info)
        !! Greedy strong-graph aggregates, local rank-revealing modified GS.
        !! The exposed setup primitive also permits independent T Bc = B tests.
        type(csr_t), intent(in) :: a
        real(dp), intent(in) :: threshold, basis(:, :)
        type(csr_t), intent(out) :: t
        real(dp), allocatable, intent(out) :: coarse_basis(:, :)
        type(fortsparse_status_t), intent(out) :: info
        type(csr_budget_t) :: budget
        integer, allocatable :: aggregate(:), ptr(:), members(:), next(:), &
                                rank(:), offset(:)
        logical, allocatable :: strong(:)
        real(dp), allocatable :: diagonal(:), q(:, :)
        real(dp) :: dot, norm, original
        integer :: n, nv, i, j, k, g, ng, v, u, r, s, row, nz
        call status_set(info, FORTSPARSE_INVALID_MATRIX, 'AMG: tentative inputs')
        if (.not. csr_valid(a)) return
        if (a%nrow /= a%ncol .or. size(basis, 1) /= a%nrow) return
        if (size(basis, 2) < 1) return
        if (.not. all(ieee_is_finite(basis))) return
        if (.not. ieee_is_finite(threshold)) return
        if (threshold < 0 .or. threshold > 1) return
        n = a%nrow
        nv = size(basis, 2)
        allocate (diagonal(n), strong(size(a%val)), aggregate(n), q(n, nv))
        do i = 1, n
            k = csr_find(a, i, i)
            if (k == 0) return
            diagonal(i) = a%val(k)
        end do
        if (any(diagonal <= 0)) return
        !$omp parallel do private(k,j)
        do i = 1, n
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                j = a%col(k)
                strong(k) = j /= i .and. abs(a%val(k)) > tiny(1.0_dp)
                strong(k) = strong(k) .and. abs(a%val(k))/sqrt(diagonal(i)) &
                            /sqrt(diagonal(j)) >= threshold
            end do
        end do
        aggregate = 0
        ng = 0
        ! First seed only where all strong neighbors are unassigned.
        do i = 1, n
            if (aggregate(i) /= 0) cycle
            s = 0
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                if (.not. strong(k)) cycle
                if (aggregate(a%col(k)) /= 0) s = 1
            end do
            if (s /= 0) cycle
            ng = ng + 1
            aggregate(i) = ng
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                if (strong(k)) aggregate(a%col(k)) = ng
            end do
        end do
        ! Attach leftover vertices, or seed an isolated remaining component.
        do i = 1, n
            if (aggregate(i) /= 0) cycle
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                if (.not. strong(k)) cycle
                g = aggregate(a%col(k))
                if (g == 0) cycle
                aggregate(i) = g
                exit
            end do
            if (aggregate(i) == 0) then
                ng = ng + 1
                aggregate(i) = ng
            end if
        end do
        allocate (ptr(ng + 1), next(ng), members(n), rank(ng), offset(ng + 1))
        next = 0
        do i = 1, n
            next(aggregate(i)) = next(aggregate(i)) + 1
        end do
        ptr(1) = 1
        do g = 1, ng
            ptr(g + 1) = ptr(g) + next(g)
        end do
        next = ptr(:ng)
        do i = 1, n
            g = aggregate(i)
            members(next(g)) = i
            next(g) = next(g) + 1
        end do
        q = 0
        rank = 0
        !$omp parallel do private(v,r,k,i,j,norm,original,u,dot) schedule(dynamic,32)
        do g = 1, ng
            r = 0
            do v = 1, nv
                original = 0
                do k = ptr(g), ptr(g + 1) - 1
                    i = members(k)
                    q(i, r + 1) = basis(i, v)
                    original = original + basis(i, v)**2
                end do
                ! Reorthogonalize to retain local rank on nearly dependent modes.
                do u = 1, 2
                    do j = 1, r
                        dot = 0
                        do k = ptr(g), ptr(g + 1) - 1
                            i = members(k)
                            dot = dot + q(i, j)*q(i, r + 1)
                        end do
                        do k = ptr(g), ptr(g + 1) - 1
                            i = members(k)
                            q(i, r + 1) = q(i, r + 1) - dot*q(i, j)
                        end do
                    end do
                end do
                norm = 0
                do k = ptr(g), ptr(g + 1) - 1
                    norm = norm + q(members(k), r + 1)**2
                end do
                if (norm <= 1.e-24_dp*original .or. norm <= tiny(1.0_dp)) cycle
                r = r + 1
                do k = ptr(g), ptr(g + 1) - 1
                    q(members(k), r) = q(members(k), r)/sqrt(norm)
                end do
                if (r == nv) exit
            end do
            if (r == 0) then
                ! Zero local modes still need a nonempty aggregate basis.
                r = 1
                do k = ptr(g), ptr(g + 1) - 1
                    q(members(k), 1) = 1/sqrt(real(ptr(g + 1) - ptr(g), dp))
                end do
            end if
            rank (g) = r
        end do
        offset(1) = 1
        do g = 1, ng
            offset(g + 1) = offset(g) + rank(g)
        end do
        allocate (coarse_basis(offset(ng + 1) - 1, nv))
        coarse_basis = 0
        !$omp parallel do private(v,r,k,i,dot)
        do g = 1, ng
            do v = 1, nv
                do r = 1, rank(g)
                    dot = 0
                    do k = ptr(g), ptr(g + 1) - 1
                        i = members(k)
                        dot = dot + q(i, r)*basis(i, v)
                    end do
                    coarse_basis(offset(g) + r - 1, v) = dot
                end do
            end do
        end do
        nz = sum(rank(aggregate))
        call csr_allocate(t, n, size(coarse_basis, 1), int(nz, i8), budget, info)
        if (info%code /= 0) return
        t%row_ptr(1) = 1
        do i = 1, n
            t%row_ptr(i + 1) = t%row_ptr(i) + rank(aggregate(i))
        end do
        !$omp parallel do private(g,r,row)
        do i = 1, n
            g = aggregate(i)
            do r = 1, rank(g)
                row = t%row_ptr(i) + r - 1
                t%col(row) = offset(g) + r - 1
                t%val(row) = q(i, r)
            end do
        end do
    end subroutine

    subroutine amg_apply(handle, r, z, info)
        !! Fixed linear preconditioner, no allocation and no procedure pointers.
        !! One handle owns mutable workspace: concurrent callers need separate handles.
        type(amg_handle_t), intent(inout) :: handle
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        type(fortsparse_status_t), intent(out), optional :: info
        type(fortsparse_status_t) :: status
        integer :: i
        z = 0
        call status_set(status, FORTSPARSE_NOT_FACTORED, 'AMG: setup required')
        if (handle%ready) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, 'AMG: vector shape')
            if (size(r) == handle%level(1)%a%nrow .and. size(z) == size(r)) then
                !$omp parallel do if(size(r) >= 512)
                do i = 1, size(r)
                    handle%level(1)%b(i) = r(i)
                    handle%level(1)%x(i) = 0
                end do
                call cycle_level(handle, 1, status)
                if (status%code == 0) z = handle%level(1)%x
            end if
        end if
        if (present(info)) info = status
    end subroutine

    recursive subroutine cycle_level(handle, l, info)
        type(amg_handle_t), intent(inout) :: handle
        integer, intent(in) :: l
        type(fortsparse_status_t), intent(out) :: info
        integer :: visit, i
        if (l == handle%count) then
            call handle%coarse%solve(handle%level(l)%b, handle%level(l)%x, info)
            return
        end if
        call smoother_apply(handle%level(l)%a, handle%level(l)%smoother, &
                            handle%level(l)%b, handle%level(l)%x, info)
        if (info%code /= 0) return
        call csr_matvec(handle%level(l)%a, handle%level(l)%x, &
                        handle%level(l)%residual)
        !$omp parallel do if(size(handle%level(l)%b) >= 512)
        do i = 1, size(handle%level(l)%b)
            handle%level(l)%residual(i) = handle%level(l)%b(i) &
                                          - handle%level(l)%residual(i)
        end do
        call csr_matvec(handle%level(l)%r, handle%level(l)%residual, &
                        handle%level(l + 1)%b)
        handle%level(l + 1)%x = 0
        do visit = 1, handle%visits
            call cycle_level(handle, l + 1, info)
            if (info%code /= 0) return
        end do
        call csr_matvec(handle%level(l)%p, handle%level(l + 1)%x, &
                        handle%level(l)%residual)
        !$omp parallel do if(size(handle%level(l)%b) >= 512)
        do i = 1, size(handle%level(l)%b)
            handle%level(l)%x(i) = handle%level(l)%x(i) + handle%level(l)%residual(i)
        end do
        call smoother_apply(handle%level(l)%a, handle%level(l)%smoother, &
                            handle%level(l)%b, handle%level(l)%x, info)
    end subroutine

    subroutine amg_destroy(handle)
        type(amg_handle_t), intent(inout) :: handle
        if (allocated(handle%level)) deallocate (handle%level)
        call handle%coarse%free()
        handle%count = 0
        handle%ready = .false.
    end subroutine

    integer function level_count(self) result(n)
        class(amg_handle_t), intent(in) :: self
        n = self%count
    end function

    integer function level_rows(self, level) result(n)
        class(amg_handle_t), intent(in) :: self
        integer, intent(in) :: level
        n = 0
        if (level < 1 .or. level > self%count) return
        n = self%level(level)%a%nrow
    end function

    real(dp) function operator_complexity(self) result(value)
        class(amg_handle_t), intent(in) :: self
        integer :: l
        value = 0
        if (self%count == 0) return
        do l = 1, self%count
            value = value + real(size(self%level(l)%a%val), dp)
        end do
        value = value/size(self%level(1)%a%val)
    end function
end module
