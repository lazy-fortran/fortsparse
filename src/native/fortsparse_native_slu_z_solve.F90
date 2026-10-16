submodule (fortsparse_native_slu) fortsparse_native_slu_z_solve
    ! Complex solves with the supernodal LU factors: operation 0 solves
    ! A x = b, 1 solves A^T x = b and 2 solves A^H x = b. The scales are real,
    ! so A^-H b = conj(A^-T conj(b)): the adjoint conjugates the right-hand
    ! side before and the solution after the transpose solve, which is exact
    ! and needs no conjugated copy of the factors. Allocation free.
    use fortsparse_slu_z_leaves, only: slu_z_scale, slu_z_apply_swaps, &
        slu_z_block_lower, slu_z_block_upper, slu_z_block_upper_t, &
        slu_z_block_lower_t, slu_z_scatter_update, slu_z_gather_update
    implicit none
contains
    module subroutine solve_complex(self, x, status, operation)
        class(slu_factor_t), intent(inout) :: self
        complex(dp), intent(inout) :: x(:)
        type(fortsparse_status_t), intent(out) :: status
        integer, optional, intent(in) :: operation
        integer :: op, k
        op = 0
        if (present(operation)) op = operation
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            'supernodal LU: invalid operation')
        if (op < 0 .or. op > 2) return
        call status_set(status, FORTSPARSE_NOT_FACTORED, 'supernodal LU: no factors')
        if (.not. self%factored) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            'supernodal LU: real factors need a real RHS')
        if (.not. self%complex_values) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            'supernodal LU: RHS dimension mismatch')
        if (size(x) /= self%n) return
        if (op == 0) then
            do k = 1, self%n
                self%zsolve_work(k) = slu_z_scale(self%drp(k), x(self%row_old(k)), &
                    1.0_dp)
            end do
            call forward_normal_z(self)
            call backward_normal_z(self)
            do k = 1, self%n
                x(self%col_old(k)) = slu_z_scale(self%dcp(k), self%zsolve_work(k), &
                    1.0_dp)
            end do
        else
            if (op == 2) x = conjg(x)
            do k = 1, self%n
                self%zsolve_work(k) = slu_z_scale(self%dcp(k), x(self%col_old(k)), &
                    1.0_dp)
            end do
            call forward_transpose_z(self)
            call backward_transpose_z(self)
            do k = 1, self%n
                x(self%row_old(k)) = slu_z_scale(self%drp(k), self%zsolve_work(k), &
                    1.0_dp)
            end do
            if (op == 2) x = conjg(x)
        end if
        self%zsolve_work(1:self%n) = (0.0_dp, 0.0_dp)
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine solve_complex

    module subroutine solve_block_complex(self, b, status, operation)
        class(slu_factor_t), intent(inout) :: self
        complex(dp), intent(inout) :: b(:,:)
        type(fortsparse_status_t), intent(out) :: status
        integer, optional, intent(in) :: operation
        integer :: j, op
        op = 0
        if (present(operation)) op = operation
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            'supernodal LU: invalid operation')
        if (op < 0 .or. op > 2) return
        call status_set(status, FORTSPARSE_NOT_FACTORED, 'supernodal LU: no factors')
        if (.not. self%factored) return
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            'supernodal LU: RHS dimension mismatch')
        if (size(b,1) /= self%n) return
        call status_set(status, FORTSPARSE_OK, '')
        do j = 1, size(b,2)
            call self%solve(b(:,j), status, op)
            if (.not. status_ok(status)) return
        end do
    end subroutine solve_block_complex

    subroutine forward_normal_z(self)
        ! Per supernode: P_d on the accumulated block, unit-lower block
        ! solve, then subtract L21 times the block from the later rows.
        class(slu_factor_t), intent(inout) :: self
        integer :: s, j, m, w, first
        integer(i8) :: l0
        complex(dp) :: value
        do s = 1, self%node_count
            m = self%panels(s)%m
            w = self%panels(s)%width
            first = self%panels(s)%first
            l0 = self%panels(s)%offset_l
            call slu_z_apply_swaps(self%zsolve_work(first), self%pivots(first), w, &
                .false.)
            call slu_z_block_lower(self%zvalues(l0+1), m, w, self%zsolve_work(first))
            if (m == w) cycle
            do j = 1, w
                value = self%zsolve_work(first+j-1)
                call slu_z_scatter_update(self%zsolve_work, self%panels(s)%rows(w+1:), &
                    self%zvalues(l0+int(j-1,i8)*m+w+1), m-w, value)
            end do
        end do
    end subroutine forward_normal_z

    subroutine backward_normal_z(self)
        class(slu_factor_t), intent(inout) :: self
        integer :: s, c, m, w, first
        integer(i8) :: l0, u0
        complex(dp) :: value
        do s = self%node_count, 1, -1
            m = self%panels(s)%m
            w = self%panels(s)%width
            first = self%panels(s)%first
            l0 = self%panels(s)%offset_l
            u0 = self%panels(s)%offset_u
            if (m > w) then
                do c = 1, w
                    value = self%zsolve_work(first+c-1)
                    call slu_z_gather_update(value, self%panels(s)%rows(w+1:), &
                        self%zvalues(u0+int(c-1,i8)*(m-w)+1), self%zsolve_work, m-w)
                    self%zsolve_work(first+c-1) = value
                end do
            end if
            call slu_z_block_upper(self%zvalues(l0+1), m, w, self%zsolve_work(first))
        end do
    end subroutine backward_normal_z

    subroutine forward_transpose_z(self)
        ! U^T is lower triangular: block solve with U_dd^T, then subtract
        ! Ut times the block from the later rows.
        class(slu_factor_t), intent(inout) :: self
        integer :: s, c, m, w, first
        integer(i8) :: l0, u0
        complex(dp) :: value
        do s = 1, self%node_count
            m = self%panels(s)%m
            w = self%panels(s)%width
            first = self%panels(s)%first
            l0 = self%panels(s)%offset_l
            u0 = self%panels(s)%offset_u
            call slu_z_block_upper_t(self%zvalues(l0+1), m, w, self%zsolve_work(first))
            if (m == w) cycle
            do c = 1, w
                value = self%zsolve_work(first+c-1)
                call slu_z_scatter_update(self%zsolve_work, self%panels(s)%rows(w+1:), &
                    self%zvalues(u0+int(c-1,i8)*(m-w)+1), m-w, value)
            end do
        end do
    end subroutine forward_transpose_z

    subroutine backward_transpose_z(self)
        ! Reverse transposed forward operator: gather L21^T, unit-upper block
        ! solve with L_dd^T, then P_d^T.
        class(slu_factor_t), intent(inout) :: self
        integer :: s, j, m, w, first
        integer(i8) :: l0
        complex(dp) :: value
        do s = self%node_count, 1, -1
            m = self%panels(s)%m
            w = self%panels(s)%width
            first = self%panels(s)%first
            l0 = self%panels(s)%offset_l
            if (m > w) then
                do j = 1, w
                    value = self%zsolve_work(first+j-1)
                    call slu_z_gather_update(value, self%panels(s)%rows(w+1:), &
                        self%zvalues(l0+int(j-1,i8)*m+w+1), self%zsolve_work, m-w)
                    self%zsolve_work(first+j-1) = value
                end do
            end if
            call slu_z_block_lower_t(self%zvalues(l0+1), m, w, self%zsolve_work(first))
            call slu_z_apply_swaps(self%zsolve_work(first), self%pivots(first), w, &
                .true.)
        end do
    end subroutine backward_transpose_z
end submodule fortsparse_native_slu_z_solve
