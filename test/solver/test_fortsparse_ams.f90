program test_fortsparse_ams
    ! AMS-preconditioned CG on lowest-order Nedelec edges of a structured cube
    ! grid: A = C^T C + beta I with C the face-edge incidence (curl), G the
    ! edge-node incidence (gradient) and Pi the nodal vector interpolation
    ! (half the edge tangent at each end node). Oracles: C G = 0, the true
    ! residual, and iteration counts against IC(0) under refinement, for the
    ! definite (beta = 1) and the singular, consistent (beta = 0) problem.
    use fortsparse, only: dp, csr_t, csr_matvec, fortsparse_status_t, pcg_params_t, &
        pcg_handle_t, pcg_setup, pcg_setup_ams, pcg_solve, pcg_destroy, PCG_IC0, &
        FORTSPARSE_OK, AMS_MULTIPLICATIVE, AMS_ADDITIVE
    use amg_test_support, only: check
    implicit none
    type(csr_t) :: a, g, pi, c
    type(pcg_params_t) :: params
    type(pcg_handle_t) :: handle
    type(fortsparse_status_t) :: info
    real(dp), allocatable :: x(:), b(:), exact(:), r(:), cg(:)
    integer :: level, side, beta_case, variant, it_ams(2, 2), it_ic0(2, 2), i, it
    ! Grid of the current Nedelec problem, for node, edge and face.
    integer :: np, ne(3), cur_side
    real(dp) :: beta, relative

    params%rtol = 1.0e-8_dp
    params%max_iterations = 2000
    do beta_case = 1, 2
        beta = merge(1.0_dp, 0.0_dp, beta_case == 1)
        do level = 1, 2
            side = 4*level
            call nedelec(side, beta, a, g, pi, c)
            if (level == 1 .and. beta_case == 1) then
                allocate (cg(c%nrow))
                call csr_matvec(c, g_times(g, [(real(i, dp)**2, i=1, g%ncol)]), cg)
                call check(maxval(abs(cg)) < 1.0e-12_dp, 'curl grad = 0')
                deallocate (cg)
            end if
            allocate (x(a%nrow), b(a%nrow), exact(a%nrow), r(a%nrow))
            do i = 1, a%nrow
                exact(i) = sin(0.31_dp*i) + cos(0.17_dp*i)
            end do
            call csr_matvec(a, exact, b)
            do variant = 1, 2
                params%ams%cycle = merge(AMS_MULTIPLICATIVE, AMS_ADDITIVE, variant == 1)
                call pcg_setup_ams(a, g, pi, params, handle, info)
                call check(info%code == FORTSPARSE_OK, 'AMS setup: '//trim(info%msg))
                x = 0
                call pcg_solve(handle, params, b, x, info, it, relative)
                call check(info%code == FORTSPARSE_OK, 'AMS-PCG solve: '//trim(info%msg))
                call csr_matvec(a, x, r)
                call check(norm2(b - r) <= 1.0e-8_dp*norm2(b)*(1 + 1.0e-6_dp), &
                    'AMS-PCG true residual')
                if (variant == 1) it_ams(beta_case, level) = it
                call pcg_destroy(handle)
            end do
            params%preconditioner = PCG_IC0
            call pcg_setup(a, params, handle, info)
            x = 0
            call pcg_solve(handle, params, b, x, info, it_ic0(beta_case, level))
            call check(info%code == FORTSPARSE_OK, 'IC(0) solve')
            call pcg_destroy(handle)
            deallocate (x, b, exact, r)
        end do
        print '(a,f3.1,a,2i5,a,2i5)', 'beta ', beta, ': AMS iterations', &
            it_ams(beta_case, :), ', IC(0)', it_ic0(beta_case, :)
        call check(it_ams(beta_case, 2) < it_ic0(beta_case, 2), 'AMS beats IC(0)')
        call check(it_ams(beta_case, 2) <= it_ams(beta_case, 1) + 8, &
            'AMS iterations stay nearly flat under refinement')
    end do
    ! AMS handles local to a routine, destroyed or not, are released when it
    ! returns.
    call nedelec(4, 1.0_dp, a, g, pi, c)
    do i = 1, 3
        call scoped_ams(i == 2)
    end do
    print '(a)', 'PASS'

contains

    subroutine scoped_ams(destroy)
        logical, intent(in) :: destroy
        type(pcg_handle_t) :: local
        real(dp), allocatable :: xs(:), bs(:)
        allocate (xs(a%nrow), bs(a%nrow))
        bs = 1
        xs = 0
        call pcg_setup_ams(a, g, pi, params, local, info)
        call check(info%code == FORTSPARSE_OK, 'scoped AMS setup')
        call pcg_solve(local, params, bs, xs, info)
        call check(info%code == FORTSPARSE_OK, 'scoped AMS solve')
        if (destroy) call pcg_destroy(local)
    end subroutine scoped_ams

    function g_times(g, v) result(y)
        type(csr_t), intent(in) :: g
        real(dp), intent(in) :: v(:)
        real(dp), allocatable :: y(:)
        allocate (y(g%nrow))
        call csr_matvec(g, v, y)
    end function g_times

    subroutine nedelec(side, beta, a, g, pi, c)
        ! side^3 cells, unit spacing; edges in x, y, z families.
        integer, intent(in) :: side
        real(dp), intent(in) :: beta
        type(csr_t), intent(out) :: a, g, pi, c
        integer :: nf(3), n0, n1, nfaces, ix, iy, iz, d, e, f, k, m
        integer, allocatable :: gi(:, :)
        real(dp), allocatable :: dense_ctc(:)
        integer :: o(3)
        cur_side = side
        np = side + 1
        n0 = np**3
        ! Edges of direction d at nodes with index < side along d.
        ne = side*np*np
        n1 = sum(ne)
        allocate (gi(n1, 2))
        e = 0
        do d = 1, 3
            do iz = 0, merge(side - 1, side, d == 3)
                do iy = 0, merge(side - 1, side, d == 2)
                    do ix = 0, merge(side - 1, side, d == 1)
                        e = e + 1
                        o = 0
                        o(d) = 1
                        gi(e, 1) = node(ix, iy, iz)
                        gi(e, 2) = node(ix + o(1), iy + o(2), iz + o(3))
                    end do
                end do
            end do
        end do
        ! G: -1 at the start node, +1 at the end node (columns sorted).
        g%nrow = n1
        g%ncol = n0
        allocate (g%row_ptr(n1 + 1), g%col(2*n1), g%val(2*n1))
        do e = 1, n1
            g%row_ptr(e) = 2*e - 1
            g%col(2*e - 1) = gi(e, 1)
            g%val(2*e - 1) = -1
            g%col(2*e) = gi(e, 2)
            g%val(2*e) = 1
        end do
        g%row_ptr(n1 + 1) = 2*n1 + 1
        ! Pi: edge e of direction d takes 1/2 of the d-component at each end.
        pi%nrow = n1
        pi%ncol = 3*n0
        allocate (pi%row_ptr(n1 + 1), pi%col(2*n1), pi%val(2*n1))
        e = 0
        do d = 1, 3
            do k = 1, ne(d)
                e = e + 1
                pi%row_ptr(e) = 2*e - 1
                pi%col(2*e - 1) = (d - 1)*n0 + gi(e, 1)
                pi%col(2*e) = (d - 1)*n0 + gi(e, 2)
                pi%val(2*e - 1:2*e) = 0.5_dp
            end do
        end do
        pi%row_ptr(n1 + 1) = 2*n1 + 1
        ! C: each face, normal d, circulates its four edges.
        nf = side*side*np
        nfaces = sum(nf)
        c%nrow = nfaces
        c%ncol = n1
        allocate (c%row_ptr(nfaces + 1), c%col(4*nfaces), c%val(4*nfaces))
        f = 0
        m = 0
        do d = 1, 3
            do iz = 0, merge(side, side - 1, d == 3)
                do iy = 0, merge(side, side - 1, d == 2)
                    do ix = 0, merge(side, side - 1, d == 1)
                        f = f + 1
                        c%row_ptr(f) = m + 1
                        call face(d, ix, iy, iz, c%col(m + 1:m + 4), c%val(m + 1:m + 4))
                        call sort4(c%col(m + 1:m + 4), c%val(m + 1:m + 4))
                        m = m + 4
                    end do
                end do
            end do
        end do
        c%row_ptr(nfaces + 1) = m + 1
        ! A = C^T C + beta I, dense accumulation per row (small grids only).
        a%nrow = n1
        a%ncol = n1
        allocate (a%row_ptr(n1 + 1), dense_ctc(n1))
        block
            integer, allocatable :: cols(:)
            real(dp), allocatable :: vals(:)
            integer :: r, p, q, nnz
            type(csr_t) :: ct
            ! A curl-curl row of lowest-order edges couples at most 33 edges.
            allocate (cols(33*n1), vals(33*n1))
            call transpose_csr(c, ct)
            nnz = 0
            a%row_ptr(1) = 1
            do r = 1, n1
                dense_ctc = 0
                do p = ct%row_ptr(r), ct%row_ptr(r + 1) - 1
                    f = ct%col(p)
                    do q = c%row_ptr(f), c%row_ptr(f + 1) - 1
                        dense_ctc(c%col(q)) = dense_ctc(c%col(q)) + ct%val(p)*c%val(q)
                    end do
                end do
                dense_ctc(r) = dense_ctc(r) + beta
                do q = 1, n1
                    if (dense_ctc(q) /= 0 .or. q == r) then
                        nnz = nnz + 1
                        cols(nnz) = q
                        vals(nnz) = dense_ctc(q)
                    end if
                end do
                a%row_ptr(r + 1) = nnz + 1
            end do
            a%col = cols(:nnz)
            a%val = vals(:nnz)
        end block
    end subroutine nedelec

    integer function node(ix, iy, iz)
        integer, intent(in) :: ix, iy, iz
        node = 1 + ix + np*(iy + np*iz)
    end function node
    integer function edge(d, ix, iy, iz)
        ! Edge of direction d starting at node (ix, iy, iz).
        integer, intent(in) :: d, ix, iy, iz
        integer :: lx, ly
        lx = merge(cur_side, np, d == 1)
        ly = merge(cur_side, np, d == 2)
        edge = sum(ne(:d - 1)) + 1 + ix + lx*(iy + ly*iz)
    end function edge
    subroutine face(d, ix, iy, iz, cols, vals)
        integer, intent(in) :: d, ix, iy, iz
        integer, intent(out) :: cols(4)
        real(dp), intent(out) :: vals(4)
        integer :: a1, a2, o1(3), o2(3)
        a1 = 1 + mod(d, 3)
        a2 = 1 + mod(d + 1, 3)
        o1 = 0
        o1(a1) = 1
        o2 = 0
        o2(a2) = 1
        ! Boundary loop of the face spanned by a1, a2: +a1, +a2 at the
        ! end of a1, -a1 at the end of a2, -a2.
        cols(1) = edge(a1, ix, iy, iz)
        cols(2) = edge(a2, ix + o1(1), iy + o1(2), iz + o1(3))
        cols(3) = edge(a1, ix + o2(1), iy + o2(2), iz + o2(3))
        cols(4) = edge(a2, ix, iy, iz)
        vals = [1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp]
    end subroutine face

    subroutine sort4(cols, vals)
        integer, intent(inout) :: cols(4)
        real(dp), intent(inout) :: vals(4)
        integer :: i, j, tc
        real(dp) :: tv
        do i = 2, 4
            tc = cols(i)
            tv = vals(i)
            j = i - 1
            do while (j >= 1)
                if (cols(j) <= tc) exit
                cols(j + 1) = cols(j)
                vals(j + 1) = vals(j)
                j = j - 1
            end do
            cols(j + 1) = tc
            vals(j + 1) = tv
        end do
    end subroutine sort4

    subroutine transpose_csr(a, t)
        type(csr_t), intent(in) :: a
        type(csr_t), intent(out) :: t
        integer :: i, k, p
        integer, allocatable :: next(:)
        t%nrow = a%ncol
        t%ncol = a%nrow
        allocate (t%row_ptr(a%ncol + 1), t%col(size(a%col)), t%val(size(a%col)), next(a%ncol))
        t%row_ptr = 0
        do k = 1, size(a%col)
            t%row_ptr(a%col(k) + 1) = t%row_ptr(a%col(k) + 1) + 1
        end do
        t%row_ptr(1) = 1
        do i = 1, a%ncol
            t%row_ptr(i + 1) = t%row_ptr(i + 1) + t%row_ptr(i)
        end do
        next = t%row_ptr(1:a%ncol)
        do i = 1, a%nrow
            do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
                p = next(a%col(k))
                t%col(p) = i
                t%val(p) = a%val(k)
                next(a%col(k)) = p + 1
            end do
        end do
    end subroutine transpose_csr
end program test_fortsparse_ams
