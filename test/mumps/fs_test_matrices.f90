module fs_test_matrices
    ! SPD finite-element test matrices: the 5-point 2D Laplacian and the P3
    ! (cubic Lagrange) Poisson stiffness matrix on a right-triangle mesh of the
    ! unit square, both with Dirichlet rows replaced by identity.
    use fortsparse, only: dp, csc_t, csc_from_triplet, fortsparse_status_t, &
        status_ok
    implicit none
    private

    public :: laplace2d, p3_poisson, inf_residual

contains

    ! 5-point Laplacian on an m-by-m interior grid (n = m*m), scaled to unity
    ! grid spacing.
    subroutine laplace2d(m, A)
        integer,     intent(in)  :: m
        type(csc_t), intent(out) :: A

        integer, allocatable  :: r(:), c(:)
        real(dp), allocatable :: v(:)
        type(fortsparse_status_t) :: status
        integer :: i, j, k, id

        allocate (r(5*m*m), c(5*m*m), v(5*m*m))
        k = 0
        do j = 1, m
            do i = 1, m
                id = (j - 1)*m + i
                call add(id, id, 4.0_dp)
                if (i > 1) call add(id, id - 1, -1.0_dp)
                if (i < m) call add(id, id + 1, -1.0_dp)
                if (j > 1) call add(id, id - m, -1.0_dp)
                if (j < m) call add(id, id + m, -1.0_dp)
            end do
        end do
        call csc_from_triplet(m*m, m*m, r(1:k), c(1:k), v(1:k), A, status)
        if (.not. status_ok(status)) error stop "laplace2d: assembly failed"
    contains
        subroutine add(ii, jj, x)
            integer,  intent(in) :: ii, jj
            real(dp), intent(in) :: x
            k = k + 1
            r(k) = ii
            c(k) = jj
            v(k) = x
        end subroutine add
    end subroutine laplace2d

    ! P3 Poisson stiffness matrix on an nc-by-nc grid of squares, each split
    ! along its diagonal into two triangles. Returns the SPD matrix A with
    ! n = (nc+1)^2 + 2*(3 nc^2 + 2 nc) + 2 nc^2 degrees of freedom.
    subroutine p3_poisson(nc, A, xyd)
        integer,     intent(in)  :: nc
        type(csc_t), intent(out) :: A
        real(dp), allocatable, intent(out), optional :: xyd(:, :) ! dof coordinates

        ! Dunavant degree-4 rule, six points, weights summing to one.
        real(dp), parameter :: pa = 0.445948490915965_dp
        real(dp), parameter :: pb = 0.091576213509771_dp
        real(dp), parameter :: wa = 0.223381589678011_dp
        real(dp), parameter :: wb = 0.109951743655322_dp
        real(dp) :: qp(3, 6), qw(6)
        real(dp) :: xy(2, 3), g(2, 3), dphi(3, 10, 6), ke(10, 10)
        integer  :: nvert, nh, nv, nd, nedge, ndof
        integer  :: gdof(10), tri(3)
        integer, allocatable :: r(:), c(:)
        real(dp), allocatable :: v(:)
        logical, allocatable :: bnd(:)
        type(fortsparse_status_t) :: status
        integer :: i, j, k, q, a1, b1, t, cnt
        real(dp) :: area, h

        h = 1.0_dp/real(nc, dp)
        nvert = (nc + 1)**2
        nh = nc*(nc + 1)   ! horizontal edges: nc per row, nc+1 rows
        nv = nc*(nc + 1)
        nd = nc*nc
        nedge = nh + nv + nd
        ndof = nvert + 2*nedge + 2*nd
        qp = 0.0_dp
        qw = 0.0_dp
        qp(:, 1) = [1.0_dp - 2.0_dp*pa, pa, pa]
        qp(:, 2) = [pa, 1.0_dp - 2.0_dp*pa, pa]
        qp(:, 3) = [pa, pa, 1.0_dp - 2.0_dp*pa]
        qp(:, 4) = [1.0_dp - 2.0_dp*pb, pb, pb]
        qp(:, 5) = [pb, 1.0_dp - 2.0_dp*pb, pb]
        qp(:, 6) = [pb, pb, 1.0_dp - 2.0_dp*pb]
        qw = [wa, wa, wa, wb, wb, wb]

        ! Reference-triangle local basis gradients with respect to barycentric
        ! coordinates; the element stiffness uses physical gradients of the
        ! barycentric coordinates, identical for every triangle of one kind.
        do q = 1, 6
            call p3_dphi(qp(:, q), dphi(:, :, q))
        end do

        allocate (bnd(ndof), source=.false.)
        allocate (r(100*2*nd + ndof), c(100*2*nd + ndof), v(100*2*nd + ndof))
        call mark_boundary()
        if (present(xyd)) call dof_coordinates()
        cnt = 0
        do j = 0, nc - 1
            do i = 0, nc - 1
                do t = 1, 2
                    if (t == 1) then
                        tri = [vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)]
                        xy(:, 1) = [i*h, j*h]
                        xy(:, 2) = [(i + 1)*h, j*h]
                        xy(:, 3) = [(i + 1)*h, (j + 1)*h]
                    else
                        tri = [vid(i, j), vid(i + 1, j + 1), vid(i, j + 1)]
                        xy(:, 1) = [i*h, j*h]
                        xy(:, 2) = [(i + 1)*h, (j + 1)*h]
                        xy(:, 3) = [i*h, (j + 1)*h]
                    end if
                    call element_dofs(i, j, t, tri, gdof)
                    call element_stiffness(xy, dphi, qw, ke, area, g)
                    do a1 = 1, 10
                        if (bnd(gdof(a1))) cycle
                        do b1 = 1, 10
                            if (bnd(gdof(b1))) cycle
                            cnt = cnt + 1
                            r(cnt) = gdof(a1)
                            c(cnt) = gdof(b1)
                            v(cnt) = ke(a1, b1)
                        end do
                    end do
                end do
            end do
        end do
        do k = 1, ndof
            if (bnd(k)) then
                cnt = cnt + 1
                r(cnt) = k
                c(cnt) = k
                v(cnt) = 1.0_dp
            end if
        end do
        call csc_from_triplet(ndof, ndof, r(1:cnt), c(1:cnt), v(1:cnt), A, status)
        if (.not. status_ok(status)) error stop "p3_poisson: assembly failed"

    contains

        integer function vid(ii, jj)
            integer, intent(in) :: ii, jj
            vid = jj*(nc + 1) + ii + 1
        end function vid

        ! edge ids (1-based within the edge numbering)
        integer function eid_h(ii, jj)
            integer, intent(in) :: ii, jj
            eid_h = jj*nc + ii + 1
        end function eid_h

        integer function eid_v(ii, jj)
            integer, intent(in) :: ii, jj
            eid_v = nh + ii*nc + jj + 1
        end function eid_v

        integer function eid_d(ii, jj)
            integer, intent(in) :: ii, jj
            eid_d = nh + nv + jj*nc + ii + 1
        end function eid_d

        ! Dof on edge e near its lower-id vertex (side=0) or higher (side=1).
        integer function edof(e, side)
            integer, intent(in) :: e, side
            edof = nvert + 2*(e - 1) + side + 1
        end function edof

        integer function idof(ii, jj, tt)
            integer, intent(in) :: ii, jj, tt
            idof = nvert + 2*nedge + 2*(jj*nc + ii) + tt
        end function idof

        ! Local ordering: 3 vertices; for each of the edges (1,2), (2,3),
        ! (3,1) the node near the first then the second vertex; interior.
        subroutine element_dofs(ii, jj, tt, tv, gd)
            integer, intent(in)  :: ii, jj, tt, tv(3)
            integer, intent(out) :: gd(10)
            integer :: e12, e23, e31, s

            if (tt == 1) then
                e12 = eid_h(ii, jj); e23 = eid_v(ii + 1, jj); e31 = eid_d(ii, jj)
            else
                e12 = eid_d(ii, jj); e23 = eid_h(ii, jj + 1); e31 = eid_v(ii, jj)
            end if
            gd(1:3) = tv
            ! Edge e12: tv(1) < tv(2) holds for both triangle kinds; likewise
            ! e23 with tv(2) vs tv(3) (kind 1: increasing; kind 2: tv(2) > tv(3)).
            gd(4) = edof(e12, 0); gd(5) = edof(e12, 1)
            if (tv(2) < tv(3)) then
                s = 0
            else
                s = 1
            end if
            gd(6) = edof(e23, s); gd(7) = edof(e23, 1 - s)
            ! Edge e31 runs 3 -> 1; lower vertex is tv(1) (a) in both kinds.
            gd(8) = edof(e31, 1); gd(9) = edof(e31, 0)
            gd(10) = idof(ii, jj, tt)
        end subroutine element_dofs

        subroutine dof_coordinates()
            integer :: ii, jj, s
            real(dp) :: f
            allocate (xyd(2, ndof))
            do jj = 0, nc
                do ii = 0, nc
                    xyd(:, vid(ii, jj)) = [ii*h, jj*h]
                end do
            end do
            do s = 0, 1
                f = (1.0_dp + s)/3.0_dp
                do jj = 0, nc
                    do ii = 0, nc - 1
                        xyd(:, edof(eid_h(ii, jj), s)) = [(ii + f)*h, jj*h]
                        xyd(:, edof(eid_v(jj, ii), s)) = [jj*h, (ii + f)*h]
                    end do
                end do
                do jj = 0, nc - 1
                    do ii = 0, nc - 1
                        xyd(:, edof(eid_d(ii, jj), s)) = [(ii + f)*h, (jj + f)*h]
                    end do
                end do
            end do
            do jj = 0, nc - 1
                do ii = 0, nc - 1
                    xyd(:, idof(ii, jj, 1)) = [(ii + 2.0_dp/3.0_dp)*h, &
                        (jj + 1.0_dp/3.0_dp)*h]
                    xyd(:, idof(ii, jj, 2)) = [(ii + 1.0_dp/3.0_dp)*h, &
                        (jj + 2.0_dp/3.0_dp)*h]
                end do
            end do
        end subroutine dof_coordinates

        subroutine mark_boundary()
            integer :: ii
            do ii = 0, nc
                bnd(vid(ii, 0)) = .true.; bnd(vid(ii, nc)) = .true.
                bnd(vid(0, ii)) = .true.; bnd(vid(nc, ii)) = .true.
            end do
            do ii = 0, nc - 1
                bnd(edof(eid_h(ii, 0), 0)) = .true.
                bnd(edof(eid_h(ii, 0), 1)) = .true.
                bnd(edof(eid_h(ii, nc), 0)) = .true.
                bnd(edof(eid_h(ii, nc), 1)) = .true.
                bnd(edof(eid_v(0, ii), 0)) = .true.
                bnd(edof(eid_v(0, ii), 1)) = .true.
                bnd(edof(eid_v(nc, ii), 0)) = .true.
                bnd(edof(eid_v(nc, ii), 1)) = .true.
            end do
        end subroutine mark_boundary

    end subroutine p3_poisson

    ! Partial derivatives of the ten P3 basis functions with respect to the
    ! three barycentric coordinates l at one point. Local order matches
    ! element_dofs: vertices, then per edge (1,2), (2,3), (3,1) the node near
    ! the first and second vertex, then the interior bubble.
    subroutine p3_dphi(l, d)
        real(dp), intent(in)  :: l(3)
        real(dp), intent(out) :: d(3, 10)

        integer :: i, j, e
        integer, parameter :: ei(3) = [1, 2, 3], ej(3) = [2, 3, 1]

        d = 0.0_dp
        do i = 1, 3
            d(i, i) = 0.5_dp*(27.0_dp*l(i)**2 - 18.0_dp*l(i) + 2.0_dp)
        end do
        do e = 1, 3
            i = ei(e)
            j = ej(e)
            ! node near i: 9/2 l_i l_j (3 l_i - 1)
            d(i, 3 + 2*e - 1) = 4.5_dp*l(j)*(6.0_dp*l(i) - 1.0_dp)
            d(j, 3 + 2*e - 1) = 4.5_dp*(3.0_dp*l(i)**2 - l(i))
            ! node near j: 9/2 l_i l_j (3 l_j - 1)
            d(j, 3 + 2*e) = 4.5_dp*l(i)*(6.0_dp*l(j) - 1.0_dp)
            d(i, 3 + 2*e) = 4.5_dp*(3.0_dp*l(j)**2 - l(j))
        end do
        d(1, 10) = 27.0_dp*l(2)*l(3)
        d(2, 10) = 27.0_dp*l(1)*l(3)
        d(3, 10) = 27.0_dp*l(1)*l(2)
    end subroutine p3_dphi

    subroutine element_stiffness(xy, dphi, qw, ke, area, g)
        real(dp), intent(in)  :: xy(2, 3), dphi(3, 10, 6), qw(6)
        real(dp), intent(out) :: ke(10, 10), area, g(2, 3)

        real(dp) :: det, gram(3, 3), t(10, 3)
        integer  :: q, k, l, a, b

        det = (xy(1, 2) - xy(1, 1))*(xy(2, 3) - xy(2, 1)) &
            - (xy(1, 3) - xy(1, 1))*(xy(2, 2) - xy(2, 1))
        area = 0.5_dp*abs(det)
        ! gradients of the barycentric coordinates
        g(:, 1) = [xy(2, 2) - xy(2, 3), xy(1, 3) - xy(1, 2)]/det
        g(:, 2) = [xy(2, 3) - xy(2, 1), xy(1, 1) - xy(1, 3)]/det
        g(:, 3) = [xy(2, 1) - xy(2, 2), xy(1, 2) - xy(1, 1)]/det
        do k = 1, 3
            do l = 1, 3
                gram(k, l) = dot_product(g(:, k), g(:, l))
            end do
        end do
        ke = 0.0_dp
        do q = 1, 6
            do a = 1, 10
                do l = 1, 3
                    t(a, l) = 0.0_dp
                    do k = 1, 3
                        t(a, l) = t(a, l) + dphi(k, a, q)*gram(k, l)
                    end do
                end do
            end do
            do b = 1, 10
                do a = 1, 10
                    ke(a, b) = ke(a, b) + area*qw(q)*dot_product(t(a, :), dphi(:, b, q))
                end do
            end do
        end do
    end subroutine element_stiffness

    ! Scaled backward error ||b - A x||_inf / (||A||_inf ||x||_inf + ||b||_inf).
    function inf_residual(A, x, b) result(res)
        use fortsparse, only: csc_matvec
        type(csc_t), intent(in) :: A
        real(dp),    intent(in) :: x(:), b(:)
        real(dp) :: res
        real(dp), allocatable :: r(:), rs(:)
        integer :: e

        r = b - csc_matvec(A, x)
        allocate (rs(A%nrow), source=0.0_dp)
        do e = 1, A%nnz
            rs(A%row_idx(e)) = rs(A%row_idx(e)) + abs(A%val(e))
        end do
        res = maxval(abs(r))/(maxval(rs)*maxval(abs(x)) + maxval(abs(b)))
    end function inf_residual

end module fs_test_matrices
