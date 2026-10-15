module fs_bench_matrices
    ! Finite-element-shaped test matrices for the backend benchmark. The
    ! sparsity is that of a conforming assembly (all local degrees of freedom
    ! of an element couple); the values are an element-wise weighted graph
    ! Laplacian, optionally Kronecker-coupled across nf fields, so the
    ! structure, fill and pivoting behaviour follow the named discretisation:
    !   p1, p2, p3  2D Lagrange P_k on an nc x nc square split into triangles,
    !               real SPD (Grad-Shafranov / Laplace stiffness);
    !   hcurl       2D lowest-order edge + nodal unknowns per triangle (the P2
    !               graph), complex symmetric indefinite (per-n curl systems);
    !   blk3        2D P2 with three coupled fields, real nonsymmetric
    !               (coupled Newton blocks);
    !   ideal       2D P2 with three coupled fields, real symmetric indefinite
    !               (damped ideal-map LDL^T systems);
    !   tet1        3D P1 on an nc^3 cube split into six tetrahedra, real SPD;
    !   tet3        3D P1 with three coupled fields, real nonsymmetric.
    use fortsparse, only: dp, csc_t, csc_z_t, csc_from_triplet, &
        fortsparse_status_t, status_ok
    implicit none
    private

    public :: bench_matrix_real, bench_matrix_complex, bench_kind_is_complex

contains

    logical function bench_kind_is_complex(kind)
        character(*), intent(in) :: kind

        bench_kind_is_complex = (trim(kind) == "hcurl")
    end function bench_kind_is_complex

    ! Real matrix of the named kind; size is cells per side. shift scales the
    ! values (refactorization benchmarks pass a different shift, same pattern).
    subroutine bench_matrix_real(kind, size_, A, shift)
        character(*), intent(in) :: kind
        integer, intent(in) :: size_
        type(csc_t), intent(out) :: A
        real(dp), intent(in), optional :: shift

        integer, allocatable :: elems(:, :)
        integer :: ndof, nf
        real(dp) :: sigma, skew, s0

        s0 = 0.0_dp
        if (present(shift)) s0 = shift
        call kind_layout(kind, size_, elems, ndof, nf, sigma, skew)
        call assemble_real(elems, ndof, nf, sigma + s0, skew, A)
    end subroutine bench_matrix_real

    subroutine bench_matrix_complex(kind, size_, A, shift)
        character(*), intent(in) :: kind
        integer, intent(in) :: size_
        type(csc_z_t), intent(out) :: A
        real(dp), intent(in), optional :: shift

        type(csc_t) :: R, S
        integer, allocatable :: elems(:, :)
        integer :: ndof, nf
        real(dp) :: sigma, skew, s0

        s0 = 0.0_dp
        if (present(shift)) s0 = shift
        call kind_layout(kind, size_, elems, ndof, nf, sigma, skew)
        ! A = (L - sigma I) + i (0.2 L + 0.05 I): complex symmetric, indefinite
        ! real part, positive definite imaginary part, hence nonsingular.
        call assemble_real(elems, ndof, nf, sigma + s0, 0.0_dp, R)
        call assemble_real(elems, ndof, nf, 0.0_dp, 0.0_dp, S)
        A%nrow = R%nrow
        A%ncol = R%ncol
        A%nnz = R%nnz
        A%col_ptr = R%col_ptr
        A%row_idx = R%row_idx
        A%val = cmplx(R%val, 0.2_dp*S%val, dp)
        call add_diag_imag(A, 0.05_dp)
    end subroutine bench_matrix_complex

    subroutine kind_layout(kind, size_, elems, ndof, nf, sigma, skew)
        character(*), intent(in) :: kind
        integer, intent(in) :: size_
        integer, allocatable, intent(out) :: elems(:, :)
        integer, intent(out) :: ndof, nf
        real(dp), intent(out) :: sigma, skew

        nf = 1
        sigma = -1.0e-3_dp ! negative shift: SPD
        skew = 0.0_dp
        select case (trim(kind))
        case ("p1")
            call tri_mesh(size_, 1, elems, ndof)
        case ("p2")
            call tri_mesh(size_, 2, elems, ndof)
        case ("p3")
            call tri_mesh(size_, 3, elems, ndof)
        case ("hcurl")
            call tri_mesh(size_, 2, elems, ndof)
            sigma = 0.7_dp
        case ("blk3")
            call tri_mesh(size_, 2, elems, ndof)
            nf = 3
            skew = 0.3_dp
        case ("ideal")
            call tri_mesh(size_, 2, elems, ndof)
            nf = 3
            sigma = 0.7_dp
        case ("tet1")
            call tet_mesh(size_, elems, ndof)
        case ("tet3")
            call tet_mesh(size_, elems, ndof)
            nf = 3
            skew = 0.3_dp
        case default
            error stop "fs_bench_matrices: unknown kind"
        end select
    end subroutine kind_layout

    ! P_k Lagrange degrees of freedom on an nc x nc grid of squares, each split
    ! along its (i,j)-(i+1,j+1) diagonal: (k-1) per edge, (k-1)(k-2)/2 per
    ! triangle. elems(:, t) lists the dofs of triangle t.
    subroutine tri_mesh(nc, k, elems, ndof)
        integer, intent(in) :: nc, k
        integer, allocatable, intent(out) :: elems(:, :)
        integer, intent(out) :: ndof

        integer :: nv, nh, nvert, ndiag, ne, nloc, nint, i, j, t, base_e, base_i
        integer :: tv(3, 2), te(3, 2), tt, m, q

        nv = (nc + 1)**2
        nh = nc*(nc + 1)
        nvert = nc*(nc + 1)
        ndiag = nc*nc
        ne = nh + nvert + ndiag
        nint = (k - 1)*(k - 2)/2
        nloc = 3 + 3*(k - 1) + nint
        base_e = nv
        base_i = nv + (k - 1)*ne
        ndof = base_i + nint*2*nc*nc
        allocate (elems(nloc, 2*nc*nc))
        t = 0
        do j = 0, nc - 1
            do i = 0, nc - 1
                tv(:, 1) = [vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)]
                te(:, 1) = [eh(i, j), ev(i + 1, j), ed(i, j)]
                tv(:, 2) = [vid(i, j), vid(i + 1, j + 1), vid(i, j + 1)]
                te(:, 2) = [ed(i, j), eh(i, j + 1), ev(i, j)]
                do tt = 1, 2
                    t = t + 1
                    elems(1:3, t) = tv(:, tt)
                    m = 3
                    do q = 1, 3
                        call edge_dofs(te(q, tt), m, t)
                    end do
                    do q = 1, nint
                        elems(m + q, t) = base_i + (t - 1)*nint + q
                    end do
                end do
            end do
        end do

    contains

        subroutine edge_dofs(e, m, t)
            integer, intent(in) :: e, t
            integer, intent(inout) :: m
            integer :: r

            do r = 1, k - 1
                elems(m + r, t) = base_e + (e - 1)*(k - 1) + r
            end do
            m = m + k - 1
        end subroutine edge_dofs

        integer function vid(ii, jj)
            integer, intent(in) :: ii, jj
            vid = jj*(nc + 1) + ii + 1
        end function vid

        integer function eh(ii, jj) ! (ii,jj)-(ii+1,jj)
            integer, intent(in) :: ii, jj
            eh = jj*nc + ii + 1
        end function eh

        integer function ev(ii, jj) ! (ii,jj)-(ii,jj+1)
            integer, intent(in) :: ii, jj
            ev = nh + jj*(nc + 1) + ii + 1
        end function ev

        integer function ed(ii, jj) ! (ii,jj)-(ii+1,jj+1)
            integer, intent(in) :: ii, jj
            ed = nh + nvert + jj*nc + ii + 1
        end function ed
    end subroutine tri_mesh

    ! P1 on an nc^3 grid of cubes, each split into the six Kuhn tetrahedra
    ! around its main diagonal.
    subroutine tet_mesh(nc, elems, ndof)
        integer, intent(in) :: nc
        integer, allocatable, intent(out) :: elems(:, :)
        integer, intent(out) :: ndof

        integer, parameter :: perm(3, 6) = reshape([1, 2, 3, 1, 3, 2, 2, 1, 3, &
            2, 3, 1, 3, 1, 2, 3, 2, 1], [3, 6])
        integer :: i, j, l, p, t, c(3), d

        ndof = (nc + 1)**3
        allocate (elems(4, 6*nc**3))
        t = 0
        do l = 0, nc - 1
            do j = 0, nc - 1
                do i = 0, nc - 1
                    do p = 1, 6
                        t = t + 1
                        c = [i, j, l]
                        elems(1, t) = gid(c)
                        do d = 1, 3
                            c(perm(d, p)) = c(perm(d, p)) + 1
                            elems(d + 1, t) = gid(c)
                        end do
                    end do
                end do
            end do
        end do

    contains

        integer function gid(cc)
            integer, intent(in) :: cc(3)
            gid = (cc(3)*(nc + 1) + cc(2))*(nc + 1) + cc(1) + 1
        end function gid
    end subroutine tet_mesh

    ! Element-wise weighted graph Laplacian Kronecker field coupling F, minus
    ! sigma on the diagonal, plus an antisymmetric part skew * w_ab for a /= b.
    subroutine assemble_real(elems, ndof, nf, sigma, skew, A)
        integer, intent(in) :: elems(:, :), ndof, nf
        real(dp), intent(in) :: sigma, skew
        type(csc_t), intent(out) :: A

        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        type(fortsparse_status_t) :: st
        integer :: nloc, t, ia, ib, f, g, c, n
        real(dp) :: w, fld

        nloc = size(elems, 1)
        n = ndof*nf
        c = 2*size(elems, 2)*(nloc*nf)**2 + n
        allocate (rows(c), cols(c), vals(c))
        c = 0
        do t = 1, size(elems, 2)
            do ia = 1, nloc
                do ib = 1, nloc
                    if (ia == ib) cycle
                    w = weight(elems(ia, t), elems(ib, t))
                    do f = 1, nf
                        do g = 1, nf
                            fld = merge(2.0_dp, 0.5_dp, f == g)
                            ! -w off the diagonal, +w on both diagonals
                            call put(gix(elems(ia, t), f), gix(elems(ib, t), g), &
                                -w*fld + merge(skew*w, -skew*w, ia < ib))
                            call put(gix(elems(ia, t), f), gix(elems(ia, t), g), w*fld)
                        end do
                    end do
                end do
            end do
        end do
        do ia = 1, n
            call put(ia, ia, -sigma)
        end do
        call csc_from_triplet(n, n, rows(:c), cols(:c), vals(:c), A, st)
        if (.not. status_ok(st)) error stop "fs_bench_matrices: assembly failed"

    contains

        integer function gix(dof, field)
            integer, intent(in) :: dof, field
            gix = (dof - 1)*nf + field
        end function gix

        subroutine put(r, s, v)
            integer, intent(in) :: r, s
            real(dp), intent(in) :: v
            c = c + 1
            rows(c) = r
            cols(c) = s
            vals(c) = v
        end subroutine put
    end subroutine assemble_real

    ! Symmetric pseudo-random edge weight in [0.5, 1.5).
    pure real(dp) function weight(i, j)
        integer, intent(in) :: i, j
        real(dp) :: s

        s = sin(real(min(i, j), dp)*12.9898_dp + real(max(i, j), dp)*78.233_dp) &
            *43758.5453_dp
        weight = 0.5_dp + (s - floor(s))
    end function weight

    subroutine add_diag_imag(A, d)
        type(csc_z_t), intent(inout) :: A
        real(dp), intent(in) :: d
        integer :: j, k

        do j = 1, A%ncol
            do k = A%col_ptr(j), A%col_ptr(j + 1) - 1
                if (A%row_idx(k) == j) A%val(k) = A%val(k) + cmplx(0.0_dp, d, dp)
            end do
        end do
    end subroutine add_diag_imag

end module fs_bench_matrices
