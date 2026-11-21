program test_fortsparse_nd
    ! The native nested dissection (fortsparse_nd) returns the orderings of
    ! METIS 5.2.1 METIS_NodeND (default options, seed 0) on synthetic graphs
    ! that reach each of its paths: one-level and five-run bisection (5000
    ! vertices and more), compression with two separators, disconnected
    ! graphs, isolated vertices, a star and a path. Oracle: the permutations
    ! METIS produced, recorded in test/data/nd_metis_5.2.1.txt, and with
    ! FORTSPARSE_TEST_METIS the linked METIS itself. With that macro and the
    ! environment variable FORTSPARSE_ND_WRITE set, the program rewrites the
    ! recording from METIS.
    use, intrinsic :: iso_c_binding, only: c_int
    use fortsparse_nd, only: nd_order
    implicit none
#ifdef FORTSPARSE_TEST_METIS
    interface
        integer(c_int) function metis_order(n, nnz, ptr, row, perm, inverse) &
                bind(C, name='fsparse_metis_order')
            import c_int
            integer(c_int), value :: n, nnz
            integer(c_int), intent(in) :: ptr(*), row(*)
            integer(c_int), intent(out) :: perm(*), inverse(*)
        end function metis_order
    end interface
#endif
    character(len=:), allocatable :: path
    character(len=512) :: arg
    integer :: unit, failures, cases
    logical :: write_mode

    path = 'test/data/nd_metis_5.2.1.txt'
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg)
        path = trim(arg)
    end if
    write_mode = .false.
#ifdef FORTSPARSE_TEST_METIS
    call get_environment_variable('FORTSPARSE_ND_WRITE', arg)
    write_mode = len_trim(arg) > 0
#endif
    if (write_mode) then
        open (newunit=unit, file=path, status='replace', action='write')
    else
        open (newunit=unit, file=path, status='old', action='read')
    end if
    failures = 0
    cases = 0
    call grid(40, 40, 1, 1, 'grid 40x40')
    call grid(12, 12, 12, 1, 'grid 12x12x12')
    call grid(40, 40, 1, 2, 'grid 40x40, 2 coupled unknowns per node')
    call grid(80, 70, 1, 1, 'grid 80x70')
    call random_graph(1500, 'random 1500')
    call two_grids(30, 20, 'two disjoint grids')
    call islands(1500, 'sparse with isolated vertices')
    call star(600, 'star 600')
    call path_graph(500, 'path 500')
    close (unit)
    if (failures > 0) error stop 1
    print '(a,i0,a)', 'PASS ', cases, ' graphs ordered as METIS 5.2.1'

contains

    !> Orders (xadj, adj) and compares with the recording (and METIS).
    subroutine compare(n, xadj, adj, name)
        integer, intent(in) :: n, xadj(0:), adj(0:)
        character(*), intent(in) :: name
        integer, allocatable :: perm(:), iperm(:), want(:)
        integer :: i, m
        character(len=128) :: label
        allocate (perm(0:n - 1), iperm(0:n - 1), want(0:n - 1))
        call nd_order(n, xadj, adj, perm, iperm)
        cases = cases + 1
        do i = 0, n - 1
            if (iperm(perm(i)) /= i) call fail('not a permutation: '//name)
        end do
#ifdef FORTSPARSE_TEST_METIS
        call metis_reference(n, xadj, adj, want)
        if (any(perm /= want)) call fail('differs from linked METIS: '//name)
        if (write_mode) then
            write (unit, '(a)') name
            write (unit, '(i0)') n
            write (unit, '(20(i0,:,1x))') want
            return
        end if
#endif
        read (unit, '(a)') label
        read (unit, *) m
        if (trim(label) /= name .or. m /= n) then
            call fail('recording out of step at '//name)
            return
        end if
        read (unit, *) want
        if (any(perm /= want)) call fail('differs from recorded METIS: '//name)
    end subroutine compare

#ifdef FORTSPARSE_TEST_METIS
    !> METIS through the reference adapter: the graph as one-based CSC
    !> (each list is ascending).
    subroutine metis_reference(n, xadj, adj, perm)
        integer, intent(in) :: n, xadj(0:), adj(0:)
        integer, intent(out) :: perm(0:)
        integer(c_int), allocatable :: p(:), ip(:)
        integer :: status
        allocate (p(n), ip(n))
        status = metis_order(n, xadj(n), xadj(0:n) + 1, adj(0:max(0, xadj(n) - 1)) + 1, p, ip)
        if (status /= 0) call fail('METIS reference failed')
        perm = p - 1
    end subroutine metis_reference
#endif

    subroutine fail(what)
        character(*), intent(in) :: what
        failures = failures + 1
        print '(a)', 'FAIL: '//what
    end subroutine fail

    !> nx x ny x nz grid, 6-neighbour coupling, dof unknowns per node all
    !> coupled to those of the node and of its neighbours (sorted lists).
    subroutine grid(nx, ny, nz, dof, name)
        integer, intent(in) :: nx, ny, nz, dof
        character(*), intent(in) :: name
        integer, allocatable :: xadj(:), adj(:)
        integer :: i, j, k, a, b, v, w, m, n, d
        integer, parameter :: ox(7) = [0, 0, -1, 0, 1, 0, 0], oy(7) = [0, -1, 0, 0, 0, 1, 0]
        integer, parameter :: oz(7) = [-1, 0, 0, 0, 0, 0, 1]
        n = nx*ny*nz*dof
        allocate (xadj(0:n), adj(0:7*dof*n))
        m = 0
        do k = 0, nz - 1
            do j = 0, ny - 1
                do i = 0, nx - 1
                    do a = 0, dof - 1
                        v = (i + nx*(j + ny*k))*dof + a
                        xadj(v) = m
                        do d = 1, 7
                            if (i + ox(d) < 0 .or. i + ox(d) >= nx .or. j + oy(d) < 0 .or. &
                                j + oy(d) >= ny .or. k + oz(d) < 0 .or. k + oz(d) >= nz) cycle
                            do b = 0, dof - 1
                                w = (i + ox(d) + nx*(j + oy(d) + ny*(k + oz(d))))*dof + b
                                if (w == v) cycle
                                adj(m) = w
                                m = m + 1
                            end do
                        end do
                    end do
                end do
            end do
        end do
        xadj(n) = m
        call compare(n, xadj, adj, name)
    end subroutine grid

    subroutine random_graph(n, name)
        integer, intent(in) :: n
        character(*), intent(in) :: name
        integer, allocatable :: eu(:), ev(:)
        integer :: i, k, e
        integer(8) :: state
        allocate (eu(5*n), ev(5*n))
        state = 12345
        e = 0
        do i = 0, n - 1
            do k = 1, 5
                ! Park-Miller: no integer overflow, the same graph at any optimization.
                state = mod(48271_8*state, 2147483647_8)
                e = e + 1
                eu(e) = i
                ev(e) = mod(i + 1 + int(mod(state, 50_8)), n)
                if (ev(e) == i) e = e - 1
            end do
        end do
        call from_edges(n, e, eu, ev, name)
    end subroutine random_graph

    subroutine two_grids(nx, ny, name)
        integer, intent(in) :: nx, ny
        character(*), intent(in) :: name
        integer, allocatable :: eu(:), ev(:)
        integer :: i, j, g, e, v
        allocate (eu(4*nx*ny), ev(4*nx*ny))
        e = 0
        do g = 0, 1
            do j = 0, ny - 1
                do i = 0, nx - 1
                    v = 2*(i + nx*j) + g
                    if (i < nx - 1) then
                        e = e + 1
                        eu(e) = v
                        ev(e) = v + 2
                    end if
                    if (j < ny - 1) then
                        e = e + 1
                        eu(e) = v
                        ev(e) = v + 2*nx
                    end if
                end do
            end do
        end do
        call from_edges(2*nx*ny, e, eu, ev, name)
    end subroutine two_grids

    subroutine islands(n, name)
        integer, intent(in) :: n
        character(*), intent(in) :: name
        integer, allocatable :: eu(:), ev(:)
        integer :: i, e
        allocate (eu(n), ev(n))
        e = 0
        do i = 0, n - 2
            if (mod(7*i, 10) < 6) then
                e = e + 1
                eu(e) = i
                ev(e) = mod(i + 1 + mod(13*i, 29), n)
                if (ev(e) == i) e = e - 1
            end if
        end do
        call from_edges(n, e, eu, ev, name)
    end subroutine islands

    subroutine star(n, name)
        integer, intent(in) :: n
        character(*), intent(in) :: name
        integer :: eu(n), ev(n), i
        do i = 1, n - 1
            eu(i) = 0
            ev(i) = i
        end do
        call from_edges(n, n - 1, eu, ev, name)
    end subroutine star

    subroutine path_graph(n, name)
        integer, intent(in) :: n
        character(*), intent(in) :: name
        integer :: eu(n), ev(n), i
        do i = 1, n - 1
            eu(i) = i - 1
            ev(i) = i
        end do
        call from_edges(n, n - 1, eu, ev, name)
    end subroutine path_graph

    !> The undirected graph of the edges, each list ascending, no repeats.
    subroutine from_edges(n, e, eu, ev, name)
        integer, intent(in) :: n, e, eu(:), ev(:)
        character(*), intent(in) :: name
        integer, allocatable :: xadj(:), adj(:), mark(:)
        logical, allocatable :: link(:, :)
        integer :: i, j, k, m
        ! Small enough for a dense adjacency, which sorts the lists for free.
        allocate (link(0:n - 1, 0:n - 1), xadj(0:n), adj(0:2*e), mark(0:n - 1))
        link = .false.
        do k = 1, e
            link(eu(k), ev(k)) = .true.
            link(ev(k), eu(k)) = .true.
        end do
        m = 0
        do i = 0, n - 1
            xadj(i) = m
            do j = 0, n - 1
                if (link(j, i) .and. j /= i) then
                    adj(m) = j
                    m = m + 1
                end if
            end do
        end do
        xadj(n) = m
        call compare(n, xadj, adj, name)
    end subroutine from_edges
end program test_fortsparse_nd
