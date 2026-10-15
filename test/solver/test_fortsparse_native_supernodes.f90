program test_fortsparse_native_supernodes
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csc, only: csc_t, csc_from_triplet, csc_matvec
    use fortsparse_status, only: fortsparse_status_t, FORTSPARSE_OK, FORTSPARSE_NOT_FACTORED
    use fortsparse_native_supernodes, only: supernodes_t
    implicit none
    integer, parameter :: sizes(*) = [1, 15, 16, 17, 31, 32, 33, 65, 193]
    integer :: k
    do k = 1, size(sizes)
        call run_case(sizes(k), k == size(sizes))
    end do
    call run_grid(17, 17)
    print '(a)', 'Native supernode block boundaries, task execution, refactor and errors pass'
contains
    subroutine run_case(n, dense)
        integer, intent(in) :: n
        logical, intent(in) :: dense
        real(dp), allocatable :: matrix(:,:)
        integer :: i, j
        allocate(matrix(n,n))
        matrix = 0.0_dp
        do j = 1, n
            do i = j+1, n
                ! Sparse disconnected chains deliberately cross block boundaries.
                if (.not. dense) then
                    if (i /= j+1 .or. mod(j,11) == 0) cycle
                end if
                matrix(i,j) = -0.25_dp/(i-j)
                matrix(j,i) = matrix(i,j)
            end do
        end do
        do j = 1, n
            matrix(j,j) = 1.0_dp+sum(abs(matrix(:,j)))
        end do
        call verify(matrix, dense)
    end subroutine run_case

    subroutine run_grid(nx, ny)
        integer, intent(in) :: nx, ny
        real(dp), allocatable :: matrix(:,:)
        integer :: i, j, p, n
        n = nx*ny
        allocate(matrix(n,n))
        matrix = 0.0_dp
        do j = 1, ny
            do i = 1, nx
                p = i+(j-1)*nx
                matrix(p,p) = 6.0_dp
                if (i < nx) then
                    matrix(p,p+1) = -1.0_dp
                    matrix(p+1,p) = -1.0_dp
                end if
                if (j < ny) then
                    matrix(p,p+nx) = -1.5_dp
                    matrix(p+nx,p) = -1.5_dp
                end if
            end do
        end do
        call verify(matrix, .false.)
    end subroutine run_grid

    subroutine verify(matrix, parallel_expected)
        real(dp), intent(in) :: matrix(:,:)
        logical, intent(in) :: parallel_expected
        type(csc_t) :: a
        type(supernodes_t) :: solver
        type(fortsparse_status_t) :: status
        logical, allocatable :: graph(:,:)
        integer, allocatable :: rows(:), columns(:), lp(:), li(:)
        real(dp), allocatable :: vals(:), x(:), b(:), expected(:)
        integer :: n, i, j, p, q, nz

        n = size(matrix,1)
        allocate(graph(n,n), lp(n+1), x(n), b(n), expected(n))
        graph = matrix /= 0.0_dp
        ! Independent graph-elimination oracle; no native symbolic routines.
        do j = 1, n
            do i = j+1, n
                if (.not. graph(i,j)) cycle
                do q = i+1, n
                    if (graph(q,j)) graph(q,i) = .true.
                end do
            end do
        end do
        lp(1) = 1
        do j = 1, n
            lp(j+1) = lp(j)+count(graph(j:n,j))
        end do
        allocate(li(lp(n+1)-1))
        p = 1
        do j = 1, n
            do i = j, n
                if (.not. graph(i,j)) cycle
                li(p) = i
                p = p+1
            end do
        end do
        nz = count(matrix /= 0.0_dp)
        allocate(rows(nz), columns(nz), vals(nz))
        p = 0
        do j = 1, n
            do i = 1, n
                if (matrix(i,j) == 0.0_dp) cycle
                p = p+1
                rows(p) = i
                columns(p) = j
                vals(p) = matrix(i,j)
            end do
        end do
        call csc_from_triplet(n, n, rows, columns, vals, a, status)
        call checked(status)
        call solver%build(n, lp, li, 4, status)
        call checked(status)
        if (solver%memory_bytes() < 8_i8*solver%factor_nnz) error stop 'factor memory accounting'
        !$ if (parallel_expected .and. solver%active_threads /= 4) error stop 'missing parallel plan'
        do i = 1, n
            expected(i) = sin(real(i,dp)/7.0_dp)
        end do
        call solver%factor(a, status)
        call checked(status)
        b = csc_matvec(a, expected)
        x = b
        call solver%solve(x, status)
        call checked(status)
        if (maxval(abs(x-expected)) > 3.0e-12_dp) error stop 'supernode known-solution error'
        if (maxval(abs(matmul(matrix,x)-b)) > 3.0e-12_dp) error stop 'supernode residual'
        a%val = 0.75_dp*a%val
        call solver%factor(a, status)
        call checked(status)
        x = 0.75_dp*b
        call solver%solve(x, status)
        call checked(status)
        if (maxval(abs(x-expected)) > 3.0e-12_dp) error stop 'supernode refactor error'
        a%val(1) = -1.0_dp
        call solver%factor(a, status)
        if (status%code /= 8) error stop 'missing supernode failed pivot'
        call solver%solve(x, status)
        if (status%code /= FORTSPARSE_NOT_FACTORED) error stop 'stale supernode factor reused'
        call solver%free()
        if (solver%memory_bytes() /= 0_i8) error stop 'supernode workspace retained after free'
    end subroutine verify

    subroutine checked(status)
        type(fortsparse_status_t), intent(in) :: status
        if (status%code /= FORTSPARSE_OK) then
            print '(i0,2a)', status%code, ': ', trim(status%msg)
            error stop 'supernode operation failed'
        end if
    end subroutine checked
end program test_fortsparse_native_supernodes
