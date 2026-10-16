program test_fortsparse_native_supernodes
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csc, only: csc_t, csc_from_triplet, csc_matvec
    use fortsparse_status, only: fortsparse_status_t, FORTSPARSE_OK, &
        FORTSPARSE_NOT_FACTORED, FORTSPARSE_NOT_POSITIVE_DEFINITE
    use fortsparse_native_supernodes, only: supernodes_t
    use fortsparse_native_cuda, only: native_cuda_available
    implicit none
    integer, parameter :: sizes(*) = [1, 15, 16, 17, 31, 32, 33, 65, 193, 257, 513]
    integer :: k, device
    character(16) :: argument
    call get_command_argument(1, argument)
    device = -1
    if (trim(argument) == '--cuda-levels') then
        if (.not. native_cuda_available) error stop 'CUDA fixture needs CUDA build'
        device = 0
        ! Held-out rectangular graph and disconnected chains exercise the
        ! dependency schedule through the actual sparse analysis/solve path.
        call run_grid(19, 11)
        call run_case(71, .false.)
    else if (trim(argument) == '--cuda') then
        if (.not. native_cuda_available) error stop 'CUDA fixture needs CUDA build'
        device = 0
        call run_case(17, .false.)
        call run_case(65, .false.)
        call run_case(193, .true.)
    else
        do k = 1, size(sizes)
            if (sizes(k) == 513) then
                call run_interleaved()
            else
                call run_case(sizes(k), sizes(k) >= 193)
            end if
        end do
    end if
    if (trim(argument) /= '--cuda-levels') then
        call run_grid(17, 17)
        call run_scattered()
    end if
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

    subroutine run_scattered()
        integer, parameter :: n = 97
        real(dp) :: matrix(n,n)
        integer :: i, j
        matrix = 0.0_dp
        do j = 1, n
            if (mod(j,7) == 0) cycle
            do i = j+1, n
                if (mod(i,7) == 0) cycle
                if (i /= j+2 .and. i /= j+9) cycle
                matrix(i,j) = -0.125_dp
                matrix(j,i) = matrix(i,j)
            end do
        end do
        do j = 1, n
            matrix(j,j) = 1.0_dp+sum(abs(matrix(:,j)))
        end do
        call verify(matrix, .false.)
    end subroutine run_scattered

    subroutine run_interleaved()
        integer, parameter :: n = 513
        real(dp) :: matrix(n,n)
        integer :: i, j
        ! Two independent dense components interleaved in panel order. Small
        ! complete subtrees therefore have noncontiguous group membership.
        matrix = 0.0_dp
        do j = 1, n
            do i = j+2, n, 2
                matrix(i,j) = -0.01_dp/(1+i-j)
                matrix(j,i) = matrix(i,j)
            end do
        end do
        ! Set diagonal dominance after both triangles are complete.
        do j = 1, n
            matrix(j,j) = 1.0_dp+sum(abs(matrix(:,j)))
        end do
        call verify(matrix, .true.)
    end subroutine run_interleaved

    subroutine verify(matrix, parallel_expected)
        real(dp), intent(in) :: matrix(:,:)
        logical, intent(in) :: parallel_expected
        type(csc_t) :: a
        type(supernodes_t) :: solver
        type(fortsparse_status_t) :: status
        logical, allocatable :: graph(:,:)
        integer, allocatable :: rows(:), columns(:), lp(:), li(:)
        real(dp), allocatable :: vals(:), x(:), b(:), expected(:), block(:,:), empty(:,:)
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
        call solver%build(n, lp, li, 4, status, device)
        call checked(status)
        if (solver%memory_bytes() < 8_i8*solver%factor_nnz) error stop 'factor memory accounting'
        if (solver%symbolic_nnz /= size(li,kind=i8)) error stop 'symbolic factor diagnostics'
        if (4_i8*solver%factor_nnz > 5_i8*solver%symbolic_nnz) &
            error stop 'amalgamation exceeds the numeric factor padding budget'
        if (parallel_expected) then
            if (solver%factor_nnz /= solver%symbolic_nnz) &
                error stop 'complete dense pattern gained structural zeros'
        end if
        if (solver%factor_storage_bytes+solver%workspace_bytes /= solver%memory_bytes()) &
            error stop 'factor/workspace attribution does not cover owned payload'
        if (device < 0) then
            if (solver%active_threads < 1 .or. solver%active_threads > 4) &
                error stop 'invalid useful team limit'
            if (n == 513) then
                !$ if (solver%active_threads < 2) error stop 'missing grouped parallel path'
            end if
        else
            if (solver%device_bytes < 8_i8*solver%factor_nnz) error stop 'missing device accounting'
        end if
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
        if (n == 33 .or. n == 193) then
            allocate(block(n,9), empty(n,0))
            do j = 1, 9
                block(:,j) = real(j-5,dp)*b
            end do
            call solver%solve_block(block, status)
            call checked(status)
            do j = 1, 9
                if (maxval(abs(block(:,j)-real(j-5,dp)*expected)) > 1.2e-11_dp) &
                    error stop 'supernode block solve known solution'
            end do
            call solver%solve_block(empty, status)
            call checked(status)
        end if
        a%val = 0.75_dp*a%val
        call solver%factor(a, status)
        call checked(status)
        x = 0.75_dp*b
        call solver%solve(x, status)
        call checked(status)
        if (maxval(abs(x-expected)) > 3.0e-12_dp) then
            print '(a,2i8,es16.7)', 'refactor device, n, error: ', device, n, &
                maxval(abs(x-expected))
            error stop 'supernode refactor error'
        end if
        ! A nonuniform positive diagonal update keeps the graph and SPD property.
        ! Manufacture the new RHS from the original CSC, independent of panels.
        do j = 1, n
            do i = a%col_ptr(j), a%col_ptr(j+1)-1
                if (a%row_idx(i) == j) a%val(i) = a%val(i)+0.01_dp*real(1+mod(j,7),dp)
            end do
        end do
        call solver%factor(a, status)
        call checked(status)
        x = csc_matvec(a, expected)
        call solver%solve(x, status)
        call checked(status)
        if (maxval(abs(x-expected)) > 3.0e-12_dp) error stop 'nonuniform refactor solve'
        ! Restore the original scaled values before the failed-pivot/recovery test.
        do j = 1, n
            do i = a%col_ptr(j), a%col_ptr(j+1)-1
                a%val(i) = 0.75_dp*matrix(a%row_idx(i),j)
            end do
        end do
        a%val(1) = -1.0_dp
        call solver%factor(a, status)
        if (status%code /= FORTSPARSE_NOT_POSITIVE_DEFINITE) error stop 'missing supernode failed pivot'
        call solver%solve(x, status)
        if (status%code /= FORTSPARSE_NOT_FACTORED) error stop 'stale supernode factor reused'
        a%val(1) = 0.75_dp*matrix(1,1)
        call solver%factor(a, status)
        call checked(status)
        x = 0.75_dp*b
        call solver%solve(x, status)
        call checked(status)
        if (maxval(abs(x-expected)) > 3.0e-12_dp) error stop 'recovery after failed factor'
        call solver%free()
        if (solver%memory_bytes() /= 0_i8) error stop 'supernode workspace retained after free'
        if (solver%factor_storage_bytes /= 0_i8 .or. solver%workspace_bytes /= 0_i8) &
            error stop 'stale supernode memory diagnostics after free'
    end subroutine verify

    subroutine checked(status)
        type(fortsparse_status_t), intent(in) :: status
        if (status%code /= FORTSPARSE_OK) then
            print '(i0,2a)', status%code, ': ', trim(status%msg)
            error stop 'supernode operation failed'
        end if
    end subroutine checked
end program test_fortsparse_native_supernodes
