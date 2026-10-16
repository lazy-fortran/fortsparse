program test_fortsparse_supernode_partition
    ! supernode_partition(parent, count) must reproduce the partition formed by
    ! scanning the explicit factor row structure (the pre-extraction rule).
    ! The structure is an independent dense graph elimination; parent/count
    ! come from native_symbolic_counts on the same permuted pattern.
    use, intrinsic :: iso_c_binding, only: c_int
    use fortsparse_kinds, only: i8
    use fortsparse_status, only: fortsparse_status_t, FORTSPARSE_OK
    use fortsparse_native_symbolic, only: native_symbolic_counts
    use fortsparse_native_supernode_schedule, only: supernode_partition
    implicit none
    interface
        integer(c_int) function amd_available() &
            bind(C, name='fsparse_native_order_amd_available')
            import :: c_int
        end function amd_available
        integer(c_int) function order_amd(n, nnz, ptr, row, perm, inverse) &
            bind(C, name='fsparse_native_order_amd')
            import :: c_int
            integer(c_int), value :: n, nnz
            integer(c_int), intent(in) :: ptr(*), row(*)
            integer(c_int), intent(out) :: perm(*), inverse(*)
        end function order_amd
    end interface
    integer, parameter :: sizes(*) = [1, 2, 9, 40, 97, 160, 260]
    real, parameter :: densities(*) = [0.0, 0.02, 0.08, 0.3, 1.0]
    integer :: k, d, cases, seed
    logical :: amd

    amd = amd_available() /= 0
    cases = 0
    seed = 12345
    do k = 1, size(sizes)
        do d = 1, size(densities)
            call check_pattern(sizes(k), densities(d), seed, cases)
        end do
    end do
    ! Chains with periodic gaps cross every width limit and padding bound.
    call check_banded(300, 1, cases)
    call check_banded(300, 5, cases)
    call check_banded(140, 139, cases)
    if (cases < 4*(size(sizes)*size(densities)+3)) &
        error stop 'partition oracle skipped cases'
    if (amd) then
        print '(a,i0,a)', 'Supernode partition matches the row-scan rule on ', &
            cases, ' natural/AMD strict/relaxed patterns'
    else
        print '(a,i0,a)', 'Supernode partition matches the row-scan rule on ', &
            cases, ' patterns (AMD unavailable: natural order repeated)'
    end if
contains
    subroutine check_pattern(n, density, seed, cases)
        integer, intent(in) :: n
        real, intent(in) :: density
        integer, intent(inout) :: seed, cases
        logical, allocatable :: graph(:,:)
        integer :: i, j
        allocate(graph(n,n))
        graph = .false.
        do j = 1, n
            graph(j,j) = .true.
            do i = j+1, n
                if (uniform(seed) < density) then
                    graph(i,j) = .true.
                    graph(j,i) = .true.
                end if
            end do
        end do
        call check_orders(graph, cases)
    end subroutine check_pattern

    subroutine check_banded(n, band, cases)
        integer, intent(in) :: n, band
        integer, intent(inout) :: cases
        logical, allocatable :: graph(:,:)
        integer :: i, j
        allocate(graph(n,n))
        graph = .false.
        do j = 1, n
            graph(j,j) = .true.
            if (mod(j,37) == 0) cycle
            do i = j+1, min(n,j+band)
                graph(i,j) = .true.
                graph(j,i) = .true.
            end do
        end do
        call check_orders(graph, cases)
    end subroutine check_banded

    subroutine check_orders(graph, cases)
        logical, intent(in) :: graph(:,:)
        integer, intent(inout) :: cases
        integer, allocatable :: ptr(:), row(:), perm(:), inverse(:)
        integer :: n, i, j, code
        n = size(graph,1)
        allocate(ptr(n+1), row(count(graph)), perm(n), inverse(n))
        call to_csc(graph, ptr, row)
        perm = [(i, i=1, n)]
        if (amd) then
            code = order_amd(n, size(row), ptr, row, perm, inverse)
            if (code /= 0) error stop 'AMD ordering failed'
            ! perm(new) = old must be a permutation.
            inverse = 0
            do j = 1, n
                if (perm(j) < 1 .or. perm(j) > n) error stop 'invalid AMD permutation'
                inverse(perm(j)) = j
            end do
            if (any(inverse == 0)) error stop 'invalid AMD permutation'
        end if
        call compare(graph, [(i, i=1, n)], cases)
        call compare(graph, perm, cases)
    end subroutine check_orders

    subroutine compare(graph, perm, cases)
        logical, intent(in) :: graph(:,:)
        integer, intent(in) :: perm(:)
        integer, intent(inout) :: cases
        logical, allocatable :: permuted(:,:), filled(:,:)
        integer, allocatable :: ptr(:), row(:), parent(:), count(:), lp(:), li(:)
        integer, allocatable :: starts(:), ends(:), m(:), ref_starts(:), ref_ends(:), ref_m(:)
        type(fortsparse_status_t) :: status
        integer :: n, i, j, q, nodes, ref_nodes, variant
        logical :: strict
        n = size(graph,1)
        allocate(permuted(n,n), parent(n), count(n), lp(n+1), starts(n), ends(n), m(n), &
            ref_starts(n), ref_ends(n), ref_m(n), ptr(n+1))
        do j = 1, n
            do i = 1, n
                permuted(i,j) = graph(perm(i), perm(j))
            end do
        end do
        allocate(row(count_true(permuted)))
        call to_csc(permuted, ptr, row)
        call native_symbolic_counts(n, ptr, row, parent, count, status)
        if (status%code /= FORTSPARSE_OK) error stop 'symbolic counts failed'
        ! Independent graph elimination for the lower factor structure.
        filled = permuted
        do j = 1, n
            do i = j+1, n
                if (.not. filled(i,j)) cycle
                do q = i+1, n
                    if (filled(q,j)) filled(q,i) = .true.
                end do
            end do
        end do
        lp(1) = 1
        do j = 1, n
            lp(j+1) = lp(j)+count_true(filled(j:n,j:j))
            if (lp(j+1)-lp(j)-1 /= count(j)) error stop 'counts disagree with elimination'
        end do
        allocate(li(lp(n+1)-1))
        q = 0
        do j = 1, n
            do i = j, n
                if (.not. filled(i,j)) cycle
                q = q+1
                li(q) = i
            end do
        end do
        do variant = 1, 2
            strict = variant == 2
            call reference_partition(n, lp, li, strict, ref_starts, ref_ends, ref_m, ref_nodes)
            call supernode_partition(n, parent, count, strict, starts, ends, m, nodes)
            if (nodes /= ref_nodes) error stop 'supernode count differs from row scan'
            if (any(starts(:nodes) /= ref_starts(:nodes))) error stop 'supernode starts differ'
            if (any(ends(:nodes) /= ref_ends(:nodes))) error stop 'supernode ends differ'
            if (any(m(:nodes) /= ref_m(:nodes))) error stop 'supernode row counts differ'
            cases = cases+1
        end do
    end subroutine compare

    subroutine reference_partition(n, lp, li, strict, starts, ends, counts, nodes)
        ! Verbatim rule of the former supernodes_t row scan (commit 600448b).
        integer, intent(in) :: n, lp(:), li(:)
        logical, intent(in) :: strict
        integer, intent(out) :: starts(:), ends(:), counts(:), nodes
        integer, allocatable :: mark(:)
        integer :: s, p, row, m, first, last, candidate_m, candidate_w, exact_limit
        integer(i8) :: original_count, candidate_count, padded
        exact_limit = 64
        if (.not. strict) then
            if (size(li,kind=i8) > 128_i8*n) exact_limit = 128
        end if
        allocate(mark(n))
        mark = 0
        nodes = 0
        first = 1
        do while (first <= n)
            nodes = nodes+1
            s = nodes
            last = first
            m = lp(first+1)-lp(first)
            original_count = m
            mark(li(lp(first):lp(first+1)-1)) = s
            do while (last < n)
                if (lp(last+1)-lp(last) < 2) exit
                if (li(lp(last)+1) /= last+1) exit
                candidate_w = last-first+2
                if (candidate_w > exact_limit) exit
                candidate_m = m
                do p = lp(last+1), lp(last+2)-1
                    row = li(p)
                    if (mark(row) == s) cycle
                    mark(row) = s
                    candidate_m = candidate_m+1
                end do
                candidate_count = original_count+lp(last+2)-lp(last+1)
                padded = int(candidate_m,i8)*candidate_w &
                    -int(candidate_w,i8)*(candidate_w-1)/2
                if (padded /= candidate_count) then
                    if (candidate_w > 32) exit
                    if (4_i8*padded > 5_i8*candidate_count) exit
                end if
                last = last+1
                m = candidate_m
                original_count = candidate_count
            end do
            starts(s) = first
            ends(s) = last
            counts(s) = m
            first = last+1
        end do
    end subroutine reference_partition

    subroutine to_csc(graph, ptr, row)
        logical, intent(in) :: graph(:,:)
        integer, intent(out) :: ptr(:), row(:)
        integer :: i, j, p
        p = 0
        ptr(1) = 1
        do j = 1, size(graph,2)
            do i = 1, size(graph,1)
                if (.not. graph(i,j)) cycle
                p = p+1
                row(p) = i
            end do
            ptr(j+1) = p+1
        end do
    end subroutine to_csc

    integer function count_true(graph)
        logical, intent(in) :: graph(:,:)
        count_true = count(graph)
    end function count_true

    real function uniform(seed)
        integer, intent(inout) :: seed
        ! Park-Miller minimal standard generator; reproducible across compilers.
        seed = int(mod(16807_i8*seed, 2147483647_i8))
        uniform = real(seed)/2147483647.0
    end function uniform
end program test_fortsparse_supernode_partition
