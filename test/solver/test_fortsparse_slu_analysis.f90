program test_fortsparse_slu_analysis
    ! Static analysis of the supernodal LU against independent oracles:
    ! dense boolean elimination of the symmetrized matched pattern, the
    ! elimination-tree postorder property, the matching/scaling contract,
    ! the exact layout plan, and refusal at every ledger phase.
    use fortsparse_kinds, only: dp, i8
    use fortsparse_memory, only: memory_budget_t, sparse_memory_plan_t, &
        MEMORY_PLAN_EXACT_LAYOUT
    use fortsparse_status, only: fortsparse_status_t, FORTSPARSE_OK, &
        FORTSPARSE_RESOURCE_LIMIT, FORTSPARSE_SINGULAR, FORTSPARSE_INVALID_MATRIX, &
        FORTSPARSE_BACKEND_UNAVAILABLE
    use fortsparse_native_amd, only: bounded_amd_available
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
    use fortsparse_native_slu, only: slu_factor_t, SLU_ORDER_AMD, SLU_ORDER_ND, &
        SLU_PHASE_COUNT
    implicit none
    integer :: seed, trial, cases
    integer, allocatable :: ptr(:), rows(:)
    real(dp), allocatable :: values(:)

    seed = 20261016
    cases = 0
    do trial = 1, 40
        call random_case(2+mod(trial*7, 39), trial, ptr, rows, values)
        call check_structure(ptr, rows, values, SLU_ORDER_AMD, .false.)
        if (mod(trial, 4) == 0) call check_structure(ptr, rows, values, SLU_ORDER_AMD, .true.)
        if (mod(trial, 5) == 0) call check_structure(ptr, rows, values, SLU_ORDER_ND, .false.)
        cases = cases+1
    end do
    call random_case(37, 3, ptr, rows, values)
    call check_injection(ptr, rows, values)
    call check_cap(ptr, rows, values)
    call check_refusals()
    print '(a,i0,a)', 'supernodal LU analysis: ', cases, ' random patterns passed'
contains
    real(dp) function uniform()
        seed = modulo(1103515245_i8*seed+12345_i8, 2147483648_i8)
        uniform = real(seed, dp)/2147483648.0_dp
    end function uniform

    subroutine random_case(n, trial, ptr, rows, values)
        ! A random permutation guarantees a perfect matching; about a third of
        ! the cases have a zero (absent) diagonal. Off-diagonal entries are
        ! mostly mirrored (near-symmetric pattern) with a few one-sided ones.
        integer, intent(in) :: n, trial
        integer, allocatable, intent(out) :: ptr(:), rows(:)
        real(dp), allocatable, intent(out) :: values(:)
        logical :: pattern(n,n)
        real(dp) :: dense(n,n)
        integer :: perm(n), i, j, k, swap
        pattern = .false.
        perm = [(i, i=1, n)]
        if (mod(trial, 3) == 0) then
            do i = n, 2, -1
                j = 1+int(uniform()*i)
                swap = perm(i)
                perm(i) = perm(j)
                perm(j) = swap
            end do
        end if
        do j = 1, n
            pattern(perm(j), j) = .true.
        end do
        do k = 1, 2*n
            i = 1+int(uniform()*n)
            j = 1+int(uniform()*n)
            pattern(i,j) = .true.
            if (uniform() < 0.85_dp) pattern(j,i) = .true.
        end do
        if (mod(trial, 3) == 0) then
            do j = 1, n
                if (perm(j) /= j) pattern(j,j) = .false.
            end do
        end if
        dense = 0.0_dp
        do j = 1, n
            do i = 1, n
                if (pattern(i,j)) dense(i,j) = (uniform()-0.5_dp)*10.0_dp**(4*uniform()-2)
            end do
            dense(perm(j),j) = 3.0_dp+uniform()
        end do
        call to_csc(dense, pattern, ptr, rows, values)
    end subroutine random_case

    subroutine to_csc(dense, pattern, ptr, rows, values)
        real(dp), intent(in) :: dense(:,:)
        logical, intent(in) :: pattern(:,:)
        integer, allocatable, intent(out) :: ptr(:), rows(:)
        real(dp), allocatable, intent(out) :: values(:)
        integer :: n, i, j, k
        n = size(dense,1)
        allocate(ptr(n+1), rows(count(pattern)), values(count(pattern)))
        k = 0
        ptr(1) = 1
        do j = 1, n
            do i = 1, n
                if (.not. pattern(i,j)) cycle
                k = k+1
                rows(k) = i
                values(k) = dense(i,j)
            end do
            ptr(j+1) = k+1
        end do
    end subroutine to_csc

    subroutine check_structure(ptr, rows, values, order, low)
        integer, intent(in) :: ptr(:), rows(:), order
        real(dp), intent(in) :: values(:)
        logical, intent(in) :: low
        type(slu_factor_t) :: lu
        type(fortsparse_status_t) :: status
        type(memory_budget_t) :: budget
        type(sparse_memory_plan_t) :: layout
        integer, allocatable :: row_old(:), col_old(:), row_new(:), col_new(:)
        integer, allocatable :: nrows(:), parent(:), size_tree(:)
        real(dp), allocatable :: rs(:), cs(:)
        logical, allocatable :: s(:,:)
        integer :: n, i, j, k, l, e, first, width, m, node, covered, a
        integer(i8) :: nnz_l
        real(dp) :: scaled

        n = size(ptr)-1
        call lu%analyze(n, ptr, rows, values, status, ordering=order, &
            low_memory=low, budget=budget)
        if (order == SLU_ORDER_ND .and. status%code == FORTSPARSE_BACKEND_UNAVAILABLE) return
        call check(status%code == FORTSPARSE_OK, 'analysis succeeds: '//trim(status%msg))
        call check(lu%plan_bytes == lu%memory_bytes(), 'exact plan equals retained layout')
        call check(budget%current == lu%memory_bytes(), 'ledger equals retained layout')
        layout = lu%plan()
        call check(layout%prediction_quality == MEMORY_PLAN_EXACT_LAYOUT, 'plan is exact layout')
        call check(layout%factor_host_bytes+layout%persistent_host_bytes &
            +layout%workspace_host_bytes == lu%memory_bytes(), 'plan classes add up')
        call check(budget%peak >= budget%current, 'ledger peak covers retained payload')
        call lu%ordering(row_old, col_old, rs, cs)
        allocate(row_new(n), col_new(n), s(n,n))
        call check(all(row_old >= 1) .and. all(row_old <= n), 'row labels in range')
        row_new = 0
        col_new = 0
        do k = 1, n
            row_new(row_old(k)) = k
            col_new(col_old(k)) = k
        end do
        call check(all(row_new > 0) .and. all(col_new > 0), 'labels are permutations')
        ! Scaled matched matrix: |M| <= 1 with a unit diagonal on matched entries.
        s = .false.
        do k = 1, n
            s(k,k) = .true.
        end do
        do j = 1, n
            l = col_new(j)
            do e = ptr(j), ptr(j+1)-1
                k = row_new(rows(e))
                s(k,l) = .true.
                s(l,k) = .true.
                scaled = rs(k)*values(e)*cs(l)
                call check(abs(scaled) <= 1.0_dp+1.0e-12_dp, 'scaled entries bounded by one')
                if (k == l) call check(abs(abs(scaled)-1.0_dp) <= 1.0e-12_dp, &
                    'matched diagonal scaled to one')
            end do
        end do
        do k = 1, n
            a = 0
            do e = ptr(col_old(k)), ptr(col_old(k)+1)-1
                if (rows(e) == row_old(k)) a = e
            end do
            call check(a > 0, 'matched entry is structurally present')
            call check(values(a) /= 0.0_dp, 'matched entry is numerically nonzero')
        end do
        ! Dense boolean elimination of S = M + M^T + I in the final order.
        do j = 1, n
            do i = j+1, n
                if (.not. s(i,j)) cycle
                do k = j+1, n
                    if (s(k,j)) then
                        s(k,i) = .true.
                        s(i,k) = .true.
                    end if
                end do
            end do
        end do
        nnz_l = 0
        allocate(parent(n), size_tree(n))
        do j = 1, n
            parent(j) = 0
            do i = n, j, -1
                if (.not. s(i,j)) cycle
                nnz_l = nnz_l+1
                if (i > j) parent(j) = i
            end do
        end do
        call check(lu%symbolic_nnz == 2_i8*nnz_l-n, 'symbolic nnz(L+U) equals dense elimination')
        ! Postorder: every subtree occupies the contiguous range ending at its root.
        size_tree = 1
        do j = 1, n
            if (parent(j) > 0) size_tree(parent(j)) = size_tree(parent(j))+size_tree(j)
        end do
        do j = 1, n
            do i = j-size_tree(j)+1, j
                call check(i >= 1, 'subtree range inside matrix')
                k = i
                do while (k /= j .and. k /= 0)
                    k = parent(k)
                end do
                call check(k == j, 'elimination tree is postordered')
            end do
        end do
        ! Supernode rows equal the column union cols(s) U struct(L(:,last)).
        covered = 0
        do node = 1, lu%node_count
            call lu%node_structure(node, first, width, m, nrows)
            call check(first == covered+1, 'supernodes are contiguous')
            covered = first+width-1
            do i = 1, width
                call check(nrows(i) == first+i-1, 'leading rows are the supernode columns')
            end do
            k = width
            do i = first+width, n
                if (.not. s(i,first+width-1)) cycle
                k = k+1
                call check(k <= m, 'row list holds the last column structure')
                if (k <= m) call check(nrows(k) == i, 'row list equals last column structure')
            end do
            call check(k == m, 'row list has no extra rows')
            do j = first, first+width-1
                do i = j, n
                    if (s(i,j)) call check(any(nrows == i), 'row list covers each column')
                end do
            end do
        end do
        call check(covered == n, 'supernodes cover all columns')
        call lu%free(budget)
        call check(budget%current == 0, 'free releases the ledger')
        call check(lu%memory_bytes() == 0, 'free releases all arrays')
    end subroutine check_structure

    subroutine check_injection(ptr, rows, values)
        integer, intent(in) :: ptr(:), rows(:)
        real(dp), intent(in) :: values(:)
        type(slu_factor_t) :: lu
        type(fortsparse_status_t) :: status
        type(memory_budget_t) :: budget
        integer :: phase, n
        n = size(ptr)-1
        do phase = 1, SLU_PHASE_COUNT
            budget = memory_budget_t()
            budget%current = 64
            budget%peak = 64
            call lu%analyze(n, ptr, rows, values, status, budget=budget, &
                inject_phase=phase)
            call check(status%code == FORTSPARSE_RESOURCE_LIMIT, 'injected phase refused')
            call check(lu%failed_phase == phase, 'refusal names the injected phase')
            call check(lu%failed_required_bytes > 0, 'required bytes reported')
            call check(lu%failed_available_bytes < lu%failed_required_bytes, &
                'available below required')
            call check(budget%current == 64, 'ledger restored after refusal')
            call check(lu%memory_bytes() == 0, 'no array retained after refusal')
            call check(.not. lu%is_analyzed(), 'refused analysis is not usable')
        end do
    end subroutine check_injection

    subroutine check_cap(ptr, rows, values)
        ! Peak of a strict run is an exact admission threshold.
        integer, intent(in) :: ptr(:), rows(:)
        real(dp), intent(in) :: values(:)
        type(slu_factor_t) :: lu
        type(fortsparse_status_t) :: status
        type(memory_budget_t) :: budget
        integer(i8) :: peak
        integer :: n
        if (bounded_amd_available() /= 1) return
        n = size(ptr)-1
        budget%limit = huge(1_i8)/4
        call lu%analyze(n, ptr, rows, values, status, budget=budget)
        call check(status%code == FORTSPARSE_OK, 'strict reference analysis')
        peak = budget%peak
        call lu%free(budget)
        budget = memory_budget_t()
        budget%limit = peak
        call lu%analyze(n, ptr, rows, values, status, budget=budget)
        call check(status%code == FORTSPARSE_OK, 'analysis fits its exact peak')
        call check(budget%peak == peak, 'capped peak reproduces the reference')
        call lu%free(budget)
        budget = memory_budget_t()
        budget%limit = peak-1
        call lu%analyze(n, ptr, rows, values, status, budget=budget)
        call check(status%code == FORTSPARSE_RESOURCE_LIMIT, 'one byte short is refused')
        call check(budget%current == 0, 'refused cap leaves no ledger residue')
        call check(lu%failed_required_bytes > lu%failed_available_bytes, &
            'cap refusal reports required/available')
        budget = memory_budget_t()
        budget%limit = peak
        call lu%analyze(n, ptr, rows, values, status, ordering=SLU_ORDER_ND, &
            budget=budget)
        call check(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, 'ND refused under a cap')
        call check(budget%current == 0, 'ND refusal leaves no ledger residue')
    end subroutine check_cap

    subroutine check_refusals()
        type(slu_factor_t) :: lu
        type(fortsparse_status_t) :: status
        integer :: ptr(4), rows(4)
        real(dp) :: values(4)
        ! Columns 1 and 2 both have entries only in row 1: structurally singular.
        ptr = [1, 2, 3, 5]
        rows = [1, 1, 2, 3]
        values = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
        call lu%analyze(3, ptr, rows, values, status)
        call check(status%code == FORTSPARSE_SINGULAR, 'structurally singular is SINGULAR')
        call check(lu%memory_bytes() == 0, 'singular analysis retains nothing')
        rows = [1, 2, 2, 3]
        values(2) = ieee_value(1.0_dp, ieee_positive_inf)
        call lu%analyze(3, ptr, rows, values, status)
        call check(status%code == FORTSPARSE_INVALID_MATRIX, 'non-finite input refused')
        values(2) = 1.0_dp
        rows = [1, 2, 3, 2]
        call lu%analyze(3, ptr, rows, values, status)
        call check(status%code == FORTSPARSE_INVALID_MATRIX, 'unsorted rows refused')
    end subroutine check_refusals

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (condition) return
        print '(a)', 'FAIL: '//message
        error stop 1
    end subroutine check
end program test_fortsparse_slu_analysis
