program test_fortsparse_memory
    use fortsparse_kinds, only: dp, i8
    use fortsparse_status, only: fortsparse_status_t, FORTSPARSE_OK, &
        FORTSPARSE_RESOURCE_LIMIT, FORTSPARSE_INDEX_OVERFLOW, &
        FORTSPARSE_INVALID_MATRIX
    use fortsparse_memory, only: memory_budget_t, sparse_memory_plan_t, &
        memory_checked_add, memory_checked_mul, plan_natural_spd_memory, &
        MEMORY_PLAN_EXACT_LAYOUT
    use fortsparse_csc, only: csc_t
    use fortsparse_native, only: native_backend_t
    use fortsparse_backend, only: FORTSPARSE_SYM_SPD
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    type(memory_budget_t) :: budget
    type(sparse_memory_plan_t) :: plan
    type(fortsparse_status_t) :: status
    type(csc_t) :: empty, bad
    logical :: graph(17,17)
    integer(i8) :: value
    integer :: i, x, y

    budget%limit = 100_i8
    call budget%reserve(60_i8, status, 'factor')
    call expect(status%code == FORTSPARSE_OK, 'reservation succeeds')
    call budget%reserve(41_i8, status, 'second phase')
    call expect(status%code == FORTSPARSE_RESOURCE_LIMIT, 'cap refuses increment')
    call expect(budget%current == 60 .and. budget%peak == 60, 'failed live ledger')
    call expect(budget%failed_required == 41 .and. budget%failed_available == 40, &
        'required and available diagnostics')
    call expect(index(status%msg, 'second phase') > 0, 'phase diagnostic')
    call budget%release(61_i8, status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX, 'release underflow')
    call expect(budget%current == 60 .and. budget%peak == 60, 'release transactional')
    call budget%release(60_i8, status)
    call expect(status%code == FORTSPARSE_OK, 'release succeeds')
    call expect(budget%current == 0 .and. budget%peak == 60, 'peak survives release')
    call budget%reserve(-1_i8, status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX, 'negative reservation')
    budget%limit = -1_i8
    call budget%reserve(0_i8, status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX, 'negative limit')
    budget = memory_budget_t()
    call budget%reserve(huge(0_i8), status)
    call expect(status%code == FORTSPARSE_OK, 'unlimited supports byte range')
    call budget%reserve(1_i8, status)
    call expect(status%code == FORTSPARSE_INDEX_OVERFLOW, 'ledger addition overflow')
    call expect(budget%current == huge(0_i8), 'overflow ledger unchanged')
    value = 37_i8
    call memory_checked_add(huge(value), 1_i8, value, status)
    call expect(status%code == FORTSPARSE_INDEX_OVERFLOW .and. value == 37, &
        'checked addition leaves output on overflow')
    call memory_checked_mul(huge(value), 2_i8, value, status)
    call expect(status%code == FORTSPARSE_INDEX_OVERFLOW .and. value == 37, &
        'checked product leaves output on overflow')
    call memory_checked_mul(huge(value), 0_i8, value, status)
    call expect(status%code == FORTSPARSE_OK .and. value == 0, 'zero product')
    call memory_checked_add(-1_i8, 1_i8, value, status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX, 'negative byte arithmetic')

    ! Independent dense graph elimination oracle, including natural-order fill.
    graph = .false.
    call graph_case(graph(:1,:1))
    do i = 2, 9
        graph(1,i) = .true.
        graph(i,1) = .true.
    end do
    call graph_case(graph(:9,:9)) ! Center first: star fills a complete tail.
    graph = .true.
    do i = 1, 17
        graph(i,i) = .false.
    end do
    call graph_case(graph)
    graph = .false.
    do y = 1, 4
        do x = 1, 4
            i = x+4*(y-1)
            if (x < 4) then
                graph(i,i+1) = .true.
                graph(i+1,i) = .true.
            end if
            if (y < 4) then
                graph(i,i+4) = .true.
                graph(i+4,i) = .true.
            end if
        end do
    end do
    call graph_case(graph) ! Includes isolated node17 and nonzero graph fill.
    graph = .false.
    do i = 1, 13, 2
        graph(i,i+2) = .true.
        graph(i+2,i) = .true.
    end do
    call graph_case(graph) ! Interleaved disconnected components.

    empty%nrow = 0
    empty%ncol = 0
    empty%nnz = 0
    allocate (empty%col_ptr(1), empty%row_idx(0), empty%val(0))
    empty%col_ptr = 1
    call plan_natural_spd_memory(empty, plan, status, nrhs=0)
    call expect(status%code == FORTSPARSE_OK .and. plan%fits, 'empty CSC plan')
    call expect(plan%peak_host_bytes == 8 .and. plan%max_rhs_batch == 0, &
        'empty reference endpoints')
    bad = empty
    bad%col_ptr = 2
    call plan_natural_spd_memory(bad, plan, status)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX, 'invalid empty endpoint')
    call plan_natural_spd_memory(empty, plan, status, nrhs=-1)
    call expect(status%code == FORTSPARSE_INVALID_MATRIX, 'negative RHS count')
    bad = empty
    bad%ncol = huge(i)
    bad%nrow = bad%ncol
    call plan_natural_spd_memory(bad, plan, status)
    call expect(status%code == FORTSPARSE_INDEX_OVERFLOW, 'dimension overflow guard')
    print '(a)', 'Memory ledger, independent fill/layout plans and cap failures pass'
contains
    subroutine expect(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print '(a)', message
            error stop 1
        end if
    end subroutine expect

    subroutine graph_case(edges)
        logical, intent(in) :: edges(:,:)
        logical :: filled(size(edges,1),size(edges,1))
        type(csc_t) :: A, invalid
        type(native_backend_t) :: reference
        type(sparse_memory_plan_t) :: result
        type(memory_budget_t) :: cap
        type(fortsparse_status_t) :: s
        integer :: n, col, row, p, k
        integer(i8) :: nz, payload
        n = size(edges,1)
        A%nrow = n
        A%ncol = n
        A%nnz = n+count(edges)
        allocate (A%col_ptr(n+1), A%row_idx(A%nnz), A%val(A%nnz))
        p = 1
        do col = 1, n
            A%col_ptr(col) = p
            do row = 1, n
                if (.not. edges(row,col) .and. row /= col) cycle
                A%row_idx(p) = row
                A%val(p) = -1.0_dp
                if (row == col) A%val(p) = 2.0_dp+count(edges(:,col))
                p = p+1
            end do
        end do
        A%col_ptr(n+1) = p
        filled = edges
        nz = int(n,i8)
        do col = 1, n
            do row = col+1, n
                if (.not. filled(row,col)) cycle
                nz = nz+1
                do k = row+1, n
                    if (.not. filled(k,col)) cycle
                    filled(k,row) = .true.
                    filled(row,k) = .true.
                end do
            end do
        end do
        call plan_natural_spd_memory(A, result, s, nrhs=8)
        call expect(s%code == FORTSPARSE_OK .and. result%fits, 'graph plan succeeds')
        call expect(result%prediction_quality == MEMORY_PLAN_EXACT_LAYOUT, &
            'exact reference layout quality')
        call expect(result%symbolic_nnz == nz, 'independent graph elimination fill')
        call expect(result%factor_host_bytes == 8*nz, 'factor payload')
        call expect(result%persistent_host_bytes == 4*nz+4*A%nnz+36*n+8, &
            'retained reference arrays')
        call expect(result%peak_host_bytes == 12*nz+4*A%nnz+36*n+8, &
            'reference peak payload')
        call expect(result%workspace_host_bytes == 0 .and. result%peak_device_bytes == 0, &
            'no transient reference RHS or device workspace')
        call expect(result%max_rhs_batch == 8, 'reference batched solve reuses vector')
        reference%reference_only = .true.
        reference%options%symmetry = FORTSPARSE_SYM_SPD
        call reference%factor_real(A, .false., s)
        call expect(s%code == FORTSPARSE_OK, 'actual reference factor')
        payload = 8_i8*(size(reference%lv,kind=i8)+size(reference%work,kind=i8))
        payload = payload+4_i8*(size(reference%li,kind=i8)+ &
            size(reference%ap,kind=i8)+size(reference%ai,kind=i8)+ &
            size(reference%lp,kind=i8)+size(reference%parent,kind=i8)+ &
            size(reference%count,kind=i8)+size(reference%mark,kind=i8)+ &
            size(reference%reach,kind=i8)+size(reference%path,kind=i8))
        call expect(payload == result%peak_host_bytes, 'actual allocated reference layout')
        call reference%free()
        cap%current = 23
        cap%peak = 23
        cap%limit = 23+result%peak_host_bytes
        call plan_natural_spd_memory(A, result, s, cap, nrhs=3)
        call expect(s%code == FORTSPARSE_OK .and. result%fits, 'exact remaining cap')
        call expect(cap%current == 23, 'planner releases actual scratch')
        call expect(cap%peak == 23+12*n, 'prediction does not count as allocation')
        cap%limit = cap%limit-1
        call plan_natural_spd_memory(A, result, s, cap, nrhs=3)
        call expect(s%code == FORTSPARSE_RESOURCE_LIMIT .and. .not. result%fits, &
            'predicted peak includes existing reservations')
        call expect(result%max_rhs_batch == 0 .and. cap%current == 23, &
            'failed plan retains caller ledger')
        call expect(cap%failed_required == result%peak_host_bytes .and. &
            cap%failed_available == cap%limit-23, 'plan failure byte diagnostic')
        cap = memory_budget_t(limit=12_i8*n-1)
        call plan_natural_spd_memory(A, result, s, cap)
        call expect(s%code == FORTSPARSE_RESOURCE_LIMIT .and. cap%current == 0, &
            'analysis scratch itself obeys cap')
        invalid = A
        invalid%col_ptr(2) = A%nnz+2
        call plan_natural_spd_memory(invalid, result, s)
        call expect(s%code == FORTSPARSE_INVALID_MATRIX, 'invalid intermediate pointer')
        invalid = A
        invalid%val(1) = ieee_value(0.0_dp, ieee_quiet_nan)
        call plan_natural_spd_memory(invalid, result, s)
        call expect(s%code == FORTSPARSE_INVALID_MATRIX, 'nonfinite values rejected')
        if (n > 1) then
            if (.not. edges(2,1)) return
            invalid = A
            invalid%val(2) = invalid%val(2)+.5_dp
            call plan_natural_spd_memory(invalid, result, s)
            call expect(s%code == FORTSPARSE_INVALID_MATRIX, 'asymmetric values rejected')
            invalid = A
            invalid%row_idx(2) = 1
            call plan_natural_spd_memory(invalid, result, s)
            call expect(s%code == FORTSPARSE_INVALID_MATRIX, 'duplicate rows rejected')
        end if
    end subroutine graph_case
end program test_fortsparse_memory
