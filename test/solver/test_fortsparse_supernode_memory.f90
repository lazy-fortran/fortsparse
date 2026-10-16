program test_fortsparse_supernode_memory
    use fortsparse_kinds, only: dp, i8
    use fortsparse_status, only: fortsparse_status_t, FORTSPARSE_OK, &
        FORTSPARSE_RESOURCE_LIMIT, FORTSPARSE_BACKEND_UNAVAILABLE, &
        FORTSPARSE_NOT_POSITIVE_DEFINITE, FORTSPARSE_NOT_FACTORED
    use fortsparse_memory, only: memory_budget_t
    use fortsparse_csc, only: csc_t
    use fortsparse_native_supernodes, only: supernodes_t
    implicit none
    logical :: openmp_enabled
    openmp_enabled = .false.
    !$ openmp_enabled = .true.
    call check_layout(1, 1, .true.)
    call check_layout(1, 65, .true.)
    call check_layout(1, 257, .true.)
    call check_layout(8, 130, .true.)
    print '(a)', 'Budgeted CPU panels, allocation lifetimes and numerical oracles pass'
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print '(a)', message
            error stop 1
        end if
    end subroutine require

    subroutine check_layout(blocks, width, numerical)
        integer, intent(in) :: blocks, width
        logical, intent(in) :: numerical
        type(csc_t) :: A
        type(supernodes_t) :: panels, uncapped
        type(memory_budget_t) :: budget, trial
        type(fortsparse_status_t) :: status
        integer, allocatable :: lp(:), li(:)
        real(dp), allocatable :: exact(:), rhs(:)
        integer :: n, col, row, last, first, p, q, threads, capped_threads
        integer(i8) :: retained, peak, baseline, cap, one_thread_peak
        logical :: reduced

        n = blocks*width
        allocate(lp(n+1), li(blocks*width*(width+1)/2), exact(n), rhs(n))
        A%nrow = n
        A%ncol = n
        A%nnz = blocks*width*width
        allocate(A%col_ptr(n+1), A%row_idx(A%nnz), A%val(A%nnz))
        p = 1
        q = 1
        do col = 1, n
            first = 1+((col-1)/width)*width
            last = first+width-1
            lp(col) = p
            do row = col, last
                li(p) = row
                p = p+1
            end do
            A%col_ptr(col) = q
            do row = first, last
                A%row_idx(q) = row
                A%val(q) = .01_dp
                if (row == col) A%val(q) = 2.01_dp
                q = q+1
            end do
            exact(col) = cos(real(col,dp))
        end do
        lp(n+1) = p
        A%col_ptr(n+1) = q
        baseline = 73_i8
        budget = memory_budget_t(current=baseline, peak=baseline)
        call panels%build(n, lp, li, 4, status, budget=budget)
        call require(status%code == FORTSPARSE_OK, 'unlimited accounted build')
        retained = panels%memory_bytes()
        peak = budget%peak
        threads = panels%active_threads
        if (blocks > 1 .and. openmp_enabled) call require(threads > 1, &
            'fixture exercises useful threaded scratch')
        call require(retained == budget%current-baseline, 'retained exact ledger')
        call require(retained == panels%factor_storage_bytes+panels%workspace_bytes, &
            'disjoint factor and other payload')
        call require(peak >= budget%current, 'setup overlap included in peak')
        call require(retained > panels%factor_storage_bytes+4_i8*n, &
            'descriptor and index payload included')
        if (numerical) then
            do first = 1, n, width
                last = first+width-1
                rhs(first:last) = 2.0_dp*exact(first:last)+.01_dp*sum(exact(first:last))
            end do
            call panels%factor(A, status)
            call require(status%code == FORTSPARSE_OK, 'budgeted dense factor oracle')
            call panels%solve(rhs, status)
            call require(status%code == FORTSPARSE_OK, 'budgeted solve')
            call require(maxval(abs(rhs-exact)) < 2.e-13_dp, &
                'independent known solution')
            A%val = 1.1_dp*A%val
            do first = 1, n, width
                last = first+width-1
                rhs(first:last) = 1.1_dp*(2.0_dp*exact(first:last) &
                    +.01_dp*sum(exact(first:last)))
            end do
            call panels%factor(A, status)
            call require(status%code == FORTSPARSE_OK, 'budgeted retained refactor')
            call panels%solve(rhs, status)
            call require(status%code == FORTSPARSE_OK, 'refactor solve')
            call require(maxval(abs(rhs-exact)) < 2.e-13_dp, 'refactor known solution')
            call require(budget%current == baseline+retained &
                .and. budget%peak == peak, &
                'factor and solve allocate no tracked payload')
            A%val(1) = -1.0_dp
            call panels%factor(A, status)
            call require(status%code == FORTSPARSE_NOT_POSITIVE_DEFINITE, &
                'strict factor reports failed pivot')
            call panels%solve(rhs, status)
            call require(status%code == FORTSPARSE_NOT_FACTORED, &
                'failed factor not reusable')
            A%val(1) = 1.1_dp*2.01_dp
            do first = 1, n, width
                last = first+width-1
                rhs(first:last) = 1.1_dp*(2.0_dp*exact(first:last) &
                    +.01_dp*sum(exact(first:last)))
            end do
            call panels%factor(A, status)
            call require(status%code == FORTSPARSE_OK, 'strict pivot recovery')
            call panels%solve(rhs, status)
            call require(status%code == FORTSPARSE_OK, 'recovered solve')
            call require(maxval(abs(rhs-exact)) < 2.e-13_dp, 'recovered known solution')
        end if
        call panels%free(budget, status)
        call require(status%code == FORTSPARSE_OK .and. budget%current == baseline, &
            'free releases only owned budget bytes')
        call panels%free(budget, status)
        call require(status%code == FORTSPARSE_OK .and. budget%current == baseline, &
            'idempotent budgeted free')
        budget%limit = peak
        call panels%build(n, lp, li, 4, status, budget=budget)
        call require(status%code == FORTSPARSE_OK, 'exact observed peak cap fits')
        call require(budget%peak <= budget%limit, 'admitted setup never exceeds cap')
        ! Replacing the same object frees its previous arrays before allocating.
        call panels%build(n, lp, li, 4, status, budget=budget)
        call require(status%code == FORTSPARSE_OK, 'same-ledger replacement')
        call require(budget%current == baseline+panels%memory_bytes(), &
            'replacement does not duplicate retained charge')
        call panels%free(budget, status)
        trial = memory_budget_t(limit=baseline+20_i8*n-1, &
            current=baseline, peak=baseline)
        call panels%build(n, lp, li, 4, status, budget=trial)
        call require(status%code == FORTSPARSE_RESOURCE_LIMIT, &
            'graph scratch preallocation cap')
        call require(trial%current == baseline .and. panels%memory_bytes() == 0, &
            'failed builder restores entry ledger')
        call require(trial%failed_required == 20_i8*n, &
            'failed phase requirement retained')
        ! Exercise failures after some graph, descriptor, row and value allocations.
        do cap = baseline+20_i8*n, peak-1_i8, max(1_i8,(peak-baseline)/7_i8)
            trial = memory_budget_t(limit=cap, current=baseline, peak=baseline)
            call panels%build(n, lp, li, 4, status, budget=trial)
            if (status%code == FORTSPARSE_OK) then
                call require(trial%peak <= cap, 'reduced schedule respects cap')
                call panels%free(trial, status)
            else
                call require(status%code == FORTSPARSE_RESOURCE_LIMIT, &
                    'bounded setup failure')
                call require(trial%current == baseline &
                    .and. panels%memory_bytes() == 0, &
                    'partial setup cleanup is transactional')
            end if
        end do
        if (threads > 1) then
            trial = memory_budget_t(current=baseline, peak=baseline)
            call panels%build(n, lp, li, 1, status, budget=trial)
            call require(status%code == FORTSPARSE_OK, 'single-thread bound')
            one_thread_peak = trial%peak
            call panels%free(trial, status)
            trial = memory_budget_t(limit=one_thread_peak, &
                current=baseline, peak=baseline)
            call panels%build(n, lp, li, 4, status, budget=trial)
            call require(status%code == FORTSPARSE_OK, &
                'workspace cap can reduce useful team')
            capped_threads = panels%active_threads
            reduced = panels%budget_reduced_threads
            call require(capped_threads < threads .and. reduced, &
                'actual reduced team reported')
            call panels%free(trial, status)
        end if
        call uncapped%build(n, lp, li, 4, status)
        call require(status%code == FORTSPARSE_OK, &
            'legacy unbudgeted build remains available')
        call require(uncapped%memory_bytes() > 0, 'unbudgeted descriptors attributed')
        call uncapped%free()
        trial = memory_budget_t(limit=peak, current=baseline, peak=baseline)
        call panels%build(n, lp, li, 4, status, device=0, budget=trial)
        call require(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
            'CPU budget rejects GPU gap')
        call require(trial%current == baseline, 'unsupported device preserves ledger')
    end subroutine check_layout
end program test_fortsparse_supernode_memory
