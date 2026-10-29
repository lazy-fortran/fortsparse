program ldlt_probe
    use fortsparse
    implicit none
    integer, parameter :: m = 70
    real(dp) :: shifts(6) = [0.0_dp, 0.05_dp, 0.3_dp, 1.0_dp, 2.0137_dp, 3.5_dp]
    integer :: s, o, ix, iy, k, c, ords(2)
    integer, allocatable :: rows(:), cols(:)
    real(dp), allocatable :: vals(:), b(:), x(:)
    type(csc_t) :: A
    type(sparse_solver_t) :: solver
    type(fortsparse_status_t) :: status
    type(sparse_factor_stats_t) :: stats
    ords = [FORTSPARSE_ORDER_NATURAL, FORTSPARSE_ORDER_AMD]
    allocate(rows(5*m*m), cols(5*m*m), vals(5*m*m), b(m*m), x(m*m))
    b = 1.0_dp
    do s = 1, size(shifts)
        k = 0
        do iy = 1, m
            do ix = 1, m
                c = (iy-1)*m + ix
                k = k+1; rows(k)=c; cols(k)=c; vals(k)=4.0_dp-shifts(s)
                if (ix>1) then; k=k+1; rows(k)=c; cols(k)=c-1; vals(k)=-1.0_dp; end if
                if (ix<m) then; k=k+1; rows(k)=c; cols(k)=c+1; vals(k)=-1.0_dp; end if
                if (iy>1) then; k=k+1; rows(k)=c; cols(k)=c-m; vals(k)=-1.0_dp; end if
                if (iy<m) then; k=k+1; rows(k)=c; cols(k)=c+m; vals(k)=-1.0_dp; end if
            end do
        end do
        call csc_from_triplet(m*m, m*m, rows(:k), cols(:k), vals(:k), A, status)
        do o = 1, 2
            solver%backend_id = FORTSPARSE_BACKEND_NATIVE
            solver%options%symmetry = FORTSPARSE_SYM_INDEFINITE
            solver%options%ordering = ords(o)
            call sparse_factor(solver, A, status)
            call sparse_factor_stats(solver, stats, status)
            print '(a,f7.4,a,i2,a,i9,a,f8.3)', 'shift ', shifts(s), ' order ', ords(o), &
                ' factor_nnz ', stats%factor_nnz, ' factor_s ', stats%factor_s
            call sparse_free(solver)
        end do
    end do
end program ldlt_probe
