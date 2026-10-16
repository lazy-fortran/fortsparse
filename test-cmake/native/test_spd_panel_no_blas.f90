program test_spd_panel_no_blas
    use fortsparse_kinds, only: dp, i8
    use fortsparse_csc, only: csc_t
    use fortsparse_memory, only: memory_budget_t
    use fortsparse_status, only: fortsparse_status_t, status_ok
    use fortsparse_spd_leaves, only: spd_panel_solve
    use fortsparse_native_supernodes, only: supernodes_t
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    integer, parameter :: n=96, m=35, pad=3, full=257
    real(dp) :: l(n+pad,n), b(m+pad,n), truth(m,n), original(m,n)
    real(dp) :: dense(full,full), rhs(full), exact(full)
    integer :: i,j,k,before,lp(full+1),li(full*(full+1)/2)
    type(csc_t) :: a
    type(supernodes_t) :: panels
    type(memory_budget_t) :: budget
    type(fortsparse_status_t) :: status
    interface
        integer(c_int) function calls() bind(C,name='fsparse_test_dtrsm_calls')
            import c_int
        end function calls
    end interface
    l=-99.0_dp
    do j=1,n
        l(j,j)=2.0_dp+0.01_dp*j
        do i=j+1,n
            l(i,j)=0.001_dp*(i-j)
        end do
        do i=1,m
            truth(i,j)=sin(real(i+2*j,dp))
        end do
    end do
    ! Independent triangular product, using only the declared lower factor.
    original=0.0_dp
    do j=1,n
        do k=1,j
            original(:,j)=original(:,j)+truth(:,k)*l(j,k)
        end do
    end do
    b=-77.0_dp; b(:m,:)=original
    before=calls()
    call spd_panel_solve(l,n,n+pad,b,m,m+pad)
    call expect(calls()==before+1,'default wide leaf did not lower to DTRSM')
    call expect(maxval(abs(b(:m,:)-truth))<1e-12_dp,'default leaf values')
    call expect(all(b(m+1:,:)==-77.0_dp),'default leaf padding')
    b=-77.0_dp; b(:m,:)=original
    before=calls()
    call spd_panel_solve(l,n,n+pad,b,m,m+pad,use_blas=.false.)
    call expect(calls()==before,'disabled leaf called DTRSM')
    call expect(maxval(abs(b(:m,:)-truth))<1e-12_dp,'generated leaf values')
    call expect(all(b(m+1:,:)==-77.0_dp),'generated leaf padding')
    ! Dense width-128 fronts trigger the same lowering in a real sparse factor.
    a%nrow=full; a%ncol=full; a%nnz=full*full
    allocate(a%col_ptr(full+1),a%row_idx(a%nnz),a%val(a%nnz))
    dense=0.01_dp
    k=0
    do j=1,full
        dense(j,j)=2.0_dp
        exact(j)=sin(real(j,dp))
        a%col_ptr(j)=(j-1)*full+1
        a%row_idx((j-1)*full+1:j*full)=[(i,i=1,full)]
        a%val((j-1)*full+1:j*full)=dense(:,j)
        lp(j)=k+1
        do i=j,full
            k=k+1; li(k)=i
        end do
    end do
    a%col_ptr(full+1)=a%nnz+1; lp(full+1)=k+1
    budget%limit=10000000_i8
    call panels%build(full,lp,li,1,status,budget=budget)
    call expect(status_ok(status),'bounded sparse panel build')
    before=calls()
    call panels%factor(a,status)
    call expect(status_ok(status),'bounded sparse factor')
    call expect(calls()==before,'strict sparse factor called DTRSM')
    rhs=matmul(dense,exact)
    call panels%solve(rhs,status)
    call expect(status_ok(status),'bounded sparse solve')
    call expect(maxval(abs(rhs-exact))<1e-12_dp,'bounded original-matrix oracle')
    call panels%free(budget,status)
    call expect(status_ok(status).and.budget%current==0,'bounded free accounting')
    call panels%build(full,lp,li,1,status)
    call expect(status_ok(status),'ordinary sparse panel build')
    before=calls()
    call panels%factor(a,status)
    call expect(status_ok(status),'ordinary sparse factor')
    call expect(calls()>before,'ordinary wide sparse factor lost DTRSM lowering')
    call panels%free()
    print '(a)', 'Wide panel generated values, padding and strict/ordinary DTRSM dispatch pass'
contains
    subroutine expect(condition,message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not.condition) then
            print *,message
            error stop 'panel dispatch oracle'
        end if
    end subroutine expect
end program test_spd_panel_no_blas
