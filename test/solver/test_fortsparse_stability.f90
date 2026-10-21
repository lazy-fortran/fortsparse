program test_fortsparse_stability
    use fortsparse
    use stability_oracle, only: stability_pencil
    implicit none
    type(sparse_solver_t) :: solver
    type(fortsparse_status_t) :: status
    type(csc_t) :: matrix
    real(dp) :: a(3,3),m(3,3),p(3,3),d(3),sigma,x(3),b(3),known(3)
    integer :: trial,j,ordering,orders(2),rows(9),cols(9),k,i,state
    integer(i8) :: negative,zero,positive
    orders=[FORTSPARSE_ORDER_DEFAULT,FORTSPARSE_ORDER_AMD]
    k=0
    do j=1,3
        do i=1,3
            k=k+1
            rows(k)=i
            cols(k)=j
        end do
    end do
    known=[1.0_dp,-2.0_dp,3.0_dp]
    do ordering=1,2
        solver%backend_id=FORTSPARSE_BACKEND_NATIVE
        solver%options%symmetry=FORTSPARSE_SYM_INDEFINITE
        solver%options%ordering=orders(ordering)
        do trial=1,40
            ! Integer unit triangular congruences: eigenvalues of (A,M) are
            ! exactly d, proved symbolically by the fixture generator.
            d=[-real(trial+2,dp),real(trial,dp),real(2*trial+3,dp)]
            call stability_pencil(real(mod(trial,5)-2,dp), &
                real(mod(3*trial,7)-3,dp),real(mod(2*trial,5)-2,dp), &
                d(1),d(2),d(3),a(1,1),a(2,1),a(3,1),a(1,2),a(2,2), &
                a(3,2),a(1,3),a(2,3),a(3,3),m(1,1),m(2,1),m(3,1), &
                m(1,2),m(2,2),m(3,2),m(1,3),m(2,3),m(3,3))
            do j=0,3
                sigma=d(1)-1
                if (j==1) sigma=0.5_dp*(d(1)+d(2))
                if (j==2) sigma=0.5_dp*(d(2)+d(3))
                if (j==3) sigma=d(3)+1
                p=a-sigma*m
                call csc_from_triplet(3,3,rows,cols,reshape(p,[9]),matrix,status)
                call sparse_factor(solver,matrix,status)
                if (status%code==FORTSPARSE_BACKEND_UNAVAILABLE.and.ordering==2) exit
                call require(status%code==FORTSPARSE_OK,'pencil factor')
                call sparse_inertia(solver,negative,zero,positive,state,status)
                call require(status%code==FORTSPARSE_OK,'inertia query')
                call require(state==FORTSPARSE_INERTIA_RELIABLE,'separated inertia')
                call require(negative==int(j,i8).and.zero==0.and.positive==3-j,'exact pencil count')
                b=matmul(p,known)
                call sparse_solve(solver,b,x,status)
                call require(status%code==FORTSPARSE_OK,'retained pencil solve')
                call require(maxval(abs(x-known))<1e-10_dp,'manufactured pencil solution')
            end do
        end do
        call sparse_destroy(solver)
    end do
    call random_congruences()
    print *, 'FortSym-derived exact pencil inertia and manufactured solves PASS'
contains
    subroutine random_congruences()
        integer, parameter :: n=43
        real(dp) :: l(n,n),ld(n,n),a(n,n),permuted(n,n),diagonal(n),r
        integer :: i,j,trial,k,seed_size
        integer, allocatable :: seed(:)
        integer :: rows(n*n),cols(n*n),permutation(n)
        type(csc_t) :: matrix
        type(sparse_solver_t) :: factor
        type(fortsparse_status_t) :: status
        integer(i8) :: negative,zero,positive
        integer :: state
        call random_seed(size=seed_size)
        allocate(seed(seed_size))
        seed=20261021
        call random_seed(put=seed)
        k=0
        do j=1,n
            permutation(j)=mod(7*(j-1),n)+1
            do i=1,n
                k=k+1
                rows(k)=i
                cols(k)=j
            end do
        end do
        do trial=1,12
            l=0
            do j=1,n
                l(j,j)=1
                call random_number(r)
                diagonal(j)=real(1+int(8*r),dp)
                call random_number(r)
                if (trial>1.and.r<0.5_dp) diagonal(j)=-diagonal(j)
                do i=j+1,min(j+2,n)
                    call random_number(r)
                    l(i,j)=real(int(5*r)-2,dp)/4
                end do
                ld(:,j)=l(:,j)*diagonal(j)
            end do
            ! All entries are dyadic rationals, so assembly is exact here.
            ! det L=1 gives exactly the inertia of diagonal D by congruence.
            a=matmul(ld,transpose(l))
            permuted=a(permutation,permutation)
            call csc_from_triplet(n,n,rows,cols,reshape(permuted,[n*n]),matrix,status)
            factor%backend_id=FORTSPARSE_BACKEND_NATIVE
            factor%options%symmetry=FORTSPARSE_SYM_INDEFINITE
            if (trial==1) factor%options%symmetry=FORTSPARSE_SYM_SPD
            call sparse_factor(factor,matrix,status)
            call require(status%code==FORTSPARSE_OK,'random congruence factor')
            call sparse_inertia(factor,negative,zero,positive,state,status)
            call require(status%code==FORTSPARSE_OK,'random congruence query')
            call require(state==FORTSPARSE_INERTIA_RELIABLE,'random congruence separation')
            call require(negative==count(diagonal<0).and.zero==0.and. &
                positive==count(diagonal>0),'exact constructed random inertia')
            call sparse_destroy(factor)
        end do
    end subroutine random_congruences

    subroutine require(ok,message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not.ok) then
            print *, 'FAIL: ',message
            error stop 1
        end if
    end subroutine
end program
