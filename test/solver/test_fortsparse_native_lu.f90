program test_fortsparse_native_lu
    use fortsparse_kinds, only: dp,i8
    use fortsparse_csc, only: csc_t,csc_z_t,csc_from_triplet,csc_matvec
    use fortsparse_native_lu, only: lu_factor_t,LU_NORMAL,LU_TRANSPOSE,LU_ADJOINT
    use fortsparse_native_lu_z, only: lu_factor_z_t
    use fortsparse_status
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
    implicit none
    call dense_real()
    call dense_complex()
    call branched_graphs()
    call scaling(20000)
    call sparse_complex(4096)
    call raw_errors()
    call limits()
    print '(a)', 'Native sparse real/complex LU values, operations, pivoting, reuse and memory pass'
contains
    subroutine checked(status)
        type(fortsparse_status_t), intent(in) :: status
        if(status%code/=FORTSPARSE_OK) then
            print *,status%code,trim(status%msg)
            error stop 'unexpected LU status'
        end if
    end subroutine checked
    subroutine real_csc(matrix,a)
        real(dp), intent(in) :: matrix(:,:)
        type(csc_t), intent(out) :: a
        integer :: i,j,p,n
        integer, allocatable :: rows(:),cols(:)
        real(dp), allocatable :: vals(:)
        type(fortsparse_status_t) :: status
        n=size(matrix,1)
        allocate(rows(count(matrix/=0)),cols(count(matrix/=0)),vals(count(matrix/=0)))
        p=0
        do j=1,n
            do i=1,n
                if(matrix(i,j)==0) cycle
                p=p+1;rows(p)=i;cols(p)=j;vals(p)=matrix(i,j)
            end do
        end do
        call csc_from_triplet(n,n,rows,cols,vals,a,status)
        call checked(status)
    end subroutine real_csc
    subroutine complex_csc(matrix,a)
        complex(dp), intent(in) :: matrix(:,:)
        type(csc_z_t), intent(out) :: a
        integer :: i,j,p,n
        integer, allocatable :: rows(:),cols(:)
        complex(dp), allocatable :: vals(:)
        type(fortsparse_status_t) :: status
        n=size(matrix,1)
        allocate(rows(count(matrix/=0)),cols(count(matrix/=0)),vals(count(matrix/=0)))
        p=0
        do j=1,n
            do i=1,n
                if(matrix(i,j)==0) cycle
                p=p+1;rows(p)=i;cols(p)=j;vals(p)=matrix(i,j)
            end do
        end do
        call csc_from_triplet(n,n,rows,cols,vals,a,status)
        call checked(status)
    end subroutine complex_csc
    subroutine dense_real()
        integer, parameter :: n=7, q(n)=[3,1,7,2,5,6,4]
        real(dp) :: l(n,n),u(n,n),matrix(n,n),base(n,n),reference(n,4),rhs(n,4),b(n,4)
        type(csc_t) :: a
        type(lu_factor_t) :: solver,small
        type(fortsparse_status_t) :: status
        integer :: i,j,op,repeat
        integer(i8) :: bytes
        l=0;u=0
        do j=1,n
            l(j,j)=1
            u(j,j)=2+.125_dp*j
            do i=j+1,n
                l(i,j)=.125_dp*sin(real(i+3*j,dp))
                u(j,i)=.25_dp*cos(real(i-2*j,dp))
            end do
        end do
        base=matmul(l,u)
        do i=1,n
            matrix(i,:)=base(mod(i,n)+1,:)
        end do
        matrix(1,1)=0 ! Exact zero diagonal requires numerical row pivoting.
        do j=1,4
            do i=1,n
                reference(i,j)=sin(real(i,dp))+.25_dp*j
            end do
        end do
        do repeat=1,2
            if(repeat==2) matrix(2,1)=matrix(2,1)+5
            call real_csc(matrix,a)
            call solver%factor(a,status,q,pivot_threshold=.5_dp)
            call checked(status)
            if(repeat==2.and..not.solver%last_reused_pattern) error stop 'LU pattern reuse lost'
            do op=LU_NORMAL,LU_ADJOINT
                if(op==LU_NORMAL) then
                    b=matmul(matrix,reference)
                else
                    b=matmul(transpose(matrix),reference)
                end if
                rhs=b
                bytes=solver%memory_bytes()
                call solver%solve(rhs,status,op)
                call checked(status)
                if(solver%memory_bytes()/=bytes) error stop 'LU solve allocated owned storage'
                if(maxval(abs(rhs-reference))>2e-11_dp) error stop 'real LU known block solution'
                if(op==LU_NORMAL) then
                    if(maxval(abs(matmul(matrix,rhs)-b))>2e-11_dp) error stop 'real LU residual'
                else
                    if(maxval(abs(matmul(transpose(matrix),rhs)-b))>2e-11_dp) error stop 'real LU transpose residual'
                end if
            end do
        end do
        call small%factor(a,status,q,low_memory=.true.)
        call checked(status)
        if(small%memory_bytes()>solver%memory_bytes()) error stop 'low-memory LU grew more storage'
        if(small%peak_memory_bytes<small%memory_bytes()) error stop 'LU peak below current'
        a%val=0
        call solver%factor(a,status,q)
        if(status%code/=FORTSPARSE_SINGULAR) error stop 'real LU singular failure'
        call solver%solve(rhs,status)
        if(status%code/=FORTSPARSE_NOT_FACTORED) error stop 'LU stale failed factor'
        call real_csc(matrix,a)
        call solver%factor(a,status,q)
        call checked(status)
        rhs=matmul(matrix,reference)
        call solver%solve(rhs,status)
        call checked(status)
        if(maxval(abs(rhs-reference))>2e-11_dp) error stop 'LU failure recovery'
        call solver%factor(a,status,[1,1,3,4,5,6,7])
        if(status%code/=FORTSPARSE_INVALID_MATRIX) error stop 'LU duplicate permutation'
        a%val(1)=ieee_value(1.0_dp,ieee_quiet_nan)
        call solver%factor(a,status)
        if(status%code/=FORTSPARSE_INVALID_MATRIX) error stop 'LU nonfinite input'
        call solver%free()
        call small%free()
        if(solver%memory_bytes()/=0) error stop 'LU free retained payload'
    end subroutine dense_real
    subroutine dense_complex()
        integer, parameter :: n=5,q(n)=[5,2,4,1,3]
        complex(dp) :: matrix(n,n),reference(n,3),rhs(n,3),b(n,3)
        type(csc_z_t) :: a
        type(lu_factor_z_t) :: solver
        type(fortsparse_status_t) :: status
        integer :: i,j,op,repeat
        do j=1,n
            do i=1,n
                matrix(i,j)=cmplx(.125_dp*sin(real(i+2*j,dp)),.25_dp*cos(real(2*i-j,dp)),dp)
            end do
            matrix(mod(j,n)+1,j)=matrix(mod(j,n)+1,j)+3
        end do
        matrix(1,1)=0
        do j=1,3
            do i=1,n
                reference(i,j)=cmplx(sin(real(i+j,dp)),cos(real(i-2*j,dp)),dp)
            end do
        end do
        do repeat=1,2
            if(repeat==2) matrix(2,1)=matrix(2,1)+cmplx(-4,.5_dp,dp)
            call complex_csc(matrix,a)
            call solver%factor(a,status,q,pivot_threshold=.1_dp,low_memory=.true.)
            call checked(status)
            do op=LU_NORMAL,LU_ADJOINT
                select case(op)
                case(LU_NORMAL)
                    b=matmul(matrix,reference)
                case(LU_TRANSPOSE)
                    b=matmul(transpose(matrix),reference)
                case(LU_ADJOINT)
                    b=matmul(conjg(transpose(matrix)),reference)
                end select
                rhs=b
                call solver%solve(rhs,status,op)
                call checked(status)
                if(maxval(abs(rhs-reference))>2e-11_dp) error stop 'complex LU known solution'
                select case(op)
                case(LU_NORMAL)
                    b=matmul(matrix,rhs)-b
                case(LU_TRANSPOSE)
                    b=matmul(transpose(matrix),rhs)-b
                case(LU_ADJOINT)
                    b=matmul(conjg(transpose(matrix)),rhs)-b
                end select
                if(maxval(abs(b))>2e-11_dp) error stop 'complex LU independent residual'
            end do
        end do
        ! Strongly unequal finite row scales must preserve transpose/adjoint semantics.
        matrix(1,:)=matrix(1,:)*1e-180_dp
        matrix(2,:)=matrix(2,:)*1e180_dp
        call complex_csc(matrix,a)
        call solver%factor(a,status,q)
        call checked(status)
        rhs=matmul(matrix,reference)
        call solver%solve(rhs,status)
        call checked(status)
        if(maxval(abs(rhs-reference))>2e-11_dp) error stop 'complex scaled LU'
        call solver%free()
    end subroutine dense_complex
    subroutine branched_graphs()
        integer, parameter :: n=13
        real(dp) :: base(n,n),matrix(n,n),expected(n),rhs(n)
        complex(dp) :: zmatrix(n,n),zexpected(n),zrhs(n)
        type(csc_t) :: a
        type(csc_z_t) :: az
        type(lu_factor_t) :: solver
        type(lu_factor_z_t) :: zsolver
        type(fortsparse_status_t) :: status
        integer :: trial,i,j,op,q(n)
        do trial=1,6
            base=0
            do j=1,n
                base(j,j)=6+.1_dp*trial
                do i=1,n
                    if(i==j.or.mod(3*i+5*j+trial,7)>1) cycle
                    base(i,j)=.2_dp*sin(real(i+2*j+trial,dp))
                end do
                expected(j)=sin(real(j,dp))
                zexpected(j)=cmplx(expected(j),cos(real(j,dp)),dp)
                q(j)=mod(5*j+trial,n)+1
            end do
            do i=1,n
                matrix(i,:)=base(mod(3*i+trial,n)+1,:)
                do j=1,n
                    zmatrix(i,j)=cmplx(matrix(i,j),.1_dp*matrix(i,j)*sin(real(i+j,dp)),dp)
                end do
            end do
            call real_csc(matrix,a)
            call complex_csc(zmatrix,az)
            call solver%factor(a,status,q)
            call checked(status)
            call zsolver%factor(az,status,q)
            call checked(status)
            do op=LU_NORMAL,LU_ADJOINT
                if(op==LU_NORMAL) then
                    rhs=matmul(matrix,expected)
                    zrhs=matmul(zmatrix,zexpected)
                else
                    rhs=matmul(transpose(matrix),expected)
                    if(op==LU_TRANSPOSE) then
                        zrhs=matmul(transpose(zmatrix),zexpected)
                    else
                        zrhs=matmul(conjg(transpose(zmatrix)),zexpected)
                    end if
                end if
                call solver%solve(rhs,status,op)
                call checked(status)
                call zsolver%solve(zrhs,status,op)
                call checked(status)
                if(maxval(abs(rhs-expected))>1e-11_dp) error stop 'branched real sparse LU'
                if(maxval(abs(zrhs-zexpected))>1e-11_dp) error stop 'branched complex sparse LU'
            end do
        end do
        call solver%free()
        call zsolver%free()
    end subroutine branched_graphs
    subroutine scaling(n)
        integer, intent(in) :: n
        type(csc_t) :: a
        type(lu_factor_t) :: solver
        type(fortsparse_status_t) :: status
        integer, allocatable :: row(:),column(:)
        real(dp), allocatable :: value(:),expected(:),rhs(:),b(:)
        integer :: j
        allocate(row(2*n),column(2*n),value(2*n),expected(n),rhs(n),b(n))
        do j=1,n
            row(2*j-1)=j;column(2*j-1)=j;value(2*j-1)=4
            row(2*j)=mod(j,n)+1;column(2*j)=j;value(2*j)=-.25_dp
            expected(j)=sin(.01_dp*j)
        end do
        call csc_from_triplet(n,n,row,column,value,a,status)
        call checked(status)
        call solver%factor_raw(n,a%col_ptr,a%row_idx,a%val,status)
        call checked(status)
        if(solver%memory_bytes()>400_i8*n) error stop 'LU sparse memory is not linear'
        if(solver%factor_nnz>5_i8*n) error stop 'LU cycle unexpectedly dense factors'
        b=csc_matvec(a,expected)
        rhs=b
        call solver%solve(rhs,status)
        call checked(status)
        if(maxval(abs(rhs-expected))>2e-12_dp) error stop 'LU sparse scaling known solution'
        if(maxval(abs(csc_matvec(a,rhs)-b))>2e-12_dp) error stop 'LU sparse scaling residual'
        call solver%free()
    end subroutine scaling
    subroutine sparse_complex(n)
        integer, intent(in) :: n
        type(csc_z_t) :: a
        type(lu_factor_z_t) :: solver
        type(fortsparse_status_t) :: status
        integer, allocatable :: row(:),column(:)
        complex(dp), allocatable :: value(:),expected(:),rhs(:),b(:)
        integer :: j
        allocate(row(2*n),column(2*n),value(2*n),expected(n),rhs(n),b(n))
        do j=1,n
            row(2*j-1)=j;column(2*j-1)=j;value(2*j-1)=cmplx(4,.5_dp,dp)
            row(2*j)=mod(j,n)+1;column(2*j)=j;value(2*j)=cmplx(-.25_dp,.1_dp,dp)
            expected(j)=cmplx(sin(.01_dp*j),cos(.03_dp*j),dp)
        end do
        call csc_from_triplet(n,n,row,column,value,a,status)
        call checked(status)
        call solver%factor_raw(n,a%col_ptr,a%row_idx,a%val,status,low_memory=.true.)
        call checked(status)
        if(solver%memory_bytes()>600_i8*n) error stop 'complex LU sparse memory not linear'
        b=(0.0_dp,0.0_dp)
        ! Independent original-A adjoint product, avoiding the CSC helper.
        do j=1,n
            b(j)=conjg(cmplx(4,.5_dp,dp))*expected(j) &
                +conjg(cmplx(-.25_dp,.1_dp,dp))*expected(mod(j,n)+1)
        end do
        rhs=b
        call solver%solve(rhs,status,LU_ADJOINT)
        call checked(status)
        if(maxval(abs(rhs-expected))>3e-12_dp) error stop 'complex raw sparse adjoint'
        call solver%free()
    end subroutine sparse_complex
    subroutine raw_errors()
        type(lu_factor_t) :: solver
        type(lu_factor_z_t) :: zsolver
        type(fortsparse_status_t) :: status
        real(dp) :: rhs(2)
        complex(dp) :: zrhs(3),zdiag(3),expected(3)
        call solver%factor_raw(2,[1,2,3],[2,1],[1._dp,1._dp],status,pivot_threshold=0._dp)
        call checked(status)
        rhs=[3._dp,2._dp]
        call solver%solve(rhs,status)
        call checked(status)
        if(maxval(abs(rhs-[2._dp,3._dp]))>1e-14_dp) error stop 'LU zero diagonal adversary'
        call solver%factor_raw(2,[1,3,5],[1,2,1,2],[1._dp,2._dp,2._dp,4._dp],status)
        if(status%code/=FORTSPARSE_SINGULAR) error stop 'LU elimination singularity'
        call solver%factor_raw(2,[1,2,4],[1,2],[1._dp,1._dp],status)
        if(status%code/=FORTSPARSE_INVALID_MATRIX) error stop 'LU raw pointer endpoint'
        call solver%factor_raw(2,[1,2,3],[1,3],[1._dp,1._dp],status)
        if(status%code/=FORTSPARSE_INVALID_MATRIX) error stop 'LU raw row bounds'
        call solver%factor_raw(2,[1,3,3],[1,1],[1._dp,1._dp],status)
        if(status%code/=FORTSPARSE_INVALID_MATRIX) error stop 'LU raw duplicate row'
        call solver%factor_raw(2,[1,3,5],[1,2,1,2], &
            [1e-300_dp,1._dp,1e300_dp,1._dp],status)
        if(status%code/=FORTSPARSE_SINGULAR) error stop 'LU scaling underflow status'
        call solver%factor_raw(0,[1],[integer ::],[real(dp) ::],status)
        call checked(status)
        call solver%solve(rhs(:0),status)
        call checked(status)
        zdiag=[cmplx(1e-280_dp,2e-280_dp,dp),cmplx(2e280_dp,-1e280_dp,dp),cmplx(.5_dp,2,dp)]
        expected=[cmplx(1,.5_dp,dp),cmplx(-.25_dp,.75_dp,dp),cmplx(.5_dp,-.25_dp,dp)]
        call zsolver%factor_raw(3,[1,2,3,4],[1,2,3],zdiag,status)
        call checked(status)
        zrhs=conjg(zdiag)*expected
        call zsolver%solve(zrhs,status,LU_ADJOINT)
        call checked(status)
        if(maxval(abs(zrhs-expected))>1e-13_dp) error stop 'scaled complex adjoint LU'
        call solver%free()
        call zsolver%free()
    end subroutine raw_errors
    subroutine limits()
        type(csc_t) :: a
        type(lu_factor_t) :: solver
        type(fortsparse_status_t) :: status
        integer(i8) :: fixed,budget
        real(dp) :: rhs(3)
        call csc_from_triplet(3,3,[1,2,3],[1,2,3],[2._dp,3._dp,4._dp],a,status)
        call checked(status)
        fixed=int(storage_size(0)/8,i8)*(10_i8*3+3+3)+16_i8*3
        budget=fixed+int(storage_size(0)/8,i8)+8_i8-1
        call solver%factor(a,status,memory_limit_bytes=budget)
        if(status%code/=FORTSPARSE_RESOURCE_LIMIT) error stop 'LU growth cap not enforced'
        if(solver%memory_bytes()>budget.or.solver%peak_memory_bytes>budget) error stop 'LU cap exceeded'
        if(solver%failed_required_bytes/=int(storage_size(0)/8,i8)+8_i8) error stop 'LU growth request bytes'
        if(solver%failed_available_bytes/=budget-fixed) error stop 'LU growth available bytes'
        call solver%factor(a,status,memory_limit_bytes=1_i8)
        if(status%code/=FORTSPARSE_RESOURCE_LIMIT) error stop 'LU fixed cap not enforced'
        if(solver%failed_required_bytes/=fixed.or.solver%failed_available_bytes/=1_i8) &
            error stop 'LU fixed budget diagnostic'

        call solver%factor(a,status,memory_limit_bytes=-1_i8)
        if(status%code/=FORTSPARSE_INVALID_MATRIX) error stop 'LU negative cap'
        if(solver%failed_required_bytes/=0.or.solver%failed_available_bytes/=0) &
            error stop 'LU stale budget diagnostic'

        call solver%factor(a,status,pivot_threshold=ieee_value(1._dp,ieee_quiet_nan))
        if(status%code/=FORTSPARSE_INVALID_MATRIX) error stop 'LU NaN threshold'
        call solver%factor(a,status,pivot_threshold=0._dp)
        call checked(status)
        rhs=[2._dp,6._dp,12._dp]
        call solver%solve(rhs,status)
        call checked(status)
        if(maxval(abs(rhs-[1._dp,2._dp,3._dp]))>1e-13_dp) error stop 'LU threshold zero'
        call solver%solve(rhs,status,99)
        if(status%code/=FORTSPARSE_INVALID_MATRIX) error stop 'LU invalid operation'
        rhs(1)=ieee_value(1._dp,ieee_quiet_nan)
        call solver%solve(rhs,status)
        if(status%code/=FORTSPARSE_INVALID_MATRIX) error stop 'LU nonfinite RHS'
        call solver%free()
    end subroutine limits
end program test_fortsparse_native_lu
