program test_fortsparse_native_refinement
    use fortsparse
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
    implicit none
    call real_lu()
    call correction_limit()
    call structured()
    call complex_indefinite()
    call memory_and_reuse()
    call low_cap_retry()
    print '(a)', 'Native CPU refinement original residuals, operations, reuse and memory pass'
contains
    subroutine check(ok,message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not.ok) error stop message
    end subroutine check
    subroutine checked(status)
        type(fortsparse_status_t), intent(in) :: status
        if (.not.status_ok(status)) then
            print *,status%code,trim(status%msg)
            error stop 'unexpected refinement status'
        end if
    end subroutine checked
    real(dp) function error_real(a,b,x) result(value)
        real(dp), intent(in) :: a(:,:),b(:),x(:)
        value=maxval(abs(b-matmul(a,x)))/(maxval(sum(abs(a),dim=2))*maxval(abs(x))+maxval(abs(b)))
    end function error_real
    real(dp) function error_complex(a,b,x) result(value)
        complex(dp), intent(in) :: a(:,:),b(:),x(:)
        value=maxval(abs(b-matmul(a,x)))/(maxval(sum(abs(a),dim=2))*maxval(abs(x))+maxval(abs(b)))
    end function error_complex
    subroutine real_lu()
        type(sparse_solver_t) :: solver
        type(csc_t) :: a,adot
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: stats,before,after
        real(dp) :: matrix(2,2),exact(2),rhs(2),x(2),poor(2),block(2,4),result(2,4)
        real(dp) :: dotx(2),zero(2),seed(2),bbar(2),abar(4)
        integer :: op,j
        matrix=reshape([1e-12_dp,1._dp,1._dp,1._dp],[2,2])
        exact=[1._dp,1._dp]
        call csc_from_triplet(2,2,[1,2,1,2],[1,1,2,2],reshape(matrix,[4]),a,status)
        call checked(status)
        solver%backend_id=FORTSPARSE_BACKEND_NATIVE
        solver%options%pivot_tolerance=0.0_dp
        call sparse_factor(solver,a,status)
        call checked(status)
        rhs=matmul(matrix,exact)
        call sparse_solve(solver,rhs,poor,status)
        call checked(status)
        call check(error_real(matrix,rhs,poor)>32*epsilon(1._dp),'adversary must require correction')
        solver%refine=.true.
        call sparse_factor(solver,a,status)
        call checked(status)
        do op=0,2
            if(op==0) then
                rhs=matmul(matrix,exact)
            else
                rhs=matmul(transpose(matrix),exact)
            end if
            call sparse_factor_stats(solver,before,status)
            call sparse_solve_op(solver,rhs,x,op,status)
            call checked(status)
            call sparse_factor_stats(solver,stats,status)
            call check(stats%refinement_steps>=1.and.stats%refinement_steps<=3,'actual real corrections')
            call check(stats%backward_error>=0.and.stats%backward_error<=32*epsilon(1._dp),'reported real backward error')
            call check(error_real(matrix,rhs,x)<=32*epsilon(1._dp),'independent real backward error')
            call check(maxval(abs(x-exact))<1e-11_dp,'corrected real known solution')
            do j=1,4
                block(:,j)=j*rhs
            end do
            call sparse_solve_block_op(solver,block,result,op,status)
            call checked(status)
            call check(maxval(abs(result(:,4)-4*exact))<1e-11_dp,'refined block RHS')
            call sparse_factor_stats(solver,after,status)
            call check(before%host_current_bytes==after%host_current_bytes,'solve retains fixed owned payload')
            call check(before%host_peak_bytes==after%host_peak_bytes,'solve adds no payload peak')
        end do
        rhs=matmul(matrix,exact)
        x=rhs
        call sparse_solve(solver,x,status)
        call checked(status)
        call check(maxval(abs(x-exact))<1e-11_dp,'refined inplace RHS saved before solve')
        adot=a
        adot%val=.125_dp*a%val
        zero=0
        call sparse_solve_jvp(solver,adot,x,zero,dotx,status)
        call checked(status)
        call check(maxval(abs(dotx+.125_dp*exact))<1e-11_dp,'refined implicit JVP')
        seed=[.3_dp,-.7_dp]
        call sparse_solve_vjp(solver,a,x,seed,bbar,abar,status,use_primal_factors=.true.)
        call checked(status)
        call check(error_real(transpose(matrix),seed,bbar)<=32*epsilon(1._dp),'refined retained transpose VJP')
        rhs(1)=ieee_value(1._dp,ieee_quiet_nan)
        call sparse_solve(solver,rhs,x,status)
        call check(.not.status_ok(status),'nonfinite refined RHS fails explicitly')
        call sparse_free(solver)
    end subroutine real_lu
    subroutine correction_limit()
        type(sparse_solver_t) :: solver
        type(csc_t) :: a
        type(fortsparse_status_t) :: status,info
        type(sparse_factor_stats_t) :: stats
        integer, parameter :: n=6
        integer :: i,j,p,rows(n*n),cols(n*n)
        real(dp) :: matrix(n,n),exact(n),rhs(n),x(n)
        p=0
        do j=1,n
            do i=1,n
                p=p+1;rows(p)=i;cols(p)=j
                matrix(i,j)=real(mod(i*17+j*31+13*13,19)-9,dp)/7
            end do
            matrix(j,j)=matrix(j,j)+4
        end do
        matrix(1,1)=2._dp**(-40)
        exact=[1._dp,-.5_dp,.25_dp,.75_dp,-1._dp,.5_dp]
        call csc_from_triplet(n,n,rows,cols,reshape(matrix,[n*n]),a,status)
        call checked(status)
        solver%backend_id=FORTSPARSE_BACKEND_NATIVE
        solver%options%pivot_tolerance=0
        solver%refine=.true.
        call sparse_factor(solver,a,status)
        call checked(status)
        rhs=matmul(matrix,exact)
        call sparse_solve(solver,rhs,x,status)
        call sparse_factor_stats(solver,stats,info)
        call checked(info)
        ! Some arithmetic schedules can converge faster. Every successful
        ! result must pass the independent original residual; a finite failure
        ! must expose the three-correction limit rather than accept the result.
        if (status_ok(status)) then
            call check(error_real(matrix,rhs,x)<=32*epsilon(1._dp),'limit fixture cannot accept poor result')
        else
            call check(status%code==FORTSPARSE_SINGULAR,'explicit correction-limit failure')
            call check(stats%refinement_steps==3,'three corrections maximum')
            call check(stats%backward_error>32*epsilon(1._dp),'failed target recorded')
        end if
        print '(a,i4,es14.5)', 'correction-limit steps/backward-error:',stats%refinement_steps,stats%backward_error
        call sparse_free(solver)
    end subroutine correction_limit
    subroutine structured()
        type(sparse_solver_t) :: solver
        type(csc_t) :: a
        type(csc_z_t) :: z
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: stats
        complex(dp) :: matrix(3,3),exact(3),rhs(3),x(3),block(3,2),result(3,2)
        complex(dp) :: dotx(3),zero(3)
        real(dp) :: rm(3,3),re(3),rb(3),rx(3)
        integer :: rows(9),cols(9),i,j,k,op,family
        k=0
        do j=1,3
            do i=1,3
                k=k+1;rows(k)=i;cols(k)=j
            end do
        end do
        exact=[cmplx(1,.2_dp,dp),cmplx(-.5_dp,.8_dp,dp),cmplx(.4_dp,-.1_dp,dp)]
        do family=1,2
            matrix=reshape([cmplx(5,0,dp),cmplx(.2_dp,.3_dp,dp),cmplx(-.1_dp,.2_dp,dp), &
                cmplx(.2_dp,-.3_dp,dp),cmplx(4,0,dp),cmplx(.4_dp,-.2_dp,dp), &
                cmplx(-.1_dp,-.2_dp,dp),cmplx(.4_dp,.2_dp,dp),cmplx(3,0,dp)],[3,3])
            solver%options%symmetry=FORTSPARSE_SYM_SPD
            if(family==1) then
                matrix(1,1)=cmplx(1e-12_dp,1e-12_dp,dp)
                matrix(2,1)=cmplx(1,.1_dp,dp)
                solver%options%symmetry=FORTSPARSE_SYM_GENERAL
                solver%options%pivot_tolerance=0
            end if
            call csc_from_triplet(3,3,rows,cols,reshape(matrix,[9]),z,status)
            call checked(status)
            solver%backend_id=FORTSPARSE_BACKEND_NATIVE
            solver%refine=.true.
            call sparse_factor(solver,z,status)
            call checked(status)
            do op=0,2
                select case(op)
                case(0)
                    rhs=matmul(matrix,exact)
                case(1)
                    rhs=matmul(transpose(matrix),exact)
                case(2)
                    rhs=matmul(conjg(transpose(matrix)),exact)
                end select
                call sparse_solve_op(solver,rhs,x,op,status)
                call checked(status)
                call sparse_factor_stats(solver,stats,status)
                call check(stats%backward_error>=0.and.stats%backward_error<=32*epsilon(1._dp),'complex error report')
                select case(op)
                case(0)
                    call check(error_complex(matrix,rhs,x)<=32*epsilon(1._dp),'complex original residual')
                case(1)
                    call check(error_complex(transpose(matrix),rhs,x)<=32*epsilon(1._dp),'complex transpose residual')
                case(2)
                    call check(error_complex(conjg(transpose(matrix)),rhs,x)<=32*epsilon(1._dp),'complex adjoint residual')
                end select
                call check(maxval(abs(x-exact))<2e-11_dp,'complex refined known solution')
                block(:,1)=rhs;block(:,2)=2*rhs
                call sparse_solve_block_op(solver,block,result,op,status)
                call checked(status)
                call check(maxval(abs(result(:,2)-2*exact))<3e-11_dp,'complex refined block')
            end do
            rhs=matmul(matrix,exact)
            call sparse_solve(solver,rhs,status)
            call checked(status)
            call check(maxval(abs(rhs-exact))<2e-11_dp,'complex refined inplace')
            zero=(0._dp,0._dp)
            block(:,1)=matmul(matrix,exact)
            call sparse_solve_jvp(solver,z,exact,zero,dotx,status)
            call checked(status)
            call check(maxval(abs(dotx+exact))<2e-11_dp,'complex refined JVP')
            call sparse_free(solver)
        end do
        rm=reshape([0._dp,1._dp,.2_dp,1._dp,0._dp,-.1_dp,.2_dp,-.1_dp,2._dp],[3,3])
        re=[1._dp,-.5_dp,.2_dp]
        call csc_from_triplet(3,3,rows,cols,reshape(rm,[9]),a,status)
        call checked(status)
        solver%options%symmetry=FORTSPARSE_SYM_INDEFINITE
        call sparse_factor(solver,a,status)
        call checked(status)
        do op=0,2
            rb=matmul(rm,re)
            call sparse_solve_op(solver,rb,rx,op,status)
            call checked(status)
            call sparse_factor_stats(solver,stats,status)
            call check(stats%backward_error>=0.and.stats%backward_error<=32*epsilon(1._dp),'LDLT refined error')
            call check(error_real(rm,rb,rx)<=32*epsilon(1._dp),'LDLT original residual')
        end do
        call sparse_free(solver)
    end subroutine structured
    subroutine complex_indefinite()
        type(sparse_solver_t) :: solver
        type(csc_z_t) :: a,changed
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: stats
        complex(dp) :: matrix(3,3),exact(3),rhs(3),x(3),block(3,2),result(3,2)
        integer :: i,j,p,op,family,rows(9),cols(9)
        integer(i8) :: cap
        logical :: reused,mask(9)
        character(16) :: expected_family
        p=0
        do j=1,3
            do i=1,3
                p=p+1;rows(p)=i;cols(p)=j
            end do
        end do
        exact=[cmplx(1,.3_dp,dp),cmplx(-.5_dp,.7_dp,dp),cmplx(.2_dp,-.1_dp,dp)]
        solver%backend_id=FORTSPARSE_BACKEND_NATIVE
        solver%refine=.true.
        do family=1,2
            matrix=(0._dp,0._dp)
            matrix(1,2)=cmplx(1,.5_dp,dp)
            matrix(1,3)=cmplx(.1_dp,-.2_dp,dp)
            matrix(2,3)=cmplx(-.2_dp,.1_dp,dp)
            matrix(3,3)=2
            do j=1,3
                do i=j+1,3
                    matrix(i,j)=matrix(j,i)
                    if(family==1) matrix(i,j)=conjg(matrix(j,i))
                end do
            end do
            if(family==1) then
                solver%options%symmetry=FORTSPARSE_SYM_HERMITIAN
                expected_family='LDLH_COMPLEX'
            else
                solver%options%symmetry=FORTSPARSE_SYM_INDEFINITE
                matrix(3,3)=cmplx(2,.4_dp,dp)
                expected_family='LDLT_COMPLEX'
            end if
            call csc_from_triplet(3,3,rows,cols,reshape(matrix,[9]),a,status)
            call checked(status)
            call sparse_factor(solver,a,status)
            call checked(status)
            call sparse_factor_stats(solver,stats,status)
            call checked(status)
            call check(stats%native_factor_family==expected_family,'complex indefinite native dispatch')
            call check(stats%factor_nnz>0.and.stats%host_current_bytes>0,'complex indefinite owned memory')
            cap=stats%host_peak_bytes
            solver%options%host_memory_limit_bytes=cap
            matrix=1.125_dp*matrix
            a%val=reshape(matrix,[9])
            call sparse_refactor(solver,a,status,reused)
            if (.not.status_ok(status)) then
                call check(status%code==FORTSPARSE_RESOURCE_LIMIT,'tight complex refactor budget fails explicitly')
                call check(.not.reused,'tight cap cannot report completed reuse')
            end if
            ! Original/checker storage is retained during a refactor's growth.
            ! Add its exact quota to the independently observed fresh peak.
            cap=cap+int(storage_size(0)/8,i8)*(4+9)+16_i8*(9+2*3)
            solver%options%host_memory_limit_bytes=cap
            call sparse_refactor(solver,a,status,reused)
            call checked(status)
            call check(reused,'complex indefinite retained pattern refactor')
            do op=0,2
                select case(op)
                case(0)
                    rhs=matmul(matrix,exact)
                case(1)
                    rhs=matmul(transpose(matrix),exact)
                case(2)
                    rhs=matmul(conjg(transpose(matrix)),exact)
                end select
                call sparse_solve_op(solver,rhs,x,op,status)
                call checked(status)
                call check(maxval(abs(x-exact))<2e-12_dp,'complex indefinite N/T/H oracle')
                call sparse_factor_stats(solver,stats,status)
                call check(stats%backward_error>=0.and.stats%backward_error<=32*epsilon(1._dp),'complex LDL refinement')
                call check(stats%host_peak_bytes<=cap,'complex indefinite refactor/growth budget')
                block(:,1)=rhs;block(:,2)=2*rhs
                call sparse_solve_block_op(solver,block,result,op,status)
                call checked(status)
                call check(maxval(abs(result(:,2)-2*exact))<3e-12_dp,'complex indefinite block oracle')
            end do
            solver%options%host_memory_limit_bytes=1_i8
            call sparse_refactor(solver,a,status,reused)
            call check(status%code==FORTSPARSE_RESOURCE_LIMIT,'complex indefinite capped refinement')
            call check(.not.reused,'failed complex refactor cannot report successful reuse')
            solver%options%host_memory_limit_bytes=cap+2000_i8
            mask=.not.((rows==1.and.cols==3).or.(rows==3.and.cols==1))
            matrix(1,3)=0;matrix(3,1)=0
            call csc_from_triplet(3,3,pack(rows,mask),pack(cols,mask), &
                pack(reshape(matrix,[9]),mask),changed,status)
            call checked(status)
            call sparse_refactor(solver,changed,status,reused)
            call checked(status)
            call check(.not.reused,'complex indefinite changed pattern full fallback')
            rhs=matmul(matrix,exact)
            call sparse_solve(solver,rhs,x,status)
            call checked(status)
            call check(maxval(abs(x-exact))<2e-12_dp,'complex changed-pattern residual refresh')
            if(family==2) then
                solver%options%symmetry=FORTSPARSE_SYM_HERMITIAN
                call sparse_factor(solver,a,status)
                call check(status%code==FORTSPARSE_INVALID_MATRIX,'Hermitian class cannot silently use general LU')
            end if
            call sparse_free(solver)
            solver%options%host_memory_limit_bytes=0_i8
        end do
    end subroutine complex_indefinite
    subroutine low_cap_retry()
        integer, parameter :: n=80
        type(sparse_solver_t) :: solver
        type(csc_t) :: a,tiny
        type(csc_z_t) :: z,ztiny
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: stats
        integer :: rows(n*n),cols(n*n),i,j,k,kind,op
        integer(i8) :: cap,fixed,quota,entry,value_bytes,lower,upper
        real(dp) :: matrix(n,n),reference(n),rhs(n),x(n)
        real(dp), pointer :: owned(:)
        complex(dp) :: zmatrix(n,n),zreference(n),zrhs(n),zx(n)
        logical :: reused
        k=0
        do j=1,n
            do i=1,n
                k=k+1;rows(k)=i;cols(k)=j
                matrix(i,j)=.01_dp*sin(real(i+2*j,dp))
                if(i==j) matrix(i,j)=real(n,dp)
                zmatrix(i,j)=cmplx(matrix(i,j),.003_dp*cos(real(2*i+j,dp)),dp)
            end do
        end do
        reference=[(sin(real(i,dp)),i=1,n)]
        zreference=cmplx(reference,.25_dp*reference,dp)
        call csc_from_triplet(n,n,rows,cols,reshape(matrix,[n*n]),a,status)
        call checked(status)
        call csc_from_triplet(n,n,rows,cols,reshape(zmatrix,[n*n]),z,status)
        call checked(status)
        call csc_from_triplet(1,1,[1],[1],[2._dp],tiny,status)
        call checked(status)
        call csc_from_triplet(1,1,[1],[1],[cmplx(2._dp,.1_dp,dp)],ztiny,status)
        call checked(status)
        lower=int(n,i8)*(n-1)/2;upper=int(n,i8)*(n+1)/2
        do kind=1,2
            solver%backend_id=FORTSPARSE_BACKEND_NATIVE
            solver%options%memory_policy=FORTSPARSE_MEMORY_LOW
            if(kind==1) then
                call sparse_factor(solver,tiny,status)
            else
                call sparse_factor(solver,ztiny,status)
            end if
            call checked(status)
            owned=>sparse_vector(solver,3)
            call check(associated(owned),'LOW retry retains user vector')
            owned=[7._dp,9._dp,11._dp]
            solver%refine=.true.
            value_bytes=8_i8*kind
            entry=int(storage_size(0)/8,i8)+value_bytes
            fixed=int(storage_size(0)/8,i8)*(11_i8*n+3_i8+n*n)+(2_i8*value_bytes+8_i8)*n
            quota=int(storage_size(0)/8,i8)*(n+1_i8+n*n)+value_bytes*(n*n+2_i8*n)
            cap=fixed+(lower+2*upper-n)*entry+quota+24_i8
            ! The checker is reserved before LU; one byte below the old exact
            ! factor schedule still fails, with the live user vector preserved.
            solver%options%host_memory_limit_bytes=cap-1
            if(kind==1) then
                call sparse_refactor(solver,a,status,reused)
            else
                call sparse_refactor(solver,z,status,reused)
            end if
            call check(status%code==FORTSPARSE_RESOURCE_LIMIT,'refined LOW retry one-byte-under cap')
            call check(all(owned==[7._dp,9._dp,11._dp]),'failed LOW retry preserves user vector')
            solver%options%host_memory_limit_bytes=cap
            if(kind==1) then
                call sparse_refactor(solver,a,status,reused)
            else
                call sparse_refactor(solver,z,status,reused)
            end if
            call checked(status)
            do op=FORTSPARSE_OP_NORMAL,FORTSPARSE_OP_ADJOINT
                if(kind==1) then
                    if(op==FORTSPARSE_OP_NORMAL) then
                        rhs=matmul(matrix,reference)
                    else
                        rhs=matmul(transpose(matrix),reference)
                    end if
                    call sparse_solve_op(solver,rhs,x,op,status)
                    call checked(status)
                    call check(maxval(abs(x-reference))<3e-13_dp,'refined LOW real independent solution')
                    if(op==FORTSPARSE_OP_NORMAL) then
                        call check(error_real(matrix,rhs,x)<=32*epsilon(1._dp),'LOW real original residual')
                    else
                        call check(error_real(transpose(matrix),rhs,x)<=32*epsilon(1._dp),'LOW real transpose residual')
                    end if
                else
                    select case(op)
                    case(FORTSPARSE_OP_NORMAL)
                        zrhs=matmul(zmatrix,zreference)
                    case(FORTSPARSE_OP_TRANSPOSE)
                        zrhs=matmul(transpose(zmatrix),zreference)
                    case(FORTSPARSE_OP_ADJOINT)
                        zrhs=matmul(conjg(transpose(zmatrix)),zreference)
                    end select
                    call sparse_solve_op(solver,zrhs,zx,op,status)
                    call checked(status)
                    call check(maxval(abs(zx-zreference))<3e-13_dp,'refined LOW complex independent solution')
                    select case(op)
                    case(FORTSPARSE_OP_NORMAL)
                        call check(error_complex(zmatrix,zrhs,zx)<=32*epsilon(1._dp),'LOW complex original residual')
                    case(FORTSPARSE_OP_TRANSPOSE)
                        call check(error_complex(transpose(zmatrix),zrhs,zx)<=32*epsilon(1._dp), &
                            'LOW complex transpose residual')
                    case(FORTSPARSE_OP_ADJOINT)
                        call check(error_complex(conjg(transpose(zmatrix)),zrhs,zx)<=32*epsilon(1._dp), &
                            'LOW complex adjoint residual')
                    end select
                end if
                call sparse_factor_stats(solver,stats,status)
                call checked(status)
                call check(stats%host_peak_bytes<=cap,'refined LOW retry peak honors total cap')
                call check(stats%backward_error>=0.and.stats%backward_error<=32*epsilon(1._dp), &
                    'refined LOW retry original-matrix error')
            end do
            call check(all(owned==[7._dp,9._dp,11._dp]),'successful LOW retry preserves user vector')
            call sparse_free(solver)
            nullify(owned)
        end do
    end subroutine low_cap_retry
    subroutine memory_and_reuse()
        type(sparse_solver_t) :: solver
        type(csc_t) :: a,changed
        type(fortsparse_status_t) :: status
        type(sparse_factor_stats_t) :: plain,refined,stats
        real(dp) :: x(2),b(2)
        real(dp), pointer :: owned(:)
        integer(i8) :: extra,cap
        logical :: reused
        call csc_from_triplet(2,2,[1,2,1,2],[1,1,2,2],[3._dp,1._dp,.5_dp,2._dp],a,status)
        call checked(status)
        solver%backend_id=FORTSPARSE_BACKEND_NATIVE
        call sparse_factor(solver,a,status)
        call checked(status)
        call sparse_factor_stats(solver,plain,status)
        solver%refine=.true.
        call sparse_factor(solver,a,status)
        call checked(status)
        call sparse_factor_stats(solver,refined,status)
        extra=int(storage_size(0)/8,i8)*(3+4)+8_i8*(4+2*2)
        call check(refined%host_current_bytes-plain%host_current_bytes==extra,'all original CSC/scratch payload counted')
        owned=>sparse_vector(solver,2)
        call check(associated(owned),'retained user vector')
        owned=[7._dp,9._dp]
        call sparse_factor_stats(solver,stats,status)
        cap=stats%host_peak_bytes
        solver%options%host_memory_limit_bytes=cap
        a%val=a%val*1.25_dp
        call sparse_refactor(solver,a,status,reused)
        call checked(status)
        call check(reused,'refinement preserves same-pattern factor reuse')
        call check(all(owned==[7._dp,9._dp]),'refactor preserves live user vector')
        call sparse_factor_stats(solver,stats,status)
        call check(stats%host_peak_bytes<=cap,'refactor peak includes retained copies under cap')
        b=[4.375_dp,3.75_dp]
        call sparse_solve(solver,b,x,status)
        call checked(status)
        call check(maxval(abs(x-1._dp))<1e-12_dp,'refactor refreshes original values')
        ! The original matrix is retained, independently of the caller's values.
        a%val=0
        call sparse_solve(solver,b,x,status)
        call checked(status)
        call check(maxval(abs(x-1._dp))<1e-12_dp,'original values survive caller mutation')
        solver%options%host_memory_limit_bytes=extra
        call sparse_refactor(solver,a,status,reused)
        call check(status%code==FORTSPARSE_RESOURCE_LIMIT,'changed cap reserves refinement before factor')
        call check(all(owned==[7._dp,9._dp]),'failed cap keeps user vector lifetime')
        solver%options%host_memory_limit_bytes=cap+1000_i8
        call csc_from_triplet(2,2,[1,2],[1,2],[2._dp,3._dp],changed,status)
        call checked(status)
        call sparse_refactor(solver,changed,status,reused)
        call checked(status)
        call check(.not.reused,'changed pattern rebuilt under cap')
        b=[2._dp,3._dp]
        call sparse_solve(solver,b,x,status)
        call checked(status)
        call check(maxval(abs(x-1._dp))<1e-12_dp,'changed-pattern original refresh')
        call check(all(owned==[7._dp,9._dp]),'changed pattern keeps user vector lifetime')
        call sparse_factor_stats(solver,stats,status)
        call check(stats%host_peak_bytes<=cap+1000_i8,'changed-pattern peak honors cap')
        print '(a,3i12)', 'refinement payload baseline/extra/with-user:', &
            plain%host_current_bytes,extra,refined%host_current_bytes+16_i8
        call sparse_free(solver)
    end subroutine memory_and_reuse
end program test_fortsparse_native_refinement
