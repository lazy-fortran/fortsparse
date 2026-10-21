program test_fortsparse_precision
    use iso_fortran_env, only: real64
    use fortnum_precision
    use fortsparse_precision_csr, only: precision_csr_t, csr_gmres_ir
    implicit none
    type(precision_csr_t) :: a, low
    type(precision_policy_t) :: policy
    integer, parameter :: n=33, nnz=3*n-2
    integer, parameter :: formats(4)=[FP64,FP32,FP16,BF16]
    integer :: row(n+1), column(nnz), i,k,f,info,steps,c
    real(real64) :: values(nnz), x(n), b(n), y(n), reference(n), residual
    logical :: ok
    k=1
    do i=1,n
        row(i)=k
        if(i>1) then
            column(k)=i-1
            values(k)=-0.37_real64
            k=k+1
        end if
        column(k)=i
        values(k)=2.1_real64
        k=k+1
        if(i<n) then
            column(k)=i+1
            values(k)=-0.21_real64
            k=k+1
        end if
    end do
    row(n+1)=k
    reference = [(sin(real(i,real64)),i=1,n)]
    ! Independent tridiagonal oracle, avoiding the CSR multiply under test.
    b = 2.1_real64*reference
    b(2:) = b(2:) - 0.37_real64*reference(:n-1)
    b(:n-1) = b(:n-1) - 0.21_real64*reference(2:)
    call a%store(n,row,column,values,FP64,ok)
    if(.not.ok) error stop 'exact CSR'
    do f=1,4
        call low%store(n,row,column,values,formats(f),ok)
        if(.not.ok) error stop 'compressed CSR'
        do c=FP64,FP32
            call low%matvec(reference,y,c,ok)
            if(.not.ok) error stop 'matvec'
            if(maxval(abs(y-b)) > &
                12*max(unit_roundoff(formats(f)),unit_roundoff(c))) &
                error stop 'SpMV roundoff'
        end do
        policy = precision_preset('cpu-mixed',ok)
        call policy%set(STAGE_PRECONDITION,formats(f),FP64,ok)
        x=0
        call csr_gmres_ir(a,b,x,policy,64*epsilon(1.0_real64),20,33,16, &
            info,steps,residual)
        if(info/=0) then
            print *, f, info, steps, residual
            error stop 'GMRES-IR convergence'
        end if
        if(maxval(abs(x-reference))>128*epsilon(1.0_real64)) &
            error stop 'GMRES-IR solution'
        call policy%set(STAGE_PRECONDITION,formats(f), &
            merge(FP64,FP32,formats(f)==FP64),ok)
        x=0
        call csr_gmres_ir(a,b,x,policy,64*epsilon(1.0_real64),20,33,16, &
            info,steps,residual)
        if(info/=0) error stop 'fp32 operator GMRES-IR'
    end do
    ! Large coefficients overflow raw binary16; scaled storage must not.
    call low%store(n,row,column,values*1.0e100_real64,FP16,ok)
    call low%matvec(reference,y,FP32,ok)
    if(.not.ok) error stop 'scaled half matvec'
    if(maxval(abs(y/1.0e100_real64-b)) > 12*unit_roundoff(FP16)) &
        error stop 'scaled half accuracy'
    row(2)=0
    call low%store(n,row,column,values,FP16,ok)
    if(ok) error stop 'invalid CSR accepted'
    print *, 'CSR storage/compute, original fp64 residual GMRES-IR PASS'
end program
