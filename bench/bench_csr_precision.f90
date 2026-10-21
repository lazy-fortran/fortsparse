program bench_csr_precision
    use fortsparse_kinds, only: dp,i8
    use fortsparse_precision_csr, only: precision_csr_t
    use fortnum_precision, only: FP64,FP32
    implicit none
    integer, parameter :: n=262144
    type(precision_csr_t) :: a
    real(dp), allocatable :: values(:),x(:),y(:),reference(:)
    integer, allocatable :: row(:),column(:)
    integer :: i,k,storage,sample,replay
    integer(i8) :: start,finish,rate
    real(dp) :: cpu0,cpu1,error
    logical :: ok
    allocate(row(n+1),column(3*n-2),values(3*n-2),x(n),y(n),reference(n))
    x=[(real(mod(i,257),dp)/129,i=1,n)]
    k=1
    do i=1,n
        row(i)=k
        if(i>1) then
            column(k)=i-1
            values(k)=-0.37_dp
            k=k+1
        end if
        column(k)=i
        values(k)=2.1_dp
        k=k+1
        if(i<n) then
            column(k)=i+1
            values(k)=-0.21_dp
            k=k+1
        end if
    end do
    row(n+1)=k
    reference=2.1_dp*x
    reference(2:)=reference(2:)-0.37_dp*x(:n-1)
    reference(:n-1)=reference(:n-1)-0.21_dp*x(2:)
    call system_clock(count_rate=rate)
    print '(a)', 'storage,sample,wall_s,cpu_s,max_error,payload_bytes'
    do storage=FP64,FP32
        call a%store(n,row,column,values,storage,ok)
        if(.not.ok) error stop 'store'
        do sample=-1,7
            call system_clock(start)
            call cpu_time(cpu0)
            do replay=1,4
                call a%matvec(x,y,storage,ok)
                if(.not.ok) error stop 'matvec'
            end do
            call cpu_time(cpu1)
            call system_clock(finish)
            error=maxval(abs(y-reference))
            if(error>1e-6_dp) error stop 'independent SpMV oracle'
            if(sample<1) cycle
            print '(i0,a,i0,3(a,es16.8),a,i0)',storage,',',sample, &
                ',',real(finish-start,dp)/real(rate,dp)/4,',',(cpu1-cpu0)/4, &
                ',',error,',',a%bytes()
        end do
    end do
end program
