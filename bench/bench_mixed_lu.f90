program bench_mixed_lu
    use fortsparse
    use fortnum_precision, only: FP64,FP32
    implicit none
    integer, parameter :: n=600,band=9,nrhs=16
    type(csc_t) :: a
    type(mixed_lu_t) :: factor
    type(fortsparse_status_t) :: status,statuses(nrhs)
    type(mixed_solve_report_t) :: reports(nrhs)
    real(dp), allocatable :: values(:),b(:,:),x(:,:),known(:,:),dense_rhs(:)
    integer, allocatable :: rows(:),columns(:)
    integer :: i,j,k,storage,r
    integer(i8) :: start,finish,rate
    real(dp) :: cpu0,cpu1,wall,cpu,error
    allocate(rows(n*(2*band+1)),columns(n*(2*band+1)),values(n*(2*band+1)))
    allocate(b(n,nrhs),x(n,nrhs),known(n,nrhs),dense_rhs(n))
    do r=1,nrhs
        known(:,r)=[(sin(real(i+r,dp)),i=1,n)]
    end do
    b=0
    k=0
    do j=1,n
        do i=max(1,j-band),min(n,j+band)
            k=k+1
            rows(k)=i
            columns(k)=j
            values(k)=real(mod(7*i+3*j,11)-5,dp)/71
            if (i==j) values(k)=4+real(mod(i,31),dp)/31
            b(i,:)=b(i,:)+values(k)*known(j,:)
        end do
    end do
    call csc_from_triplet(n,n,rows(:k),columns(:k),values(:k),a,status)
    if(status%code/=FORTSPARSE_OK) error stop 'assemble'
    call system_clock(count_rate=rate)
    print '(a)', 'storage,phase,wall_s,cpu_s,max_error,payload_bytes,steps,rcond'
    do storage=FP64,FP32
        call cpu_time(cpu0)
        call system_clock(start)
        call factor%factor(a,status,storage)
        call system_clock(finish)
        call cpu_time(cpu1)
        if(status%code/=FORTSPARSE_OK) error stop 'factor'
        wall=real(finish-start,dp)/real(rate,dp)
        cpu=cpu1-cpu0
        call factor%solve(b(:,1),dense_rhs,reports(1),status)
        error=maxval(abs(dense_rhs-known(:,1)))
        if(status%code/=FORTSPARSE_OK.or.error>1e-12_dp) error stop 'warm solution'
        print '(i0,a,3(a,es16.8),a,i0,a,i0,a,es16.8)',reports(1)%storage,',factor', &
            ',',wall,',',cpu,',',error,',',factor%bytes(),',',reports(1)%steps, &
            ',',reports(1)%rcond_estimate
        call cpu_time(cpu0)
        call system_clock(start)
        call factor%solve_block(b,x,reports,statuses)
        call system_clock(finish)
        call cpu_time(cpu1)
        error=maxval(abs(x-known))
        if(any(statuses%code/=FORTSPARSE_OK).or.error>1e-12_dp) error stop 'block solution'
        print '(i0,a,3(a,es16.8),a,i0,a,i0,a,es16.8)',reports(1)%storage,',solve16', &
            ',',real(finish-start,dp)/real(rate,dp),',',cpu1-cpu0,',',error, &
            ',',factor%bytes(),',',maxval(reports%steps),',',reports(1)%rcond_estimate
        call factor%free()
    end do
end program
