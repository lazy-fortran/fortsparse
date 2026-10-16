program test_fortsparse_lu_leaves
    use fortsparse_kinds, only: dp
    use fortsparse_lu_leaves
    implicit none
    real(dp) :: x(6),dx(6),v(2),dv(2),seed(2),bar(6),plus(2),minus(2),h
    complex(dp) :: a,b,c,expected,actual
    integer :: i
    x=[.7_dp,-.2_dp,1.3_dp,.4_dp,-.8_dp,.3_dp]
    dx=[.2_dp,-.4_dp,.1_dp,.3_dp,-.2_dp,.5_dp]
    seed=[.9_dp,-.6_dp]
    h=1e-6_dp
    call lu_z_update_jvp(x(1),x(2),x(3),x(4),x(5),x(6), &
        dx(1),dx(2),dx(3),dx(4),dx(5),dx(6),v(1),v(2),dv(1),dv(2))
    a=cmplx(x(1),x(2),dp);b=cmplx(x(3),x(4),dp);c=cmplx(x(5),x(6),dp)
    expected=a-b*c
    if(abs(cmplx(v(1),v(2),dp)-expected)>1e-14_dp) error stop 'LU complex update value'
    call lu_z_update_vjp(x(1),x(2),x(3),x(4),x(5),x(6),seed(1),seed(2), &
        v(1),v(2),bar(1),bar(2),bar(3),bar(4),bar(5),bar(6))
    if(abs(dot_product(seed,dv)-dot_product(bar,dx))>1e-13_dp) error stop 'LU update pairing'
    plus=update_probe(x+h*dx);minus=update_probe(x-h*dx)
    if(maxval(abs((plus-minus)/(2*h)-dv))>2e-9_dp) error stop 'LU update finite difference'
    call lu_z_quotient_jvp(x(1),x(2),x(3),x(4),dx(1),dx(2),dx(3),dx(4), &
        v(1),v(2),dv(1),dv(2))
    if(abs(cmplx(v(1),v(2),dp)-a/b)>1e-14_dp) error stop 'LU quotient value'
    call lu_z_quotient_vjp(x(1),x(2),x(3),x(4),seed(1),seed(2), &
        v(1),v(2),bar(1),bar(2),bar(3),bar(4))
    if(abs(dot_product(seed,dv)-dot_product(bar(:4),dx(:4)))>1e-13_dp) error stop 'LU quotient pairing'
    plus=quotient_probe(x+h*dx);minus=quotient_probe(x-h*dx)
    if(maxval(abs((plus-minus)/(2*h)-dv))>2e-9_dp) error stop 'LU quotient finite difference'
    do i=-1,1,2
        a=cmplx(10._dp**(300*i),10._dp**(300*i),dp)
        b=conjg(a)
        actual=lu_z_divide(a,b)
        if(abs(actual-cmplx(0,1,dp))>2e-15_dp) error stop 'LU scaled complex division'
    end do
    call indexed_updates()
    print '(a)', 'Generated complex LU arithmetic and real-pairing JVP/VJP pass'
contains
    subroutine indexed_updates()
        integer, parameter :: rows(5)=[8,2,6,3,5]
        complex(dp) :: base(10),direction(10),weights(5),dweights(5),shared,dshared
        complex(dp) :: result(10),expected_vector(10),xp(10),xm(10),tangent(10)
        complex(dp) :: seedz,value,derivative,wp,wm,coefficient,dcoefficient
        real(dp) :: primal(2),dprimal(2),bars(6),lhs,rhs_pair
        integer :: j,k,adj
        logical :: conjugate
        do j=1,10
            base(j)=cmplx(real(j-4,dp)/8,real(3-j,dp)/16,dp)
            direction(j)=cmplx(real(mod(j,3)-1,dp)/8,real(j,dp)/32,dp)
        end do
        do j=1,5
            weights(j)=cmplx(real(j-3,dp)/4,real(2-j,dp)/8,dp)
            dweights(j)=cmplx(real(3-j,dp)/16,real(j,dp)/32,dp)
        end do
        shared=cmplx(-.75_dp,.375_dp,dp)
        dshared=cmplx(.125_dp,-.25_dp,dp)
        seedz=cmplx(.875_dp,-.625_dp,dp)
        result=base
        call lu_z_scatter_update(result,rows,weights,5,shared)
        expected_vector=base
        do j=1,5
            expected_vector(rows(j))=base(rows(j))-weights(j)*shared
        end do
        if(maxval(abs(result-expected_vector))>1e-14_dp) &
            error stop 'complex indexed scatter independent value/padding'
        call lu_z_scatter_update(result,rows,weights,0,shared)
        if(any(result/=expected_vector)) error stop 'zero-length complex scatter'
        xp=base+h*direction;xm=base-h*direction
        call lu_z_scatter_update(xp,rows,weights+h*dweights,5,shared+h*dshared)
        call lu_z_scatter_update(xm,rows,weights-h*dweights,5,shared-h*dshared)
        tangent=direction
        rhs_pair=0.0_dp
        do j=1,5
            k=rows(j)
            call lu_z_update_jvp(real(base(k),dp),aimag(base(k)), &
                real(weights(j),dp),aimag(weights(j)),real(shared,dp),aimag(shared), &
                real(direction(k),dp),aimag(direction(k)),real(dweights(j),dp), &
                aimag(dweights(j)),real(dshared,dp),aimag(dshared), &
                primal(1),primal(2),dprimal(1),dprimal(2))
            tangent(k)=cmplx(dprimal(1),dprimal(2),dp)
            call lu_z_update_vjp(real(base(k),dp),aimag(base(k)), &
                real(weights(j),dp),aimag(weights(j)),real(shared,dp),aimag(shared), &
                real(seedz,dp),aimag(seedz),primal(1),primal(2), &
                bars(1),bars(2),bars(3),bars(4),bars(5),bars(6))
            rhs_pair=rhs_pair+dot_product(bars,[real(direction(k),dp), &
                aimag(direction(k)),real(dweights(j),dp),aimag(dweights(j)), &
                real(dshared,dp),aimag(dshared)])
        end do
        if(maxval(abs((xp-xm)/(2*h)-tangent))>2e-9_dp) &
            error stop 'complex indexed scatter versus symbolic JVP'
        lhs=0.0_dp
        do j=1,5
            lhs=lhs+real(conjg(seedz)*tangent(rows(j)),dp)
        end do
        if(abs(lhs-rhs_pair)>1e-13_dp) error stop 'complex scatter real-pairing VJP'
        do adj=0,1
            conjugate=adj==1
            value=shared
            call lu_z_gather_update(value,rows,weights,base,5,conjugate)
            coefficient=shared
            do j=1,5
                if(conjugate) then
                    coefficient=coefficient-conjg(weights(j))*base(rows(j))
                else
                    coefficient=coefficient-weights(j)*base(rows(j))
                end if
            end do
            if(abs(value-coefficient)>1e-14_dp) error stop 'complex gather dot oracle'
            call lu_z_gather_update(value,rows,weights,base,0,conjugate)
            if(value/=coefficient) error stop 'zero-length complex gather'
            wp=shared+h*dshared;wm=shared-h*dshared
            call lu_z_gather_update(wp,rows,weights+h*dweights, &
                base+h*direction,5,conjugate)
            call lu_z_gather_update(wm,rows,weights-h*dweights, &
                base-h*direction,5,conjugate)
            value=shared;derivative=dshared
            rhs_pair=real(conjg(seedz)*dshared,dp)
            do j=1,5
                k=rows(j)
                coefficient=weights(j);dcoefficient=dweights(j)
                if(conjugate) then
                    coefficient=conjg(coefficient);dcoefficient=conjg(dcoefficient)
                end if
                call lu_z_update_jvp(real(value,dp),aimag(value), &
                    real(coefficient,dp),aimag(coefficient),real(base(k),dp), &
                    aimag(base(k)),real(derivative,dp),aimag(derivative), &
                    real(dcoefficient,dp),aimag(dcoefficient), &
                    real(direction(k),dp),aimag(direction(k)), &
                    primal(1),primal(2),dprimal(1),dprimal(2))
                call lu_z_update_vjp(real(value,dp),aimag(value), &
                    real(coefficient,dp),aimag(coefficient),real(base(k),dp), &
                    aimag(base(k)),real(seedz,dp),aimag(seedz),primal(1),primal(2), &
                    bars(1),bars(2),bars(3),bars(4),bars(5),bars(6))
                value=cmplx(primal(1),primal(2),dp)
                derivative=cmplx(dprimal(1),dprimal(2),dp)
                rhs_pair=rhs_pair+dot_product(bars(3:6),[real(dcoefficient,dp), &
                    aimag(dcoefficient),real(direction(k),dp),aimag(direction(k))])
            end do
            if(abs((wp-wm)/(2*h)-derivative)>2e-9_dp) &
                error stop 'complex gather versus symbolic JVP'
            if(abs(real(conjg(seedz)*derivative,dp)-rhs_pair)>1e-13_dp) &
                error stop 'complex gather real-pairing VJP'
        end do
    end subroutine indexed_updates
    function update_probe(y) result(value)
        real(dp), intent(in) :: y(6)
        real(dp) :: value(2)
        complex(dp) :: z
        z=cmplx(y(1),y(2),dp)-cmplx(y(3),y(4),dp)*cmplx(y(5),y(6),dp)
        value=[real(z,dp),aimag(z)]
    end function update_probe
    function quotient_probe(y) result(value)
        real(dp), intent(in) :: y(6)
        real(dp) :: value(2)
        complex(dp) :: z
        z=cmplx(y(1),y(2),dp)/cmplx(y(3),y(4),dp)
        value=[real(z,dp),aimag(z)]
    end function quotient_probe
end program test_fortsparse_lu_leaves
