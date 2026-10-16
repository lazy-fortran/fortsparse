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
    print '(a)', 'Generated complex LU arithmetic and real-pairing JVP/VJP pass'
contains
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
