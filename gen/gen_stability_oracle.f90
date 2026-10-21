program gen_stability_oracle
    ! Unit lower L: det L=1, M=L L^T>0, A=L diag(d) L^T.
    ! det(A-sigma M)=product(d-sigma); the exact generalized eigenvalues
    ! are d. Check the polynomial identity before emitting the fixture.
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t,sym,num,operator(+),operator(-),operator(*)
    use fortsym_engine, only: engine_result_t,VERDICT_TRUE
    use fortsym_engine_native, only: native_engine_t,make_native_engine
    use fortsym_kernel, only: kernel_spec_t,emit_kernel,TARGET_FORTRAN_CPU
    use fortsym_string, only: str,chars
    implicit none
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(engine_result_t) :: check
    type(expr_t) :: l(3,3),a(3,3),m(3,3),p(3,3),d(3),sigma,roots(18)
    type(kernel_spec_t) :: spec
    integer :: i,j,k,unit
    character(3) :: name
    call arena%init()
    engine=make_native_engine(arena)
    do i=1,3
        do j=1,3
            l(i,j)=num(arena,0)
        end do
        l(i,i)=num(arena,1)
    end do
    l(2,1)=sym(arena,'l21')
    l(3,1)=sym(arena,'l31')
    l(3,2)=sym(arena,'l32')
    d=[sym(arena,'d1'),sym(arena,'d2'),sym(arena,'d3')]
    sigma=sym(arena,'sigma')
    do j=1,3
        do i=1,3
            a(i,j)=num(arena,0)
            m(i,j)=num(arena,0)
            do k=1,3
                a(i,j)=a(i,j)+l(i,k)*d(k)*l(j,k)
                m(i,j)=m(i,j)+l(i,k)*l(j,k)
            end do
            p(i,j)=a(i,j)-sigma*m(i,j)
        end do
    end do
    check=engine%zero_test(determinant(p)-(d(1)-sigma)*(d(2)-sigma)*(d(3)-sigma))
    if (.not.check%ok.or.check%verdict/=VERDICT_TRUE) &
        error stop 'FortSym did not prove the congruence pencil identity'
    roots(:9)=reshape(a,[9])
    roots(10:)=reshape(m,[9])
    spec%target=TARGET_FORTRAN_CPU
    spec%name=str('stability_pencil')
    spec%generator=str('gen/gen_stability_oracle.f90')
    spec%generator_revision=str('7e8ef7abdd90c49bb5f75e2c119e7d181656bcc6')
    spec%regenerate_command=str('gen/regenerate_stability.sh /path/to/fortsym')
    spec%pure_procedure=.true.
    spec%args=[str('l21'),str('l31'),str('l32'),str('d1'),str('d2'),str('d3')]
    allocate(spec%outputs(18))
    do k=1,18
        j=mod((k-1)/3,3)+1
        i=mod(k-1,3)+1
        write(name,'(a,i1,i1)') 'a',i,j
        if(k>9) name(1:1)='m'
        spec%outputs(k)=str(name)
    end do
    open(newunit=unit,file='test/solver/stability_oracle.f90',status='replace')
    write(unit,'(a)') '! Generated after FortSym proves det(A-sigma M)=product(d-sigma).'
    write(unit,'(a)') 'module stability_oracle'
    write(unit,'(a)') 'implicit none'
    write(unit,'(a)') 'contains'
    write(unit,'(a)') chars(emit_kernel(roots,spec))
    write(unit,'(a)') 'end module'
    close(unit)
    print '(a)', 'FortSym exact congruence pencil determinant identity: PROVED'
contains
    function determinant(p) result(det)
        type(expr_t), intent(in) :: p(3,3)
        type(expr_t) :: det
        det=p(1,1)*(p(2,2)*p(3,3)-p(2,3)*p(3,2))- &
            p(1,2)*(p(2,1)*p(3,3)-p(2,3)*p(3,1))+ &
            p(1,3)*(p(2,1)*p(3,2)-p(2,2)*p(3,1))
    end function
end program
