program native_only_smoke
    use fortsparse
    use fortsparse_native, only: native_backend_t
    use fortsparse_native_general, only: native_general_backend_t
    implicit none
    type(sparse_solver_t) :: solver
    type(sparse_factor_stats_t) :: stats
    type(fortsparse_status_t) :: status
    type(csc_t) :: A
    type(csc_z_t) :: Z
    real(dp) :: dense(3,3), exact(3), rhs(3), x(3)
    complex(dp) :: zdense(3,3), zexact(3), zrhs(3), zx(3)
    integer, parameter :: ptr(4) = [1,4,7,10]
    integer, parameter :: rows(9) = [1,2,3,1,2,3,1,2,3]

    call require(sparse_auto_backend(9,.false.) == FORTSPARSE_BACKEND_NATIVE, &
        'native-only AUTO real capability')
    call require(sparse_auto_backend(huge(1),.true.) == FORTSPARSE_BACKEND_NATIVE, &
        'native-only AUTO complex capability at large size')
    dense = reshape([0.0_dp,1.0_dp,4.0_dp, 2.0_dp,0.0_dp,1.0_dp, &
        1.0_dp,3.0_dp,2.0_dp], [3,3])
    A%nrow = 3
    A%ncol = 3
    A%nnz = 9
    A%col_ptr = ptr
    A%row_idx = rows
    A%val = reshape(dense,[9])
    exact = [1.0_dp,-2.0_dp,.5_dp]
    rhs = matmul(dense,exact)
    solver%options%host_memory_limit_bytes = 1_i8
    call sparse_factor(solver,A,status)
    call require(status%code == FORTSPARSE_RESOURCE_LIMIT, 'AUTO real cap refuses growth')
    call require(.not. solver%factored, 'cap failure leaves unfactored state')
    solver%options%host_memory_limit_bytes = 1024_i8*1024_i8
    call sparse_factor(solver,A,status)
    call checked(status)
    select type(backend => solver%backend)
    type is(native_general_backend_t)
    class default
        error stop 'AUTO did not construct native general backend'
    end select
    call sparse_solve(solver,rhs,x,status)
    call checked(status)
    call require(maxval(abs(x-exact)) < 2.e-13_dp, 'independent pivoted real solve')
    call sparse_factor_stats(solver,stats,status)
    call checked(status)
    call require(stats%strict_memory_supported, 'native real cap contract')
    call require(stats%host_peak_bytes <= solver%options%host_memory_limit_bytes, &
        'real payload remains below cap')
    call sparse_destroy(solver)

    Z%nrow = 3
    Z%ncol = 3
    Z%nnz = 9
    Z%col_ptr = ptr
    Z%row_idx = rows
    zdense = cmplx(dense,reshape([.1_dp,-.2_dp,0.0_dp,.3_dp,.1_dp, &
        -.1_dp,.4_dp,-.2_dp,.2_dp],[3,3]),dp)
    Z%val = reshape(zdense,[9])
    zexact = cmplx(exact,[.2_dp,.1_dp,-.3_dp],dp)
    zrhs = matmul(zdense,zexact)
    solver%options = sparse_options_t()
    solver%options%host_memory_limit_bytes = 1024_i8*1024_i8
    call sparse_factor(solver,Z,status)
    call checked(status)
    select type(backend => solver%backend)
    type is(native_general_backend_t)
    class default
        error stop 'AUTO did not construct native complex backend'
    end select
    call sparse_solve(solver,zrhs,zx,status)
    call checked(status)
    call require(maxval(abs(zx-zexact)) < 2.e-13_dp, 'independent complex solve')
    zrhs = matmul(conjg(transpose(zdense)),zexact)
    call sparse_solve_op(solver,zrhs,zx,FORTSPARSE_OP_ADJOINT,status)
    call checked(status)
    call require(maxval(abs(zx-zexact)) < 2.e-13_dp, 'native complex adjoint solve')
    call sparse_destroy(solver)

    dense = reshape([4.0_dp,-1.0_dp,0.0_dp, -1.0_dp,4.0_dp,-1.0_dp, &
        0.0_dp,-1.0_dp,4.0_dp], [3,3])
    A%val = reshape(dense,[9])
    rhs = matmul(dense,exact)
    solver%options = sparse_options_t()
    solver%options%symmetry = FORTSPARSE_SYM_SPD
    solver%options%host_memory_limit_bytes = 4096_i8
    call sparse_factor(solver,A,status)
    call checked(status)
    select type(backend => solver%backend)
    type is(native_backend_t)
    class default
        error stop 'AUTO SPD did not construct native SPD backend'
    end select
    call sparse_solve(solver,rhs,x,status)
    call checked(status)
    call require(maxval(abs(x-exact)) < 2.e-13_dp, 'independent capped SPD solve')
    call sparse_destroy(solver)

    solver%options = sparse_options_t()
    solver%backend_id = FORTSPARSE_BACKEND_SUPERLU
    call sparse_factor(solver,A,status)
    call require(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
        'explicit unavailable SuperLU status')
    solver%backend_id = FORTSPARSE_BACKEND_SUPERLU_ILU
    call sparse_factor(solver,A,status)
    call require(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
        'explicit unavailable SuperLU ILU status')
    solver%backend_id = FORTSPARSE_BACKEND_UMFPACK_IPC
    call sparse_factor(solver,A,status)
    call require(status%code == FORTSPARSE_BACKEND_UNAVAILABLE, &
        'IPC without helper fails explicitly')
    call sparse_destroy(solver)
    print '(a)', 'Native-only AUTO real/complex/SPD, caps and unavailable providers pass'
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print '(a)', message
            error stop 1
        end if
    end subroutine require
    subroutine checked(s)
        type(fortsparse_status_t), intent(in) :: s
        if (.not. status_ok(s)) then
            print '(i0,2a)', s%code, ': ', s%msg
            error stop 1
        end if
    end subroutine checked
end program native_only_smoke
