! Pivoted-LU reference for the static-structure diagnostic.
! Usage: lu_probe matrix.mtx backend matching ordering threshold [pref.bin]
!   backend 1: native scalar LU (lu_factor_t via observer copy exposing p/q)
!   backend 2: SuperLU through the public API (default SuperLU options)
!   matching 0 none, 1 native structural, 2 weighted preference file
!   ordering 1 AMD, 2 METIS ND (production *_general wrappers on P*A)
! Prints factor nnz, pivot departures, original-COO real128 backward errors.
program lu_probe
    use, intrinsic :: iso_c_binding, only: c_int
    use, intrinsic :: iso_fortran_env, only: real64, int64, int32
    use fortsparse_csc, only: csc_t, csc_from_triplet
    use fortsparse_status, only: fortsparse_status_t
    use fortsparse_native_matching, only: structural_row_matching, matching_graph_rows
    use fortsparse_native_lu, only: lu_factor_t
    use fortsparse, only: sparse_solver_t, sparse_factor_stats_t, sparse_factor, &
        sparse_factor_stats, sparse_solve, FORTSPARSE_BACKEND_SUPERLU
    implicit none
    integer, parameter :: dp = real64, i8 = int64, qp = selected_real_kind(33)
    interface
        integer(c_int) function order_amd_general(n,nnz,ptr,row,perm,inverse) &
            bind(C,name='fsparse_native_order_amd_general')
            import :: c_int
            integer(c_int), value :: n,nnz
            integer(c_int), intent(in) :: ptr(*),row(*)
            integer(c_int), intent(out) :: perm(*),inverse(*)
        end function
        integer(c_int) function order_nd_general(n,nnz,ptr,row,perm,inverse) &
            bind(C,name='fsparse_native_order_general')
            import :: c_int
            integer(c_int), value :: n,nnz
            integer(c_int), intent(in) :: ptr(*),row(*)
            integer(c_int), intent(out) :: perm(*),inverse(*)
        end function
    end interface
    type(csc_t) :: a
    type(fortsparse_status_t) :: status
    type(lu_factor_t) :: lu
    type(sparse_solver_t) :: slu
    type(sparse_factor_stats_t) :: stats
    integer :: n, m, nz, e, u, backend, matching, ordering, i, j, code, departures, row
    integer, allocatable :: rowc(:), colc(:), pref(:), inv(:), work(:), graph(:), perm(:), iperm(:)
    real(dp), allocatable :: val(:), known(:), b(:), x(:)
    real(qp), allocatable :: rq(:), denom(:), anorm(:)
    real(dp) :: threshold, t
    integer(i8) :: t0, t1, rate, fnnz
    character(512) :: path, arg, wpath

    call get_command_argument(1, path)
    call get_command_argument(2, arg); read(arg,*) backend
    call get_command_argument(3, arg); read(arg,*) matching
    call get_command_argument(4, arg); read(arg,*) ordering
    call get_command_argument(5, arg); read(arg,*) threshold
    if (matching == 2) call get_command_argument(6, wpath)
    open(newunit=u, file=trim(path), status='old')
    block
        character(512) :: line
        do
            read(u,'(a)') line
            if (line(1:1) /= '%') exit
        end do
        read(line,*) n, m, nz
    end block
    allocate(rowc(nz), colc(nz), val(nz), known(n), b(n), x(n), rq(n), denom(n), anorm(n))
    do e = 1, nz
        read(u,*) rowc(e), colc(e), val(e)
    end do
    close(u)
    call csc_from_triplet(n, n, rowc, colc, val, a, status)
    if (status%code /= 0) error stop 'invalid matrix'
    do i = 1, n
        known(i) = sin(real(i,dp)) + 0.125_dp*cos(real(3*i,dp))
    end do
    rq = 0; anorm = 0
    do e = 1, nz
        rq(rowc(e)) = rq(rowc(e)) + real(val(e),qp)*real(known(colc(e)),qp)
        anorm(rowc(e)) = anorm(rowc(e)) + abs(real(val(e),qp))
    end do
    b = real(rq, dp)
    x = b
    departures = -1
    call system_clock(t0, rate)
    if (backend == 1) then
        allocate(pref(n), inv(n), work(5*n), graph(size(a%row_idx)), perm(n), iperm(n))
        select case (matching)
        case (0)
            graph = a%row_idx
        case (1)
            call structural_row_matching(n, a%col_ptr, a%row_idx, a%val, pref, inv, work, status)
            if (status%code /= 0) error stop 'matching failed'
            call matching_graph_rows(a%col_ptr, a%row_idx, inv, graph)
        case (2)
            open(newunit=u, file=trim(wpath), access='stream', status='old')
            read(u) pref
            close(u)
            do i = 1, n
                inv(pref(i)) = i
            end do
            call matching_graph_rows(a%col_ptr, a%row_idx, inv, graph)
        end select
        if (ordering == 1) then
            code = order_amd_general(n, size(graph), a%col_ptr, graph, perm, iperm)
        else
            code = order_nd_general(n, size(graph), a%col_ptr, graph, perm, iperm)
        end if
        if (code /= 0) error stop 'ordering failed'
        write(*,'(a)') 'READY'; flush(6)
        call system_clock(t0)
        if (matching == 0) then
            call lu%factor(a, status, column_permutation=perm, pivot_threshold=threshold, pruning=.true.)
        else
            call lu%factor(a, status, column_permutation=perm, pivot_threshold=threshold, &
                pruning=.true., row_preference=pref)
        end if
        call system_clock(t1)
        if (status%code /= 0) then
            write(*,*) 'FAIL', status%code, trim(status%msg)
            stop
        end if
        fnnz = lu%factor_nnz
        departures = 0
        do j = 1, n
            if (matching == 0) then
                row = lu%q(j)
            else
                row = pref(lu%q(j))
            end if
            if (lu%p(j) /= row) departures = departures + 1
        end do
        call lu%solve(x, status)
        if (status%code /= 0) error stop 'solve failed'
        write(*,'(a,i0,a,i0,a,i0)') 'NATIVE_LU ln=', lu%ln, ' un=', lu%un, ' peak_bytes=', lu%peak_memory_bytes
    else
        slu%backend_id = FORTSPARSE_BACKEND_SUPERLU
        write(*,'(a)') 'READY'; flush(6)
        call system_clock(t0)
        call sparse_factor(slu, a, status)
        call system_clock(t1)
        if (status%code /= 0) then
            write(*,*) 'FAIL', status%code, trim(status%msg)
            stop
        end if
        call sparse_factor_stats(slu, stats, status)
        fnnz = stats%factor_nnz
        call sparse_solve(slu, b, x, status)
        if (status%code /= 0) error stop 'solve failed'
    end if
    t = real(t1-t0,dp)/rate
    rq = -real(b,qp); denom = abs(real(b,qp))
    do e = 1, nz
        rq(rowc(e)) = rq(rowc(e)) + real(val(e),qp)*real(x(colc(e)),qp)
        denom(rowc(e)) = denom(rowc(e)) + abs(real(val(e),qp)*real(x(colc(e)),qp))
    end do
    write(*,'(a,*(g0,:,","))') 'RESULT,', backend, matching, ordering, threshold, fnnz, departures, &
        real(maxval(abs(rq))/(maxval(anorm)*maxval(abs(real(x,qp)))+maxval(abs(real(b,qp)))),dp), &
        real(maxval(abs(rq)/max(denom,tiny(1._qp))),dp), maxval(abs(x-known))/maxval(abs(known)), t
end program lu_probe
