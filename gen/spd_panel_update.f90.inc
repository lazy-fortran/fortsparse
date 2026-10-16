    ! Full L spans have distinct original row IDs. Active lane bits are graph
    ! reach, independent of numerical cancellation. No allocation or dispatch.
    pure subroutine spd_panel_scatter(x, masks, indices, values, n, coeff, lanes, active)
        real(dp), intent(inout) :: x(8,*)
        integer, intent(inout) :: masks(*)
        integer, intent(in) :: indices(*), n, lanes, active
        real(dp), intent(in) :: values(*), coeff(*)
        real(dp) :: a, b, c, y
        integer :: i, lane, row, valid, first_lane, last_lane, range_mask
        ! Uniform mask permits contiguous lane SIMD; arbitrary masks retain
        ! the original guarded schedule. Caller contract: 1 <= lanes <= 8.
        if (active == 2**lanes-1) then
            do i = 1, n
                row = indices(i)
                masks(row) = ior(masks(row), active)
                b = values(i)
                !$omp simd private(a,c,y)
                do lane = 1, lanes
                    a = x(lane,row)
                    c = coeff(lane)
@PANEL_VALUE@
                    x(lane,row) = y
                end do
            end do
        else
            valid = ibits(active,0,lanes)
            first_lane = 1
            last_lane = 0
            range_mask = 0
            if (valid /= 0) then
                first_lane = trailz(valid)+1
                last_lane = bit_size(valid)-leadz(valid)
                range_mask = ishft(2**(last_lane-first_lane+1)-1,first_lane-1)
            end if
            ! An empty run evaluates no numerical entries. Gapped masks retain
            ! the original guarded path; no inactive lane is evaluated.
            if (valid == range_mask) then
                do i = 1, n
                    row = indices(i)
                    masks(row) = ior(masks(row), active)
                    b = values(i)
                    !$omp simd private(a,c,y)
                    do lane = first_lane, last_lane
                        a = x(lane,row)
                        c = coeff(lane)
@PANEL_VALUE@
                        x(lane,row) = y
                    end do
                end do
            else
                do i = 1, n
                    row = indices(i)
                    masks(row) = ior(masks(row), active)
                    b = values(i)
                    !$omp simd private(a,c,y)
                    do lane = 1, lanes
                        if (.not. btest(active,lane-1)) cycle
                        a = x(lane,row)
                        c = coeff(lane)
@PANEL_VALUE@
                        x(lane,row) = y
                    end do
                end do
            end if
        end if
    end subroutine spd_panel_scatter

    pure subroutine spd_panel_scatter_jvp(x, dx, masks, indices, values, &
            dvalues, n, coeff, dcoeff, lanes, active)
        real(dp), intent(inout) :: x(8,*), dx(8,*)
        integer, intent(inout) :: masks(*)
        integer, intent(in) :: indices(*), n, lanes, active
        real(dp), intent(in) :: values(*), dvalues(*), coeff(*), dcoeff(*)
        real(dp) :: a, b, c, da, db, dc, y, dy
        integer :: i, lane, row, valid, first_lane, last_lane, range_mask
        ! Uniform mask permits contiguous lane SIMD; arbitrary masks retain
        ! the original guarded schedule. Caller contract: 1 <= lanes <= 8.
        if (active == 2**lanes-1) then
            do i = 1, n
                row = indices(i)
                masks(row) = ior(masks(row), active)
                b = values(i)
                db = dvalues(i)
                !$omp simd private(a,c,da,dc,y,dy)
                do lane = 1, lanes
                    a = x(lane,row)
                    c = coeff(lane)
                    da = dx(lane,row)
                    dc = dcoeff(lane)
@PANEL_JVP@
                    x(lane,row) = y
                    dx(lane,row) = dy
                end do
            end do
        else
            valid = ibits(active,0,lanes)
            first_lane = 1
            last_lane = 0
            range_mask = 0
            if (valid /= 0) then
                first_lane = trailz(valid)+1
                last_lane = bit_size(valid)-leadz(valid)
                range_mask = ishft(2**(last_lane-first_lane+1)-1,first_lane-1)
            end if
            ! An empty run evaluates no numerical entries. Gapped masks retain
            ! the original guarded path; no inactive lane is evaluated.
            if (valid == range_mask) then
                do i = 1, n
                    row = indices(i)
                    masks(row) = ior(masks(row), active)
                    b = values(i)
                    db = dvalues(i)
                    !$omp simd private(a,c,da,dc,y,dy)
                    do lane = first_lane, last_lane
                        a = x(lane,row)
                        c = coeff(lane)
                        da = dx(lane,row)
                        dc = dcoeff(lane)
@PANEL_JVP@
                        x(lane,row) = y
                        dx(lane,row) = dy
                    end do
                end do
            else
                do i = 1, n
                    row = indices(i)
                    masks(row) = ior(masks(row), active)
                    b = values(i)
                    db = dvalues(i)
                    !$omp simd private(a,c,da,dc,y,dy)
                    do lane = 1, lanes
                        if (.not. btest(active,lane-1)) cycle
                        a = x(lane,row)
                        c = coeff(lane)
                        da = dx(lane,row)
                        dc = dcoeff(lane)
@PANEL_JVP@
                        x(lane,row) = y
                        dx(lane,row) = dy
                    end do
                end do
            end if
        end if
    end subroutine spd_panel_scatter_jvp

    pure subroutine spd_panel_scatter_vjp(x, xbar, indices, values, n, &
            coeff, lanes, active, values_bar, coeff_bar)
        real(dp), intent(in) :: x(8,*), values(*), coeff(*)
        real(dp), intent(inout) :: xbar(8,*), values_bar(*), coeff_bar(*)
        integer, intent(in) :: indices(*), n, lanes, active
        real(dp) :: a, b, c, seed, bacc, cacc, y, abar, bnext, cnext
        integer :: i, lane, row
        do i = 1, n
            row = indices(i)
            b = values(i)
            do lane = 1, lanes
                if (.not. btest(active,lane-1)) cycle
                a = x(lane,row)
                c = coeff(lane)
                seed = xbar(lane,row)
                bacc = values_bar(i)
                cacc = coeff_bar(lane)
@PANEL_VJP@
                xbar(lane,row) = abar
                values_bar(i) = bnext
                coeff_bar(lane) = cnext
            end do
        end do
    end subroutine spd_panel_scatter_vjp
