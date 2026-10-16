module fortsparse_ilutp_selection
    ! Integer sparse scheduling and bounded selection, independently written.
    use fortsparse_kinds, only: dp
    implicit none
    private
    public :: lower_push, lower_pop, select_largest
contains
    subroutine lower_push(heap, count, column)
        integer, intent(inout) :: heap(:), count
        integer, intent(in) :: column
        integer :: position, parent
        count = count + 1
        position = count
        do while (position > 1)
            parent = position/2
            if (heap(parent) <= column) exit
            heap(position) = heap(parent)
            position = parent
        end do
        heap(position) = column
    end subroutine

    integer function lower_pop(heap, count) result(column)
        integer, intent(inout) :: heap(:), count
        integer :: position, child, last
        column = heap(1)
        last = heap(count)
        count = count - 1
        if (count == 0) return
        position = 1
        do while (position <= count/2)
            child = 2*position
            if (child < count) then
                if (heap(child + 1) < heap(child)) child = child + 1
            end if
            if (last <= heap(child)) exit
            heap(position) = heap(child)
            position = child
        end do
        heap(position) = last
    end function

    subroutine select_largest(values, touched, count, row, lower, cap, &
                              selected, kept, eligible)
        real(dp), intent(in) :: values(:)
        integer, intent(in) :: touched(:), count, row, cap
        logical, intent(in) :: lower
        integer, intent(out) :: selected(:), kept, eligible
        integer :: k, column, position, parent, child, candidate
        ! A bounded heap keeps its worst selected entry at the root. Equal
        ! magnitudes prefer the smaller current elimination-column number.
        kept = 0
        eligible = 0
        do k = 1, count
            column = touched(k)
            if (lower) then
                if (column >= row) cycle
            else
                if (column <= row) cycle
            end if
            if (abs(values(column)) <= 0.0_dp) cycle
            eligible = eligible + 1
            if (cap == 0) cycle
            if (kept < cap) then
                kept = kept + 1
                position = kept
                do while (position > 1)
                    parent = position/2
                    if (.not. worse(column, selected(parent))) exit
                    selected(position) = selected(parent)
                    position = parent
                end do
                selected(position) = column
            else
                if (.not. worse(selected(1), column)) cycle
                position = 1
                do while (position <= kept/2)
                    child = 2*position
                    if (child < kept) then
                        candidate = selected(child + 1)
                        if (worse(candidate, &
                                  selected(child))) child = child + 1
                    end if
                    if (.not. worse(selected(child), column)) exit
                    selected(position) = selected(child)
                    position = child
                end do
                selected(position) = column
            end if
        end do
    contains
        logical function worse(left, right)
            integer, intent(in) :: left, right
            if (abs(values(left)) < abs(values(right))) then
                worse = .true.
            else if (abs(values(left)) > abs(values(right))) then
                worse = .false.
            else
                worse = left > right
            end if
        end function
    end subroutine
end module fortsparse_ilutp_selection
