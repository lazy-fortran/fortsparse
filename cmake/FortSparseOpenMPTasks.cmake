# Probe the compiler/runtime pair instead of assuming OpenMP 5 iterator
# dependencies from a compiler version. Older tools retain the correct
# individual-update scheduler in fortsparse_native_supernodes.F90.
include(CheckFortranSourceCompiles)
function(fortsparse_check_task_iterator)
    set(CMAKE_REQUIRED_LIBRARIES OpenMP::OpenMP_Fortran)
    check_fortran_source_compiles([=[
program task_iterator
    implicit none
    integer :: token(4), sources(3), values(4), i, first, last
    token = 0
    sources = [1,2,3]
    values = 0
    first = 1
    last = 4
    !$omp parallel num_threads(2)
    !$omp single
    do i = 1, 3
        !$omp task firstprivate(i) depend(out:token(i))
        values(i) = i
        !$omp end task
    end do
    !$omp task depend(iterator(i=first:last), in:token(sources(i))) &
    !$omp depend(out:token(4))
    values(4) = sum(values(1:3))
    !$omp end task
    !$omp end single
    !$omp end parallel
end program task_iterator
]=] FORTSPARSE_HAVE_OMP_DEPEND_ITERATOR SRC_EXT F90)
    if(FORTSPARSE_HAVE_OMP_DEPEND_ITERATOR)
        target_compile_definitions(fortsparse PRIVATE FORTSPARSE_HAVE_OMP_DEPEND_ITERATOR)
    endif()
endfunction()
fortsparse_check_task_iterator()
