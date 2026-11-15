# Dense kernels are fortnum's; fortsparse needs no system BLAS or LAPACK.
# fortsparse calls the fortnum_ routines directly, and SuperLU and the UMFPACK
# helper link fortnum's BLAS and LAPACK under their standard names
# (fortnum_lapack77). A consumer that defines BLAS::BLAS and LAPACK::LAPACK
# before adding fortsparse keeps them for these two libraries.
if(NOT TARGET BLAS::BLAS)
    add_library(BLAS::BLAS INTERFACE IMPORTED GLOBAL)
    set_target_properties(BLAS::BLAS PROPERTIES
        INTERFACE_LINK_LIBRARIES fortnum_lapack77)
endif()
if(NOT TARGET LAPACK::LAPACK)
    add_library(LAPACK::LAPACK INTERFACE IMPORTED GLOBAL)
    set_target_properties(LAPACK::LAPACK PROPERTIES
        INTERFACE_LINK_LIBRARIES fortnum_lapack77)
endif()

# Imported targets that link these (the SuiteSparse libraries) inherit the
# module directories of fortnum, which CMake requires to exist when generating.
foreach(_target fortnum_lapack77 fortnum fortgen_runtime)
    get_target_property(_dirs ${_target} INTERFACE_INCLUDE_DIRECTORIES)
    foreach(_dir IN LISTS _dirs)
        string(REGEX REPLACE "^[$]<BUILD_INTERFACE:(.*)>$" "\\1" _dir "${_dir}")
        if(_dir AND NOT _dir MATCHES "[$]<")
            file(MAKE_DIRECTORY "${_dir}")
        endif()
    endforeach()
endforeach()
