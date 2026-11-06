# Prefer optimized system kernels without overriding a consumer's BLAS policy.
# Load before backend packages put their own FindLAPACK on CMAKE_MODULE_PATH.
if(TARGET BLAS::BLAS)
    set(BLAS_FOUND TRUE)
    if(NOT BLAS_LIBRARIES)
        set(BLAS_LIBRARIES BLAS::BLAS)
    endif()
endif()
if(TARGET LAPACK::LAPACK)
    set(LAPACK_FOUND TRUE)
    if(NOT LAPACK_LIBRARIES)
        set(LAPACK_LIBRARIES LAPACK::LAPACK)
    endif()
endif()
if(NOT BLAS_FOUND AND NOT TARGET BLAS::BLAS)
    if(NOT DEFINED BLA_VENDOR AND NOT DEFINED ENV{BLA_VENDOR})
        foreach(_fortsparse_vendor OpenBLAS Intel10_64lp_seq)
            set(BLA_VENDOR "${_fortsparse_vendor}")
            find_package(BLAS QUIET MODULE)
            if(BLAS_FOUND)
                break()
            endif()
        endforeach()
        if(NOT BLAS_FOUND)
            set(BLA_VENDOR All)
        endif()
    endif()
    find_package(BLAS REQUIRED MODULE)
endif()
if(NOT LAPACK_FOUND AND NOT TARGET LAPACK::LAPACK)
    find_package(LAPACK REQUIRED MODULE)
endif()
set(_fortsparse_blas_libraries "${BLAS_LIBRARIES}")
set(_fortsparse_lapack_libraries "${LAPACK_LIBRARIES}")

# LAPACK < 3.7 (Apple Accelerate among them) has no DSYTRF_RK. The indefinite
# panel factor then uses DSYTRF and converts its storage. Recheck whenever the
# selected libraries change.
include(CheckFortranFunctionExists)
include(CMakePushCheckState)
if(NOT "${_fortsparse_lapack_libraries};${_fortsparse_blas_libraries}" STREQUAL
        "${_FORTSPARSE_DSYTRF_RK_CHECKED_WITH}")
    unset(FORTSPARSE_LAPACK_HAS_DSYTRF_RK CACHE)
    set(_FORTSPARSE_DSYTRF_RK_CHECKED_WITH
        "${_fortsparse_lapack_libraries};${_fortsparse_blas_libraries}"
        CACHE INTERNAL "")
endif()
cmake_push_check_state(RESET)
set(CMAKE_REQUIRED_LIBRARIES ${_fortsparse_lapack_libraries}
    ${_fortsparse_blas_libraries})
set(CMAKE_REQUIRED_QUIET ON)
check_fortran_function_exists(dsytrf_rk FORTSPARSE_LAPACK_HAS_DSYTRF_RK)
cmake_pop_check_state()
if(NOT FORTSPARSE_LAPACK_HAS_DSYTRF_RK)
    message(STATUS "fortsparse: LAPACK lacks DSYTRF_RK (LAPACK < 3.7); "
        "indefinite supernodal panels use DSYTRF (Bunch-Kaufman)")
endif()
