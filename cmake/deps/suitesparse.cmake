# Fetch and build SuiteSparse/UMFPACK (GPL) from source via ExternalProject,
# linked against fortnum's BLAS and LAPACK (fortnum_lapack77). Builds only the packages UMFPACK requires without CHOLMOD:
# SuiteSparse_config + AMD + CAMD + CCOLAMD + COLAMD + UMFPACK.
#
# Exposes IMPORTED targets SuiteSparse::umfpack, SuiteSparse::amd and
# SuiteSparse::suitesparse_config. Only the GPL helper executable links these;
# libfortsparse never does.
#
# Pattern adapted from sparse_draft/cmake/FetchSuiteSparse.cmake: ExternalProject
# (not FetchContent) so we can restrict the project set and pass the dense
# libraries out of tree. The top-level CMake must have added fortnum and
# included fortsparse_blas.cmake before this file.

include(ExternalProject)

# The out-of-tree SuiteSparse build sees files, not targets: fortnum's
# archives in link order and the Fortran runtime they need.
set(FORTSPARSE_DENSE_TARGETS fortnum_lapack77 fortnum fortgen_runtime)
set(FORTSPARSE_DENSE_LIBRARIES "")
foreach(_target IN LISTS FORTSPARSE_DENSE_TARGETS)
    get_target_property(_dir ${_target} ARCHIVE_OUTPUT_DIRECTORY)
    if(NOT _dir)
        get_target_property(_dir ${_target} BINARY_DIR)
    endif()
    list(APPEND FORTSPARSE_DENSE_LIBRARIES
        "${_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}${_target}${CMAKE_STATIC_LIBRARY_SUFFIX}")
endforeach()
if(OpenMP_Fortran_FOUND)
    list(APPEND FORTSPARSE_DENSE_LIBRARIES ${OpenMP_Fortran_LIBRARIES})
endif()
foreach(_dir IN LISTS CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES)
    list(APPEND FORTSPARSE_DENSE_LIBRARIES "-L${_dir}")
endforeach()
foreach(_lib IN LISTS CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)
    if(IS_ABSOLUTE "${_lib}" OR _lib MATCHES "^-")
        list(APPEND FORTSPARSE_DENSE_LIBRARIES "${_lib}")
    else()
        list(APPEND FORTSPARSE_DENSE_LIBRARIES "-l${_lib}")
    endif()
endforeach()

# Join the list with the pipe LIST_SEPARATOR so it reaches the sub-build as a
# ";"-separated list without tripping ExternalProject's genexp evaluation.
string(REPLACE ";" "|" FORTSPARSE_DENSE_LIBRARIES_ARG "${FORTSPARSE_DENSE_LIBRARIES}")

set(SUITESPARSE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/suitesparse-install)
set(SUITESPARSE_INCLUDE_DIR ${SUITESPARSE_INSTALL_PREFIX}/include)

set(UMFPACK_LIBRARY_PATH
    ${SUITESPARSE_INSTALL_PREFIX}/lib/libumfpack${CMAKE_STATIC_LIBRARY_SUFFIX})
set(AMD_LIBRARY_PATH
    ${SUITESPARSE_INSTALL_PREFIX}/lib/libamd${CMAKE_STATIC_LIBRARY_SUFFIX})
set(CAMD_LIBRARY_PATH
    ${SUITESPARSE_INSTALL_PREFIX}/lib/libcamd${CMAKE_STATIC_LIBRARY_SUFFIX})
set(CCOLAMD_LIBRARY_PATH
    ${SUITESPARSE_INSTALL_PREFIX}/lib/libccolamd${CMAKE_STATIC_LIBRARY_SUFFIX})
set(COLAMD_LIBRARY_PATH
    ${SUITESPARSE_INSTALL_PREFIX}/lib/libcolamd${CMAKE_STATIC_LIBRARY_SUFFIX})
set(SUITESPARSE_CONFIG_LIBRARY_PATH
    ${SUITESPARSE_INSTALL_PREFIX}/lib/libsuitesparseconfig${CMAKE_STATIC_LIBRARY_SUFFIX})

# The helper includes <suitesparse/umfpack.h>, so the include dir must be the
# parent of the suitesparse/ header directory the install creates.
file(MAKE_DIRECTORY ${SUITESPARSE_INCLUDE_DIR})

ExternalProject_Add(
    SuiteSparse
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    URL https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v7.8.3.tar.gz
    URL_HASH MD5=242e38ecfc8a3e3aa6b7d8d44849c5cf
    DEPENDS ${FORTSPARSE_DENSE_TARGETS}
    # Pass list-valued args pipe-separated and tell ExternalProject to translate
    # "|" back into ";" inside the sub-build. A literal ";" (or a $<SEMICOLON>
    # genexp) in CMAKE_ARGS makes ExternalProject emit an "add_custom_command
    # EVAL" error when fortsparse is configured as a nested FetchContent
    # sub-project, which silently drops the source build.
    LIST_SEPARATOR |
    CMAKE_ARGS
        -G Ninja
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
        -DCMAKE_BUILD_TYPE=Release
        -DBUILD_SHARED_LIBS=OFF
        -DSUITESPARSE_ENABLE_PROJECTS=suitesparse_config|amd|camd|ccolamd|colamd|umfpack
        -DUMFPACK_USE_CHOLMOD=OFF
        -DSUITESPARSE_USE_OPENMP=OFF
        -DSUITESPARSE_DEMOS=OFF
        -DBLAS_LIBRARIES=${FORTSPARSE_DENSE_LIBRARIES_ARG}
        -DLAPACK_LIBRARIES=${FORTSPARSE_DENSE_LIBRARIES_ARG}
        -DCMAKE_INSTALL_PREFIX=${SUITESPARSE_INSTALL_PREFIX}
        -DCMAKE_INSTALL_LIBDIR=lib
    BUILD_BYPRODUCTS
        ${UMFPACK_LIBRARY_PATH}
        ${AMD_LIBRARY_PATH}
        ${CAMD_LIBRARY_PATH}
        ${CCOLAMD_LIBRARY_PATH}
        ${COLAMD_LIBRARY_PATH}
        ${SUITESPARSE_CONFIG_LIBRARY_PATH})

# BLAS::BLAS and LAPACK::LAPACK carry the link interface the GPL helper
# needs; attach them to the lowest IMPORTED target so the helper resolves the
# dense kernels against fortnum.
set(_fortsparse_blas_iface "")
if(TARGET BLAS::BLAS)
    list(APPEND _fortsparse_blas_iface BLAS::BLAS)
elseif(BLAS_LIBRARIES)
    list(APPEND _fortsparse_blas_iface ${BLAS_LIBRARIES})
endif()
if(TARGET LAPACK::LAPACK)
    list(APPEND _fortsparse_blas_iface LAPACK::LAPACK)
elseif(LAPACK_LIBRARIES)
    list(APPEND _fortsparse_blas_iface ${LAPACK_LIBRARIES})
endif()

add_library(SuiteSparse::suitesparse_config STATIC IMPORTED)
set_target_properties(SuiteSparse::suitesparse_config PROPERTIES
    IMPORTED_LOCATION ${SUITESPARSE_CONFIG_LIBRARY_PATH}
    INTERFACE_INCLUDE_DIRECTORIES ${SUITESPARSE_INCLUDE_DIR}
    INTERFACE_LINK_LIBRARIES "${_fortsparse_blas_iface}")
add_dependencies(SuiteSparse::suitesparse_config SuiteSparse)

add_library(SuiteSparse::amd STATIC IMPORTED)
set_target_properties(SuiteSparse::amd PROPERTIES
    IMPORTED_LOCATION ${AMD_LIBRARY_PATH}
    INTERFACE_INCLUDE_DIRECTORIES ${SUITESPARSE_INCLUDE_DIR}
    INTERFACE_LINK_LIBRARIES SuiteSparse::suitesparse_config)
add_dependencies(SuiteSparse::amd SuiteSparse)

add_library(SuiteSparse::umfpack STATIC IMPORTED)
set_target_properties(SuiteSparse::umfpack PROPERTIES
    IMPORTED_LOCATION ${UMFPACK_LIBRARY_PATH}
    INTERFACE_INCLUDE_DIRECTORIES ${SUITESPARSE_INCLUDE_DIR}
    INTERFACE_LINK_LIBRARIES SuiteSparse::amd)
add_dependencies(SuiteSparse::umfpack SuiteSparse)
