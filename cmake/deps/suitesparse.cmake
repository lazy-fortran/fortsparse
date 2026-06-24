# Fetch and build SuiteSparse/UMFPACK (GPL) from source via ExternalProject,
# pointed at the fetched reference BLAS/LAPACK so the whole chain needs no system
# numerical libraries. Builds only the packages UMFPACK requires without CHOLMOD:
# SuiteSparse_config + AMD + CAMD + CCOLAMD + COLAMD + UMFPACK.
#
# Exposes IMPORTED targets SuiteSparse::umfpack, SuiteSparse::amd and
# SuiteSparse::suitesparse_config. Only the GPL helper executable links these;
# libfortsparse never does.
#
# Pattern adapted from sparse_draft/cmake/FetchSuiteSparse.cmake: ExternalProject
# (not FetchContent) so we can restrict the project set and pass an out-of-tree
# BLAS. include() cmake/deps/reflapack.cmake before this file so the `blas` /
# `lapack` targets and their archive paths exist.

include(ExternalProject)

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
    # Pass the SUITESPARSE_ENABLE_PROJECTS list pipe-separated and tell
    # ExternalProject to translate "|" back into ";" inside the sub-build. A
    # literal ";" (or a $<SEMICOLON> genexp) in CMAKE_ARGS makes ExternalProject
    # emit an "add_custom_command EVAL" error when fortsparse is configured as a
    # nested FetchContent sub-project, which silently drops the source build.
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
        -DBLAS_LIBRARIES=${FORTSPARSE_REFLAPACK_BLAS_LIB}
        -DLAPACK_LIBRARIES=${FORTSPARSE_REFLAPACK_LAPACK_LIB}
        -DCMAKE_INSTALL_PREFIX=${SUITESPARSE_INSTALL_PREFIX}
        -DCMAKE_INSTALL_LIBDIR=lib
    BUILD_BYPRODUCTS
        ${UMFPACK_LIBRARY_PATH}
        ${AMD_LIBRARY_PATH}
        ${CAMD_LIBRARY_PATH}
        ${CCOLAMD_LIBRARY_PATH}
        ${COLAMD_LIBRARY_PATH}
        ${SUITESPARSE_CONFIG_LIBRARY_PATH})

# The reference BLAS/LAPACK archives must exist before SuiteSparse configures.
add_dependencies(SuiteSparse blas lapack)

add_library(SuiteSparse::suitesparse_config STATIC IMPORTED)
set_target_properties(SuiteSparse::suitesparse_config PROPERTIES
    IMPORTED_LOCATION ${SUITESPARSE_CONFIG_LIBRARY_PATH}
    INTERFACE_INCLUDE_DIRECTORIES ${SUITESPARSE_INCLUDE_DIR}
    INTERFACE_LINK_LIBRARIES "blas;lapack")
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
