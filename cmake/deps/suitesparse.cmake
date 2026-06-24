# Fetch and build SuiteSparse/UMFPACK (GPL) from source via ExternalProject,
# linked against the SYSTEM BLAS/LAPACK so the niche solver is the only thing
# built from source. Builds only the packages UMFPACK requires without CHOLMOD:
# SuiteSparse_config + AMD + CAMD + CCOLAMD + COLAMD + UMFPACK.
#
# Exposes IMPORTED targets SuiteSparse::umfpack, SuiteSparse::amd and
# SuiteSparse::suitesparse_config. Only the GPL helper executable links these;
# libfortsparse never does.
#
# Pattern adapted from sparse_draft/cmake/FetchSuiteSparse.cmake: ExternalProject
# (not FetchContent) so we can restrict the project set and pass the system
# BLAS/LAPACK out of tree. The top-level CMake must have run find_package(BLAS)
# and find_package(LAPACK) before including this file.

include(ExternalProject)

# Resolve the find_package BLAS/LAPACK results to concrete library file paths.
# ExternalProject configures out of tree and cannot see CMake imported targets,
# so pass it real file paths. find_package may yield either a list of file paths
# in BLAS_LIBRARIES / LAPACK_LIBRARIES or only the BLAS::BLAS / LAPACK::LAPACK
# imported targets (newer FindBLAS/FindLAPACK); in the latter case extract their
# IMPORTED_LOCATION. Concrete paths also survive the nested-FetchContent scope
# re-evaluation that a $<TARGET_FILE:...> genexp would abort on.
function(_fortsparse_resolve_libs out_var raw_libs imported_target)
    set(_resolved "")
    if(raw_libs)
        foreach(_lib IN LISTS raw_libs)
            if(TARGET ${_lib})
                get_target_property(_loc ${_lib} IMPORTED_LOCATION)
                if(_loc)
                    list(APPEND _resolved "${_loc}")
                endif()
            else()
                list(APPEND _resolved "${_lib}")
            endif()
        endforeach()
    endif()
    if(NOT _resolved AND TARGET ${imported_target})
        get_target_property(_loc ${imported_target} IMPORTED_LOCATION)
        if(_loc)
            list(APPEND _resolved "${_loc}")
        endif()
        get_target_property(_iface ${imported_target} INTERFACE_LINK_LIBRARIES)
        if(_iface)
            foreach(_lib IN LISTS _iface)
                if(TARGET ${_lib})
                    get_target_property(_loc ${_lib} IMPORTED_LOCATION)
                    if(_loc)
                        list(APPEND _resolved "${_loc}")
                    endif()
                else()
                    list(APPEND _resolved "${_lib}")
                endif()
            endforeach()
        endif()
    endif()
    set(${out_var} "${_resolved}" PARENT_SCOPE)
endfunction()

_fortsparse_resolve_libs(FORTSPARSE_BLAS_LIBS "${BLAS_LIBRARIES}" BLAS::BLAS)
_fortsparse_resolve_libs(FORTSPARSE_LAPACK_LIBS "${LAPACK_LIBRARIES}" LAPACK::LAPACK)

# Join multi-entry lists with the pipe LIST_SEPARATOR so they reach the sub-build
# as ";"-separated lists without tripping ExternalProject's genexp evaluation.
string(REPLACE ";" "|" FORTSPARSE_BLAS_LIBS_ARG "${FORTSPARSE_BLAS_LIBS}")
string(REPLACE ";" "|" FORTSPARSE_LAPACK_LIBS_ARG "${FORTSPARSE_LAPACK_LIBS}")

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
        -DBLAS_LIBRARIES=${FORTSPARSE_BLAS_LIBS_ARG}
        -DLAPACK_LIBRARIES=${FORTSPARSE_LAPACK_LIBS_ARG}
        -DCMAKE_INSTALL_PREFIX=${SUITESPARSE_INSTALL_PREFIX}
        -DCMAKE_INSTALL_LIBDIR=lib
    BUILD_BYPRODUCTS
        ${UMFPACK_LIBRARY_PATH}
        ${AMD_LIBRARY_PATH}
        ${CAMD_LIBRARY_PATH}
        ${CCOLAMD_LIBRARY_PATH}
        ${COLAMD_LIBRARY_PATH}
        ${SUITESPARSE_CONFIG_LIBRARY_PATH})

# The system BLAS/LAPACK targets carry the link interface the GPL helper needs;
# attach them to the lowest IMPORTED target so the helper resolves dense kernels
# against the system OpenBLAS/LAPACK.
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
