# Fetch and build Reference-LAPACK from source (static), giving the targets
# `blas` and `lapack`. This is the zero-system-dependency BLAS/LAPACK path: it
# needs only a Fortran compiler. Correctness-first reference implementation, not
# tuned for speed; consumers that want a fast BLAS use FORTSPARSE_USE_SYSTEM_DEPS=ON
# and supply their own.
#
# Idempotent: only fetches once per configure tree. Other deps (SuperLU's
# external BLAS path, SuiteSparse) consume the `blas`/`lapack` targets and the
# static archive paths exported here.

include(FetchContent)

if(NOT TARGET lapack)
    # Reference-LAPACK honours the standard BUILD_SHARED_LIBS and its own CBLAS
    # / LAPACKE switches. Keep everything static and skip the C wrappers we do
    # not use so the build stays Fortran-only and quick.
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    set(CBLAS OFF CACHE BOOL "" FORCE)
    set(LAPACKE OFF CACHE BOOL "" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON CACHE BOOL "" FORCE)

    FetchContent_Declare(reflapack
        GIT_REPOSITORY https://github.com/Reference-LAPACK/lapack
        GIT_TAG v3.12.1
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(reflapack)
endif()

# Absolute paths to the built static archives, consumed by the SuiteSparse
# ExternalProject (which configures out of tree and cannot see the `blas` /
# `lapack` CMake targets directly). Use the concrete archive paths under the
# Reference-LAPACK build tree rather than a $<TARGET_FILE:...> generator
# expression: ExternalProject_Add re-evaluates CMAKE_ARGS genexps in an
# add_custom_command scope that cannot resolve the `blas` / `lapack` targets
# when fortsparse is a nested FetchContent sub-project, which aborts the
# configure with "No target blas". Reference-LAPACK installs its static
# archives into <build>/lib regardless of nesting.
set(FORTSPARSE_REFLAPACK_BLAS_LIB
    "${reflapack_BINARY_DIR}/lib/libblas${CMAKE_STATIC_LIBRARY_SUFFIX}"
    CACHE INTERNAL "reference BLAS static archive")
set(FORTSPARSE_REFLAPACK_LAPACK_LIB
    "${reflapack_BINARY_DIR}/lib/liblapack${CMAKE_STATIC_LIBRARY_SUFFIX}"
    CACHE INTERNAL "reference LAPACK static archive")
