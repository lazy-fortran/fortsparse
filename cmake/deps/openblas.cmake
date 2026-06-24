# Fetch and build OpenBLAS from source (static), giving the target `openblas`
# (alias of `openblas_static`) and a single static archive that carries BOTH the
# BLAS and the bundled LAPACK symbols. This is the perf-critical dense-kernel
# provider for the zero-system-dependency chain: UMFPACK's dense factorisation
# runs on OpenBLAS, where the netlib reference BLAS was both too slow and
# numerically off versus the tuned references.
#
# Build configuration:
#   DYNAMIC_ARCH=1     runtime CPU dispatch -- one binary picks the right kernel
#                      on old Sandy Bridge and new nodes alike, no SIGILL.
#   USE_THREAD=0       single-threaded: deterministic, and no oversubscription
#                      when many per-rank UMFPACK helpers run on one node.
#   BUILD_SHARED_LIBS  OFF -- one static libopenblas.a, nothing to ship/relink.
#   NOFORTRAN=0 +      a Fortran compiler is present, so OpenBLAS compiles its
#   NO_LAPACK=0        bundled LAPACK into the same archive -- both BLAS and
#                      LAPACK symbols resolve from one libopenblas.a.
#
# Idempotent: only fetches once per configure tree. The SuiteSparse
# ExternalProject and the GPL helper consume the concrete archive path exported
# here; libfortsparse never links it.

include(FetchContent)

if(NOT TARGET openblas)
    # A parent project that already ran find_package(BLAS/LAPACK) (NEO-2 does)
    # leaves BLAS_LIBRARIES / LAPACK_LIBRARIES / *_FOUND set in this scope. Left
    # in place they let downstream consumers relink the system BLAS, defeating
    # the zero-system-dependency guarantee. Clear them for the sub-build scope so
    # the OpenBLAS build actually happens and nothing reaches for system OpenBLAS.
    set(BLAS_LIBRARIES "")
    set(LAPACK_LIBRARIES "")
    set(BLAS_FOUND FALSE)
    set(LAPACK_FOUND FALSE)

    # OpenBLAS build switches. DYNAMIC_ARCH and USE_THREAD are OpenBLAS-native
    # CMake options; NO_LAPACK stays off (default) so the bundled LAPACK is built.
    set(DYNAMIC_ARCH ON CACHE BOOL "" FORCE)
    set(USE_THREAD OFF CACHE BOOL "" FORCE)
    set(USE_LOCKING ON CACHE BOOL "" FORCE)
    set(NO_LAPACK OFF CACHE BOOL "" FORCE)
    set(NO_LAPACKE ON CACHE BOOL "" FORCE)
    set(NO_CBLAS ON CACHE BOOL "" FORCE)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(BUILD_WITHOUT_LAPACK OFF CACHE BOOL "" FORCE)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON CACHE BOOL "" FORCE)

    # v0.3.29 raises OpenBLAS's own cmake_minimum_required to 3.16, which keeps
    # it configurable under CMake 4.x; v0.3.28 still declared a pre-3.5 minimum
    # that CMake 4 rejects.
    FetchContent_Declare(openblas
        GIT_REPOSITORY https://github.com/OpenMathLib/OpenBLAS
        GIT_TAG v0.3.29
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(openblas)

    # The bundled LAPACK is Fortran, so the static libopenblas.a carries
    # unresolved Fortran-runtime symbols. The GPL helper is a C executable and
    # links OpenBLAS through this target, so propagate the Fortran runtime and
    # libm as interface link libraries; otherwise the C link driver leaves the
    # LAPACK references unresolved. OpenBLAS does not export these itself.
    # OpenBLAS links openblas_static with the plain signature, so match it (a
    # keyword call here errors). Plain libs propagate to consumers as link
    # interface, which is what the helper needs.
    target_link_libraries(openblas_static gfortran m)
endif()

# Absolute path to the built static archive, consumed by the SuiteSparse
# ExternalProject (which configures out of tree and cannot see the `openblas`
# CMake target) and by the GPL helper. Use the concrete archive path under the
# OpenBLAS build tree rather than a $<TARGET_FILE:openblas> generator
# expression: ExternalProject_Add re-evaluates CMAKE_ARGS genexps in an
# add_custom_command scope that cannot resolve the `openblas` target when
# fortsparse is a nested FetchContent sub-project, which aborts the configure.
# The single archive holds both BLAS and LAPACK symbols, so both variables point
# at it. CACHE INTERNAL so the path survives nested-scope re-evaluation.
set(FORTSPARSE_OPENBLAS_LIB
    "${openblas_BINARY_DIR}/lib/libopenblas${CMAKE_STATIC_LIBRARY_SUFFIX}"
    CACHE INTERNAL "OpenBLAS static archive (BLAS + bundled LAPACK)")
