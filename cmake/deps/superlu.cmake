# Fetch and build sequential SuperLU (BSD) from source, static. Provides the
# `superlu` target. This is the zero-system-dependency replacement for
# find_package(SuperLU): it needs no system BLAS.
#
# SuperLU's own bundled CBLAS also creates a target literally named `blas`,
# which collides with the reference-BLAS `blas` target the SuiteSparse chain
# needs. Both must live in the same CMake tree, so we disable the bundled CBLAS
# and point SuperLU at the fetched reference BLAS instead (still source-built,
# still no system dependency). include(deps/reflapack) must run before this file
# so the `blas` target exists.
#
# The fortsparse SuperLU C shim includes its headers with a `superlu/` prefix
# (e.g. <superlu/slu_ddefs.h>), matching the usual distro layout. The upstream
# build exposes the bare SRC directory on the `superlu` target, so we stage a
# `superlu/` symlink to those headers and attach it to the target's interface
# includes.

include(FetchContent)
include(deps/reflapack)

if(NOT TARGET superlu)
    # Use the fetched reference BLAS rather than SuperLU's bundled CBLAS to keep
    # a single `blas` target in the tree.
    set(enable_internal_blaslib OFF CACHE BOOL "" FORCE)
    set(TPL_BLAS_LIBRARIES "$<TARGET_FILE:blas>" CACHE STRING "" FORCE)
    set(enable_tests OFF CACHE BOOL "" FORCE)
    set(enable_doc OFF CACHE BOOL "" FORCE)
    set(enable_examples OFF CACHE BOOL "" FORCE)
    set(enable_fortran OFF CACHE BOOL "" FORCE)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON CACHE BOOL "" FORCE)

    FetchContent_Declare(superlu
        GIT_REPOSITORY https://github.com/xiaoyeli/superlu
        GIT_TAG v7.0.0
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(superlu)

    # Stage a `superlu/` prefix directory so the shim's <superlu/slu_*defs.h>
    # includes resolve against the upstream SRC headers (which also carry the
    # configured superlu_config.h next to them).
    set(_fortsparse_superlu_incstage "${superlu_BINARY_DIR}/fortsparse_include")
    file(MAKE_DIRECTORY "${_fortsparse_superlu_incstage}")
    if(NOT EXISTS "${_fortsparse_superlu_incstage}/superlu")
        file(CREATE_LINK "${superlu_SOURCE_DIR}/SRC"
            "${_fortsparse_superlu_incstage}/superlu" SYMBOLIC)
    endif()
    target_include_directories(superlu INTERFACE
        "$<BUILD_INTERFACE:${_fortsparse_superlu_incstage}>")

    # SuperLU records the reference BLAS only as a path string, so add the
    # `blas` target as an explicit link dependency. This both resolves its BLAS
    # symbols for downstream linkers and orders `blas` to build first.
    target_link_libraries(superlu PUBLIC blas)
endif()
