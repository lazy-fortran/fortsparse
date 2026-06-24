# Fetch and build sequential SuperLU (BSD) from source, static. Provides the
# `superlu` target. This is the zero-system-dependency replacement for
# find_package(SuperLU): it needs no system BLAS.
#
# SuperLU builds its own bundled CBLAS (enable_internal_blaslib=ON) so it needs
# no external BLAS at all. That keeps SuperLU fully self-contained and decoupled
# from the OpenBLAS the SuiteSparse/UMFPACK chain links: the two backends share
# no BLAS target, so there is no name collision and SuperLU stays buildable even
# when the UMFPACK chain is off.
#
# The fortsparse SuperLU C shim includes its headers with a `superlu/` prefix
# (e.g. <superlu/slu_ddefs.h>), matching the usual distro layout. The upstream
# build exposes the bare SRC directory on the `superlu` target, so we stage a
# `superlu/` symlink to those headers and attach it to the target's interface
# includes.

include(FetchContent)

if(NOT TARGET superlu)
    # Build SuperLU's bundled CBLAS so the backend carries its own BLAS and
    # needs no external library.
    set(enable_internal_blaslib ON CACHE BOOL "" FORCE)
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
endif()
