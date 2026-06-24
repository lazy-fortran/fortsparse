cmake_minimum_required(VERSION 3.22)

project(fortsparse LANGUAGES Fortran C)

# Export compile_commands.json so clangd resolves the C shim include paths.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Detect standalone vs subproject build so testing defaults sanely.
set(FORTSPARSE_IS_STANDALONE OFF)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(FORTSPARSE_IS_STANDALONE ON)
endif()

option(FORTSPARSE_BUILD_TESTING "Build tests for fortsparse"
    ${FORTSPARSE_IS_STANDALONE})

# SuperLU (BSD) is the default in-process backend. A consumer that only uses
# the out-of-process UMFPACK-IPC backend can switch it off and then needs no
# SuperLU at all: the shim sources are not compiled, SuperLU is not searched
# for or linked, and the FORTSPARSE_HAVE_SUPERLU macro stays undefined so the
# solver factory reports the backend unavailable instead of failing to link.
option(FORTSPARSE_ENABLE_SUPERLU "Enable the in-process SuperLU backend" ON)

# Zero-dependency build for the niche libraries only. OFF (default) fetches and
# builds SuperLU and the SuiteSparse/UMFPACK chain from source; ON locates them
# via find_package, for distros that ship them. Either way BLAS/LAPACK always
# come from the system: they are universally present and optimized, and reusing
# the same system OpenBLAS as the rest of the stack keeps numerics consistent.
# The two paths are orthogonal to FORTSPARSE_ENABLE_SUPERLU: when SuperLU is
# disabled nothing SuperLU-related is fetched, found, or linked under either
# setting.
option(FORTSPARSE_USE_SYSTEM_DEPS
    "Use system SuperLU/UMFPACK via find_package instead of fetching them"
    OFF)

# Roadmap switches kept inert in v0.1. Enzyme reverse-mode autodiff and device
# offloading land in later releases; the flags document the direction without
# pulling in any toolchain today.
option(FORTSPARSE_ENABLE_ENZYME "Enable Enzyme-based reverse-mode autodiff" OFF)
option(FORTSPARSE_ENABLE_DEVICE_OFFLOAD "Enable GPU/device offloading" OFF)

# Shared library inclusion needs PIC-clean objects.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Pin a modern Fortran standard for free-form, full-feature builds.
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# SuperLU (BSD) is the default in-process backend; the library links only it.
# Provisioned only when the backend is enabled: from the system via
# find_package, or fetched and built from source for the zero-dependency build.
# Either path yields the SuperLU::SuperLU / superlu target src/ links against.
if(FORTSPARSE_ENABLE_SUPERLU)
    if(FORTSPARSE_USE_SYSTEM_DEPS)
        find_package(SuperLU REQUIRED)
    else()
        include(deps/superlu)
    endif()
endif()
find_package(Threads REQUIRED)

# BLAS/LAPACK always come from the system: optimized, universally present, and
# the same OpenBLAS the rest of the stack uses, so numerics stay consistent. The
# fetched UMFPACK and SuperLU link against this, not a from-source BLAS. A parent
# project that already ran find_package(BLAS/LAPACK) (NEO-2 does) leaves the
# result variables set; reuse them instead of searching again, and honour
# BLA_VENDOR so a user can select MKL. REQUIRED because the GPL helper cannot
# link UMFPACK without a dense BLAS.
if(NOT BLAS_FOUND)
    find_package(BLAS REQUIRED)
endif()
if(NOT LAPACK_FOUND)
    find_package(LAPACK REQUIRED)
endif()

add_subdirectory(src)

# The UMFPACK (GPL) backend lives entirely in a separate helper executable. It
# is the only target that links UMFPACK; libfortsparse never does. The library
# reaches it out-of-process over shared memory and named semaphores.
#
# System path: built only when find_package locates UMFPACK, so a host without
# it still gets the SuperLU-only library. Fetched path: fetch and build the
# SuiteSparse/UMFPACK chain from source, linked against the system BLAS/LAPACK,
# and always build the helper.
if(FORTSPARSE_USE_SYSTEM_DEPS)
    find_package(UMFPACK)
    set(FORTSPARSE_BUILD_HELPER ${UMFPACK_FOUND})
    set(FORTSPARSE_UMFPACK_TARGET UMFPACK::UMFPACK)
else()
    include(deps/suitesparse)
    set(FORTSPARSE_BUILD_HELPER ON)
    set(FORTSPARSE_UMFPACK_TARGET SuiteSparse::umfpack)
endif()

if(FORTSPARSE_BUILD_HELPER)
    add_executable(fortsparse_umfpack_helper
        helper/fortsparse_umfpack_helper.c)
    target_include_directories(fortsparse_umfpack_helper PRIVATE
        ${PROJECT_SOURCE_DIR}/include)
    target_link_libraries(fortsparse_umfpack_helper PRIVATE
        ${FORTSPARSE_UMFPACK_TARGET} Threads::Threads)
    # SuiteSparse_config references libm (hypot). With system BLAS on the link
    # line libm is no longer pulled in transitively the way the fetched OpenBLAS
    # dragged it in, so the static SuiteSparse archives leave hypot unresolved.
    # Request libm explicitly; it must follow the archives on the link line,
    # which a separate PRIVATE entry guarantees. Skip on Windows, where the math
    # routines live in the C runtime and no libm exists.
    if(NOT WIN32)
        target_link_libraries(fortsparse_umfpack_helper PRIVATE m)
    endif()
endif()

# Place the UMFPACK helper next to a consumer's executable so fortsparse's
# exe-dir discovery finds it at run time with no env var or PATH change, and
# optionally install it alongside that executable. Call once per executable
# that may issue sparse solves:
#
#   fortsparse_colocate_helper(my_app)                 # build tree only
#   fortsparse_colocate_helper(my_app ${install_dir})  # also install there
#
# A no-op when the helper was not built (a SuperLU-only build links nothing
# GPL and needs no separate process), so consumers can call it unconditionally.
# Available to FetchContent/add_subdirectory consumers once fortsparse has been
# made available; the parent owns its own target names, so it must name each
# executable here rather than fortsparse attaching to them implicitly.
function(fortsparse_colocate_helper target)
    if(NOT TARGET fortsparse_umfpack_helper)
        return()
    endif()
    add_custom_command(TARGET ${target} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:fortsparse_umfpack_helper>
            $<TARGET_FILE_DIR:${target}>
        VERBATIM)
    if(ARGN)
        install(PROGRAMS $<TARGET_FILE:fortsparse_umfpack_helper>
            DESTINATION ${ARGN})
    endif()
endfunction()

if(FORTSPARSE_BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()
