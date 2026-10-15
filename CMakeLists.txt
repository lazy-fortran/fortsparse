cmake_minimum_required(VERSION 3.22)

project(fortsparse LANGUAGES Fortran C)

# Standard install directory variables (CMAKE_INSTALL_BINDIR etc.).
include(GNUInstallDirs)

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
option(FORTSPARSE_TEST_OOM
    "Also run the address-limit SuperLU allocation-failure test (manual)" OFF)
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

# Install the UMFPACK helper to CMAKE_INSTALL_BINDIR. ON makes a consumer that
# installs its own executables to the standard bin dir get the helper
# co-located on `make install` with no extra call; the self-locating library
# then finds it next to the installed binary. A consumer that installs
# elsewhere uses fortsparse_colocate_helper(<exe> <dir>) instead and can set
# this OFF to keep fortsparse from contributing a rule to its install.
option(FORTSPARSE_INSTALL_HELPER
    "Install the UMFPACK helper to CMAKE_INSTALL_BINDIR" ON)

# MUMPS (CeCILL-C) in-process backend: threaded multifrontal Cholesky/LDL^T/LU
# with METIS nested dissection, factor-once/solve-many and block solves. It is
# linked, never copied, and stays optional: OFF builds nothing MUMPS-related.
# Needs MUMPS built with OpenMP and METIS, found via find_package(MUMPS CONFIG)
# (the scivision/mumps CMake build, installed to a prefix on CMAKE_PREFIX_PATH or
# MUMPS_DIR) or else by cmake/FindMUMPS.cmake (distribution packages, sequential
# or MPI-parallel).
option(FORTSPARSE_ENABLE_MUMPS "Enable the in-process MUMPS backend" OFF)
option(FORTSPARSE_BUILD_BENCHMARKS "Build bounded SPD comparison drivers" OFF)
option(FORTSPARSE_ENABLE_NATIVE_CUDA "Build native CUDA numerical kernels" OFF)

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

# Select dense kernels before any backend package can discover its own BLAS
# or prepend a FindLAPACK module. Reuse the parent's choice when present.
include(cmake/fortsparse_blas.cmake)

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
if(FORTSPARSE_ENABLE_MUMPS)
    # One load only (a second find_package(MUMPS) redefines imported targets).
    # MUMPS::ZMUMPS exists when the installation was built with the complex
    # library; src/ then links it, else the complex entry points are stubs
    # returning FORTSPARSE_BACKEND_UNAVAILABLE. A sequential scivision/mumps
    # config defines MPI::MPI_C and MPI::MPI_Fortran itself (as its mpiseq
    # shim), which fails when a parent project already ran find_package(MPI);
    # the find module, which reuses existing MPI targets, is used then.
    if(NOT TARGET MPI::MPI_C AND NOT TARGET MPI::MPI_Fortran)
        find_package(MUMPS CONFIG QUIET COMPONENTS d)
    endif()
    if(NOT MUMPS_FOUND)
        find_package(MUMPS MODULE REQUIRED)
    endif()
endif()
# Some MUMPS configs rerun a private FindLAPACK even when the imported target
# already exists. Keep the selected libraries for targets using the variables.
set(BLAS_LIBRARIES "${_fortsparse_blas_libraries}")
set(LAPACK_LIBRARIES "${_fortsparse_lapack_libraries}")
# OpenMP threads the CSC matvec; without it the matvec is serial.
option(FORTSPARSE_ENABLE_OPENMP "Thread the CSC matvec with OpenMP" ON)
if(FORTSPARSE_ENABLE_OPENMP OR FORTSPARSE_ENABLE_MUMPS)
    find_package(OpenMP COMPONENTS Fortran)
endif()

if(FORTSPARSE_ENABLE_NATIVE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
endif()
add_subdirectory(src)
include(cmake/fortsparse_native_ordering.cmake)
fortsparse_configure_native_ordering(fortsparse)
if(FORTSPARSE_ENABLE_NATIVE_CUDA)
    include(cmake/FortSparseNativeCUDA.cmake)
endif()

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
    if(FORTSPARSE_INSTALL_HELPER)
        install(TARGETS fortsparse_umfpack_helper
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
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
    if(FORTSPARSE_ENABLE_NATIVE_CUDA)
        fortsparse_add_cuda_tests()
    endif()
endif()

if(FORTSPARSE_BUILD_BENCHMARKS)
    add_subdirectory(bench)
endif()
