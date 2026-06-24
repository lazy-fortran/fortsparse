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

# Zero-system-dependency build. OFF (default) fetches and builds every
# numerical dependency from source -- reference BLAS/LAPACK, SuperLU, and the
# SuiteSparse/UMFPACK chain -- so a consumer needs only a C/Fortran compiler,
# CMake and Ninja. ON uses find_package for SuperLU/UMFPACK/BLAS as before, for
# distros that ship them. The two paths are orthogonal to
# FORTSPARSE_ENABLE_SUPERLU: when SuperLU is disabled nothing SuperLU-related is
# fetched, found, or linked under either setting.
option(FORTSPARSE_USE_SYSTEM_DEPS
    "Use system SuperLU/UMFPACK/BLAS via find_package instead of fetching them"
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

add_subdirectory(src)

# The UMFPACK (GPL) backend lives entirely in a separate helper executable. It
# is the only target that links UMFPACK; libfortsparse never does. The library
# reaches it out-of-process over shared memory and named semaphores.
#
# System path: built only when find_package locates UMFPACK, so a host without
# it still gets the SuperLU-only library. Zero-dependency path: fetch and build
# the SuiteSparse/UMFPACK chain (pointed at the fetched reference BLAS/LAPACK)
# and always build the helper.
if(FORTSPARSE_USE_SYSTEM_DEPS)
    find_package(UMFPACK)
    set(FORTSPARSE_BUILD_HELPER ${UMFPACK_FOUND})
    set(FORTSPARSE_UMFPACK_TARGET UMFPACK::UMFPACK)
else()
    include(deps/reflapack)
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
endif()

if(FORTSPARSE_BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()
