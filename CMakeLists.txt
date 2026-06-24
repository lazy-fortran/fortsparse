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
find_package(SuperLU REQUIRED)
find_package(Threads REQUIRED)

add_subdirectory(src)

# The UMFPACK (GPL) backend lives entirely in a separate helper executable. It
# is the only target that links UMFPACK; libfortsparse never does. The library
# reaches it out-of-process over shared memory and named semaphores. Built only
# when UMFPACK is available so the SuperLU-only build stays self-contained.
find_package(UMFPACK)
if(UMFPACK_FOUND)
    add_executable(fortsparse_umfpack_helper
        helper/fortsparse_umfpack_helper.c)
    target_include_directories(fortsparse_umfpack_helper PRIVATE
        ${PROJECT_SOURCE_DIR}/include)
    target_link_libraries(fortsparse_umfpack_helper PRIVATE
        UMFPACK::UMFPACK Threads::Threads)
endif()

if(FORTSPARSE_BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()
