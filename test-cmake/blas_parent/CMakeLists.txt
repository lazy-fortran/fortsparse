cmake_minimum_required(VERSION 3.22)
project(fortsparse_parent_blas LANGUAGES Fortran C)

option(PARENT_TARGETS_ONLY "Supply only imported dense-kernel targets" ON)
option(PARENT_SYSTEM_DEPS "Find system sparse solvers" ON)
if(PARENT_TARGETS_ONLY)
    # A consumer may supply targets without find_package result variables.
    find_library(parent_blas NAMES openblas blas REQUIRED)
    find_library(parent_lapack NAMES openblas lapack REQUIRED)
    add_library(BLAS::BLAS INTERFACE IMPORTED)
    set_property(TARGET BLAS::BLAS PROPERTY INTERFACE_LINK_LIBRARIES "${parent_blas}")
    add_library(LAPACK::LAPACK INTERFACE IMPORTED)
    set_property(TARGET LAPACK::LAPACK PROPERTY INTERFACE_LINK_LIBRARIES
        "${parent_lapack};BLAS::BLAS")
    # Rediscovery must not replace the parent's targets or require a vendor.
    set(BLA_VENDOR ParentTargetsOnly)
else()
    find_package(LAPACK REQUIRED)
endif()
set(FORTSPARSE_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(FORTSPARSE_USE_SYSTEM_DEPS ${PARENT_SYSTEM_DEPS} CACHE BOOL "" FORCE)
get_filename_component(fortsparse_source "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
add_subdirectory("${fortsparse_source}" fortsparse)

add_executable(parent_blas_solve
    "${fortsparse_source}/test/solver/test_fortsparse_solve_real.f90")
target_link_libraries(parent_blas_solve PRIVATE fortsparse)
enable_testing()
add_test(NAME parent_blas_solve COMMAND parent_blas_solve)
