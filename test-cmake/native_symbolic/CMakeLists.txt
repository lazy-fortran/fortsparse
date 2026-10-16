cmake_minimum_required(VERSION 3.20)
project(fortsparse_symbolic_oracle LANGUAGES Fortran)
enable_testing()
add_library(symbolic_test_library STATIC
    ../../src/fortsparse_kinds.f90
    ../../src/fortsparse_status.f90)
include(../../cmake/fortsparse_native_symbolic.cmake)
fortsparse_configure_native_symbolic(symbolic_test_library)
fortsparse_add_native_symbolic_test(symbolic_test_library)
