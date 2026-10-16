cmake_minimum_required(VERSION 3.22)
project(fortsparse_native_only LANGUAGES Fortran C)

# Exercise the default fetching mode with every sparse dependency disabled.
# Disconnected FetchContent and disabled discovery make accidental dependence
# fail at configure time instead of succeeding because this host has packages.
set(FETCHCONTENT_FULLY_DISCONNECTED ON CACHE BOOL "" FORCE)
set(CMAKE_DISABLE_FIND_PACKAGE_SuperLU ON)
set(CMAKE_DISABLE_FIND_PACKAGE_UMFPACK ON)
set(CMAKE_DISABLE_FIND_PACKAGE_AMD ON)
set(CMAKE_DISABLE_FIND_PACKAGE_METIS ON)
set(FORTSPARSE_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(FORTSPARSE_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
set(FORTSPARSE_USE_SYSTEM_DEPS OFF CACHE BOOL "" FORCE)
foreach(option IN ITEMS SUPERLU UMFPACK_IPC MUMPS AMD METIS NATIVE_CUDA)
    set(FORTSPARSE_ENABLE_${option} OFF CACHE BOOL "" FORCE)
endforeach()
get_filename_component(fortsparse_source "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
add_subdirectory("${fortsparse_source}" fortsparse)
foreach(provider IN ITEMS fortsparse_umfpack_helper SuperLU::SuperLU superlu
    UMFPACK::UMFPACK SuiteSparse::umfpack MUMPS::DMUMPS MUMPS::ZMUMPS
    fortsparse_amd SuiteSparse::AMD SuiteSparse::AMD_static AMD::AMD amd
    fortsparse_metis METIS::METIS METIS::metis metis METIS fortsym fortgen)
    if(TARGET ${provider})
        message(FATAL_ERROR "Native-only build unexpectedly provides ${provider}")
    endif()
endforeach()

add_executable(native_only_smoke smoke.f90)
target_link_libraries(native_only_smoke PRIVATE fortsparse)
enable_testing()
add_test(NAME native_only_smoke COMMAND native_only_smoke)
set_tests_properties(native_only_smoke PROPERTIES TIMEOUT 1
    ENVIRONMENT "FORTSPARSE_UMFPACK_HELPER=${CMAKE_CURRENT_BINARY_DIR}/no-helper")
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    add_test(NAME native_only_links COMMAND ${CMAKE_COMMAND}
        -DEXECUTABLE=$<TARGET_FILE:native_only_smoke>
        -P ${CMAKE_CURRENT_SOURCE_DIR}/check_links.cmake)
    set_tests_properties(native_only_links PROPERTIES TIMEOUT 1)
endif()
