# Fixtures/records link neither a sparse solver nor BLAS: also used by MKL-only.
add_library(fortsparse_family_fixtures STATIC fixtures.f90 records.f90
    ../process_metrics.c)
set_target_properties(fortsparse_family_fixtures PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
target_include_directories(fortsparse_family_fixtures PUBLIC
    "${CMAKE_CURRENT_BINARY_DIR}/modules")
add_executable(bench_native_families bench_native_families.f90)
target_link_libraries(bench_native_families PRIVATE fortsparse
    fortsparse_family_fixtures)
fortsparse_colocate_helper(bench_native_families)
set(_family_mini5 "${CMAKE_CURRENT_SOURCE_DIR}/../fixtures/neo2-mini5/matrix.mtx")
if(FORTSPARSE_BENCH_MKL_ROOT)
    add_executable(bench_native_families_pardiso bench_native_families_pardiso.f90)
    target_link_libraries(bench_native_families_pardiso PRIVATE
        fortsparse_family_fixtures "${FORTSPARSE_BENCH_MKL_RT}"
        ${CMAKE_DL_LIBS} Threads::Threads m)
endif()
if(FORTSPARSE_BUILD_TESTING)
    if(Python3_Interpreter_FOUND)
        add_test(NAME native_family_selector COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/test_selector.py")
        set_tests_properties(native_family_selector PROPERTIES
            LABELS "benchmark;fast" TIMEOUT 1)
    endif()
    add_executable(test_native_family_fixtures test_fixtures.f90)
    target_link_libraries(test_native_family_fixtures PRIVATE
        fortsparse_family_fixtures)
    add_test(NAME native_family_fixtures COMMAND test_native_family_fixtures
        "${_family_mini5}")
    set_tests_properties(native_family_fixtures PROPERTIES
        LABELS "benchmark;fast" TIMEOUT 1)
    foreach(_case neo2-mini5 drift-real drift-complex hpd-chain ldlt-chain)
        set(_level 3)
        if(_case STREQUAL "neo2-mini5")
            set(_level 5)
        endif()
        add_test(NAME "native_family_${_case}" COMMAND bench_native_families
            native "${_case}" "${_level}" 1 rhs8 H heldout "${_family_mini5}")
        set_tests_properties("native_family_${_case}" PROPERTIES
            LABELS "native;benchmark;fast" TIMEOUT 1
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    endforeach()
    if(TARGET bench_native_families_pardiso)
        add_test(NAME native_family_pardiso_transpose
            COMMAND bench_native_families_pardiso hpd-chain 3 1 rhs8 T heldout
                "${_family_mini5}")
        set_tests_properties(native_family_pardiso_transpose PROPERTIES
            LABELS "benchmark;pardiso;fast" TIMEOUT 1
            ENVIRONMENT "MKL_NUM_THREADS=1;OMP_NUM_THREADS=1;MKL_THREADING_LAYER=GNU")
        foreach(_case neo2-mini5 drift-complex ldlt-chain)
            add_test(NAME "native_family_pardiso_${_case}"
                COMMAND bench_native_families_pardiso "${_case}" 3 1 refactor H
                    training "${_family_mini5}")
            set_tests_properties("native_family_pardiso_${_case}" PROPERTIES
                LABELS "benchmark;pardiso;fast" TIMEOUT 1
                ENVIRONMENT "MKL_NUM_THREADS=1;OMP_NUM_THREADS=1;MKL_THREADING_LAYER=GNU")
        endforeach()
    endif()
endif()
