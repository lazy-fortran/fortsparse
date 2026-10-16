# Fixtures/records link neither a sparse solver nor BLAS: also used by MKL-only.
add_library(fortsparse_family_fixtures STATIC fixtures.f90 records.f90
    real_fixture_io.f90
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
# Optional GPL benchmark executable; fixtures/records remain permissive and
# libfortsparse never links this target or UMFPACK. Independent of IPC enablement.
set(_family_umf_target "")
if(TARGET UMFPACK::UMFPACK)
    set(_family_umf_target UMFPACK::UMFPACK)
elseif(TARGET SuiteSparse::umfpack)
    set(_family_umf_target SuiteSparse::umfpack)
else()
    find_package(UMFPACK QUIET)
    if(TARGET UMFPACK::UMFPACK)
        set(_family_umf_target UMFPACK::UMFPACK)
    endif()
endif()
if(_family_umf_target)
    add_executable(bench_native_families_umfpack
        umfpack_direct/bench_native_families_umfpack.f90
        umfpack_direct/umfpack_direct.c)
    target_link_libraries(bench_native_families_umfpack PRIVATE
        fortsparse_family_fixtures ${_family_umf_target})
    if(NOT WIN32)
        target_link_libraries(bench_native_families_umfpack PRIVATE m)
    endif()
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
    add_executable(test_native_family_real_io test_real_fixture_io.f90)
    target_link_libraries(test_native_family_real_io PRIVATE fortsparse_family_fixtures)
    add_test(NAME native_family_real_io COMMAND test_native_family_real_io)
    set_tests_properties(native_family_real_io PROPERTIES
        LABELS "benchmark;fast" TIMEOUT 1)
    if(Python3_Interpreter_FOUND)
        add_test(NAME native_family_real_audit COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/test_real_fixture.py")
        set_tests_properties(native_family_real_audit PROPERTIES
            LABELS "benchmark;fast" TIMEOUT 1)
        foreach(_peer superlu umfpack mumps umfpack-direct umfpack-direct-noir intel-pardiso)
            set(_real_driver bench_native_families)
            if(_peer STREQUAL "superlu" AND NOT FORTSPARSE_ENABLE_SUPERLU)
                continue()
            elseif(_peer STREQUAL "umfpack" AND NOT FORTSPARSE_BUILD_HELPER)
                continue()
            elseif(_peer STREQUAL "mumps" AND NOT FORTSPARSE_ENABLE_MUMPS)
                continue()
            elseif(_peer MATCHES "^umfpack-direct")
                set(_real_driver bench_native_families_umfpack)
            elseif(_peer STREQUAL "intel-pardiso")
                set(_real_driver bench_native_families_pardiso)
            endif()
            if(TARGET ${_real_driver})
                add_test(NAME "native_family_real_${_peer}"
                    COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/test_real_fixture.py"
                    $<TARGET_FILE:${_real_driver}> "${_peer}")
                set_tests_properties("native_family_real_${_peer}" PROPERTIES
                    LABELS "benchmark;fast" TIMEOUT 1
                    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;MKL_THREADING_LAYER=GNU")
            endif()
        endforeach()
    endif()
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

if(FORTSPARSE_BUILD_TESTING AND TARGET bench_native_families_umfpack)
    foreach(_mode default no-refinement)
        set(_suffix "")
        if(_mode STREQUAL "no-refinement")
            set(_suffix "noir_")
        endif()
        foreach(_operation N T H)
            foreach(_case neo2-mini5 drift-real drift-complex hpd-chain ldlt-chain)
                set(_direct_level 3)
                if(_case STREQUAL "neo2-mini5")
                    set(_direct_level 5)
                endif()
                add_test(NAME "native_family_umfpack_direct_${_suffix}${_case}_${_operation}"
                    COMMAND bench_native_families_umfpack "${_case}" ${_direct_level} 1 rhs8
                        "${_operation}" heldout "${_family_mini5}" "${_mode}")
                set_tests_properties("native_family_umfpack_direct_${_suffix}${_case}_${_operation}"
                    PROPERTIES LABELS "benchmark;umfpack;fast" TIMEOUT 1
                    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
            endforeach()
        endforeach()
        foreach(_workload rhs1 warm refactor)
            foreach(_case drift-real drift-complex)
                add_test(NAME "native_family_umfpack_direct_${_suffix}${_case}_${_workload}"
                    COMMAND bench_native_families_umfpack "${_case}" 3 1
                        "${_workload}" H training "${_family_mini5}" "${_mode}")
                set_tests_properties("native_family_umfpack_direct_${_suffix}${_case}_${_workload}"
                    PROPERTIES LABELS "benchmark;umfpack;fast" TIMEOUT 1
                    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
            endforeach()
        endforeach()
    endforeach()
endif()

if(FORTSPARSE_BUILD_TESTING AND TARGET bench_native_families_umfpack
        AND Python3_Interpreter_FOUND AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
    add_test(NAME native_family_umfpack_direct_records
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/test_umfpack_direct.py"
            $<TARGET_FILE:bench_native_families_umfpack> "${_family_mini5}")
    set_tests_properties(native_family_umfpack_direct_records PROPERTIES
        LABELS "benchmark;umfpack;memory;fast" TIMEOUT 1
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()

if(FORTSPARSE_BUILD_TESTING AND Python3_Interpreter_FOUND)
    add_test(NAME native_family_memory_aliases COMMAND ${Python3_EXECUTABLE}
        "${CMAKE_CURRENT_SOURCE_DIR}/test_memory_aliases.py"
        $<TARGET_FILE:bench_native_families> "${_family_mini5}")
    set_tests_properties(native_family_memory_aliases PROPERTIES
        LABELS "native;benchmark;fast" TIMEOUT 1
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()
