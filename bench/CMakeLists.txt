# Included explicitly by the owning root when benchmark builds are requested.
find_package(Threads REQUIRED)
add_library(fortsparse_spd_fixtures STATIC spd_fixtures.f90 process_metrics.c)
set_target_properties(fortsparse_spd_fixtures PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
target_include_directories(fortsparse_spd_fixtures PUBLIC
    "${CMAKE_CURRENT_BINARY_DIR}/modules")
add_executable(bench_spd bench_spd.f90)
target_link_libraries(bench_spd PRIVATE fortsparse fortsparse_spd_fixtures)
fortsparse_colocate_helper(bench_spd)
include(gpu/native.cmake)
include(gpu/cudss.cmake)
if(FORTSPARSE_BUILD_TESTING)
    add_executable(test_spd_fixtures test_spd_fixtures.f90)
    target_link_libraries(test_spd_fixtures PRIVATE fortsparse_spd_fixtures)
    add_test(NAME spd_fixtures COMMAND test_spd_fixtures)
    set_tests_properties(spd_fixtures PROPERTIES LABELS "benchmark;fast" TIMEOUT 1)
    if(FORTSPARSE_ENABLE_NATIVE_CUDA)
        add_test(NAME spd_native_cuda
            COMMAND bench_spd native-cuda grid2 1 1 rhs1 training)
        set_tests_properties(spd_native_cuda PROPERTIES
            LABELS "native;cuda;benchmark;fast" TIMEOUT 1
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    endif()
    if(FORTSPARSE_ENABLE_MUMPS)
        add_test(NAME spd_mumps_ordering
            COMMAND bench_spd mumps grid3 6 1 rhs1 training)
        set_tests_properties(spd_mumps_ordering PROPERTIES
            LABELS "benchmark;mumps;fast" TIMEOUT 1
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    endif()
    add_executable(test_spd_matrix_market test_spd_matrix_market.f90)
    target_link_libraries(test_spd_matrix_market PRIVATE fortsparse_spd_fixtures)
    add_test(NAME spd_matrix_market COMMAND test_spd_matrix_market)
    set_tests_properties(spd_matrix_market PROPERTIES LABELS "benchmark;fast" TIMEOUT 1)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        if(TARGET fortsparse_umfpack_helper)
            add_test(NAME spd_helper_postexec_memory COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/native_families/test_helper_memory.py
                $<TARGET_FILE:fortsparse_umfpack_helper>)
            set_tests_properties(spd_helper_postexec_memory PROPERTIES
                LABELS "ipc;benchmark;fast" TIMEOUT 1
                ENVIRONMENT "OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1")
        endif()
        add_test(NAME spd_benchmark_memory COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_memory_benchmark.py
            $<TARGET_FILE:bench_spd>)
        set_tests_properties(spd_benchmark_memory PROPERTIES
            LABELS "native;benchmark;memory;fast" TIMEOUT 1
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        add_test(NAME spd_benchmark_analysis COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_benchmark_analysis.py)
        set_tests_properties(spd_benchmark_analysis PROPERTIES
            LABELS "benchmark;fast" TIMEOUT 1)
        add_test(NAME spd_benchmark_report COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_report_spd.py)
        set_tests_properties(spd_benchmark_report PROPERTIES
            LABELS "benchmark;fast" TIMEOUT 1)
    endif()
endif()

set(FORTSPARSE_BENCH_MKL_ROOT "$ENV{MKLROOT}" CACHE PATH
    "Optional oneMKL prefix for the isolated PARDISO benchmark")
if(FORTSPARSE_BENCH_MKL_ROOT)
    find_library(FORTSPARSE_BENCH_MKL_RT NAMES mkl_rt
        PATHS "${FORTSPARSE_BENCH_MKL_ROOT}/lib" NO_DEFAULT_PATH REQUIRED)
    add_executable(bench_spd_pardiso bench_spd_pardiso.f90)
    # No fortsparse, BLAS or LAPACK target: keep Intel/OpenBLAS isolated.
    target_link_libraries(bench_spd_pardiso PRIVATE fortsparse_spd_fixtures
        "${FORTSPARSE_BENCH_MKL_RT}" ${CMAKE_DL_LIBS} Threads::Threads m)
endif()

add_subdirectory(native_families)

add_executable(bench_csr_precision bench_csr_precision.f90)
target_link_libraries(bench_csr_precision PRIVATE fortsparse)
add_executable(bench_mixed_lu bench_mixed_lu.f90)
target_link_libraries(bench_mixed_lu PRIVATE fortsparse)

add_executable(bench_amg bench_amg.f90 process_metrics.c)
target_link_libraries(bench_amg PRIVATE fortsparse_amg_test_support)
