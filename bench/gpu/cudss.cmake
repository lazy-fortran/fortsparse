include_guard(GLOBAL)
set(FORTSPARSE_BENCH_CUDSS_ROOT "" CACHE PATH "External cuDSS CUDA13 SDK prefix")
if(NOT FORTSPARSE_BENCH_CUDSS_ROOT)
    return()
endif()
find_package(CUDAToolkit 13 REQUIRED)
find_path(FORTSPARSE_BENCH_CUDSS_INCLUDE cudss.h
    PATHS "${FORTSPARSE_BENCH_CUDSS_ROOT}/include" NO_DEFAULT_PATH REQUIRED)
find_library(FORTSPARSE_BENCH_CUDSS_LIBRARY cudss
    PATHS "${FORTSPARSE_BENCH_CUDSS_ROOT}/lib" NO_DEFAULT_PATH REQUIRED)
add_executable(bench_spd_cudss
    "${CMAKE_CURRENT_LIST_DIR}/cudss_driver.f90"
    "${CMAKE_CURRENT_LIST_DIR}/cudss_shim.c")
target_compile_features(bench_spd_cudss PRIVATE c_std_11)
target_include_directories(bench_spd_cudss PRIVATE
    "${FORTSPARSE_BENCH_CUDSS_INCLUDE}")
target_compile_definitions(bench_spd_cudss PRIVATE
    FSBENCH_CUDSS_THREADING_LAYER="${FORTSPARSE_BENCH_CUDSS_ROOT}/lib/libcudss_mtlayer_gomp.so")
target_link_libraries(bench_spd_cudss PRIVATE fortsparse_spd_fixtures
    "${FORTSPARSE_BENCH_CUDSS_LIBRARY}" CUDA::cudart ${CMAKE_DL_LIBS})
include("${CMAKE_CURRENT_LIST_DIR}/process_device_memory.cmake")
fsbench_add_device_sampler(bench_spd_cudss)
# Vendor SDK and CUDA runtime are confined to this isolated executable.
set_property(TARGET bench_spd_cudss PROPERTY BUILD_RPATH
    "${FORTSPARSE_BENCH_CUDSS_ROOT}/lib;${CUDAToolkit_LIBRARY_DIR}")
if(FORTSPARSE_BUILD_TESTING)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        add_test(NAME gpu_benchmark_runner COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_LIST_DIR}/test_run_gpu_spd.py")
        set_tests_properties(gpu_benchmark_runner PROPERTIES
            LABELS "benchmark;fast" TIMEOUT 1)
    endif()
    add_test(NAME spd_cudss_host
        COMMAND bench_spd_cudss grid2 1 1 rhs1 training 0 host tuned)
    add_test(NAME spd_cudss_refactor
        COMMAND bench_spd_cudss grid2 1 1 refactor training 0 host tuned)
    add_test(NAME spd_cudss_resident
        COMMAND bench_spd_cudss grid2 1 1 rhs8 training 0 resident tuned)
    # Tiny lower-symmetric MatrixMarket oracles are retained in the build tree.
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/cudss-common.mtx"
        "%%MatrixMarket matrix coordinate real symmetric\n4 4 8\n1 1 4\n2 1 -1\n4 1 -0.25\n2 2 5\n3 2 -0.5\n3 3 6\n4 3 -1\n4 4 7\n")
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/cudss-common.perm" "4\n2 0 3 1\n")
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/cudss-nonpositive.mtx"
        "%%MatrixMarket matrix coordinate real symmetric\n3 3 3\n1 1 2\n2 2 -1\n3 3 3\n")
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/cudss-identity.perm" "3\n0 1 2\n")
    add_test(NAME spd_cudss_common
        COMMAND bench_spd_cudss matrix-market 1 1 rhs1 training 0 host common
            "${CMAKE_CURRENT_BINARY_DIR}/cudss-common.mtx" 1
            "${CMAKE_CURRENT_BINARY_DIR}/cudss-common.perm")
    add_test(NAME spd_cudss_nonpositive
        COMMAND bench_spd_cudss matrix-market 1 1 rhs1 training 0 host common
            "${CMAKE_CURRENT_BINARY_DIR}/cudss-nonpositive.mtx" 1
            "${CMAKE_CURRENT_BINARY_DIR}/cudss-identity.perm" 2)
    set_tests_properties(spd_cudss_nonpositive PROPERTIES
        PASS_REGULAR_EXPRESSION "PASS cuDSS rejected expected nonpositive minor")
    add_test(NAME spd_cudss_warm
        COMMAND bench_spd_cudss matrix-market 1 1 warm training 0 host common
            "${CMAKE_CURRENT_BINARY_DIR}/cudss-common.mtx" 1
            "${CMAKE_CURRENT_BINARY_DIR}/cudss-common.perm")
    set_tests_properties(spd_cudss_host spd_cudss_refactor spd_cudss_resident
        spd_cudss_common spd_cudss_nonpositive spd_cudss_warm
        PROPERTIES LABELS "benchmark;cuda;fast" TIMEOUT 1
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()
