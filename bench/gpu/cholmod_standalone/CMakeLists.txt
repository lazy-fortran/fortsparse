cmake_minimum_required(VERSION 3.23)
project(FortSparseCHOLMODComparison LANGUAGES C Fortran)
# This independent build never configures or links the FortSparse library.
set(bench_root "${CMAKE_CURRENT_LIST_DIR}/../..")
add_library(fortsparse_spd_fixtures STATIC
    "${bench_root}/spd_fixtures.f90" "${bench_root}/process_metrics.c")
set_target_properties(fortsparse_spd_fixtures PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
target_include_directories(fortsparse_spd_fixtures PUBLIC
    "${CMAKE_CURRENT_BINARY_DIR}/modules")
include("${CMAKE_CURRENT_LIST_DIR}/../cholmod.cmake")

enable_testing()
foreach(workload rhs1 rhs8 refactor warm repeated-rhs)
    add_test(NAME gpu_cholmod_${workload}
        COMMAND bench_spd_cholmod grid2 1 1 ${workload} training 0 host tuned - 3)
    set_tests_properties(gpu_cholmod_${workload} PROPERTIES
        LABELS "benchmark;cuda;fast" TIMEOUT 1
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;CHOLMOD_GPU_POOL_BYTES=268435456")
endforeach()
set(common_order "${CMAKE_CURRENT_BINARY_DIR}/identity-order.txt")
file(WRITE "${common_order}" "128\n")
foreach(i RANGE 0 127)
    file(APPEND "${common_order}" "${i}\n")
endforeach()
add_test(NAME gpu_cholmod_common_order
    COMMAND bench_spd_cholmod band 1 1 rhs1 training 0 host common - 1
        "${common_order}")
set_tests_properties(gpu_cholmod_common_order PROPERTIES
    LABELS "benchmark;cuda;fast" TIMEOUT 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;CHOLMOD_GPU_POOL_BYTES=268435456")

find_package(Python3 COMPONENTS Interpreter QUIET)
set(FORTSPARSE_BENCH_GPU_RUNNER
    "${CMAKE_CURRENT_LIST_DIR}/../run_gpu_spd.py" CACHE FILEPATH
    "Shared runner used to validate actual CHOLMOD CSV integration")
if(Python3_Interpreter_FOUND AND EXISTS "${FORTSPARSE_BENCH_GPU_RUNNER}")
    foreach(workload refactor warm)
        add_test(NAME gpu_cholmod_csv_${workload}
            COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_LIST_DIR}/../cholmod_csv_test.py"
                --driver $<TARGET_FILE:bench_spd_cholmod>
                --runner "${FORTSPARSE_BENCH_GPU_RUNNER}"
                --workload ${workload})
        set_tests_properties(gpu_cholmod_csv_${workload} PROPERTIES
            LABELS "benchmark;cuda;fast" TIMEOUT 1
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;CHOLMOD_GPU_POOL_BYTES=268435456")
    endforeach()
endif()
