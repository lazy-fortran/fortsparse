add_executable(bench_spd_native_gpu ${CMAKE_CURRENT_LIST_DIR}/native_driver.f90)
target_link_libraries(bench_spd_native_gpu PRIVATE
    fortsparse fortsparse_spd_fixtures)
include(${CMAKE_CURRENT_LIST_DIR}/process_device_memory.cmake)
fsbench_add_device_sampler(bench_spd_native_gpu)
if(FORTSPARSE_BUILD_TESTING)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        add_test(NAME gpu_spd_front_analysis COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_LIST_DIR}/test_analyze_gpu_spd.py)
        set_tests_properties(gpu_spd_front_analysis PROPERTIES
            LABELS "benchmark;fast" TIMEOUT 1)
    endif()
    foreach(workload IN ITEMS rhs1 rhs8 warm refactor repeated-rhs)
        add_test(NAME gpu_contract_native_${workload}
            COMMAND bench_spd_native_gpu grid2 1 1 ${workload}
                heldout -1 host tuned)
        set_tests_properties(gpu_contract_native_${workload} PROPERTIES
            LABELS "native;benchmark;fast" TIMEOUT 1
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        if(FORTSPARSE_ENABLE_NATIVE_CUDA)
            add_test(NAME gpu_contract_native_cuda_${workload}
                COMMAND bench_spd_native_gpu grid2 1 1 ${workload}
                    heldout 0 host tuned)
            set_tests_properties(gpu_contract_native_cuda_${workload} PROPERTIES
                LABELS "native;cuda;benchmark;fast" TIMEOUT 1
                ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        endif()
    endforeach()
endif()
