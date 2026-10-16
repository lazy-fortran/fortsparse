add_executable(bench_spd_native_gpu ${CMAKE_CURRENT_LIST_DIR}/native_driver.f90)
target_link_libraries(bench_spd_native_gpu PRIVATE
    fortsparse fortsparse_spd_fixtures)
if(FORTSPARSE_BUILD_TESTING)
    foreach(workload IN ITEMS rhs1 rhs8 warm refactor repeated-rhs)
        add_test(NAME gpu_contract_native_${workload}
            COMMAND bench_spd_native_gpu grid2 1 1 ${workload}
                heldout -1 host tuned)
        set_tests_properties(gpu_contract_native_${workload} PROPERTIES
            LABELS "native;benchmark;fast" TIMEOUT 1
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    endforeach()
endif()
