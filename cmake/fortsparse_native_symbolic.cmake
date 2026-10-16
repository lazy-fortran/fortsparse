include_guard(GLOBAL)

function(fortsparse_configure_native_symbolic target)
    target_sources(${target} PRIVATE
        "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../src/native/fortsparse_native_symbolic.f90")
endfunction()

function(fortsparse_add_native_symbolic_test target)
    add_executable(test_native_symbolic_counts
        "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../test/solver/test_native_symbolic_counts.f90")
    target_link_libraries(test_native_symbolic_counts PRIVATE ${target})
    add_test(NAME test_native_symbolic_counts COMMAND test_native_symbolic_counts)
    set_tests_properties(test_native_symbolic_counts PROPERTIES
        LABELS "native;fast" TIMEOUT 1)
endfunction()
