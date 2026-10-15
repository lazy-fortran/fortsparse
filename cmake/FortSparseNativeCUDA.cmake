# Included after creating the fortsparse target, only when CUDA is requested.
# Ordinary CPU builds never need CUDA. Architectures follow the caller's
# CMAKE_CUDA_ARCHITECTURES, so installed binaries can target declared runners.
# CUDA is enabled in the root before creating src/ so that directory inherits
# the language's compile rules. The optional target is configured here.
target_sources(fortsparse PRIVATE
    ${PROJECT_SOURCE_DIR}/src/native/fsparse_spd_cuda.cu
    ${PROJECT_SOURCE_DIR}/src/native/fortsparse_native_cuda.F90)
target_compile_definitions(fortsparse PRIVATE FORTSPARSE_HAVE_NATIVE_CUDA)
target_include_directories(fortsparse PRIVATE ${PROJECT_SOURCE_DIR}/src/generated)
target_link_libraries(fortsparse PUBLIC CUDA::cudart CUDA::cublas)
set_target_properties(fortsparse PROPERTIES CUDA_STANDARD 14 CUDA_STANDARD_REQUIRED ON)

function(fortsparse_add_cuda_tests)
    foreach(test IN ITEMS test_fortsparse_spd_cuda_leaves test_fortsparse_spd_cuda)
        add_executable(${test} ${PROJECT_SOURCE_DIR}/test-cmake/cuda/${test}.cu)
        target_link_libraries(${test} PRIVATE fortsparse)
        target_include_directories(${test} PRIVATE
            ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src/generated)
        set_target_properties(${test} PROPERTIES
            CUDA_STANDARD 14 CUDA_STANDARD_REQUIRED ON)
        add_test(NAME ${test} COMMAND ${test})
        set_tests_properties(${test} PROPERTIES LABELS "native;cuda;fast" TIMEOUT 1)
    endforeach()
    add_executable(test_fortsparse_native_cuda_binding
        ${PROJECT_SOURCE_DIR}/test-cmake/cuda/test_fortsparse_native_cuda_binding.f90)
    target_link_libraries(test_fortsparse_native_cuda_binding PRIVATE fortsparse)
    add_test(NAME test_fortsparse_native_cuda_binding
        COMMAND test_fortsparse_native_cuda_binding)
    set_tests_properties(test_fortsparse_native_cuda_binding PROPERTIES
        LABELS "native;cuda;fast" TIMEOUT 1)
endfunction()
