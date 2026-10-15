# Included after creating the fortsparse target, only when CUDA is requested.
# Ordinary CPU builds never need CUDA. Architectures follow the caller's
# CMAKE_CUDA_ARCHITECTURES, so installed binaries can target declared runners.
include(CheckLanguage)
check_language(CUDA)
if(NOT CMAKE_CUDA_COMPILER)
    message(FATAL_ERROR "FortSparse native CUDA requested but no CUDA compiler found")
endif()
enable_language(CUDA)
find_package(CUDAToolkit REQUIRED)
target_sources(fortsparse PRIVATE
    ${PROJECT_SOURCE_DIR}/src/native/fsparse_spd_cuda.cu)
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
endfunction()
