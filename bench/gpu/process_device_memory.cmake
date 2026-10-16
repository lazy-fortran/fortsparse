include_guard(GLOBAL)
function(fsbench_add_device_sampler target)
    target_sources(${target} PRIVATE
        "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/process_device_memory.c")
    if(TARGET CUDA::cudart)
        find_path(FSBENCH_NVML_INCLUDE_DIR nvml.h
            HINTS ${CUDAToolkit_INCLUDE_DIRS})
        if(FSBENCH_NVML_INCLUDE_DIR)
            target_compile_definitions(${target} PRIVATE FSBENCH_HAVE_NVML)
            target_include_directories(${target} PRIVATE "${FSBENCH_NVML_INCLUDE_DIR}")
            target_link_libraries(${target} PRIVATE CUDA::cudart ${CMAKE_DL_LIBS})
        endif()
    endif()
endfunction()
