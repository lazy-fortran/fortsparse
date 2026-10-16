# Explicit opt-in standalone GPL CHOLMOD/GPU comparison; no core dependency.
set(FORTSPARSE_BENCH_CHOLMOD_ROOT "" CACHE PATH
    "Private CUDA-enabled SuiteSparse prefix for isolated CHOLMOD benchmark")
if(FORTSPARSE_BENCH_CHOLMOD_ROOT)
    set(_cholmod_saved_prefix "${CMAKE_PREFIX_PATH}")
    list(PREPEND CMAKE_PREFIX_PATH "${FORTSPARSE_BENCH_CHOLMOD_ROOT}")
    find_package(CHOLMOD CONFIG REQUIRED PATHS
        "${FORTSPARSE_BENCH_CHOLMOD_ROOT}/lib/cmake/CHOLMOD"
        NO_DEFAULT_PATH)
    set(CMAKE_PREFIX_PATH "${_cholmod_saved_prefix}")
    find_package(CUDAToolkit REQUIRED)
    add_executable(bench_spd_cholmod
        "${CMAKE_CURRENT_LIST_DIR}/cholmod_driver.f90"
        "${CMAKE_CURRENT_LIST_DIR}/cholmod_shim.c")
    target_link_libraries(bench_spd_cholmod PRIVATE
        fortsparse_spd_fixtures SuiteSparse::CHOLMOD CUDA::cudart)
    set_target_properties(bench_spd_cholmod PROPERTIES C_STANDARD 11)
endif()
