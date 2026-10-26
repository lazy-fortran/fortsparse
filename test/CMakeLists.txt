# fortsparse test suite. Each test is registered only when its source file is
# present, so the library configures and builds before the test programs are
# written and every test activates automatically once its source lands.
# Tests that fpm cannot build (MUMPS/MPI, the C IPC driver) live in
# test-cmake/, which fpm does not scan. Their CTest registrations are unchanged.

function(fortsparse_add_test name source)
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
        return()
    endif()
    add_executable(${name} ${source})
    target_link_libraries(${name} PRIVATE fortsparse)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

fortsparse_add_native_symbolic_test(fortsparse)

fortsparse_add_test(test_fortsparse_smoke test_fortsparse_smoke.f90)
fortsparse_add_test(test_fortsparse_kinds test_fortsparse_kinds.f90)
fortsparse_add_test(test_fortsparse_status test_fortsparse_status.f90)
fortsparse_add_test(test_fortsparse_memory solver/test_fortsparse_memory.f90)
set_tests_properties(test_fortsparse_memory PROPERTIES LABELS "native;fast" TIMEOUT 1)
fortsparse_add_test(test_fortsparse_spd_leaves
    solver/test_fortsparse_spd_leaves.f90)
fortsparse_add_test(test_fortsparse_native_lu
    solver/test_fortsparse_native_lu.f90)
fortsparse_add_test(test_fortsparse_ilutp solver/test_fortsparse_ilutp.f90)
set_tests_properties(test_fortsparse_ilutp PROPERTIES
    LABELS "native;preconditioner;fast" TIMEOUT 1)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(test_ilutp_allocations
        ../test-cmake/native/test_ilutp_allocations.f90
        ../test-cmake/native/ilutp_alloc_probe.c)
    target_link_libraries(test_ilutp_allocations PRIVATE fortsparse)
    target_link_options(test_ilutp_allocations PRIVATE "LINKER:--wrap=malloc"
        "LINKER:--wrap=calloc" "LINKER:--wrap=realloc")
    add_test(NAME test_ilutp_allocations COMMAND test_ilutp_allocations)
    set_tests_properties(test_ilutp_allocations PROPERTIES
        LABELS "native;preconditioner;memory;fast" TIMEOUT 1)
endif()
fortsparse_add_test(test_fortsparse_native_lu_lifetime
    ../test-cmake/native/test_native_lu_lifetime.f90)
target_sources(test_fortsparse_native_lu_lifetime PRIVATE ../test-cmake/native/lu_discard_counter.c)
set_tests_properties(test_fortsparse_native_lu_lifetime PROPERTIES
    LABELS "native;lu;fast" TIMEOUT 1 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_definitions(test_fortsparse_native_lu_lifetime PRIVATE FS_TEST_WRAP_ADVICE)
    target_link_options(test_fortsparse_native_lu_lifetime PRIVATE "LINKER:--wrap=madvise" "LINKER:--wrap=malloc"
        "LINKER:--wrap=calloc" "LINKER:--wrap=realloc" "LINKER:--wrap=free")
    add_test(NAME test_fortsparse_native_lu_advice_failure COMMAND test_fortsparse_native_lu_lifetime 1)
    set_tests_properties(test_fortsparse_native_lu_advice_failure PROPERTIES
        LABELS "native;lu;fast" TIMEOUT 1 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()
# Always compile the no-advice fallback. Kernel-call wrapping is Linux-only.
add_executable(test_native_lu_discard_portable
    ../test-cmake/native/test_native_lu_discard.c ../src/native/fsparse_native_lu_memory.c)
target_compile_definitions(test_native_lu_discard_portable PRIVATE
    FORTSPARSE_LU_DISCARD_DISABLED FS_TEST_DISCARD_DISABLED)
add_test(NAME test_native_lu_discard_portable COMMAND test_native_lu_discard_portable)
set_tests_properties(test_native_lu_discard_portable PROPERTIES LABELS "native;lu;fast" TIMEOUT 1)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(test_native_lu_discard
        ../test-cmake/native/test_native_lu_discard.c ../src/native/fsparse_native_lu_memory.c)
    target_compile_definitions(test_native_lu_discard PRIVATE FS_TEST_WRAP_ADVICE)
    target_link_options(test_native_lu_discard PRIVATE "LINKER:--wrap=madvise")
    add_test(NAME test_native_lu_discard COMMAND test_native_lu_discard)
    set_tests_properties(test_native_lu_discard PROPERTIES LABELS "native;lu;fast" TIMEOUT 1)
endif()
fortsparse_add_test(test_fortsparse_lu_leaves
    solver/test_fortsparse_lu_leaves.f90)
fortsparse_add_test(test_fortsparse_lu_pruning
    solver/test_fortsparse_lu_pruning.f90)
set_tests_properties(test_fortsparse_lu_pruning PROPERTIES
    LABELS "native;lu;fast" TIMEOUT 1)
set_tests_properties(test_fortsparse_native_lu test_fortsparse_lu_leaves PROPERTIES
    LABELS "native;lu;fast" TIMEOUT 1)
fortsparse_add_test(test_fortsparse_hpd_leaves
    solver/test_fortsparse_hpd_leaves.f90)
set_tests_properties(test_fortsparse_hpd_leaves PROPERTIES
    LABELS "native;fast" TIMEOUT 1)
fortsparse_add_test(test_fortsparse_native_hpd
    solver/test_fortsparse_native_hpd.f90)
set_tests_properties(test_fortsparse_native_hpd PROPERTIES
    LABELS "native;fast" TIMEOUT 1 PROCESSORS 2
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
fortsparse_add_test(test_fortsparse_native
    solver/test_fortsparse_native.f90)
fortsparse_add_test(test_fortsparse_native_supernodes
    solver/test_fortsparse_native_supernodes.f90)
fortsparse_add_test(test_fortsparse_supernode_partition
    solver/test_fortsparse_supernode_partition.f90)
set_tests_properties(test_fortsparse_supernode_partition PROPERTIES
    LABELS "native;fast" TIMEOUT 1)
fortsparse_add_test(test_fortsparse_slu_analysis
    solver/test_fortsparse_slu_analysis.f90)
set_tests_properties(test_fortsparse_slu_analysis PROPERTIES
    LABELS "native;lu;fast" TIMEOUT 1)
fortsparse_add_test(test_fortsparse_slu_leaves
    solver/test_fortsparse_slu_leaves.f90)
set_tests_properties(test_fortsparse_slu_leaves PROPERTIES
    LABELS "native;lu;fast" TIMEOUT 1)
fortsparse_add_test(test_fortsparse_slu_factor
    solver/test_fortsparse_slu_factor.f90)
set_tests_properties(test_fortsparse_slu_factor PROPERTIES
    LABELS "native;lu;fast" TIMEOUT 1)
fortsparse_add_test(test_fortsparse_slu_z_leaves
    solver/test_fortsparse_slu_z_leaves.f90)
set_tests_properties(test_fortsparse_slu_z_leaves PROPERTIES
    LABELS "native;lu;fast" TIMEOUT 1)
fortsparse_add_test(test_fortsparse_slu_z_factor
    solver/test_fortsparse_slu_z_factor.f90)
set_tests_properties(test_fortsparse_slu_z_factor PROPERTIES
    LABELS "native;lu;fast" TIMEOUT 1)
fortsparse_add_test(test_fortsparse_slu_z_api
    solver/test_fortsparse_slu_z_api.f90)
set_tests_properties(test_fortsparse_slu_z_api PROPERTIES
    LABELS "native;lu;fast" TIMEOUT 1
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
fortsparse_add_test(test_fortsparse_slu_api
    solver/test_fortsparse_slu_api.f90)
set_tests_properties(test_fortsparse_slu_api PROPERTIES
    LABELS "native;lu;fast" TIMEOUT 1)
fortsparse_add_test(test_fortsparse_slu_fallback
    solver/test_fortsparse_slu_fallback.f90)
set_tests_properties(test_fortsparse_slu_fallback PROPERTIES
    LABELS "native;lu;fast" TIMEOUT 1)
fortsparse_add_test(test_fortsparse_slu_threads
    solver/test_fortsparse_slu_threads.f90)
set_tests_properties(test_fortsparse_slu_threads PROPERTIES
    LABELS "native;lu;fast" TIMEOUT 1
    ENVIRONMENT "OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1")
if(TARGET test_fortsparse_slu_threads)
    # Exercise parallel analysis independently of the audited AMD cap ABI.
    # The original entry retains every cap and LOW_MEMORY assertion.
    add_test(NAME test_fortsparse_slu_threads_uncapped
        COMMAND test_fortsparse_slu_threads uncapped)
    set_tests_properties(test_fortsparse_slu_threads_uncapped PROPERTIES
        LABELS "native;lu;fast" TIMEOUT 1
        ENVIRONMENT "OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=8")
    add_test(NAME test_fortsparse_slu_threads_repeat
        COMMAND test_fortsparse_slu_threads repeat)
    set_tests_properties(test_fortsparse_slu_threads_repeat PROPERTIES
        LABELS "native;lu;fast" TIMEOUT 1
        ENVIRONMENT "OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1")
endif()
fortsparse_add_test(test_fortsparse_supernode_memory
    solver/test_fortsparse_supernode_memory.f90)
set_tests_properties(test_fortsparse_supernode_memory PROPERTIES
    LABELS "native;fast" TIMEOUT 1)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_sources(test_fortsparse_supernode_memory PRIVATE
        ../test-cmake/native/forbid_supernode_blas.c)
    target_link_options(test_fortsparse_supernode_memory PRIVATE
        "LINKER:--wrap=dgemm_" "LINKER:--wrap=dsyrk_" "LINKER:--wrap=dgemv_"
        "LINKER:--wrap=dtrsv_" "LINKER:--wrap=dtrsm_")
endif()
set_tests_properties(test_fortsparse_spd_leaves test_fortsparse_native
    test_fortsparse_native_supernodes PROPERTIES LABELS "native;fast" TIMEOUT 1)
add_executable(test_native_symbolic
    ${PROJECT_SOURCE_DIR}/cmake/tests/test_native_symbolic.c)
target_link_libraries(test_native_symbolic PRIVATE fortsparse)
add_test(NAME test_native_symbolic COMMAND test_native_symbolic)
set_tests_properties(test_native_symbolic PROPERTIES LABELS "native;fast" TIMEOUT 1)
add_executable(test_native_symbolic_portable
    ${PROJECT_SOURCE_DIR}/cmake/tests/test_native_symbolic.c
    ../src/native/fsparse_native_symbolic.c)
target_compile_definitions(test_native_symbolic_portable PRIVATE FORTSPARSE_SYMBOLIC_MALLOC_ONLY)
add_test(NAME test_native_symbolic_portable COMMAND test_native_symbolic_portable)
set_tests_properties(test_native_symbolic_portable PROPERTIES LABELS "native;fast" TIMEOUT 1)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(test_native_symbolic_failure ../test-cmake/native/test_native_symbolic_failure.f90
        ../test-cmake/native/native_symbolic_fail.c)
    target_link_libraries(test_native_symbolic_failure PRIVATE fortsparse)
    target_link_options(test_native_symbolic_failure PRIVATE
        "LINKER:--wrap=fsparse_native_symbolic_alloc")
    add_test(NAME test_native_symbolic_failure COMMAND test_native_symbolic_failure)
    set_tests_properties(test_native_symbolic_failure PROPERTIES LABELS "native;fast" TIMEOUT 1)
endif()
add_executable(test_native_order
    ${PROJECT_SOURCE_DIR}/cmake/tests/test_native_order.c)
target_link_libraries(test_native_order PRIVATE fortsparse)
if(FORTSPARSE_NATIVE_METIS_FOUND)
    target_compile_definitions(test_native_order PRIVATE FORTSPARSE_TEST_METIS)
endif()
add_test(NAME test_native_order COMMAND test_native_order)
set_tests_properties(test_native_order PROPERTIES LABELS "native;fast" TIMEOUT 1)
add_executable(test_native_order_amd
    ${PROJECT_SOURCE_DIR}/cmake/tests/test_native_order_amd.c)
target_link_libraries(test_native_order_amd PRIVATE fortsparse)
if(FORTSPARSE_NATIVE_AMD_FOUND)
    target_compile_definitions(test_native_order_amd PRIVATE FORTSPARSE_TEST_AMD)
endif()
add_test(NAME test_native_order_amd COMMAND test_native_order_amd)
set_tests_properties(test_native_order_amd PROPERTIES LABELS "native;fast" TIMEOUT 1)
fortsparse_add_test(test_fortsparse_csc csc/test_fortsparse_csc.f90)
fortsparse_add_test(test_fortsparse_matvec_omp csc/test_fortsparse_matvec_omp.f90)

if(FORTSPARSE_USE_SYSTEM_DEPS AND FORTSPARSE_ENABLE_SUPERLU)
    add_test(NAME test_fortsparse_parent_blas
        COMMAND ${CMAKE_CTEST_COMMAND} --build-and-test
            ${PROJECT_SOURCE_DIR}/test-cmake/blas_parent
            ${PROJECT_BINARY_DIR}/test-parent-blas
            --build-generator "${CMAKE_GENERATOR}"
            --build-options
                -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            --test-command ${CMAKE_CTEST_COMMAND} -R ^parent_blas_solve$
                --output-on-failure)
    set_tests_properties(test_fortsparse_parent_blas PROPERTIES
        ENVIRONMENT "CMAKE_BUILD_PARALLEL_LEVEL=2;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1"
        PROCESSORS 2)
endif()

if(NOT FORTSPARSE_USE_SYSTEM_DEPS AND FORTSPARSE_ENABLE_SUPERLU)
    add_test(NAME test_fortsparse_fetched_blas_reconfigure
        COMMAND ${CMAKE_CTEST_COMMAND} --build-and-test
            ${PROJECT_SOURCE_DIR}/test-cmake/blas_parent
            ${PROJECT_BINARY_DIR}/test-fetched-blas
            --build-generator "${CMAKE_GENERATOR}"
            --build-two-config --build-target parent_blas_solve
            --build-options
                -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                -DPARENT_TARGETS_ONLY=OFF -DPARENT_SYSTEM_DEPS=OFF
            --test-command ${CMAKE_CTEST_COMMAND} -R ^parent_blas_solve$
                --output-on-failure)
    set_tests_properties(test_fortsparse_fetched_blas_reconfigure PROPERTIES
        ENVIRONMENT "CMAKE_BUILD_PARALLEL_LEVEL=2;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1"
        PROCESSORS 2)
endif()

# Strong-scaling benchmark of csc_matvec; built, not registered as a test.
add_executable(bench_fortsparse_matvec ../bench/bench_csc_matvec.f90)
target_link_libraries(bench_fortsparse_matvec PRIVATE fortsparse)

if(UNIX)
    add_executable(test_fortsparse_ipc_wait
        ../test-cmake/ipc/test_fortsparse_ipc_wait.c)
    target_include_directories(test_fortsparse_ipc_wait PRIVATE
        ${PROJECT_SOURCE_DIR}/include)
    target_link_libraries(test_fortsparse_ipc_wait PRIVATE fortsparse)
    add_test(NAME test_fortsparse_ipc_wait COMMAND test_fortsparse_ipc_wait)
    set_tests_properties(test_fortsparse_ipc_wait PROPERTIES TIMEOUT 10)
endif()

# These cases factor through the default backend (MUMPS in a sequential MUMPS
# build with ZMUMPS, SuperLU otherwise); the SuperLU-specific ones pin SuperLU.
# They are registered when the SuperLU backend is compiled in. The IPC, selflocate, kinds,
# status, and csc tests do not depend on it and stay registered always.
if(FORTSPARSE_ENABLE_SUPERLU)
    fortsparse_add_test(test_fortsparse_solve_real
        solver/test_fortsparse_solve_real.f90)
    fortsparse_add_test(test_fortsparse_solve_derivatives
        solver/test_fortsparse_solve_derivatives.f90)
    fortsparse_add_test(test_fortsparse_solve_complex
        solver/test_fortsparse_solve_complex.f90)
    fortsparse_add_test(test_fortsparse_complex_solve_derivatives
        solver/test_fortsparse_complex_solve_derivatives.f90)
    fortsparse_add_test(test_fortsparse_poisson
        solver/test_fortsparse_poisson.f90)
    fortsparse_add_test(test_fortsparse_superlu_symmetric
        solver/test_fortsparse_superlu_symmetric.f90)
    fortsparse_add_test(test_fortsparse_superlu_refactor
        solver/test_fortsparse_superlu_refactor.f90)
    fortsparse_add_test(test_fortsparse_errors solver/test_fortsparse_errors.f90)
    # Forces a real SuperLU allocation failure under an address-space limit.
    fortsparse_add_test(test_fortsparse_superlu_status
        solver/test_fortsparse_superlu_status.F90)
    target_compile_definitions(test_fortsparse_superlu_status PRIVATE
        FORTSPARSE_HAVE_SUPERLU)
    if(TARGET test_fortsparse_superlu_status AND UNIX AND FORTSPARSE_TEST_OOM)
        # Opt-in (-DFORTSPARSE_TEST_OOM=ON): the window of working address
        # limits is toolchain dependent and SuperLU can abort uncatchably just
        # above it, so this is not part of the default suite. Reruns the test
        # under an address-space limit (the default registration above covers
        # the mapping cases only). Single-threaded BLAS/OpenMP so the
        # limit is not eaten by per-thread buffers.
        add_test(NAME test_fortsparse_superlu_status_oom
            COMMAND sh -c "ulimit -v 250000 && exec \"$0\""
            $<TARGET_FILE:test_fortsparse_superlu_status>)
        set_tests_properties(test_fortsparse_superlu_status_oom PROPERTIES
            LABELS oom
            ENVIRONMENT
            "FORTSPARSE_TEST_ULIMIT_V=1;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1")
    endif()
    # The no-helper error path must run where no helper sits next to it,
    # otherwise exe-directory discovery would resolve one and the backend would
    # be available. Keep this executable out of the bin dir that holds
    # fortsparse_umfpack_helper.
    if(TARGET test_fortsparse_errors)
        set_target_properties(test_fortsparse_errors PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/no_helper")
    endif()
endif()

# IPC backend round-trip through the GPL helper. Registered only when the
# helper target exists; the test needs the helper path in the environment.
if(TARGET fortsparse_umfpack_helper AND
   EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/ipc/test_fortsparse_umfpack_ipc.f90")
    add_executable(test_fortsparse_umfpack_ipc ipc/test_fortsparse_umfpack_ipc.f90)
    target_link_libraries(test_fortsparse_umfpack_ipc PRIVATE fortsparse)
    add_test(NAME test_fortsparse_umfpack_ipc
        COMMAND test_fortsparse_umfpack_ipc)
    set_tests_properties(test_fortsparse_umfpack_ipc PROPERTIES
        ENVIRONMENT
        "FORTSPARSE_UMFPACK_HELPER=$<TARGET_FILE:fortsparse_umfpack_helper>")
endif()

# External context copies keep options but never factors or handles.
fortsparse_add_test(test_fortsparse_external_copy ipc/test_fortsparse_external_copy.f90)
set_tests_properties(test_fortsparse_external_copy PROPERTIES
    LABELS "external;fast" TIMEOUT 1)
if(TARGET fortsparse_umfpack_helper)
    set_tests_properties(test_fortsparse_external_copy PROPERTIES ENVIRONMENT
        "FORTSPARSE_UMFPACK_HELPER=$<TARGET_FILE:fortsparse_umfpack_helper>")
endif()

# Explicit pivot control and a pre-numeric symbolic-memory estimate gate.
if(TARGET fortsparse_umfpack_helper)
    add_executable(test_fortsparse_umfpack_controls
        ipc/test_fortsparse_umfpack_controls.f90)
    target_link_libraries(test_fortsparse_umfpack_controls PRIVATE fortsparse)
    add_executable(test_legacy_umfpack_helper
        ../test-cmake/ipc/test_legacy_umfpack_helper.c)
    target_include_directories(test_legacy_umfpack_helper PRIVATE
        ${PROJECT_SOURCE_DIR}/include)
    target_link_libraries(test_legacy_umfpack_helper PRIVATE
        ${FORTSPARSE_UMFPACK_TARGET} Threads::Threads)
    if(NOT WIN32)
        target_link_libraries(test_legacy_umfpack_helper PRIVATE m)
    endif()
    add_dependencies(test_fortsparse_umfpack_controls
        fortsparse_umfpack_helper test_legacy_umfpack_helper)
    add_test(NAME test_fortsparse_umfpack_controls COMMAND ${CMAKE_COMMAND}
        -DTEST=$<TARGET_FILE:test_fortsparse_umfpack_controls>
        -DHELPER=$<TARGET_FILE:fortsparse_umfpack_helper>
        -DLEGACY=$<TARGET_FILE:test_legacy_umfpack_helper>
        -DRECEIPT=${CMAKE_CURRENT_BINARY_DIR}/legacy-factor-requests.txt
        -P ${PROJECT_SOURCE_DIR}/test-cmake/ipc/run_umfpack_controls.cmake)
endif()

# Exe-directory discovery of the helper, with no environment variable and no
# PATH entry. The test executable is placed in the same directory as the helper
# binary so find_helper resolves it next to the running program. No ENVIRONMENT
# property is set: that is the point of this test.
if(TARGET fortsparse_umfpack_helper AND
   EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/ipc/test_fortsparse_ipc_selflocate.f90")
    add_executable(test_fortsparse_ipc_selflocate
        ipc/test_fortsparse_ipc_selflocate.f90)
    target_link_libraries(test_fortsparse_ipc_selflocate PRIVATE fortsparse)
    set_target_properties(test_fortsparse_ipc_selflocate PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY
        "$<TARGET_FILE_DIR:fortsparse_umfpack_helper>")
    add_test(NAME test_fortsparse_ipc_selflocate
        COMMAND test_fortsparse_ipc_selflocate)
endif()

# MUMPS backend tests and benchmark. The benchmark is built but not registered
# as a test: it takes problem size and thread count on the command line.
if(FORTSPARSE_ENABLE_MUMPS)
    add_library(fs_test_mpi STATIC ../test-cmake/mumps/fs_test_mpi.F90)
    if(MUMPS_parallel)
        target_compile_definitions(fs_test_mpi PRIVATE FORTSPARSE_MUMPS_MPI)
        target_link_libraries(fs_test_mpi PRIVATE MPI::MPI_Fortran)
        fortsparse_add_test(test_fortsparse_mumps_mpi_lifecycle
            ../test-cmake/mumps/test_fortsparse_mumps_mpi_lifecycle.f90)
        target_link_libraries(test_fortsparse_mumps_mpi_lifecycle PRIVATE fs_test_mpi)
        add_executable(test_fortsparse_mumps_mpi_ranks
            ../test-cmake/mumps/test_fortsparse_mumps_mpi_ranks.f90)
        target_link_libraries(test_fortsparse_mumps_mpi_ranks PRIVATE
            fortsparse fs_test_mpi MPI::MPI_Fortran)
        add_test(NAME test_fortsparse_mumps_mpi_ranks
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_fortsparse_mumps_mpi_ranks>)
        set_tests_properties(test_fortsparse_mumps_mpi_ranks PROPERTIES
            PASS_REGULAR_EXPRESSION "PASS" FAIL_REGULAR_EXPRESSION "FAIL"
            TIMEOUT 60 PROCESSORS 2)
    endif()
    add_library(fs_test_matrices STATIC ../test-cmake/mumps/fs_test_matrices.f90)
    target_link_libraries(fs_test_matrices PUBLIC fortsparse fs_test_mpi)
    add_executable(test_fortsparse_mumps ../test-cmake/mumps/test_fortsparse_mumps.f90)
    target_link_libraries(test_fortsparse_mumps PRIVATE fs_test_matrices)
    add_test(NAME test_fortsparse_mumps COMMAND test_fortsparse_mumps)
    add_executable(test_fortsparse_mumps_refactor
        ../test-cmake/mumps/test_fortsparse_mumps_refactor.f90)
    target_link_libraries(test_fortsparse_mumps_refactor PRIVATE
        fs_test_matrices)
    add_test(NAME test_fortsparse_mumps_refactor
        COMMAND test_fortsparse_mumps_refactor)
    if(TARGET MUMPS::ZMUMPS)
        add_executable(test_fortsparse_mumps_complex
            ../test-cmake/mumps/test_fortsparse_mumps_complex.f90)
        target_link_libraries(test_fortsparse_mumps_complex PRIVATE
            fortsparse fs_test_mpi ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
        if(OpenMP_Fortran_FOUND)
            target_link_libraries(test_fortsparse_mumps_complex PRIVATE
                OpenMP::OpenMP_Fortran)
        endif()
        add_test(NAME test_fortsparse_mumps_complex
            COMMAND test_fortsparse_mumps_complex)
        # The default AUTO backend picks SuperLU or MUMPS per matrix.
        if(FORTSPARSE_ENABLE_SUPERLU)
            add_executable(test_fortsparse_auto
                ../test-cmake/mumps/test_fortsparse_auto.f90)
            target_link_libraries(test_fortsparse_auto PRIVATE fs_test_matrices)
            add_test(NAME test_fortsparse_auto COMMAND test_fortsparse_auto)
        endif()
    endif()
    add_executable(bench_fortsparse_mumps
        ../test-cmake/mumps/bench_fortsparse_mumps.f90)
    target_link_libraries(bench_fortsparse_mumps PRIVATE fs_test_matrices)
    # SuperLU vs MUMPS vs AUTO on FE-shaped matrices (README, "Backend
    # selection"); built, not registered.
    if(FORTSPARSE_ENABLE_SUPERLU)
        add_executable(bench_fortsparse_backends
            ../test-cmake/mumps/fs_bench_matrices.f90
            ../test-cmake/mumps/bench_fortsparse_backends.f90)
        target_link_libraries(bench_fortsparse_backends PRIVATE fortsparse)
    endif()
endif()

fortsparse_add_test(test_fortsparse_memory_api solver/test_fortsparse_memory_api.f90)
set_tests_properties(test_fortsparse_memory_api PROPERTIES LABELS "native;fast" TIMEOUT 1)
fortsparse_add_test(test_fortsparse_ldlt_blocked solver/test_fortsparse_ldlt_blocked.f90)
set_tests_properties(test_fortsparse_ldlt_blocked PROPERTIES LABELS "native;fast" TIMEOUT 1)
fortsparse_add_test(test_fortsparse_ldlt_supernodal solver/test_fortsparse_ldlt_supernodal.f90)
set_tests_properties(test_fortsparse_ldlt_supernodal PROPERTIES
    LABELS "native;fast" TIMEOUT 1 PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
fortsparse_add_test(test_fortsparse_ldlt_merge solver/test_fortsparse_ldlt_merge.f90)
set_tests_properties(test_fortsparse_ldlt_merge PROPERTIES
    LABELS "native;performance;fast" TIMEOUT 1 PROCESSORS 2
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
fortsparse_add_test(test_fortsparse_native_ldlt
    solver/test_fortsparse_native_ldlt.f90)
set_tests_properties(test_fortsparse_native_ldlt PROPERTIES
    LABELS "native;fast" TIMEOUT 1)

fortsparse_add_test(test_fortsparse_native_general_api solver/test_fortsparse_native_general_api.f90)
set_tests_properties(test_fortsparse_native_general_api PROPERTIES
    LABELS "native;lu;fast" TIMEOUT 1 WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")

fortsparse_add_test(test_fortsparse_context_copy
    solver/test_fortsparse_context_copy.f90)
set_tests_properties(test_fortsparse_context_copy PROPERTIES
    LABELS "native;lu;fast" TIMEOUT 1)

fortsparse_add_test(test_fortsparse_native_ldlt_z
    solver/test_fortsparse_native_ldlt_z.f90)
set_tests_properties(test_fortsparse_native_ldlt_z PROPERTIES
    LABELS "native;fast" TIMEOUT 1)
fortsparse_add_test(test_fortsparse_inertia solver/test_fortsparse_inertia.f90)
set_tests_properties(test_fortsparse_inertia PROPERTIES LABELS "native;fast" TIMEOUT 1)
# This configure/build check is distinct from ordinary sub-second numerics.
# It disables sparse dependency discovery/fetching and inspects actual linkage.
get_target_property(_precision_fortnum_src fortnum SOURCE_DIR)
get_filename_component(_precision_fortnum_root "${_precision_fortnum_src}" DIRECTORY)
get_target_property(_precision_runtime_src fortgen_runtime SOURCE_DIR)
get_filename_component(_precision_fortgen_root "${_precision_runtime_src}" DIRECTORY)
add_test(NAME test_fortsparse_native_only_build
    COMMAND ${CMAKE_CTEST_COMMAND} --build-and-test
        ${PROJECT_SOURCE_DIR}/test-cmake/native_only
        ${PROJECT_BINARY_DIR}/test-native-only
        --build-generator "${CMAKE_GENERATOR}"
        --build-options
            -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_BUILD_TYPE=Release
            -DFETCHCONTENT_SOURCE_DIR_FORTNUM=${_precision_fortnum_root}
            -DFETCHCONTENT_SOURCE_DIR_FORTGEN_RUNTIME=${_precision_fortgen_root}
        --test-command ${CMAKE_CTEST_COMMAND} -R ^native_only_ --output-on-failure)
set_tests_properties(test_fortsparse_native_only_build PROPERTIES
    LABELS "build;configuration" TIMEOUT 120 PROCESSORS 2
    ENVIRONMENT "CMAKE_BUILD_PARALLEL_LEVEL=2;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1")
if(NOT FORTSPARSE_ENABLE_SUPERLU AND NOT FORTSPARSE_ENABLE_MUMPS)
    fortsparse_add_test(test_fortsparse_native_only ../test-cmake/native_only/smoke.f90)
    set_tests_properties(test_fortsparse_native_only PROPERTIES
        LABELS "native;fast" TIMEOUT 1
        ENVIRONMENT "FORTSPARSE_UMFPACK_HELPER=${CMAKE_CURRENT_BINARY_DIR}/no-helper")
endif()
fortsparse_add_test(test_fortsparse_native_spd_refactor
    solver/test_fortsparse_native_spd_refactor.f90)
set_tests_properties(test_fortsparse_native_spd_refactor PROPERTIES
    LABELS "native;spd;fast" TIMEOUT 1)
fortsparse_add_test(test_fortsparse_native_refinement solver/test_fortsparse_native_refinement.f90)
set_tests_properties(test_fortsparse_native_refinement PROPERTIES
    LABELS "native;refinement;fast" TIMEOUT 1)
fortsparse_add_test(test_fortsparse_componentwise_refinement
    solver/test_fortsparse_componentwise_refinement.f90)
set_tests_properties(test_fortsparse_componentwise_refinement PROPERTIES
    LABELS "native;refinement;lu;fast" TIMEOUT 1)
fortsparse_add_test(test_fortsparse_refined_pivot_default
    solver/test_fortsparse_refined_pivot_default.f90)
set_tests_properties(test_fortsparse_refined_pivot_default PROPERTIES
    LABELS "native;refinement;lu;fast" TIMEOUT 1)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(test_native_refinement_allocations
        ../test-cmake/native/test_native_refinement_allocations.f90
        ../test-cmake/native/refinement_allocation_counter.c)
    target_link_libraries(test_native_refinement_allocations PRIVATE fortsparse)
    target_link_options(test_native_refinement_allocations PRIVATE
        "LINKER:--wrap=malloc" "LINKER:--wrap=calloc" "LINKER:--wrap=realloc")
    add_test(NAME test_native_refinement_allocations COMMAND test_native_refinement_allocations)
    set_tests_properties(test_native_refinement_allocations PROPERTIES LABELS "native;refinement;fast" TIMEOUT 1)
endif()

# Independent dense manufactured RHS tests for retained N/T/H on external factors.
add_executable(test_external_retained_ops ../test-cmake/test_external_retained_ops.F90)
target_link_libraries(test_external_retained_ops PRIVATE fortsparse)
if(FORTSPARSE_ENABLE_MUMPS)
    target_compile_definitions(test_external_retained_ops PRIVATE FS_TEST_MUMPS)
    target_link_libraries(test_external_retained_ops PRIVATE fs_test_mpi)
endif()
foreach(_backend IN ITEMS superlu umfpack mumps)
    if(_backend STREQUAL "superlu")
        set(_id 1)
    elseif(_backend STREQUAL "umfpack")
        set(_id 2)
    else()
        set(_id 3)
    endif()
    add_test(NAME test_external_retained_ops_${_backend}
        COMMAND test_external_retained_ops ${_id})
    set_tests_properties(test_external_retained_ops_${_backend} PROPERTIES
        SKIP_RETURN_CODE 77 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endforeach()
if(TARGET fortsparse_umfpack_helper)
    add_dependencies(test_external_retained_ops fortsparse_umfpack_helper)
    set_tests_properties(test_external_retained_ops_umfpack PROPERTIES
        ENVIRONMENT "FORTSPARSE_UMFPACK_HELPER=$<TARGET_FILE:fortsparse_umfpack_helper>;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    add_executable(test_factor_only_umfpack_helper ../test-cmake/ipc/test_factor_only_umfpack_helper.c)
    target_include_directories(test_factor_only_umfpack_helper PRIVATE ${PROJECT_SOURCE_DIR}/include)
    target_link_libraries(test_factor_only_umfpack_helper PRIVATE ${FORTSPARSE_UMFPACK_TARGET} Threads::Threads)
    if(NOT WIN32)
        target_link_libraries(test_factor_only_umfpack_helper PRIVATE m)
    endif()
    add_test(NAME test_external_retained_ops_stale COMMAND test_external_retained_ops 2 stale)
    set_tests_properties(test_external_retained_ops_stale PROPERTIES
        ENVIRONMENT "FORTSPARSE_UMFPACK_HELPER=$<TARGET_FILE:test_factor_only_umfpack_helper>;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()

fortsparse_add_test(test_fortsparse_native_spd_refinement solver/test_fortsparse_native_spd_refinement.f90)
set_tests_properties(test_fortsparse_native_spd_refinement PROPERTIES
    LABELS "native;refinement;fast" TIMEOUT 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(test_native_spd_refinement_allocations
        ../test-cmake/native/test_native_spd_refinement_allocations.f90
        ../test-cmake/native/refinement_allocation_counter.c)
    target_link_libraries(test_native_spd_refinement_allocations PRIVATE fortsparse)
    target_link_options(test_native_spd_refinement_allocations PRIVATE
        "LINKER:--wrap=malloc" "LINKER:--wrap=calloc" "LINKER:--wrap=realloc")
    add_test(NAME test_native_spd_refinement_allocations COMMAND test_native_spd_refinement_allocations)
    set_tests_properties(test_native_spd_refinement_allocations PROPERTIES
        LABELS "native;refinement;fast" TIMEOUT 1
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()

fortsparse_add_test(test_fortsparse_refined_memory_plan solver/test_fortsparse_refined_memory_plan.f90)
set_tests_properties(test_fortsparse_refined_memory_plan PROPERTIES
    LABELS "native;refinement;fast" TIMEOUT 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(test_general_amd
        ../test-cmake/native/test_general_amd.f90
        ../test-cmake/native/count_general_amd.c)
    target_link_libraries(test_general_amd PRIVATE fortsparse)
    target_link_options(test_general_amd PRIVATE
        "LINKER:--wrap=fsparse_native_order_amd_general")
    if(FORTSPARSE_NATIVE_AMD_FOUND)
        target_compile_definitions(test_general_amd PRIVATE FORTSPARSE_TEST_AMD)
    endif()
    add_test(NAME test_general_amd COMMAND test_general_amd)
    set_tests_properties(test_general_amd PROPERTIES LABELS "native;lu;fast" TIMEOUT 1)
endif()

add_executable(test_native_order_amd_bounded ../test-cmake/native/test_amd_bounded.c)
target_link_libraries(test_native_order_amd_bounded PRIVATE fortsparse ${CMAKE_DL_LIBS})
target_include_directories(test_native_order_amd_bounded PRIVATE ../src/native)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND FORTSPARSE_NATIVE_BOUNDED_AMD_FOUND
        AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_definitions(test_native_order_amd_bounded PRIVATE FS_TEST_BOUNDED_AMD)
    target_link_libraries(test_native_order_amd_bounded PRIVATE ${FORTSPARSE_NATIVE_AMD_TARGET})
    target_link_options(test_native_order_amd_bounded PRIVATE "LINKER:--wrap=amd_2" "LINKER:--wrap=amd_version")
endif()
add_test(NAME test_native_order_amd_bounded COMMAND test_native_order_amd_bounded)
set_tests_properties(test_native_order_amd_bounded PROPERTIES LABELS "native;lu;fast" TIMEOUT 1)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(test_general_amd_bounded
        ../test-cmake/native/test_general_amd_bounded.f90
        ../test-cmake/native/count_bounded_amd.c)
    target_link_libraries(test_general_amd_bounded PRIVATE fortsparse)
    target_include_directories(test_general_amd_bounded PRIVATE ../src/native)
    target_link_options(test_general_amd_bounded PRIVATE "LINKER:--wrap=fsparse_native_amd_bounded")
    add_test(NAME test_general_amd_bounded COMMAND test_general_amd_bounded)
    set_tests_properties(test_general_amd_bounded PROPERTIES LABELS "native;lu;refinement;fast" TIMEOUT 1
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(test_spd_panel_no_blas
        ../test-cmake/native/test_spd_panel_no_blas.f90
        ../test-cmake/native/count_spd_dtrsm.c)
    target_link_libraries(test_spd_panel_no_blas PRIVATE fortsparse)
    target_link_options(test_spd_panel_no_blas PRIVATE "LINKER:--wrap=dtrsm_")
    add_test(NAME test_spd_panel_no_blas COMMAND test_spd_panel_no_blas)
    set_tests_properties(test_spd_panel_no_blas PROPERTIES LABELS "native;memory;fast"
        TIMEOUT 1 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()

# Exercise the public ABI from actual C and C++ consumers, including C99 complex.
enable_language(CXX)
add_executable(test_native_c_api ../test-cmake/native/test_native_c_api.c)
target_link_libraries(test_native_c_api PRIVATE fortsparse)
add_test(NAME test_native_c_api COMMAND test_native_c_api)
add_executable(test_native_cxx_api ../test-cmake/native/test_native_cxx_api.cpp)
target_compile_features(test_native_cxx_api PRIVATE cxx_std_11)
target_link_libraries(test_native_cxx_api PRIVATE fortsparse)
add_test(NAME test_native_cxx_api COMMAND test_native_cxx_api)
if(FORTSPARSE_NATIVE_AMD_FOUND)
    target_compile_definitions(test_native_c_api PRIVATE FS_TEST_EXPECT_AMD)
    target_compile_definitions(test_native_cxx_api PRIVATE FS_TEST_EXPECT_AMD)
endif()
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    foreach(_c_api_target test_native_c_api test_native_cxx_api)
        target_sources(${_c_api_target} PRIVATE ../test-cmake/native/count_general_amd.c)
        target_compile_definitions(${_c_api_target} PRIVATE FS_TEST_COUNT_AMD)
        target_link_options(${_c_api_target} PRIVATE "LINKER:--wrap=fsparse_native_order_amd_general")
    endforeach()
endif()
set_tests_properties(test_native_c_api test_native_cxx_api PROPERTIES
    LABELS "native;c-api;fast" TIMEOUT 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(test_spd_amd_bounded
        ../test-cmake/native/test_spd_amd_bounded.f90
        ../test-cmake/native/probe_spd_amd_allocations.c)
    target_link_libraries(test_spd_amd_bounded PRIVATE fortsparse)
    target_include_directories(test_spd_amd_bounded PRIVATE ../src/native)
    target_link_options(test_spd_amd_bounded PRIVATE
        "LINKER:--wrap=fsparse_native_amd_bounded" "LINKER:--wrap=fsparse_native_order_amd"
        "LINKER:--wrap=malloc" "LINKER:--wrap=calloc" "LINKER:--wrap=realloc" "LINKER:--wrap=free")
    add_test(NAME test_spd_amd_bounded COMMAND test_spd_amd_bounded)
    set_tests_properties(test_spd_amd_bounded PROPERTIES LABELS "native;memory;refinement;fast"
        TIMEOUT 1 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()

# Public synthetic NEO-2 spline system: preserve the original observable gate,
# not just a tiny normwise backward error in its ill-conditioned coefficient basis.
add_executable(test_native_lu_spline ../test-cmake/native/test_native_lu_spline.f90)
target_link_libraries(test_native_lu_spline PRIVATE fortsparse)
add_test(NAME test_native_lu_spline COMMAND test_native_lu_spline
    "${CMAKE_CURRENT_SOURCE_DIR}/../test-cmake/native/data")
set_tests_properties(test_native_lu_spline PROPERTIES
    LABELS "native;lu;fast" TIMEOUT 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
fortsparse_add_test(test_native_lu_singletons ../test-cmake/native/test_native_lu_singletons.f90)
set_tests_properties(test_native_lu_singletons PROPERTIES
    LABELS "native;lu;fast" TIMEOUT 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

fortsparse_add_test(test_fortsparse_native_raw_ldlt solver/test_fortsparse_native_raw_ldlt.f90)
set_tests_properties(test_fortsparse_native_raw_ldlt PROPERTIES
    LABELS "native;fast" TIMEOUT 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(test_ldlt_amd ../test-cmake/native/test_ldlt_amd.f90
        ../test-cmake/native/probe_spd_amd_allocations.c)
    target_link_libraries(test_ldlt_amd PRIVATE fortsparse)
    target_include_directories(test_ldlt_amd PRIVATE ../src/native)
    target_link_options(test_ldlt_amd PRIVATE
        "LINKER:--wrap=fsparse_native_amd_bounded" "LINKER:--wrap=fsparse_native_order_amd"
        "LINKER:--wrap=malloc" "LINKER:--wrap=calloc" "LINKER:--wrap=realloc" "LINKER:--wrap=free")
    add_test(NAME test_ldlt_amd COMMAND test_ldlt_amd)
    set_tests_properties(test_ldlt_amd PROPERTIES LABELS "native;memory;refinement;fast"
        TIMEOUT 1 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(test_csc_normalize_symmetric ../test-cmake/native/test_csc_normalize_symmetric.f90
        ../test-cmake/native/refinement_allocation_counter.c)
    target_link_libraries(test_csc_normalize_symmetric PRIVATE fortsparse)
    target_link_options(test_csc_normalize_symmetric PRIVATE
        "LINKER:--wrap=malloc" "LINKER:--wrap=calloc" "LINKER:--wrap=realloc")
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        target_compile_options(test_csc_normalize_symmetric PRIVATE
            $<$<COMPILE_LANGUAGE:Fortran>:-ffpe-trap=invalid,zero,overflow>)
    endif()
    add_test(NAME test_csc_normalize_symmetric COMMAND test_csc_normalize_symmetric)
    set_tests_properties(test_csc_normalize_symmetric PROPERTIES LABELS "native;csc;fast"
        TIMEOUT 1 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()

fortsparse_add_test(test_lu_row_preference ../test-cmake/native/test_lu_row_preference.f90)
set_tests_properties(test_lu_row_preference PROPERTIES LABELS "native;lu;fast"
    TIMEOUT 1 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

fortsparse_add_test(test_fortsparse_structural_matching solver/test_fortsparse_structural_matching.f90)
fortsparse_add_test(test_fortsparse_matching_api solver/test_fortsparse_matching_api.f90)
set_tests_properties(test_fortsparse_structural_matching test_fortsparse_matching_api
    PROPERTIES LABELS "native;lu;memory;fast" TIMEOUT 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
fortsparse_add_test(test_fortsparse_mpm solver/test_fortsparse_mpm.f90)
set_tests_properties(test_fortsparse_mpm PROPERTIES LABELS "native;lu;fast" TIMEOUT 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(test_mpm_allocations
        ../test-cmake/native/test_mpm_allocations.f90
        ../test-cmake/native/probe_spd_amd_allocations.c)
    target_include_directories(test_mpm_allocations PRIVATE ../src/native)
    target_link_libraries(test_mpm_allocations PRIVATE fortsparse)
    target_link_options(test_mpm_allocations PRIVATE
        "LINKER:--wrap=malloc" "LINKER:--wrap=calloc" "LINKER:--wrap=realloc" "LINKER:--wrap=free"
        "LINKER:--wrap=fsparse_native_amd_bounded" "LINKER:--wrap=fsparse_native_order_amd")
    add_test(NAME test_mpm_allocations COMMAND test_mpm_allocations)
    set_tests_properties(test_mpm_allocations PROPERTIES LABELS "native;lu;memory;fast" TIMEOUT 1
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(test_slu_allocations
        ../test-cmake/native/test_slu_allocations.f90
        ../test-cmake/native/probe_spd_amd_allocations.c)
    target_include_directories(test_slu_allocations PRIVATE ../src/native)
    target_link_libraries(test_slu_allocations PRIVATE fortsparse)
    target_link_options(test_slu_allocations PRIVATE
        "LINKER:--wrap=malloc" "LINKER:--wrap=calloc" "LINKER:--wrap=realloc" "LINKER:--wrap=free"
        "LINKER:--wrap=fsparse_native_amd_bounded" "LINKER:--wrap=fsparse_native_order_amd")
    add_test(NAME test_slu_allocations COMMAND test_slu_allocations)
    set_tests_properties(test_slu_allocations PROPERTIES LABELS "native;lu;memory;fast" TIMEOUT 1
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(test_matching_allocations
        ../test-cmake/native/test_matching_allocations.f90
        ../test-cmake/native/probe_spd_amd_allocations.c)
    target_include_directories(test_matching_allocations PRIVATE ../src/native)
    target_link_libraries(test_matching_allocations PRIVATE fortsparse)
    target_link_options(test_matching_allocations PRIVATE
        "LINKER:--wrap=malloc" "LINKER:--wrap=calloc" "LINKER:--wrap=realloc" "LINKER:--wrap=free"
        "LINKER:--wrap=fsparse_native_amd_bounded" "LINKER:--wrap=fsparse_native_order_amd")
    add_test(NAME test_matching_allocations COMMAND test_matching_allocations)
    set_tests_properties(test_matching_allocations PROPERTIES LABELS "native;lu;memory;fast" TIMEOUT 1
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(test_native_c_api_memory ../test-cmake/native/test_native_c_api_memory.c
        ../test-cmake/native/probe_spd_amd_allocations.c)
    target_link_libraries(test_native_c_api_memory PRIVATE fortsparse)
    target_include_directories(test_native_c_api_memory PRIVATE ../src/native)
    target_link_options(test_native_c_api_memory PRIVATE
        "LINKER:--wrap=fsparse_native_amd_bounded" "LINKER:--wrap=fsparse_native_order_amd"
        "LINKER:--wrap=malloc" "LINKER:--wrap=calloc" "LINKER:--wrap=realloc" "LINKER:--wrap=free")
    add_test(NAME test_native_c_api_memory COMMAND test_native_c_api_memory)
    set_tests_properties(test_native_c_api_memory PROPERTIES LABELS "native;c-api;memory;fast"
        TIMEOUT 1 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(test_structured_raw_allocations
        ../test-cmake/native/test_structured_raw_allocations.f90
        ../test-cmake/native/probe_spd_amd_allocations.c)
    target_include_directories(test_structured_raw_allocations PRIVATE ../src/native)
    target_link_libraries(test_structured_raw_allocations PRIVATE fortsparse)
    target_link_options(test_structured_raw_allocations PRIVATE
        "LINKER:--wrap=malloc" "LINKER:--wrap=calloc" "LINKER:--wrap=realloc" "LINKER:--wrap=free"
        "LINKER:--wrap=fsparse_native_amd_bounded" "LINKER:--wrap=fsparse_native_order_amd")
    add_test(NAME test_structured_raw_allocations COMMAND test_structured_raw_allocations)
    set_tests_properties(test_structured_raw_allocations PROPERTIES LABELS "native;memory;fast"
        TIMEOUT 1 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()

fortsparse_add_test(test_fortsparse_structured_raw solver/test_fortsparse_structured_raw.f90)
set_tests_properties(test_fortsparse_structured_raw PROPERTIES LABELS "native;memory;fast"
    TIMEOUT 1 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

if(FORTSPARSE_NATIVE_METIS_FOUND)
    foreach(_c_api_nd_target test_native_c_api test_native_cxx_api test_native_c_api_memory)
        if(TARGET ${_c_api_nd_target})
            target_compile_definitions(${_c_api_nd_target} PRIVATE FS_TEST_EXPECT_METIS)
        endif()
    endforeach()
endif()
if(TARGET test_native_c_api_memory)
    target_link_options(test_native_c_api_memory PRIVATE "LINKER:--wrap=fsparse_native_order_general")
endif()

add_executable(test_native_lu_compensated ../test-cmake/native/test_native_lu_compensated.F90)
target_link_libraries(test_native_lu_compensated PRIVATE fortsparse)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_sources(test_native_lu_compensated PRIVATE ../test-cmake/native/refinement_allocation_counter.c)
    target_compile_definitions(test_native_lu_compensated PRIVATE FS_HAVE_ALLOC_PROBE)
    target_link_options(test_native_lu_compensated PRIVATE
        "LINKER:--wrap=malloc" "LINKER:--wrap=calloc" "LINKER:--wrap=realloc")
endif()
add_test(NAME test_native_lu_compensated COMMAND test_native_lu_compensated)
set_tests_properties(test_native_lu_compensated PROPERTIES
    LABELS "native;lu;memory;fast" TIMEOUT 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

# Experimental low-level real LU target panels; public defaults stay scalar.
fortsparse_add_test(test_fortsparse_lu_panel solver/test_fortsparse_lu_panel.f90)
set_tests_properties(test_fortsparse_lu_panel PROPERTIES LABELS "native;lu;fast" TIMEOUT 1)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(test_lu_panel_allocations
        ../test-cmake/native/test_lu_panel_allocations.f90
        ../test-cmake/native/lu_panel_alloc_probe.c)
    target_link_libraries(test_lu_panel_allocations PRIVATE fortsparse)
    target_link_options(test_lu_panel_allocations PRIVATE "LINKER:--wrap=malloc"
        "LINKER:--wrap=calloc" "LINKER:--wrap=realloc" "LINKER:--wrap=free")
    add_test(NAME test_lu_panel_allocations COMMAND test_lu_panel_allocations)
    set_tests_properties(test_lu_panel_allocations PROPERTIES
        LABELS "native;lu;memory;fast" TIMEOUT 1)
endif()

fortsparse_add_test(test_fortsparse_precision solver/test_fortsparse_precision.f90)

fortsparse_add_test(test_fortsparse_mixed_lu solver/test_fortsparse_mixed_lu.f90)
set_tests_properties(test_fortsparse_mixed_lu PROPERTIES
    LABELS "native;precision;fast" TIMEOUT 1 PROCESSORS 2
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
fortsparse_add_test(test_fortsparse_stability solver/test_fortsparse_stability.f90)
target_sources(test_fortsparse_stability PRIVATE solver/stability_oracle.f90)
set_tests_properties(test_fortsparse_stability PROPERTIES
    LABELS "native;stability;fast" TIMEOUT 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(test_mixed_allocations ../test-cmake/native/test_mixed_allocations.f90
        ../test-cmake/native/ilutp_alloc_probe.c)
    target_link_libraries(test_mixed_allocations PRIVATE fortsparse)
    target_link_options(test_mixed_allocations PRIVATE "LINKER:--wrap=malloc"
        "LINKER:--wrap=calloc" "LINKER:--wrap=realloc")
    add_test(NAME test_mixed_allocations COMMAND test_mixed_allocations)
    set_tests_properties(test_mixed_allocations PROPERTIES LABELS "native;precision;memory;fast"
        TIMEOUT 1 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()

fortsparse_add_test(test_fortsparse_amg solver/test_fortsparse_amg.f90)
fortsparse_add_test(test_fortsparse_csr_smoothers solver/test_fortsparse_csr_smoothers.f90)
foreach(amg_test test_fortsparse_amg test_fortsparse_csr_smoothers)
    target_link_libraries(${amg_test} PRIVATE fortsparse_amg_test_support)
    set_tests_properties(${amg_test} PROPERTIES LABELS "native;preconditioner;fast"
        TIMEOUT 1 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
endforeach()

if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_executable(test_amg_allocations
        ../test-cmake/native/test_amg_allocations.f90
        ../test-cmake/native/ilutp_alloc_probe.c)
    target_link_libraries(test_amg_allocations PRIVATE fortsparse_amg_test_support)
    target_link_options(test_amg_allocations PRIVATE "LINKER:--wrap=malloc"
        "LINKER:--wrap=calloc" "LINKER:--wrap=realloc")
    add_test(NAME test_amg_allocations COMMAND test_amg_allocations)
    set_tests_properties(test_amg_allocations PROPERTIES
        LABELS "native;preconditioner;memory;fast" TIMEOUT 1
        ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
endif()
