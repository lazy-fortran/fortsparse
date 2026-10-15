# fortsparse test suite. Each test is registered only when its source file is
# present, so the library configures and builds before the test programs are
# written and every test activates automatically once its source lands.

function(fortsparse_add_test name source)
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
        return()
    endif()
    add_executable(${name} ${source})
    target_link_libraries(${name} PRIVATE fortsparse)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

fortsparse_add_test(test_fortsparse_smoke test_fortsparse_smoke.f90)
fortsparse_add_test(test_fortsparse_kinds test_fortsparse_kinds.f90)
fortsparse_add_test(test_fortsparse_status test_fortsparse_status.f90)
fortsparse_add_test(test_fortsparse_csc csc/test_fortsparse_csc.f90)
fortsparse_add_test(test_fortsparse_matvec_omp csc/test_fortsparse_matvec_omp.f90)

# Strong-scaling benchmark of csc_matvec; built, not registered as a test.
add_executable(bench_fortsparse_matvec ../bench/bench_csc_matvec.f90)
target_link_libraries(bench_fortsparse_matvec PRIVATE fortsparse)

if(UNIX)
    add_executable(test_fortsparse_ipc_wait ipc/test_fortsparse_ipc_wait.c)
    target_include_directories(test_fortsparse_ipc_wait PRIVATE
        ${PROJECT_SOURCE_DIR}/include)
    target_link_libraries(test_fortsparse_ipc_wait PRIVATE fortsparse)
    add_test(NAME test_fortsparse_ipc_wait COMMAND test_fortsparse_ipc_wait)
    set_tests_properties(test_fortsparse_ipc_wait PROPERTIES TIMEOUT 10)
endif()

# These cases factor through the default SUPERLU backend, so they only make
# sense when the SuperLU backend is compiled in. The IPC, selflocate, kinds,
# status, and csc tests do not depend on it and stay registered always.
if(FORTSPARSE_ENABLE_SUPERLU)
    fortsparse_add_test(test_fortsparse_solve_real
        solver/test_fortsparse_solve_real.f90)
    fortsparse_add_test(test_fortsparse_solve_derivatives
        solver/test_fortsparse_solve_derivatives.f90)
    fortsparse_add_test(test_fortsparse_solve_complex
        solver/test_fortsparse_solve_complex.f90)
    fortsparse_add_test(test_fortsparse_poisson
        solver/test_fortsparse_poisson.f90)
    fortsparse_add_test(test_fortsparse_errors solver/test_fortsparse_errors.f90)
    # Forces a real SuperLU allocation failure under an address-space limit.
    fortsparse_add_test(test_fortsparse_superlu_status
        solver/test_fortsparse_superlu_status.f90)
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
    add_library(fs_test_mpi STATIC mumps/fs_test_mpi.F90)
    if(MUMPS_parallel)
        target_compile_definitions(fs_test_mpi PRIVATE FORTSPARSE_MUMPS_MPI)
        target_link_libraries(fs_test_mpi PRIVATE MPI::MPI_Fortran)
        fortsparse_add_test(test_fortsparse_mumps_mpi_lifecycle
            mumps/test_fortsparse_mumps_mpi_lifecycle.f90)
        target_link_libraries(test_fortsparse_mumps_mpi_lifecycle PRIVATE fs_test_mpi)
        add_executable(test_fortsparse_mumps_mpi_ranks
            mumps/test_fortsparse_mumps_mpi_ranks.f90)
        target_link_libraries(test_fortsparse_mumps_mpi_ranks PRIVATE
            fortsparse fs_test_mpi MPI::MPI_Fortran)
        add_test(NAME test_fortsparse_mumps_mpi_ranks
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_fortsparse_mumps_mpi_ranks>)
        set_tests_properties(test_fortsparse_mumps_mpi_ranks PROPERTIES
            PASS_REGULAR_EXPRESSION "PASS" FAIL_REGULAR_EXPRESSION "FAIL"
            TIMEOUT 60 PROCESSORS 2)
    endif()
    add_library(fs_test_matrices STATIC mumps/fs_test_matrices.f90)
    target_link_libraries(fs_test_matrices PUBLIC fortsparse fs_test_mpi)
    add_executable(test_fortsparse_mumps mumps/test_fortsparse_mumps.f90)
    target_link_libraries(test_fortsparse_mumps PRIVATE fs_test_matrices)
    add_test(NAME test_fortsparse_mumps COMMAND test_fortsparse_mumps)
    add_executable(test_fortsparse_mumps_refactor
        mumps/test_fortsparse_mumps_refactor.f90)
    target_link_libraries(test_fortsparse_mumps_refactor PRIVATE
        fs_test_matrices)
    add_test(NAME test_fortsparse_mumps_refactor
        COMMAND test_fortsparse_mumps_refactor)
    if(TARGET MUMPS::ZMUMPS)
        add_executable(test_fortsparse_mumps_complex
            mumps/test_fortsparse_mumps_complex.f90)
        target_link_libraries(test_fortsparse_mumps_complex PRIVATE
            fortsparse fs_test_mpi ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
        if(OpenMP_Fortran_FOUND)
            target_link_libraries(test_fortsparse_mumps_complex PRIVATE
                OpenMP::OpenMP_Fortran)
        endif()
        add_test(NAME test_fortsparse_mumps_complex
            COMMAND test_fortsparse_mumps_complex)
    endif()
    add_executable(bench_fortsparse_mumps mumps/bench_fortsparse_mumps.f90)
    target_link_libraries(bench_fortsparse_mumps PRIVATE fs_test_matrices)
endif()
