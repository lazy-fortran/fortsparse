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
fortsparse_add_test(test_fortsparse_solve_real solver/test_fortsparse_solve_real.f90)
fortsparse_add_test(test_fortsparse_solve_complex
    solver/test_fortsparse_solve_complex.f90)
fortsparse_add_test(test_fortsparse_poisson solver/test_fortsparse_poisson.f90)
fortsparse_add_test(test_fortsparse_errors solver/test_fortsparse_errors.f90)
# The no-helper error path must run where no helper sits next to it, otherwise
# exe-directory discovery would resolve one and the backend would be available.
# Keep this executable out of the bin dir that holds fortsparse_umfpack_helper.
if(TARGET test_fortsparse_errors)
    set_target_properties(test_fortsparse_errors PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/no_helper")
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
