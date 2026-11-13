include(CheckFortranSourceRuns)

file(READ "${CMAKE_CURRENT_LIST_DIR}/../src/fortsparse_compiler.f90" _compiler_probe)
file(READ "${CMAKE_CURRENT_LIST_DIR}/../test/solver/test_fortsparse_compiler_self_assignment.f90"
    _compiler_probe_driver)
check_fortran_source_runs("${_compiler_probe}\n${_compiler_probe_driver}"
    FORTSPARSE_COMPILER_SELF_ASSIGNMENT_SAFE SRC_EXT F90)
if(NOT FORTSPARSE_COMPILER_SELF_ASSIGNMENT_SAFE)
    message(FATAL_ERROR
        "The Fortran compiler failed the solver self-assignment lifetime probe. "
        "Use a compiler that passes test_fortsparse_compiler_self_assignment. "
        "Cross builds need a working CMAKE_CROSSCOMPILING_EMULATOR.")
endif()
