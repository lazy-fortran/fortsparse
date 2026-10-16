# Query the built program's actual loader resolution. Static sparse targets
# are independently excluded during this fixture's configure. No solver runs.
find_program(loader_inspector NAMES ldd REQUIRED)
execute_process(COMMAND "${loader_inspector}" "${EXECUTABLE}"
    RESULT_VARIABLE result OUTPUT_VARIABLE resolved ERROR_VARIABLE diagnostic)
if(NOT result EQUAL 0)
    if(resolved MATCHES "not a dynamic executable|statically linked" OR
       diagnostic MATCHES "not a dynamic executable|statically linked")
        return()
    endif()
    message(FATAL_ERROR "Native-only dependency inspection failed: ${diagnostic}")
endif()
if(resolved MATCHES "not found")
    message(FATAL_ERROR "Native-only unresolved runtime dependency: ${resolved}")
endif()
string(REGEX MATCHALL "[A-Za-z0-9_+.-]+[.]so([.][0-9]+)*" libraries "${resolved}")
foreach(name IN LISTS libraries)
    string(TOLOWER "${name}" name)
    if(name MATCHES "^(lib)?(superlu|umfpack|[sdcz]?mumps|cholmod|klu|spqr|btf|camd|ccolamd|colamd|amd|metis|suitesparseconfig|pastix|pardiso|cudss|fortsym|fortgen|symengine|yacas|flint)([._-]|$)")
        message(FATAL_ERROR "Native-only execution links forbidden dependency ${name}")
    endif()
endforeach()
message(STATUS "Native-only runtime resolution contains no sparse solver or CAS")
