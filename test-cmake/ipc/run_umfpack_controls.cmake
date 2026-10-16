# One focused behavioral check, with the actual helper then a stale protocol.
execute_process(COMMAND "${CMAKE_COMMAND}" -E env
    "FORTSPARSE_UMFPACK_HELPER=${HELPER}" "${TEST}"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "UMFPACK control/budget regression failed: ${result}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E env
    "FORTSPARSE_UMFPACK_HELPER=${LEGACY}"
    "FORTSPARSE_LEGACY_RECEIPT=${RECEIPT}" "${TEST}" legacy
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Legacy-helper rejection failed: ${result}")
endif()
