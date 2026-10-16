execute_process(COMMAND "${FREEFEM}" -nw -nowait -v 0 "${SCRIPT}"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 0.8)
if(result EQUAL 0 OR NOT "${output}\n${error}" MATCHES "FortSparse native FreeFem failure 1")
    message(FATAL_ERROR "Expected singular native failure without fallback: ${result}\n${output}\n${error}")
endif()
message(STATUS "PASS FreeFem native singular failure without fallback")
