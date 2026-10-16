foreach(_case real complex hrot_complex)
    add_test(NAME freefem_native_${_case}
        COMMAND "${FREEFEM_EXECUTABLE}" -nw -nowait -v 0 "${CMAKE_CURRENT_SOURCE_DIR}/${_case}.edp")
endforeach()
foreach(_case real complex)
    add_test(NAME freefem_native_singular_${_case} COMMAND "${CMAKE_COMMAND}"
        "-DFREEFEM=${FREEFEM_EXECUTABLE}" "-DSCRIPT=${CMAKE_CURRENT_SOURCE_DIR}/singular_${_case}.edp"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/expect_failure.cmake")
endforeach()
set_tests_properties(freefem_native_real freefem_native_complex freefem_native_hrot_complex
    freefem_native_singular_real freefem_native_singular_complex PROPERTIES
    LABELS "native;freefem;fast" TIMEOUT 1
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    ENVIRONMENT "FF_LOADPATH=$<TARGET_FILE_DIR:fortsparse_freefem>;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    find_program(FORTSPARSE_LDD NAMES ldd REQUIRED)
    find_program(FORTSPARSE_NM NAMES nm REQUIRED)
    add_test(NAME freefem_native_module_linkage COMMAND "${CMAKE_COMMAND}"
        "-DLIBRARY=$<TARGET_FILE:fortsparse_freefem>" "-DLDD=${FORTSPARSE_LDD}"
        "-DNM=${FORTSPARSE_NM}" -P "${CMAKE_CURRENT_SOURCE_DIR}/module_linkage.cmake")
    set_tests_properties(freefem_native_module_linkage PROPERTIES LABELS "native;freefem;fast" TIMEOUT 1)
endif()
