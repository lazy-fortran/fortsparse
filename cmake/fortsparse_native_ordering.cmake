include_guard(GLOBAL)
option(FORTSPARSE_ENABLE_METIS "Use optional METIS native fill-reducing ordering" ON)

function(fortsparse_configure_native_ordering target)
    set(FORTSPARSE_NATIVE_METIS_FOUND OFF PARENT_SCOPE)
    target_sources(${target} PRIVATE
        "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../src/native/fsparse_native_order.c")
    if(NOT FORTSPARSE_ENABLE_METIS)
        return()
    endif()

    # Reuse a package target, including one created by an enabled MUMPS build.
    set(_metis_target)
    foreach(_candidate fortsparse_metis METIS::METIS METIS::metis metis METIS)
        if(TARGET ${_candidate})
            set(_metis_target ${_candidate})
            break()
        endif()
    endforeach()
    if(NOT _metis_target)
        find_path(FORTSPARSE_METIS_INCLUDE_DIR NAMES metis.h
            HINTS ${METIS_ROOT} ENV METIS_ROOT PATH_SUFFIXES include)
        find_library(FORTSPARSE_METIS_LIBRARY NAMES metis
            HINTS ${METIS_ROOT} ENV METIS_ROOT PATH_SUFFIXES lib lib64)
        mark_as_advanced(FORTSPARSE_METIS_INCLUDE_DIR FORTSPARSE_METIS_LIBRARY)
        if(FORTSPARSE_METIS_INCLUDE_DIR AND FORTSPARSE_METIS_LIBRARY)
            add_library(fortsparse_metis UNKNOWN IMPORTED GLOBAL)
            set_target_properties(fortsparse_metis PROPERTIES
                IMPORTED_LOCATION "${FORTSPARSE_METIS_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${FORTSPARSE_METIS_INCLUDE_DIR}")
            set(_metis_target fortsparse_metis)
        endif()
    endif()
    if(_metis_target)
        target_link_libraries(${target} PRIVATE ${_metis_target})
        target_compile_definitions(${target} PRIVATE FORTSPARSE_HAVE_METIS)
        set(FORTSPARSE_NATIVE_METIS_FOUND ON PARENT_SCOPE)
    else()
        message(STATUS "METIS unavailable: native nested dissection disabled")
        set(FORTSPARSE_NATIVE_METIS_FOUND OFF PARENT_SCOPE)
    endif()
endfunction()
