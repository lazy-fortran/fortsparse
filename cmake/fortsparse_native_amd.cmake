include_guard(GLOBAL)
option(FORTSPARSE_ENABLE_AMD "Use optional BSD-3-Clause AMD native ordering" ON)

function(fortsparse_configure_native_amd target)
    set(FORTSPARSE_NATIVE_AMD_FOUND OFF PARENT_SCOPE)
    set(FORTSPARSE_NATIVE_BOUNDED_AMD_FOUND OFF PARENT_SCOPE)
    target_sources(${target} PRIVATE
        "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../src/native/fsparse_native_order_amd.c"
        "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../src/native/fsparse_native_amd_bounded.c"
        "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../src/native/fsparse_native_amd_core.c")
    if(NOT FORTSPARSE_ENABLE_AMD)
        return()
    endif()
    if(NOT TARGET SuiteSparse::AMD AND NOT TARGET SuiteSparse::AMD_static
       AND NOT TARGET AMD::AMD AND NOT TARGET amd)
        find_package(AMD CONFIG QUIET)
    endif()
    set(_amd_target)
    foreach(_candidate fortsparse_amd SuiteSparse::AMD SuiteSparse::AMD_static AMD::AMD amd)
        if(TARGET ${_candidate})
            set(_amd_target ${_candidate})
            break()
        endif()
    endforeach()
    if(NOT _amd_target)
        find_path(FORTSPARSE_AMD_INCLUDE_DIR NAMES amd.h
            HINTS ${AMD_ROOT} ENV AMD_ROOT PATH_SUFFIXES include include/suitesparse suitesparse)
        find_library(FORTSPARSE_AMD_LIBRARY NAMES amd
            HINTS ${AMD_ROOT} ENV AMD_ROOT PATH_SUFFIXES lib lib64)
        mark_as_advanced(FORTSPARSE_AMD_INCLUDE_DIR FORTSPARSE_AMD_LIBRARY)
        if(FORTSPARSE_AMD_INCLUDE_DIR AND FORTSPARSE_AMD_LIBRARY)
            if(FORTSPARSE_AMD_LIBRARY MATCHES "[.]a$")
                find_library(FORTSPARSE_AMD_CONFIG_LIBRARY NAMES suitesparseconfig
                    HINTS ${AMD_ROOT} ENV AMD_ROOT PATH_SUFFIXES lib lib64)
                mark_as_advanced(FORTSPARSE_AMD_CONFIG_LIBRARY)
                if(NOT FORTSPARSE_AMD_CONFIG_LIBRARY)
                    message(STATUS "Static AMD needs its SuiteSparse_config dependency")
                    return()
                endif()
            endif()
            add_library(fortsparse_amd UNKNOWN IMPORTED GLOBAL)
            set_target_properties(fortsparse_amd PROPERTIES
                IMPORTED_LOCATION "${FORTSPARSE_AMD_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${FORTSPARSE_AMD_INCLUDE_DIR}")
            if(FORTSPARSE_AMD_LIBRARY MATCHES "[.]a$")
                set_property(TARGET fortsparse_amd PROPERTY INTERFACE_LINK_LIBRARIES
                    "${FORTSPARSE_AMD_CONFIG_LIBRARY}")
                if(NOT WIN32)
                    set_property(TARGET fortsparse_amd APPEND PROPERTY
                        INTERFACE_LINK_LIBRARIES m)
                endif()
            endif()
            set(_amd_target fortsparse_amd)
        endif()
    endif()
    if(_amd_target)
        target_link_libraries(${target} PRIVATE ${_amd_target})
        target_compile_definitions(${target} PRIVATE FORTSPARSE_HAVE_AMD)
        set(FORTSPARSE_NATIVE_AMD_FOUND ON PARENT_SCOPE)
        set(FORTSPARSE_NATIVE_AMD_TARGET ${_amd_target} PARENT_SCOPE)
        # amd_2 is public but deliberately unchecked/destructive. Only the
        # reviewed upstream 3.3.4 allocation closure enables the bounded ABI.
        include(CheckCSourceCompiles)
        set(CMAKE_REQUIRED_LIBRARIES ${_amd_target})
        unset(FORTSPARSE_AMD_BOUNDED_ABI CACHE)
        check_c_source_compiles([=[
#include <stdint.h>
#include <amd.h>
#if AMD_MAIN_VERSION != 3 || AMD_SUB_VERSION != 3 || AMD_SUBSUB_VERSION != 4
#error Unreviewed AMD release
#endif
_Static_assert(sizeof(int)==sizeof(int32_t), "32-bit C int required");
typedef void (*amd2_t)(int32_t,int32_t*,int32_t*,int32_t*,int32_t,int32_t,
    int32_t*,int32_t*,int32_t*,int32_t*,int32_t*,int32_t*,int32_t*,double*,double*);
_Static_assert(_Generic(&amd_2, amd2_t:1, default:0), "amd_2 signature mismatch");
int main(void) {
  int version[3]; int32_t w[1]={0};
  amd_version(version);
  amd_2(1,w,w,w,1,0,w,w,w,w,w,w,w,0,0);
  return version[0];
}]=] FORTSPARSE_AMD_BOUNDED_ABI)
        if(FORTSPARSE_AMD_BOUNDED_ABI)
            target_compile_definitions(${target} PRIVATE FORTSPARSE_HAVE_BOUNDED_AMD)
            set(FORTSPARSE_NATIVE_BOUNDED_AMD_FOUND ON PARENT_SCOPE)
        else()
            message(STATUS "Audited AMD caller-workspace ABI unavailable; project-owned bounded AMD used")
        endif()
    else()
        message(STATUS "SuiteSparse AMD unavailable: project-owned AMD ordering used")
    endif()
endfunction()
