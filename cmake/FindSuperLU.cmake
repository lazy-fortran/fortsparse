# Locate sequential SuperLU (BSD). Provides the imported target
# SuperLU::SuperLU with its include directory and library. Headers usually live
# under a "superlu" suffix (e.g. /usr/include/superlu/slu_ddefs.h).

find_path(SuperLU_INCLUDE_DIR
    NAMES slu_ddefs.h
    PATH_SUFFIXES superlu SuperLU
)

find_library(SuperLU_LIBRARY
    NAMES superlu
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SuperLU
    REQUIRED_VARS SuperLU_LIBRARY SuperLU_INCLUDE_DIR
)

if(SuperLU_FOUND AND NOT TARGET SuperLU::SuperLU)
    add_library(SuperLU::SuperLU UNKNOWN IMPORTED)
    set_target_properties(SuperLU::SuperLU PROPERTIES
        IMPORTED_LOCATION "${SuperLU_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${SuperLU_INCLUDE_DIR}"
    )
endif()

mark_as_advanced(SuperLU_INCLUDE_DIR SuperLU_LIBRARY)
