# FindUMFPACK.cmake -- locate SuiteSparse/UMFPACK and its config library.
#
# Defines the imported target UMFPACK::UMFPACK carrying the include directory
# (umfpack.h, usually under a suitesparse/ suffix) and both link libraries
# (umfpack and suitesparseconfig). Sets UMFPACK_FOUND on success.

find_path(UMFPACK_INCLUDE_DIR
    NAMES umfpack.h
    PATH_SUFFIXES suitesparse)

find_library(UMFPACK_LIBRARY
    NAMES umfpack)

find_library(SUITESPARSECONFIG_LIBRARY
    NAMES suitesparseconfig)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(UMFPACK
    REQUIRED_VARS UMFPACK_LIBRARY SUITESPARSECONFIG_LIBRARY UMFPACK_INCLUDE_DIR)

if(UMFPACK_FOUND AND NOT TARGET UMFPACK::UMFPACK)
    add_library(UMFPACK::UMFPACK UNKNOWN IMPORTED)
    set_target_properties(UMFPACK::UMFPACK PROPERTIES
        IMPORTED_LOCATION "${UMFPACK_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${UMFPACK_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "${SUITESPARSECONFIG_LIBRARY}")
endif()

mark_as_advanced(UMFPACK_INCLUDE_DIR UMFPACK_LIBRARY SUITESPARSECONFIG_LIBRARY)
