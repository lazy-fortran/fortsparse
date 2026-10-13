# FindMUMPS.cmake -- locate a MUMPS installation that ships no CMake package,
# such as the distribution packages (Arch, Debian, Fedora). The top level tries
# the scivision/mumps package config first and falls back to this module.
#
# Distribution MUMPS libraries do not record all their dependencies, so the
# targets carry them explicitly: METIS, SCOTCH/esmumps, LAPACK/BLAS, and for an
# MPI build ScaLAPACK and MPI. A sequential build is recognised by the libmpiseq
# shim installed beside libdmumps.
#
# Hints: MUMPS_ROOT, or MUMPS_DIR pointing at the prefix or its cmake/ folder.
#
# Defines:
#   MUMPS::DMUMPS     double-precision MUMPS with its dependencies
#   MUMPS::ZMUMPS     double-complex MUMPS, when libzmumps is installed
#   MUMPS_parallel    ON for an MPI build, OFF for a sequential one
#   MUMPS_FOUND, MUMPS_INCLUDE_DIR

set(_mumps_hints)
if(MUMPS_DIR)
    list(APPEND _mumps_hints "${MUMPS_DIR}" "${MUMPS_DIR}/..")
endif()

find_path(MUMPS_INCLUDE_DIR NAMES dmumps_c.h
    HINTS ${_mumps_hints} PATH_SUFFIXES include include/mumps include/MUMPS)
find_library(MUMPS_DMUMPS_LIBRARY NAMES dmumps
    HINTS ${_mumps_hints} PATH_SUFFIXES lib lib64)
find_library(MUMPS_COMMON_LIBRARY NAMES mumps_common
    HINTS ${_mumps_hints} PATH_SUFFIXES lib lib64)

set(_mumps_deps)
if(MUMPS_DMUMPS_LIBRARY)
    get_filename_component(_mumps_libdir "${MUMPS_DMUMPS_LIBRARY}" DIRECTORY)
    find_library(MUMPS_ZMUMPS_LIBRARY NAMES zmumps HINTS "${_mumps_libdir}")
    find_library(MUMPS_PORD_LIBRARY NAMES pord HINTS "${_mumps_libdir}")
    # Only a shim next to libdmumps marks that installation as sequential.
    find_library(MUMPS_MPISEQ_LIBRARY NAMES mpiseq
        PATHS "${_mumps_libdir}" NO_DEFAULT_PATH)
    find_library(MUMPS_MPISEQ_C_LIBRARY NAMES mpiseq_c
        PATHS "${_mumps_libdir}" NO_DEFAULT_PATH)
    find_library(MUMPS_MPISEQ_FORTRAN_LIBRARY NAMES mpiseq_fortran
        PATHS "${_mumps_libdir}" NO_DEFAULT_PATH)
    list(APPEND _mumps_deps "${MUMPS_COMMON_LIBRARY}")
    if(MUMPS_PORD_LIBRARY)
        list(APPEND _mumps_deps "${MUMPS_PORD_LIBRARY}")
    endif()

    # Orderings MUMPS may have been built with; harmless when unused.
    find_library(MUMPS_METIS_LIBRARY NAMES metis HINTS "${_mumps_libdir}")
    foreach(_lib esmumps scotch scotcherr)
        find_library(MUMPS_${_lib}_LIBRARY NAMES ${_lib} HINTS "${_mumps_libdir}")
        mark_as_advanced(MUMPS_${_lib}_LIBRARY)
    endforeach()
    foreach(_var METIS esmumps scotch scotcherr)
        if(MUMPS_${_var}_LIBRARY)
            list(APPEND _mumps_deps "${MUMPS_${_var}_LIBRARY}")
        endif()
    endforeach()
endif()

set(_mumps_parallel_vars)
if(MUMPS_MPISEQ_LIBRARY OR MUMPS_MPISEQ_C_LIBRARY)
    set(MUMPS_parallel OFF)
    foreach(_var MPISEQ MPISEQ_C MPISEQ_FORTRAN)
        if(MUMPS_${_var}_LIBRARY)
            list(APPEND _mumps_deps "${MUMPS_${_var}_LIBRARY}")
        endif()
    endforeach()
elseif(MUMPS_DMUMPS_LIBRARY)
    set(MUMPS_parallel ON)
    # FindMPI reuses MPI targets a parent project already created.
    find_package(MPI COMPONENTS C Fortran)
    find_library(MUMPS_SCALAPACK_LIBRARY
        NAMES scalapack scalapack-openmpi scalapack-mpich
        HINTS "${_mumps_libdir}")
    set(_mumps_parallel_vars MUMPS_SCALAPACK_LIBRARY MPI_C_FOUND MPI_Fortran_FOUND)
endif()

if(MUMPS_DMUMPS_LIBRARY AND NOT LAPACK_FOUND)
    find_package(LAPACK)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS
    REQUIRED_VARS MUMPS_DMUMPS_LIBRARY MUMPS_COMMON_LIBRARY MUMPS_INCLUDE_DIR
        LAPACK_FOUND ${_mumps_parallel_vars})

if(MUMPS_FOUND)
    if(MUMPS_parallel)
        list(APPEND _mumps_deps "${MUMPS_SCALAPACK_LIBRARY}")
    endif()
    list(APPEND _mumps_deps ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    if(MUMPS_parallel)
        list(APPEND _mumps_deps MPI::MPI_Fortran MPI::MPI_C)
    endif()
    foreach(_arith DMUMPS ZMUMPS)
        if(MUMPS_${_arith}_LIBRARY AND NOT TARGET MUMPS::${_arith})
            add_library(MUMPS::${_arith} UNKNOWN IMPORTED)
            set_target_properties(MUMPS::${_arith} PROPERTIES
                IMPORTED_LOCATION "${MUMPS_${_arith}_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_INCLUDE_DIR}"
                INTERFACE_LINK_LIBRARIES "${_mumps_deps}")
        endif()
    endforeach()
endif()

mark_as_advanced(MUMPS_INCLUDE_DIR MUMPS_DMUMPS_LIBRARY MUMPS_ZMUMPS_LIBRARY
    MUMPS_COMMON_LIBRARY MUMPS_PORD_LIBRARY MUMPS_MPISEQ_LIBRARY
    MUMPS_MPISEQ_C_LIBRARY MUMPS_MPISEQ_FORTRAN_LIBRARY MUMPS_METIS_LIBRARY
    MUMPS_SCALAPACK_LIBRARY)
