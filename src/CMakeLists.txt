# fortsparse core library (MIT). Links only permissive code: SuperLU (BSD) and
# the MIT C shims. UMFPACK (GPL) is never linked here; it is reached through a
# separate helper executable. Module sources are kept alphabetical and one per
# line to minimise merge conflicts when backends are appended.
add_library(fortsparse STATIC)

target_sources(fortsparse PRIVATE
    # Fortran module sources (alphabetical, one per line; append here)
    fortsparse.f90
    fortsparse_backend.f90
    fortsparse_csc.f90
    fortsparse_kinds.f90
    fortsparse_solver.F90
    fortsparse_status.f90
    fortsparse_version.f90
    ipc/fortsparse_ipc_proto.f90
    ipc/fortsparse_umfpack_ipc.f90
    # MIT C shims
    ipc/fsparse_ipc.c
)

# The SuperLU backend is an optional component. When enabled, compile its
# Fortran binding and C shim, define FORTSPARSE_HAVE_SUPERLU so the solver
# factory references the backend, and link SuperLU. When disabled, none of
# this is compiled or linked and the library needs no SuperLU at all.
if(FORTSPARSE_ENABLE_SUPERLU)
    target_sources(fortsparse PRIVATE
        superlu/fortsparse_superlu.f90
        superlu/fsparse_superlu.c)
    target_compile_definitions(fortsparse PRIVATE FORTSPARSE_HAVE_SUPERLU)
endif()

# Consumers find .mod files in the module output directory; the C shims include
# the shared protocol and IPC headers from the project include/ directory.
target_include_directories(fortsparse
    PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY}
    PRIVATE ${PROJECT_SOURCE_DIR}/include)

# Threads are always needed; SuperLU (BSD) only when its backend is enabled.
# The SuperLU target differs by provisioning: the system FindSuperLU module
# exports SuperLU::SuperLU, the fetched in-tree build exports superlu.
target_link_libraries(fortsparse PUBLIC Threads::Threads)
if(FORTSPARSE_ENABLE_SUPERLU)
    if(FORTSPARSE_USE_SYSTEM_DEPS)
        target_link_libraries(fortsparse PUBLIC SuperLU::SuperLU)
    else()
        target_link_libraries(fortsparse PUBLIC superlu)
    endif()
endif()

# The MUMPS backend (CeCILL-C library, MIT binding). Optional; see the top-level
# FORTSPARSE_ENABLE_MUMPS. The module is compiled with OpenMP when available so
# the refinement residual runs threaded.
if(FORTSPARSE_ENABLE_MUMPS)
    target_sources(fortsparse PRIVATE
        ../backends/mumps/fortsparse_mumps.f90
        ../backends/mumps/fsparse_mumps.c)
    target_compile_definitions(fortsparse PRIVATE FORTSPARSE_HAVE_MUMPS)
    target_link_libraries(fortsparse PUBLIC MUMPS::DMUMPS)
    if(TARGET MUMPS::ZMUMPS)
        target_sources(fortsparse PRIVATE ../backends/mumps/fsparse_zmumps.c)
        target_link_libraries(fortsparse PUBLIC MUMPS::ZMUMPS)
        set(FORTSPARSE_MUMPS_COMPLEX ON PARENT_SCOPE)
    else()
        target_sources(fortsparse PRIVATE
            ../backends/mumps/fsparse_zmumps_stub.c)
    endif()
    # A static METIS built together with MUMPS needs its GKlib after it on the
    # link line; the MUMPS-provided FindMETIS module does not add it.
    find_library(FORTSPARSE_GKLIB NAMES GKlib
        HINTS "${MUMPS_DIR}/.." "${MUMPS_DIR}/../lib")
    if(FORTSPARSE_GKLIB)
        target_link_libraries(fortsparse PUBLIC ${FORTSPARSE_GKLIB})
    endif()
    if(OpenMP_Fortran_FOUND)
        target_link_libraries(fortsparse PUBLIC OpenMP::OpenMP_Fortran)
    endif()
endif()
