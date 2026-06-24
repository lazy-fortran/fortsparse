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
    fortsparse_solver.f90
    fortsparse_status.f90
    fortsparse_version.f90
    ipc/fortsparse_ipc_proto.f90
    ipc/fortsparse_umfpack_ipc.f90
    superlu/fortsparse_superlu.f90
    # MIT C shims
    ipc/fsparse_ipc.c
    superlu/fsparse_superlu.c
)

# Consumers find .mod files in the module output directory; the C shims include
# the shared protocol and IPC headers from the project include/ directory.
target_include_directories(fortsparse
    PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY}
    PRIVATE ${PROJECT_SOURCE_DIR}/include)

# Permissive backends only: SuperLU (BSD) plus threads for portability.
target_link_libraries(fortsparse PUBLIC SuperLU::SuperLU Threads::Threads)
