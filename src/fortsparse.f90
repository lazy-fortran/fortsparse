module fortsparse
    ! Umbrella module. Clients do `use fortsparse` to reach the full public
    ! surface: kinds, status, version, CSC types and builders, and the
    ! backend-pluggable solver API.
    use fortsparse_kinds, only: dp, sp, i4, i8
    use fortsparse_status, only: fortsparse_status_t, status_ok, status_set, &
        FORTSPARSE_OK, FORTSPARSE_SINGULAR, FORTSPARSE_INVALID_MATRIX, &
        FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_INTERNAL_ERROR
    use fortsparse_version, only: fortsparse_version_string
    use fortsparse_csc, only: csc_t, csc_z_t, csc_from_triplet, csc_is_valid, &
        csc_matvec
    use fortsparse_solver, only: sparse_solver_t, sparse_factor, sparse_solve, &
        sparse_free, sparse_solve_once, FORTSPARSE_BACKEND_SUPERLU, &
        FORTSPARSE_BACKEND_UMFPACK_IPC
    implicit none
    private

    ! Kinds
    public :: dp, sp, i4, i8

    ! Status
    public :: fortsparse_status_t, status_ok, status_set
    public :: FORTSPARSE_OK, FORTSPARSE_SINGULAR, FORTSPARSE_INVALID_MATRIX
    public :: FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_NOT_FACTORED
    public :: FORTSPARSE_INTERNAL_ERROR

    ! Version
    public :: fortsparse_version_string

    ! Storage
    public :: csc_t, csc_z_t, csc_from_triplet, csc_is_valid, csc_matvec

    ! Solver
    public :: sparse_solver_t, sparse_factor, sparse_solve, sparse_free
    public :: sparse_solve_once, FORTSPARSE_BACKEND_SUPERLU
    public :: FORTSPARSE_BACKEND_UMFPACK_IPC

end module fortsparse
