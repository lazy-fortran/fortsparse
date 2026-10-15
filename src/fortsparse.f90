module fortsparse
    ! Umbrella module. Clients do `use fortsparse` to reach the full public
    ! surface: kinds, status, version, CSC types and builders, and the
    ! backend-pluggable solver API.
    use fortsparse_kinds, only: dp, sp, i4, i8
    use fortsparse_status, only: fortsparse_status_t, status_ok, status_set, &
        FORTSPARSE_OK, FORTSPARSE_SINGULAR, FORTSPARSE_INVALID_MATRIX, &
        FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_INTERNAL_ERROR, FORTSPARSE_OUT_OF_MEMORY, &
        FORTSPARSE_INDEX_OVERFLOW
    use fortsparse_version, only: fortsparse_version_string
    use fortsparse_csc, only: csc_t, csc_z_t, csc_from_triplet, csc_is_valid, &
        csc_conjugate_transpose, csc_matmul, csc_matvec, csc_transpose
    use fortsparse_solver, only: sparse_solver_t, sparse_factor, sparse_refactor, &
        sparse_solve, &
        sparse_free, sparse_destroy, sparse_solve_once, sparse_vector, &
        sparse_solve_jvp, sparse_solve_vjp, sparse_solve_block, &
        sparse_factor_stats, sparse_options_t, sparse_factor_stats_t, &
        FORTSPARSE_BACKEND_SUPERLU, FORTSPARSE_BACKEND_UMFPACK_IPC, &
        FORTSPARSE_BACKEND_MUMPS, FORTSPARSE_SYM_GENERAL, FORTSPARSE_SYM_SPD, &
        FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_SYM_HERMITIAN, &
        FORTSPARSE_ORDER_DEFAULT, FORTSPARSE_ORDER_AMD, FORTSPARSE_ORDER_ND
    implicit none
    private

    ! Kinds
    public :: dp, sp, i4, i8

    ! Status
    public :: fortsparse_status_t, status_ok, status_set
    public :: FORTSPARSE_OK, FORTSPARSE_SINGULAR, FORTSPARSE_INVALID_MATRIX
    public :: FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_NOT_FACTORED
    public :: FORTSPARSE_INTERNAL_ERROR
    public :: FORTSPARSE_OUT_OF_MEMORY, FORTSPARSE_INDEX_OVERFLOW

    ! Version
    public :: fortsparse_version_string

    ! Storage
    public :: csc_t, csc_z_t, csc_from_triplet, csc_is_valid, csc_matmul
    public :: csc_matvec, csc_transpose
    public :: csc_conjugate_transpose

    ! Solver
    public :: sparse_solver_t, sparse_factor, sparse_refactor, sparse_solve, &
        sparse_free
    public :: sparse_destroy, sparse_solve_once, sparse_vector
    public :: sparse_solve_jvp, sparse_solve_vjp
    public :: FORTSPARSE_BACKEND_SUPERLU
    public :: FORTSPARSE_BACKEND_UMFPACK_IPC, FORTSPARSE_BACKEND_MUMPS
    public :: sparse_solve_block, sparse_factor_stats
    public :: sparse_options_t, sparse_factor_stats_t
    public :: FORTSPARSE_SYM_GENERAL, FORTSPARSE_SYM_SPD
    public :: FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_SYM_HERMITIAN
    public :: FORTSPARSE_ORDER_DEFAULT, FORTSPARSE_ORDER_AMD, FORTSPARSE_ORDER_ND

end module fortsparse
