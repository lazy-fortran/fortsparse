module fortsparse
    ! Umbrella module. Clients do `use fortsparse` to reach the full public
    ! surface: kinds, status, version, CSC types and builders, and the
    ! backend-pluggable solver API.
    use fortsparse_kinds, only: dp, sp, i4, i8
    use fortsparse_status, only: fortsparse_status_t, status_ok, status_set, &
        FORTSPARSE_OK, FORTSPARSE_SINGULAR, FORTSPARSE_INVALID_MATRIX, &
        FORTSPARSE_BACKEND_UNAVAILABLE, FORTSPARSE_NOT_FACTORED, &
        FORTSPARSE_INTERNAL_ERROR, FORTSPARSE_OUT_OF_MEMORY, &
        FORTSPARSE_INDEX_OVERFLOW, FORTSPARSE_NOT_POSITIVE_DEFINITE, &
        FORTSPARSE_RESOURCE_LIMIT
    use fortsparse_backend, only: FORTSPARSE_MEMORY_AUTO, FORTSPARSE_MEMORY_LOW, &
        FORTSPARSE_OP_NORMAL, FORTSPARSE_OP_TRANSPOSE, FORTSPARSE_OP_ADJOINT
    use fortsparse_memory, only: MEMORY_PLAN_EXACT_LAYOUT, MEMORY_PLAN_BOUND, &
        MEMORY_PLAN_ESTIMATE
    use fortsparse_version, only: fortsparse_version_string
    use fortsparse_ilutp, only: sparse_ilutp_t, sparse_ilutp_stats_t
    use fortsparse_csc, only: csc_t, csc_z_t, csc_from_triplet, csc_is_valid, &
        csc_conjugate_transpose, csc_matmul, csc_matvec, csc_transpose, &
        csc_normalize_symmetric_raw
    use fortsparse_solver, only: sparse_solver_t, sparse_factor, sparse_refactor, &
        sparse_solve, sparse_solve_op, sparse_solve_block_op, sparse_analyze, sparse_memory_plan_t, &
        sparse_free, sparse_destroy, sparse_copy, sparse_solve_once, sparse_vector, &
        sparse_solve_jvp, sparse_solve_vjp, sparse_solve_block, &
        sparse_factor_stats, sparse_options_t, sparse_factor_stats_t, &
        FORTSPARSE_BACKEND_SUPERLU, FORTSPARSE_BACKEND_UMFPACK_IPC, &
        FORTSPARSE_BACKEND_MUMPS, FORTSPARSE_BACKEND_DEFAULT, &
        FORTSPARSE_BACKEND_AUTO, FORTSPARSE_BACKEND_NATIVE, &
        FORTSPARSE_BACKEND_NATIVE_REFERENCE, sparse_auto_backend, &
        FORTSPARSE_SYM_GENERAL, FORTSPARSE_SYM_SPD, &
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
    public :: FORTSPARSE_NOT_POSITIVE_DEFINITE, &
        FORTSPARSE_RESOURCE_LIMIT

    ! Version
    public :: fortsparse_version_string

    ! Storage
    public :: csc_t, csc_z_t, csc_from_triplet, csc_is_valid, csc_matmul
    public :: csc_matvec, csc_transpose
    public :: csc_conjugate_transpose, csc_normalize_symmetric_raw

    ! Solver
    public :: sparse_ilutp_t, sparse_ilutp_stats_t
    public :: sparse_solver_t, sparse_factor, sparse_refactor, sparse_solve, &
        sparse_free
    public :: sparse_destroy, sparse_copy, sparse_solve_once, sparse_vector
    public :: sparse_solve_jvp, sparse_solve_vjp
    public :: FORTSPARSE_BACKEND_SUPERLU
    public :: FORTSPARSE_BACKEND_UMFPACK_IPC, FORTSPARSE_BACKEND_MUMPS
    public :: FORTSPARSE_BACKEND_DEFAULT, FORTSPARSE_BACKEND_AUTO
    public :: FORTSPARSE_BACKEND_NATIVE, FORTSPARSE_BACKEND_NATIVE_REFERENCE
    public :: sparse_auto_backend
    public :: sparse_solve_block, sparse_factor_stats
    public :: FORTSPARSE_MEMORY_AUTO, FORTSPARSE_MEMORY_LOW
    public :: sparse_solve_op, sparse_solve_block_op
    public :: FORTSPARSE_OP_NORMAL, FORTSPARSE_OP_TRANSPOSE, FORTSPARSE_OP_ADJOINT
    public :: sparse_analyze, sparse_memory_plan_t
    public :: MEMORY_PLAN_EXACT_LAYOUT, MEMORY_PLAN_BOUND, MEMORY_PLAN_ESTIMATE
    public :: sparse_options_t, sparse_factor_stats_t
    public :: FORTSPARSE_SYM_GENERAL, FORTSPARSE_SYM_SPD
    public :: FORTSPARSE_SYM_INDEFINITE, FORTSPARSE_SYM_HERMITIAN
    public :: FORTSPARSE_ORDER_DEFAULT, FORTSPARSE_ORDER_AMD, FORTSPARSE_ORDER_ND
    public :: native_ordering_available

contains

    logical function native_ordering_available(ordering) result(available)
        ! Built native graph capability only; matrix class, device and memory
        ! policy are still checked by factorization. No graph is allocated here.
        use, intrinsic :: iso_c_binding, only: c_int
        integer, intent(in) :: ordering
        interface
            integer(c_int) function amd_available() &
                bind(C, name='fsparse_native_order_amd_available')
                import c_int
            end function
            integer(c_int) function nd_available() &
                bind(C, name='fsparse_native_order_available')
                import c_int
            end function
        end interface
        select case(ordering)
        case(FORTSPARSE_ORDER_DEFAULT)
            available = .true.
        case(FORTSPARSE_ORDER_AMD)
            available = amd_available() == 1
        case(FORTSPARSE_ORDER_ND)
            available = nd_available() == 1
        case default
            available = .false.
        end select
    end function native_ordering_available

end module fortsparse
