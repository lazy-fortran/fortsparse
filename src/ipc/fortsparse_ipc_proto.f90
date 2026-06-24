module fortsparse_ipc_proto
    ! Shared-memory protocol mirror for the out-of-process backend (MIT). The
    ! bind(c) header type matches struct fsparse_shm_header in
    ! include/fsparse_proto.h field for field, so the Fortran side and the GPL
    ! helper agree on the layout. Also resolves the helper binary from the
    ! environment or PATH. No UMFPACK symbol appears here.
    use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_char, &
        c_size_t, c_null_char
    implicit none
    private

    ! Directory of the running executable, from the MIT C shim. Lets a program
    ! find the helper next to it with no PATH entry and no environment variable.
    interface
        function fsparse_self_dir(buf, n) result(length) bind(c, &
                name="fsparse_self_dir")
            import :: c_char, c_size_t
            character(kind=c_char), intent(out) :: buf(*)
            integer(c_size_t), value            :: n
            integer(c_size_t)                   :: length
        end function fsparse_self_dir
    end interface

    public :: shm_header_t
    public :: OP_FACTOR_REAL, OP_FACTOR_COMPLEX, OP_SOLVE_REAL, &
        OP_SOLVE_COMPLEX, OP_FREE, OP_SHUTDOWN
    public :: ST_OK, ST_SINGULAR, ST_ERROR
    public :: find_helper

    ! Opcodes and normalized status codes, matching fsparse_proto.h.
    integer(c_int32_t), parameter :: OP_FACTOR_REAL = 1
    integer(c_int32_t), parameter :: OP_FACTOR_COMPLEX = 2
    integer(c_int32_t), parameter :: OP_SOLVE_REAL = 3
    integer(c_int32_t), parameter :: OP_SOLVE_COMPLEX = 4
    integer(c_int32_t), parameter :: OP_FREE = 5
    integer(c_int32_t), parameter :: OP_SHUTDOWN = 6
    integer(c_int32_t), parameter :: ST_OK = 0
    integer(c_int32_t), parameter :: ST_SINGULAR = 1
    integer(c_int32_t), parameter :: ST_ERROR = 2

    character(*), parameter :: HELPER_NAME = "fortsparse_umfpack_helper"

    ! Header at offset 0 of the shared mapping. Field order and kinds mirror
    ! struct fsparse_shm_header exactly so c_f_pointer aliases it safely.
    type, bind(c) :: shm_header_t
        integer(c_int32_t) :: opcode
        integer(c_int32_t) :: status
        integer(c_int64_t) :: n
        integer(c_int64_t) :: nnz
        integer(c_int32_t) :: refine
        integer(c_int32_t) :: is_complex
        integer(c_int64_t) :: off_colptr
        integer(c_int64_t) :: off_rowidx
        integer(c_int64_t) :: off_ax
        integer(c_int64_t) :: off_az
        integer(c_int64_t) :: off_b
        integer(c_int64_t) :: off_bz
        integer(c_int64_t) :: off_x
        integer(c_int64_t) :: off_xz
        integer(c_int64_t) :: data_bytes
    end type shm_header_t

contains

    ! Resolve the helper path. Try, in order, FORTSPARSE_UMFPACK_HELPER, the
    ! directory of the running executable, then the PATH directories. Returns an
    ! empty string when no helper is found.
    function find_helper() result(path)
        character(:), allocatable :: path

        path = helper_from_env()
        if (len(path) > 0) return
        path = helper_from_exe_dir()
        if (len(path) > 0) return
        path = helper_from_path()
    end function find_helper

    ! Build <exe_dir>/fortsparse_umfpack_helper and return it if it exists.
    function helper_from_exe_dir() result(path)
        character(:), allocatable :: path
        character(kind=c_char)    :: buf(4096)
        integer(c_size_t)         :: n
        integer                   :: i
        logical                   :: ok

        path = ""
        n = fsparse_self_dir(buf, int(size(buf), c_size_t))
        if (n <= 0_c_size_t) return
        do i = 1, int(n)
            path = path//buf(i)
        end do
        path = path//"/"//HELPER_NAME
        inquire (file=path, exist=ok)
        if (.not. ok) path = ""
    end function helper_from_exe_dir

    ! Read FORTSPARSE_UMFPACK_HELPER if it names an existing file.
    function helper_from_env() result(path)
        character(:), allocatable :: path
        character(4096)           :: buf
        integer                   :: n
        logical                   :: ok

        path = ""
        call get_environment_variable("FORTSPARSE_UMFPACK_HELPER", buf, n)
        if (n <= 0) return
        inquire (file=buf(1:n), exist=ok)
        if (ok) path = buf(1:n)
    end function helper_from_env

    ! Search the PATH directories for an executable helper of the canonical
    ! name and return the first existing match.
    function helper_from_path() result(path)
        character(:), allocatable :: path
        character(:), allocatable :: env
        character(4096)           :: buf
        integer                   :: n, lo, hi, sep
        logical                   :: ok

        path = ""
        call get_environment_variable("PATH", buf, n)
        if (n <= 0) return
        env = buf(1:n)
        lo = 1
        do while (lo <= len(env))
            sep = index(env(lo:), ":")
            if (sep == 0) then
                hi = len(env)
            else
                hi = lo + sep - 2
            end if
            if (hi >= lo) then
                path = env(lo:hi)//"/"//HELPER_NAME
                inquire (file=path, exist=ok)
                if (ok) return
            end if
            if (sep == 0) exit
            lo = hi + 2
        end do
        path = ""
    end function helper_from_path

end module fortsparse_ipc_proto
