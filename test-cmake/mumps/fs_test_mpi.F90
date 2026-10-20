module fs_test_mpi
#ifdef FORTSPARSE_MUMPS_MPI
    use mpi, only: MPI_Init_thread, MPI_Finalize, MPI_THREAD_SERIALIZED, MPI_SUCCESS
#endif
    implicit none
    private
    public :: test_mpi_init, test_mpi_finalize, test_mpi_parallel
#ifdef FORTSPARSE_MUMPS_MPI
    logical, parameter :: test_mpi_parallel = .true.
#else
    logical, parameter :: test_mpi_parallel = .false.
#endif
contains
    subroutine test_mpi_init()
#ifdef FORTSPARSE_MUMPS_MPI
        integer :: ierr, provided
        call MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ierr)
        if (ierr /= MPI_SUCCESS .or. provided < MPI_THREAD_SERIALIZED) &
            error stop "MPI initialization failed"
#endif
    end subroutine test_mpi_init

    subroutine test_mpi_finalize()
#ifdef FORTSPARSE_MUMPS_MPI
        integer :: ierr
        call MPI_Finalize(ierr)
        if (ierr /= MPI_SUCCESS) error stop "MPI finalization failed"
#endif
    end subroutine test_mpi_finalize
end module fs_test_mpi
