/* Link stub used when the MUMPS installation has no complex (ZMUMPS) library.
 * The complex entry points report "unavailable" (info1 = -2000) so the Fortran
 * layer returns FORTSPARSE_BACKEND_UNAVAILABLE instead of failing to link.
 * MIT (our code). */

#include <stddef.h>

void *fsparse_mumps_factor_z(long n, long nnz, const int *colptr1,
                             const int *rowidx1, const void *val, int mode,
                             int ordering, int nthreads, int verbose,
                             int *info1, int *info2)
{
    (void) n; (void) nnz; (void) colptr1; (void) rowidx1; (void) val;
    (void) mode; (void) ordering; (void) nthreads; (void) verbose;
    *info1 = -2000;
    *info2 = 0;
    return NULL;
}

int fsparse_mumps_refactor_z(void *handle, long nnz, const int *colptr1,
                            const int *rowidx1, const void *val, int *info1,
                            int *info2)
{
    (void) handle; (void) nnz; (void) colptr1; (void) rowidx1; (void) val;
    *info1 = -2000;
    *info2 = 0;
    return -1;
}

int fsparse_mumps_solve_z(void *handle, int nrhs, void *rhs, long ld,
                          int *info2)
{
    (void) handle; (void) nrhs; (void) rhs; (void) ld;
    *info2 = 0;
    return -2000;
}

void fsparse_mumps_stats_z(void *handle, long *istat, double *rstat)
{
    (void) handle; (void) istat; (void) rstat;
}

void fsparse_mumps_free_z(void *handle) { (void) handle; }
