/* Complex double (ZMUMPS) twin of fsparse_mumps.c: same handle life cycle,
 * same lock, same 32-bit limits. MIT (our code); MUMPS is linked, never copied.
 *
 * Modes (argument `mode`):
 *   0  general: full matrix, LU (SYM=0)
 *   1  Hermitian: A = A^H, lower triangle read, upper mirrored by conjugation,
 *      diagonal taken real, LU on the expanded matrix (SYM=0). MUMPS has no
 *      Hermitian factorization, so this is a general LU that exploits no
 *      symmetry; it is robust for definite and indefinite Hermitian systems.
 *   2  complex symmetric: A = A^T (not A^H), lower triangle, LDL^T (SYM=2).
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "zmumps_c.h"

#define JOB_INIT (-1)
#define JOB_END (-2)
#define USE_COMM_WORLD (-987654)
#define ICN(k) icntl[(k) - 1]
#define INFOG(k) infog[(k) - 1]
#define RINFOG(k) rinfog[(k) - 1]

void fsparse_mumps_lock(void);
void fsparse_mumps_unlock(void);

void *fsparse_mumps_factor_z(long n, long nnz, const long *colptr1,
                             const long *rowidx1,
                             const mumps_double_complex *val, int mode,
                             int ordering, int nthreads, int verbose,
                             int *info1, int *info2);
int fsparse_mumps_solve_z(void *handle, int nrhs, mumps_double_complex *rhs,
                          long ld, int *info2);
void fsparse_mumps_stats_z(void *handle, long *istat, double *rstat);
void fsparse_mumps_free_z(void *handle);

typedef struct {
    ZMUMPS_STRUC_C id;
    MUMPS_INT *irn;
    MUMPS_INT *jcn;
    mumps_double_complex *a;
    long n;
    double t_analysis;
    double t_factor;
    int ok;
} fsparse_zmumps_t;

static double wall(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double) ts.tv_sec + 1.0e-9 * (double) ts.tv_nsec;
}

static void mumps_call(fsparse_zmumps_t *h, int job)
{
    h->id.job = job;
    fsparse_mumps_lock();
    zmumps_c(&h->id);
    fsparse_mumps_unlock();
}

static void destroy(fsparse_zmumps_t *h)
{
    if (h == NULL) return;
    if (h->ok) mumps_call(h, JOB_END);
    free(h->irn);
    free(h->jcn);
    free(h->a);
    free(h);
}

void *fsparse_mumps_factor_z(long n, long nnz, const long *colptr1,
                             const long *rowidx1,
                             const mumps_double_complex *val, int mode,
                             int ordering, int nthreads, int verbose,
                             int *info1, int *info2)
{
    fsparse_zmumps_t *h;
    long col, e, cnt = 0, cap;
    int attempt, pct, sym;
    double t0;

    *info1 = 0;
    *info2 = 0;
    if (mode < 0 || mode > 2) {
        *info1 = -1002;
        return NULL;
    }
    cap = (mode == 1) ? 2 * nnz : nnz;
    if (n > 2147483647L || cap > 2147483647L) {
        *info1 = -1000; /* beyond 32-bit MUMPS integers */
        return NULL;
    }
    h = (fsparse_zmumps_t *) calloc(1, sizeof(*h));
    if (h == NULL) {
        *info1 = -1001;
        return NULL;
    }
    h->n = n;
    h->irn = (MUMPS_INT *) malloc((size_t) (cap > 0 ? cap : 1) * sizeof(MUMPS_INT));
    h->jcn = (MUMPS_INT *) malloc((size_t) (cap > 0 ? cap : 1) * sizeof(MUMPS_INT));
    h->a = (mumps_double_complex *) malloc((size_t) (cap > 0 ? cap : 1) *
                                           sizeof(mumps_double_complex));
    if (h->irn == NULL || h->jcn == NULL || h->a == NULL) {
        *info1 = -1001;
        destroy(h);
        return NULL;
    }
    for (col = 0; col < n; ++col) {
        for (e = colptr1[col] - 1; e < colptr1[col + 1] - 1; ++e) {
            long row = rowidx1[e]; /* 1-based */
            if (mode != 0 && row < col + 1) continue; /* lower triangle only */
            h->irn[cnt] = (MUMPS_INT) row;
            h->jcn[cnt] = (MUMPS_INT) (col + 1);
            h->a[cnt] = val[e];
            if (mode == 1 && row == col + 1) h->a[cnt].i = 0.0;
            ++cnt;
            if (mode == 1 && row != col + 1) { /* conjugate mirror */
                h->irn[cnt] = (MUMPS_INT) (col + 1);
                h->jcn[cnt] = (MUMPS_INT) row;
                h->a[cnt].r = val[e].r;
                h->a[cnt].i = -val[e].i;
                ++cnt;
            }
        }
    }

    sym = (mode == 2) ? 2 : 0;
    h->id.comm_fortran = USE_COMM_WORLD;
    h->id.par = 1;
    h->id.sym = sym;
    mumps_call(h, JOB_INIT);
    if (h->id.INFOG(1) < 0) {
        *info1 = h->id.INFOG(1);
        *info2 = h->id.INFOG(2);
        free(h->irn);
        free(h->jcn);
        free(h->a);
        free(h);
        return NULL;
    }
    h->ok = 1;

    h->id.ICN(1) = verbose ? 6 : -1;
    h->id.ICN(2) = verbose ? 6 : -1;
    h->id.ICN(3) = verbose ? 6 : -1;
    h->id.ICN(4) = verbose ? 2 : 0;
    h->id.ICN(7) = ordering;
    h->id.ICN(10) = 0; /* refinement is done by the Fortran layer */
    if (nthreads > 0) h->id.ICN(16) = nthreads;

    h->id.n = (MUMPS_INT) n;
    h->id.nnz = (MUMPS_INT8) cnt;
    h->id.irn = h->irn;
    h->id.jcn = h->jcn;
    h->id.a = h->a;

    t0 = wall();
    mumps_call(h, 1);
    h->t_analysis = wall() - t0;
    if (h->id.INFOG(1) < 0) {
        *info1 = h->id.INFOG(1);
        *info2 = h->id.INFOG(2);
        destroy(h);
        return NULL;
    }
    pct = 20;
    t0 = wall();
    for (attempt = 0; attempt < 6; ++attempt) {
        h->id.ICN(14) = pct;
        mumps_call(h, 2);
        if (h->id.INFOG(1) != -9 && h->id.INFOG(1) != -8 &&
            h->id.INFOG(1) != -17 && h->id.INFOG(1) != -20)
            break;
        pct *= 2;
    }
    h->t_factor = wall() - t0;
    if (h->id.INFOG(1) < 0) {
        *info1 = h->id.INFOG(1);
        *info2 = h->id.INFOG(2);
        destroy(h);
        return NULL;
    }
    free(h->irn);
    free(h->jcn);
    free(h->a);
    h->irn = NULL;
    h->jcn = NULL;
    h->a = NULL;
    h->id.irn = NULL;
    h->id.jcn = NULL;
    h->id.a = NULL;
    return h;
}

int fsparse_mumps_solve_z(void *handle, int nrhs, mumps_double_complex *rhs,
                          long ld, int *info2)
{
    fsparse_zmumps_t *h = (fsparse_zmumps_t *) handle;

    h->id.nrhs = nrhs;
    h->id.lrhs = (MUMPS_INT) ld;
    h->id.rhs = rhs;
    mumps_call(h, 3);
    *info2 = h->id.INFOG(2);
    return h->id.INFOG(1);
}

/* Same layout as fsparse_mumps_stats_d. */
void fsparse_mumps_stats_z(void *handle, long *istat, double *rstat)
{
    fsparse_zmumps_t *h = (fsparse_zmumps_t *) handle;
    long v;

    v = h->id.INFOG(29);
    istat[0] = v >= 0 ? v : -v * 1000000L;
    istat[1] = h->id.INFOG(16);
    istat[2] = h->id.INFOG(22);
    v = h->id.INFOG(9);
    istat[3] = v >= 0 ? v : -v * 1000000L;
    istat[4] = h->id.ICN(16);
    istat[5] = h->id.INFOG(7);
    rstat[0] = h->id.RINFOG(1);
    rstat[1] = h->id.RINFOG(2);
    rstat[2] = h->id.RINFOG(3);
    rstat[3] = h->t_analysis;
    rstat[4] = h->t_factor;
}

void fsparse_mumps_free_z(void *handle)
{
    destroy((fsparse_zmumps_t *) handle);
}
