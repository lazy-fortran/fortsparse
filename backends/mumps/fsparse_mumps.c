/* Thin C shim over sequential, OpenMP-threaded MUMPS (CeCILL-C). Exposes a
 * factor/solve/free C ABI that Fortran binds through iso_c_binding. The
 * analysis and numerical factorization live in an opaque handle so repeated
 * (and multi-column) solves reuse them.
 *
 * License of this file: MIT (our code). MUMPS is linked, never copied; its
 * headers are only included here. CeCILL-C is a weak-copyleft licence that
 * permits linking from MIT code; the library stays an optional backend.
 *
 * This file binds real double precision; fsparse_zmumps.c binds
 * complex. The matrix is handed to MUMPS in assembled coordinate format (ICNTL(5)=0, ICNTL(18)=0); for a symmetric
 * option only the lower triangle is passed. */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

#include "dmumps_c.h"

#define JOB_INIT (-1)
#define JOB_END (-2)
#define USE_COMM_WORLD (-987654)
#define ICN(k) icntl[(k) - 1]
#define INFOG(k) infog[(k) - 1]
#define RINFOG(k) rinfog[(k) - 1]

void *fsparse_mumps_factor_d(long n, long nnz, const long *colptr1,
                             const long *rowidx1, const double *val, int sym,
                             int ordering, int nthreads, int verbose,
                             int *info1, int *info2);
int fsparse_mumps_solve_d(void *handle, int nrhs, double *rhs, long ld,
                          int *info2);
void fsparse_mumps_stats_d(void *handle, long *istat, double *rstat);
void fsparse_mumps_free_d(void *handle);

typedef struct {
    DMUMPS_STRUC_C id;
    MUMPS_INT *irn;
    MUMPS_INT *jcn;
    double *a;
    long n;
    double t_analysis; /* wall seconds of the symbolic analysis */
    double t_factor;   /* wall seconds of the numerical factorization */
    int ok; /* instance initialised (JOB=-1 done, JOB=-2 still owed) */
} fsparse_mumps_t;

/* MUMPS keeps process-global state, so two threads must never be inside it at
 * once, even on different instances. Every MUMPS call (real and complex) runs
 * under this lock; the parallelism comes from MUMPS' own OpenMP threads. */
static pthread_mutex_t mumps_lock = PTHREAD_MUTEX_INITIALIZER;

void fsparse_mumps_lock(void) { pthread_mutex_lock(&mumps_lock); }
void fsparse_mumps_unlock(void) { pthread_mutex_unlock(&mumps_lock); }

static double wall(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double) ts.tv_sec + 1.0e-9 * (double) ts.tv_nsec;
}

static void mumps_call(fsparse_mumps_t *h, int job)
{
    h->id.job = job;
    fsparse_mumps_lock();
    dmumps_c(&h->id);
    fsparse_mumps_unlock();
}

static void destroy(fsparse_mumps_t *h)
{
    if (h == NULL) return;
    if (h->ok) mumps_call(h, JOB_END);
    free(h->irn);
    free(h->jcn);
    free(h->a);
    free(h);
}

void *fsparse_mumps_factor_d(long n, long nnz, const long *colptr1,
                             const long *rowidx1, const double *val, int sym,
                             int ordering, int nthreads, int verbose,
                             int *info1, int *info2)
{
    fsparse_mumps_t *h;
    long col, e, cnt = 0;
    int attempt, pct;
    double t0;

    *info1 = 0;
    *info2 = 0;
    if (n > 2147483647L || nnz > 2147483647L) {
        *info1 = -1000; /* beyond 32-bit MUMPS integers */
        return NULL;
    }
    h = (fsparse_mumps_t *) calloc(1, sizeof(*h));
    if (h == NULL) {
        *info1 = -1001;
        return NULL;
    }
    h->n = n;
    h->irn = (MUMPS_INT *) malloc((size_t) (nnz > 0 ? nnz : 1) * sizeof(MUMPS_INT));
    h->jcn = (MUMPS_INT *) malloc((size_t) (nnz > 0 ? nnz : 1) * sizeof(MUMPS_INT));
    h->a = (double *) malloc((size_t) (nnz > 0 ? nnz : 1) * sizeof(double));
    if (h->irn == NULL || h->jcn == NULL || h->a == NULL) {
        *info1 = -1001;
        destroy(h);
        return NULL;
    }
    for (col = 0; col < n; ++col) {
        for (e = colptr1[col] - 1; e < colptr1[col + 1] - 1; ++e) {
            long row = rowidx1[e]; /* 1-based */
            if (sym != 0 && row < col + 1) continue; /* lower triangle only */
            h->irn[cnt] = (MUMPS_INT) row;
            h->jcn[cnt] = (MUMPS_INT) (col + 1);
            h->a[cnt] = val[e];
            ++cnt;
        }
    }

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

    /* Quiet by default: errors are reported through INFOG, not stdout. */
    h->id.ICN(1) = verbose ? 6 : -1;
    h->id.ICN(2) = verbose ? 6 : -1;
    h->id.ICN(3) = verbose ? 6 : -1;
    h->id.ICN(4) = verbose ? 2 : 0;
    h->id.ICN(7) = ordering; /* 5 = METIS nested dissection */
    h->id.ICN(10) = 0;       /* refinement is done by the Fortran layer */
    if (nthreads > 0) h->id.ICN(16) = nthreads;

    h->id.n = (MUMPS_INT) n;
    h->id.nnz = (MUMPS_INT8) cnt;
    h->id.irn = h->irn;
    h->id.jcn = h->jcn;
    h->id.a = h->a;

    /* Workspace slack ICNTL(14) grows until MUMPS stops asking for more
     * (INFOG(1) = -9 / -8 / -17 / -20 are memory-estimate failures). */
    t0 = wall();
    mumps_call(h, 1); /* analysis */
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
        mumps_call(h, 2); /* numerical factorization */
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
    /* SYM=1 factors without pivoting and does not itself reject a matrix that
     * is not positive definite; INFOG(12) counts the negative pivots. */
    if (sym == 1 && h->id.INFOG(12) > 0) {
        *info1 = -40;
        *info2 = h->id.INFOG(12);
        destroy(h);
        return NULL;
    }
    /* The coordinate arrays are no longer needed for the solve. */
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

/* Solve with nrhs columns stored contiguously with leading dimension ld.
 * rhs holds the right-hand sides on entry and the solutions on return. */
int fsparse_mumps_solve_d(void *handle, int nrhs, double *rhs, long ld,
                          int *info2)
{
    fsparse_mumps_t *h = (fsparse_mumps_t *) handle;

    h->id.nrhs = nrhs;
    h->id.lrhs = (MUMPS_INT) ld;
    h->id.rhs = rhs;
    mumps_call(h, 3);
    *info2 = h->id.INFOG(2);
    return h->id.INFOG(1);
}

/* istat: [1] factor entries (effective), [2] peak factorization memory (MB,
 * max over processes), [3] effective factorization memory (MB), [4] estimated
 * factor storage (entries), [5] OpenMP threads requested, [6] ordering used.
 * rstat: [1] estimated flops in assembly, [2] elimination flops,
 * [3] factorization flops, [4] analysis wall time, [5] factor wall time. */
void fsparse_mumps_stats_d(void *handle, long *istat, double *rstat)
{
    fsparse_mumps_t *h = (fsparse_mumps_t *) handle;
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

void fsparse_mumps_free_d(void *handle)
{
    destroy((fsparse_mumps_t *) handle);
}
