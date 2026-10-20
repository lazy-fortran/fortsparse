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
#define ICN(k) icntl[(k) - 1]
/* MUMPS 5.9 L0-OMP tree parallelism (ICNTL(48)=1, its default) parks the
 * threads that finish their subtrees in usleep(1000) (dfac_omp_m.F,
 * DMUMPS_PERFORM_COPIES), a fixed ~1 ms per factorization with two or more
 * threads. It pays off only from this order on. */
#define FSPARSE_MUMPS_L0_MIN_N 5000
/* Automatic ordering (ordering < 0): AMD first; when the AMD analysis
 * estimates more than this many elimination flops per stored entry (a
 * three-dimensional or strongly coupled pattern), METIS nested dissection
 * is analysed instead, whose fill is then far smaller than its cost. */
#define FSPARSE_MUMPS_ND_FLOPS_PER_NNZ 3000.0
#define INFOG(k) infog[(k) - 1]
#define RINFOG(k) rinfog[(k) - 1]

void fsparse_mumps_lock(void);
void fsparse_mumps_unlock(void);
int fsparse_mumps_mpi_ready(void);
MUMPS_INT fsparse_mumps_comm(void);

void *fsparse_mumps_factor_z(long n, long nnz, const int *colptr1,
                             const int *rowidx1,
                             const mumps_double_complex *val, int mode,
                             int ordering, int nthreads, int verbose,
                             int *info1, int *info2);
int fsparse_mumps_solve_z(void *handle, int nrhs, mumps_double_complex *rhs,
                          long ld, int *info2);
void fsparse_mumps_stats_z(void *handle, long *istat, double *rstat);
int fsparse_mumps_refactor_z(void *handle, long nnz, const int *colptr1,
                            const int *rowidx1,
                            const mumps_double_complex *val, int *info1,
                            int *info2);
void fsparse_mumps_free_z(void *handle);

typedef struct {
    ZMUMPS_STRUC_C id;
    MUMPS_INT *irn;
    MUMPS_INT *jcn;
    mumps_double_complex *a;
    long n;
    double t_analysis;
    double t_factor;
    int mode;  /* fsparse_mumps_factor_z mode */
    int pct;   /* MUMPS ICNTL(14) slack that last worked */
    int sym;   /* MUMPS SYM used */
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
    if (h->ok && fsparse_mumps_mpi_ready()) mumps_call(h, JOB_END);
    free(h->irn);
    free(h->jcn);
    free(h->a);
    free(h);
}

/* Numerical factorization with ICNTL(14) slack that grows until MUMPS stops
 * asking for more workspace (INFOG(1) = -9 / -8 / -17 / -20). Returns 0 or -1
 * with info1/info2 set. The handle stays valid; the caller decides whether to
 * destroy it. */
static int numeric(fsparse_zmumps_t *h, int *info1, int *info2)
{
    int attempt, pct = h->pct;
    double t0 = wall();

    for (attempt = 0; attempt < 6; ++attempt) {
        h->id.ICN(14) = pct;
        mumps_call(h, 2); /* numerical factorization */
        if (h->id.INFOG(1) != -9 && h->id.INFOG(1) != -8 &&
            h->id.INFOG(1) != -17 && h->id.INFOG(1) != -20)
            break;
        pct *= 2;
    }
    h->pct = pct;
    h->t_factor = wall() - t0;
    if (h->id.INFOG(1) < 0) {
        *info1 = h->id.INFOG(1);
        *info2 = h->id.INFOG(2);
        return -1;
    }
    /* SYM=1 factors without pivoting and does not itself reject a matrix that
     * is not positive definite; INFOG(12) counts the negative pivots. */
    if (h->sym == 1 && h->id.INFOG(12) > 0) {
        *info1 = -40;
        *info2 = h->id.INFOG(12);
        return -1;
    }
    return 0;
}

void *fsparse_mumps_factor_z(long n, long nnz, const int *colptr1,
                             const int *rowidx1,
                             const mumps_double_complex *val, int mode,
                             int ordering, int nthreads, int verbose,
                             int *info1, int *info2)
{
    fsparse_zmumps_t *h;
    long col, e, cnt = 0, cap;
    int sym;
    double t0;

    *info1 = 0;
    *info2 = 0;
    if (!fsparse_mumps_mpi_ready()) {
        *info1 = -1004;
        return NULL;
    }
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
    h->id.comm_fortran = fsparse_mumps_comm();
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
    h->id.ICN(7) = ordering < 0 ? 0 : ordering; /* 0 AMD, 5 METIS */
    h->id.ICN(10) = 0; /* refinement is done by the Fortran layer */
    if (nthreads > 0) h->id.ICN(16) = nthreads;
    /* JOB_INIT disables L0-OMP when MUMPS was built without OpenMP. */
    if (n < FSPARSE_MUMPS_L0_MIN_N) h->id.ICN(48) = 0;

    h->id.n = (MUMPS_INT) n;
    h->id.nnz = (MUMPS_INT8) cnt;
    h->id.irn = h->irn;
    h->id.jcn = h->jcn;
    h->id.a = h->a;

    t0 = wall();
    mumps_call(h, 1); /* analysis */
    if (ordering < 0 && h->id.INFOG(1) >= 0 &&
        h->id.RINFOG(1) > FSPARSE_MUMPS_ND_FLOPS_PER_NNZ * (double) nnz) {
        h->id.ICN(7) = 5;
        mumps_call(h, 1);
    }
    h->t_analysis = wall() - t0;
    if (h->id.INFOG(1) < 0) {
        *info1 = h->id.INFOG(1);
        *info2 = h->id.INFOG(2);
        destroy(h);
        return NULL;
    }
    h->pct = 20;
    h->sym = sym;
    h->mode = mode;
    if (numeric(h, info1, info2) != 0) {
        destroy(h);
        return NULL;
    }
    /* The coordinate arrays stay: fsparse_mumps_refactor_z refills the values
     * and reuses this analysis. */
    return h;
}

/* Numerical refactorization of the same sparsity pattern with new values: the
 * symbolic analysis is kept. The caller guarantees colptr1/rowidx1 equal the
 * pattern of the original factorization (the Fortran layer checks it); a
 * differing entry count is still caught here. Returns 0, or the MUMPS error in
 * info1 (the handle then holds no usable factors and must be freed). */
int fsparse_mumps_refactor_z(void *handle, long nnz, const int *colptr1,
                            const int *rowidx1, const mumps_double_complex *val,
                            int *info1, int *info2)
{
    fsparse_zmumps_t *h = (fsparse_zmumps_t *) handle;
    long col, e, cnt = 0;

    *info1 = 0;
    *info2 = 0;
    if (!fsparse_mumps_mpi_ready()) {
        *info1 = -1004;
        return -1;
    }
    (void) nnz;
    if (h == NULL || h->a == NULL) {
        *info1 = -1003;
        return -1;
    }
    for (col = 0; col < h->n; ++col) {
        for (e = colptr1[col] - 1; e < colptr1[col + 1] - 1; ++e) {
            long row = rowidx1[e];
            if (h->mode != 0 && row < col + 1) continue;
            if (cnt >= h->id.nnz) goto mismatch;
            h->a[cnt] = val[e];
            if (h->mode == 1 && row == col + 1) h->a[cnt].i = 0.0;
            ++cnt;
            if (h->mode == 1 && row != col + 1) {
                if (cnt >= h->id.nnz) goto mismatch;
                h->a[cnt].r = val[e].r;
                h->a[cnt].i = -val[e].i;
                ++cnt;
            }
        }
    }
    if (cnt != h->id.nnz) goto mismatch;
    h->id.a = h->a;
    h->t_analysis = 0.0; /* statistics describe this call: no analysis */
    return numeric(h, info1, info2);
mismatch:
    *info1 = -1003; /* pattern does not match the analysed one */
    return -1;
}

int fsparse_mumps_solve_op_z(void *handle, int nrhs, mumps_double_complex *rhs,
                          long ld, int operation, int *info2)
{
    fsparse_zmumps_t *h = (fsparse_zmumps_t *) handle;

    if (!fsparse_mumps_mpi_ready()) {
        *info2 = 0;
        return -1004;
    }
    if (operation < 0 || operation > 2) { *info2 = 0; return -1005; }
    h->id.nrhs = nrhs;
    h->id.lrhs = (MUMPS_INT) ld;
    h->id.rhs = rhs;
    /* MUMPS zsol_driver.F: ICNTL(9)=1 is A, other values are A^T.
     * A^H x=b is conj(A^T\conj(b)); symmetric mode ignores ICNTL(9). */
    h->id.ICN(9) = operation == 0 ? 1 : 0;
    if (operation == 2)
        for (int j = 0; j < nrhs; ++j)
            for (long i = 0; i < h->id.n; ++i) rhs[j*ld+i].i = -rhs[j*ld+i].i;
    mumps_call(h, 3);
    if (operation == 2)
        for (int j = 0; j < nrhs; ++j)
            for (long i = 0; i < h->id.n; ++i) rhs[j*ld+i].i = -rhs[j*ld+i].i;
    h->id.ICN(9) = 1;
    *info2 = h->id.INFOG(2);
    return h->id.INFOG(1);
}

int fsparse_mumps_solve_z(void *handle, int nrhs, mumps_double_complex *rhs,
                          long ld, int *info2)
{
    return fsparse_mumps_solve_op_z(handle, nrhs, rhs, ld, 0, info2);
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
