/* Out-of-process UMFPACK helper (GPL). SEPARATE standalone executable.
 *
 * This is the ONLY file that includes umfpack.h and links -lumfpack. It never
 * becomes part of libfortsparse. The MIT library drives it at arm's length
 * over the shared memory and named semaphores created by the parent process
 * (see include/fsparse_ipc.h, include/fsparse_proto.h). FSF treats a separate
 * process behind an IPC boundary as a GPL-compatible separation, so the GPL of
 * UMFPACK stays confined to this binary.
 *
 * argv: [shm_name, sem_req_name, sem_done_name, bytes]. The helper attaches to
 * those objects, then loops: wait on req, dispatch on header->opcode, post
 * done. The Symbolic/Numeric factors stay resident across solves.
 *
 * License of this file: GPL (it links GPL UMFPACK). */

#include "fsparse_proto.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <suitesparse/umfpack.h>

#if defined(_WIN32)

#include <windows.h>

typedef struct {
    HANDLE map;
    HANDLE sem_req;
    HANDLE sem_done;
    void *region;
} helper_ctx;

static int attach(helper_ctx *c, char **argv, int64_t bytes)
{
    c->map = OpenFileMappingA(FILE_MAP_ALL_ACCESS, FALSE, argv[1]);
    if (c->map == NULL) return 1;
    c->region = MapViewOfFile(c->map, FILE_MAP_ALL_ACCESS, 0, 0,
                              (SIZE_T) bytes);
    if (c->region == NULL) return 1;
    c->sem_req = OpenSemaphoreA(SEMAPHORE_ALL_ACCESS, FALSE, argv[2]);
    c->sem_done = OpenSemaphoreA(SEMAPHORE_ALL_ACCESS, FALSE, argv[3]);
    if (c->sem_req == NULL || c->sem_done == NULL) return 1;
    return 0;
}

static void wait_req(helper_ctx *c)
{
    WaitForSingleObject(c->sem_req, INFINITE);
}

static void post_done(helper_ctx *c)
{
    ReleaseSemaphore(c->sem_done, 1, NULL);
}

static void detach(helper_ctx *c)
{
    if (c->region != NULL) UnmapViewOfFile(c->region);
    if (c->map != NULL) CloseHandle(c->map);
    if (c->sem_req != NULL) CloseHandle(c->sem_req);
    if (c->sem_done != NULL) CloseHandle(c->sem_done);
}

#else

#include <fcntl.h>
#include <semaphore.h>
#include <sys/mman.h>
#include <unistd.h>

typedef struct {
    int fd;
    sem_t *sem_req;
    sem_t *sem_done;
    void *region;
    int64_t bytes;
} helper_ctx;

static int attach(helper_ctx *c, char **argv, int64_t bytes)
{
    /* Set the failure sentinels first: detach() tests against SEM_FAILED and
     * fd < 0, which a zeroed struct does not portably provide (SEM_FAILED is
     * not guaranteed to be NULL). */
    c->bytes = bytes;
    c->region = NULL;
    c->fd = -1;
    c->sem_req = SEM_FAILED;
    c->sem_done = SEM_FAILED;
    c->fd = shm_open(argv[1], O_RDWR, 0600);
    if (c->fd < 0) return 1;
    c->region = mmap(NULL, (size_t) bytes, PROT_READ | PROT_WRITE, MAP_SHARED,
                     c->fd, 0);
    if (c->region == MAP_FAILED) { c->region = NULL; return 1; }
    c->sem_req = sem_open(argv[2], 0);
    c->sem_done = sem_open(argv[3], 0);
    if (c->sem_req == SEM_FAILED || c->sem_done == SEM_FAILED) return 1;
    return 0;
}

static void wait_req(helper_ctx *c)
{
    while (sem_wait(c->sem_req) != 0)
        ; /* retry on EINTR */
}

static void post_done(helper_ctx *c)
{
    sem_post(c->sem_done);
}

static void detach(helper_ctx *c)
{
    if (c->region != NULL) munmap(c->region, (size_t) c->bytes);
    if (c->fd >= 0) close(c->fd);
    if (c->sem_req != SEM_FAILED) sem_close(c->sem_req);
    if (c->sem_done != SEM_FAILED) sem_close(c->sem_done);
}

#endif

/* Resident UMFPACK factorization, real or complex. */
typedef struct {
    void *symbolic;
    void *numeric;
    int is_complex;
} factors_t;

static char *base(void *region, int64_t off)
{
    return (char *) region + off;
}

/* Map a UMFPACK status onto the normalized protocol status. */
static int32_t norm_status(int st)
{
    if (st == UMFPACK_OK) return FSPARSE_ST_OK;
    if (st == UMFPACK_WARNING_singular_matrix) return FSPARSE_ST_SINGULAR;
    return FSPARSE_ST_ERROR;
}

static void free_factors(factors_t *f)
{
    if (f->numeric != NULL) {
        if (f->is_complex) umfpack_zl_free_numeric(&f->numeric);
        else umfpack_dl_free_numeric(&f->numeric);
        f->numeric = NULL;
    }
    if (f->symbolic != NULL) {
        if (f->is_complex) umfpack_zl_free_symbolic(&f->symbolic);
        else umfpack_dl_free_symbolic(&f->symbolic);
        f->symbolic = NULL;
    }
}

static int32_t do_factor_real(fsparse_shm_header *h, void *region,
                              factors_t *f)
{
    int64_t *colptr = (int64_t *) base(region, h->off_colptr);
    int64_t *rowidx = (int64_t *) base(region, h->off_rowidx);
    double *ax = (double *) base(region, h->off_ax);
    int st;

    free_factors(f);
    f->is_complex = 0;
    st = umfpack_dl_symbolic(h->n, h->n, colptr, rowidx, ax, &f->symbolic,
                             NULL, NULL);
    if (st != UMFPACK_OK) return norm_status(st);
    st = umfpack_dl_numeric(colptr, rowidx, ax, f->symbolic, &f->numeric,
                            NULL, NULL);
    return norm_status(st);
}

static int32_t do_factor_complex(fsparse_shm_header *h, void *region,
                                 factors_t *f)
{
    int64_t *colptr = (int64_t *) base(region, h->off_colptr);
    int64_t *rowidx = (int64_t *) base(region, h->off_rowidx);
    double *ax = (double *) base(region, h->off_ax);
    double *az = (double *) base(region, h->off_az);
    int st;

    free_factors(f);
    f->is_complex = 1;
    st = umfpack_zl_symbolic(h->n, h->n, colptr, rowidx, ax, az, &f->symbolic,
                             NULL, NULL);
    if (st != UMFPACK_OK) return norm_status(st);
    st = umfpack_zl_numeric(colptr, rowidx, ax, az, f->symbolic, &f->numeric,
                            NULL, NULL);
    return norm_status(st);
}

static int32_t do_solve_real(fsparse_shm_header *h, void *region, factors_t *f)
{
    int64_t *colptr = (int64_t *) base(region, h->off_colptr);
    int64_t *rowidx = (int64_t *) base(region, h->off_rowidx);
    double *ax = (double *) base(region, h->off_ax);
    double *b = (double *) base(region, h->off_b);
    double *x = (double *) base(region, h->off_x);
    int st;

    if (f->numeric == NULL) return FSPARSE_ST_ERROR;
    st = umfpack_dl_solve(UMFPACK_A, colptr, rowidx, ax, x, b, f->numeric,
                          NULL, NULL);
    return norm_status(st);
}

static int32_t do_solve_complex(fsparse_shm_header *h, void *region,
                                factors_t *f)
{
    int64_t *colptr = (int64_t *) base(region, h->off_colptr);
    int64_t *rowidx = (int64_t *) base(region, h->off_rowidx);
    double *ax = (double *) base(region, h->off_ax);
    double *az = (double *) base(region, h->off_az);
    double *bx = (double *) base(region, h->off_b);
    double *bz = (double *) base(region, h->off_bz);
    double *xx = (double *) base(region, h->off_x);
    double *xz = (double *) base(region, h->off_xz);
    int st;

    if (f->numeric == NULL) return FSPARSE_ST_ERROR;
    st = umfpack_zl_solve(UMFPACK_A, colptr, rowidx, ax, az, xx, xz, bx, bz,
                          f->numeric, NULL, NULL);
    return norm_status(st);
}

static int dispatch(fsparse_shm_header *h, void *region, factors_t *f)
{
    switch (h->opcode) {
    case FSPARSE_OP_FACTOR_REAL:
        h->status = do_factor_real(h, region, f);
        return 0;
    case FSPARSE_OP_FACTOR_COMPLEX:
        h->status = do_factor_complex(h, region, f);
        return 0;
    case FSPARSE_OP_SOLVE_REAL:
        h->status = do_solve_real(h, region, f);
        return 0;
    case FSPARSE_OP_SOLVE_COMPLEX:
        h->status = do_solve_complex(h, region, f);
        return 0;
    case FSPARSE_OP_FREE:
        free_factors(f);
        h->status = FSPARSE_ST_OK;
        return 0;
    case FSPARSE_OP_SHUTDOWN:
        return 1;
    default:
        /* Unknown opcode: complete the doorbell with an error status so the
         * parent unblocks instead of hanging. Only shutdown ends the loop. */
        h->status = FSPARSE_ST_ERROR;
        return 0;
    }
}

int main(int argc, char **argv)
{
    helper_ctx c;
    factors_t f;
    fsparse_shm_header *h;
    int64_t bytes;

    if (argc < 5) {
        fprintf(stderr, "usage: %s shm req done bytes\n", argv[0]);
        return 2;
    }
    memset(&c, 0, sizeof(c));
    memset(&f, 0, sizeof(f));
    bytes = (int64_t) strtoll(argv[4], NULL, 10);
    if (attach(&c, argv, bytes) != 0) {
        detach(&c);
        return 2;
    }
    h = (fsparse_shm_header *) c.region;

    for (;;) {
        wait_req(&c);
        if (dispatch(h, c.region, &f) != 0) break;
        post_done(&c);
    }

    free_factors(&f);
    detach(&c);
    return 0;
}
