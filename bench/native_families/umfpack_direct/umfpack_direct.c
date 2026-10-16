/* SPDX-License-Identifier: GPL-2.0-or-later
 * Benchmark-only UMFPACK client. Never compiled or linked into libfortsparse.
 * Uses the installed official I4 di/zi API and unchanged native defaults.
 */
#include <umfpack.h>
#include <complex.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

_Static_assert(sizeof(int) == sizeof(int32_t) && INT_MAX == INT32_MAX,
               "UMFPACK di/zi comparator requires 32-bit C int");
_Static_assert(sizeof(double complex) == 2 * sizeof(double),
               "UMFPACK packed complex requires adjacent real/imaginary doubles");

typedef struct {
    void *symbolic, *numeric;
    int is_complex;
    int32_t n;
    double control[UMFPACK_CONTROL], info[UMFPACK_INFO];
    int64_t factor_nnz;
    double factor_peak_mb;
} benchmark_umf_t;

void *fsbench_umf_create(int is_complex)
{
    benchmark_umf_t *ctx = calloc(1, sizeof(*ctx));
    if (!ctx) return NULL;
    ctx->is_complex = is_complex;
    ctx->factor_nnz = -1;
    ctx->factor_peak_mb = -1;
    if (is_complex) umfpack_zi_defaults(ctx->control);
    else umfpack_di_defaults(ctx->control);
    return ctx;
}

int fsbench_umf_symbolic(void *owner, int32_t n, const int32_t *ap,
                        const int32_t *ai, const void *ax)
{
    benchmark_umf_t *ctx = owner;
    ctx->n = n;
    if (ctx->is_complex)
        return umfpack_zi_symbolic(n, n, ap, ai, ax, NULL, &ctx->symbolic,
                                  ctx->control, ctx->info);
    return umfpack_di_symbolic(n, n, ap, ai, ax, &ctx->symbolic,
                              ctx->control, ctx->info);
}

int fsbench_umf_numeric(void *owner, const int32_t *ap, const int32_t *ai,
                       const void *ax)
{
    benchmark_umf_t *ctx = owner;
    int status;
    /* Reuse Symbolic but drop the previous Numeric before recomputing it.
     * This excludes no live overlap: the old numeric allocation is gone. */
    if (ctx->is_complex) {
        umfpack_zi_free_numeric(&ctx->numeric);
        status = umfpack_zi_numeric(ap, ai, ax, NULL, ctx->symbolic,
                                   &ctx->numeric, ctx->control, ctx->info);
    } else {
        umfpack_di_free_numeric(&ctx->numeric);
        status = umfpack_di_numeric(ap, ai, ax, ctx->symbolic,
                                   &ctx->numeric, ctx->control, ctx->info);
    }
    if (status == UMFPACK_OK) {
        double fill = ctx->info[UMFPACK_LNZ] + ctx->info[UMFPACK_UNZ] - ctx->n;
        double units = ctx->info[UMFPACK_PEAK_MEMORY];
        double unit_bytes = ctx->info[UMFPACK_SIZE_OF_UNIT];
        ctx->factor_nnz = -1;
        if (isfinite(fill) && fill >= 0 && fill < (double) INT64_MAX)
            ctx->factor_nnz = (int64_t) fill;
        /* Info memory counts are Units, not bytes. Preserve these vendor
         * diagnostics before solve overwrites Info. This is not process RSS
         * or exact allocator-owned payload; the runner uses postexec VmHWM. */
        if (isfinite(units) && isfinite(unit_bytes) && units >= 0 && unit_bytes > 0) {
            double mb = units * unit_bytes / 1.e6;
            if (isfinite(mb) && mb > ctx->factor_peak_mb) ctx->factor_peak_mb = mb;
        }
    }
    return status;
}

int fsbench_umf_solve(void *owner, int operation, const int32_t *ap,
                     const int32_t *ai, const void *ax, const void *b, void *x)
{
    benchmark_umf_t *ctx = owner;
    int system;
    if (operation < 0 || operation > 2) return UMFPACK_ERROR_invalid_system;
    system = UMFPACK_A;
    if (operation != 0) system = UMFPACK_At;
    if (ctx->is_complex && operation == 1) system = UMFPACK_Aat;
    if (ctx->is_complex)
        return umfpack_zi_solve(system, ap, ai, ax, NULL, x, NULL, b, NULL,
                               ctx->numeric, ctx->control, ctx->info);
    return umfpack_di_solve(system, ap, ai, ax, x, b, ctx->numeric,
                           ctx->control, ctx->info);
}

void fsbench_umf_stats(void *owner, int64_t *factor_nnz, double *factor_peak_mb)
{
    benchmark_umf_t *ctx = owner;
    *factor_nnz = ctx->factor_nnz;
    *factor_peak_mb = ctx->factor_peak_mb;
}

void fsbench_umf_destroy(void *owner)
{
    benchmark_umf_t *ctx = owner;
    if (!ctx) return;
    if (ctx->is_complex) {
        umfpack_zi_free_numeric(&ctx->numeric);
        umfpack_zi_free_symbolic(&ctx->symbolic);
    } else {
        umfpack_di_free_numeric(&ctx->numeric);
        umfpack_di_free_symbolic(&ctx->symbolic);
    }
    free(ctx);
}
