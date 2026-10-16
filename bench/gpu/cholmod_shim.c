/* Standalone benchmark adapter; never linked into the FortSparse library.
 * The executable uses CHOLMOD GPL-2.0+ Supernodal/GPU modules.
 */
#include <cholmod.h>
#include <cuda_runtime_api.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double wall_seconds(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

struct cholmod_bench {
    cholmod_common common;
    cholmod_sparse *matrix;
    cholmod_factor *factor;
    cholmod_dense *solution, *work_y, *work_e;
    int64_t *permutation;
    int *source_positions;
    size_t lower_nnz, wrapper_host_bytes, device_peak, pinned_peak;
    size_t gpu_calls[4];
    int device;
};

void fsbench_cholmod_free(void *opaque);

void *fsbench_cholmod_create(int n, int nnz, const int *ptr,
    const int *row, const double *values, int threads, int device,
    int ordering, const int64_t *permutation, int64_t pool_bytes,
    double *setup_times, int *status)
{
    double begin = wall_seconds();
    *status = 0;
    setup_times[0] = setup_times[1] = -1;
    struct cholmod_bench *c = calloc(1, sizeof(*c));
    if (!c) { *status = -100; return NULL; }
    if (!cholmod_l_start(&c->common)) {
        free(c); *status = -101; return NULL;
    }
    c->device = device;
    c->common.print = 0;
    c->common.nthreads_max = threads;
    c->common.supernodal = CHOLMOD_SUPERNODAL;
    c->common.useGPU = 1;
    c->common.maxGpuMemBytes = (size_t)pool_bytes;
    c->common.maxGpuMemFraction = 0;
    if (cudaSetDevice(device) != cudaSuccess || cudaFree(NULL) != cudaSuccess) {
        *status = -102; fsbench_cholmod_free(c); return NULL;
    }
    setup_times[0] = wall_seconds() - begin;
    begin = wall_seconds();
    if (ordering != 0) {
        c->common.nmethods = 1;
        c->common.method[0].ordering = ordering == 1 ? CHOLMOD_AMD :
            ordering == 2 ? CHOLMOD_METIS : CHOLMOD_GIVEN;
        if (ordering == 3) c->common.postorder = 0;
    }
    for (int j = 0; j < n; ++j)
        for (int p = ptr[j] - 1; p < ptr[j+1] - 1; ++p)
            if (row[p] - 1 >= j) ++c->lower_nnz;
    c->matrix = cholmod_l_allocate_sparse(n, n, c->lower_nnz, 1, 1, -1,
        CHOLMOD_REAL, &c->common);
    c->source_positions = malloc(c->lower_nnz * sizeof(int));
    c->permutation = ordering == 3 ? malloc((size_t)n * sizeof(int64_t)) : NULL;
    if (!c->matrix || !c->source_positions ||
        (ordering == 3 && !c->permutation)) {
        *status = -103; fsbench_cholmod_free(c); return NULL;
    }
    if (c->permutation)
        memcpy(c->permutation, permutation, (size_t)n * sizeof(int64_t));
    c->wrapper_host_bytes = sizeof(*c) + c->lower_nnz * sizeof(int) +
        (c->permutation ? (size_t)n * sizeof(int64_t) : 0);
    int64_t *ap = c->matrix->p, *ai = c->matrix->i;
    double *ax = c->matrix->x;
    size_t k = 0;
    for (int j = 0; j < n; ++j) {
        ap[j] = k;
        for (int p = ptr[j] - 1; p < ptr[j+1] - 1; ++p) {
            if (row[p] - 1 < j) continue;
            ai[k] = row[p] - 1;
            ax[k] = values[p];
            c->source_positions[k++] = p;
        }
    }
    ap[n] = k;
    setup_times[1] = wall_seconds() - begin;
    (void)nnz;
    return c;
}

void fsbench_cholmod_release_factor(void *opaque)
{
    struct cholmod_bench *c = opaque;
    cholmod_l_free_factor(&c->factor, &c->common);
}

int fsbench_cholmod_analyze(void *opaque)
{
    struct cholmod_bench *c = opaque;
    cholmod_l_free_factor(&c->factor, &c->common);
    c->factor = c->permutation ? cholmod_l_analyze_p(c->matrix,
        c->permutation, NULL, 0, &c->common) :
        cholmod_l_analyze(c->matrix, &c->common);
    return c->factor ? c->common.status : -104;
}

int fsbench_cholmod_factor(void *opaque, const double *values)
{
    struct cholmod_bench *c = opaque;
    double *ax = c->matrix->x;
    for (size_t k = 0; k < c->lower_nnz; ++k)
        ax[k] = values[c->source_positions[k]];
    int ok = cholmod_l_factorize(c->matrix, c->factor, &c->common);
    if (cudaDeviceSynchronize() != cudaSuccess) return -105;
    c->gpu_calls[0] += c->common.cholmod_gpu_gemm_calls;
    c->gpu_calls[1] += c->common.cholmod_gpu_syrk_calls;
    c->gpu_calls[2] += c->common.cholmod_gpu_trsm_calls;
    c->gpu_calls[3] += c->common.cholmod_gpu_potrf_calls;
    if (c->common.dev_mempool_size > c->device_peak)
        c->device_peak = c->common.dev_mempool_size;
    if (c->common.host_pinned_mempool_size > c->pinned_peak)
        c->pinned_peak = c->common.host_pinned_mempool_size;
    if (!ok || c->factor->minor != c->matrix->nrow) return -106;
    return c->common.status;
}

int fsbench_cholmod_solve(void *opaque, int nrhs, const double *rhs, double *x)
{
    struct cholmod_bench *c = opaque;
    size_t n = c->matrix->nrow;
    cholmod_dense b = {0};
    b.nrow = n; b.ncol = nrhs; b.nzmax = n * nrhs; b.d = n;
    b.x = (void *)rhs; b.xtype = CHOLMOD_REAL; b.dtype = CHOLMOD_DOUBLE;
    int ok = cholmod_l_solve2(CHOLMOD_A, c->factor, &b, NULL, &c->solution,
        NULL, &c->work_y, &c->work_e, &c->common);
    if (!ok || !c->solution) return -107;
    memcpy(x, c->solution->x, n * (size_t)nrhs * sizeof(double));
    return c->common.status;
}

/* stats: CHOLMOD tracked host peak; wrapper host; device pool peak; pinned
 * pool peak; four GPU BLAS call counts; numerical factor stored entries;
 * device (-1 without positive counter evidence); factor permutation. */
void fsbench_cholmod_stats(void *opaque, double *stats, int64_t *permutation)
{
    struct cholmod_bench *c = opaque;
    stats[0] = c->common.memory_usage;
    stats[1] = c->wrapper_host_bytes;
    stats[2] = c->device_peak;
    stats[3] = c->pinned_peak;
    size_t calls = 0;
    for (int j = 0; j < 4; ++j) { stats[4+j] = c->gpu_calls[j]; calls += c->gpu_calls[j]; }
    stats[8] = c->factor->is_super ? c->factor->xsize : c->factor->nzmax;
    stats[9] = calls ? c->device : -1;
    memcpy(permutation, c->factor->Perm, c->matrix->nrow * sizeof(int64_t));
}

void fsbench_cholmod_free(void *opaque)
{
    struct cholmod_bench *c = opaque;
    if (!c) return;
    cholmod_l_free_dense(&c->solution, &c->common);
    cholmod_l_free_dense(&c->work_y, &c->common);
    cholmod_l_free_dense(&c->work_e, &c->common);
    cholmod_l_free_factor(&c->factor, &c->common);
    cholmod_l_free_sparse(&c->matrix, &c->common);
    cholmod_l_finish(&c->common);
    free(c->source_positions);
    free(c->permutation);
    free(c);
}

int64_t fsbench_cholmod_ordering_hash(void *opaque)
{
    struct cholmod_bench *c = opaque;
    const int64_t *p = c->factor->Perm;
    /* FNV-1a over each zero-based new->old int64 in little-endian byte order. */
    uint64_t hash = UINT64_C(14695981039346656037);
    for (size_t j = 0; j < c->matrix->nrow; ++j) {
        uint64_t value = (uint64_t)p[j];
        for (int byte = 0; byte < 8; ++byte) {
            hash ^= (value >> (8 * byte)) & 255;
            hash *= UINT64_C(1099511628211);
        }
    }
    int64_t result;
    memcpy(&result, &hash, sizeof(result));
    return result;
}

void fsbench_cholmod_version(int *version)
{
    cholmod_version(version);
}
