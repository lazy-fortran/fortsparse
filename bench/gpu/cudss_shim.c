/* Project-owned benchmark adapter. NVIDIA cuDSS remains an external SDK.
 * The required SDK disclaimer and Government End Users notice are reproduced
 * below and in CUDSS-NOTICE.md; no NVIDIA implementation is copied here.
 * NOTWITHSTANDING ANY TERMS OR CONDITIONS TO THE CONTRARY IN THE
 * LICENSE AGREEMENT, NVIDIA MAKES NO REPRESENTATION ABOUT THE
 * SUITABILITY OF THESE LICENSED DELIVERABLES FOR ANY PURPOSE. IT IS
 * PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY OF ANY KIND.
 * NVIDIA DISCLAIMS ALL WARRANTIES WITH REGARD TO THESE LICENSED
 * DELIVERABLES, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY,
 * NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
 * NOTWITHSTANDING ANY TERMS OR CONDITIONS TO THE CONTRARY IN THE
 * LICENSE AGREEMENT, IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY
 * SPECIAL, INDIRECT, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, OR ANY
 * DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
 * WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
 * ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
 * OF THESE LICENSED DELIVERABLES.
 * U.S. Government End Users. These Licensed Deliverables are a
 * "commercial item" as that term is defined at 48 C.F.R. 2.101 (OCT
 * 1995), consisting of "commercial computer software" and "commercial
 * computer software documentation" as such terms are used in 48
 * C.F.R. 12.212 (SEPT 1995) and is provided to the U.S. Government
 * only as a commercial end item. Consistent with 48 C.F.R.12.212 and
 * 48 C.F.R. 227.7202-1 through 227.7202-4 (JUNE 1995), all
 * U.S. Government End Users acquire the Licensed Deliverables with
 * only those rights set forth herein.
 */
#define _POSIX_C_SOURCE 200809L
#include <cuda_runtime_api.h>
#include <cudss.h>
#include "process_device_memory.h"
#include <stdint.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

enum { RUNTIME, CONVERSION, ANALYSIS, NUMERIC, REFACTOR, SOLVE, H2D, D2H,
       SYNC, CLEANUP, FACTOR, COLD_FACTOR, SAMPLING, NTIMES };
enum { REQUESTED_PEAK, USED_PEAK, RESERVED_PEAK, WRAPPER_BYTES, PROCESS_SAMPLE,
       TOTAL_PEAK, PINNED_BYTES, FACTOR_NNZ, ACTUAL_DEVICE, ACTIVE_THREADS,
       ORDER_HASH, VERSION, RELEASE_THRESHOLD, DEVICE_ESTIMATE, HOST_ESTIMATE,
       POOL_LIMIT, SAMPLE_COUNT, NCOUNTS };
typedef struct {
    cudaStream_t stream;
    cudaMemPool_t pool;
    cudssHandle_t handle;
    cudssConfig_t config;
    cudssData_t data;
    cudssMatrix_t a, b, x;
    int n, nnz, nrhs, *ptr, *row, *source, *dptr, *drow;
    double *values, *dvalues, *db, *dx;
    _Atomic uint64_t requested, high;
    double times[NTIMES];
    int64_t counts[NCOUNTS];
} context;

static double seconds(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + 1.e-9 * (double)t.tv_nsec;
}
static int allocate(void *opaque, void **pointer, size_t bytes, cudaStream_t stream) {
    context *c = opaque;
    cudaError_t code = cudaMallocFromPoolAsync(pointer, bytes, c->pool, stream);
    if (code != cudaSuccess) return (int)code;
    uint64_t live = atomic_fetch_add(&c->requested, bytes) + bytes;
    uint64_t high = atomic_load(&c->high);
    while (high < live && !atomic_compare_exchange_weak(&c->high, &high, live)) {}
    return 0;
}
static int release(void *opaque, void *pointer, size_t bytes, cudaStream_t stream) {
    context *c = opaque;
    cudaError_t code = cudaFreeAsync(pointer, stream);
    if (code == cudaSuccess) atomic_fetch_sub(&c->requested, bytes);
    return (int)code;
}
static void snapshot(context *c) {
    double start = seconds();
    uint64_t value;
    if (c->pool) {
        if (cudaMemPoolGetAttribute(c->pool, cudaMemPoolAttrUsedMemHigh, &value) == cudaSuccess)
            if ((int64_t)value > c->counts[USED_PEAK]) c->counts[USED_PEAK] = (int64_t)value;
        if (cudaMemPoolGetAttribute(c->pool, cudaMemPoolAttrReservedMemHigh, &value) == cudaSuccess)
            if ((int64_t)value > c->counts[RESERVED_PEAK]) c->counts[RESERVED_PEAK] = (int64_t)value;
    }
    c->counts[REQUESTED_PEAK] = (int64_t)atomic_load(&c->high);
    int64_t process = fsbench_gpu_process_memory((int)c->counts[ACTUAL_DEVICE]);
    if (c->counts[SAMPLE_COUNT] < 0) c->counts[SAMPLE_COUNT] = 0;
    if (process >= 0) ++c->counts[SAMPLE_COUNT];
    if (process > c->counts[PROCESS_SAMPLE]) c->counts[PROCESS_SAMPLE] = process;
    c->times[SAMPLING] += seconds() - start;
}
static int synchronize(context *c) {
    double start = seconds();
    cudaError_t code = cudaStreamSynchronize(c->stream);
    c->times[SYNC] += seconds() - start;
    return (int)code;
}
static int copy(context *c, void *to, const void *from, size_t bytes, enum cudaMemcpyKind kind) {
    double start = seconds();
    int code = (int)cudaMemcpyAsync(to, from, bytes, kind, c->stream);
    if (!code) code = synchronize(c);
    c->times[kind == cudaMemcpyHostToDevice ? H2D : D2H] += seconds() - start;
    return code;
}
static int execute(context *c, int phase, int timer) {
    double start = seconds();
    int code = (int)cudssExecute(c->handle, phase, c->config, c->data, c->a, c->x, c->b);
    if (!code) code = synchronize(c);
    c->times[timer] += seconds() - start;
    if (!code) {
        int info = 0;
        size_t written;
        code = (int)cudssDataGet(c->handle, c->data, CUDSS_DATA_INFO,
                               &info, sizeof(info), &written);
        if (!code && info) code = 1000 + info;
    }
    return code;
}
static void destroy(context *c) {
    double start = seconds();
    if (c->a) cudssMatrixDestroy(c->a);
    if (c->b) cudssMatrixDestroy(c->b);
    if (c->x) cudssMatrixDestroy(c->x);
    if (c->data) cudssDataDestroy(c->handle, c->data);
    if (c->config) cudssConfigDestroy(c->config);
    if (c->handle) cudssDestroy(c->handle);
    size_t dense = (size_t)c->n * c->nrhs * sizeof(double);
    if (c->dptr) release(c, c->dptr, ((size_t)c->n + 1) * sizeof(int), c->stream);
    if (c->drow) release(c, c->drow, (size_t)c->nnz * sizeof(int), c->stream);
    if (c->dvalues) release(c, c->dvalues, (size_t)c->nnz * sizeof(double), c->stream);
    if (c->db) release(c, c->db, dense, c->stream);
    if (c->dx) release(c, c->dx, dense, c->stream);
    if (c->stream) synchronize(c);
    snapshot(c);
    if (c->pool) cudaMemPoolDestroy(c->pool);
    if (c->stream) cudaStreamDestroy(c->stream);
    free(c->ptr); free(c->row); free(c->values); free(c->source);
    c->times[CLEANUP] += seconds() - start;
}

#define CUDA(call) do { int result = (int)(call); if (result) { \
    fprintf(stderr, "%s: CUDA %d\n", #call, result); return 2000 + result; } } while (0)
#define CUDSS(call) do { int result = (int)(call); if (result) { \
    fprintf(stderr, "%s: cuDSS %d\n", #call, result); return result; } } while (0)

static int initialize(context *c, int n, const int *ptr, const int *row,
                      const double *values, int nrhs, int threads, int device,
                      int ordering, const int *permutation) {
    double start = seconds();
    c->n = n; c->nrhs = nrhs;
    c->ptr = malloc(((size_t)n + 1) * sizeof(int));
    c->source = malloc((size_t)(ptr[n] - 1) * sizeof(int));
    c->row = malloc((size_t)(ptr[n] - 1) * sizeof(int));
    c->values = malloc((size_t)(ptr[n] - 1) * sizeof(double));
    if (!c->ptr || !c->source || !c->row || !c->values) return 6;
    /* Symmetric lower CSC is upper CSR; preserve exact original values. */
    for (int j = 0; j < n; ++j) {
        c->ptr[j] = c->nnz;
        for (int p = ptr[j] - 1; p < ptr[j + 1] - 1; ++p) {
            if (row[p] < j + 1) continue;
            c->row[c->nnz] = row[p] - 1;
            c->source[c->nnz] = p;
            c->values[c->nnz++] = values[p];
        }
    }
    c->ptr[n] = c->nnz;
    c->times[CONVERSION] = seconds() - start;
    start = seconds();
    CUDA(cudaSetDevice(device));
    CUDA(cudaFree(NULL));
    int actual;
    CUDA(cudaGetDevice(&actual));
    c->counts[ACTUAL_DEVICE] = actual;
    CUDA(cudaStreamCreateWithFlags(&c->stream, cudaStreamNonBlocking));
    struct cudaMemPoolProps properties = {0};
    properties.allocType = cudaMemAllocationTypePinned;
    properties.location.type = cudaMemLocationTypeDevice;
    properties.location.id = device;
    properties.maxSize = 256u * 1024u * 1024u;
    const char *limit = getenv("FSBENCH_CUDSS_POOL_LIMIT_BYTES");
    if (limit) {
        char *end;
        unsigned long long parsed = strtoull(limit, &end, 10);
        if (!parsed || *end || parsed > SIZE_MAX) return 2;
        properties.maxSize = (size_t)parsed;
    }
    c->counts[POOL_LIMIT] = (int64_t)properties.maxSize;
    CUDA(cudaMemPoolCreate(&c->pool, &properties));
    uint64_t threshold = UINT64_MAX;
    CUDA(cudaMemPoolSetAttribute(c->pool, cudaMemPoolAttrReleaseThreshold, &threshold));
    c->counts[RELEASE_THRESHOLD] = -1; /* UINT64_MAX, documented as retain-all. */
    CUDSS(cudssCreate(&c->handle));
    CUDSS(cudssSetStream(c->handle, c->stream));
    cudssDeviceMemHandler_t handler = {0};
    handler.ctx = c; handler.device_alloc = allocate; handler.device_free = release;
    strcpy(handler.name, "fortsparse-benchmark-private-pool");
    CUDSS(cudssSetDeviceMemHandler(c->handle, &handler));
    CUDSS(cudssSetThreadingLayer(c->handle, FSBENCH_CUDSS_THREADING_LAYER));
    CUDSS(cudssConfigCreate(&c->config));
    CUDSS(cudssDataCreate(c->handle, &c->data));
    CUDSS(cudssConfigSet(c->config, CUDSS_CONFIG_HOST_NTHREADS, &threads, sizeof(threads)));
    int disabled = 0;
    CUDSS(cudssConfigSet(c->config, CUDSS_CONFIG_HYBRID_EXECUTE_MODE, &disabled, sizeof(disabled)));
    CUDSS(cudssConfigSet(c->config, CUDSS_CONFIG_HYBRID_MEMORY_MODE, &disabled, sizeof(disabled)));
    CUDSS(cudssConfigSet(c->config, CUDSS_CONFIG_IR_N_STEPS, &disabled, sizeof(disabled)));
    cudssReorderingAlg_t algorithm = ordering == 1 ? CUDSS_REORDERING_ALG_AMD :
        ordering == 2 ? CUDSS_REORDERING_ALG_NESTED_DISSECTION :
        ordering == 3 ? CUDSS_REORDERING_ALG_NONE : CUDSS_REORDERING_ALG_DEFAULT;
    CUDSS(cudssConfigSet(c->config, CUDSS_CONFIG_REORDERING_ALG, &algorithm, sizeof(algorithm)));
    if (ordering == 3)
        CUDSS(cudssDataSet(c->handle, c->data, CUDSS_DATA_USER_PERM,
                          permutation, (size_t)n * sizeof(int)));
    CUDA(allocate(c, (void **)&c->dptr, ((size_t)n + 1) * sizeof(int), c->stream));
    CUDA(allocate(c, (void **)&c->drow, (size_t)c->nnz * sizeof(int), c->stream));
    CUDA(allocate(c, (void **)&c->dvalues, (size_t)c->nnz * sizeof(double), c->stream));
    CUDA(allocate(c, (void **)&c->db, (size_t)n * nrhs * sizeof(double), c->stream));
    CUDA(allocate(c, (void **)&c->dx, (size_t)n * nrhs * sizeof(double), c->stream));
    c->counts[WRAPPER_BYTES] = ((int64_t)n + 1 + c->nnz) * sizeof(int) +
        ((int64_t)c->nnz + 2 * (int64_t)n * nrhs) * sizeof(double);
    c->times[RUNTIME] = seconds() - start;
    CUDA(copy(c, c->dptr, c->ptr, ((size_t)n + 1) * sizeof(int), cudaMemcpyHostToDevice));
    CUDA(copy(c, c->drow, c->row, (size_t)c->nnz * sizeof(int), cudaMemcpyHostToDevice));
    CUDA(copy(c, c->dvalues, c->values, (size_t)c->nnz * sizeof(double), cudaMemcpyHostToDevice));
    CUDSS(cudssMatrixCreateCsr(&c->a, n, n, c->nnz, c->dptr, NULL, c->drow,
        c->dvalues, CUDSS_R_32I, CUDSS_R_32I, CUDSS_R_64F, CUDSS_MTYPE_SPD,
        CUDSS_MVIEW_UPPER, CUDSS_BASE_ZERO));
    CUDSS(cudssMatrixCreateDn(&c->b, n, nrhs, n, c->db, CUDSS_R_64F, CUDSS_LAYOUT_COL_MAJOR));
    CUDSS(cudssMatrixCreateDn(&c->x, n, nrhs, n, c->dx, CUDSS_R_64F, CUDSS_LAYOUT_COL_MAJOR));
    int result = execute(c, CUDSS_PHASE_ANALYSIS, ANALYSIS);
    if (result) return result;
    int *actual_permutation = malloc((size_t)n * sizeof(int));
    if (!actual_permutation) return 6;
    size_t written;
    result = (int)cudssDataGet(c->handle, c->data, CUDSS_DATA_PERM_REORDER_ROW,
                              actual_permutation, (size_t)n * sizeof(int), &written);
    if (!result) {
        uint64_t hash = UINT64_C(14695981039346656037);
        for (int j = 0; j < n; ++j) {
            if (ordering == 3 && actual_permutation[j] != permutation[j]) result = 2;
            hash = (hash ^ (uint32_t)actual_permutation[j]) * UINT64_C(1099511628211);
        }
        c->counts[ORDER_HASH] = (int64_t)hash;
    }
    free(actual_permutation);
    if (result) return result;
    int64_t estimate[16];
    if (cudssDataGet(c->handle, c->data, CUDSS_DATA_MEMORY_ESTIMATES,
                     estimate, sizeof(estimate), &written) == CUDSS_STATUS_SUCCESS) {
        c->counts[DEVICE_ESTIMATE] = estimate[1]; c->counts[HOST_ESTIMATE] = estimate[3];
    }
    if (cudssDataGet(c->handle, c->data, CUDSS_DATA_LU_NNZ, &c->counts[FACTOR_NNZ],
                     sizeof(int64_t), &written) != CUDSS_STATUS_SUCCESS)
        c->counts[FACTOR_NNZ] = -1;
    return execute(c, CUDSS_PHASE_FACTORIZATION, NUMERIC);
}

int fsbench_cudss_run(int n, const int *ptr, const int *row, const double *values,
                     const double *updated, const double *rhs, double *solution,
                     int nrhs, int threads, int device, int workload, int resident,
                     int ordering, const int *permutation, int repeats,
                     double times[NTIMES], int64_t counts[NCOUNTS]) {
    context c = {0};
    for (int i = 0; i < NCOUNTS; ++i) c.counts[i] = -1;
    c.counts[PINNED_BYTES] = -1; /* Internal vendor host registrations unknown. */
    c.counts[VERSION] = CUDSS_VERSION;
    c.times[REFACTOR] = -1.;
    double warm_sampling = 0.;
    int64_t warm_counts[NCOUNTS];
    for (int i = 0; i < NCOUNTS; ++i) warm_counts[i] = -1;
    double start = seconds();
    int result = initialize(&c, n, ptr, row, values, nrhs, threads, device, ordering, permutation);
    c.times[FACTOR] = resident ? c.times[ANALYSIS] + c.times[NUMERIC] : seconds() - start;
    c.times[COLD_FACTOR] = c.times[FACTOR];
    snapshot(&c);
    if (!result && workload == 2) { /* Warm runtime, fresh factors/analysis. */
        double cold = c.times[COLD_FACTOR];
        destroy(&c);
        warm_sampling = c.times[SAMPLING];
        memcpy(warm_counts, c.counts, sizeof(warm_counts));
        memset(&c, 0, sizeof(c));
        for (int i = 0; i < NCOUNTS; ++i) c.counts[i] = -1;
        c.counts[VERSION] = CUDSS_VERSION;
        c.times[REFACTOR] = -1.;
        start = seconds();
        result = initialize(&c, n, ptr, row, values, nrhs, threads, device, ordering, permutation);
        c.times[FACTOR] = resident ? c.times[ANALYSIS] + c.times[NUMERIC] : seconds() - start;
        c.times[COLD_FACTOR] = cold;
        snapshot(&c);
    }
    if (!result && workload == 3) {
        start = seconds();
        for (int p = 0; p < c.nnz; ++p) c.values[p] = updated[c.source[p]];
        result = copy(&c, c.dvalues, c.values, (size_t)c.nnz * sizeof(double), cudaMemcpyHostToDevice);
        c.times[REFACTOR] = 0.;
        if (!result) result = execute(&c, CUDSS_PHASE_REFACTORIZATION, REFACTOR);
        if (!resident) c.times[REFACTOR] = seconds() - start;
        snapshot(&c);
    }
    size_t bytes = (size_t)n * nrhs * sizeof(double);
    for (int repeat = 0; !result && repeat < repeats; ++repeat) {
        double before_solve = c.times[SOLVE];
        start = seconds();
        result = copy(&c, c.db, rhs + (size_t)repeat * n * nrhs, bytes, cudaMemcpyHostToDevice);
        if (!result) result = execute(&c, CUDSS_PHASE_SOLVE, SOLVE);
        if (!result) result = copy(&c, solution + (size_t)repeat * n * nrhs, c.dx,
                                   bytes, cudaMemcpyDeviceToHost);
        if (!resident) c.times[SOLVE] = before_solve + seconds() - start;
        snapshot(&c);
    }
    destroy(&c);
    c.times[SAMPLING] += warm_sampling;
    for (int i = REQUESTED_PEAK; i <= PROCESS_SAMPLE; ++i)
        if (warm_counts[i] > c.counts[i]) c.counts[i] = warm_counts[i];
    if (warm_counts[SAMPLE_COUNT] > 0) c.counts[SAMPLE_COUNT] += warm_counts[SAMPLE_COUNT];
    memcpy(times, c.times, sizeof(c.times));
    memcpy(counts, c.counts, sizeof(c.counts));
    return result;
}
