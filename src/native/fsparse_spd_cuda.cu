#include "fsparse_spd_cuda.h"
#include "fortsparse_spd_leaves.cuh"
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <limits>
#include <new>

namespace {
constexpr size_t workspace_bytes = 1024 * 1024;
struct Context {
    int device;
    size_t capacity;
    double *arena = nullptr;
    int *info = nullptr;
    void *workspace = nullptr;
    cudaStream_t stream = nullptr;
    cublasHandle_t blas = nullptr;
    int *maps = nullptr;
    size_t map_count = 0;
    bool factor_batch = false;
};

int activate(Context *ctx)
{
    return ctx ? static_cast<int>(cudaSetDevice(ctx->device)) : -1;
}

bool fits(Context *ctx, size_t offset, size_t count)
{
    return ctx && offset <= ctx->capacity && count <= ctx->capacity-offset;
}

bool matrix_fits(Context *ctx, size_t offset, int ld, int rows, int columns)
{
    if (rows < 0 || columns < 0 || ld < (rows > 0 ? rows : 1)) return false;
    if (rows == 0 || columns == 0) return fits(ctx, offset, 0);
    const size_t limit = std::numeric_limits<size_t>::max();
    if (static_cast<size_t>(columns-1) > (limit-rows)/ld) return false;
    return fits(ctx, offset, static_cast<size_t>(columns-1)*ld+rows);
}

/* Small diagonal panel schedule, one cooperative block. Numerical leaves
 * are generated; topology, synchronization and ownership are explicit here.
 */
__global__ void factor_panel(double *a, int n, int lda, int *info, bool retain_info)
{
    if (!retain_info && threadIdx.x == 0) *info = 0;
    __syncthreads();
    if (*info != 0) return;
    for (int j = 0; j < n; ++j) {
        if (threadIdx.x == 0) {
            double pivot = a[j+static_cast<size_t>(j)*lda];
            for (int p = 0; p < j; ++p)
                spd_update_value(pivot, a[j+static_cast<size_t>(p)*lda],
                                 a[j+static_cast<size_t>(p)*lda], &pivot);
            if (!(pivot > 0.0) || !isfinite(pivot)) *info = j+1;
            else spd_root_value(pivot, a+j+static_cast<size_t>(j)*lda);
        }
        __syncthreads();
        if (*info != 0) return;
        for (int i = j+1+threadIdx.x; i < n; i += blockDim.x) {
            double value = a[i+static_cast<size_t>(j)*lda];
            for (int p = 0; p < j; ++p)
                spd_update_value(value, a[i+static_cast<size_t>(p)*lda],
                                 a[j+static_cast<size_t>(p)*lda], &value);
            spd_divide_value(value, a[j+static_cast<size_t>(j)*lda],
                             a+i+static_cast<size_t>(j)*lda);
        }
        __syncthreads();
    }
}

__global__ void scatter_update(double *destination, int ldd, int columns,
                               const double *update, int m, int n,
                               const int *rows, const int *cols, const int *info)
{
    if (*info != 0) return;
    const size_t index = threadIdx.x+static_cast<size_t>(blockIdx.x)*blockDim.x;
    if (index >= static_cast<size_t>(m)*n) return;
    const int i = index % m, j = index / m;
    if (i < j) return;
    const int row = rows[i], column = cols[j];
    if (row <= 0 || row > ldd || column <= 0 || column > columns) return;
    double *entry = destination+(row-1)+static_cast<size_t>(column-1)*ldd;
    spd_update_value(*entry, -1.0, update[index], entry);
}
} // namespace

extern "C" int fsparse_spd_cuda_create(void **context, int device, size_t capacity)
{
    if (!context) return -1;
    *context = nullptr;
    if (capacity == 0 || capacity > std::numeric_limits<size_t>::max()/sizeof(double))
        return -1;
    Context *ctx = new (std::nothrow) Context;
    if (!ctx) return static_cast<int>(cudaErrorMemoryAllocation);
    ctx->device = device;
    ctx->capacity = capacity;
    int result = activate(ctx);
    if (result) {
        delete ctx;
        return result;
    }
    if (!result) result = cudaStreamCreateWithFlags(&ctx->stream, cudaStreamNonBlocking);
    if (!result) result = cudaMalloc(&ctx->arena, capacity*sizeof(double));
    if (!result) result = cudaMalloc(&ctx->info, sizeof(int));
    if (!result) result = cudaMemsetAsync(ctx->info, 0, sizeof(int), ctx->stream);
    if (!result) result = cudaMalloc(&ctx->workspace, workspace_bytes);
    if (!result) {
        cublasStatus_t code = cublasCreate(&ctx->blas);
        if (code == CUBLAS_STATUS_SUCCESS) code = cublasSetStream(ctx->blas, ctx->stream);
        if (code == CUBLAS_STATUS_SUCCESS)
            code = cublasSetWorkspace(ctx->blas, ctx->workspace, workspace_bytes);
        if (code != CUBLAS_STATUS_SUCCESS) result = 1000+code;
    }
    if (result) {
        fsparse_spd_cuda_destroy(ctx);
        return result;
    }
    *context = ctx;
    return 0;
}

extern "C" int fsparse_spd_cuda_destroy(void *context)
{
    auto *ctx = static_cast<Context*>(context);
    if (!ctx) return 0;
    int result = activate(ctx);
    if (result) return result;
    if (ctx->stream) result = cudaStreamSynchronize(ctx->stream);
    if (ctx->blas) cublasDestroy(ctx->blas);
    if (ctx->maps) cudaFree(ctx->maps);
    if (ctx->workspace) cudaFree(ctx->workspace);
    if (ctx->info) cudaFree(ctx->info);
    if (ctx->arena) cudaFree(ctx->arena);
    if (ctx->stream) cudaStreamDestroy(ctx->stream);
    delete ctx;
    return result;
}

extern "C" int fsparse_spd_cuda_upload(void *context, size_t offset,
                                       const double *source, size_t count)
{
    auto *ctx = static_cast<Context*>(context);
    if (!fits(ctx, offset, count) || (count && !source)) return -1;
    int result = activate(ctx);
    if (result || count == 0) return result;
    return cudaMemcpyAsync(ctx->arena+offset, source, count*sizeof(double),
                           cudaMemcpyHostToDevice, ctx->stream);
}

extern "C" int fsparse_spd_cuda_download(void *context, size_t offset,
                                         double *destination, size_t count)
{
    auto *ctx = static_cast<Context*>(context);
    if (!fits(ctx, offset, count) || (count && !destination)) return -1;
    int result = activate(ctx);
    if (result || count == 0) return result;
    return cudaMemcpyAsync(destination, ctx->arena+offset, count*sizeof(double),
                           cudaMemcpyDeviceToHost, ctx->stream);
}

extern "C" int fsparse_spd_cuda_synchronize(void *context)
{
    auto *ctx = static_cast<Context*>(context);
    int result = activate(ctx);
    return result ? result : static_cast<int>(cudaStreamSynchronize(ctx->stream));
}

extern "C" int fsparse_spd_cuda_panel_factor(void *context, size_t offset,
                                             int n, int lda)
{
    auto *ctx = static_cast<Context*>(context);
    if (!matrix_fits(ctx, offset, lda, n, n)) return -1;
    int result = activate(ctx);
    if (result) return result;
    factor_panel<<<1,256,0,ctx->stream>>>(ctx->arena+offset, n, lda, ctx->info,
                                        ctx->factor_batch);
    return cudaGetLastError();
}

extern "C" int fsparse_spd_cuda_factor_info(void *context, int *info)
{
    auto *ctx = static_cast<Context*>(context);
    if (!info) return -1;
    int result = activate(ctx);
    if (result) return result;
    result = cudaMemcpyAsync(info, ctx->info, sizeof(int),
                             cudaMemcpyDeviceToHost, ctx->stream);
    if (!result) result = cudaStreamSynchronize(ctx->stream);
    ctx->factor_batch = false;
    return result;
}

extern "C" int fsparse_spd_cuda_set_maps(void *context, const int *maps, size_t count)
{
    auto *ctx = static_cast<Context*>(context);
    if (!ctx || (count && !maps) || count > std::numeric_limits<size_t>::max()/sizeof(int))
        return -1;
    int result = activate(ctx);
    if (result) return result;
    result = cudaStreamSynchronize(ctx->stream);
    if (result) return result;
    if (ctx->maps) cudaFree(ctx->maps);
    ctx->maps = nullptr;
    ctx->map_count = 0;
    if (!count) return 0;
    result = cudaMalloc(&ctx->maps, count*sizeof(int));
    if (result) return result;
    ctx->map_count = count;
    result = cudaMemcpyAsync(ctx->maps, maps, count*sizeof(int),
                             cudaMemcpyHostToDevice, ctx->stream);
    return result ? result : static_cast<int>(cudaStreamSynchronize(ctx->stream));
}

extern "C" int fsparse_spd_cuda_begin_factor(void *context)
{
    auto *ctx = static_cast<Context*>(context);
    int result = activate(ctx);
    if (result) return result;
    ctx->factor_batch = true;
    return cudaMemsetAsync(ctx->info, 0, sizeof(int), ctx->stream);
}

extern "C" int fsparse_spd_cuda_scatter_update(void *context, size_t update_offset,
                                              int m, int n, size_t destination_offset,
                                              int destination_rows, int destination_columns,
                                              size_t row_map_offset, size_t column_map_offset)
{
    auto *ctx = static_cast<Context*>(context);
    if (!matrix_fits(ctx, update_offset, m, m, n) ||
        !matrix_fits(ctx, destination_offset, destination_rows,
                     destination_rows, destination_columns)) return -1;
    if (row_map_offset > ctx->map_count || column_map_offset > ctx->map_count) return -1;
    if (static_cast<size_t>(m) > ctx->map_count-row_map_offset ||
        static_cast<size_t>(n) > ctx->map_count-column_map_offset) return -1;
    int result = activate(ctx);
    if (result || m == 0 || n == 0) return result;
    const size_t count = static_cast<size_t>(m)*n;
    if ((count+255)/256 > std::numeric_limits<unsigned>::max()) return -1;
    scatter_update<<<static_cast<unsigned>((count+255)/256),256,0,ctx->stream>>>(
        ctx->arena+destination_offset, destination_rows, destination_columns,
        ctx->arena+update_offset, m, n, ctx->maps+row_map_offset,
        ctx->maps+column_map_offset, ctx->info);
    return cudaGetLastError();
}

extern "C" int fsparse_spd_cuda_panel_solve(void *context, size_t l_offset,
                                            int n, int ldl, size_t b_offset,
                                            int m, int ldb)
{
    auto *ctx = static_cast<Context*>(context);
    if (!matrix_fits(ctx, l_offset, ldl, n, n) ||
        !matrix_fits(ctx, b_offset, ldb, m, n)) return -1;
    int result = activate(ctx);
    if (result || m == 0 || n == 0) return result;
    const double one = 1.0;
    const cublasStatus_t code = cublasDtrsm(ctx->blas, CUBLAS_SIDE_RIGHT,
        CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT,
        m, n, &one, ctx->arena+l_offset, ldl, ctx->arena+b_offset, ldb);
    return code == CUBLAS_STATUS_SUCCESS ? 0 : 1000+code;
}

static int schur_gemm(void *context, size_t c_offset, int ldc,
                       size_t left_offset, int ldl, size_t right_offset, int ldr,
                       int m, int n, int k, double beta)
{
    auto *ctx = static_cast<Context*>(context);
    if (!matrix_fits(ctx, c_offset, ldc, m, n) ||
        !matrix_fits(ctx, left_offset, ldl, m, k) ||
        !matrix_fits(ctx, right_offset, ldr, n, k)) return -1;
    int result = activate(ctx);
    if (result || m == 0 || n == 0) return result;
    const double minus_one = -1.0;
    const cublasStatus_t code = cublasDgemm(ctx->blas, CUBLAS_OP_N, CUBLAS_OP_T,
        m, n, k, &minus_one, ctx->arena+left_offset, ldl,
        ctx->arena+right_offset, ldr, &beta, ctx->arena+c_offset, ldc);
    return code == CUBLAS_STATUS_SUCCESS ? 0 : 1000+code;
}

extern "C" int fsparse_spd_cuda_schur_update(void *context, size_t c_offset,
                                             int ldc, size_t left_offset, int ldl,
                                             size_t right_offset, int ldr,
                                             int m, int n, int k)
{
    return schur_gemm(context, c_offset, ldc, left_offset, ldl,
                       right_offset, ldr, m, n, k, 1.0);
}

extern "C" int fsparse_spd_cuda_schur_product(void *context, size_t c_offset,
                                              int ldc, size_t left_offset, int ldl,
                                              size_t right_offset, int ldr,
                                              int m, int n, int k)
{
    return schur_gemm(context, c_offset, ldc, left_offset, ldl,
                       right_offset, ldr, m, n, k, 0.0);
}

extern "C" size_t fsparse_spd_cuda_allocated_bytes(void *context)
{
    auto *ctx = static_cast<Context*>(context);
    return ctx ? ctx->capacity*sizeof(double)+sizeof(int)+workspace_bytes
                 +ctx->map_count*sizeof(int) : 0;
}
