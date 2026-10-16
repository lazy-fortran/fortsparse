#include "fsparse_spd_cuda.h"
#include "fortsparse_spd_leaves.cuh"
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <limits>
#include <new>
#include <vector>
#include <utility>
#include <climits>

namespace {
constexpr size_t workspace_bytes = 1024 * 1024;
struct DeviceTask {
    fsparse_spd_cuda_schur_task task;
    size_t block_end;
};
struct Batch {
    size_t first, count, block_first, blocks;
};
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
    DeviceTask *tasks = nullptr;
    size_t task_count = 0;
    std::vector<Batch> batches;
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
    spd_factor_panel_schedule(a, n, n, lda, info, retain_info);
}

__global__ void factor_front(double *a, int width, int rows, int lda,
                             int *info, bool retain_info)
{
    spd_factor_panel_schedule(a, width, rows, lda, info, retain_info);
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

__global__ void schur_scatter(const double *left, int ldl,
                              const double *right, int ldr, int m, int n, int k,
                              double *destination, int ldd, int columns,
                              const int *rows, const int *cols, const int *info)
{
    if (*info != 0) return;
    const size_t index = threadIdx.x+static_cast<size_t>(blockIdx.x)*blockDim.x;
    spd_schur_scatter_entry(left, ldl, right, ldr, m, n, k,
                            destination, ldd, columns, rows, cols, index);
}
// Each block finds its retained task in a prefix schedule. A large task
// keeps the same 256-entry tiling as the single-task kernel. Arithmetic and
// mapped lower-triangle semantics come from the retained generated schedule.
__global__ void schur_batch(double *arena, const int *maps,
                            const DeviceTask *tasks, Batch batch, const int *info)
{
    if (*info != 0) return;
    const size_t block = batch.block_first+blockIdx.x;
    size_t lo = batch.first, hi = batch.first+batch.count;
    while (lo < hi) {
        const size_t middle = lo+(hi-lo)/2;
        if (tasks[middle].block_end <= block) lo = middle+1;
        else hi = middle;
    }
    const DeviceTask &item = tasks[lo];
    const auto &t = item.task;
    const size_t begin = lo == 0 ? 0 : tasks[lo-1].block_end;
    const size_t index = (block-begin)*blockDim.x+threadIdx.x;
    spd_schur_scatter_entry(arena+t.left_offset, t.ldl,
        arena+t.right_offset, t.ldr, t.m, t.n, t.k,
        arena+t.destination_offset, t.destination_rows, t.destination_columns,
        maps+t.row_map_offset, maps+t.column_map_offset, index);
}

void clear_plan(Context *ctx)
{
    if (ctx->tasks) cudaFree(ctx->tasks);
    ctx->tasks = nullptr;
    ctx->task_count = 0;
    std::vector<Batch>().swap(ctx->batches);
}
} // namespace

extern "C" int fsparse_spd_cuda_create_deferred(void **context, int device, size_t capacity)
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
    if (result) {
        fsparse_spd_cuda_destroy(ctx);
        return result;
    }
    *context = ctx;
    return 0;
}

extern "C" int fsparse_spd_cuda_prepare_blas(void *context)
{
    auto *ctx = static_cast<Context*>(context);
    int result = activate(ctx);
    if (result) return result;
    // This is setup, never a hidden allocation in a numerical factor batch.
    if (ctx->factor_batch) return -1;
    if (ctx->blas) return 0;
    result = cudaMalloc(&ctx->workspace, workspace_bytes);
    if (result) return result;
    cublasStatus_t code = cublasCreate(&ctx->blas);
    if (code == CUBLAS_STATUS_SUCCESS) code = cublasSetStream(ctx->blas, ctx->stream);
    if (code == CUBLAS_STATUS_SUCCESS)
        code = cublasSetWorkspace(ctx->blas, ctx->workspace, workspace_bytes);
    if (code == CUBLAS_STATUS_SUCCESS) return 0;
    if (ctx->blas) cublasDestroy(ctx->blas);
    ctx->blas = nullptr;
    cudaFree(ctx->workspace);
    ctx->workspace = nullptr;
    return 1000+code;
}

extern "C" int fsparse_spd_cuda_create(void **context, int device, size_t capacity)
{
    int result = fsparse_spd_cuda_create_deferred(context, device, capacity);
    if (!result) result = fsparse_spd_cuda_prepare_blas(*context);
    if (result && context && *context) {
        fsparse_spd_cuda_destroy(*context);
        *context = nullptr;
    }
    return result;
}

extern "C" int fsparse_spd_cuda_destroy(void *context)
{
    auto *ctx = static_cast<Context*>(context);
    if (!ctx) return 0;
    int result = activate(ctx);
    if (result) return result;
    if (ctx->stream) result = cudaStreamSynchronize(ctx->stream);
    if (ctx->blas) cublasDestroy(ctx->blas);
    clear_plan(ctx);
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

extern "C" int fsparse_spd_cuda_front_factor(void *context, size_t offset,
                                             int width, int rows, int lda)
{
    auto *ctx = static_cast<Context*>(context);
    if (width < 0 || width > 8 || rows < width ||
        !matrix_fits(ctx, offset, lda, rows, width)) return -1;
    int result = activate(ctx);
    if (result || width == 0) return result;
    factor_front<<<1,256,0,ctx->stream>>>(ctx->arena+offset, width, rows, lda,
                                        ctx->info, ctx->factor_batch);
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
    clear_plan(ctx);
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

extern "C" int fsparse_spd_cuda_schur_scatter(void *context,
                                             size_t left_offset, int ldl,
                                             size_t right_offset, int ldr,
                                             int m, int n, int k,
                                             size_t destination_offset,
                                             int destination_rows, int destination_columns,
                                             size_t row_map_offset, size_t column_map_offset)
{
    auto *ctx = static_cast<Context*>(context);
    if (!matrix_fits(ctx, left_offset, ldl, m, k) ||
        !matrix_fits(ctx, right_offset, ldr, n, k) ||
        !matrix_fits(ctx, destination_offset, destination_rows,
                     destination_rows, destination_columns)) return -1;
    if (row_map_offset > ctx->map_count || column_map_offset > ctx->map_count) return -1;
    if (static_cast<size_t>(m) > ctx->map_count-row_map_offset ||
        static_cast<size_t>(n) > ctx->map_count-column_map_offset) return -1;
    int result = activate(ctx);
    if (result || m == 0 || n == 0 || k == 0) return result;
    const size_t limit = std::numeric_limits<size_t>::max();
    if (static_cast<size_t>(m) > limit/static_cast<size_t>(n)) return -1;
    const size_t count = static_cast<size_t>(m)*n;
    const size_t blocks = count/256+(count%256 != 0);
    if (blocks > std::numeric_limits<unsigned>::max()) return -1;
    schur_scatter<<<static_cast<unsigned>(blocks),256,0,ctx->stream>>>(
        ctx->arena+left_offset, ldl, ctx->arena+right_offset, ldr, m, n, k,
        ctx->arena+destination_offset, destination_rows, destination_columns,
        ctx->maps+row_map_offset, ctx->maps+column_map_offset, ctx->info);
    return cudaGetLastError();
}

extern "C" int fsparse_spd_cuda_set_schur_batches(void *context,
    const fsparse_spd_cuda_schur_task *tasks, size_t task_count,
    const size_t *batch_ptr, size_t batch_count)
{
    auto *ctx = static_cast<Context*>(context);
    const size_t limit = std::numeric_limits<size_t>::max();
    if (!ctx || !batch_ptr || (task_count && !tasks) ||
        task_count > limit/sizeof(DeviceTask) || batch_count > limit/sizeof(Batch) ||
        batch_ptr[0] != 0 || batch_ptr[batch_count] != task_count) return -1;
    try {
        std::vector<DeviceTask> device_tasks(task_count);
        std::vector<Batch> batches(batch_count);
        size_t blocks = 0;
        for (size_t b = 0; b < batch_count; ++b) {
            const size_t first = batch_ptr[b], end = batch_ptr[b+1];
            if (first > end || end > task_count) return -1;
            Batch &batch = batches[b];
            batch = {first, end-first, blocks, 0};
            size_t destination_end = 0, first_destination = ctx->capacity;
            size_t source_end = 0;
            for (size_t i = first; i < end; ++i) {
                const auto &t = tasks[i];
                if (t.k < 0 || t.k > 8 ||
                    !matrix_fits(ctx, t.left_offset, t.ldl, t.m, t.k) ||
                    !matrix_fits(ctx, t.right_offset, t.ldr, t.n, t.k) ||
                    !matrix_fits(ctx, t.destination_offset, t.destination_rows,
                                 t.destination_rows, t.destination_columns) ||
                    t.row_map_offset > ctx->map_count ||
                    t.column_map_offset > ctx->map_count ||
                    size_t(t.m) > ctx->map_count-t.row_map_offset ||
                    size_t(t.n) > ctx->map_count-t.column_map_offset) return -1;
                if (t.m && t.n && t.k) {
                    // Matrix bounds above also prove these products cannot overflow.
                    if (size_t(t.m) > limit/size_t(t.n)) return -1;
                    const size_t entries = size_t(t.m)*t.n;
                    const size_t added = entries/256+(entries%256 != 0);
                    if (added > limit-blocks) return -1;
                    blocks += added;
                    if (blocks-batch.block_first > INT_MAX) return -1;
                    if (t.destination_offset < destination_end) return -1;
                    if (first_destination == ctx->capacity)
                        first_destination = t.destination_offset;
                    destination_end = t.destination_offset+
                        size_t(t.destination_rows)*t.destination_columns;
                    const size_t left_end = t.left_offset+size_t(t.k-1)*t.ldl+t.m;
                    const size_t right_end = t.right_offset+size_t(t.k-1)*t.ldr+t.n;
                    if (left_end > source_end) source_end = left_end;
                    if (right_end > source_end) source_end = right_end;
                }
                device_tasks[i] = {t, blocks};
            }
            if (source_end > first_destination) return -1;
            batch.blocks = blocks-batch.block_first;
        }
        int result = activate(ctx);
        if (result) return result;
        DeviceTask *uploaded = nullptr;
        if (task_count) {
            result = cudaMalloc(&uploaded, task_count*sizeof(DeviceTask));
            if (!result) result = cudaMemcpyAsync(uploaded, device_tasks.data(),
                task_count*sizeof(DeviceTask), cudaMemcpyHostToDevice, ctx->stream);
        }
        if (!result) result = cudaStreamSynchronize(ctx->stream);
        if (result) { if (uploaded) cudaFree(uploaded); return result; }
        clear_plan(ctx);
        ctx->tasks = uploaded;
        ctx->task_count = task_count;
        ctx->batches = std::move(batches);
        return 0;
    } catch (const std::bad_alloc &) {
        return static_cast<int>(cudaErrorMemoryAllocation);
    }
}

extern "C" int fsparse_spd_cuda_schur_batch(void *context, size_t batch)
{
    auto *ctx = static_cast<Context*>(context);
    if (!ctx || batch >= ctx->batches.size()) return -1;
    const Batch &item = ctx->batches[batch];
    if (!item.blocks) return 0;
    int result = activate(ctx);
    if (result) return result;
    schur_batch<<<static_cast<unsigned>(item.blocks),256,0,ctx->stream>>>(
        ctx->arena, ctx->maps, ctx->tasks, item, ctx->info);
    return cudaGetLastError();
}

extern "C" size_t fsparse_spd_cuda_host_plan_bytes(void *context)
{
    auto *ctx = static_cast<Context*>(context);
    return ctx ? sizeof(Context)+ctx->batches.capacity()*sizeof(Batch) : 0;
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
    if (!ctx->blas) return -1; // Call prepare_blas during setup first.
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
    if (!ctx->blas) return -1; // Never initialize vendor state in this hot call.
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
    return ctx ? ctx->capacity*sizeof(double)+sizeof(int)
                 +(ctx->workspace ? workspace_bytes : 0)
                 +ctx->map_count*sizeof(int)+ctx->task_count*sizeof(DeviceTask) : 0;
}
