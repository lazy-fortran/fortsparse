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
struct Level {
    size_t panel_first = 0, panel_count = 0;
    Batch updates = {};
};
struct SolvePanel {
    fsparse_spd_cuda_level_panel panel;
    size_t row_first;
    int first;
};
struct SolveLevel { size_t first, count; };
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
    bool factor_ready = false;
    SolvePanel *solve_panels = nullptr;
    int *solve_rows = nullptr;
    int solve_n = 0;
    size_t solve_panel_count = 0, solve_row_count = 0, solve_factor_end = 0;
    std::vector<SolveLevel> solve_levels;
    DeviceTask *tasks = nullptr;
    size_t task_count = 0;
    std::vector<Batch> batches;
    fsparse_spd_cuda_level_panel *level_panels = nullptr;
    size_t level_panel_count = 0;
    DeviceTask *level_tasks = nullptr;
    size_t level_task_count = 0;
    std::vector<Level> levels;
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

// A global failure can change while other CTAs execute. Each CTA snapshots it
// once into shared storage, then all threads take identical barrier branches.
__global__ void factor_level(double *arena,
    const fsparse_spd_cuda_level_panel *panels, size_t first, int *info)
{
    __shared__ int local_info;
    if (threadIdx.x == 0) local_info = atomicAdd(info, 0);
    __syncthreads();
    const auto &p = panels[first+blockIdx.x];
    spd_factor_panel_schedule(arena+p.offset, p.width, p.rows, p.lda,
                              &local_info, true);
    if (threadIdx.x == 0 && local_info != 0) atomicCAS(info, 0, local_info);
}

__global__ void update_level(double *arena, const int *maps,
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
    const auto &item = tasks[lo];
    const auto &t = item.task;
    const size_t begin = lo == 0 ? 0 : tasks[lo-1].block_end;
    spd_schur_scatter_entry(arena+t.left_offset, t.ldl,
        arena+t.right_offset, t.ldr, t.m, t.n, t.k,
        arena+t.destination_offset, t.destination_rows, t.destination_columns,
        maps+t.row_map_offset, maps+t.column_map_offset,
        (block-begin)*blockDim.x+threadIdx.x, true);
}

template<bool Forward>
__global__ void solve_level(double *arena, const SolvePanel *panels, size_t first,
    const int *rows, double *rhs, int ld, int nrhs, const int *info)
{
    const auto &s = panels[first+blockIdx.x];
    const auto &p = s.panel;
    const int *tail = rows ? rows+s.row_first : nullptr;
    if (Forward)
        spd_forward_solve_schedule(arena+p.offset,p.width,p.rows,p.lda,
                                  tail,s.first,rhs,ld,nrhs,info);
    else
        spd_backward_solve_schedule(arena+p.offset,p.width,p.rows,p.lda,
                                   tail,s.first,rhs,ld,nrhs,info);
}

void clear_solve(Context *ctx)
{
    if (ctx->solve_panels) cudaFree(ctx->solve_panels);
    if (ctx->solve_rows) cudaFree(ctx->solve_rows);
    ctx->solve_panels = nullptr;
    ctx->solve_rows = nullptr;
    ctx->solve_n = 0;
    ctx->solve_panel_count = ctx->solve_row_count = ctx->solve_factor_end = 0;
    ctx->factor_ready = false;
    std::vector<SolveLevel>().swap(ctx->solve_levels);
}

void clear_levels(Context *ctx)
{
    if (ctx->level_panels) cudaFree(ctx->level_panels);
    if (ctx->level_tasks) cudaFree(ctx->level_tasks);
    ctx->level_panels = nullptr;
    ctx->level_tasks = nullptr;
    ctx->level_panel_count = ctx->level_task_count = 0;
    std::vector<Level>().swap(ctx->levels);
}

void clear_schur(Context *ctx)
{
    if (ctx->tasks) cudaFree(ctx->tasks);
    ctx->tasks = nullptr;
    ctx->task_count = 0;
    std::vector<Batch>().swap(ctx->batches);
}
void clear_plan(Context *ctx)
{
    clear_schur(ctx);
    clear_levels(ctx);
    clear_solve(ctx);
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
    if (count && offset < ctx->solve_factor_end) ctx->factor_ready = false;
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
    ctx->factor_ready = false;
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
    ctx->factor_ready = false;
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
    ctx->factor_ready = result == 0 && *info == 0;
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
    ctx->factor_ready = false;
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
        clear_schur(ctx);
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

extern "C" int fsparse_spd_cuda_set_levels(void *context,
    const fsparse_spd_cuda_level_panel *panels, size_t panel_count,
    const fsparse_spd_cuda_level_task *tasks, size_t task_count, size_t level_count)
{
    auto *ctx = static_cast<Context*>(context);
    const size_t limit = std::numeric_limits<size_t>::max();
    if (!ctx || ctx->factor_batch || (panel_count && !panels) ||
        (task_count && !tasks) || level_count > INT_MAX || panel_count > INT_MAX ||
        (panel_count && !level_count) || (!panel_count && (level_count || task_count)) ||
        panel_count > limit/sizeof(*panels) || task_count > limit/sizeof(DeviceTask) ||
        level_count > limit/sizeof(Level) || level_count+1 > limit/sizeof(size_t))
        return -1;
    try {
        std::vector<Level> levels(level_count);
        std::vector<size_t> panel_ptr(level_count+1,0), task_ptr(level_count+1,0);
        size_t previous_end = 0, thin_count = 0;
        for (size_t i = 0; i < panel_count; ++i) {
            const auto &p = panels[i];
            if (p.width < 0 || p.rows < p.width || p.level < 0 ||
                size_t(p.level) >= level_count ||
                !matrix_fits(ctx,p.offset,p.lda,p.rows,p.width) ||
                size_t(p.width) > limit/size_t(p.lda) ||
                !fits(ctx,p.offset,size_t(p.lda)*p.width) || p.offset < previous_end)
                return -1;
            previous_end = p.offset+size_t(p.lda)*p.width;
            if (p.width <= 8 && p.width) { ++panel_ptr[p.level+1]; ++thin_count; }
        }
        for (size_t i = 0; i < task_count; ++i) {
            const auto &item = tasks[i];
            const auto &t = item.task;
            if (item.source >= panel_count || item.destination >= panel_count ||
                item.destination <= item.source) return -1;
            const auto &s = panels[item.source], &d = panels[item.destination];
            if (s.width <= 0 || s.width > 8 || s.level >= d.level ||
                t.k != s.width || t.ldl != s.lda || t.ldr != s.lda ||
                t.m < 0 || t.n < 0 ||
                !matrix_fits(ctx,t.left_offset,t.ldl,t.m,t.k) ||
                !matrix_fits(ctx,t.right_offset,t.ldr,t.n,t.k) ||
                t.left_offset < s.offset || t.right_offset < s.offset ||
                t.left_offset-s.offset > size_t(s.rows) ||
                t.right_offset-s.offset > size_t(s.rows) ||
                size_t(t.m) > size_t(s.rows)-(t.left_offset-s.offset) ||
                size_t(t.n) > size_t(s.rows)-(t.right_offset-s.offset) ||
                t.destination_offset != d.offset || t.destination_rows != d.lda ||
                t.destination_columns != d.width ||
                t.row_map_offset > ctx->map_count || t.column_map_offset > ctx->map_count ||
                size_t(t.m) > ctx->map_count-t.row_map_offset ||
                size_t(t.n) > ctx->map_count-t.column_map_offset ||
                (t.n && size_t(t.m) > limit/size_t(t.n))) return -1;
            ++task_ptr[s.level+1];
        }
        for (size_t l = 0; l < level_count; ++l) {
            panel_ptr[l+1] += panel_ptr[l]; task_ptr[l+1] += task_ptr[l];
            levels[l].panel_first = panel_ptr[l];
            levels[l].panel_count = panel_ptr[l+1]-panel_ptr[l];
        }
        std::vector<fsparse_spd_cuda_level_panel> device_panels(thin_count);
        std::vector<DeviceTask> device_tasks(task_count);
        auto panel_next = panel_ptr, task_next = task_ptr;
        for (size_t i = 0; i < panel_count; ++i)
            if (panels[i].width <= 8 && panels[i].width)
                device_panels[panel_next[panels[i].level]++] = panels[i];
        for (size_t i = 0; i < task_count; ++i)
            device_tasks[task_next[panels[tasks[i].source].level]++] = {tasks[i].task,0};
        size_t blocks = 0;
        for (size_t l = 0; l < level_count; ++l) {
            auto &batch = levels[l].updates;
            batch = {task_ptr[l],task_ptr[l+1]-task_ptr[l],blocks,0};
            for (size_t i = task_ptr[l]; i < task_ptr[l+1]; ++i) {
                const auto &t = device_tasks[i].task;
                const size_t entries = size_t(t.m)*t.n;
                const size_t added = entries/256+(entries%256 != 0);
                if (added > limit-blocks) return -1;
                blocks += added;
                if (blocks-batch.block_first > INT_MAX) return -1;
                device_tasks[i].block_end = blocks;
            }
            batch.blocks = blocks-batch.block_first;
        }
        int result = activate(ctx);
        fsparse_spd_cuda_level_panel *uploaded_panels = nullptr;
        DeviceTask *uploaded_tasks = nullptr;
        if (!result && thin_count)
            result = cudaMalloc(&uploaded_panels,thin_count*sizeof(*uploaded_panels));
        if (!result && thin_count) result = cudaMemcpyAsync(uploaded_panels,
            device_panels.data(),thin_count*sizeof(*uploaded_panels),
            cudaMemcpyHostToDevice,ctx->stream);
        if (!result && task_count)
            result = cudaMalloc(&uploaded_tasks,task_count*sizeof(*uploaded_tasks));
        if (!result && task_count) result = cudaMemcpyAsync(uploaded_tasks,
            device_tasks.data(),task_count*sizeof(*uploaded_tasks),
            cudaMemcpyHostToDevice,ctx->stream);
        if (!result) result = cudaStreamSynchronize(ctx->stream);
        if (result) {
            if (uploaded_panels) cudaFree(uploaded_panels);
            if (uploaded_tasks) cudaFree(uploaded_tasks);
            return result;
        }
        clear_levels(ctx);
        ctx->level_panels = uploaded_panels; ctx->level_panel_count = thin_count;
        ctx->level_tasks = uploaded_tasks; ctx->level_task_count = task_count;
        ctx->levels = std::move(levels);
        return 0;
    } catch (const std::bad_alloc &) {
        return static_cast<int>(cudaErrorMemoryAllocation);
    }
}

extern "C" int fsparse_spd_cuda_factor_level(void *context, size_t level)
{
    auto *ctx = static_cast<Context*>(context);
    if (!ctx || !ctx->factor_batch || level >= ctx->levels.size()) return -1;
    const auto &item = ctx->levels[level];
    if (!item.panel_count) return 0;
    int result = activate(ctx);
    if (result) return result;
    factor_level<<<static_cast<unsigned>(item.panel_count),256,0,ctx->stream>>>(
        ctx->arena,ctx->level_panels,item.panel_first,ctx->info);
    return cudaGetLastError();
}

extern "C" int fsparse_spd_cuda_update_level(void *context, size_t level)
{
    auto *ctx = static_cast<Context*>(context);
    if (!ctx || !ctx->factor_batch || level >= ctx->levels.size()) return -1;
    const auto &item = ctx->levels[level].updates;
    if (!item.blocks) return 0;
    int result = activate(ctx);
    if (result) return result;
    update_level<<<static_cast<unsigned>(item.blocks),256,0,ctx->stream>>>(
        ctx->arena,ctx->maps,ctx->level_tasks,item,ctx->info);
    return cudaGetLastError();
}

extern "C" size_t fsparse_spd_cuda_host_plan_bytes(void *context)
{
    auto *ctx = static_cast<Context*>(context);
    return ctx ? sizeof(Context)+ctx->batches.capacity()*sizeof(Batch)
                 +ctx->levels.capacity()*sizeof(Level)
                 +ctx->solve_levels.capacity()*sizeof(SolveLevel) : 0;
}

extern "C" int fsparse_spd_cuda_set_solve_plan(void *context, int n,
    const fsparse_spd_cuda_level_panel *panels, size_t panel_count,
    const int *tail_rows, size_t row_count, size_t level_count)
{
    auto *ctx = static_cast<Context*>(context);
    const size_t limit = std::numeric_limits<size_t>::max();
    if (!ctx || n <= 0 || !panels || !panel_count || !level_count ||
        panel_count > size_t(n) || level_count > panel_count ||
        (row_count && !tail_rows) || row_count > limit/sizeof(int)) return -1;
    try {
        std::vector<SolvePanel> ordered(panel_count);
        std::vector<int> owner(n);
        std::vector<size_t> ptr(level_count+1,0);
        std::vector<SolveLevel> levels(level_count);
        size_t row = 0, factor_end = 0;
        int first = 0;
        for (size_t i = 0; i < panel_count; ++i) {
            const auto &p = panels[i];
            if (p.width <= 0 || p.width > 128 || p.width > n-first ||
                p.rows < p.width || p.level < 0 || size_t(p.level) >= level_count ||
                !matrix_fits(ctx,p.offset,p.lda,p.rows,p.width) ||
                !fits(ctx,p.offset,size_t(p.width)*p.lda) || p.offset < factor_end ||
                size_t(p.rows-p.width) > row_count-row) return -1;
            factor_end = p.offset+size_t(p.width)*p.lda;
            ordered[i] = {p,row,first};
            for (int j = 0; j < p.width; ++j) owner[first+j] = i;
            first += p.width;
            row += p.rows-p.width;
            ++ptr[p.level+1];
        }
        if (first != n || row != row_count) return -1;
        for (const auto &item : ordered) {
            const auto &p = item.panel;
            int previous = item.first+p.width;
            for (int j = 0; j < p.rows-p.width; ++j) {
                const int target = tail_rows[item.row_first+j];
                if (target <= previous || target > n ||
                    panels[owner[target-1]].level <= p.level) return -1;
                previous = target;
            }
        }
        for (size_t l = 0; l < level_count; ++l) {
            ptr[l+1] += ptr[l];
            levels[l] = {ptr[l],ptr[l+1]-ptr[l]};
        }
        auto next = ptr;
        std::vector<SolvePanel> upload(panel_count);
        for (const auto &item : ordered) upload[next[item.panel.level]++] = item;
        int result = activate(ctx);
        if (result) return result;
        SolvePanel *device_panels = nullptr;
        int *device_rows = nullptr;
        result = cudaMalloc(&device_panels,panel_count*sizeof(SolvePanel));
        if (!result && row_count) result = cudaMalloc(&device_rows,row_count*sizeof(int));
        if (!result) result = cudaMemcpyAsync(device_panels,upload.data(),
            panel_count*sizeof(SolvePanel),cudaMemcpyHostToDevice,ctx->stream);
        if (!result && row_count) result = cudaMemcpyAsync(device_rows,tail_rows,
            row_count*sizeof(int),cudaMemcpyHostToDevice,ctx->stream);
        if (!result) result = cudaStreamSynchronize(ctx->stream);
        if (result) {
            if (device_panels) cudaFree(device_panels);
            if (device_rows) cudaFree(device_rows);
            return result;
        }
        clear_solve(ctx);
        ctx->solve_panels = device_panels; ctx->solve_rows = device_rows;
        ctx->solve_panel_count = panel_count; ctx->solve_row_count = row_count;
        ctx->solve_factor_end = factor_end; ctx->solve_n = n;
        ctx->solve_levels = std::move(levels);
        return 0;
    } catch (const std::bad_alloc &) { return int(cudaErrorMemoryAllocation); }
}

extern "C" int fsparse_spd_cuda_solve(void *context, size_t rhs_offset, int ld, int nrhs)
{
    auto *ctx = static_cast<Context*>(context);
    if (!ctx || !ctx->factor_ready || !ctx->solve_panels || nrhs < 0 || nrhs > 8 ||
        rhs_offset < ctx->solve_factor_end ||
        !matrix_fits(ctx,rhs_offset,ld,ctx->solve_n,nrhs)) return -1;
    if (!nrhs) return 0;
    int result = activate(ctx);
    if (result) return result;
    for (const auto &level : ctx->solve_levels) {
        if (!level.count) continue;
        solve_level<true><<<unsigned(level.count),256,0,ctx->stream>>>(ctx->arena,
            ctx->solve_panels,level.first,ctx->solve_rows,ctx->arena+rhs_offset,ld,nrhs,ctx->info);
        result = cudaGetLastError(); if (result) return result;
    }
    for (size_t l = ctx->solve_levels.size(); l-- > 0;) {
        const auto &level = ctx->solve_levels[l];
        if (!level.count) continue;
        solve_level<false><<<unsigned(level.count),256,0,ctx->stream>>>(ctx->arena,
            ctx->solve_panels,level.first,ctx->solve_rows,ctx->arena+rhs_offset,ld,nrhs,ctx->info);
        result = cudaGetLastError(); if (result) return result;
    }
    return 0;
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
                 +ctx->map_count*sizeof(int)+ctx->task_count*sizeof(DeviceTask)
                 +ctx->level_task_count*sizeof(DeviceTask)
                 +ctx->level_panel_count*sizeof(fsparse_spd_cuda_level_panel)
                 +ctx->solve_panel_count*sizeof(SolvePanel)+ctx->solve_row_count*sizeof(int) : 0;
}
