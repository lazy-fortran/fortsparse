// Per-output contraction and mapped lower scatter. Hardware launch geometry
// belongs to the native adapter; arithmetic uses the retained generated leaf.
#include <cstddef>
// Cooperative lower-panel factorization with optional trailing rows. The
// diagonal and every trailing row follow the same generated scalar leaves.
// All block threads must participate; rows may exceed the block width.
__device__ __forceinline__ void spd_factor_panel_schedule(
    double *a, int width, int rows, int lda, int *info, bool retain_info)
{
    if (!retain_info && threadIdx.x == 0) *info = 0;
    __syncthreads();
    if (*info != 0) return;
    for (int j = 0; j < width; ++j) {
        if (threadIdx.x == 0) {
            double pivot = a[j+static_cast<std::size_t>(j)*lda];
            for (int p = 0; p < j; ++p)
                spd_update_value(pivot, a[j+static_cast<std::size_t>(p)*lda],
                                 a[j+static_cast<std::size_t>(p)*lda], &pivot);
            if (!(pivot > 0.0) || !isfinite(pivot)) *info = j+1;
            else spd_root_value(pivot, a+j+static_cast<std::size_t>(j)*lda);
        }
        __syncthreads();
        if (*info != 0) return;
        for (int i = j+1+threadIdx.x; i < rows; i += blockDim.x) {
            double value = a[i+static_cast<std::size_t>(j)*lda];
            for (int p = 0; p < j; ++p)
                spd_update_value(value, a[i+static_cast<std::size_t>(p)*lda],
                                 a[j+static_cast<std::size_t>(p)*lda], &value);
            spd_divide_value(value, a[j+static_cast<std::size_t>(j)*lda],
                             a+i+static_cast<std::size_t>(j)*lda);
        }
        __syncthreads();
    }
}

__device__ __forceinline__ void spd_schur_scatter_entry(
    const double *left, int ldl, const double *right, int ldr,
    int m, int n, int k, double *destination, int ldd, int columns,
    const int *rows, const int *cols, std::size_t index, bool atomic = false)
{
    if (index >= static_cast<std::size_t>(m)*n) return;
    const int i = index % m, j = index / m;
    if (i < j) return;
    const int row = rows[i], column = cols[j];
    if (row <= 0 || row > ldd || column <= 0 || column > columns) return;
    double product = 0.0;
    for (int p = 0; p < k; ++p)
        spd_update_value(product, left[i+static_cast<std::size_t>(p)*ldl],
                         right[j+static_cast<std::size_t>(p)*ldr], &product);
    double *entry = destination+(row-1)+static_cast<std::size_t>(column-1)*ldd;
    // The same generated negative-dot increment feeds serial and parallel
    // hardware schedules; only the accumulation ordering changes.
    if (atomic) atomicAdd(entry, product);
    else spd_update_value(*entry, -1.0, product, entry);
}

// Bounded host-workflow solve tile: one cooperative block owns one panel.
// Tail rows belong to strictly later levels; only forward scatters need atomics.
// Arithmetic uses the same generated update/divide leaves as the CPU reference.
__device__ __forceinline__ void spd_forward_solve_schedule(
    const double *factor, int width, int rows, int lda, const int *tail_rows,
    int first, double *rhs, int ldrhs, int nrhs, const int *info)
{
    __shared__ double values[128*8];
    if (*info != 0) return;
    for (int index = threadIdx.x; index < width*nrhs; index += blockDim.x) {
        const int i = index%width, r = index/width;
        values[i+128*r] = rhs[first+i+static_cast<std::size_t>(r)*ldrhs];
    }
    __syncthreads();
    for (int j = 0; j < width; ++j) {
        for (int r = threadIdx.x; r < nrhs; r += blockDim.x)
            spd_divide_value(values[j+128*r], factor[j+static_cast<std::size_t>(j)*lda],
                             values+j+128*r);
        __syncthreads();
        const int remaining = width-j-1;
        for (int index = threadIdx.x; index < remaining*nrhs; index += blockDim.x) {
            const int i = j+1+index%remaining, r = index/remaining;
            spd_update_value(values[i+128*r], factor[i+static_cast<std::size_t>(j)*lda],
                             values[j+128*r], values+i+128*r);
        }
        __syncthreads();
    }
    for (int index = threadIdx.x; index < width*nrhs; index += blockDim.x) {
        const int i = index%width, r = index/width;
        rhs[first+i+static_cast<std::size_t>(r)*ldrhs] = values[i+128*r];
    }
    const int tail = rows-width;
    for (std::size_t index = threadIdx.x; index < static_cast<std::size_t>(tail)*nrhs;
         index += blockDim.x) {
        const int i = index%tail, r = index/tail;
        double contribution = 0.0;
        for (int j = 0; j < width; ++j)
            spd_update_value(contribution, factor[width+i+static_cast<std::size_t>(j)*lda],
                             values[j+128*r], &contribution);
        atomicAdd(rhs+tail_rows[i]-1+static_cast<std::size_t>(r)*ldrhs, contribution);
    }
}

__device__ __forceinline__ void spd_backward_solve_schedule(
    const double *factor, int width, int rows, int lda, const int *tail_rows,
    int first, double *rhs, int ldrhs, int nrhs, const int *info)
{
    __shared__ double values[128*8];
    if (*info != 0) return;
    for (int index = threadIdx.x; index < width*nrhs; index += blockDim.x) {
        const int j = index%width, r = index/width;
        double value = rhs[first+j+static_cast<std::size_t>(r)*ldrhs];
        for (int i = 0; i < rows-width; ++i)
            spd_update_value(value, factor[width+i+static_cast<std::size_t>(j)*lda],
                             rhs[tail_rows[i]-1+static_cast<std::size_t>(r)*ldrhs], &value);
        values[j+128*r] = value;
    }
    __syncthreads();
    for (int j = width-1; j >= 0; --j) {
        for (int r = threadIdx.x; r < nrhs; r += blockDim.x)
            spd_divide_value(values[j+128*r], factor[j+static_cast<std::size_t>(j)*lda],
                             values+j+128*r);
        __syncthreads();
        for (int index = threadIdx.x; index < j*nrhs; index += blockDim.x) {
            const int i = index%j, r = index/j;
            spd_update_value(values[i+128*r], factor[j+static_cast<std::size_t>(i)*lda],
                             values[j+128*r], values+i+128*r);
        }
        __syncthreads();
    }
    for (int index = threadIdx.x; index < width*nrhs; index += blockDim.x) {
        const int i = index%width, r = index/width;
        rhs[first+i+static_cast<std::size_t>(r)*ldrhs] = values[i+128*r];
    }
}
