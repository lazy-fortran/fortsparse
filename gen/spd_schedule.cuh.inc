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
    const int *rows, const int *cols, std::size_t index)
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
    spd_update_value(*entry, -1.0, product, entry);
}
