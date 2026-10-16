// Per-output contraction and mapped lower scatter. Hardware launch geometry
// belongs to the native adapter; arithmetic uses the retained generated leaf.
#include <cstddef>
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
