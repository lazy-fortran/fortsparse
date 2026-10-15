/* MIT: common validation for the native ordering adapters. */
#ifndef FSPARSE_ORDER_GRAPH_H
#define FSPARSE_ORDER_GRAPH_H
#include <limits.h>

static int fsparse_graph_has_row(const int *ptr, const int *row,
                                int column, int wanted) {
    int lo = ptr[column] - 1, hi = ptr[column + 1] - 1;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (row[mid] < wanted) lo = mid + 1;
        else hi = mid;
    }
    return lo < ptr[column + 1] - 1 && row[lo] == wanted;
}

/* Full symmetric canonical CSC, with one-based indices. */
static int fsparse_graph_valid(int n, int nnz, const int *ptr, const int *row) {
    if (n <= 0 || nnz < 0 || nnz == INT_MAX || !ptr) return 0;
    if (nnz && !row) return 0;
    if (ptr[0] != 1 || ptr[n] != nnz + 1) return 0;
    for (int col = 0; col < n; ++col) {
        if (ptr[col] < 1 || ptr[col] > ptr[col + 1] || ptr[col + 1] > nnz + 1)
            return 0;
        int previous = 0;
        for (int entry = ptr[col] - 1; entry < ptr[col + 1] - 1; ++entry) {
            if (row[entry] < 1 || row[entry] > n || row[entry] <= previous)
                return 0;
            previous = row[entry];
        }
    }
    for (int col = 0; col < n; ++col) {
        for (int entry = ptr[col] - 1; entry < ptr[col + 1] - 1; ++entry) {
            if (row[entry] != col + 1 &&
                !fsparse_graph_has_row(ptr, row, row[entry] - 1, col + 1))
                return 0;
        }
    }
    return 1;
}
#endif
