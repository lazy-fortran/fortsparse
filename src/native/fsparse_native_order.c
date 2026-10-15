/* MIT wrapper around the optional, separately linked Apache-2.0 METIS.
 * No ordering implementation or third-party headers are vendored here. */
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>

#ifdef FORTSPARSE_HAVE_METIS
#include <metis.h>

static int has_row(const int *ptr, const int *row, int column, int wanted) {
    int lo = ptr[column] - 1, hi = ptr[column + 1] - 1;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (row[mid] < wanted) lo = mid + 1;
        else hi = mid;
    }
    return lo < ptr[column + 1] - 1 && row[lo] == wanted;
}
#endif

/* Input: full symmetric canonical CSC, Fortran one-based indices.
 * Output: perm[new] = old, inverse[old] = new, both one-based.
 * Status: 0 OK, 2 invalid matrix, 3 unavailable, 5 METIS failure, 6 OOM. */
int fsparse_native_order(int n, int nnz, const int *ptr, const int *row,
                        int *perm, int *inverse) {
#ifndef FORTSPARSE_HAVE_METIS
    (void)n; (void)nnz; (void)ptr; (void)row; (void)perm; (void)inverse;
    return 3;
#else
    idx_t *xadj = NULL, *adj = NULL, *p = NULL, *ip = NULL;
    idx_t options[METIS_NOPTIONS], vertices, count = 0;
    int result = 0;
    if (n <= 0 || nnz < 0 || nnz == INT_MAX || !ptr || !perm || !inverse)
        return 2;
    if (nnz && !row) return 2;
    if (ptr[0] != 1 || ptr[n] != nnz + 1) return 2;
    for (int col = 0; col < n; ++col) {
        if (ptr[col] < 1 || ptr[col] > ptr[col + 1] || ptr[col + 1] > nnz + 1)
            return 2;
        int previous = 0;
        for (int entry = ptr[col] - 1; entry < ptr[col + 1] - 1; ++entry) {
            if (row[entry] < 1 || row[entry] > n || row[entry] <= previous)
                return 2;
            previous = row[entry];
        }
    }
    for (int col = 0; col < n; ++col) {
        for (int entry = ptr[col] - 1; entry < ptr[col + 1] - 1; ++entry) {
            if (row[entry] != col + 1 &&
                !has_row(ptr, row, row[entry] - 1, col + 1)) return 2;
        }
    }
    if ((size_t)n + 1 > SIZE_MAX / sizeof(idx_t) ||
        (size_t)nnz > SIZE_MAX / sizeof(idx_t)) return 6;
    xadj = malloc(((size_t)n + 1) * sizeof(*xadj));
    adj = malloc((nnz ? (size_t)nnz : 1) * sizeof(*adj));
    p = malloc((size_t)n * sizeof(*p));
    ip = malloc((size_t)n * sizeof(*ip));
    if (!xadj || !adj || !p || !ip) {
        result = 6;
        goto done;
    }
    for (int col = 0; col < n; ++col) {
        xadj[col] = count;
        for (int entry = ptr[col] - 1; entry < ptr[col + 1] - 1; ++entry) {
            if (row[entry] != col + 1) adj[count++] = (idx_t)(row[entry] - 1);
        }
    }
    xadj[n] = count;
    vertices = (idx_t)n;
    int status = METIS_SetDefaultOptions(options);
    if (status != METIS_OK) {
        result = 5;
        goto done;
    }
    options[METIS_OPTION_SEED] = 0;
    options[METIS_OPTION_NUMBERING] = 0;
    status = METIS_NodeND(&vertices, xadj, adj, NULL, options, p, ip);
    if (status != METIS_OK) {
        result = status == METIS_ERROR_MEMORY ? 6 :
                 status == METIS_ERROR_INPUT ? 2 : 5;
        goto done;
    }
    for (int i = 0; i < n; ++i) {
        if (p[i] < 0 || p[i] >= n || ip[i] < 0 || ip[i] >= n) {
            result = 5;
            goto done;
        }
    }
    for (int i = 0; i < n; ++i) {
        if (ip[p[i]] != i) {
            result = 5;
            goto done;
        }
    }
    for (int i = 0; i < n; ++i) {
        perm[i] = (int)p[i] + 1;
        inverse[i] = (int)ip[i] + 1;
    }
done:
    free(xadj); free(adj); free(p); free(ip);
    return result;
#endif
}
