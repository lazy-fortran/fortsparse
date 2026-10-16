/* MIT wrapper around the optional, separately linked Apache-2.0 METIS.
 * No ordering implementation or third-party headers are vendored here. */
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>

#ifdef FORTSPARSE_HAVE_METIS
#include <metis.h>
#include "fsparse_order_graph.h"
#endif

int fsparse_native_order_available(void) {
#ifdef FORTSPARSE_HAVE_METIS
    return 1;
#else
    return 0;
#endif
}

#ifdef FORTSPARSE_HAVE_METIS
static int compare_index(const void *left, const void *right) {
    const idx_t a = *(const idx_t *)left, b = *(const idx_t *)right;
    return (a > b) - (a < b);
}
#endif

/* Input: canonical CSC, Fortran one-based indices. General input orders the
 * undirected union of A and its transpose, with no repeated or self edges.
 * Output: perm[new] = old, inverse[old] = new, both one-based.
 * Status: 0 OK, 2 invalid matrix, 3 unavailable, 5 METIS failure, 6 OOM. */
static int order_metis(int n, int nnz, const int *ptr, const int *row,
                       int *perm, int *inverse, int general) {
#ifndef FORTSPARSE_HAVE_METIS
    (void)n; (void)nnz; (void)ptr; (void)row; (void)perm; (void)inverse;
    (void)general;
    return 3;
#else
    idx_t *xadj = NULL, *adj = NULL, *p = NULL, *ip = NULL, *next = NULL;
    idx_t options[METIS_NOPTIONS], vertices, count = 0;
    int result = 0;
    if (!perm || !inverse) return 2;
    if (general ? !fsparse_csc_graph_valid(n, nnz, ptr, row) :
                  !fsparse_graph_valid(n, nnz, ptr, row)) return 2;
    const uint64_t capacity = (uint64_t)nnz * (general ? 2u : 1u);
    if (capacity > IDX_MAX) return 6;
    if ((size_t)n + 1 > SIZE_MAX / sizeof(idx_t) ||
        capacity > SIZE_MAX / sizeof(idx_t)) return 6;
    xadj = malloc(((size_t)n + 1) * sizeof(*xadj));
    adj = malloc((capacity ? (size_t)capacity : 1) * sizeof(*adj));
    p = malloc((size_t)n * sizeof(*p));
    ip = malloc((size_t)n * sizeof(*ip));
    if (general) next = calloc((size_t)n, sizeof(*next));
    if (!xadj || !adj || !p || !ip || (general && !next)) {
        result = 6;
        goto done;
    }
    if (general) {
        for (int col = 0; col < n; ++col) {
            for (int entry = ptr[col] - 1; entry < ptr[col + 1] - 1; ++entry) {
                const int r = row[entry] - 1;
                if (r == col) continue;
                ++next[col];
                ++next[r];
            }
        }
        xadj[0] = 0;
        for (int col = 0; col < n; ++col) {
            xadj[col + 1] = xadj[col] + next[col];
            next[col] = xadj[col];
        }
        for (int col = 0; col < n; ++col) {
            for (int entry = ptr[col] - 1; entry < ptr[col + 1] - 1; ++entry) {
                const int r = row[entry] - 1;
                if (r == col) continue;
                adj[next[col]++] = (idx_t)r;
                adj[next[r]++] = (idx_t)col;
            }
        }
        for (int col = 0; col < n; ++col) {
            const idx_t begin = xadj[col], end = xadj[col + 1];
            qsort(adj + begin, (size_t)(end - begin), sizeof(*adj), compare_index);
            xadj[col] = count;
            idx_t previous = -1;
            for (idx_t entry = begin; entry < end; ++entry) {
                const idx_t r = adj[entry];
                if (r != previous) adj[count++] = r;
                previous = r;
            }
        }
    } else {
        for (int col = 0; col < n; ++col) {
            xadj[col] = count;
            for (int entry = ptr[col] - 1; entry < ptr[col + 1] - 1; ++entry) {
                if (row[entry] != col + 1) adj[count++] = (idx_t)(row[entry] - 1);
            }
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
    free(xadj); free(adj); free(p); free(ip); free(next);
    return result;
#endif
}

/* Preserve the original symmetric-input contract used by SPD factors. */
int fsparse_native_order(int n, int nnz, const int *ptr, const int *row,
                        int *perm, int *inverse) {
    return order_metis(n, nnz, ptr, row, perm, inverse, 0);
}

int fsparse_native_order_general(int n, int nnz, const int *ptr, const int *row,
                                int *perm, int *inverse) {
    return order_metis(n, nnz, ptr, row, perm, inverse, 1);
}
