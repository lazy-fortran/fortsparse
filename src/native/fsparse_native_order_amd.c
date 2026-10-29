/* MIT wrapper; AMD's BSD-3-Clause implementation remains an external library.
 * Without it, the built-in AMD core is used through the bounded adapter. */
#include <stdint.h>
#include <stdlib.h>
#include "fsparse_order_graph.h"
#ifdef FORTSPARSE_HAVE_AMD
#include <amd.h>
#else
#include "fsparse_native_amd_bounded.h"
#endif

/* Capability query without allocating a graph or invoking the vendor. */
int fsparse_native_order_amd_available(void) {
    return 1;
}

/* Same public index/status conventions as fsparse_native_order (METIS).
 * perm[new] = old; inverse[old] = new; all are one-based. */
static int order_amd(int n, int nnz, const int *ptr, const int *row,
                     int *perm, int *inverse, int symmetric) {
#ifndef FORTSPARSE_HAVE_AMD
    int64_t arena_bytes = 0, output_bytes = 0;
    if (!perm || !inverse) return 2;
    if (symmetric && !fsparse_graph_valid(n, nnz, ptr, row)) return 2;
    int code = fsparse_native_amd_bounded_size(n, nnz, ptr, row, 0,
                                               &arena_bytes, &output_bytes);
    if (code != 0) return code;
    if ((uint64_t)arena_bytes > SIZE_MAX) return 6;
    int32_t *arena = malloc(arena_bytes ? (size_t)arena_bytes : 1);
    if (!arena) return 6;
    code = fsparse_native_amd_bounded(n, nnz, ptr, row, 0, arena, arena_bytes,
                                      perm, inverse);
    free(arena);
    return code;
#else
    int32_t *ap = NULL, *ai = NULL, *p = NULL, *ip = NULL;
    int result = 0;
    if (!perm || !inverse) return 2;
    if (symmetric ? !fsparse_graph_valid(n, nnz, ptr, row) :
                    !fsparse_csc_graph_valid(n, nnz, ptr, row)) return 2;
    if ((size_t)n + 1 > SIZE_MAX / sizeof(int32_t) ||
        (size_t)nnz > SIZE_MAX / sizeof(int32_t)) return 6;
    ap = malloc(((size_t)n + 1) * sizeof(*ap));
    ai = malloc((nnz ? (size_t)nnz : 1) * sizeof(*ai));
    p = malloc((size_t)n * sizeof(*p));
    ip = calloc((size_t)n, sizeof(*ip));
    if (!ap || !ai || !p || !ip) {
        result = 6;
        goto done;
    }
    for (int col = 0; col <= n; ++col) ap[col] = (int32_t)(ptr[col] - 1);
    for (int entry = 0; entry < nnz; ++entry) ai[entry] = (int32_t)(row[entry] - 1);
    int status = amd_order((int32_t)n, ap, ai, p, NULL, NULL);
    if (status != AMD_OK && status != AMD_OK_BUT_JUMBLED) {
        result = status == AMD_OUT_OF_MEMORY ? 6 :
                 status == AMD_INVALID ? 2 : 5;
        goto done;
    }
    for (int i = 0; i < n; ++i) {
        if (p[i] < 0 || p[i] >= n || ip[p[i]] != 0) {
            result = 5;
            goto done;
        }
        ip[p[i]] = (int32_t)i + 1;
    }
    for (int i = 0; i < n; ++i) {
        perm[i] = (int)p[i] + 1;
        inverse[i] = (int)ip[i];
    }
done:
    free(ap); free(ai); free(p); free(ip);
    return result;
#endif
}

/* The original ABI retains structural-symmetry validation for SPD callers. */
int fsparse_native_order_amd(int n, int nnz, const int *ptr, const int *row,
                            int *perm, int *inverse) {
    return order_amd(n, nnz, ptr, row, perm, inverse, 1);
}

/* AMD accepts unsymmetric input and internally orders the A+A^T pattern.
 * Source: AMD/Include/amd.h (BSD-3-Clause), input contract of amd_order.
 * LU consumes this as a column permutation; numerical row pivots stay native. */
int fsparse_native_order_amd_general(int n, int nnz, const int *ptr, const int *row,
                                    int *perm, int *inverse) {
    return order_amd(n, nnz, ptr, row, perm, inverse, 0);
}
