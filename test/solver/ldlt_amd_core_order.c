/* The built-in AMD core on a full symmetric one-based CSC pattern, also in
 * builds whose AMD ordering is SuiteSparse's: perm[new] = old, one-based.
 * Returns 0, or nonzero when the core or an allocation fails. */
#include <stdint.h>
#include <stdlib.h>
/* relative to this file, so fpm builds find it without an include path */
#include "../../src/native/fsparse_native_amd_core.h"

int test_amd_core_order(int n, const int *ptr, const int *row, int *perm) {
    int32_t m = 0;
    for (int j = 0; j < n; ++j)
        for (int p = ptr[j] - 1; p < ptr[j + 1] - 1; ++p)
            if (row[p] - 1 != j) ++m;
    int32_t iwlen = m + m / 5 + n;
    int32_t *work = malloc(sizeof(int32_t) * (9 * (size_t)n + (size_t)iwlen));
    if (!work) return 1;
    int32_t *pe = work, *len = pe + n, *nv = len + n, *next = nv + n,
            *last = next + n, *head = last + n, *elen = head + n,
            *degree = elen + n, *w = degree + n, *iw = w + n;
    int32_t used = 0;
    for (int j = 0; j < n; ++j) {
        pe[j] = used;
        for (int p = ptr[j] - 1; p < ptr[j + 1] - 1; ++p)
            if (row[p] - 1 != j) iw[used++] = row[p] - 1;
        len[j] = used - pe[j];
    }
    int code = fsparse_amd_core(n, pe, iw, len, iwlen, m, nv, next, last,
                                head, elen, degree, w);
    for (int k = 0; k < n && code == 0; ++k) perm[k] = last[k] + 1;
    free(work);
    return code;
}
