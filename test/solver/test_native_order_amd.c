#include <stdio.h>
#include <string.h>

int fsparse_native_order_amd(int, int, const int *, const int *, int *, int *);
int fsparse_native_order_amd_general(int, int, const int *, const int *, int *, int *);

static int check_inverse(const int *perm, const int *inverse, int n) {
    for (int i = 0; i < n; ++i) {
        if (perm[i] < 1 || perm[i] > n) return 0;
        if (inverse[perm[i] - 1] != i + 1) return 0;
    }
    return 1;
}

int main(void) {
    /* A cycle with a chord, one attached vertex, and one isolated vertex. */
    const int ptr[7] = {1, 5, 8, 12, 16, 18, 19};
    const int row[18] = {1, 2, 3, 4, 1, 2, 3, 1, 2, 3, 4,
                         1, 3, 4, 5, 4, 5, 6};
    int perm[6], inverse[6], again[6], again_inverse[6];
    int status = fsparse_native_order_amd(6, 18, ptr, row, perm, inverse);
#ifdef FORTSPARSE_TEST_AMD
    if (status != 0 || !check_inverse(perm, inverse, 6)) return 1;
    if (fsparse_native_order_amd(6, 18, ptr, row, again, again_inverse) != 0)
        return 2;
    if (memcmp(perm, again, sizeof perm) ||
        memcmp(inverse, again_inverse, sizeof inverse)) return 3;
    const int bad_ptr[3] = {1, 4, 3}, bad_row[2] = {1, 3};
    const int two_ptr[3] = {1, 2, 3}, asymmetric_row[2] = {1, 1};
    if (fsparse_native_order_amd(2, 2, bad_ptr, row, perm, inverse) != 2)
        return 4;
    if (fsparse_native_order_amd(2, 2, two_ptr, bad_row, perm, inverse) != 2)
        return 5;
    if (fsparse_native_order_amd(2, 2, two_ptr, asymmetric_row,
                                perm, inverse) != 2) return 6;
    const int duplicate_ptr[3] = {1, 3, 3}, duplicate_row[2] = {1, 1};
    if (fsparse_native_order_amd(2, 2, duplicate_ptr, duplicate_row,
                                perm, inverse) != 2) return 7;
    const int diagonal_ptr[4] = {1, 2, 3, 4}, diagonal_row[3] = {1, 2, 3};
    if (fsparse_native_order_amd(3, 3, diagonal_ptr, diagonal_row,
                                perm, inverse) != 0 ||
        !check_inverse(perm, inverse, 3)) return 8;
    const int empty_ptr[4] = {1, 1, 1, 1};
    if (fsparse_native_order_amd(3, 0, empty_ptr, NULL, perm, inverse) != 0 ||
        !check_inverse(perm, inverse, 3)) return 9;
    /* Genuinely one-way edges: AMD forms its own A+A^T ordering graph. */
    const int general_ptr[6] = {1, 3, 5, 7, 9, 10};
    const int general_row[9] = {1, 3, 2, 5, 3, 4, 1, 4, 5};
    if (fsparse_native_order_amd_general(5, 9, general_ptr, general_row,
                                        perm, inverse) != 0 ||
        !check_inverse(perm, inverse, 5)) return 11;
    if (fsparse_native_order_amd(5, 9, general_ptr, general_row, perm, inverse) != 2)
        return 12;
    if (fsparse_native_order_amd_general(2, 2, bad_ptr, row, perm, inverse) != 2 ||
        fsparse_native_order_amd_general(2, 2, two_ptr, bad_row, perm, inverse) != 2 ||
        fsparse_native_order_amd_general(2, 2, duplicate_ptr, duplicate_row,
                                        perm, inverse) != 2) return 13;
#else
    (void)check_inverse; (void)again; (void)again_inverse;
    if (status != 3) return 10;
    if (fsparse_native_order_amd_general(6, 18, ptr, row, perm, inverse) != 3)
        return 14;
#endif
    puts("PASS AMD permutation, inverse, deterministic ordering and status");
    return 0;
}
