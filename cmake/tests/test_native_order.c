#include <stdio.h>
#include <string.h>

int fsparse_native_order(int, int, const int *, const int *, int *, int *);
int fsparse_native_order_general(int, int, const int *, const int *, int *, int *);
int fsparse_native_order_available(void);

static int check_inverse(const int *perm, const int *inverse, int n) {
    for (int i = 0; i < n; ++i) {
        if (perm[i] < 1 || perm[i] > n) return 0;
        if (inverse[perm[i] - 1] != i + 1) return 0;
    }
    return 1;
}

int main(void) {
    /* A path plus diagonal and an isolated vertex: connected-component
     * handling and removal of self edges must preserve a bijection. */
    const int ptr[6] = {1, 3, 6, 9, 11, 12};
    const int row[11] = {1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5};
    int perm[5], inverse[5], again[5], again_inverse[5];
    int status = fsparse_native_order(5, 11, ptr, row, perm, inverse);
#ifdef FORTSPARSE_TEST_METIS
    if (fsparse_native_order_available() != 1) return 10;
    if (status != 0 || !check_inverse(perm, inverse, 5)) return 1;
    if (fsparse_native_order(5, 11, ptr, row, again, again_inverse) != 0)
        return 2;
    if (memcmp(perm, again, sizeof perm) ||
        memcmp(inverse, again_inverse, sizeof inverse)) return 3;
    /* Partial mirrors, one-way edges, self edges and an isolated vertex.
     * The independent full undirected graph above must give the same order. */
    const int general_ptr[6] = {1, 3, 6, 8, 9, 10};
    const int general_row[9] = {1, 2, 1, 2, 3, 3, 4, 4, 5};
    if (fsparse_native_order_general(5, 9, general_ptr, general_row,
                                    again, again_inverse) != 0) return 11;
    if (memcmp(perm, again, sizeof perm) ||
        memcmp(inverse, again_inverse, sizeof inverse)) return 12;
    const int bad_ptr[3] = {1, 4, 3}, bad_row[2] = {1, 3};
    const int two_ptr[3] = {1, 2, 3}, asymmetric_row[2] = {1, 1};
    if (fsparse_native_order(2, 2, bad_ptr, row, perm, inverse) != 2)
        return 4;
    if (fsparse_native_order(2, 2, two_ptr, bad_row, perm, inverse) != 2)
        return 5;
    if (fsparse_native_order(2, 2, two_ptr, asymmetric_row, perm, inverse) != 2)
        return 6; /* asymmetric graph */
    const int duplicate_ptr[3] = {1, 3, 3}, duplicate_row[2] = {1, 1};
    if (fsparse_native_order(2, 2, duplicate_ptr, duplicate_row,
                            perm, inverse) != 2) return 9;
    if (fsparse_native_order_general(2, 2, duplicate_ptr, duplicate_row,
                                    perm, inverse) != 2) return 13;
    const int diagonal_ptr[4] = {1, 2, 3, 4}, diagonal_row[3] = {1, 2, 3};
    if (fsparse_native_order(3, 3, diagonal_ptr, diagonal_row,
                            perm, inverse) != 0 ||
        !check_inverse(perm, inverse, 3)) return 7;
    const int empty_ptr[4] = {1, 1, 1, 1};
    if (fsparse_native_order_general(3, 0, empty_ptr, NULL,
                                    perm, inverse) != 0 ||
        !check_inverse(perm, inverse, 3)) return 14;
#else
    (void)check_inverse; (void)again; (void)again_inverse;
    if (status != 3) return 8;
    if (fsparse_native_order_available() != 0 ||
        fsparse_native_order_general(5, 11, ptr, row, perm, inverse) != 3)
        return 15;
#endif
    puts("PASS native ordering permutation, inverse, determinism and status");
    return 0;
}
