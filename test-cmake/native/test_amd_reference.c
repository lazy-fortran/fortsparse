/* The native AMD ordering equals SuiteSparse amd_order's permutation, which
 * the test loads at run time where it is installed (exit 77 skips otherwise).
 * The bounded low-memory and 64-bit paths order alike. With the argument
 * "bench" it times both orderings on larger grids. */
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "fsparse_native_amd_bounded.h"

int fsparse_native_order_amd_general(int, int, const int *, const int *, int *, int *);
int fsparse_native_order_amd_general_wide(int, int, const int *, const int *, int *, int *);
typedef int (*amd_order_t)(int32_t, const int32_t *, const int32_t *, int32_t *,
                           double *, double *);

static int failures;
static void require(int ok, const char *what, const char *name) {
    if (!ok) { fprintf(stderr, "FAIL %s: %s\n", name, what); ++failures; }
}
static double now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

struct csc { int n, nnz, *ptr, *row; };
static uint64_t seed = 88172645463325252ull;
static uint64_t xorshift(void) {
    seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17;
    return seed;
}
/* One-based CSC from a dense-free generator: column j holds the rows for
 * which keep(i, j) is true, in increasing order. */
static struct csc stencil(int nx, int ny, int nz, int r, int full) {
    struct csc a;
    int n = nx * ny * nz, cap = (2 * r + 1) * (2 * r + 1) * (2 * r + 1), k = 0;
    a.n = n;
    a.ptr = malloc((n + 1) * sizeof(int));
    a.row = malloc((size_t)n * cap * sizeof(int));
    for (int z = 0; z < nz; z++) for (int y = 0; y < ny; y++) for (int x = 0; x < nx; x++) {
        int j = x + nx * (y + ny * z);
        a.ptr[j] = k + 1;
        for (int dz = -r; dz <= r; dz++) for (int dy = -r; dy <= r; dy++) for (int dx = -r; dx <= r; dx++) {
            if (!full && abs(dx) + abs(dy) + abs(dz) > 1) continue;
            int X = x + dx, Y = y + dy, Z = z + dz;
            if (X < 0 || Y < 0 || Z < 0 || X >= nx || Y >= ny || Z >= nz) continue;
            a.row[k++] = X + nx * (Y + ny * Z) + 1;
        }
    }
    a.ptr[n] = k + 1;
    a.nnz = k;
    return a;
}
/* Unsymmetric random pattern with per-column density, optional dense
 * columns and rows, and empty columns. */
static struct csc random_pattern(int n, int per_column, int dense, int empty) {
    struct csc a;
    char *mark = calloc(n, 1);
    a.n = n;
    a.ptr = malloc((n + 1) * sizeof(int));
    a.row = malloc((size_t)n * (per_column + 2) * sizeof(int) + (size_t)dense * n * sizeof(int) * 2);
    int k = 0;
    for (int j = 0; j < n; j++) {
        a.ptr[j] = k + 1;
        if (j % 97 == 3 && empty) continue;
        int count = (j < dense) ? n : per_column;
        for (int c = 0; c < count; c++) {
            int i = (j < dense) ? c : (int)(xorshift() % (uint64_t)n);
            if (count < n) mark[i] = 1;
        }
        if (count < n) {
            mark[j] = 1;
            for (int d = 0; d < dense; d++) mark[n - 1 - d] = 1;
            for (int i = 0; i < n; i++) if (mark[i]) { a.row[k++] = i + 1; mark[i] = 0; }
        } else {
            for (int i = 0; i < n; i++) a.row[k++] = i + 1;
        }
    }
    a.ptr[n] = k + 1;
    a.nnz = k;
    free(mark);
    return a;
}

static void compare(amd_order_t amd_order, struct csc a, const char *name, int bench) {
    int n = a.n;
    int *perm = malloc(n * sizeof(int)), *inverse = malloc(n * sizeof(int));
    int *perm2 = malloc(n * sizeof(int)), *inverse2 = malloc(n * sizeof(int));
    int32_t *ap = malloc((n + 1) * sizeof(int32_t)), *ai = malloc((a.nnz ? a.nnz : 1) * sizeof(int32_t));
    int32_t *p = malloc(n * sizeof(int32_t));
    for (int j = 0; j <= n; j++) ap[j] = a.ptr[j] - 1;
    for (int q = 0; q < a.nnz; q++) ai[q] = a.row[q] - 1;
    double t_native = 1e30, t_ref = 1e30;
    int reps = bench ? 5 : 1;
    for (int rep = 0; rep < reps; rep++) {
        double t = now();
        require(fsparse_native_order_amd_general(n, a.nnz, a.ptr, a.row, perm, inverse) == 0,
                "native order", name);
        t = now() - t;
        if (t < t_native) t_native = t;
        t = now();
        require(amd_order(n, ap, ai, p, NULL, NULL) >= 0, "SuiteSparse order", name);
        t = now() - t;
        if (t < t_ref) t_ref = t;
    }
    int same = 1;
    for (int k = 0; k < n; k++) same &= perm[k] == p[k] + 1;
    require(same, "permutation equals amd_order", name);
    require(fsparse_native_order_amd_general_wide(n, a.nnz, a.ptr, a.row, perm2, inverse2) == 0,
            "64-bit order", name);
    same = 1;
    for (int k = 0; k < n; k++) same &= perm2[k] == perm[k];
    require(same, "64-bit permutation equals 32-bit", name);
    int64_t bytes = 0, outputs = 0;
    require(fsparse_native_amd_bounded_size(n, a.nnz, a.ptr, a.row, 1, &bytes, &outputs) == 0,
            "low-memory size", name);
    int32_t *arena = malloc((size_t)bytes);
    require(fsparse_native_amd_bounded(n, a.nnz, a.ptr, a.row, 1, arena, bytes, perm2, inverse2) == 0,
            "low-memory order", name);
    same = 1;
    for (int k = 0; k < n; k++) same &= perm2[k] == perm[k];
    require(same, "low-memory permutation equals default", name);
    if (bench)
        printf("%-28s n=%8d nnz=%9d native %.4f s SuiteSparse %.4f s ratio %.2f\n",
               name, n, a.nnz, t_native, t_ref, t_native / t_ref);
    free(arena); free(perm); free(inverse); free(perm2); free(inverse2);
    free(ap); free(ai); free(p);
    free(a.ptr); free(a.row);
}

int main(int argc, char **argv) {
    const char *names[] = {"libamd.so.3", "libamd.so", "libamd.dylib", NULL};
    void *library = NULL;
    for (int i = 0; names[i] && !library; i++) library = dlopen(names[i], RTLD_NOW);
    if (!library) {
        printf("SKIP: SuiteSparse AMD is not installed\n");
        return 77;
    }
    amd_order_t amd_order = (amd_order_t)dlsym(library, "amd_order");
    if (!amd_order) {
        printf("SKIP: libamd has no amd_order\n");
        return 77;
    }
    int bench = argc > 1;
    if (bench) {
        compare(amd_order, stencil(400, 400, 1, 1, 0), "2D 5-point 400^2", 1);
        compare(amd_order, stencil(400, 400, 1, 1, 1), "2D 9-point 400^2", 1);
        compare(amd_order, stencil(50, 50, 50, 1, 0), "3D 7-point 50^3", 1);
        compare(amd_order, stencil(36, 36, 36, 1, 1), "3D 27-point 36^3", 1);
        compare(amd_order, stencil(250, 250, 1, 2, 1), "2D 25-point 250^2", 1);
        compare(amd_order, random_pattern(100000, 6, 0, 1), "random unsymmetric 1e5", 1);
    } else {
        compare(amd_order, stencil(1, 1, 1, 1, 0), "single unknown", 0);
        compare(amd_order, stencil(60, 60, 1, 1, 0), "2D 5-point 60^2", 0);
        compare(amd_order, stencil(60, 60, 1, 1, 1), "2D 9-point 60^2", 0);
        compare(amd_order, stencil(14, 14, 14, 1, 0), "3D 7-point 14^3", 0);
        compare(amd_order, stencil(12, 12, 12, 1, 1), "3D 27-point 12^3", 0);
        compare(amd_order, stencil(40, 40, 1, 2, 1), "2D 25-point 40^2", 0);
        compare(amd_order, random_pattern(3000, 5, 0, 1), "random with empty columns", 0);
        compare(amd_order, random_pattern(3000, 4, 3, 0), "random with dense rows", 0);
        compare(amd_order, random_pattern(400, 60, 0, 0), "random dense-ish", 0);
    }
    if (failures) return 1;
    printf("PASS: native AMD equals SuiteSparse amd_order\n");
    return 0;
}
