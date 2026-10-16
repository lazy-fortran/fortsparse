#include "fsparse_spd_cuda.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <initializer_list>
#include <limits>

namespace {
void check(int code)
{
    if (code) {
        std::fprintf(stderr, "native CUDA operation failed: %d\n", code);
        std::exit(1);
    }
}
void close(double value, double expected)
{
    if (!std::isfinite(value) || std::abs(value-expected) > 2.0e-13) {
        std::fprintf(stderr, "native CUDA mismatch: %.17g vs %.17g\n", value, expected);
        std::exit(1);
    }
}

void fused_oracle(void *ctx)
{
    // Different operand strides, nontrivial one-based maps and unused cells.
    // The independent oracle uses a direct host dot product, not emitted leaves.
    const int maps[8] = {1,3,0,8,7, 1,3,5};
    double left[56], right[32], initial[35], actual[35], expected[35];
    for (int p = 0; p < 8; ++p) {
        for (int i = 0; i < 7; ++i) left[i+7*p] = (i < 5) ? (i-2*p+1)*.125 : 999;
        for (int j = 0; j < 4; ++j) right[j+4*p] = (j < 3) ? (2*j+p-1)*.25 : 777;
    }
    for (int i = 0; i < 35; ++i) initial[i] = 10+i*.0625;
    check(fsparse_spd_cuda_set_maps(ctx, maps, 8));
    check(fsparse_spd_cuda_upload(ctx, 32, left, 56));
    check(fsparse_spd_cuda_upload(ctx, 96, right, 32));
    const size_t retained = fsparse_spd_cuda_allocated_bytes(ctx);
    for (int width : {0,1,3,8}) {
        check(fsparse_spd_cuda_begin_factor(ctx));
        check(fsparse_spd_cuda_upload(ctx, 160, initial, 35));
        // Reusing the destination applies two contractions; no reset/scratch.
        for (int repeat = 0; repeat < 2; ++repeat)
            check(fsparse_spd_cuda_schur_scatter(ctx, 32, 7, 96, 4, 5, 3,
                                                width, 160, 7, 5, 0, 5));
        check(fsparse_spd_cuda_download(ctx, 160, actual, 35));
        check(fsparse_spd_cuda_synchronize(ctx));
        for (int i = 0; i < 35; ++i) expected[i] = initial[i];
        for (int j = 0; j < 3; ++j)
            for (int i = j; i < 5; ++i) {
                if (maps[i] <= 0 || maps[i] > 7) continue;
                double product = 0;
                for (int p = 0; p < width; ++p) product += left[i+7*p]*right[j+4*p];
                expected[maps[i]-1+7*(maps[5+j]-1)] -= 2*product;
            }
        for (int i = 0; i < 35; ++i) close(actual[i], expected[i]);
        if (fsparse_spd_cuda_allocated_bytes(ctx) != retained) std::exit(1);
    }
    // Missing/invalid column entries also leave the entire destination intact.
    const int missing[8] = {1,3,0,-1,7, 0,-1,6};
    check(fsparse_spd_cuda_set_maps(ctx, missing, 8));
    check(fsparse_spd_cuda_upload(ctx, 160, initial, 35));
    check(fsparse_spd_cuda_schur_scatter(ctx, 32, 7, 96, 4, 5, 3, 8, 160, 7, 5, 0, 5));
    check(fsparse_spd_cuda_download(ctx, 160, actual, 35));
    check(fsparse_spd_cuda_synchronize(ctx));
    for (int i = 0; i < 35; ++i) close(actual[i], initial[i]);
    check(fsparse_spd_cuda_set_maps(ctx, maps, 8));
    if (fsparse_spd_cuda_schur_scatter(ctx, 32, 4, 96, 4, 5, 3, 8, 160, 7, 5, 0, 5) != -1 ||
        fsparse_spd_cuda_schur_scatter(ctx, 32, 7, 96, 4, 5, 3, -1, 160, 7, 5, 0, 5) != -1 ||
        fsparse_spd_cuda_schur_scatter(ctx, 32, 7, 96, 4, -1, 3, 8, 160, 7, 5, 0, 5) != -1 ||
        fsparse_spd_cuda_schur_scatter(ctx, 32, 7, 96, 4, 5, 3, 8, 240, 7, 5, 0, 5) != -1 ||
        fsparse_spd_cuda_schur_scatter(ctx, 32, 7, 96, 4, 5, 3, 8, 160, 7, 5, 4, 5) != -1 ||
        fsparse_spd_cuda_schur_scatter(ctx, 32, 7, 96, 4, 5, 3, 8, 160, 7, 5, 0, 6) != -1 ||
        fsparse_spd_cuda_schur_scatter(ctx, std::numeric_limits<size_t>::max(),
                                      7, 96, 4, 5, 3, 8, 160, 7, 5, 0, 5) != -1 ||
        fsparse_spd_cuda_schur_scatter(nullptr, 32, 7, 96, 4, 5, 3, 8, 160, 7, 5, 0, 5) != -1)
        std::exit(1);
    // A failed earlier panel suppresses the fused write in the retained stream.
    const double bad_pivot = -1;
    check(fsparse_spd_cuda_begin_factor(ctx));
    check(fsparse_spd_cuda_upload(ctx, 0, &bad_pivot, 1));
    check(fsparse_spd_cuda_panel_factor(ctx, 0, 1, 1));
    check(fsparse_spd_cuda_schur_scatter(ctx, 32, 7, 96, 4, 5, 3, 8, 160, 7, 5, 0, 5));
    int info = 0;
    check(fsparse_spd_cuda_factor_info(ctx, &info));
    if (info != 1) std::exit(1);
    check(fsparse_spd_cuda_download(ctx, 160, actual, 35));
    check(fsparse_spd_cuda_synchronize(ctx));
    for (int i = 0; i < 35; ++i) close(actual[i], initial[i]);
}
}

int main()
{
    void *ctx = nullptr;
    check(fsparse_spd_cuda_create(&ctx, 0, 256));
    if (fsparse_spd_cuda_allocated_bytes(ctx) < 256*sizeof(double)) return 1;
    /* Lower triangular oracle; arbitrary upper values are ignored. */
    const double l[16] = {2,.5,.25,.75, 0,3,.5,.25, 0,0,4,.5, 0,0,0,5};
    double matrix[16] = {}, factor[16];
    for (int j = 0; j < 4; ++j)
        for (int i = 0; i < 4; ++i)
            for (int k = 0; k < 4; ++k) matrix[i+4*j] += l[i+4*k]*l[j+4*k];
    check(fsparse_spd_cuda_upload(ctx, 0, matrix, 16));
    check(fsparse_spd_cuda_panel_factor(ctx, 0, 4, 4));
    int info = -1;
    check(fsparse_spd_cuda_factor_info(ctx, &info));
    if (info != 0) return 1;
    check(fsparse_spd_cuda_download(ctx, 0, factor, 16));
    check(fsparse_spd_cuda_synchronize(ctx));
    for (int j = 0; j < 4; ++j)
        for (int i = j; i < 4; ++i) close(factor[i+4*j], l[i+4*j]);

    const double rhs[12] = {1,2,3, -1,-2,-3, .5,1,2, 2,1,0};
    double solved[12];
    check(fsparse_spd_cuda_upload(ctx, 32, rhs, 12));
    check(fsparse_spd_cuda_panel_solve(ctx, 0, 4, 4, 32, 3, 3));
    check(fsparse_spd_cuda_download(ctx, 32, solved, 12));
    check(fsparse_spd_cuda_synchronize(ctx));
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 3; ++i) {
            double residual = -rhs[i+3*j];
            for (int k = 0; k < 4; ++k) residual += solved[i+3*k]*l[j+4*k];
            close(residual, 0.0);
        }
    }
    const double left[6] = {1,2,3, -1,.5,2};
    const double right[8] = {2,1,0,-1, .5,-2,1,3};
    double updated[12];
    check(fsparse_spd_cuda_upload(ctx, 64, left, 6));
    check(fsparse_spd_cuda_upload(ctx, 80, right, 8));
    check(fsparse_spd_cuda_upload(ctx, 96, rhs, 12));
    check(fsparse_spd_cuda_schur_update(ctx, 96, 3, 64, 3, 80, 4, 3, 4, 2));
    check(fsparse_spd_cuda_download(ctx, 96, updated, 12));
    check(fsparse_spd_cuda_synchronize(ctx));
    for (int j = 0; j < 4; ++j)
        for (int i = 0; i < 3; ++i) {
            double expected = rhs[i+3*j];
            for (int k = 0; k < 2; ++k) expected -= left[i+3*k]*right[j+4*k];
            close(updated[i+3*j], expected);
        }
    /* Reused product scratch overwrites its previous contents, unlike update. */
    for (int repeat = 0; repeat < 2; ++repeat) {
        check(fsparse_spd_cuda_schur_product(ctx, 96, 3, 64, 3, 80, 4, 3, 4, 2));
        check(fsparse_spd_cuda_download(ctx, 96, updated, 12));
        check(fsparse_spd_cuda_synchronize(ctx));
        for (int j = 0; j < 4; ++j)
            for (int i = 0; i < 3; ++i) {
                double expected = 0;
                for (int k = 0; k < 2; ++k) expected -= left[i+3*k]*right[j+4*k];
                close(updated[i+3*j], expected);
            }
    }

    fused_oracle(ctx);
    matrix[0] = -1;
    check(fsparse_spd_cuda_upload(ctx, 0, matrix, 16));
    check(fsparse_spd_cuda_panel_factor(ctx, 0, 4, 4));
    check(fsparse_spd_cuda_factor_info(ctx, &info));
    if (info != 1) return 1;
    if (fsparse_spd_cuda_upload(ctx, 250, matrix, 16) != -1) return 1;
    if (fsparse_spd_cuda_panel_factor(ctx, 240, 4, 8) != -1) return 1;
    if (fsparse_spd_cuda_schur_update(ctx, 96, 2, 64, 3, 80, 4, 3, 4, 2) != -1)
        return 1;
    check(fsparse_spd_cuda_destroy(ctx));
    std::puts("CUDA resident factor, triangular solve, Schur/scatter updates and errors pass");
    return 0;
}
