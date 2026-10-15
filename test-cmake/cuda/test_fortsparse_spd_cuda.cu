#include "fsparse_spd_cuda.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>

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
    std::puts("CUDA resident factor, triangular solve, Schur update and errors pass");
    return 0;
}
