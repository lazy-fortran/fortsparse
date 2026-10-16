#include "fsparse_spd_cuda.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>

namespace {
void check(int code)
{
    if (!code) return;
    std::fprintf(stderr, "thin-panel CUDA operation failed: %d\n", code);
    std::exit(1);
}
void close(double actual, double expected)
{
    if (std::isfinite(actual) &&
        std::abs(actual-expected) <= 3e-13*std::fmax(1.0, std::abs(expected))) return;
    std::fprintf(stderr, "thin-panel mismatch: %.17g vs %.17g\n", actual, expected);
    std::exit(1);
}
void require(bool value)
{
    if (!value) std::exit(1);
}

void manufactured_front(void *ctx, int width, int rows, int padding)
{
    const int lda = rows+padding;
    const int columns = width ? width : 1;
    std::vector<double> lower(lda*columns), a(lda*columns), actual(a.size());
    const size_t retained = fsparse_spd_cuda_allocated_bytes(ctx);
    for (int repeat = 0; repeat < 2; ++repeat) {
        for (int j = 0; j < columns; ++j)
            for (int i = 0; i < lda; ++i)
                lower[i+lda*j] = i < j ? 0.0 : i == j ? 2.0+j*.25 :
                    .015625*(1+(i+3*j+repeat)%11);
        // Construct A11=L11*L11^T and A21=L21*L11^T independently.
        // Upper entries and leading-dimension padding are sentinels.
        for (int j = 0; j < columns; ++j)
            for (int i = 0; i < lda; ++i) {
                a[i+lda*j] = 9000+i+j;
                if (j >= width || i < j || i >= rows) continue;
                a[i+lda*j] = 0;
                for (int k = 0; k <= j; ++k)
                    a[i+lda*j] += lower[i+lda*k]*lower[j+lda*k];
            }
        check(fsparse_spd_cuda_upload(ctx, 0, a.data(), a.size()));
        check(fsparse_spd_cuda_begin_factor(ctx));
        // No vendor setup or device allocation is permitted in a factor batch.
        require(fsparse_spd_cuda_prepare_blas(ctx) == -1);
        check(fsparse_spd_cuda_front_factor(ctx, 0, width, rows, lda));
        int info = -1;
        check(fsparse_spd_cuda_factor_info(ctx, &info));
        require(info == 0);
        check(fsparse_spd_cuda_download(ctx, 0, actual.data(), actual.size()));
        check(fsparse_spd_cuda_synchronize(ctx));
        for (int j = 0; j < columns; ++j)
            for (int i = 0; i < lda; ++i)
                close(actual[i+lda*j], j < width && i >= j && i < rows ?
                      lower[i+lda*j] : a[i+lda*j]);
        require(fsparse_spd_cuda_allocated_bytes(ctx) == retained);
    }
}

void failed_pivot(void *ctx)
{
    const double bad[4] = {4, 0, 999, -1};
    const double untouched[4] = {9, 3, 999, 10};
    double actual[4];
    check(fsparse_spd_cuda_upload(ctx, 0, bad, 4));
    check(fsparse_spd_cuda_upload(ctx, 16, untouched, 4));
    check(fsparse_spd_cuda_begin_factor(ctx));
    check(fsparse_spd_cuda_front_factor(ctx, 0, 2, 2, 2));
    check(fsparse_spd_cuda_front_factor(ctx, 16, 1, 2, 2));
    int info = 0;
    check(fsparse_spd_cuda_factor_info(ctx, &info));
    require(info == 2);
    check(fsparse_spd_cuda_download(ctx, 16, actual, 4));
    check(fsparse_spd_cuda_synchronize(ctx));
    for (int i = 0; i < 4; ++i) close(actual[i], untouched[i]);
    // The next refactor batch clears the retained failure and solves normally.
    check(fsparse_spd_cuda_begin_factor(ctx));
    check(fsparse_spd_cuda_front_factor(ctx, 16, 1, 2, 2));
    check(fsparse_spd_cuda_factor_info(ctx, &info));
    require(info == 0);
    check(fsparse_spd_cuda_download(ctx, 16, actual, 4));
    check(fsparse_spd_cuda_synchronize(ctx));
    close(actual[0], 3);
    close(actual[1], 1);
}

void explicit_blas(void *ctx)
{
    const double lower[4] = {2,.5,999,3};
    const double rhs[4] = {1,2,3,4};
    const double left[4] = {1,2,3,4};
    double actual[4];
    check(fsparse_spd_cuda_upload(ctx, 0, lower, 4));
    check(fsparse_spd_cuda_upload(ctx, 8, rhs, 4));
    check(fsparse_spd_cuda_upload(ctx, 16, left, 4));
    const size_t before = fsparse_spd_cuda_allocated_bytes(ctx);
    require(fsparse_spd_cuda_panel_solve(ctx, 0, 2, 2, 8, 2, 2) == -1);
    require(fsparse_spd_cuda_schur_product(ctx, 24, 2, 16, 2, 16, 2, 2, 2, 2) == -1);
    require(fsparse_spd_cuda_allocated_bytes(ctx) == before);
    check(fsparse_spd_cuda_prepare_blas(ctx));
    const size_t prepared = fsparse_spd_cuda_allocated_bytes(ctx);
    require(prepared > before);
    check(fsparse_spd_cuda_prepare_blas(ctx));
    require(fsparse_spd_cuda_allocated_bytes(ctx) == prepared);
    check(fsparse_spd_cuda_panel_solve(ctx, 0, 2, 2, 8, 2, 2));
    check(fsparse_spd_cuda_download(ctx, 8, actual, 4));
    check(fsparse_spd_cuda_synchronize(ctx));
    for (int i = 0; i < 2; ++i) {
        close(actual[i]*2, rhs[i]);
        close(actual[i]*.5+actual[i+2]*3, rhs[i+2]);
    }
    check(fsparse_spd_cuda_schur_product(ctx, 24, 2, 16, 2, 16, 2, 2, 2, 2));
    check(fsparse_spd_cuda_schur_update(ctx, 24, 2, 16, 2, 16, 2, 2, 2, 2));
    check(fsparse_spd_cuda_download(ctx, 24, actual, 4));
    check(fsparse_spd_cuda_synchronize(ctx));
    for (int j = 0; j < 2; ++j)
        for (int i = 0; i < 2; ++i) {
            double product = 0;
            for (int k = 0; k < 2; ++k) product += left[i+2*k]*left[j+2*k];
            close(actual[i+2*j], -2*product);
        }
    require(fsparse_spd_cuda_allocated_bytes(ctx) == prepared);
}
}

int main()
{
    void *ctx = nullptr;
    constexpr size_t capacity = 9000;
    check(fsparse_spd_cuda_create_deferred(&ctx, 0, capacity));
    require(fsparse_spd_cuda_allocated_bytes(ctx) == capacity*sizeof(double)+sizeof(int));
    manufactured_front(ctx, 0, 300, 5);
    manufactured_front(ctx, 1, 513, 3);
    manufactured_front(ctx, 3, 300, 7);
    manufactured_front(ctx, 8, 1031, 6);
    failed_pivot(ctx);
    require(fsparse_spd_cuda_front_factor(ctx, 0, 9, 10, 10) == -1);
    require(fsparse_spd_cuda_front_factor(ctx, 0, 3, 2, 3) == -1);
    require(fsparse_spd_cuda_front_factor(ctx, 0, 3, 300, 299) == -1);
    require(fsparse_spd_cuda_front_factor(nullptr, 0, 1, 1, 1) == -1);
    require(fsparse_spd_cuda_front_factor(ctx, std::numeric_limits<size_t>::max(),
                                         1, 1, 1) == -1);
    explicit_blas(ctx);
    check(fsparse_spd_cuda_destroy(ctx));
    std::puts("CUDA thin front values, padding, refactor, pivot failure and setup pass");
}
