#include "fsparse_spd_cuda.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <initializer_list>
#include <limits>
#include <vector>

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

void batched_oracle()
{
    void *ctx = nullptr;
    check(fsparse_spd_cuda_create(&ctx, 0, 2048));
    std::vector<double> initial(2048), actual(2048), expected(2048);
    for (size_t i = 0; i < initial.size(); ++i) initial[i] = (int(i%19)-9)*.03125;
    int maps[61];
    for (int i = 0; i < 37; ++i) maps[i] = i+1;
    maps[3] = 0; maps[7] = 42;
    const int columns[9] = {1,3,4,6,7,8,9,10,12};
    for (int j = 0; j < 9; ++j) maps[37+j] = columns[j];
    for (int i = 0; i < 11; ++i) maps[46+i] = i+2;
    for (int j = 0; j < 4; ++j) maps[57+j] = 2*j+1;
    maps[59] = -1;
    check(fsparse_spd_cuda_set_maps(ctx, maps, 61));
    fsparse_spd_cuda_schur_task tasks[4] = {
        {32,200,1024,0,37,40,12,37,9,3,40,12},
        {32,200,1024,0,37,40,12,37,9,0,40,12}, // zero work in prefix schedule
        {400,460,1600,46,57,13,5,11,4,2,14,8},
        {32,200,1024,0,37,40,12,37,9,3,40,12}
    };
    const size_t ptr[4] = {0,0,3,4}; // empty batch, two destinations, then reuse
    const size_t device_before = fsparse_spd_cuda_allocated_bytes(ctx);
    const size_t host_before = fsparse_spd_cuda_host_plan_bytes(ctx);
    check(fsparse_spd_cuda_set_schur_batches(ctx, tasks, 4, ptr, 3));
    const size_t device_bytes = fsparse_spd_cuda_allocated_bytes(ctx);
    const size_t host_bytes = fsparse_spd_cuda_host_plan_bytes(ctx);
    if (device_bytes <= device_before || host_bytes == 0) std::exit(1);
    auto apply = [&](const fsparse_spd_cuda_schur_task &t) {
        for (int j = 0; j < t.n; ++j)
            for (int i = j; i < t.m; ++i) {
                int row = maps[t.row_map_offset+i], col = maps[t.column_map_offset+j];
                if (row <= 0 || row > t.destination_rows ||
                    col <= 0 || col > t.destination_columns) continue;
                double product = 0;
                for (int p = 0; p < t.k; ++p)
                    product += initial[t.left_offset+i+size_t(p)*t.ldl]*
                               initial[t.right_offset+j+size_t(p)*t.ldr];
                expected[t.destination_offset+row-1+size_t(col-1)*t.destination_rows] -= product;
            }
    };
    for (int repeat = 0; repeat < 3; ++repeat) {
        expected = initial;
        check(fsparse_spd_cuda_begin_factor(ctx));
        check(fsparse_spd_cuda_upload(ctx, 0, initial.data(), initial.size()));
        for (size_t batch : {size_t(0),size_t(1),size_t(1),size_t(2)})
            check(fsparse_spd_cuda_schur_batch(ctx, batch));
        apply(tasks[0]); apply(tasks[2]); apply(tasks[0]); apply(tasks[2]); apply(tasks[3]);
        check(fsparse_spd_cuda_download(ctx, 0, actual.data(), actual.size()));
        check(fsparse_spd_cuda_synchronize(ctx));
        for (size_t i = 0; i < actual.size(); ++i) close(actual[i], expected[i]);
        if (fsparse_spd_cuda_allocated_bytes(ctx) != device_bytes ||
            fsparse_spd_cuda_host_plan_bytes(ctx) != host_bytes) std::exit(1);
    }
    // Invalid setup is rejected before replacing the previous valid plan.
    auto saved = tasks[2];
    for (int error = 0; error < 6; ++error) {
        tasks[2] = saved;
        if (error == 0) tasks[2].destination_offset = 1024; // overlapping destinations
        if (error == 1) tasks[2].left_offset = 1700; // read/write alias across tasks
        if (error == 2) tasks[2].row_map_offset = 60;
        if (error == 3) tasks[2].left_offset = std::numeric_limits<size_t>::max();
        if (error == 4) tasks[2].ldl = 10;
        if (error == 5) tasks[2].k = 9;
        if (fsparse_spd_cuda_set_schur_batches(ctx, tasks, 4, ptr, 3) != -1) std::exit(1);
    }
    tasks[2] = saved;
    const size_t reversed[4] = {0,3,2,4};
    if (fsparse_spd_cuda_set_schur_batches(ctx, tasks, 4, reversed, 3) != -1 ||
        fsparse_spd_cuda_schur_batch(ctx, 3) != -1 ||
        fsparse_spd_cuda_schur_batch(nullptr, 0) != -1) std::exit(1);
    const double bad = -1;
    check(fsparse_spd_cuda_begin_factor(ctx));
    check(fsparse_spd_cuda_upload(ctx, 0, initial.data(), initial.size()));
    check(fsparse_spd_cuda_upload(ctx, 0, &bad, 1));
    check(fsparse_spd_cuda_panel_factor(ctx, 0, 1, 1));
    check(fsparse_spd_cuda_schur_batch(ctx, 1));
    int info = 0;
    check(fsparse_spd_cuda_factor_info(ctx, &info));
    if (info != 1) std::exit(1);
    check(fsparse_spd_cuda_download(ctx, 0, actual.data(), actual.size()));
    check(fsparse_spd_cuda_synchronize(ctx));
    for (size_t i = 1; i < actual.size(); ++i) close(actual[i], initial[i]);
    // Replacing maps clears stale plans and releases both metadata allocations.
    check(fsparse_spd_cuda_set_maps(ctx, maps, 61));
    if (fsparse_spd_cuda_schur_batch(ctx, 0) != -1 ||
        fsparse_spd_cuda_host_plan_bytes(ctx) != host_before ||
        fsparse_spd_cuda_allocated_bytes(ctx) != device_before) std::exit(1);
    check(fsparse_spd_cuda_destroy(ctx));
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
    batched_oracle();
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
