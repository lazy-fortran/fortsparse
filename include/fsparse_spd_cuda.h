#ifndef FSPARSE_SPD_CUDA_H
#define FSPARSE_SPD_CUDA_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

/* Single-device, single-stream retained arena; offsets count doubles.
 * Calls enqueue work; synchronize before reusing host upload/download memory.
 * Return 0 on success, -1 for invalid arguments, positive CUDA codes, or
 * 1000 + cuBLAS status. Device pivot info is 0 or a one-based failed pivot.
 * Layout is column-major. Overlapping operands are unsupported.
 */
int fsparse_spd_cuda_create(void **context, int device, size_t capacity);
/* Setup-only thin path: create without a cuBLAS handle/workspace. Explicit
 * BLAS operations require prepare_blas before begin_factor; they never allocate
 * it implicitly. Existing create retains its eager preparation semantics. */
int fsparse_spd_cuda_create_deferred(void **context, int device, size_t capacity);
int fsparse_spd_cuda_prepare_blas(void *context);
int fsparse_spd_cuda_destroy(void *context);
int fsparse_spd_cuda_upload(void *context, size_t offset,
                            const double *source, size_t count);
int fsparse_spd_cuda_download(void *context, size_t offset,
                              double *destination, size_t count);
int fsparse_spd_cuda_synchronize(void *context);
int fsparse_spd_cuda_panel_factor(void *context, size_t offset, int n, int lda);
/* Width 0..8 lower panel plus its trailing solve in one cooperative block.
 * rows>=width, lda>=max(1,rows); upper entries and padding are preserved. */
int fsparse_spd_cuda_front_factor(void *context, size_t offset,
                                  int width, int rows, int lda);
int fsparse_spd_cuda_factor_info(void *context, int *info);
/* Setup only: retain compact integer row/column maps on the device. */
int fsparse_spd_cuda_set_maps(void *context, const int *maps, size_t count);
/* Begin a whole factorization; pivot failure then survives every panel until
 * factor_info synchronizes and ends this batch. */
int fsparse_spd_cuda_begin_factor(void *context);
int fsparse_spd_cuda_scatter_update(void *context, size_t update_offset,
                                    int m, int n, size_t destination_offset,
                                    int destination_rows, int destination_columns,
                                    size_t row_map_offset, size_t column_map_offset);
/* Fused C(mapped rows, mapped columns) -= left right^T for source i >= j.
 * Maps are one-based; absent/out-of-range entries are skipped. Positive row
 * and column entries must be unique, as in the analyzed sparse update plan.
 * No product scratch is required. */
int fsparse_spd_cuda_schur_scatter(void *context, size_t left_offset, int ldl,
                                   size_t right_offset, int ldr, int m, int n, int k,
                                   size_t destination_offset, int destination_rows,
                                   int destination_columns, size_t row_map_offset,
                                   size_t column_map_offset);
/* Setup-only plan for thin (0 <= k <= 8) updates. Batch pointers are
 * zero-based task offsets, length batch_count+1. In each batch, nonempty
 * destinations must be ordered, disjoint full panel slices; all operand
 * reads must precede the first destination. Unique maps have the same
 * contract as schur_scatter. These conditions make one source panel safe
 * without atomics; batches remain stream ordered. set_maps invalidates this
 * plan. Empty batches/tasks enqueue no work. No factor-time allocation. */
typedef struct fsparse_spd_cuda_schur_task {
    size_t left_offset, right_offset, destination_offset;
    size_t row_map_offset, column_map_offset;
    int ldl, ldr, m, n, k, destination_rows, destination_columns;
} fsparse_spd_cuda_schur_task;
int fsparse_spd_cuda_set_schur_batches(void *context,
    const fsparse_spd_cuda_schur_task *tasks, size_t task_count,
    const size_t *batch_ptr, size_t batch_count);
int fsparse_spd_cuda_schur_batch(void *context, size_t batch);
/* Retained host context and launch metadata, separate from device bytes. */
size_t fsparse_spd_cuda_host_plan_bytes(void *context);
int fsparse_spd_cuda_panel_solve(void *context, size_t l_offset, int n, int ldl,
                                 size_t b_offset, int m, int ldb);
int fsparse_spd_cuda_schur_update(void *context, size_t c_offset, int ldc,
                                  size_t left_offset, int ldl,
                                  size_t right_offset, int ldr,
                                  int m, int n, int k);
/* Overwrite scratch C := -left right^T; retained scratch need not be zero. */
int fsparse_spd_cuda_schur_product(void *context, size_t c_offset, int ldc,
                                   size_t left_offset, int ldl,
                                   size_t right_offset, int ldr,
                                   int m, int n, int k);
size_t fsparse_spd_cuda_allocated_bytes(void *context);

#ifdef __cplusplus
}
#endif
#endif
