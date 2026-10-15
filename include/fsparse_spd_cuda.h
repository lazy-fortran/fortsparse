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
int fsparse_spd_cuda_destroy(void *context);
int fsparse_spd_cuda_upload(void *context, size_t offset,
                            const double *source, size_t count);
int fsparse_spd_cuda_download(void *context, size_t offset,
                              double *destination, size_t count);
int fsparse_spd_cuda_synchronize(void *context);
int fsparse_spd_cuda_panel_factor(void *context, size_t offset, int n, int lda);
int fsparse_spd_cuda_factor_info(void *context, int *info);
int fsparse_spd_cuda_panel_solve(void *context, size_t l_offset, int n, int ldl,
                                 size_t b_offset, int m, int ldb);
int fsparse_spd_cuda_schur_update(void *context, size_t c_offset, int ldc,
                                  size_t left_offset, int ldl,
                                  size_t right_offset, int ldr,
                                  int m, int n, int k);
size_t fsparse_spd_cuda_allocated_bytes(void *context);

#ifdef __cplusplus
}
#endif
#endif
