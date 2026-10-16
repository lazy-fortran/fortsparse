#ifndef FSPARSE_NATIVE_AMD_BOUNDED_H
#define FSPARSE_NATIVE_AMD_BOUNDED_H
#include <stdint.h>
/* MIT adapter; external AMD 3.3.4 uses caller workspace. No allocations here.
 * Inputs and outputs follow the existing one-based native ordering ABI. */
int fsparse_native_amd_bounded_available(void);
int fsparse_native_amd_bounded_size(int n, int nnz, const int *ptr, const int *row,
    int low, int64_t *workspace_bytes, int64_t *output_bytes);
int fsparse_native_amd_bounded(int n, int nnz, const int *ptr, const int *row,
    int low, int32_t *workspace, int64_t workspace_bytes, int *perm, int *inverse);
#endif
