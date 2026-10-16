/* MIT. Native single-process general LU, double precision. */
#ifndef FORTSPARSE_NATIVE_H
#define FORTSPARSE_NATIVE_H
#include <stdint.h>
#ifdef __cplusplus
#include <complex>
typedef std::complex<double> fs_native_complex;
extern "C" {
#else
#include <complex.h>
typedef double _Complex fs_native_complex;
#endif
/* Status codes match fortsparse_status: 0 success; 1 singular; 2 invalid;
 * 3 unavailable; 4 unfactored; 5 internal; 6 allocation; 7 index overflow;
 * 8 resource limit; 9 not positive definite. All objects are process-local.
 * create output must point to NULL. free accepts NULL and clears its argument.
 * Never use a copied handle after free, pass a foreign pointer, or concurrently
 * call one handle. Independent handles may be used independently.
 */
int32_t fs_native_create(int32_t scalar_kind, void **handle); /* 0 real, 1 complex */
int32_t fs_native_free(void **handle);
/* Changing refinement invalidates factors. Default false; applied on next factor.
 * Native ordering is DEFAULT (natural); no external factorization fallback.
 */
int32_t fs_native_set_refine(void *handle, int32_t refine);
enum fs_native_ordering {
    FS_NATIVE_ORDER_DEFAULT = 0,
    FS_NATIVE_ORDER_AMD = 1,
    FS_NATIVE_ORDER_ND = 2
};
/* Returns status 0 supported, 3 unavailable, 2 invalid enum. No context/error
 * mutation. ND is recognized but unavailable for this general-LU API.
 */
int32_t fs_native_ordering_available(int32_t ordering);
/* Default remains natural. Successful policy changes invalidate factors;
 * setting the same policy retains them. Invalid/unavailable requests preserve
 * policy and factors, and store the error; a valid setter clears that error.
 * AMD requires the optional ordering library, never an external factor backend.
 */
int32_t fs_native_set_ordering(void *handle, int32_t ordering);
enum fs_native_memory_policy {
    FS_NATIVE_MEMORY_AUTO = 0,
    FS_NATIVE_MEMORY_LOW = 1
};
/* Per-handle solver-owned array payload, including C conversion/output staging.
 * Zero limits mean unlimited. Negative limits/unknown policies preserve factors.
 * A successful change frees/invalidates factors; same policy retains them.
 * CPU general LU owns zero device bytes; the device cap is retained explicitly.
 * Allocator/runtime/object-descriptor overhead and caller arrays are excluded.
 */
int32_t fs_native_set_memory(void *, int64_t host_limit_bytes,
    int64_t device_limit_bytes, int32_t memory_policy);
struct fs_native_memory_stats {
    int64_t host_limit_bytes, device_limit_bytes;
    int64_t host_current_bytes, host_peak_bytes;
    int64_t device_current_bytes, device_peak_bytes;
    int64_t budget_required_bytes, budget_available_bytes;
    int32_t memory_policy, strict_memory_supported, peak_is_bound;
};
/* Pure diagnostic query: does not change last_error. Peak describes the last
 * factor/solve, not RSS or a lifetime maximum. peak_is_bound marks conservative
 * backend admission/history counters; -1 is unknown (opaque uncapped AMD).
 * Required/available describe the failed allocation phase, not final fill.
 */
int32_t fs_native_get_memory_stats(void *, struct fs_native_memory_stats *);
/* Canonical zero-based CSC: n>0, 0<=nnz<INT32_MAX, colptr[n+1] starts at 0,
 * ends at nnz, monotone; rowind[nnz] strictly increasing within each column.
 * Arrays must contain their declared number of accessible elements. Input is
 * borrowed only during the call; conversion and retained factors are owned.
 * Even nnz=0 requires non-NULL rowind/values pointers.
 * Any factor/refactor failure invalidates prior factors. A valid later call can
 * recover. refactor retains symbolic work when the backend supports that pattern.
 */
int32_t fs_native_factor_d(void *, int32_t n, int32_t nnz,
    const int32_t *colptr, const int32_t *rowind, const double *values);
int32_t fs_native_factor_z(void *, int32_t n, int32_t nnz,
    const int32_t *colptr, const int32_t *rowind, const fs_native_complex *values);
int32_t fs_native_refactor_d(void *, int32_t n, int32_t nnz,
    const int32_t *colptr, const int32_t *rowind, const double *values);
int32_t fs_native_refactor_z(void *, int32_t n, int32_t nnz,
    const int32_t *colptr, const int32_t *rowind, const fs_native_complex *values);
/* op=0 normal, 1 transpose, 2 conjugate transpose. Column-major n*nrhs, nrhs>0.
 * All input values must be finite. b/x may alias; output changes only on success.
 * Complex pointers use C double complex / C++ std::complex<double> array layout,
 * with Fortran ISO_C_BINDING c_double_complex; no complex value crosses by value.
 */
int32_t fs_native_solve_d(void *, int32_t op, int32_t n, int32_t nrhs,
    const double *b, double *x);
int32_t fs_native_solve_z(void *, int32_t op, int32_t n, int32_t nrhs,
    const fs_native_complex *b, fs_native_complex *x);
/* Copies last diagnostic, always NUL-terminates when capacity>0; returns last
 * status (or 2 for NULL handle/invalid output). Does not change stored status.
 */
int32_t fs_native_last_error(void *, char *message, int32_t capacity);
#ifdef __cplusplus
}
#endif
#endif
