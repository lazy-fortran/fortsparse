/* Behavioral interposition: actual wide-panel BLAS lowering, no timing claim. */
#include <stddef.h>
static int calls;
void __real_dtrsm_(const char *, const char *, const char *, const char *,
    const int *, const int *, const double *, const double *, const int *,
    double *, const int *, size_t, size_t, size_t, size_t);
void __wrap_dtrsm_(const char *side, const char *uplo, const char *trans,
    const char *diag, const int *m, const int *n, const double *alpha,
    const double *a, const int *lda, double *b, const int *ldb,
    size_t sl, size_t ul, size_t tl, size_t dl) {
    ++calls;
    __real_dtrsm_(side,uplo,trans,diag,m,n,alpha,a,lda,b,ldb,sl,ul,tl,dl);
}
int fsparse_test_dtrsm_calls(void) { return calls; }
