/* Behavioral interposition: actual wide-panel BLAS lowering, no timing claim.
   FORTNUM_DTRSM is the linker symbol of fortnum_blas3::fortnum_dtrsm, the
   entry point fortsparse calls; character lengths follow the reference
   arguments. */
#include <stddef.h>
#define WRAP_(name) __wrap_##name
#define REAL_(name) __real_##name
#define WRAP(name) WRAP_(name)
#define REAL(name) REAL_(name)
static int calls;
void REAL(FORTNUM_DTRSM)(const char *, const char *, const char *, const char *,
    const int *, const int *, const double *, const double *, const int *,
    double *, const int *, size_t, size_t, size_t, size_t);
void WRAP(FORTNUM_DTRSM)(const char *side, const char *uplo, const char *trans,
    const char *diag, const int *m, const int *n, const double *alpha,
    const double *a, const int *lda, double *b, const int *ldb,
    size_t sl, size_t ul, size_t tl, size_t dl) {
    __atomic_add_fetch(&calls, 1, __ATOMIC_RELAXED);
    REAL(FORTNUM_DTRSM)(side,uplo,trans,diag,m,n,alpha,a,lda,b,ldb,sl,ul,tl,dl);
}
int fsparse_test_dtrsm_calls(void) { return __atomic_load_n(&calls, __ATOMIC_RELAXED); }
