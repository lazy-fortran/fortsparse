/* Thin C shim over sequential SuperLU (BSD). Exposes a clean
 * factor/solve/free C ABI that Fortran binds via iso_c_binding. The
 * factorization (L, U, perm_c, perm_r, the column-permuted SuperMatrix and the
 * original CSC SuperMatrix) is held in an opaque handle so repeated solves
 * reuse it. License of this file: MIT (our code). It links only SuperLU. */

#include <stdlib.h>
#include <string.h>

#include <superlu/slu_ddefs.h>
#include <superlu/slu_zdefs.h>

void *fsparse_slu_factor_d(long n, long nnz, const long *colptr1,
                           const long *rowidx1, const double *val, int *info);
int fsparse_slu_solve_d(void *handle, const double *b, double *x, long n);
void fsparse_slu_free_d(void *handle);
void *fsparse_slu_factor_z(long n, long nnz, const long *colptr1,
                           const long *rowidx1, const double *val, int *info);
int fsparse_slu_solve_z(void *handle, const double *b, double *x, long n);
void fsparse_slu_free_z(void *handle);

/* Persistent real factorization. */
typedef struct {
    SuperMatrix A;  /* original CSC, SLU_NC */
    SuperMatrix AC; /* column-permuted, SLU_NCP */
    SuperMatrix L;  /* L factor, SLU_SC */
    SuperMatrix U;  /* U factor, SLU_NC */
    int *perm_c;
    int *perm_r;
    int *etree;
    int n;
    int_t *colptr0; /* 0-based CSC, owned by A */
    int_t *rowidx0;
    double *aval;
    superlu_options_t options;
    GlobalLU_t glu;
} slu_real_t;

/* Persistent complex factorization. */
typedef struct {
    SuperMatrix A;
    SuperMatrix AC;
    SuperMatrix L;
    SuperMatrix U;
    int *perm_c;
    int *perm_r;
    int *etree;
    int n;
    int_t *colptr0;
    int_t *rowidx0;
    doublecomplex *aval;
    superlu_options_t options;
    GlobalLU_t glu;
} slu_cplx_t;

/* Copy a 1-based long CSC index array into a freshly malloc'd 0-based int_t
 * array. Returns NULL on allocation failure. */
static int_t *copy_to_int_t0(const long *src, long count)
{
    int_t *dst = (int_t *) malloc((size_t) count * sizeof(int_t));
    long k;
    if (dst == NULL) return NULL;
    for (k = 0; k < count; ++k) dst[k] = (int_t) (src[k] - 1);
    return dst;
}

void *fsparse_slu_factor_d(long n, long nnz, const long *colptr1,
                           const long *rowidx1, const double *val, int *info)
{
    slu_real_t *h;
    SuperMatrix B; /* dummy, unused by dgstrf */
    SuperLUStat_t stat;
    int panel_size, relax;
    int_t lu_info = 0;

    *info = 0;
    h = (slu_real_t *) calloc(1, sizeof(slu_real_t));
    if (h == NULL) { *info = -1; return NULL; }

    h->n = (int) n;
    h->colptr0 = copy_to_int_t0(colptr1, n + 1);
    h->rowidx0 = copy_to_int_t0(rowidx1, nnz);
    h->aval = (double *) malloc((size_t) nnz * sizeof(double));
    h->perm_c = (int *) malloc((size_t) n * sizeof(int));
    h->perm_r = (int *) malloc((size_t) n * sizeof(int));
    h->etree = (int *) malloc((size_t) n * sizeof(int));
    if (h->colptr0 == NULL || h->rowidx0 == NULL || h->aval == NULL ||
        h->perm_c == NULL || h->perm_r == NULL || h->etree == NULL) {
        *info = -1;
        fsparse_slu_free_d(h);
        return NULL;
    }
    memcpy(h->aval, val, (size_t) nnz * sizeof(double));

    dCreate_CompCol_Matrix(&h->A, (int) n, (int) n, (int_t) nnz, h->aval,
                           h->rowidx0, h->colptr0, SLU_NC, SLU_D, SLU_GE);

    set_default_options(&h->options);
    h->options.ColPerm = COLAMD;
    get_perm_c(COLAMD, &h->A, h->perm_c);
    sp_preorder(&h->options, &h->A, h->perm_c, h->etree, &h->AC);

    panel_size = sp_ienv(1);
    relax = sp_ienv(2);
    StatInit(&stat);
    memset(&B, 0, sizeof(B));
    dgstrf(&h->options, &h->AC, relax, panel_size, h->etree, NULL, 0,
           h->perm_c, h->perm_r, &h->L, &h->U, &h->glu, &stat, &lu_info);
    StatFree(&stat);

    if (lu_info != 0) {
        /* lu_info > 0: zero pivot at column lu_info (singular). < 0: bad arg. */
        *info = (lu_info > 0) ? 1 : -1;
        fsparse_slu_free_d(h);
        return NULL;
    }
    return h;
}

int fsparse_slu_solve_d(void *handle, const double *b, double *x, long n)
{
    slu_real_t *h = (slu_real_t *) handle;
    SuperMatrix B;
    SuperLUStat_t stat;
    int info = 0;

    if (h == NULL) return -1;
    memcpy(x, b, (size_t) n * sizeof(double));
    dCreate_Dense_Matrix(&B, (int) n, 1, x, (int) n, SLU_DN, SLU_D, SLU_GE);
    StatInit(&stat);
    dgstrs(NOTRANS, &h->L, &h->U, h->perm_c, h->perm_r, &B, &stat, &info);
    StatFree(&stat);
    Destroy_SuperMatrix_Store(&B);
    return info;
}

void fsparse_slu_free_d(void *handle)
{
    slu_real_t *h = (slu_real_t *) handle;
    if (h == NULL) return;
    if (h->A.Store != NULL) Destroy_SuperMatrix_Store(&h->A);
    if (h->AC.Store != NULL) Destroy_CompCol_Permuted(&h->AC);
    if (h->L.Store != NULL) Destroy_SuperNode_Matrix(&h->L);
    if (h->U.Store != NULL) Destroy_CompCol_Matrix(&h->U);
    free(h->perm_c);
    free(h->perm_r);
    free(h->etree);
    free(h->colptr0);
    free(h->rowidx0);
    free(h->aval);
    free(h);
}

void *fsparse_slu_factor_z(long n, long nnz, const long *colptr1,
                           const long *rowidx1, const double *val, int *info)
{
    slu_cplx_t *h;
    SuperLUStat_t stat;
    int panel_size, relax;
    int_t lu_info = 0;

    *info = 0;
    h = (slu_cplx_t *) calloc(1, sizeof(slu_cplx_t));
    if (h == NULL) { *info = -1; return NULL; }

    h->n = (int) n;
    h->colptr0 = copy_to_int_t0(colptr1, n + 1);
    h->rowidx0 = copy_to_int_t0(rowidx1, nnz);
    h->aval = (doublecomplex *) malloc((size_t) nnz * sizeof(doublecomplex));
    h->perm_c = (int *) malloc((size_t) n * sizeof(int));
    h->perm_r = (int *) malloc((size_t) n * sizeof(int));
    h->etree = (int *) malloc((size_t) n * sizeof(int));
    if (h->colptr0 == NULL || h->rowidx0 == NULL || h->aval == NULL ||
        h->perm_c == NULL || h->perm_r == NULL || h->etree == NULL) {
        *info = -1;
        fsparse_slu_free_z(h);
        return NULL;
    }
    /* val is interleaved real(dp) pairs == doublecomplex {r, i}. */
    memcpy(h->aval, val, (size_t) nnz * sizeof(doublecomplex));

    zCreate_CompCol_Matrix(&h->A, (int) n, (int) n, (int_t) nnz, h->aval,
                           h->rowidx0, h->colptr0, SLU_NC, SLU_Z, SLU_GE);

    set_default_options(&h->options);
    h->options.ColPerm = COLAMD;
    get_perm_c(COLAMD, &h->A, h->perm_c);
    sp_preorder(&h->options, &h->A, h->perm_c, h->etree, &h->AC);

    panel_size = sp_ienv(1);
    relax = sp_ienv(2);
    StatInit(&stat);
    zgstrf(&h->options, &h->AC, relax, panel_size, h->etree, NULL, 0,
           h->perm_c, h->perm_r, &h->L, &h->U, &h->glu, &stat, &lu_info);
    StatFree(&stat);

    if (lu_info != 0) {
        *info = (lu_info > 0) ? 1 : -1;
        fsparse_slu_free_z(h);
        return NULL;
    }
    return h;
}

int fsparse_slu_solve_z(void *handle, const double *b, double *x, long n)
{
    slu_cplx_t *h = (slu_cplx_t *) handle;
    SuperMatrix B;
    SuperLUStat_t stat;
    int info = 0;

    if (h == NULL) return -1;
    memcpy(x, b, (size_t) n * sizeof(doublecomplex));
    zCreate_Dense_Matrix(&B, (int) n, 1, (doublecomplex *) x, (int) n, SLU_DN,
                         SLU_Z, SLU_GE);
    StatInit(&stat);
    zgstrs(NOTRANS, &h->L, &h->U, h->perm_c, h->perm_r, &B, &stat, &info);
    StatFree(&stat);
    Destroy_SuperMatrix_Store(&B);
    return info;
}

void fsparse_slu_free_z(void *handle)
{
    slu_cplx_t *h = (slu_cplx_t *) handle;
    if (h == NULL) return;
    if (h->A.Store != NULL) Destroy_SuperMatrix_Store(&h->A);
    if (h->AC.Store != NULL) Destroy_CompCol_Permuted(&h->AC);
    if (h->L.Store != NULL) Destroy_SuperNode_Matrix(&h->L);
    if (h->U.Store != NULL) Destroy_CompCol_Matrix(&h->U);
    free(h->perm_c);
    free(h->perm_r);
    free(h->etree);
    free(h->colptr0);
    free(h->rowidx0);
    free(h->aval);
    free(h);
}
