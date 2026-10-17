/* Thin C shim over sequential SuperLU (BSD). Exposes a clean
 * factor/solve/free C ABI that Fortran binds via iso_c_binding. The
 * factorization (L, U, perm_c, perm_r, the column-permuted SuperMatrix and the
 * original CSC SuperMatrix) is held in an opaque handle so repeated solves
 * reuse it. License of this file: MIT (our code). It links only SuperLU. */

#ifdef FORTSPARSE_HAVE_SUPERLU

#include <stdlib.h>
#include <string.h>
#include <limits.h>

#if defined(__has_include)
#if __has_include(<superlu/slu_ddefs.h>)
#include <superlu/slu_ddefs.h>
#include <superlu/slu_zdefs.h>
#else
#include <slu_ddefs.h>
#include <slu_zdefs.h>
#endif
#else
#include <superlu/slu_ddefs.h>
#include <superlu/slu_zdefs.h>
#endif

void *fsparse_slu_factor_d(int n, int nnz, const int *colptr1,
                           const int *rowidx1, const double *val, int symmetric, int ilu,
                           int *info);
int fsparse_slu_refactor_d(void *handle, int n, int nnz, const int *colptr1,
                           const int *rowidx1, const double *val, int symmetric);
int fsparse_slu_solve_d(void *handle, const double *b, double *x, long n);
void fsparse_slu_free_d(void *handle);
void *fsparse_slu_factor_z(long n, long nnz, const long *colptr1,
                           const long *rowidx1, const double *val, int symmetric,
                           int *info);
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
    double *rscale, *cscale;
    char equed;
    int incomplete;
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

/* The public real CSC uses default Fortran integers. Convert directly into
 * SuperLU storage without a matrix-sized intermediate long array. */
static int_t *copy_int_to_int_t0(const int *src, int count)
{
    int_t *dst = (int_t *) malloc((size_t) count * sizeof(int_t));
    int k;
    if (dst == NULL) return NULL;
    for (k = 0; k < count; ++k) dst[k] = (int_t) (src[k] - 1);
    return dst;
}

/* Shim factor codes: 0 ok, 1 singular, 2 allocation failure, 3 index
 * overflow (size exceeds SuperLU's int_t), -1 internal.
 * SuperLU gstrf: 0 < info <= n means U(info,info) is exactly zero; info > n
 * means a memory allocation failed (info = bytes allocated so far + n);
 * info < 0 is an invalid argument. */
static int classify_lu_info(long lu_info, long n)
{
    if (lu_info == 0) return 0;
    if (lu_info < 0) return -1;
    return (lu_info <= n) ? 1 : 2;
}

/* Unsymmetric input: COLAMD on A^T A with partial pivoting (SuperLU
 * default). Symmetric input (symmetric=1, the caller's SPD/INDEFINITE/
 * HERMITIAN hint): the SuperLU User's Guide setting for symmetric or
 * structurally symmetric matrices -- minimum degree on A^T + A, symmetric
 * mode, diagonal pivots preferred above the threshold 0.001 (threshold
 * partial pivoting otherwise). The full matrix is retained, including
 * off-diagonal pivots needed by indefinite systems. */
static void set_options(superlu_options_t *options, int symmetric)
{
    set_default_options(options);
    if (symmetric) {
        options->ColPerm = MMD_AT_PLUS_A;
        options->SymmetricMode = YES;
        options->DiagPivotThresh = 0.001;
    } else {
        options->ColPerm = COLAMD;
    }
}

void *fsparse_slu_factor_d(int n, int nnz, const int *colptr1,
                           const int *rowidx1, const double *val, int symmetric, int ilu,
                           int *info)
{
    slu_real_t *h;
    SuperMatrix B; /* dummy, unused by dgstrf */
    SuperLUStat_t stat;
    int panel_size, relax;
    int_t lu_info = 0;
    double rowcnd, colcnd, amax;
    int equ_info;

    *info = 0;
    h = (slu_real_t *) calloc(1, sizeof(slu_real_t));
    if (h == NULL) { *info = 2; return NULL; }

    h->n = (int) n;
    h->incomplete = ilu;
    h->equed = 'N';
    h->colptr0 = copy_int_to_int_t0(colptr1, n + 1);
    h->rowidx0 = copy_int_to_int_t0(rowidx1, nnz);
    h->aval = (double *) malloc((size_t) nnz * sizeof(double));
    h->perm_c = (int *) malloc((size_t) n * sizeof(int));
    h->perm_r = (int *) malloc((size_t) n * sizeof(int));
    h->etree = (int *) malloc((size_t) n * sizeof(int));
    if (h->colptr0 == NULL || h->rowidx0 == NULL || h->aval == NULL ||
        h->perm_c == NULL || h->perm_r == NULL || h->etree == NULL) {
        *info = 2;
        fsparse_slu_free_d(h);
        return NULL;
    }
    memcpy(h->aval, val, (size_t) nnz * sizeof(double));

    dCreate_CompCol_Matrix(&h->A, (int) n, (int) n, (int_t) nnz, h->aval,
                           h->rowidx0, h->colptr0, SLU_NC, SLU_D, SLU_GE);

    set_options(&h->options, symmetric);
    if (ilu) {
        ilu_set_default_options(&h->options);
        /* Equilibration without the optional MC64 dependency. ILU uses the
         * general-matrix policy even if a symmetry hint was supplied. */
        h->options.RowPerm = NOROWPERM;
        h->options.ColPerm = COLAMD;
        h->options.ILU_DropTol = 1e-4;
        h->options.ILU_FillFactor = 10.0;
        h->rscale = (double *) malloc((size_t) n * sizeof(double));
        h->cscale = (double *) malloc((size_t) n * sizeof(double));
        if (h->rscale == NULL || h->cscale == NULL) {
            *info = 2;
            fsparse_slu_free_d(h);
            return NULL;
        }
        dgsequ(&h->A, h->rscale, h->cscale, &rowcnd, &colcnd, &amax, &equ_info);
        if (equ_info != 0) {
            *info = 1;
            fsparse_slu_free_d(h);
            return NULL;
        }
        dlaqgs(&h->A, h->rscale, h->cscale, rowcnd, colcnd, amax, &h->equed);
    }
    get_perm_c(h->options.ColPerm, &h->A, h->perm_c);
    sp_preorder(&h->options, &h->A, h->perm_c, h->etree, &h->AC);

    panel_size = sp_ienv(1);
    relax = sp_ienv(2);
    StatInit(&stat);
    memset(&B, 0, sizeof(B));
    if (ilu)
        dgsitrf(&h->options, &h->AC, relax, panel_size, h->etree, NULL, 0,
                h->perm_c, h->perm_r, &h->L, &h->U, &h->glu, &stat, &lu_info);
    else
        dgstrf(&h->options, &h->AC, relax, panel_size, h->etree, NULL, 0,
                h->perm_c, h->perm_r, &h->L, &h->U, &h->glu, &stat, &lu_info);
    StatFree(&stat);

    /* dgsitrf may perturb zero pivots and still produce a usable
     * approximate inverse. A non-NULL ILU handle disambiguates this count
     * from the exact factor/error status codes. */
    if (ilu && lu_info >= 0 && lu_info <= n) {
        *info = (int) lu_info;
        return h;
    }
    if (lu_info != 0) {
        *info = classify_lu_info((long) lu_info, (long) n);
        fsparse_slu_free_d(h);
        return NULL;
    }
    return h;
}

/* Numerical refactorization of a retained real handle with new values on the
 * same CSC pattern and symmetry hint. The column ordering perm_c depends on
 * the pattern only, so it is kept and the minimum-degree ordering (the
 * dominant cost of a small factorization) is skipped; the column elimination
 * tree, row pivoting and L/U are recomputed exactly as in a fresh factor, so
 * the result is identical to fsparse_slu_factor_d on the same input.
 * Returns 4 (handle unchanged) when the pattern or hint differs, otherwise the
 * fsparse_slu_factor_d info codes; on a nonzero code other than 4 the handle
 * holds no factorization and must be freed. */
int fsparse_slu_refactor_d(void *handle, int n, int nnz, const int *colptr1,
                           const int *rowidx1, const double *val, int symmetric)
{
    slu_real_t *h = (slu_real_t *) handle;
    SuperLUStat_t stat;
    int panel_size, relax;
    int_t lu_info = 0;
    int k;

    if (h == NULL || h->n != n || h->incomplete) return 4;
    if ((h->options.SymmetricMode == YES) != (symmetric != 0)) return 4;
    if (((NCformat *) h->A.Store)->nnz != (int_t) nnz) return 4;
    for (k = 0; k <= n; ++k)
        if (h->colptr0[k] != (int_t) (colptr1[k] - 1)) return 4;
    for (k = 0; k < nnz; ++k)
        if (h->rowidx0[k] != (int_t) (rowidx1[k] - 1)) return 4;

    if (h->AC.Store != NULL) Destroy_CompCol_Permuted(&h->AC);
    if (h->L.Store != NULL) Destroy_SuperNode_Matrix(&h->L);
    if (h->U.Store != NULL) Destroy_CompCol_Matrix(&h->U);
    h->AC.Store = NULL;
    h->L.Store = NULL;
    h->U.Store = NULL;
    memcpy(h->aval, val, (size_t) nnz * sizeof(double));
    sp_preorder(&h->options, &h->A, h->perm_c, h->etree, &h->AC);

    panel_size = sp_ienv(1);
    relax = sp_ienv(2);
    StatInit(&stat);
    dgstrf(&h->options, &h->AC, relax, panel_size, h->etree, NULL, 0,
           h->perm_c, h->perm_r, &h->L, &h->U, &h->glu, &stat, &lu_info);
    StatFree(&stat);
    if (lu_info != 0) return classify_lu_info((long) lu_info, (long) n);
    return 0;
}

int fsparse_slu_solve_op_d(void *handle, const double *b, double *x, long n,
                            int operation)
{
    slu_real_t *h = (slu_real_t *) handle;
    SuperMatrix B;
    SuperLUStat_t stat;
    int info = 0;
    long k;

    if (h == NULL || n != h->n || operation < 0 || operation > 2) return -1;
    memcpy(x, b, (size_t) n * sizeof(double));
    /* E = R A C. Normal: C E^-1 R; transpose/adjoint: R E^-T C. */
    if (operation == 0) {
        if (h->equed == 'R' || h->equed == 'B')
            for (k = 0; k < n; ++k) x[k] *= h->rscale[k];
    } else {
        if (h->equed == 'C' || h->equed == 'B')
            for (k = 0; k < n; ++k) x[k] *= h->cscale[k];
    }
    dCreate_Dense_Matrix(&B, (int) n, 1, x, (int) n, SLU_DN, SLU_D, SLU_GE);
    StatInit(&stat);
    dgstrs(operation == 0 ? NOTRANS : TRANS, &h->L, &h->U,
           h->perm_c, h->perm_r, &B, &stat, &info);
    StatFree(&stat);
    Destroy_SuperMatrix_Store(&B);
    if (operation == 0) {
        if (h->equed == 'C' || h->equed == 'B')
            for (k = 0; k < n; ++k) x[k] *= h->cscale[k];
    } else {
        if (h->equed == 'R' || h->equed == 'B')
            for (k = 0; k < n; ++k) x[k] *= h->rscale[k];
    }
    return info;
}

int fsparse_slu_solve_d(void *handle, const double *b, double *x, long n)
{
    return fsparse_slu_solve_op_d(handle, b, x, n, 0);
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
    free(h->rscale);
    free(h->cscale);
    free(h);
}

void *fsparse_slu_factor_z(long n, long nnz, const long *colptr1,
                           const long *rowidx1, const double *val, int symmetric,
                           int *info)
{
    slu_cplx_t *h;
    SuperLUStat_t stat;
    int panel_size, relax;
    int_t lu_info = 0;

    *info = 0;
    /* SuperLU indices are int_t; refuse sizes that would silently truncate. */
    if (sizeof(int_t) < sizeof(long) &&
        (n > (long) INT_MAX || nnz > (long) INT_MAX)) {
        *info = 3;
        return NULL;
    }
    h = (slu_cplx_t *) calloc(1, sizeof(slu_cplx_t));
    if (h == NULL) { *info = 2; return NULL; }

    h->n = (int) n;
    h->colptr0 = copy_to_int_t0(colptr1, n + 1);
    h->rowidx0 = copy_to_int_t0(rowidx1, nnz);
    h->aval = (doublecomplex *) malloc((size_t) nnz * sizeof(doublecomplex));
    h->perm_c = (int *) malloc((size_t) n * sizeof(int));
    h->perm_r = (int *) malloc((size_t) n * sizeof(int));
    h->etree = (int *) malloc((size_t) n * sizeof(int));
    if (h->colptr0 == NULL || h->rowidx0 == NULL || h->aval == NULL ||
        h->perm_c == NULL || h->perm_r == NULL || h->etree == NULL) {
        *info = 2;
        fsparse_slu_free_z(h);
        return NULL;
    }
    /* val is interleaved real(dp) pairs == doublecomplex {r, i}. */
    memcpy(h->aval, val, (size_t) nnz * sizeof(doublecomplex));

    zCreate_CompCol_Matrix(&h->A, (int) n, (int) n, (int_t) nnz, h->aval,
                           h->rowidx0, h->colptr0, SLU_NC, SLU_Z, SLU_GE);

    set_options(&h->options, symmetric);
    get_perm_c(h->options.ColPerm, &h->A, h->perm_c);
    sp_preorder(&h->options, &h->A, h->perm_c, h->etree, &h->AC);

    panel_size = sp_ienv(1);
    relax = sp_ienv(2);
    StatInit(&stat);
    zgstrf(&h->options, &h->AC, relax, panel_size, h->etree, NULL, 0,
           h->perm_c, h->perm_r, &h->L, &h->U, &h->glu, &stat, &lu_info);
    StatFree(&stat);

    if (lu_info != 0) {
        *info = classify_lu_info((long) lu_info, (long) n);
        fsparse_slu_free_z(h);
        return NULL;
    }
    return h;
}

int fsparse_slu_solve_op_z(void *handle, const double *b, double *x, long n,
                            int operation)
{
    slu_cplx_t *h = (slu_cplx_t *) handle;
    SuperMatrix B;
    SuperLUStat_t stat;
    int info = 0;

    if (h == NULL || n != h->n || operation < 0 || operation > 2) return -1;
    memcpy(x, b, (size_t) n * sizeof(doublecomplex));
    zCreate_Dense_Matrix(&B, (int) n, 1, (doublecomplex *) x, (int) n, SLU_DN,
                         SLU_Z, SLU_GE);
    StatInit(&stat);
    zgstrs(operation == 0 ? NOTRANS : (operation == 1 ? TRANS : CONJ),
           &h->L, &h->U, h->perm_c, h->perm_r, &B, &stat, &info);
    StatFree(&stat);
    Destroy_SuperMatrix_Store(&B);
    return info;
}

int fsparse_slu_solve_z(void *handle, const double *b, double *x, long n)
{
    return fsparse_slu_solve_op_z(handle, b, x, n, 0);
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

#endif /* FORTSPARSE_HAVE_SUPERLU */
