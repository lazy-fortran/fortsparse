/* MIT: owned graph construction and bounded AMD workspace.
 * The public amd_2 interface is unchecked and destructive. This adapter creates
 * its complete canonical symmetric off-diagonal graph before the vendor call.
 * Audit: AMD 3.3.4 amd_2 -> amd_postorder -> iterative amd_post_tree; no heap
 * calls. See docs/native-amd-bounded.md. No global SuiteSparse hooks are used.
 * Without the audited vendor ABI the project-owned core takes the same arena. */
#include "fsparse_native_amd_bounded.h"
#include <limits.h>
#include <stddef.h>
#include "fsparse_order_graph.h"
_Static_assert(sizeof(int) == sizeof(int32_t), "AMD adapter needs 32-bit C int");
#ifdef FORTSPARSE_HAVE_BOUNDED_AMD
#include <amd.h>
#if AMD_MAIN_VERSION != 3 || AMD_SUB_VERSION != 3 || AMD_SUBSUB_VERSION != 4
#error "Bounded AMD requires the audited AMD 3.3.4 header"
#endif
#else
#include "fsparse_native_amd_core.h"
#endif

enum { FS_OK=0, FS_INVALID=2, FS_UNAVAILABLE=3, FS_INTERNAL=5, FS_INDEX=7 };
int fsparse_native_amd_bounded_available(void) {
#ifdef FORTSPARSE_HAVE_BOUNDED_AMD
    int version[3];
    amd_version(version);
    return version[0]==3 && version[1]==3 && version[2]==4;
#else
    return 1;
#endif
}

struct layout { int32_t k, m, iwlen; int64_t words; };
static int plan(int n, int nnz, const int *ptr, const int *row, int low,
                struct layout *out) {
    if (!fsparse_native_amd_bounded_available()) return FS_UNAVAILABLE;
    if (n<=0 || nnz<0 || (low!=0 && low!=1)) return FS_INVALID;
    /* Keep AMD's documented public front-end size restriction; additionally
     * check the symmetrized graph and every destructive workspace offset. */
    if (n>=INT32_MAX/4 || nnz>=INT32_MAX/4) return FS_INDEX;
    if (!fsparse_csc_graph_valid(n,nnz,ptr,row)) return FS_INVALID;
    int64_t k=0, pairs=0;
    for (int j=0;j<n;++j) {
        for (int p=ptr[j]-1;p<ptr[j+1]-1;++p) {
            int i=row[p]-1;
            if (i==j) continue;
            ++k;
            if (i<j || !fsparse_graph_has_row(ptr,row,i,j+1)) ++pairs;
        }
    }
    /* n,nnz guards bound these i8 products well below INT64_MAX. */
    int64_t m=2*pairs;
    int64_t iwlen=m+(low ? 0 : m/5)+n;
    int64_t words=k+8*(int64_t)n+1+iwlen;
    if (k>=INT32_MAX || m>=INT32_MAX || iwlen>=INT32_MAX ||
        words>INT64_MAX/4 || (uint64_t)words>SIZE_MAX/sizeof(int32_t)) return FS_INDEX;
    out->k=(int32_t)k; out->m=(int32_t)m;
    out->iwlen=(int32_t)iwlen; out->words=words;
    return FS_OK;
}

int fsparse_native_amd_bounded_size(int n, int nnz, const int *ptr, const int *row,
    int low, int64_t *workspace_bytes, int64_t *output_bytes) {
    if (!workspace_bytes || !output_bytes) return FS_INVALID;
    struct layout s;
    int code=plan(n,nnz,ptr,row,low,&s);
    if (code!=FS_OK) return code;
    *workspace_bytes=4*s.words;
    *output_bytes=8*(int64_t)n;
    return FS_OK;
}

int fsparse_native_amd_bounded(int n, int nnz, const int *ptr, const int *row,
    int low, int32_t *workspace, int64_t workspace_bytes, int *perm, int *inverse) {
    struct layout s;
    int code=plan(n,nnz,ptr,row,low,&s);
    if (code!=FS_OK) return code;
    if (!workspace || !perm || !inverse || workspace_bytes<4*s.words) return FS_INVALID;
    int32_t *tp=workspace, *ti=tp+n+1, *pe=ti+s.k, *len=pe+n,
        *nv=len+n, *head=nv+n, *elen=head+n, *degree=elen+n,
        *w=degree+n, *iw=w+n;
    for (int j=0;j<n;++j) len[j]=0;
    for (int j=0;j<n;++j) for (int p=ptr[j]-1;p<ptr[j+1]-1;++p)
        if (row[p]-1!=j) ++len[row[p]-1];
    tp[0]=0;
    for (int j=0;j<n;++j) { tp[j+1]=tp[j]+len[j]; nv[j]=tp[j]; }
    if (tp[n]!=s.k) return FS_INTERNAL;
    /* Ascending source columns give sorted transpose rows without sorting. */
    for (int j=0;j<n;++j) for (int p=ptr[j]-1;p<ptr[j+1]-1;++p) {
        int i=row[p]-1;
        if (i!=j) ti[nv[i]++]=j;
    }
    int used=0;
    for (int j=0;j<n;++j) {
        int a=ptr[j]-1, end=ptr[j+1]-1, t=tp[j], tend=tp[j+1];
        pe[j]=used;
        int previous=-1;
        while (a<end || t<tend) {
            int left=a<end ? row[a]-1 : INT_MAX;
            int right=t<tend ? ti[t] : INT_MAX;
            int i=left<right ? left : right;
            if (left==i) ++a;
            if (right==i) ++t;
            if (i==j) continue;
            if (i<0 || i>=n || i<=previous || used>=s.m) return FS_INTERNAL;
            iw[used++]=i; previous=i;
        }
        len[j]=used-pe[j];
    }
    if (used!=s.m) return FS_INTERNAL;
#ifdef FORTSPARSE_HAVE_BOUNDED_AMD
    amd_2((int32_t)n,pe,iw,len,s.iwlen,s.m,nv,(int32_t *)inverse,
          (int32_t *)perm,head,elen,degree,w,NULL,NULL);
#else
    if (fsparse_amd_core((int32_t)n,pe,iw,len,s.iwlen,s.m,nv,(int32_t *)inverse,
                         (int32_t *)perm,head,elen,degree,w)!=0) return FS_INTERNAL;
#endif
    for (int i=0;i<n;++i)
        if (perm[i]<0 || perm[i]>=n || inverse[i]<0 || inverse[i]>=n ||
            perm[inverse[i]]!=i) return FS_INTERNAL;
    for (int i=0;i<n;++i) { ++perm[i]; ++inverse[i]; }
    return FS_OK;
}
