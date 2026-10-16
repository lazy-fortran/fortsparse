#include "fsparse_native_amd_bounded.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#ifdef __GLIBC__
#include <dlfcn.h>
#endif
#ifdef FS_TEST_BOUNDED_AMD
#include <amd.h>
#endif
static int calls, expected_n, expected_m, expected_low, graph[16][16];
static int mismatch_version;
static volatile int tracking;
static volatile unsigned long long allocations;
/* Test-executable symbol interposition reaches AMD's DSO too. It neither
 * changes SuiteSparse hooks nor preloads instrumentation into production. */
#ifdef __GLIBC__
extern void *__libc_malloc(size_t);
extern void *__libc_calloc(size_t,size_t);
extern void *__libc_realloc(void *,size_t);
void *malloc(size_t n) { if(tracking) ++allocations; return __libc_malloc(n); }
void *calloc(size_t n,size_t s) { if(tracking) ++allocations; return __libc_calloc(n,s); }
void *realloc(void *p,size_t n) { if(tracking) ++allocations; return __libc_realloc(p,n); }
#endif
static void require(int condition,const char *message) {
    if (!condition) { fprintf(stderr,"%s\n",message); exit(1); }
}
#ifdef FS_TEST_BOUNDED_AMD
void __real_amd_version(int *);
void __wrap_amd_version(int *version) {
    __real_amd_version(version);
    if(mismatch_version) ++version[2];
}
void __real_amd_2(int32_t,int32_t*,int32_t*,int32_t*,int32_t,int32_t,
    int32_t*,int32_t*,int32_t*,int32_t*,int32_t*,int32_t*,int32_t*,double*,double*);
void __wrap_amd_2(int32_t n,int32_t *pe,int32_t *iw,int32_t *len,int32_t iwlen,int32_t pfree,
    int32_t *nv,int32_t *next,int32_t *last,int32_t *head,int32_t *elen,
    int32_t *degree,int32_t *w,double *control,double *info) {
    ++calls;
    require(n==expected_n && pfree==expected_m,"independent symmetric graph size");
    require(iwlen==pfree+n+(expected_low ? 0 : pfree/5),"LOW/default elbow workspace");
    require(!control && !info,"bounded default control without uncharged arrays");
    int used=0;
    for(int j=0;j<n;++j) {
        require(pe[j]==used,"canonical symmetric row pointer");
        int count=0;
        for(int i=0;i<n;++i) if(i!=j && (graph[i][j] || graph[j][i])) {
            require(count<len[j] && iw[pe[j]+count]==i,"independent original plus transpose graph oracle");
            ++count;
        }
        require(count==len[j],"off-diagonal row count");
        used+=count;
    }
    __real_amd_2(n,pe,iw,len,iwlen,pfree,nv,next,last,head,elen,degree,w,control,info);
}
#endif
int main(void) {
    int64_t bytes=17,outputs=19;
    int ptr[17],row[256],perm[16],inverse[16];
    if(!fsparse_native_amd_bounded_available()) {
        require(fsparse_native_amd_bounded_size(1,0,ptr,row,0,&bytes,&outputs)==3,"missing audited AMD size gate");
        require(fsparse_native_amd_bounded(1,0,ptr,row,0,NULL,0,perm,inverse)==3,"missing audited AMD execute gate");
        return 0;
    }
#if defined(__GLIBC__) && defined(FS_TEST_BOUNDED_AMD)
    void *(*probe_allocate)(size_t,size_t)=(void *(*)(size_t,size_t))dlsym(RTLD_DEFAULT,"SuiteSparse_malloc");
    void *(*probe_free)(void *)=(void *(*)(void *))dlsym(RTLD_DEFAULT,"SuiteSparse_free");
    require(probe_allocate && probe_free,"loaded SuiteSparse DSO allocator boundary");
    tracking=1;
    void *probe=probe_allocate(1,1);
    tracking=0;
    require(probe!=NULL && allocations>0,"DSO allocator interception self-check");
    probe_free(probe);
    allocations=0;
#endif
    int empty[1]={1};
    require(fsparse_native_amd_bounded_size(0,0,empty,NULL,0,&bytes,&outputs)==2,"empty dimension rejected");
    require(fsparse_native_amd_bounded_size(INT_MAX,0,empty,NULL,0,&bytes,&outputs)==7,"large dimension before pointer access");
    require(bytes==17 && outputs==19,"query failure leaves outputs unchanged");
    for(int kind=0;kind<6;++kind) {
        expected_n=kind==0 ? 1 : 9;
        memset(graph,0,sizeof(graph));
        for(int i=0;i<expected_n;++i) graph[i][i]=kind!=1;
        for(int i=1;i<expected_n;++i) {
            if(kind==2) { graph[0][i]=1; graph[i][0]=1; }
            if(kind==3) graph[i][0]=1;
            if(kind==4) { graph[0][i]=1; if(i%2) graph[i][0]=1; }
            if(kind==5 && i<expected_n-2) { graph[i][i-1]=1; graph[i-1][i]=1; }
        }
        int nz=0,k=0;
        for(int j=0;j<expected_n;++j) {
            ptr[j]=nz+1;
            for(int i=0;i<expected_n;++i) if(graph[i][j]) {
                row[nz++]=i+1;
                if(i!=j) ++k;
            }
        }
        ptr[expected_n]=nz+1;
        expected_m=0;
        for(int i=0;i<expected_n;++i) for(int j=0;j<expected_n;++j)
            if(i!=j && (graph[i][j] || graph[j][i])) ++expected_m;
        for(expected_low=0;expected_low<=1;++expected_low) {
            require(fsparse_native_amd_bounded_size(expected_n,nz,ptr,row,expected_low,&bytes,&outputs)==0,"size query");
            int64_t iwlen=expected_m+expected_n+(expected_low ? 0 : expected_m/5);
            require(bytes==4*(k+8*(int64_t)expected_n+1+iwlen) && outputs==8*expected_n,"exact arena and output bytes");
            int32_t *work=malloc((size_t)bytes);
            require(work!=NULL,"test arena allocation");
            require(fsparse_native_amd_bounded(expected_n,nz,ptr,row,expected_low,work,bytes-1,perm,inverse)==2,"one-byte-short arena");
            tracking=1;
            int status=fsparse_native_amd_bounded(expected_n,nz,ptr,row,expected_low,work,bytes,perm,inverse);
            tracking=0;
            require(status==0,"bounded AMD execute");
            require(allocations==0,"DSO-aware AMD invocation allocated heap storage");
            for(int i=0;i<expected_n;++i) require(perm[inverse[i]-1]==i+1,"one-based permutation inverse");
            free(work);
        }
    }
    ptr[0]=0;
    require(fsparse_native_amd_bounded_size(9,0,ptr,row,0,&bytes,&outputs)==2,"invalid CSC pointers");
    int dup_ptr[3]={1,3,3},dup_row[2]={2,2};
    require(fsparse_native_amd_bounded_size(2,2,dup_ptr,dup_row,0,&bytes,&outputs)==2,"duplicate CSC rows");
    dup_row[1]=3;
    require(fsparse_native_amd_bounded_size(2,2,dup_ptr,dup_row,0,&bytes,&outputs)==2,"out of range row");
#ifdef FS_TEST_BOUNDED_AMD
    require(calls==12,"all vendor invocations observed");
    mismatch_version=1;
    require(!fsparse_native_amd_bounded_available(),"unreviewed runtime version rejected");
    require(fsparse_native_amd_bounded_size(1,0,empty,NULL,0,&bytes,&outputs)==3,"version mismatch query gate");
    require(fsparse_native_amd_bounded(1,0,empty,NULL,0,NULL,0,perm,inverse)==3,"version mismatch execute gate");
#endif
#ifdef __GLIBC__
    printf("bounded AMD: canonical graph oracle, exact buffers, %d vendor calls, DSO-aware heap allocations %llu\n",calls,allocations);
#else
    printf("bounded AMD: exact buffers and permutations; DSO-aware heap oracle unavailable\n");
#endif
    return 0;
}
