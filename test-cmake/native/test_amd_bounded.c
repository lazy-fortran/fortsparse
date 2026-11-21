#include "fsparse_native_amd_bounded.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
static int expected_n, expected_m, expected_low, graph[16][16];
static volatile int tracking;
static volatile unsigned long long allocations;
/* Test-executable symbol interposition counts every heap allocation. */
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
int main(void) {
    int64_t bytes=17,outputs=19;
    int ptr[17],row[256],perm[16],inverse[16];
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
        int nz=0;
        for(int j=0;j<expected_n;++j) {
            ptr[j]=nz+1;
            for(int i=0;i<expected_n;++i) if(graph[i][j]) {
                row[nz++]=i+1;
            }
        }
        ptr[expected_n]=nz+1;
        expected_m=0;
        for(int i=0;i<expected_n;++i) for(int j=0;j<expected_n;++j)
            if(i!=j && (graph[i][j] || graph[j][i])) ++expected_m;
        for(expected_low=0;expected_low<=1;++expected_low) {
            require(fsparse_native_amd_bounded_size(expected_n,nz,ptr,row,expected_low,&bytes,&outputs)==0,"size query");
            int64_t iwlen=expected_m+expected_n+(expected_low ? 0 : expected_m/5);
            require(bytes==4*(7*(int64_t)expected_n+iwlen) && outputs==8*expected_n,"exact arena and output bytes");
            int32_t *work=malloc((size_t)bytes);
            require(work!=NULL,"test arena allocation");
            require(fsparse_native_amd_bounded(expected_n,nz,ptr,row,expected_low,work,bytes-1,perm,inverse)==2,"one-byte-short arena");
            tracking=1;
            int status=fsparse_native_amd_bounded(expected_n,nz,ptr,row,expected_low,work,bytes,perm,inverse);
            tracking=0;
            require(status==0,"bounded AMD execute");
            require(allocations==0,"bounded AMD allocated heap storage");
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
    printf("bounded AMD: exact buffers and permutations, heap allocations %llu\n",allocations);
    return 0;
}
