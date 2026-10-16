/* Independent static-link allocation and failure oracle; no production hooks. */
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include "fsparse_native_amd_bounded.h"
void *__real_malloc(size_t);
void *__real_calloc(size_t,size_t);
void *__real_realloc(void *,size_t);
void __real_free(void *);
static int enabled, amd_calls, opaque_calls;
static int64_t calls, fail_at, live, peak;
static struct { void *ptr; size_t bytes; } slots[8192];
static void remove_slot(void *ptr) {
    if (!enabled || !ptr) return;
    for (int i=0;i<8192;++i) if (slots[i].ptr==ptr) {
        live-=(int64_t)slots[i].bytes; slots[i].ptr=NULL; return;
    }
}
static void add_slot(void *ptr,size_t n) {
    if (!enabled || !ptr) return;
    /* Fortran's zero-length numeric arrays may request malloc(1): this is
       allocator minimum granularity, not nonzero array payload. */
    if (n==1) n=0;
    for (int i=0;i<8192;++i) if (!slots[i].ptr) {
        slots[i].ptr=ptr;slots[i].bytes=n;live+=(int64_t)n;
        if (live>peak) peak=live;
        return;
    }
    abort();
}
static int fail(void) { return enabled && ++calls==fail_at; }
void fsparse_spd_probe_begin(int64_t nth) {
    for (int i=0;i<8192;++i) slots[i].ptr=NULL;
    calls=live=peak=0;fail_at=nth;enabled=1;
}
void fsparse_spd_probe_end(int64_t *count,int64_t *current,int64_t *maximum) {
    enabled=0;*count=calls;*current=live;*maximum=peak;
}
void *__wrap_malloc(size_t n) {
    if (fail()) return NULL;
    void *p=__real_malloc(n);add_slot(p,n);return p;
}
void *__wrap_calloc(size_t n,size_t s) {
    if (fail()) return NULL;
    void *p=__real_calloc(n,s);add_slot(p,n*s);return p;
}
void *__wrap_realloc(void *p,size_t n) {
    if (fail()) return NULL;
    int slot=-1;
    if (enabled && p) for (int i=0;i<8192;++i)
        if (slots[i].ptr==p) { slot=i;break; }
    void *q=__real_realloc(p,n);
    if (q || n==0) {
        if (slot>=0) {live-=(int64_t)slots[slot].bytes;slots[slot].ptr=NULL;}
        add_slot(q,n);
    }
    return q;
}
void __wrap_free(void *p) { remove_slot(p);__real_free(p); }
int __real_fsparse_native_amd_bounded(int,int,const int*,const int*,int,int32_t*,int64_t,int*,int*);
int __wrap_fsparse_native_amd_bounded(int n,int nz,const int *ptr,const int *row,
    int low,int32_t *work,int64_t bytes,int *perm,int *inverse) {
    ++amd_calls;
    return __real_fsparse_native_amd_bounded(n,nz,ptr,row,low,work,bytes,perm,inverse);
}
int __real_fsparse_native_order_amd(int,int,const int*,const int*,int*,int*);
int __wrap_fsparse_native_order_amd(int n,int nz,const int *ptr,const int *row,int *p,int *ip) {
    ++opaque_calls;return __real_fsparse_native_order_amd(n,nz,ptr,row,p,ip);
}
int fsparse_spd_probe_amd_calls(void) { return amd_calls; }
int fsparse_spd_probe_opaque_calls(void) { return opaque_calls; }
