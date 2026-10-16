/* Executable-local, serial allocation and advisory oracle. No production state. */
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <errno.h>
#if defined(FS_TEST_WRAP_ADVICE)
static int count, deny_advice, active, failed_ordinal, ordinal;
static int64_t current_bytes, peak_bytes;
static struct { uintptr_t address; size_t bytes; } live[256];
void *__real_malloc(size_t);
void *__real_calloc(size_t,size_t);
void *__real_realloc(void*,size_t);
void __real_free(void*);
static void remember(void *p,size_t bytes) {
    if(!active || !p) return;
    for(size_t i=0;i<256;++i) if(!live[i].address) {
        live[i].address=(uintptr_t)p;live[i].bytes=bytes;current_bytes+=bytes;
        if(current_bytes>peak_bytes) peak_bytes=current_bytes;return;
    }
    abort();
}
static size_t forget(uintptr_t address) {
    if(!address) return 0;
    for(size_t i=0;i<256;++i) if(live[i].address==address) {
        size_t bytes=live[i].bytes;current_bytes-=bytes;
        live[i].address=0;return bytes;
    }
    return 0;
}
static int refuse(void) {
    if(!active) return 0;
    ++ordinal;
    if(ordinal!=failed_ordinal) return 0;
    errno=ENOMEM;return 1;
}
int fsparse_lu_alloc_test_begin(int nth) {
    ordinal=0;failed_ordinal=nth;peak_bytes=current_bytes;active=1;return 1;
}
int fsparse_lu_alloc_test_end(void) { active=0;return ordinal; }
int64_t fsparse_lu_alloc_test_current(void) { return current_bytes; }
int64_t fsparse_lu_alloc_test_peak(void) { return peak_bytes; }
void *__wrap_malloc(size_t bytes) {
    if(refuse()) return NULL;
    void *p=__real_malloc(bytes);remember(p,bytes);return p;
}
void *__wrap_calloc(size_t n,size_t width) {
    if(refuse()) return NULL;
    void *p=__real_calloc(n,width);remember(p,n*width);return p;
}
void *__wrap_realloc(void *old,size_t bytes) {
    if(refuse()) return NULL;
    const uintptr_t address=(uintptr_t)old;
    void *p=__real_realloc(old,bytes);
    if(p || !bytes) { forget(address);remember(p,bytes); }
    return p;
}
void __wrap_free(void *p) { forget((uintptr_t)p);__real_free(p); }
int __real_madvise(void *,size_t,int);
void fsparse_lu_discard_test_mode(int denied) { count=0;deny_advice=denied; }
int fsparse_lu_discard_test_calls(void) { return count; }
int __wrap_madvise(void *p,size_t bytes,int mode) {
    ++count;
    if(deny_advice) { errno=ENOMEM;return -1; }
    return __real_madvise(p,bytes,mode);
}
#else
int fsparse_lu_alloc_test_begin(int nth) { (void)nth;return 0; }
int fsparse_lu_alloc_test_end(void) { return 0; }
int64_t fsparse_lu_alloc_test_current(void) { return 0; }
int64_t fsparse_lu_alloc_test_peak(void) { return 0; }
void fsparse_lu_discard_test_mode(int denied) { (void)denied; }
int fsparse_lu_discard_test_calls(void) { return -1; }
#endif
