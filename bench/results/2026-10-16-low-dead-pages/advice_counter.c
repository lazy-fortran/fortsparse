/* Diagnostic ELF only. No values returned here enter ordinary timing results. */
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <errno.h>
#include <sys/mman.h>
static uint64_t helper_calls,eligible_payload_calls,advice_calls,advised_bytes,failures;
static uint64_t unrelated_calls,invalid_span;
static int depth;
static uintptr_t origin;
static size_t capacity;
static int errors[256];
void __real_fsparse_native_lu_discard_payload(void*,size_t);
int __real_madvise(void*,size_t,int);
void __wrap_fsparse_native_lu_discard_payload(void *payload,size_t bytes) {
    ++helper_calls;if(bytes>=1024u*1024u) ++eligible_payload_calls;
    ++depth;origin=(uintptr_t)payload;capacity=bytes;
    __real_fsparse_native_lu_discard_payload(payload,bytes);--depth;
}
int __wrap_madvise(void *p,size_t bytes,int mode) {
    const int result=__real_madvise(p,bytes,mode), error=errno;
    if(depth) {
        ++advice_calls;advised_bytes+=bytes;
        const uintptr_t first=(uintptr_t)p;
        if(mode!=MADV_DONTNEED || first<=origin || first>UINTPTR_MAX-bytes ||
           origin>UINTPTR_MAX-capacity || first+bytes>=origin+capacity) ++invalid_span;
        if(result) { ++failures;if(error>=0 && error<256) ++errors[error]; }
    } else ++unrelated_calls;
    errno=error;return result;
}
__attribute__((destructor)) static void report(void) {
    fprintf(stderr,"ADVICE_DIAGNOSTIC {\"helper_calls\":%llu,\"eligible_payload_calls\":%llu,"
        "\"advice_calls\":%llu,\"cumulative_advised_bytes\":%llu,\"failures\":%llu,"
        "\"unrelated_calls\":%llu,\"invalid_interior_spans\":%llu,\"errno_counts\":{",
        (unsigned long long)helper_calls,(unsigned long long)eligible_payload_calls,
        (unsigned long long)advice_calls,(unsigned long long)advised_bytes,
        (unsigned long long)failures,(unsigned long long)unrelated_calls,
        (unsigned long long)invalid_span);
    int first=1;
    for(int i=0;i<256;++i) if(errors[i]) {
        fprintf(stderr,"%s\"%d\":%d",first?"":",",i,errors[i]);first=0;
    }
    fprintf(stderr,"}}\n");
}
