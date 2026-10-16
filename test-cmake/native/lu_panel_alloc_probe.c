/* SPDX-License-Identifier: MIT
 * Test-only requested-payload accounting and exhaustive allocation failures. */
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
void *__real_malloc(size_t);
void *__real_calloc(size_t, size_t);
void *__real_realloc(void *, size_t);
void __real_free(void *);
static struct { void *ptr; size_t bytes; } owned[256];
static int enabled, calls, fail_after, failed;
static int64_t current, peak;
static int slot(void *p)
{
    int i;
    for (i=0;i<256;++i) if (owned[i].ptr==p) return i;
    return -1;
}
static int reject(void)
{
    if (!enabled) return 0;
    ++calls;
    if (fail_after==0) { fail_after=-1; failed=1; return 1; }
    if (fail_after>0) --fail_after;
    return 0;
}
static void record(void *p,size_t bytes)
{
    int i;
    if (!enabled || !p) return;
    i=slot(NULL);
    if (i<0) abort();
    owned[i].ptr=p; owned[i].bytes=bytes;
    current+=(int64_t)bytes;
    if (current>peak) peak=current;
}
static void forget(void *p)
{
    int i;
    if (!p) return;
    i=slot(p);
    if (i<0) return;
    current-=(int64_t)owned[i].bytes;
    owned[i].ptr=NULL; owned[i].bytes=0;
}
void fs_panel_arm(int after)
{
    enabled=1; fail_after=after; calls=0; failed=0; peak=current;
}
int fs_panel_stop(void) { enabled=0; return calls; }
int fs_panel_failed(void) { return failed; }
int64_t fs_panel_current(void) { return current; }
int64_t fs_panel_peak(void) { return peak; }
void *__wrap_malloc(size_t n)
{
    void *p;
    if (reject()) return NULL;
    p=__real_malloc(n); record(p,n); return p;
}
void *__wrap_calloc(size_t n,size_t s)
{
    void *p;
    if (reject()) return NULL;
    p=__real_calloc(n,s); record(p,n*s); return p;
}
void *__wrap_realloc(void *old,size_t n)
{
    void *p;
    int i=old ? slot(old) : -1;
    if (reject()) return NULL;
    p=__real_realloc(old,n);
    if (p || n==0) {
        if (i>=0) {
            current-=(int64_t)owned[i].bytes;
            owned[i].ptr=NULL; owned[i].bytes=0;
        }
        record(p,n);
    }
    return p;
}
void __wrap_free(void *p) { forget(p); __real_free(p); }
