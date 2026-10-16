/* SPDX-License-Identifier: MIT
 * Test-only allocator fault injection and allocation-free apply observation. */
#include <stddef.h>
void *__real_malloc(size_t);
void *__real_calloc(size_t, size_t);
void *__real_realloc(void *, size_t);
static int enabled, remaining, calls;
void fs_ilutp_alloc_arm(int after) { enabled=1; remaining=after; calls=0; }
int fs_ilutp_alloc_stop(void) { enabled=0; return calls; }
static int reject(void)
{
    if (!enabled) return 0;
    ++calls;
    if (remaining == 0) { remaining=-1; return 1; }
    if (remaining > 0) --remaining;
    return 0;
}
void *__wrap_malloc(size_t n) { return reject() ? NULL : __real_malloc(n); }
void *__wrap_calloc(size_t n, size_t size)
{ return reject() ? NULL : __real_calloc(n, size); }
void *__wrap_realloc(void *p, size_t n)
{ return reject() ? NULL : __real_realloc(p, n); }
