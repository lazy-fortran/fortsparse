/* Independent payload/guard oracle, plus injected advisory failure. */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <errno.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#if defined(__linux__)
#include <unistd.h>
#include <sys/mman.h>
#endif
void fsparse_native_lu_discard_payload(void *, size_t);
#if defined(FS_TEST_WRAP_ADVICE)
static int deny, calls, last_result;
static void *advised;
static size_t advised_bytes;
int __real_madvise(void *, size_t, int);
int __wrap_madvise(void *p, size_t bytes, int mode)
{
    ++calls; advised = p; advised_bytes = bytes;
    assert(mode == MADV_DONTNEED);
    if (deny) { errno = ENOMEM; last_result = -1; return -1; }
    last_result = __real_madvise(p, bytes, mode);
    return last_result;
}
#endif
int main(void)
{
    const size_t bytes = (size_t)2 * 1024 * 1024 + 137;
    unsigned char *raw = malloc(bytes + 16384);
    assert(raw != NULL);
    for (int offset = 0; offset < 4; ++offset) {
        memset(raw, 0x5a, bytes + 16384);
        unsigned char *payload = raw + 37 + offset;
        errno = EDOM;
        fsparse_native_lu_discard_payload(payload, bytes);
        assert(errno == EDOM);
#if defined(FS_TEST_WRAP_ADVICE) && !defined(FS_TEST_DISCARD_DISABLED)
        assert(calls == offset + 1 && last_result == 0);
        const uintptr_t lo = (uintptr_t)advised;
        const uintptr_t hi = lo + advised_bytes;
        const size_t page = (size_t)sysconf(_SC_PAGESIZE);
        assert(lo > (uintptr_t)payload && hi < (uintptr_t)payload + bytes);
        assert(lo % page == 0 && hi % page == 0 && hi > lo);
        for (size_t i = 0; i < bytes + 16384; ++i) {
            const uintptr_t address = (uintptr_t)(raw + i);
            assert(raw[i] == (address >= lo && address < hi ? 0 : 0x5a));
        }
#else
        /* Portable/no-advice compilation keeps every payload/guard byte. */
        for (size_t i = 0; i < bytes + 16384; ++i) assert(raw[i] == 0x5a);
#endif
    }
#if defined(FS_TEST_WRAP_ADVICE)
    const int before = calls;
    fsparse_native_lu_discard_payload(NULL, bytes);
    fsparse_native_lu_discard_payload(raw, 0);
    fsparse_native_lu_discard_payload(raw, 16384);
    fsparse_native_lu_discard_payload((void *)(UINTPTR_MAX - 1), bytes);
    assert(calls == before);
    deny = 1;
    memset(raw, 0x6b, bytes + 16384);
    errno = ERANGE;
    fsparse_native_lu_discard_payload(raw + 17, bytes);
    assert(errno == ERANGE);
#ifndef FS_TEST_DISCARD_DISABLED
    assert(calls == before + 1 && last_result == -1);
#else
    assert(calls == before);
#endif
    for (size_t i = 0; i < bytes + 16384; ++i) assert(raw[i] == 0x6b);
#endif
    free(raw);
    return 0;
}
