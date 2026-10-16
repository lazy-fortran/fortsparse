/* Release only dead, wholly owned LOW factor payload pages before free.
 * No allocator metadata, partial boundary pages or process-wide policy changes.
 * The caller must finish every read/copy before calling this advisory helper. */
#ifndef _DEFAULT_SOURCE
#define _DEFAULT_SOURCE
#endif
#include <stddef.h>
#include <stdint.h>
#if defined(__linux__) && !defined(FORTSPARSE_LU_DISCARD_DISABLED)
#include <errno.h>
#include <sys/mman.h>
#include <unistd.h>
#endif

void fsparse_native_lu_discard_payload(void *payload, size_t bytes)
{
#if defined(__linux__) && !defined(FORTSPARSE_LU_DISCARD_DISABLED)
    const size_t minimum = (size_t)1024 * 1024;
    if (payload == NULL || bytes < minimum) return;
    const int saved_errno = errno;
    const long result = sysconf(_SC_PAGESIZE);
    if (result > 0) {
        const uintptr_t page = (uintptr_t)result;
        const uintptr_t first = (uintptr_t)payload;
        if (bytes <= UINTPTR_MAX - first && page <= bytes / 2) {
            /* Strictly interior endpoints even for page-aligned payloads. */
            const uintptr_t begin = first + page - first % page;
            const uintptr_t last = first + bytes - 1;
            const uintptr_t end = last - last % page;
            if (begin < end)
                (void)madvise((void *)begin, end - begin, MADV_DONTNEED);
        }
    }
    /* Advisory failure never prevents the ordinary owner deallocation. */
    errno = saved_errno;
#else
    (void)payload;
    (void)bytes;
#endif
}
