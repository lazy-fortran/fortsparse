/* Context-owned symbolic L row scratch. Large Linux mappings are released
 * directly instead of entering the process allocator's retained heap arena.
 * This does not change any process-wide allocator policy. */
#ifndef _DEFAULT_SOURCE
#define _DEFAULT_SOURCE
#endif
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#if defined(__linux__) && !defined(FORTSPARSE_SYMBOLIC_MALLOC_ONLY)
#include <sys/mman.h>
#define FSPARSE_SYMBOLIC_MAPPED 1
#endif

#define FSPARSE_SYMBOLIC_MAP_MIN ((size_t)1024 * 1024)

void *fsparse_native_symbolic_alloc(size_t bytes)
{
    if (bytes == 0 || bytes > (size_t)PTRDIFF_MAX) return NULL;
#ifdef FSPARSE_SYMBOLIC_MAPPED
    if (bytes >= FSPARSE_SYMBOLIC_MAP_MIN) {
        void *owner = mmap(NULL, bytes, PROT_READ | PROT_WRITE,
                           MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        return owner == MAP_FAILED ? NULL : owner;
    }
#endif
    return malloc(bytes);
}

void fsparse_native_symbolic_free(void *owner, size_t bytes)
{
    if (owner == NULL) return;
#ifdef FSPARSE_SYMBOLIC_MAPPED
    if (bytes >= FSPARSE_SYMBOLIC_MAP_MIN) {
        (void)munmap(owner, bytes);
        return;
    }
#else
    (void)bytes;
#endif
    free(owner);
}
