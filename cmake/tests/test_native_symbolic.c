#ifndef _DEFAULT_SOURCE
#define _DEFAULT_SOURCE
#endif
#include <errno.h>
#include <stddef.h>
#include <stdint.h>
#if defined(__linux__) && !defined(FORTSPARSE_SYMBOLIC_MALLOC_ONLY)
#include <sys/mman.h>
#endif

void *fsparse_native_symbolic_alloc(size_t bytes);
void fsparse_native_symbolic_free(void *owner, size_t bytes);

int main(void)
{
    /* Exercise malloc and mapping-sized backing stores across reuse cycles,
       including odd byte sizes: munmap accepts the allocation's original size. */
    const size_t sizes[] = {4, 4092, 1024*1024-4, 1024*1024, 1024*1024+12};
    for (int repeat = 0; repeat < 3; ++repeat) {
        for (size_t test = 0; test < sizeof sizes / sizeof sizes[0]; ++test) {
            const size_t bytes = sizes[test];
            int32_t *rows = fsparse_native_symbolic_alloc(bytes);
            if (rows == NULL) return 1;
            for (size_t i = 0; i < bytes / sizeof *rows; ++i)
                rows[i] = (int32_t)(i ^ (size_t)repeat);
            for (size_t i = 0; i < bytes / sizeof *rows; ++i)
                if (rows[i] != (int32_t)(i ^ (size_t)repeat)) return 2;
#if defined(__linux__) && !defined(FORTSPARSE_SYMBOLIC_MALLOC_ONLY)
            const uintptr_t address = (uintptr_t)rows;
#endif
            fsparse_native_symbolic_free(rows, bytes);
#if defined(__linux__) && !defined(FORTSPARSE_SYMBOLIC_MALLOC_ONLY)
            if (bytes >= 1024*1024) {
                unsigned char residency;
                errno = 0;
                if (mincore((void *)address, 1, &residency) != -1 || errno != ENOMEM)
                    return 3; /* Released pages must no longer be mapped. */
            }
#endif
        }
    }
    if (fsparse_native_symbolic_alloc(0) != NULL) return 4;
    if (fsparse_native_symbolic_alloc(SIZE_MAX) != NULL) return 5;
    fsparse_native_symbolic_free(NULL, 0);
    fsparse_native_symbolic_free(NULL, 1024*1024);
    return 0;
}
