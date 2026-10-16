#include <stddef.h>

static int fail_allocation;
void * __real_fsparse_native_symbolic_alloc(size_t bytes);
void fsparse_test_symbolic_fail(int enabled)
{
    fail_allocation = enabled;
}
void * __wrap_fsparse_native_symbolic_alloc(size_t bytes)
{
    return fail_allocation ? NULL : __real_fsparse_native_symbolic_alloc(bytes);
}
