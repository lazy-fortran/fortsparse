#include "fsparse_native_amd_bounded.h"
static int calls;
int __real_fsparse_native_amd_bounded(int,int,const int*,const int*,int,int32_t*,int64_t,int*,int*);
int __wrap_fsparse_native_amd_bounded(int n,int nz,const int *ptr,const int *row,
    int low,int32_t *work,int64_t bytes,int *perm,int *inverse) {
    ++calls;
    return __real_fsparse_native_amd_bounded(n,nz,ptr,row,low,work,bytes,perm,inverse);
}
int fsparse_test_bounded_amd_calls(void) { return calls; }
