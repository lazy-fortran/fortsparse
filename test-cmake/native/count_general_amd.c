static int calls, nonidentity;
int __real_fsparse_native_order_amd_general(int, int, const int *, const int *, int *, int *);
int __wrap_fsparse_native_order_amd_general(int n, int nz, const int *ptr,
                                           const int *row, int *q, int *inverse) {
    ++calls;
    int code = __real_fsparse_native_order_amd_general(n, nz, ptr, row, q, inverse);
    if (!code) for (int i = 0; i < n; ++i) if (q[i] != i+1) nonidentity = 1;
    return code;
}
int fsparse_test_general_amd_calls(void) { return calls; }
int fsparse_test_general_amd_nonidentity(void) { return nonidentity; }
int fsparse_test_general_amd_available(void) {
#ifdef FORTSPARSE_TEST_AMD
    return 1;
#else
    return 0;
#endif
}
