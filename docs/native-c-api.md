# Native C and C++ boundary

`include/fortsparse_native.h` exposes single-process real/complex general native
LU. `fs_native_create(0,&handle)` selects real, `1` selects complex. Initialize
`handle` to NULL. The context owns its factors; there is no global handle table.
`fs_native_free(&handle)` releases them and clears the pointer; repeated free is
safe. Foreign/stale/copied freed handles and concurrent calls to one context are
outside the contract. Independent contexts are reentrant.

CSC is zero-based, square, with I4 dimensions, sorted unique rows in each column,
and finite double values. The adapter validates before converting to Fortran's
one-based indices. Arrays are borrowed only during a call; subsequent caller
mutation/free is safe. Even an empty pattern requires non-NULL row/value pointers.
CSC conversion and output staging use temporary adapter-owned allocations with
checked allocation status; these are additional to backend storage/statistics.
This API does not expose a memory cap or claim its staging is backend-accounted.

`factor_d/z` factor afresh; `refactor_d/z` permits the owning backend to retain
symbolic work for an unchanged pattern. Every failed factor/refactor invalidates
old factors, including invalid input. A later valid factor/refactor can recover.
Default ordering is the native DEFAULT (natural) ordering; optional AMD/METIS and
all external sparse solver backends may be disabled in this API's minimal build.
`set_refine(handle,1)` enables native original-equation residual refinement on the
next factor. Changing the flag invalidates factors. It does not promise another
vendor's iteration sequence or a rigorous forward-error certificate.

`set_ordering(handle,FS_NATIVE_ORDER_AMD)` explicitly enables optional AMD
ordering. `FS_NATIVE_ORDER_DEFAULT=0` stays the default; `AMD=1` uses only the
optional BSD ordering library, with no external sparse factorization fallback.
`ordering_available(order)` returns status 0 for supported, 3 for unavailable,
and 2 for an invalid enum; it does not change any context diagnostic. ND=2 is
recognized but currently unavailable for this general-LU boundary. The setter
checks support before changing policy: invalid/unavailable requests preserve
policy and factors, successful changes invalidate factors, and successful
same-policy calls retain them. A valid setter clears the previous diagnostic.
This does not change refinement or add a memory-cap policy.

`solve_d/z` accepts N/T/H operations 0/1/2 and a contiguous column-major
`n*nrhs` block. It validates the dimension, finite input, operation, and finite
output; output changes only on success. Input/output may alias because the
solution is staged. Setup and solve allocation are boundary costs; this is not
an allocation-free hot-loop C interface. Return codes match `fortsparse_status`;
`last_error` copies the diagnostic with bounded NUL termination.

Complex pointers are C99 `double _Complex`, C++ `std::complex<double>` and
Fortran `complex(c_double_complex)`. No complex scalar crosses the ABI by value.
The C++ array real/imaginary representation follows the
[C++ specification](https://eel.is/c++draft/complex.numbers.general#4);
Fortran uses [ISO_C_BINDING intrinsic-type interoperability](https://gcc.gnu.org/onlinedocs/gfortran/Intrinsic-Types.html).
Actual C99 and C++ complex known-solution tests independently exercise N/T/H,
including unequal imaginary off-diagonal entries. A struct pretending to be a
C complex value and GNU complex return conventions are not used.

Build/test the minimal dependency configuration:

```sh
cmake -S . -B build-c-api -G Ninja -DCMAKE_BUILD_TYPE=Debug \
  -DFORTSPARSE_ENABLE_SUPERLU=OFF -DFORTSPARSE_ENABLE_UMFPACK_IPC=OFF \
  -DFORTSPARSE_ENABLE_MUMPS=OFF -DFORTSPARSE_ENABLE_AMD=OFF \
  -DFORTSPARSE_ENABLE_METIS=OFF
cmake --build build-c-api --target test_native_c_api test_native_cxx_api -j2
ctest --test-dir build-c-api -R '^test_native_c(xx)?_api$' --output-on-failure -j2
```

The tests cover independent dense-known-x block RHS, retained input lifetimes,
N/T/H, in-place output, real/complex refinement, refactor and changed patterns,
singular/invalid/nonfinite failures, recovery, shape/overflow checks, independent
contexts and repeated create/factor/free. With AMD enabled and disabled, they
also check ordering support/status, retained factors after rejected/same-policy
requests, invalidation on policy changes and actual AMD dispatch. No performance
admission is inferred.
