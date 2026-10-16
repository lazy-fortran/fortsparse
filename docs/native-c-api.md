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
checked allocation status. The handle-scoped memory contract below includes
these buffers together with backend storage.

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
This does not change refinement or memory policy.

`set_memory(handle, host_limit_bytes, device_limit_bytes, policy)` sets 64-bit
per-handle array-payload caps. Zero means unlimited; policies are
`FS_NATIVE_MEMORY_AUTO=0` and RAM-only `FS_NATIVE_MEMORY_LOW=1`. Negative limits
or invalid policies preserve factors and policy, returning invalid status.
Successful changes free/invalidate factors; same-policy calls preserve them.
A valid setter clears the last diagnostic and allocation-failure byte fields.
This boundary is CPU native general LU only: device allocations are zero, and
retaining a device limit does not enable GPU factorization. External
factorization backends are not selected.

The host cap covers conversion arrays, retained factors, numerical growth,
refinement reservations, ordering workspace, solve output staging and their
simultaneous old/new lifetimes. For n columns and m entries, conversion payload
is `4*(n+1)+(4+value_bytes)*m`, where value_bytes is 8 real or 16 complex.
Before conversion, retained backend payload plus this entire buffer must fit;
the backend then receives the remaining positive quota. An exhausted positive
quota never becomes the backend's zero-means-unlimited setting. Solve staging
adds `value_bytes*n*nrhs` to retained payload. C pointers map to contiguous
Fortran views, avoiding a hidden RHS copy. Resource exhaustion returns status 8
before the allocation; operating-system allocation failure returns status 6.
Failed factor/refactor invalidates factors; failed solve leaves factors and
caller output intact. Raising limits and retrying can recover.

`get_memory_stats(handle,&stats)` is diagnostic-only and preserves `last_error`.
It reports configured limits/policy, current payload, the last factor/solve peak,
device bytes, failed-phase required/available bytes and strict-cap support.
The `peak_is_bound` flag distinguishes conservative backend admission/history
bounds from the exact solve-staging interval; unknown opaque AMD peaks remain
`-1`. Never sum peaks from separate operations. Under a strict cap, that enforced
limit also bounds any conservative historical counter. Required/available bytes
describe the failed allocation phase; zero fields mean no retained byte failure
diagnostic. This API does not predict pivot-dependent fill or process RSS.
Audited AMD caller-workspace execution supports positive caps/LOW; an unavailable
audited ABI fails closed. Opaque uncapped AMD stays available without a strict
memory claim. Payload excludes caller A/B/X, fixed context/backend descriptors,
allocator/page rounding, runtime/code, stacks and process/job overhead. Each
handle has its own cap; concurrent handles must be added by the application.

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

The focused `test_native_c_api_memory` independently traces static-link
allocations and injects every conversion/backend/refinement allocation failure.
Fresh real/complex peaks match reported payload; compound refactor traces are
compared with the maximum of matching operation bounds, not one unrelated peak.
Solve traces require exactly one explicit staging allocation, including aliased
N/T/H blocks. Tests also cover one-byte cap boundaries, zero/unlimited and tiny
caps, unchanged failure output, independent handles and lifecycle recovery.
These allocation fixtures establish payload accounting, not RSS predictability
or performance admission.
