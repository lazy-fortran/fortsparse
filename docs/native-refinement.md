# Native CPU refinement

Set `solver%refine=.true.` before factor/refactor to enable original-matrix
residual checking and iterative refinement for native CPU real/complex LU,
complex HPD Cholesky, real pivoted LDLT, and complex Hermitian/transpose-symmetric
pivoted LDL. This option is retained at factor
entry; changing it requires a factor/refactor call. This is numerical backward
error control, not a forward-error bound or a rigorous certificate.

The backend retains the original canonical full CSC indices and values and two
length-n vectors of the corresponding value type. Normal, transpose and adjoint
residuals use that original matrix, with conjugation only for the adjoint.
Caller mutation or deallocation of the original matrix after factorization does
not change the checker. Each solve saves its RHS before modifying the output,
including in-place solves and derivative RHS. Block solves process one RHS at a
time using the same bounded workspace.

The acceptance criterion is
`||b-op(A)x||_inf / (||op(A)||_inf ||x||_inf + ||b||_inf) <= 32 epsilon`.
A zero residual is accepted, including the zero RHS. The backend takes at most
three correction solves with the retained factors. Nonfinite iterates/residuals,
nonrepresentable matrix norms, or failure to meet the target return an explicit
failure (`SINGULAR` with a refinement diagnostic); that status does not prove
that the original matrix is mathematically singular. Failed output is not an
accepted solution. Exponent normalization avoids overflowing the product in the
backward-error denominator. Residual/correction arithmetic uses the retained
generated update/divide leaves. Factors remain untrusted candidates checked
against the original matrix.

`refinement_steps` reports completed corrections for a vector, or the maximum
across a successful block. `backward_error` is the final estimate, or the maximum
across a successful block. `-1` means unavailable, for example an unrefined solve
or a nonfinite residual. A finite three-step failure retains its failed error.
Requested/active teams, factor family and factor statistics retain their meanings;
original CSC population and norm preparation follow the existing `factor_s`
measurement; quota reservation remains part of its dispatch overhead.

For default-width integers, additional payload is
`4*(n+1+nnz) + value_bytes*(nnz+2*n)`, with value_bytes 8 real or 16 complex.
The implementation uses `storage_size` for the actual integer width. It reserves
this quota and any live `sparse_vector` payload before invoking the factors, so
factor growth/reuse copies cannot consume the refinement quota. Original storage
is populated after a successful factorization. Reserving future storage is
conservative for a fresh workflow: its observed peak can be smaller than the
cap needed to reserve that storage throughout a retained refactor's growth.
Same-shape storage is reused;
obsolete buffers are released before a changed shape/type is allocated. All
retained arrays and their actual overlap with factor growth contribute to
`host_current_bytes` and `host_peak_bytes`. These are allocation payloads;
allocator overhead and process RSS are separate. Zero residual budget remaining
never becomes the factor interface's zero-means-unlimited setting.

Refinement adds no solve-time allocations. The focused Linux linker-wrapper
oracle measures the additional retained allocation payload and verifies zero
linked allocation calls during correction-required solve, block and JVP paths.
It does not intercept private allocation calls inside shared runtime libraries.
The functional fixtures independently evaluate original dense residuals for all
operations and cover correction-required LU, a three-correction failure,
HPD/LDLT, in-place/block solves, implicit derivatives, changed caps/patterns,
caller mutation and retained user-vector lifetime. Neither those fixtures nor
payload measurements establish external-solver performance admission.
