# Native CPU refinement

Set `solver%refine=.true.` before factor/refactor to enable original-matrix
residual checking and iterative refinement for native CPU real SPD Cholesky,
real/complex LU,
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
in the general-family adapter, original CSC population and norm preparation follow
the existing `factor_s` measurement; quota reservation remains part of its dispatch
overhead.

For default-width integers, additional payload is
`4*(n+1+nnz) + value_bytes*(nnz+2*n)`, with value_bytes 8 real or 16 complex.
The implementation uses `storage_size` for the actual integer width. It reserves
this quota and any live `sparse_vector` payload before invoking the factors, so
factor growth/reuse copies cannot consume that quota. This reservation is an
allocation gate, not a measurement of simultaneously live arrays. In the general-family
adapter, a fresh workflow allocates/populates original-matrix storage after
factorization;
a same-shape refactor keeps that storage live during factor growth. Obsolete
buffers are released before a changed shape/type is allocated.

After successful factorization, `host_current_bytes` and `host_peak_bytes`
report tracked allocated-array payload, including actual overlap with factor growth. They exclude unused
reservation, allocator overhead and process RSS. Consequently, a fresh observed
`host_peak_bytes` alone does **not** guarantee that a same-pattern refactor fits.
The conservative reservation gate can also require more than a fresh run's
observed live peak. In the three-by-three complex fixture, the refinement quota
is 292 bytes: refactoring with only the fresh observed peak as the cap returned
`RESOURCE_LIMIT`; adding that exact quota allowed the tested retained workflow.
The tight-cap failure is retained as evidence, not removed from the test history.
Changed numerical values can alter pivots, fill and growth even with the same
input pattern, so that fixture's accepted cap is not a general prediction.
Zero residual budget remaining never becomes the factors' zero-means-unlimited
setting.

Refinement adds no solve-time allocations. The focused Linux linker-wrapper
oracle measures the additional retained allocation payload and verifies zero
linked allocation calls during correction-required solve, block and JVP paths.
It does not intercept private allocation calls inside shared runtime libraries.
The functional fixtures independently evaluate original dense residuals for all
operations and cover correction-required LU, a three-correction failure,
HPD/LDLT, in-place/block solves, implicit derivatives, changed caps/patterns,
caller mutation and retained user-vector lifetime. Neither those fixtures nor
payload measurements establish external-solver performance admission.

Real SPD refinement supports natural-order and audited explicit-AMD CPU
simplicial and blocked Cholesky.
`ORDER_DEFAULT` selects natural order when refinement is requested, including
uncapped inputs that might otherwise select automatic AMD. Explicit AMD uses
the shared bounded AMD 3.3.4 workspace even without a positive cap; unavailable
ABI, nested dissection and GPU refinement return `BACKEND_UNAVAILABLE`.
AMD correction RHS and results are permuted through existing numeric work,
while residuals remain in the original matrix/RHS order. The ordinary unrefined planner and
execution paths retain their existing behavior. Original CSC and two real vectors
are reserved and populated before SPD analysis, so planning scratch and blocked
construction overlap are included in the host budget. Same-pattern refactor
refreshes the retained original values without reallocating the checker; toggling
refinement or changing caps/patterns rebuilds through the ordinary factor path.
Original norm preparation contributes to factor-entry overhead; scalar and block
solves allocate no checker storage. These SPD checks use the same original-CSC
norm and residual procedures as the native general adapter.

`sparse_analyze` supports exact natural compact CPU SPD layout planning with
refinement for `LOW_MEMORY` or `NATIVE_REFERENCE`. Its persistent payload includes
the original CSC and two vectors; its peak covers their overlap with graph
planning and factor storage. The same checked quota calculation is used by the
factor adapters. Internal SPD analysis already has that quota charged and does
not add it a second time. RHS count changes do not add retained vector storage:
block refinement processes one RHS at a time.

On an already factored solver, preflight shares the current host cap with all
retained payload and verifies space for the proposed workflow, preserving the
existing overlap policy. A refused forecast leaves those factors and their live
allocation statistics unchanged. Forecast bytes are not reported as allocated
storage. The plan describes array layout, not positive definiteness, conditioning,
allocator overhead or RSS. Natural blocked `AUTO`, GPU and explicit ordering
preflight remain unavailable; unrefined exact-layout numbers are unchanged.
