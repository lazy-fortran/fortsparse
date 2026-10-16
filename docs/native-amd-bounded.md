# Bounded AMD ordering for native LU

Explicit AMD ordering supports positive host allocation caps and `MEMORY_LOW`
for native real/complex general LU when the audited AMD 3.3.4 caller-workspace
ABI is available. Natural `DEFAULT`, `NATIVE_REFERENCE`, and existing SPD
ordering validation retain their previous policies. This is CPU ordering; row
pivoting, value scaling, sparse numerical factors and retained N/T/H solves are
unchanged. These checks establish correctness and array-payload accounting,
not performance admission or completion of the broader memory roadmap.

## Library boundary

The existing uncapped `amd_order` adapter remains unchanged. It owns opaque
vendor allocation, so complete peak bytes remain unknown and
`strict_memory_supported` remains false for that path. A positive host cap or
`MEMORY_LOW` selects the separate allocation-free `amd_2` adapter. Its public
header warns that inputs are unchecked and destroyed; FortSparse therefore
constructs and validates the complete private graph before calling AMD.

Configure checks require the reviewed 3.3.4 header, exact function signature,
32-bit C integer ABI and linkable `amd_2`/`amd_version` symbols. Every bounded
operation checks the runtime version too. Other header/runtime versions or a
missing interface return `BACKEND_UNAVAILABLE`; uncapped capability remains
available independently. There are no global SuiteSparse allocator hooks,
production preload or hardcoded installation paths.

The allocation contract trusts the audited upstream **AMD 3.3.4** implementation:
`amd_2 -> amd_postorder -> amd_post_tree` uses caller arrays and iterative tree
traversal, with no heap allocation in that closure. Matching version numbers
are ABI guards, not proof that a locally modified library has identical behavior;
a patched build claiming 3.3.4 must preserve that contract. The source and
compiled-library audit must be repeated before accepting another release.
AMD remains an externally linked BSD 3-clause dependency of the MIT core.

## Owned layout and cap admission

For canonical one-based CSC, let `K` be its stored off-diagonal entry count and
`m` the directed off-diagonal count in the structural union `A + A^T`.
The query validates CSC, counts the graph using allocation-free mirrored-entry
searches, and returns checked workspace and output byte counts. Values are
borrowed; numerical zeros remain structural entries. No value copy or dense
matrix is constructed. Counting costs O(K log(max_degree+1)); graph construction
uses an off-diagonal transpose followed by sorted union merges in O(n+K+m).

The caller owns one int32 arena and two n-element permutation outputs:

| Arrays | int32 count |
| --- | --- |
| transpose pointers/indices | `n+1+K` |
| AMD vectors except its two output vectors | `7*n` |
| adjacency and compression space | `iwlen` |
| permutation and inverse outputs | `2*n` |

Thus workspace bytes are `4*(K+8*n+1+iwlen)` and output bytes are `8*n`.
The ordinary bounded path uses `iwlen=m+floor(m/5)+n`, matching upstream's elbow
room. `MEMORY_LOW` uses the documented supported minimum `m+n`, saving exactly
`4*floor(m/5)` bytes at the cost of potentially more graph compression. The
transpose stays inside the arena during AMD and is counted throughout.
Default control and null diagnostic pointers add no owned diagnostic arrays.
All counts and conversions are checked before allocation, including the
symmetrized graph and workspace index limits.

Cap admission includes all old factors, retained user-vector storage and original
refinement CSC/work vectors, including a future refinement reservation if those
buffers have not yet been populated. Both output arrays and the arena are charged
before allocating them. After ordering, the arena and inverse are released;
the `4*n` permutation input stays charged through LU factorization and growth.
A positive exhausted cap never becomes LU's zero-means-unlimited argument.

`host_current_bytes` counts retained allocated array payload. `host_peak_bytes`
takes the maximum of actual ordering and factor/growth phases, including their
live external buffers; unused future refinement reservations are admission gates,
not observed allocated bytes. Fixed call frames, object/allocator metadata,
loader/BLAS/runtime costs and RSS remain separate. A known ordering-wrapper floor
is never substituted for the opaque uncapped path's complete peak.

A successful same-pattern/same-strategy factor retains its column-order identity.
Refactor copies that order into a charged `4*n` input buffer and skips graph
construction/AMD. Refused ordering setup preserves that identity for recovery;
a numerical factor attempt invalidates it until a successful factor is available.
Changed pattern or strategy recomputes ordering as needed. The original-matrix
checker and user vector retain their existing lifetimes across same-shape reuse.

## Verification

Focused tests independently reconstruct dense boolean `A + A^T` and inspect the
actual graph passed to AMD for asymmetric/mirrored, diagonal-only, empty-edge and
isolated-vertex cases. They check canonical rows, exact buffers, both elbow-room
policies, invalid dimensions/pointers/duplicates, permutation inverses and runtime
version refusal. A Linux/glibc test executable first verifies that its allocator
interposition reaches the loaded SuiteSparse shared-library boundary, then checks
zero heap allocations during twelve real `amd_2` calls. It does not modify
SuiteSparse configuration or install an interposer in production.

Real/complex tests use independent dense manufactured N/T/H and eight-RHS
solutions with numerical row-pivot adversaries. They cover exact/one-byte-under
caps, old-factor overlap, refusal/recovery, order-copy reuse without a new AMD
call, changed pattern/strategy, original-matrix refinement, retained user vectors,
raw borrowed CSC and both bounded policies. AMD-disabled and isolated
missing-entry/unreviewed-header configurations fail closed; the ordinary SPD
wrapper and uncapped behavior keep their existing regression checks.
