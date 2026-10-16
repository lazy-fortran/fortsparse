# Elimination-tree path-union column counts

`native_symbolic_counts` is project-owned graph-analysis code. It computes the
elimination forest and off-diagonal Cholesky column counts without enumerating
the factor's filled edges. This is a composition of standard elimination-tree,
DFS-order, lowest-common-ancestor and tree-difference techniques, not a claim
of a new sparse-factorization algorithm. There is no copied external code.

The caller supplies an already validated canonical, sorted, one-based full
symmetric CSC pattern in the chosen ordering. Values and definiteness are
irrelevant to structural counts. Numerical pivoting/cancellation is excluded;
the result is the usual symbolic Cholesky envelope. Roots have parent zero,
and every nonzero parent has a larger index than its child. Counts exclude
the diagonal. The helper checks argument extents and leaves full canonical
pattern validation with its owning backend.

For column `j`, let `S` be its original upper-triangular neighbors `i < j`.
Every such `i` is a descendant of `j` in the elimination tree. The filled
predecessors in row `j` are exactly the union of tree paths from `S` towards
`j`, excluding `j` itself. Equivalently, each vertex `v` on this union receives
one off-diagonal entry in column `v` of the lower factor. The standard
ancestor-compression tree construction establishes these parents without
walking every eventual filled path for every row.

Sort `S` by DFS preorder. Add one to `delta(i)` for every `i` in `S`, subtract
one at the LCA of each adjacent pair, and, when `S` is nonempty, subtract one
at `j`. Finally accumulate `delta` from children into parents. To verify the
formula, fix a descendant `v` of `j`. Its subtree occupies a contiguous DFS
interval. If it contains `r > 0` starts, exactly `r - 1` adjacent pairs have
their LCA inside that subtree: they are the consecutive starts inside the
interval. Its subtree difference sum is therefore `r - (r - 1) = 1`. If it
contains no starts, the sum is zero. At `j` and every ancestor of `j`, all
starts and adjacent LCAs contribute a net one, canceled by the subtraction at
`j`. Other subtrees contribute zero. Summing these indicators over all `j`
therefore gives precisely each off-diagonal factor column count.

The implementation uses an iterative DFS, an in-place heapsort per original
upper column, binary-lifting LCA and 64-bit signed differences. Lifting tables
include an explicit zero vertex at every level, avoiding unguarded sentinel
indexing and dependence on Fortran short-circuit behavior. Parent indices
increase, so vertex order itself is a valid bottom-up accumulation order.
Setup is bounded by `O(nnz(A) log n + n log n)` time and `O(n log n)` temporary
storage, independent of `nnz(L)`; all helper workspaces are released on return.
This does not remove the separate filled-pattern construction required by
the numerical factorization. It is setup topology, not a continuous FortSym
leaf or an operation in a numerical hot loop.

The focused test independently eliminates dense Boolean graphs by connecting
surviving neighbors into cliques. It compares both counts and parents for all
undirected graphs through five vertices and larger paths, cliques, stars,
ladders, disconnected blocks and deterministic sparse/dense graphs, including
reversed orderings. This executable oracle supports the argument without
substituting numerical probes for it.

For isolated bounds-checked correctness verification:

```sh
cmake -S test-cmake/native_symbolic -B /var/tmp/fortsparse-symbolic-check \
  -G Ninja -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_Fortran_FLAGS='-Wall -Wextra -fcheck=all -fbacktrace'
cmake --build /var/tmp/fortsparse-symbolic-check -j 2
ctest --test-dir /var/tmp/fortsparse-symbolic-check \
  -R '^test_native_symbolic_counts$' --output-on-failure
```

The owning library includes `cmake/fortsparse_native_symbolic.cmake`, calls
`fortsparse_configure_native_symbolic(fortsparse)` once, and registers the
focused test with `fortsparse_add_native_symbolic_test(fortsparse)` when testing
is enabled. Backend integration and threshold selection are separate changes.

The backend uses this count schedule for reordered matrices with at least
4096 columns. Small and low-fill natural-order cases retain the simple reach
counter. Filled row-index construction remains a separate pass. Current
five-trial 3D analysis medians improved by roughly 0.7 ms; complete timings
remain subject to independent numeric-work variability, so no whole-solve
speedup is inferred from that setup measurement alone.
