# fortsparse roadmap

`fortsparse` owns explicit sparse storage, sparse products, and pluggable
direct solves. It is a numerical dependency for compact-support covariance and
Markov-precision models; it does not own GP semantics, autodiff source
transformation, or optimizer loops.

## Completion rules

Every item requires an independent behavioral oracle, documented licensing
boundary, focused tests, and a reproducible benchmark where performance is
relevant. A test that only inspects repository state is not an oracle. Real and
complex derivative products must satisfy finite-difference and adjoint
identities. GPL SuiteSparse code remains isolated in the existing helper
process boundary.

## Work order

- [x] Establish MIT library ownership, SuperLU backend, GPL-isolated UMFPACK
  helper, CSC construction, status API, and solver lifetime rules.
- [x] Add real and complex solve JVP/VJP products with independent solve and
  adjoint-oracle tests.
- [ ] Add an explicit CSR/ELL view for matrix-free row-owned products without
  changing the CSC direct-solve API.
- [ ] Add batched sparse MVM/matmat products with resident OpenACC correctness
  tests and storage/nonzero diagnostics.
- [ ] Add compact-support kernel assembly helpers only where the sparsity
  pattern is explicit; keep kernel formulas in `fortml`.
- [ ] Benchmark CSR/ELL products against dense PyTorch and KeOps on matched
  float64 compact-support workloads through `fortml-bench`.
- [ ] Add optional `nvfortran`/CUDA sparse backends only after a complete
  transfer-inclusive workload beats the resident CPU/CSR baseline.
- [ ] Publish a versioned API and benchmark report for downstream `fortml`.

## Dependency boundary

`fortnum` owns dense and iterative numerical primitives, `fortml` owns GP
observation/kernel semantics, and `fortopt` owns optimization. `fortsparse`
may expose products and solves to all three through public status-bearing
interfaces, but it must not depend on model or optimizer modules.

## Research record

The ignored `.provenance/` directory contains exact upstream sparse-library
revisions, solver papers, license checks, and benchmark metadata. The source
tree contains no third-party GPL headers; the UMFPACK helper remains a separate
GPL executable.
