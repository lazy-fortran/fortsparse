# Changelog

## Unreleased
- Real symmetric indefinite matrices on the native CPU backend are factored by
  a multifrontal LDL^T with threshold 1x1/2x2 and delayed pivots and a
  pair-compressed AMD order; complex Hermitian and HPD matrices use it, or the
  supernodal Cholesky, through their real embedding. The scalar factors remain
  for host caps, LOW_MEMORY and NATIVE_REFERENCE.
- Removed `options%ldlt_blocked` and `options%ldlt_supernodal` and the
  `ldlt_blocked_panels` statistic.
- Nested-dissection ordering (`FORTSPARSE_ORDER_ND`) is built in, also under
  fpm: a Fortran port of METIS 5.2.1 `METIS_NodeND` that returns METIS's
  orderings at METIS's speed. `FORTSPARSE_ENABLE_METIS` is now off by default
  and links METIS only as the test reference.

## 0.1.0 (2026-11-14)
- Retained sparse real and complex factorizations behind a common backend API.
- Native CPU factorization, optional external backends and CUDA support for
  selected positive-definite systems.
- Implicit solve derivatives, memory controls and preconditioners.
- See [implementation status](docs/status.md) for tested paths and limitations.
