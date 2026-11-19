# Changelog

## Unreleased
- Real symmetric indefinite matrices on the native CPU backend are factored by
  a multifrontal LDL^T with threshold 1x1/2x2 and delayed pivots and a
  pair-compressed AMD order; complex Hermitian and HPD matrices use it, or the
  supernodal Cholesky, through their real embedding. The scalar factors remain
  for host caps, LOW_MEMORY and NATIVE_REFERENCE.
- Removed `options%ldlt_blocked` and `options%ldlt_supernodal` and the
  `ldlt_blocked_panels` statistic.

## 0.1.0 (2026-11-14)
- Retained sparse real and complex factorizations behind a common backend API.
- Native CPU factorization, optional external backends and CUDA support for
  selected positive-definite systems.
- Implicit solve derivatives, memory controls and preconditioners.
- See [implementation status](docs/status.md) for tested paths and limitations.
