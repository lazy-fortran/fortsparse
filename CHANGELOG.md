# Changelog

## Unreleased
- The native supernodal LU no longer keeps an assembly map after its first
  factorization: the map (`4(n+1)+12 nnz` bytes, 19 % of a retained factor on
  a coupled-MHD Jacobian) is built by the first `sparse_refactor` of the same
  analysis, so a factor-once caller saves it at no cost in time. Delayed
  pivots rebuild the layout from the factored matrix and no longer need it.
- A symmetric matrix can be passed as one triangle: `options%triangle =
  FORTSPARSE_TRIANGLE_LOWER` or `_UPPER` (LAPACK's uplo; the other triangle is
  not read) with an SPD, INDEFINITE or HERMITIAN hint, real matrices. The
  triangle is expanded to an exactly symmetric full matrix, so the bitwise
  symmetry the native LDL^T demanded of full storage no longer applies, and
  `sparse_refactor` keeps the ordering and symbolic analysis of the expansion.
  `sparse_inertia` reports the inertia of such factors.
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
