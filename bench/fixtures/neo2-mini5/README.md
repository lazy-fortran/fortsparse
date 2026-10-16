# NEO-2 mini5 numeric fixture

This is the existing NEO-2 five-by-five mini-example, reconstructed as exact
integer mathematical data. It is a real general nonsymmetric system with one
all-ones RHS. It exercises neither kinetic sparsity nor production memory.

- `matrix.mtx`: 1-based coordinate real general MatrixMarket, 22 entries.
- `rhs.mtx`: real array MatrixMarket, five-by-one.
- `oracle.json`: exact determinant and rational solution; decimal values are
  convenience approximations.
- `provenance.json`: pinned source identity, reconstruction recipe and hashes.
- `LICENSE.source`: original NEO-2 MIT license and attribution.

The oracle was obtained independently with Cramer's rule, evaluating each 5x5
determinant as its 120 signed permutation products using integer arithmetic.
Exact rational substitution gives Ax-b=0. This did not invoke UMFPACK,
FortSparse, BLAS or LAPACK. A consumer can independently reproduce the oracle
with exact rational elimination, or compare a dense pivoted solve and its
scaled residual. No implementation source from UMFPACK or other GPL components
is included. No production execution time or solver tolerance is inferred.
