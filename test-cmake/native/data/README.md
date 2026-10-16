# Public synthetic spline fixture

`spline54.mtx` and `spline54-rhs.mtx` are the exact binary64 n=54, nnz=180
system assembled by the public NEO-2 `TEST/test_spline_sparse_assembly.f90`
first call to `COMMON/spline_cof.f90` (`splinecof3`, sw1=1, sw2=2).
Input nodes are [0,.08,.21,.37,.52,.70,.86,1]; y=sin(2*x)+.25*x*x,
lambda=1, zero boundary data and unit weight. No private plasma data is used.
The published-source checkout inspected is NEO-2 commit
2b8d5c0e9fea7226de1b0decc9951f7d4e769162. Fixture values were captured at the
native adapter boundary before factorization, without changing assembly.
MatrixMarket header percent signs are normalized; numeric values are unchanged.

The original observable is the cubic evaluated at the eight input nodes,
with absolute error <=1e-12. Coefficients a/b/c/d for interval i are solution
entries 7*i+[1,2,3,4] (i zero-based). The endpoint uses interval 6, h=.14.
This is a numerical interpolation assertion, not a rigorous certificate.

The matrix is badly conditioned (estimated infinity condition ~1.76e14).
Before the fix, native unrefined natural LU fails the node gate at 1.8083e-8,
although original-A residual infinity norm is 9.68e-14. Independent dense DGESV
and band DGBSV also fail (~3.90e-9); they cannot serve as exact centers here.
An independent 80-decimal-digit LU of the exact binary64 matrix/RHS and legacy
UMFPACK both pass at ~1.128e-13. Thus a tiny normwise backward error alone is
insufficient for this observable. External raw diagnostic evidence is retained
under /var/tmp/neo2-native-fortsparse/spline-fixture; the committed regression
uses only the public matrix/RHS and original analytical node data.
