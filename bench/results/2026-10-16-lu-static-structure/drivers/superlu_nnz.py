"""SuperLU fill oracle mirroring FortSparse's general SuperLU settings
(dgstrf, COLAMD, DiagPivotThresh=1.0, no equilibration, no row permutation),
via scipy's bundled SuperLU. Usage: python3 -I superlu_nnz.py matrix.mtx"""
import sys, json
import numpy as np, scipy, scipy.io
from scipy.sparse.linalg import splu
A = scipy.io.mmread(sys.argv[1]).tocsc(); A.sum_duplicates()
n = A.shape[0]
known = np.sin(np.arange(1, n + 1)) + 0.125 * np.cos(3.0 * np.arange(1, n + 1))
b = A @ known
lu = splu(A, permc_spec="COLAMD", diag_pivot_thresh=1.0, options=dict(Equil=False, RowPerm="NOROWPERM"))
x = lu.solve(b)
r = b - A @ x
out = {"scipy": scipy.__version__, "nnzL_incl_diag": int(lu.L.nnz), "nnzU_incl_diag": int(lu.U.nnz),
       "nnzLU": int(lu.L.nnz + lu.U.nnz - n),
       "row_departures": int(np.sum(lu.perm_r != lu.perm_c)),
       "normwise_be": float(np.max(np.abs(r)) / (abs(A).sum(axis=1).max() * np.max(np.abs(x)) + np.max(np.abs(b))))}
print(json.dumps(out))
