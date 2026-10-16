"""Max-product matching oracle (scipy LAPJV-style min-weight full bipartite
matching on c_ij = log(colmax_j) - log|a_ij| + 1) and diagonal statistics.
Usage: python3 -I wmatch.py matrix.mtx out_prefix [structural_preference.bin]
Writes out_prefix-wpref.bin (int32, 1-based preferred row per column)."""
import sys, json, time
import numpy as np
import scipy.io, scipy.sparse as sp
from scipy.sparse.csgraph import min_weight_full_bipartite_matching

path, out = sys.argv[1], sys.argv[2]
A = scipy.io.mmread(path).tocsc()
A.sum_duplicates()
n = A.shape[0]
absA = abs(A).tocsc()
absA.eliminate_zeros()  # numerically zero entries cannot be matched
colmax = np.asarray(absA.max(axis=0).todense()).ravel()
rowmax = np.asarray(absA.max(axis=1).todense()).ravel()

def diag_stats(pref, label):
    # pref[j] = 0-based row matched to column j
    v = np.asarray(abs(A)[pref, np.arange(n)]).ravel()
    rc = v / colmax
    rr = v / rowmax[pref]
    s = {"label": label,
         "structural_zero_or_numeric_zero": int(np.sum(v == 0)),
         "col_rel_lt_0.1": int(np.sum(rc < 0.1)),
         "col_rel_lt_1e-3": int(np.sum(rc < 1e-3)),
         "col_rel_lt_1e-8": int(np.sum(rc < 1e-8)),
         "row_rel_lt_0.1": int(np.sum(rr < 0.1)),
         "both_rel_lt_0.1": int(np.sum((rc < 0.1) & (rr < 0.1))),
         "col_rel_eq_1": int(np.sum(rc >= 1.0)),
         "min_col_rel": float(rc[v > 0].min()) if np.any(v > 0) else 0.0,
         "median_col_rel": float(np.median(rc)),
         "sum_log10_abs": float(np.sum(np.log10(np.where(v > 0, v, np.nan)))) if np.all(v > 0) else None}
    return s

res = {"n": n, "nnz": int(A.nnz), "nnz_numeric_nonzero": int(absA.nnz)}
res["identity"] = diag_stats(np.arange(n), "identity (no matching)")
C = absA.copy()
C.data = np.log(colmax[C.indices * 0 + np.repeat(np.arange(n), np.diff(C.indptr))]) - np.log(C.data) + 1.0
if len(sys.argv) > 4:  # given weighted preference (e.g. Fortran prototype, dual-certified)
    pref = np.fromfile(sys.argv[4], dtype="<i4").astype(np.int64) - 1
    res["weighted_source"] = sys.argv[4]
else:
    t0 = time.perf_counter()
    r, c = min_weight_full_bipartite_matching(C)  # rows r[k] matched to cols c[k]
    res["scipy_match_s_unqualified"] = time.perf_counter() - t0
    pref = np.empty(n, dtype=np.int64); pref[c] = r
assert len(np.unique(pref)) == n
res["weighted"] = diag_stats(pref, "max-product matching (scipy oracle)")
(pref + 1).astype("<i4").tofile(out + "-wpref.bin")
if len(sys.argv) > 3:
    sp_pref = np.fromfile(sys.argv[3], dtype="<i4").astype(np.int64) - 1
    assert len(np.unique(sp_pref)) == n
    res["structural"] = diag_stats(sp_pref, "native structural Hopcroft-Karp matching")
    res["weighted_vs_structural_same_rows"] = int(np.sum(sp_pref == pref))
print(json.dumps(res, indent=1))
