"""Dense emulation of PARDISO-style static-structure LU numerics.

Usage: python3 -I static_probe.py matrix.mtx prefix matching ordering scaling \
           mode tau [tau ...]
  matching: weighted | structural   (prefix-wpref.bin / prefix-spref.bin)
  ordering: ND | AMD                 (prefix-<matching>-<ord>-order.bin and
                                      -snodes.bin written by symb_probe)
  scaling : mc64 | none   (mc64: duals from Bellman-Ford on the optimal
                           matching's difference constraints; checked)
  mode    : static | restricted   (restricted = partial pivoting inside the
                                   supernode diagonal block only)
  tau     : perturbation, pivot |p| < tau*||M||_inf replaced by sign*tau*||M||

The symmetric order already includes the etree postorder used for supernodes.
Dense storage only emulates arithmetic; zeros outside the symmetrized
structure stay exactly zero under (restricted) diagonal pivoting.  Residuals
and backward errors use the ORIGINAL matrix in float64 (no extended
precision available in numpy on arm64); labelled accordingly.
"""
import sys, json, time
import numpy as np
import scipy.io
from scipy.linalg import solve_triangular

path, prefix, matching, ordering, scaling, mode = sys.argv[1:7]
taus = [float(t) for t in sys.argv[7:]]
A = scipy.io.mmread(path).tocsc()
A.sum_duplicates()
n = A.shape[0]
known = np.sin(np.arange(1, n + 1)) + 0.125 * np.cos(3.0 * np.arange(1, n + 1))
b = A @ known
absA = abs(A).tocsc(); absA.eliminate_zeros()
colmax = np.asarray(absA.max(axis=0).todense()).ravel()
pref = np.fromfile(f"{prefix}-{'wpref' if matching == 'weighted' else 'spref'}.bin", dtype="<i4").astype(np.int64) - 1
inv = np.empty(n, np.int64); inv[pref] = np.arange(n)
out = {"n": n, "matching": matching, "ordering": ordering, "scaling": scaling, "mode": mode}

dr = np.ones(n); dc = np.ones(n)
if scaling == "mc64":
    coo = absA.tocoo()
    i, j, a = coo.row.astype(np.int64), coo.col.astype(np.int64), coo.data
    c = np.log(colmax[j]) - np.log(a)
    # c at matched entry of row i: (i, inv[i])
    matched_c = np.full(n, np.nan)
    sel = j == inv[i]
    matched_c[i[sel]] = c[sel]
    assert not np.isnan(matched_c).any(), "matching uses an entry absent from |A|>0"
    k = inv[i]                      # edge k -> j, weight c_ij - c_{i,k}
    w = c - matched_c[i]
    v = np.zeros(n); t0 = time.perf_counter(); it = 0
    while True:
        it += 1
        cand = v[k] + w
        vn = v.copy(); np.minimum.at(vn, j, cand)
        if np.array_equal(vn, v) or np.max(v - vn) < 1e-14:
            v = vn; break
        v = vn
        if it > n + 1:
            raise RuntimeError("negative cycle: matching not optimal")
    u = matched_c - v[inv]        # u_i = c_{i,inv i} - v_{inv i}
    slack = c - u[i] - v[j]
    out["dual_iterations"] = it
    out["dual_min_slack"] = float(slack.min())
    out["dual_matched_max_abs_slack"] = float(np.abs(slack[sel]).max())
    out["dual_s_unqualified"] = time.perf_counter() - t0
    dr = np.exp(u); dc = np.exp(v) / colmax  # row i of A, column j of A
    S = (A.multiply(dr[:, None])).multiply(dc[None, :]).tocsc()
    out["scaled_max_abs"] = float(abs(S).max())
    out["scaled_matched_min_abs"] = float(np.min(np.abs(S[pref, np.arange(n)])))
    out["log10_scale_range_rows"] = float(np.ptp(np.log10(dr)))
    out["log10_scale_range_cols"] = float(np.ptp(np.log10(dc)))
order = np.fromfile(f"{prefix}-{matching}-{ordering}-order.bin", dtype="<i4").astype(np.int64) - 1
starts = np.fromfile(f"{prefix}-{matching}-{ordering}-snodes.bin", dtype="<i4").astype(np.int64) - 1
assert np.array_equal(np.sort(order), np.arange(n))
# M = Q^T (Dr P A Dc) Q; matched-row-relabelled index q -> original row pref[q], column q.
Ms = A.multiply(dr[:, None]).multiply(dc[None, :]).tocsr()[pref[order], :].tocsc()[:, order]
assert np.all(Ms.diagonal() != 0)
snend = np.empty(n, np.int64)
bounds = list(starts) + [n]
for s in range(len(starts)):
    snend[bounds[s]:bounds[s + 1]] = bounds[s + 1]
normM = float(abs(Ms).sum(axis=1).max())
F0 = np.asfortranarray(Ms.toarray())
maxM = float(np.abs(F0).max())
del Ms

def factor(F, tau_abs):
    rowperm = np.arange(n); pert = swaps = thresh_viol = 0
    b0 = 0; s = 0
    while b0 < n:
        b1 = b0
        while s < len(starts) and (b1 == b0 or bounds[s + 1] - b0 <= 128):
            b1 = bounds[s + 1]; s += 1
        for kk in range(b0, b1):
            col = F[kk:, kk]
            cmax = np.max(np.abs(col))
            if mode == "restricted":
                e = snend[kk]
                r = kk + int(np.argmax(np.abs(F[kk:e, kk])))
                if r != kk:
                    F[[kk, r], :] = F[[r, kk], :]; rowperm[[kk, r]] = rowperm[[r, kk]]; swaps += 1
            p = F[kk, kk]
            if abs(p) < 0.1 * cmax:
                thresh_viol += 1
            if abs(p) < tau_abs:
                F[kk, kk] = tau_abs if p >= 0 else -tau_abs; pert += 1
            F[kk + 1:, kk] /= F[kk, kk]
            if kk + 1 < b1:
                F[kk + 1:, kk + 1:b1] -= np.outer(F[kk + 1:, kk], F[kk, kk + 1:b1])
        if b1 < n:
            F[b0:b1, b1:] = solve_triangular(F[b0:b1, b0:b1], F[b0:b1, b1:], lower=True, unit_diagonal=True)
            F[b1:, b1:] -= F[b1:, b0:b1] @ F[b0:b1, b1:]
        b0 = b1
    return rowperm, pert, swaps, thresh_viol

def make_solve(F, rowperm):
    def solve(rhs):
        t = (dr * rhs)[pref][order][rowperm]
        t = solve_triangular(F, t, lower=True, unit_diagonal=True, check_finite=False)
        t = solve_triangular(F, t, lower=False, check_finite=False)
        y = np.empty(n); y[order] = t
        return dc * y
    return solve

def errors(x):
    r = b - A @ x
    denom = abs(A) @ np.abs(x) + np.abs(b)
    nrm = float(np.max(np.abs(r)) / (abs(A).sum(axis=1).max() * np.max(np.abs(x)) + np.max(np.abs(b))))
    cw = float(np.max(np.abs(r) / np.where(denom > 0, denom, 1)))
    fw = float(np.max(np.abs(x - known)) / np.max(np.abs(known)))
    return nrm, cw, fw

out["runs"] = []
for tau in taus:
    F = F0.copy(order="F")
    t0 = time.perf_counter()
    rowperm, pert, swaps, viol = factor(F, tau * normM)
    tf = time.perf_counter() - t0
    run = {"tau": tau, "perturbed": pert, "in_supernode_swaps": swaps,
           "threshold0.1_violations": viol, "factor_s_unqualified": tf,
           "growth_maxU_over_maxM": float(np.abs(np.triu(F)).max() / maxM),
           "max_abs_L": float(np.abs(np.tril(F, -1)).max()),
           "min_abs_Udiag": float(np.min(np.abs(np.diag(F)))), "refinement": []}
    solve = make_solve(F, rowperm)
    x = solve(b)
    for step in range(0, 6):
        nrm, cw, fw = errors(x)
        run["refinement"].append({"step": step, "normwise_be": nrm, "componentwise_be": cw, "forward_rel": fw})
        if not np.all(np.isfinite(x)):
            break
        x = x + solve(b - A @ x)
    out["runs"].append(run)
    del F
print(json.dumps(out, indent=1))
