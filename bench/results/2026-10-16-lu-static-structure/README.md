# Static symmetrized-structure LU diagnostic for the KIN6D Jacobians

Scope: correctness and structure only. This run was on faepmac1, an Apple M3 Ultra running macOS 26.6.2, which is not the timing reference, and the host was shared with my concurrent runs.

- Every time below is labelled "unqualified". None of them is performance or Pareto evidence.
- Drivers are retained in `drivers/`. Raw outputs, the scratch build and the
  matrices stay outside Git: outputs in mailuefterl
  `~/code/sparse-dev/benchmarks/2026-10-16-lu-static-structure/`, matrices at the
  hashed source paths below.

## Inputs, provenance and hashes

| File | SHA256 | Source |
|---|---|---|
| data/kin6d-5525.mtx | 707d5e280ed1f345f12581a1086be4bfe42ff00e5a1291e87b3176b4b5288078 | mailuefterl `…/nonnested-csc-probe/lu-source-structure/` |
| data/kin6d-5525-nd-order.bin | dc1e09b22da053d16bcc79698bd178282dfb7f14d453d28b83d232a3612aa25d | same directory (retained; not used, see note) |
| data/kin6d-nonnested-csc-19989.mtx | 13ffe01f42f89dcefd7a692a9d0a725e0983b84fbe7f2d4445f3ab5cb447b1e5 | mailuefterl `/var/tmp/` (path recorded in `nonnested-csc-probe/matrices.json` and `19989-diagnostic/receipt.json`) |

- Both hashes match `matrices.json` / the manifest.
- n=5,525 has 203,054 entries, 7 of them explicit zeros. n=19,989 has 801,006 entries, 6 explicit zeros.
- The pattern is 95.7% structurally symmetric: 95.67% and 95.61% of off-diagonal entries have a transpose partner. nnz(A+Aᵀ) is 211,638 and 835,386.

Build and toolchain:
- FortSparse `main` 16e80c3 (clean tree), built with gfortran 16.1.0 (Homebrew GCC) at `-O3`, configured as requested in `build/`.
- Library `build/lib/libfortsparse.a`: hash in the list at the end.
- Dependencies: METIS 5.1.0, AMD 3.3.4 (SuiteSparse), SuperLU 7.0.1, OpenBLAS 0.3.33. Python oracles use scipy 1.18.0 and numpy 2.5.1 (Accelerate BLAS).

## Drivers

All are in `src/`.

- **`symb_probe.f90`.** Matching (none, native Hopcroft–Karp, or weighted from file) × production `fsparse_native_order_amd_general` / `fsparse_native_order_general` (METIS ND) on PmA. Then:
  - symmetrized pattern;
  - `native_symbolic_counts`, cross-checked against an independent ancestor-walk etree and count oracle (all 12 pattern/order pairs per matrix agree exactly);
  - etree postorder;
  - replica of the SPD supernode merge rule (`supernodes.F90:280-322`, expressed through parent/count);
  - LU flop and payload accounting.
  - It also writes the final order and the supernode starts.
- **`mpm_probe.f90`.** Independent maximum-product matching prototype (Dijkstra shortest augmenting path on c_ij = log colmax_j − log|a_ij|, with duals). It checks the LP certificate.
- **`lu_probe.f90`** together with `obs/fortsparse_native_lu_observer.f90`. The observer is a copy of the generated `src/native/fortsparse_native_lu.f90` (sha 785070b2…) with only `p,pinv,q,ln,un` made public (2-line diff). It runs the production scalar LU exactly as `general_factor_real` does (pruning on, `row_preference` = matching), counts pivot departures `p(k) ≠ preferred_row(q(k))`, and computes original-COO real128 backward errors. Mode 2 runs SuperLU through the public API.
- **`wmatch.py`.** Diagonal statistics; scipy matching oracle for n=5,525.
- **`static_probe.py`.** MC64-dual scaling (Bellman–Ford on the optimal matching's difference constraints, with feasibility checks), followed by a **dense emulation** of the static-structure numerics in the exact final order:
  - right-looking with restricted (within-supernode diagonal block) or no pivoting;
  - PARDISO-style perturbation |p| < τ‖Ã‖∞ → sign·τ‖Ã‖∞;
  - refinement against the ORIGINAL A, in float64.
  - Dense storage only emulates the arithmetic: structural zeros of the symmetrized structure stay exact zeros under these pivots.
- **`superlu_nnz.py`.** scipy-bundled SuperLU with FortSparse's general settings (`dgstrf`, COLAMD, DiagPivotThresh=1, no equilibration or row permutation). FortSparse's SuperLU backend does not report factor nnz (`stats%factor_nnz = -1`).

Commands (from `out/`; `S` = this directory):

```
$S/mpm_probe  M  out/kNNNN                                  # weighted matching (n=19,989 used it)
python3 -I src/wmatch.py M out/kNNNN [spref.bin [wpref.bin]] # diag stats (+ scipy matching for 5,525)
$S/symb_probe M out/kNNNN out/kNNNN-wpref.bin               # symbolic tables, order/snode files
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout 60 $S/lu_probe M <backend> <matching> <ordering> 0.1 [wpref]
python3 -I src/static_probe.py M out/kNNNN <weighted|structural> <ND|AMD> <mc64|none> <restricted|static> 1e-13 [1e-8]
python3 -I src/superlu_nnz.py M
python3 src/summarize_symb.py out/symb-*.log ; python3 src/summarize_static.py out/static-*.json
```

Raw outputs: `out/symb-{5525,19989}.log`, `out/wmatch-*.json`, `out/lu-{5525,19989}.log`, `out/static-*.json` (+ `.time`), `out/superlu-*.json`.

## (d) Matching feasibility: diagonal quality

Measured as |a_{σ(j)j}| relative to column max (col) or row max (row).

| n | matching | zero diag | col<0.1 | row<0.1 | both<0.1 | min col ratio | Σlog10\|diag\| |
|---|---|---:|---:|---:|---:|---:|---:|
| 5,525 | identity | 737 | 1,474 | 1,474 | 737 | 0 | n/a |
| 5,525 | native Hopcroft–Karp | 0 | 2,002 | 2,002 | 1,056 | 1.1e-7 | 7,059.5 |
| 5,525 | max-product | 0 | 1,474 | 1,474 | **0** | 2.5e-5 | **8,667.7** |
| 19,989 | identity | 2,849 | 5,697 | 5,698 | 2,849 | 0 | n/a |
| 19,989 | native Hopcroft–Karp | 0 | 7,519 | 7,520 | 3,644 | 5.9e-10 | 29,539.1 |
| 19,989 | max-product | 0 | 5,697 | 5,698 | **0** | 2.0e-5 | **35,228.5** |

MPM prototype (`mpm_probe`):

| n | free after greedy | heap pops | Σc = Σu+Σv | min reduced cost | max matched reduced | time (unqualified) |
|---|---:|---:|---|---:|---:|---:|
| 5,525 | 1,474 | 17,476 | 12262.628979608347 / …346 | −2.4e-15 | 2.4e-15 | 7.6 ms |
| 19,989 | 5,700 | 76,254 | 50881.30152561534 / …558 | −4.9e-15 | 4.9e-15 | 31 ms |

- For n=5,525 the prototype matching is identical row-for-row to scipy's `min_weight_full_bipartite_matching` (5,525/5,525).
- For n=19,989, scipy ran for more than 3.5 minutes without finishing and was killed. The dual certificate establishes optimality instead.
- After MC64-dual scaling (`static_probe`): max|DrADc| = 1+2.7e-15 and the matched minimum is 1−2.3e-15. The log10 scale ranges are 4.6 for rows and 7.6 for columns (n=5,525), and 4.8 and 9.1 (n=19,989).
- The col<0.1 and row<0.1 sets after MPM are disjoint. This is a badly scaled mixed gauge/voltage system, not a structurally weak diagonal, so the raw ratios are scale artefacts. Scaling removes them.

## (a) Static symbolic structure of the symmetrized matched matrix

Notes on the table:
- nnz(L) includes the diagonal; nnz(L+U) sym = 2nnz(L) − n.
- LU flops = Σ(c + 2c²), with c the off-diagonal column count.
- Supernodes follow the SPD rule on the etree-postordered order (relaxed). "exact" means fundamental chains only.
- "pad" = padded/exact L+U entries. "payload" = 8·padded L+U + 4·Σm.
- Flop shares are computed on padded supernodal flops. "GEMM" is the share of external Schur updates 2(m−w)²w in diag + TRSM + GEMM; "k≥8" restricts it to width≥8 sources.

n = 5,525

| matching | order | zero diag | nnz(L) | nnz(L+U) sym | LU flops | snodes exact given/post | snodes relaxed | max w | max front | pad | payload MB | flops w≥8 | flops w≥16 | GEMM | GEMM k≥8 | widths 1,2-3,4-7,8-15,16-31,32-63,64-128 |
|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| none | AMD | 737 | 782,510 | 1,559,495 | 4.24e8 | 1933/1932 | 1783 | 128 | 477 | 1.013 | 12.8 | 0.978 | 0.945 | 0.763 | 0.742 | 1346,224,118,38,34,12,11 |
| none | ND | 737 | 421,872 | 838,219 | 9.9e7 | 1841/1822 | 1475 | 64 | 233 | 1.036 | 7.1 | 0.945 | 0.842 | 0.762 | 0.712 | 1047,98,167,86,41,29,7 |
| structural | AMD | 0 | 532,926 | 1,060,327 | 1.68e8 | 1479/1478 | 1208 | 64 | 325 | 1.015 | 8.8 | 0.930 | 0.914 | 0.794 | 0.729 | 640,182,259,41,60,17,9 |
| structural | ND | 0 | 559,194 | 1,112,863 | 1.77e8 | 1497/1497 | 1133 | 64 | 316 | 1.037 | 9.4 | 0.965 | 0.899 | 0.791 | 0.758 | 641,191,111,111,39,32,8 |
| weighted | AMD | 0 | 525,058 | 1,044,591 | 1.82e8 | 1511/1510 | 1184 | 64 | 379 | 1.010 | 8.6 | 0.965 | 0.923 | 0.811 | 0.779 | 739,57,110,187,74,8,9 |
| weighted | ND | 0 | 488,710 | **971,895** | 1.32e8 | 1536/1536 | 1114 | 64 | 266 | 1.034 | 8.2 | 0.968 | 0.864 | 0.773 | 0.744 | 748,51,82,152,45,27,9 |

n = 19,989

| matching | order | zero diag | nnz(L) | nnz(L+U) sym | LU flops | snodes exact given/post | snodes relaxed | max w | max front | pad | payload MB | flops w≥8 | flops w≥16 | GEMM | GEMM k≥8 | widths |
|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| none | AMD | 2849 | 4,874,353 | 9,728,717 | 5.07e9 | 6148/6147 | 5718 | 128 | 1065 | 1.007 | 79.3 | 0.992 | 0.973 | 0.846 | 0.838 | 4066,1013,261,172,100,59,47 |
| none | ND | 2849 | 2,073,823 | 4,127,657 | 7.4e8 | 5734/5699 | 4350 | 64 | 473 | 1.033 | 34.7 | 0.974 | 0.921 | 0.816 | 0.792 | 2799,399,517,311,185,119,20 |
| structural | AMD | 0 | 2,957,259 | 5,894,529 | 1.65e9 | 4328/4327 | 3440 | 128 | 655 | 1.009 | 48.2 | 0.968 | 0.963 | 0.790 | 0.761 | 1280,701,1062,83,220,64,30 |
| structural | ND | 0 | 2,848,224 | 5,676,459 | 1.43e9 | 4418/4413 | 3132 | 128 | 612 | 1.021 | 47.0 | 0.982 | 0.948 | 0.786 | 0.770 | 1280,719,453,378,172,108,22 |
| weighted | AMD | 0 | 2,549,913 | 5,079,837 | 1.29e9 | 4110/4109 | 3381 | 64 | 609 | 1.005 | 41.5 | 0.975 | 0.956 | 0.851 | 0.828 | 1604,295,618,505,282,49,28 |
| weighted | ND | 0 | 2,481,498 | **4,943,007** | 1.07e9 | 4631/4631 | 3056 | 64 | 531 | 1.023 | 41.0 | 0.981 | 0.928 | 0.830 | 0.813 | 1626,254,365,482,175,123,31 |

The "none" rows have structurally zero diagonals (737 and 2,849). Static pivoting is therefore impossible for them; they appear for structure reference only.

## (b) Existing native scalar LU, threshold 0.1 (current revision on this Mac; 60 s cap)

| n | matching | order | factor_nnz | ln / un | departures | normwise / cw BE (real128, no refinement) | fwd err | factor s (unqualified) | owned peak B |
|---|---|---|---:|---|---:|---|---:|---:|---:|
| 5,525 | structural | AMD | 4,449,636 | 1,566,238 / 2,877,873 | 1,191 | 9.2e-16 / 3.7e-15 | 5.1e-10 | 0.38 | 101,851,224 |
| 5,525 | structural | ND | 4,233,464 | 1,482,493 / 2,745,446 | 1,059 | 3.5e-16 / 3.5e-15 | 5.8e-10 | 0.32 | 101,851,224 |
| 5,525 | none | ND | 2,728,370 | 1,146,169 / 1,576,676 | 1,230 | 7.3e-16 / 4.5e-15 | 3.3e-10 | 0.16 | 53,092,440 |
| 5,525 | none | AMD | 7,964,640 | 3,003,177 / 4,955,938 | 1,583 | 2.5e-15 / 2.4e-14 | 1.6e-9 | 1.12 | 189,145,176 |
| 5,525 | weighted (MPM pref) | ND | 2,000,896 | 992,285 / 1,003,086 | 76 | 2.2e-16 / 2.4e-15 | 9.3e-10 | 0.11 | 43,262,040 |
| 5,525 | weighted (MPM pref) | AMD | 1,884,566 | 1,125,885 / 753,156 | 146 | 1.8e-16 / 1.3e-15 | 1.4e-10 | 0.10 | 45,228,120 |
| 19,989 | none | ND | 39,516,398 | 9,614,117 / 29,882,292 | 4,732 | 8.7e-16 / 6.1e-15 | 4.5e-9 | 5.18 | 709,206,360 |
| 19,989 | structural | ND | 35,531,571 | 10,678,713 / 24,832,869 | 4,023 | 3.9e-15 / 1.6e-14 | 2.2e-9 | 5.49 | 809,869,656 |
| 19,989 | structural | AMD | 57,455,407 | 13,291,098 / 44,144,320 | 4,480 | 1.5e-13 / 5.8e-13 | 6.0e-8 | 11.20 | 1,413,849,432 |
| 19,989 | weighted (MPM pref) | ND | 20,567,718 | 8,775,032 / 11,772,697 | 327 | 8.5e-16 / 5.3e-15 | 3.1e-9 | 2.09 | 435,528,024 |

- Nothing was censored: every run finished inside the 60 s cap.
- Earlier mailuefterl diagnostics measured these on other revisions and METIS builds: 5,525 unmatched ND 3,334,519; 19,989 unmatched ND 31,749,631; matched ND 49,039,912; matched AMD timed out at 15 s. Their factor_nnz differences are attributed to a different METIS build/revision and not investigated here. The structural conclusion does not depend on them.
- `departures` counts all steps whose pivot row differs from the preferred row, structural singletons included.
- The weighted "MPM pref" rows feed the MPM matching as a row preference into the unchanged scalar LU, which keeps its own row max-abs scaling. This is a side observation, not the proposed method.

## (c) SuperLU

- Through FortSparse's in-process backend: completed, but factor_nnz is not reported (−1).
  - n=5,525: normwise / cw BE 6.3e-16 / 9.7e-11, with no refinement.
  - n=19,989: normwise / cw BE 1.9e-15 / 7.2e-10.
- With the same settings through scipy-bundled SuperLU: nnz(L)+nnz(U)−n = **2,215,127** (n=5,525) and **15,261,749** (n=19,989). Normwise BE is 5.3e-16 and 1.7e-15.

## Numerical viability of static pivoting (dense emulation in the exact final order)

Notes on the table:
- τ is the perturbation relative to ‖Ã‖∞.
- "viol0.1" counts pivots below 0.1× the max of their remaining column, i.e. the ones global threshold pivoting would have rejected.
- Growth = max|U| / max|Ã|.
- BE uses the original A in float64: normwise / componentwise, with "ref k" the value after k refinement steps.

| n | matching / order / scaling | pivoting | τ | perturbed | in-block swaps | viol0.1 | growth | max\|L\| | BE ref0 | BE ref1 | BE ref3 | fwd ref3 |
|---|---|---|---|---:|---:|---:|---:|---:|---|---|---|---:|
| 5,525 | weighted/ND/mc64 | restricted | 1e-13, 1e-8 | 0 | 1,184 | 19 | 131 | 131 | 6.6e-16 / 1.9e-13 | 1.1e-16 / 6.5e-16 | 8.6e-17 / 4.1e-16 | 2.0e-10 |
| 5,525 | weighted/ND/mc64 | none (static) | 1e-13, 1e-8 | 0 | 0 | 86 | 562 | 7,960 | 4.2e-15 / 7.9e-14 | 8.6e-17 / 3.7e-16 | 6.4e-17 / 5.1e-16 | 4.9e-10 |
| 5,525 | weighted/AMD/mc64 | restricted | 1e-13, 1e-8 | 0 | 1,337 | 83 | 1.27e4 | 5,620 | 4.7e-15 / 2.1e-13 | 8.6e-17 / 4.9e-16 | 8.6e-17 / 6.3e-16 | 6.5e-10 |
| 5,525 | weighted/AMD/mc64 | none | 1e-13, 1e-8 | 0 | 0 | 144 | 8.8e3 | 1,380 | 2.1e-15 / 1.6e-13 | 6.4e-17 / 4.3e-16 | 8.6e-17 / 4.9e-16 | 4.9e-10 |
| 5,525 | weighted/ND/unscaled | restricted | 1e-13 | 0 | 2,025 | 823 | 0.99 | 1.2e4 | 5.0e-15 / 4.2e-11 | 1.3e-16 / 4.8e-16 | 8.6e-17 / 5.6e-16 | 4.6e-10 |
| 5,525 | weighted/ND/unscaled | restricted | 1e-8 | **182** | 2,023 | 822 | 0.99 | 1.2e4 | 2.9e-9 / 4.9e-4 | 7.8e-10 / 4.4e-4 | 4.7e-10 / 3.6e-4 | 0.8 |
| 5,525 | weighted/ND/unscaled | none | 1e-13 | 0 | 0 | 1,185 | 14.6 | 4.0e4 | 2.4e-15 / 5.0e-14 | 8.6e-17 / 4.5e-16 | 6.4e-17 / 4.4e-16 | 2.3e-10 |
| 5,525 | weighted/ND/unscaled | none | 1e-8 | **184** | 0 | 1,176 | 14.6 | 4.0e4 | 1.7e-9 / 4.9e-4 | 8.0e-10 / 4.4e-4 | 3.6e-10 / 3.6e-4 | 0.8 |
| 5,525 | structural/ND/unscaled | restricted | 1e-13 | **136** | 2,031 | 1,152 | 3.1e16 | 1.0e21 | 6e-5 / 0.89 | 6e-5 / 0.99 | 5.9e-5 / 0.97 | 7.3e4 |
| 5,525 | structural/ND/unscaled | restricted | 1e-8 | **318** | 2,020 | 1,068 | 1.6e11 | 2.6e13 | 3.7e-7 / 0.023 | 5.7e-7 / 0.012 | 2.8e-7 / 2.9e-3 | 17 |
| 5,525 | structural/ND/unscaled | none | 1e-13 | **187** | 0 | 2,368 | 4.7e12 | 8.3e15 | 4.2e-5 / 0.96 | 3.8e-5 / 0.95 | 2.8e-5 / 0.97 | 3.1e8 |
| 5,525 | structural/ND/unscaled | none | 1e-8 | **368** | 0 | 2,083 | 6.8e8 | 1.0e11 | 1e-7 / 1.5e-3 | 8.1e-10 / 4.4e-4 | 4e-10 / 3.6e-4 | 0.8 |
| 19,989 | weighted/ND/mc64 | restricted | 1e-13 | 0 | 4,695 | 145 | 93.6 | 98.2 | 2.1e-15 / 1.1e-12 | 1.9e-16 / 5.2e-16 | 1.9e-16 / 7.5e-16 | 1.1e-9 |
| 19,989 | weighted/ND/mc64 | none (static) | 1e-13 | 0 | 0 | 437 | 1.38e3 | 3,710 | 1.5e-14 / 1.2e-13 | 1.5e-16 / 6.3e-16 | 1.5e-16 / 5.4e-16 | 1.7e-9 |
| 19,989 | weighted/AMD/mc64 | restricted | 1e-13 | 0 | 5,148 | 310 | 2.82e3 | 1,640 | 3.7e-14 / 6.8e-13 | 1.9e-16 / 5.8e-16 | 1.5e-16 / 5.4e-16 | 2.0e-9 |
| 19,989 | structural/ND/unscaled | restricted | 1e-13 | **528** | 8,086 | 4,947 | 9.4e16 | 1.2e21 | 2.8e-5 / 0.98 | 1e-5 / 0.98 | 4.6e-6 / 0.99 | 72 |

How to read the table:
- Forward errors of 1e-10 to 1e-9 agree with the pivoted scalar LU and SuperLU forward errors on the same known solution, so they reflect conditioning.
- The float64 residual floor limits BE ref≥1. The scalar LU numbers use real128 residuals; the two are not directly comparable below about 1e-16.
- Emulation peak RSS was 14.5–17.5 GB for n=19,989 (dense 3.2 GB matrix plus copies). Its time (86–92 s per run) is not evidence.

## Conclusions

1. **The static symmetrized structure is viable for both KIN6D Jacobians, and much smaller than pivoted fill.** With MPM + ND, nnz(L+U) is:
   - 0.97 M versus 2.0–8.0 M for the native pivoted LU (all variants) and 2.2 M for SuperLU (n=5,525);
   - 4.94 M versus 20.6–57.5 M native and 15.3 M SuperLU (n=19,989).

   That is 2.3× smaller than SuperLU and 4.2–11.6× smaller than any production native variant. The planned payload is 8.2 MB and 41 MB, against native owned peaks of 43–189 MB and 436 MB–1.41 GB.

   The structure-symmetrization overhead is small: the pattern is 95.7% symmetric. The pivoted fill is dominated by pivoting-induced fill: 1,059–1,583 and 4,023–4,732 row departures with structural matching.
2. **Weighted matching and dual scaling are necessary and sufficient here.**
   - With MPM + MC64 scaling, restricted supernode pivoting needs 0 perturbations (τ=1e-13). Growth is 94–131 with ND. One refinement step reaches componentwise BE ≤ 6.5e-16.
   - Even no pivoting at all works (growth ≤ 1.4e3).
   - Native Hopcroft–Karp matching fails: growth 7e8–9e16, 136–528 perturbations, componentwise BE stuck at 3.6e-4 to 0.99 after refinement (both n).
   - Without scaling, τ=1e-8 perturbs about 180 pivots and fails. The perturbation must be relative to the scaled norm, and τ=1e-13 (the PARDISO default) is right.
3. **Most of the arithmetic is dense.** With the existing SPD relaxed-amalgamation rule, supernodes of width ≥8 carry 96.8–98.1% of the LU flops (≥16: 86–93%). External Schur GEMMs are 77–83% of the flops (74–81% with k≥8), TRSM 16–20%, and diagonal blocks 1.2–2.2%. Relaxed padding costs only 2.3–3.4% entries under ND.

   The largest fronts are small: 266 and 531 rows. Total LU flops are 1.3e8 and 1.07e9. For comparison, the n=19,989 pivoted LU performs over 1e10 scalar scatter entries in the earlier mailuefterl counts.

   So supernodal BLAS-3 is justified, but the case is modest in size. The win comes from fill and memory first, and only then from BLAS-3 throughput.
4. **ND beats AMD** on the matched symmetrized structure, in fill, flops and growth. Postordering changes almost nothing for AMD or METIS output; keep it as a cheap safety step.
5. **Limits.**
   - The numerical evidence is a dense emulation in float64, not the sparse implementation.
   - Two snapshots of one consumer are not a corpus.
   - Mac numbers are unqualified.
   - The `kin6d-5525-nd-order.bin` capture was retained but not used: all orders here were regenerated with this revision's wrappers.

## Hashes (scratch binaries and sources)

```
lu_probe.f90       a2cff1931da48116b683f5a4c3fd59000adc16810077da48ad6bcb6552c44cd4
mpm_probe.f90      8aeee79ed8dc6990fcca9b523f1ac71786b0775acc7677688c9543abc9bc36b7
symb_probe.f90     5dd0e6605dec8085e96897da437671679cf71ba9efb6a18b5fe500aff02d4302
static_probe.py    ea49865e81db93d512e301787a890accb80c2e0a9b37d2fcb398fd34dde046fb
wmatch.py          5c9e55e1f5fb5d4f2b0da2854a2dc0b541d44d217a2a26cf1d7d61e5115e5a99
superlu_nnz.py     29a5d347b8ca06a551c992aa4152dcee511db1d795c84d97fca48d984e218cd5
observer LU copy   f3b975ea0819b12f772d02027da3436b7d3aa5b377016ce54e17db3487b5c466 (orig 785070b2a477dac878c1bcbf4139e7e1f2cee2975bc23a16594b7ef6fc8d1459)
symb_probe (bin)   f1ce9bd12d9b99c18a0c1842bf0574c9d36eaa0475563cba70d8e979409db00d
lu_probe (bin)     95cd32419848776f33c7af54fadaf1877ff546fe68bd39a9077ff406b6e4c12e
mpm_probe (bin)    e42debebd4eb26e390cf92a258409f73df6e18e30bf7efbbb78dc2d4e93885bf
libfortsparse.a    caa9247ac35aff42cb1a9dda63259f4847d47235a83412c0cf6ab9b1f3fed95d
```
